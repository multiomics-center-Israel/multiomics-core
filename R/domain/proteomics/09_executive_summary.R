#' Proteomics Executive Summary Generation
#'
#' Generates a concise executive summary of the proteomics analysis results
#' and saves it as executive_summary.md in the output directory.
#' Adapted from R/domain/rnaseq/09_executive_summary.R.

# ==============================================================================
# MAIN FUNCTION
# ==============================================================================

#' Generate proteomics executive summary
#'
#' @param de_res       DE results list (with summary_df)
#' @param pathway_res  Pathway results (optional)
#' @param qc_pre_obj   QC pre-DE results (with objects$outlier_res etc.)
#' @param pre          Preprocessed data list
#' @param config       Full config
#' @param out_dir      Output directory
#' @param ppi_res      PPI results (optional)
#' @param adv_stats    Advanced stats results (optional)
#' @return Path to the generated executive_summary.md file
#' @export
generate_proteomics_executive_summary <- function(de_res,
                                                    pathway_res = NULL,
                                                    qc_pre_obj = NULL,
                                                    pre = NULL,
                                                    config = NULL,
                                                    out_dir = ".",
                                                    ppi_res = NULL,
                                                    adv_stats = NULL) {

    message("=== Generating Proteomics Executive Summary ===")

    cfg <- config$modes$proteomics
    summary_points <- character()

    # --- Sample & protein stats ---
    if (!is.null(pre)) {
        n_samples <- ncol(pre$expr_imp_single)
        n_proteins_raw <- if (!is.null(pre$expr_raw)) nrow(pre$expr_raw) else NA
        n_proteins_filt <- nrow(pre$expr_imp_single)

        summary_points <- c(summary_points,
            sprintf("Analyzed **%d samples** with **%s proteins** after filtering (from %s raw)",
                    n_samples,
                    format(n_proteins_filt, big.mark = ","),
                    format(n_proteins_raw, big.mark = ",")))

        # Imputation
        imp_method <- cfg$imputation$method %||% "none"
        if (tolower(imp_method) != "none") {
            summary_points <- c(summary_points,
                sprintf("Missing values imputed using **%s** method", imp_method))
        }

        # Missingness
        if (!is.null(pre$expr_raw)) {
            total_vals <- prod(dim(pre$expr_raw))
            if (total_vals > 0) {
                pct_missing <- round(100 * sum(is.na(pre$expr_raw)) / total_vals, 1)
                summary_points <- c(summary_points,
                    sprintf("Missingness in raw data: **%.1f%%**", pct_missing))
            }
        }
    }

    # --- DE results ---
    if (!is.null(de_res) && !is.null(de_res$summary_df)) {
        de_stats <- get_de_summary_stats_proteomics(de_res$summary_df, cfg)
        for (stat in de_stats) {
            summary_points <- c(summary_points, stat$text)
        }

        top_proteins <- get_top_de_proteins(de_res$summary_df, cfg, n = 10)
        if (length(top_proteins) > 0) {
            summary_points <- c(summary_points,
                sprintf("  Top DE proteins: %s", paste(top_proteins, collapse = ", ")))
        }
    }

    # --- Pathway highlights ---
    pw_highlights <- get_pathway_highlights(pathway_res, n = 5)
    if (length(pw_highlights) > 0) {
        summary_points <- c(summary_points,
            sprintf("Top enriched pathways: %s", paste(pw_highlights, collapse = "; ")))
    }

    # --- PPI summary ---
    if (!is.null(ppi_res) && !is.null(ppi_res$summary)) {
        s <- ppi_res$summary
        summary_points <- c(summary_points,
            sprintf("PPI network: **%d nodes**, **%d edges**, %d communities",
                    s$n_nodes %||% 0, s$n_edges %||% 0, s$n_communities %||% 0))
        if (!is.null(s$top_hubs) && length(s$top_hubs) > 0) {
            summary_points <- c(summary_points,
                sprintf("  Top hub proteins: %s", paste(head(s$top_hubs, 5), collapse = ", ")))
        }
    }

    # --- Advanced stats summary ---
    if (!is.null(adv_stats) && !is.null(adv_stats$summary)) {
        as_s <- adv_stats$summary
        if (!is.null(as_s$current_power)) {
            summary_points <- c(summary_points,
                sprintf("Statistical power at current sample size: **%.0f%%**",
                        as_s$current_power * 100))
        }
        if (!is.null(as_s$n_significant_effects)) {
            summary_points <- c(summary_points,
                sprintf("Proteins with significant effect sizes: **%d**",
                        as_s$n_significant_effects))
        }
    }

    # --- Outlier summary ---
    if (!is.null(qc_pre_obj) && !is.null(qc_pre_obj$objects$outlier_res)) {
        n_out <- qc_pre_obj$objects$outlier_res$n_outliers %||% 0
        if (n_out > 0) {
            summary_points <- c(summary_points,
                sprintf("**%d potential outlier sample(s)** detected", n_out))
        }
    }

    # Format and save
    formatted <- format_executive_summary(summary_points)

    summary_path <- file.path(out_dir, "executive_summary.md")
    writeLines(formatted, summary_path)

    message("Executive summary saved to: ", summary_path)
    summary_path
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Extract DE summary stats from proteomics summary_df
#' @noRd
get_de_summary_stats_proteomics <- function(summary_df, cfg) {
    # Use pipeline pass columns as the single source of truth for DE counts.
    # This ensures executive summary, report, interactive panel, and pipeline
    # summary all show the same numbers.
    pass_cols <- grep("^pass\\.imputs\\.", colnames(summary_df), value = TRUE)

    # Fallback: recompute from padj+LFC if no pass columns
    if (length(pass_cols) == 0) {
        padj_cols <- grep("^padj\\.imputs\\.", colnames(summary_df), value = TRUE)
        if (length(padj_cols) == 0) return(list())

        padj_cut <- cfg$de$p_cutoff %||% 0.05
        lfc_cut  <- cfg$de$linear_fc_cutoff %||% 1.5
        stats <- list()
        for (pcol in padj_cols) {
            cn <- sub("^padj\\.imputs\\.", "", pcol)
            fc_col <- paste0("linearFC.imputs.", cn)
            if (!(fc_col %in% colnames(summary_df))) next
            lfc_vals <- as.numeric(summary_df[[fc_col]])
            log2fc   <- signed_fc_to_log2(lfc_vals)
            is_sig   <- !is.na(as.numeric(summary_df[[pcol]])) &
                        as.numeric(summary_df[[pcol]]) <= padj_cut &
                        abs(log2fc) >= log2(lfc_cut)
            n_up   <- sum(is_sig & log2fc > 0, na.rm = TRUE)
            n_down <- sum(is_sig & log2fc < 0, na.rm = TRUE)
            stats[[cn]] <- list(contrast = cn, n_sig = n_up + n_down,
                text = sprintf("**%s**: %d significant proteins (%d up, %d down)",
                               cn, n_up + n_down, n_up, n_down))
        }
        return(stats)
    }

    stats <- list()
    for (pcol in pass_cols) {
        cn <- sub("^pass\\.imputs\\.", "", pcol)
        fc_col <- paste0("linearFC.imputs.", cn)
        if (!(fc_col %in% colnames(summary_df))) next

        is_sig <- !is.na(summary_df[[pcol]]) & summary_df[[pcol]] == 1
        lfc_vals <- as.numeric(summary_df[[fc_col]])
        n_up   <- sum(is_sig & lfc_vals > 0, na.rm = TRUE)
        n_down <- sum(is_sig & lfc_vals < 0, na.rm = TRUE)
        n_sig  <- n_up + n_down

        stats[[cn]] <- list(
            contrast = cn, n_sig = n_sig,
            text = sprintf("**%s**: %d significant proteins (%d up, %d down)",
                           cn, n_sig, n_up, n_down)
        )
    }
    stats
}

#' Get top DE proteins by significance
#' @noRd
get_top_de_proteins <- function(summary_df, cfg, n = 10) {
    id_col <- cfg$de_table$id_col %||% "FeatureID"
    if (!(id_col %in% colnames(summary_df))) return(character())

    pass_col <- intersect("pass_any_contrast", colnames(summary_df))
    if (length(pass_col) > 0) {
        sig_df <- summary_df[summary_df[[pass_col[1]]] == TRUE, , drop = FALSE]
    } else {
        padj_cols <- grep("^padj\\.imputs\\.", colnames(summary_df), value = TRUE)
        if (length(padj_cols) == 0) return(character())
        min_padj <- apply(summary_df[, padj_cols, drop = FALSE], 1, min, na.rm = TRUE)
        sig_df <- summary_df[min_padj < 0.05, , drop = FALSE]
    }

    if (nrow(sig_df) == 0) return(character())

    padj_cols <- grep("^padj\\.imputs\\.", colnames(sig_df), value = TRUE)
    if (length(padj_cols) > 0) {
        min_padj <- apply(sig_df[, padj_cols, drop = FALSE], 1, min, na.rm = TRUE)
        sig_df <- sig_df[order(min_padj), ]
    }

    head(sig_df[[id_col]], n)
}
