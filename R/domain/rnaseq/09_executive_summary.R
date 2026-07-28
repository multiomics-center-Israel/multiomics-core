#' Executive Summary Generation for RNA-seq Reports
#'
#' Generates a concise executive summary of the analysis results
#' and saves it as executive_summary.md in the output directory.

# ==============================================================================
# MAIN FUNCTION
# ==============================================================================

#' Generate RNA-seq executive summary
#'
#' @param de_res DE results list (with summary_df)
#' @param pathway_res Pathway/enrichment results
#' @param qc_pre_obj QC results from pre-DE step
#' @param pre Preprocessed data list
#' @param config Full config
#' @param out_dir Output directory
#' @return Path to the generated executive_summary.md file
#' @export
generate_rnaseq_executive_summary <- function(de_res,
                                               pathway_res = NULL,
                                               qc_pre_obj = NULL,
                                               pre = NULL,
                                               config = NULL,
                                               out_dir = ".") {

    message("=== Generating Executive Summary ===")

    summary_points <- character()

    # --- Sample & gene stats ---
    if (!is.null(pre)) {
        n_samples <- ncol(pre$expr_filt)
        n_genes_raw <- if (!is.null(pre$expr_raw)) nrow(pre$expr_raw) else NA
        n_genes_filt <- nrow(pre$expr_filt)

        summary_points <- c(summary_points,
            sprintf("Analyzed **%d samples** with **%s genes** after filtering (from %s raw)",
                    n_samples,
                    format(n_genes_filt, big.mark = ","),
                    format(n_genes_raw, big.mark = ",")))
    }

    # --- DE results ---
    de_stats <- get_de_summary_stats(de_res, rna_cfg = config$modes$rna)
    if (!is.null(de_stats)) {
        for (stat in de_stats) {
            summary_points <- c(summary_points, stat$text)
        }

        top_genes <- get_top_de_genes(de_res, n = 10)
        if (length(top_genes) > 0) {
            summary_points <- c(summary_points,
                sprintf("  Top DE genes: %s", paste(top_genes, collapse = ", ")))
        }
    }

    # --- Pathway highlights ---
    pw_highlights <- get_pathway_highlights(pathway_res, n = 5)
    if (length(pw_highlights) > 0) {
        summary_points <- c(summary_points,
            sprintf("Top enriched pathways: %s", paste(pw_highlights, collapse = "; ")))
    }

    # --- QC flags ---
    if (!is.null(qc_pre_obj) && is.list(qc_pre_obj)) {
        if (!is.null(qc_pre_obj$outliers)) {
            n_out <- sum(qc_pre_obj$outliers$is_outlier, na.rm = TRUE)
            if (n_out > 0) {
                summary_points <- c(summary_points,
                    sprintf("**%d potential outlier sample(s)** detected by PCA", n_out))
            }
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

#' Extract DE summary statistics from multiomics-core de_res
#' @noRd
get_de_summary_stats <- function(de_res, rna_cfg = NULL) {
    if (is.null(de_res) || is.null(de_res$summary_df)) return(NULL)

    df <- de_res$summary_df

    # Use pipeline pass columns as the single source of truth for DE counts.
    # RNA-seq pattern: {contrast}_pass (TRUE/FALSE)
    pass_cols <- grep("_pass$", colnames(df), value = TRUE)
    pass_cols <- setdiff(pass_cols, "pass_any_contrast")

    if (length(pass_cols) > 0) {
        stats <- list()
        for (i in seq_along(pass_cols)) {
            pcol <- pass_cols[i]
            cn <- sub("_pass$", "", pcol)
            is_sig <- !is.na(df[[pcol]]) & df[[pcol]] %in% c(TRUE, 1)

            # Find FC column for direction
            fc_col <- paste0("log2FoldChange.", cn)
            if (!(fc_col %in% colnames(df))) {
                fc_col <- paste0("linearFC.", cn)
            }
            if (fc_col %in% colnames(df)) {
                fc_vals <- as.numeric(df[[fc_col]])
                n_up   <- sum(is_sig & fc_vals > 0, na.rm = TRUE)
                n_down <- sum(is_sig & fc_vals < 0, na.rm = TRUE)
                text <- sprintf("**%s**: %d significant genes (%d up, %d down)",
                                cn, n_up + n_down, n_up, n_down)
            } else {
                n_sig <- sum(is_sig)
                text <- sprintf("**%s**: %d significant genes", cn, n_sig)
            }
            stats[[i]] <- list(contrast = cn, n_sig = sum(is_sig), text = text)
        }
        return(stats)
    }

    # Fallback: recompute from padj + LFC using config thresholds
    padj_cols <- grep("^padj\\.", colnames(df), value = TRUE)
    if (length(padj_cols) == 0) return(NULL)

    padj_cut <- if (!is.null(rna_cfg)) rna_cfg$de$p_cutoff %||% 0.05 else 0.05
    lfc_cut  <- if (!is.null(rna_cfg)) log2(rna_cfg$de$linear_fc_cutoff %||% 1.5) else log2(1.5)

    stats <- list()
    for (i in seq_along(padj_cols)) {
        pcol <- padj_cols[i]
        cn <- sub("^padj\\.", "", pcol)
        lcol <- paste0("log2FoldChange.", cn)
        if (!(lcol %in% colnames(df))) next

        sig <- !is.na(df[[pcol]]) & df[[pcol]] <= padj_cut &
               !is.na(df[[lcol]]) & abs(df[[lcol]]) >= lfc_cut
        n_up   <- sum(sig & df[[lcol]] > 0, na.rm = TRUE)
        n_down <- sum(sig & df[[lcol]] < 0, na.rm = TRUE)
        n_sig  <- n_up + n_down
        text <- sprintf("**%s**: %d significant genes (%d up, %d down)",
                        cn, n_sig, n_up, n_down)
        stats[[i]] <- list(contrast = cn, n_sig = n_sig, text = text)
    }
    stats
}

#' Get top DE genes by significance and fold change
#' @noRd
get_top_de_genes <- function(de_res, n = 10) {
    if (is.null(de_res) || is.null(de_res$summary_df)) return(character())

    df <- de_res$summary_df

    # Use pass_any_contrast if available
    pass_col <- intersect(c("pass_any_contrast"), colnames(df))
    id_col <- intersect(c("FeatureID", "gene_id", "gene_symbol"), colnames(df))

    if (length(id_col) == 0) return(character())
    id_col <- id_col[1]

    if (length(pass_col) > 0) {
        sig_df <- df[df[[pass_col[1]]] == TRUE, , drop = FALSE]
    } else {
        padj_cols <- grep("^padj\\.", colnames(df), value = TRUE)
        if (length(padj_cols) == 0) return(character())
        sig_df <- df[!is.na(df[[padj_cols[1]]]) & df[[padj_cols[1]]] < 0.05, , drop = FALSE]
    }

    if (nrow(sig_df) == 0) return(character())

    # Sort by minimum padj
    padj_cols <- grep("^padj\\.", colnames(sig_df), value = TRUE)
    if (length(padj_cols) > 0) {
        min_padj <- apply(sig_df[, padj_cols, drop = FALSE], 1, min, na.rm = TRUE)
        sig_df <- sig_df[order(min_padj), ]
    }

    head(sig_df[[id_col]], n)
}

#' Get pathway highlights from enrichment results
#' @noRd
get_pathway_highlights <- function(pathway_res, n = 5) {
    if (is.null(pathway_res)) return(character())

    pathways <- character()

    # Handle list of results (per contrast/database)
    if (is.data.frame(pathway_res)) {
        name_col <- intersect(c("pathway", "Description", "term", "name"), colnames(pathway_res))
        if (length(name_col) > 0) {
            pathways <- head(pathway_res[[name_col[1]]], n)
        }
    } else if (is.list(pathway_res)) {
        for (res in pathway_res) {
            if (is.data.frame(res)) {
                name_col <- intersect(c("pathway", "Description", "term", "name"), colnames(res))
                if (length(name_col) > 0) {
                    pathways <- c(pathways, head(res[[name_col[1]]], 3))
                }
            }
        }
        pathways <- head(unique(pathways), n)
    }

    # Clean pathway names
    pathways <- gsub("^GO:[0-9]+\\s*", "", pathways)
    pathways <- vapply(pathways, function(x) {
        if (nchar(x) > 60) paste0(substr(x, 1, 57), "...") else x
    }, character(1))

    pathways
}

#' Format summary points into markdown
#' @noRd
format_executive_summary <- function(summary_points) {
    if (length(summary_points) == 0) {
        return("No summary available.\n")
    }

    formatted <- ""

    for (point in summary_points) {
        if (startsWith(point, "  ")) {
            formatted <- paste0(formatted, "  - ", trimws(point), "\n")
        } else {
            formatted <- paste0(formatted, "- ", point, "\n")
        }
    }

    formatted
}
