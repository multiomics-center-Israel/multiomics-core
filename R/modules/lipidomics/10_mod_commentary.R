#' Lipidomics Commentary Module
#'
#' Generates AI-powered or data-driven commentary for all lipidomics figures.
#' Commentary is saved as JSON and embedded in the final HTML report.
#' Enriched with experiment-level data (top DE lipids, PCA variance, groups)
#' for richer, context-aware interpretations.
#'
#' @param de_res      Result from mod_lipidomics_de() — list(summary_df, files, method)
#' @param qc_pre_obj  Result from mod_lipidomics_qc_pre() — list with PCA objects
#' @param config      Full pipeline config
#' @param out_dir     Lipidomics output directory (e.g. .../lipidomics)
#' @param deps        Unused. Accepts the plot-producing analysis results so
#'   targets forces those PNGs to be written before this target scans out_dir;
#'   without it a parallel run can generate commentary before feature
#'   selection / class / biomarker / pathway / clustering / enhanced-QC plots
#'   exist, leaving those report sections without AI interpretation.
#' @return Path to commentary_all.json (character, format = "file")
#' @export
mod_lipidomics_commentary <- function(de_res, qc_pre_obj = NULL, config, out_dir,
                                       deps = NULL) {

    lipid_cfg    <- config$modes$lipidomics
    comment_cfg  <- lipid_cfg$commentary %||% list()

    # Skip if commentary is not enabled
    if (!isTRUE(comment_cfg$enabled)) {
        message("Lipidomics commentary disabled in config (commentary.enabled: false)")
        commentary_dir <- file.path(out_dir, "commentary")
        dir.create(commentary_dir, recursive = TRUE, showWarnings = FALSE)
        out_file <- file.path(commentary_dir, "commentary_all.json")
        jsonlite::write_json(list(), out_file, auto_unbox = TRUE)
        return(out_file)
    }

    # 1. Discover all figures
    figures_tbl <- discover_lipidomics_figures(out_dir, config)

    if (nrow(figures_tbl) == 0) {
        message("No lipidomics figures found for commentary.")
        commentary_dir <- file.path(out_dir, "commentary")
        dir.create(commentary_dir, recursive = TRUE, showWarnings = FALSE)
        out_file <- file.path(commentary_dir, "commentary_all.json")
        jsonlite::write_json(list(), out_file, auto_unbox = TRUE)
        return(out_file)
    }

    # 2. Build enriched data for richer commentary
    extra_data <- build_lipidomics_extra_data(de_res, qc_pre_obj, config)

    # 3. Override config so generate_all_commentary reads our commentary settings
    #    (it hardcodes config$modes$rna$commentary — same workaround as proteomics)
    config_for_commentary <- config
    if (is.null(config_for_commentary$modes$rna))
        config_for_commentary$modes$rna <- list()
    config_for_commentary$modes$rna$commentary <- comment_cfg

    # 4. Generate commentary
    commentary_dir <- file.path(out_dir, "commentary")
    commentary_list <- generate_all_commentary(
        figures_tbl = figures_tbl,
        config      = config_for_commentary,
        de_summary  = NULL,
        extra_data  = extra_data,
        output_dir  = commentary_dir
    )

    out_file <- file.path(commentary_dir, "commentary_all.json")
    message("Lipidomics commentary complete: ", length(commentary_list), " figures processed")
    out_file
}


#' Build enriched data for lipidomics commentary context
#'
#' Extracts top DE lipids per contrast, PCA variance explained,
#' experimental groups, sample counts, and pipeline methods from
#' the pipeline results and config.
#'
#' @param de_res      DE results object
#' @param qc_pre_obj  QC pre results (with PCA variance info)
#' @param config      Full pipeline config
#' @return Named list with top_lipids, pca_variance, groups, n_samples, etc.
build_lipidomics_extra_data <- function(de_res, qc_pre_obj = NULL, config) {
    lipid_cfg <- config$modes$lipidomics
    extra <- list()

    extra$domain <- "lipidomics"

    # --- Experimental groups and sample count ---
    tryCatch({
        group_col <- lipid_cfg$de$condition_column %||%
                     lipid_cfg$effects$color %||% "condition"
        extra$group_col <- group_col

        if (!is.null(de_res$summary_df)) {
            contrasts <- lipid_cfg$de$contrasts
            if (is.list(contrasts)) contrasts <- unlist(contrasts)
            extra$contrasts <- contrasts
        }
    }, error = function(e) {
        message("  Could not extract group info: ", e$message)
    })

    # --- PCA variance explained ---
    tryCatch({
        if (!is.null(qc_pre_obj)) {
            var_expl <- qc_pre_obj$objects$var_expl
            if (is.null(var_expl)) var_expl <- qc_pre_obj$var_expl
            if (!is.null(var_expl) && length(var_expl) >= 2) {
                extra$pca_variance <- as.numeric(var_expl[1:min(3, length(var_expl))])
                message("  PCA variance: PC1=", round(extra$pca_variance[1], 1),
                        "%, PC2=", round(extra$pca_variance[2], 1), "%")
            }
        }
    }, error = function(e) {
        message("  Could not extract PCA variance: ", e$message)
    })

    # --- Normalization and DE methods ---
    norm_info <- lipid_cfg$normalization %||% list()
    norm_parts <- character(0)
    if (!is.null(norm_info$method) && tolower(norm_info$method) != "none")
        norm_parts <- c(norm_parts, norm_info$method)
    if (!is.null(norm_info$transform) && tolower(norm_info$transform) != "none")
        norm_parts <- c(norm_parts, norm_info$transform)
    if (!is.null(norm_info$scaling) && tolower(norm_info$scaling) != "none")
        norm_parts <- c(norm_parts, norm_info$scaling)
    extra$norm_method <- if (length(norm_parts) > 0) paste(norm_parts, collapse = " + ") else "none"

    extra$de_method <- lipid_cfg$de$method %||% "limma"
    extra$de_cutoffs <- list(
        p_cutoff  = lipid_cfg$de$p_cutoff %||% 0.05,
        fc_cutoff = lipid_cfg$de$linear_fc_cutoff %||% 1.5
    )

    # --- Top DE lipids per contrast (up to 5 up + 5 down) ---
    extra$top_genes <- list()
    tryCatch({
        if (!is.null(de_res$summary_df)) {
            sdf <- de_res$summary_df

            # Detect lipid name column
            name_col <- NULL
            for (candidate in c("display_name", "Name", "feature_id")) {
                if (candidate %in% names(sdf)) {
                    name_col <- candidate
                    break
                }
            }
            if (is.null(name_col)) name_col <- names(sdf)[1]

            # Get all contrasts from column names
            pval_cols <- grep("^P\\.Value", names(sdf), value = TRUE)
            fc_cols   <- grep("^logFC", names(sdf), value = TRUE)

            for (pc in pval_cols) {
                contrast_name <- sub("^P\\.Value\\.?", "", pc)
                if (!nzchar(contrast_name)) contrast_name <- "main"

                # Find matching FC column
                fc_col <- fc_cols[grepl(contrast_name, fc_cols, fixed = TRUE)]
                if (length(fc_col) == 0 && length(fc_cols) == 1) fc_col <- fc_cols[1]
                if (length(fc_col) == 0) next
                fc_col <- fc_col[1]

                pval_vals <- as.numeric(sdf[[pc]])
                fc_vals   <- as.numeric(sdf[[fc_col]])
                lipids    <- as.character(sdf[[name_col]])

                # Filter significant lipids
                p_cut <- extra$de_cutoffs$p_cutoff %||% 0.05
                sig_mask <- !is.na(pval_vals) & pval_vals <= p_cut
                if (sum(sig_mask) == 0) next

                sig_df <- data.frame(
                    gene   = lipids[sig_mask],
                    log2FC = fc_vals[sig_mask],
                    padj   = pval_vals[sig_mask],
                    stringsAsFactors = FALSE
                )

                # Top 5 up + top 5 down by absolute log2FC
                up_df <- sig_df[sig_df$log2FC > 0, ]
                dn_df <- sig_df[sig_df$log2FC < 0, ]
                up_df <- up_df[order(-abs(up_df$log2FC)), ]
                dn_df <- dn_df[order(-abs(dn_df$log2FC)), ]
                top_up <- head(up_df, 5)
                top_dn <- head(dn_df, 5)

                top_combined <- rbind(
                    if (nrow(top_up) > 0) cbind(top_up, direction = "up") else NULL,
                    if (nrow(top_dn) > 0) cbind(top_dn, direction = "down") else NULL
                )

                if (!is.null(top_combined) && nrow(top_combined) > 0) {
                    extra$top_genes[[contrast_name]] <- top_combined
                    message("  Top DE lipids for ", contrast_name, ": ",
                            nrow(top_up), " up + ", nrow(top_dn), " down")
                }
            }
        }
    }, error = function(e) {
        message("  Could not extract top DE lipids: ", e$message)
    })

    extra
}
