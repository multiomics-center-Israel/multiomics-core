#' RNA-seq Commentary Module
#'
#' Generates AI-powered or data-driven commentary for all pipeline figures.
#' Commentary is saved as JSON and embedded in the final HTML report.
#' Enriched with experiment-level data (top DE genes, PCA variance, groups)
#' for richer, context-aware interpretations.
#'
#' @param de_res      Result from mod_rnaseq_de() — list(dds, tables)
#' @param qc_pre_obj  Result from mod_rnaseq_qc_pre() — list with PCA objects
#' @param config      Full pipeline config
#' @param out_dir     RNA output directory (e.g. .../rna)
#' @return Path to commentary_all.json (character, format = "file")
#' @export
mod_rnaseq_commentary <- function(de_res, qc_pre_obj = NULL, config, out_dir) {

    rna_cfg      <- config$modes$rna
    comment_cfg  <- rna_cfg$commentary %||% list()

    # Skip if commentary is not enabled
    if (!isTRUE(comment_cfg$enabled)) {
        message("Commentary disabled in config (commentary.enabled: false)")
        commentary_dir <- file.path(out_dir, "commentary")
        dir.create(commentary_dir, recursive = TRUE, showWarnings = FALSE)
        out_file <- file.path(commentary_dir, "commentary_all.json")
        jsonlite::write_json(list(), out_file, auto_unbox = TRUE)
        return(out_file)
    }

    # 1. Discover all figures
    figures_tbl <- discover_figures(out_dir, config)

    if (nrow(figures_tbl) == 0) {
        message("No figures found for commentary.")
        commentary_dir <- file.path(out_dir, "commentary")
        dir.create(commentary_dir, recursive = TRUE, showWarnings = FALSE)
        out_file <- file.path(commentary_dir, "commentary_all.json")
        jsonlite::write_json(list(), out_file, auto_unbox = TRUE)
        return(out_file)
    }

    # 2. Load DE summary counts for context
    de_summary <- NULL
    de_counts_file <- file.path(out_dir, "Datasets", "de_summary_counts.tsv")
    if (file.exists(de_counts_file)) {
        de_summary <- read.delim(de_counts_file, stringsAsFactors = FALSE)
    }

    # 3. Build enriched data for richer commentary
    extra_data <- build_extra_data(de_res, qc_pre_obj, config)

    # 4. Generate commentary
    commentary_dir <- file.path(out_dir, "commentary")
    commentary_list <- generate_all_commentary(
        figures_tbl = figures_tbl,
        config      = config,
        de_summary  = de_summary,
        extra_data  = extra_data,
        output_dir  = commentary_dir
    )

    out_file <- file.path(commentary_dir, "commentary_all.json")
    message("Commentary module complete: ", length(commentary_list), " figures processed")
    out_file
}


#' Build enriched data for commentary context
#'
#' Extracts top DE genes per contrast, PCA variance explained,
#' experimental groups, sample counts, and pipeline methods from
#' the pipeline results and config.
#'
#' @param de_res      DE results object
#' @param qc_pre_obj  QC pre results (with PCA variance info)
#' @param config      Full pipeline config
#' @return Named list with top_genes, pca_variance, groups, n_samples, etc.
build_extra_data <- function(de_res, qc_pre_obj = NULL, config) {
    rna_cfg <- config$modes$rna
    extra <- list()

    # --- Experimental groups and sample count ---
    tryCatch({
        group_col <- rna_cfg$filtering$group_col %||% rna_cfg$effects$color
        if (!is.null(de_res$tables) && length(de_res$tables) > 0) {
            # Try to get groups from contrast names
            extra$groups <- names(de_res$tables)
        }
        # Try to get from metadata if available
        if (!is.null(de_res$dds)) {
            meta <- tryCatch(
                as.data.frame(SummarizedExperiment::colData(de_res$dds)),
                error = function(e) NULL
            )
            if (!is.null(meta) && !is.null(group_col) && group_col %in% names(meta)) {
                extra$groups <- as.character(unique(meta[[group_col]]))
                extra$n_samples <- nrow(meta)
            }
        }
    }, error = function(e) {
        message("  Could not extract group info: ", e$message)
    })

    # --- PCA variance explained ---
    tryCatch({
        if (!is.null(qc_pre_obj)) {
            # qc_pre_obj$objects$var_expl contains variance explained per PC
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
    extra$norm_method <- rna_cfg$normalization$method %||% "TMM"
    extra$de_method <- paste0("DESeq2 (", rna_cfg$de$deseq_mode %||% "Wald", ")")
    extra$de_cutoffs <- list(
        p_cutoff  = rna_cfg$de$p_cutoff %||% 0.05,
        fc_cutoff = rna_cfg$de$linear_fc_cutoff %||% 1.5
    )

    # --- Top DE genes per contrast (up to 5 up + 5 down) ---
    extra$top_genes <- list()
    tryCatch({
        if (!is.null(de_res$summary_df)) {
            sdf <- de_res$summary_df

            # Detect gene name column
            gene_col <- NULL
            for (candidate in c("GeneName", "gene_name", "Symbol", "FeatureID")) {
                if (candidate %in% names(sdf)) {
                    gene_col <- candidate
                    break
                }
            }
            if (is.null(gene_col)) gene_col <- names(sdf)[1]

            # Get all contrasts from column names
            padj_cols <- grep("^padj\\.", names(sdf), value = TRUE)
            fc_cols   <- grep("^(log2FoldChange|linearFC)\\.", names(sdf), value = TRUE)

            for (pc in padj_cols) {
                contrast_name <- sub("^padj\\.", "", pc)

                # Find matching FC column
                fc_col <- fc_cols[grepl(contrast_name, fc_cols, fixed = TRUE)]
                if (length(fc_col) == 0) next
                fc_col <- fc_col[1]

                padj_vals <- as.numeric(sdf[[pc]])
                fc_vals   <- as.numeric(sdf[[fc_col]])
                genes     <- as.character(sdf[[gene_col]])

                # Convert linear FC to log2 if needed
                is_linear <- grepl("^linearFC\\.", fc_col)
                log2fc <- if (is_linear) {
                    sign(fc_vals) * log2(pmax(abs(fc_vals), 1e-10))
                } else {
                    fc_vals
                }

                # Filter significant genes
                sig_mask <- !is.na(padj_vals) & padj_vals <= (extra$de_cutoffs$p_cutoff %||% 0.05)
                if (sum(sig_mask) == 0) next

                sig_df <- data.frame(
                    gene   = genes[sig_mask],
                    log2FC = log2fc[sig_mask],
                    padj   = padj_vals[sig_mask],
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
                    message("  Top genes for ", contrast_name, ": ",
                            nrow(top_up), " up + ", nrow(top_dn), " down")
                }
            }
        }
    }, error = function(e) {
        message("  Could not extract top DE genes: ", e$message)
    })

    extra
}
