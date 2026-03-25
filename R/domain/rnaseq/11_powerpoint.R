# R/domain/rnaseq/11_powerpoint.R
#
# Generate a PowerPoint presentation summarizing RNA-seq pipeline results.
# Each slide includes an AI-generated bottom line via Claude Code CLI,
# following the same pattern as the metabolomics 11_powerpoint.R.
#
# Uses shared helpers from R/core/16_pptx_helpers.R:
#   pptx_enabled(), pptx_add_title_slide(), pptx_add_design_slide(),
#   pptx_add_figure_slide(), pptx_add_table_slide(),
#   pptx_fp_*() font functions, PPTX_COL_* constants,
#   generate_slide_bottom_line()
#
# Requires: officer, flextable (optional, for tables)


# ==============================================================================
# POWERPOINT GENERATION
# ==============================================================================

#' Generate an RNA-seq summary PowerPoint
#'
#' @param pre            Preprocessed data list (expr_work, expr_filt, meta, row_data).
#' @param qc_pre_obj     QC pre-DE results list (with $plots, $files).
#' @param de_res         DE results list (with $summary_df, $tables).
#' @param clustering_obj Clustering results list or NULL.
#' @param pathway_res    Pathway/enrichment results list or NULL.
#' @param config         Full pipeline config.
#' @param out_dir        Output directory for this mode.
#' @return Character path to the generated .pptx file, or character(0) if skipped.
#' @export
generate_rnaseq_pptx <- function(pre, qc_pre_obj, de_res, clustering_obj,
                                  pathway_res, config, out_dir) {

    # ---- Guard: check if PowerPoint generation is enabled ----
    if (!pptx_enabled(config, "rna")) {
        message("RNA-seq PowerPoint generation disabled in config.")
        return(character(0))
    }

    if (!requireNamespace("officer", quietly = TRUE)) {
        warning("officer package not available -- skipping PowerPoint generation")
        return(character(0))
    }

    message("Generating RNA-seq PowerPoint presentation...")

    # ---- Config extraction ----
    rna_cfg  <- config$modes$rna
    de_cfg   <- rna_cfg$de %||% list()
    norm_cfg <- rna_cfg$normalization %||% list()
    diag_dir <- file.path(out_dir, "Diagnostic_plots")
    qc_post_dir <- file.path(diag_dir, "QC_post")

    # Project info
    project_name <- config$project$name %||% "RNA-seq Analysis"
    analyst      <- gsub("_", " ", config$project$analyst %||% "")
    date_str     <- format(Sys.Date(), "%B %d, %Y")
    organism     <- rna_cfg$annotation$organism %||% "unknown"

    # Condition / group column
    cond_col <- de_cfg$condition_column %||% rna_cfg$effects$color %||% "Group"
    # Handle array config for multi-color PCA
    if (length(cond_col) > 1) cond_col <- as.character(cond_col[[1]])

    # Group info
    groups <- if (!is.null(pre$meta) && cond_col %in% colnames(pre$meta)) {
        table(pre$meta[[cond_col]])
    } else NULL

    # DE stats
    n_genes   <- if (!is.null(pre$expr_filt)) nrow(pre$expr_filt) else NA
    n_samples <- if (!is.null(pre$expr_filt)) ncol(pre$expr_filt) else NA
    n_de      <- if (!is.null(de_res$summary_df) &&
                      "pass_any_contrast" %in% colnames(de_res$summary_df)) {
        sum(de_res$summary_df$pass_any_contrast %in% c(TRUE, 1), na.rm = TRUE)
    } else 0
    de_method <- de_res$method %||% "DESeq2"

    # Thresholds
    p_cut   <- de_cfg$p_cutoff %||% 0.05
    lfc_cut <- log2(de_cfg$linear_fc_cutoff %||% 1.5)

    # Contrast names from de_res$tables
    contrasts <- names(de_res$tables)

    # Groups string
    groups_str <- if (!is.null(groups)) {
        paste(sprintf("%s (n=%d)", names(groups), as.integer(groups)), collapse = "  vs  ")
    } else "N/A"

    # Gene name resolver: use row_data if available (annotation-based names)
    name_lut <- NULL
    if (!is.null(pre$row_data)) {
        # Try several common gene name columns in order of preference
        name_col <- NULL
        candidates <- c("Gene_Symbol", "gene_symbol", "symbol",
                         "gene_name", "Name", "Gene_Name")
        for (cc in candidates) {
            if (cc %in% colnames(pre$row_data)) {
                name_col <- cc
                break
            }
        }
        if (!is.null(name_col)) {
            # Build ID column
            id_col <- NULL
            id_candidates <- c("FeatureID", "gene_id", "Ensembl_ID", "feature_id")
            for (ic in id_candidates) {
                if (ic %in% colnames(pre$row_data)) {
                    id_col <- ic
                    break
                }
            }
            if (is.null(id_col)) {
                # Fallback: use rownames
                ids <- rownames(pre$row_data)
            } else {
                ids <- as.character(pre$row_data[[id_col]])
            }
            nms <- as.character(pre$row_data[[name_col]])
            name_lut <- stats::setNames(nms, ids)
        }
    }

    resolve <- function(ids) {
        if (!is.null(name_lut)) {
            resolved <- name_lut[ids]
            ifelse(is.na(resolved) | resolved == "", ids, resolved)
        } else {
            ids
        }
    }

    # Footer
    footer_txt <- paste0("multiomics-core | ", analyst, " | ", date_str,
                          "   |   Generated by AI-assisted pipeline")

    # ---- Create presentation ----
    pptx <- officer::read_pptx()

    # ==================================================================
    # SLIDE 1: Title
    # ==================================================================
    title_context <- sprintf(
        "Title slide for RNA-seq analysis. Project: %s. Analyst: %s. Organism: %s. ",
        project_name, analyst, organism)
    if (!is.null(groups)) {
        title_context <- paste0(title_context,
            sprintf("Groups: %s. ", paste(sprintf("%s (n=%d)", names(groups),
                                                   as.integer(groups)), collapse = ", ")))
    }
    title_context <- paste0(title_context,
        sprintf("%s genes, %d samples, %d DE genes (DESeq2).",
                format(n_genes, big.mark = ","), n_samples, n_de))
    title_bl <- generate_slide_bottom_line(title_context, config = config, mode = "rna")

    stats_text <- sprintf("%s genes  |  %d samples  |  %s  |  %d DE genes",
                           format(n_genes, big.mark = ","), n_samples, groups_str, n_de)

    pptx <- pptx_add_title_slide(
        pptx       = pptx,
        main_title = "RNA-seq Analysis",
        project_name = project_name,
        analyst    = analyst,
        date_str   = date_str,
        stats_text = stats_text,
        bl         = title_bl,
        footer     = footer_txt
    )

    # ==================================================================
    # SLIDE 2: Study Design
    # ==================================================================
    norm_method <- norm_cfg$method %||% "TMMlogCPM"

    design_context <- sprintf(
        "Study design. Organism: %s. %d samples, %s genes. Groups: %s. "
        , organism, n_samples, format(n_genes, big.mark = ","), groups_str)
    design_context <- paste0(design_context,
        sprintf("DE method: %s. Normalization: %s.", de_method, norm_method))
    design_bl <- generate_slide_bottom_line(design_context, config = config, mode = "rna")

    design_items <- list(
        list(label = "Organism",        value = organism),
        list(label = "Samples",         value = as.character(n_samples)),
        list(label = "Genes (filtered)", value = format(n_genes, big.mark = ",")),
        list(label = "Groups",          value = groups_str),
        list(label = "DE Method",       value = de_method),
        list(label = "Normalization",   value = norm_method),
        list(label = "P-value cutoff (adj)", value = as.character(p_cut)),
        list(label = "log2(FC) cutoff", value = sprintf("%.2f (linear FC = %.1f)",
                                                         lfc_cut, 2^lfc_cut)),
        list(label = "DE genes",        value = as.character(n_de))
    )

    pptx <- pptx_add_design_slide(
        pptx   = pptx,
        items  = design_items,
        bl     = design_bl,
        footer = footer_txt
    )

    # ==================================================================
    # SLIDE 3: PCA PC1 vs PC2
    # ==================================================================
    pca_png <- file.path(diag_dir, "PCA_PC1.vs.PC2.png")
    if (file.exists(pca_png)) {
        pca_context <- sprintf(
            "PCA plot (PC1 vs PC2) of RNA-seq normalized expression. %d samples, groups: %s.",
            n_samples, groups_str)
        pca_bl <- generate_slide_bottom_line(pca_context, image_path = pca_png,
                                              config = config, mode = "rna")
        pptx <- pptx_add_figure_slide(
            pptx     = pptx,
            title    = "PCA: PC1 vs PC2",
            img      = pca_png,
            subtitle = sprintf("%d samples | %s | %s normalization",
                               n_samples, groups_str, norm_method),
            bl       = pca_bl,
            footer   = footer_txt
        )
    }

    # ==================================================================
    # SLIDE 4: Sample Correlation Heatmap
    # ==================================================================
    cor_png <- file.path(diag_dir, "sample_correlation_heatmap.png")
    if (file.exists(cor_png)) {
        cor_context <- "Sample-to-sample Pearson correlation heatmap after normalization."
        cor_bl <- generate_slide_bottom_line(cor_context, image_path = cor_png,
                                              config = config, mode = "rna")
        pptx <- pptx_add_figure_slide(
            pptx   = pptx,
            title  = "Sample Correlation Heatmap",
            img    = cor_png,
            bl     = cor_bl,
            footer = footer_txt
        )
    }

    # ==================================================================
    # SLIDES: Volcano plot (per contrast)
    # ==================================================================
    for (ctr in contrasts) {
        volcano_png <- file.path(qc_post_dir, paste0("volcano_", ctr, ".png"))
        if (!file.exists(volcano_png)) next

        ctr_tbl <- de_res$tables[[ctr]]
        if (is.null(ctr_tbl) || nrow(ctr_tbl) == 0) next

        # Count up/down regulated genes using log2FoldChange and padj
        n_up <- sum(!is.na(ctr_tbl$padj) & ctr_tbl$padj < p_cut &
                     !is.na(ctr_tbl$log2FoldChange) & ctr_tbl$log2FoldChange >= lfc_cut,
                    na.rm = TRUE)
        n_down <- sum(!is.na(ctr_tbl$padj) & ctr_tbl$padj < p_cut &
                       !is.na(ctr_tbl$log2FoldChange) & ctr_tbl$log2FoldChange <= -lfc_cut,
                      na.rm = TRUE)

        # Top hits by absolute fold change
        sig_idx <- which(!is.na(ctr_tbl$padj) & ctr_tbl$padj < p_cut &
                          abs(ctr_tbl$log2FoldChange) >= lfc_cut)
        top_names <- ""
        if (length(sig_idx) > 0) {
            ordered_idx <- sig_idx[order(abs(ctr_tbl$log2FoldChange[sig_idx]),
                                          decreasing = TRUE)]
            top_ids <- ctr_tbl$FeatureID[utils::head(ordered_idx, 8)]
            top_names <- paste(resolve(top_ids), collapse = ", ")
        }

        vol_stats <- sprintf(
            '{"contrast":"%s","up":%d,"down":%d,"total_tested":%d,"top":"%s"}',
            gsub("_vs_", " vs ", ctr), n_up, n_down, nrow(ctr_tbl), top_names)
        vol_context <- sprintf(
            "Volcano plot: %s. %d up-regulated, %d down-regulated (padj < %.2g, |log2FC| >= %.2f). Top genes: %s.",
            gsub("_vs_", " vs ", ctr), n_up, n_down, p_cut, lfc_cut, top_names)
        vol_bl <- generate_slide_bottom_line(vol_context, stats_json = vol_stats,
                                              image_path = volcano_png, config = config, mode = "rna")

        subtitle <- sprintf("Up: %d  |  Down: %d  |  Total DE: %d  (padj < %s, |log2FC| >= %.2f)",
                             n_up, n_down, n_up + n_down, p_cut, lfc_cut)
        pptx <- pptx_add_figure_slide(
            pptx     = pptx,
            title    = paste("Volcano Plot:", gsub("_vs_", " vs ", ctr)),
            img      = volcano_png,
            subtitle = subtitle,
            bl       = vol_bl,
            footer   = footer_txt
        )
    }

    # ==================================================================
    # SLIDES: MA plot (per contrast)
    # ==================================================================
    for (ctr in contrasts) {
        ma_png <- file.path(qc_post_dir, paste0("ma_", ctr, ".png"))
        if (!file.exists(ma_png)) next

        ma_context <- sprintf("MA plot: %s. log2(FC) vs average expression for RNA-seq.",
                               gsub("_vs_", " vs ", ctr))
        ma_bl <- generate_slide_bottom_line(ma_context, image_path = ma_png,
                                             config = config, mode = "rna")
        pptx <- pptx_add_figure_slide(
            pptx   = pptx,
            title  = paste("MA Plot:", gsub("_vs_", " vs ", ctr)),
            img    = ma_png,
            bl     = ma_bl,
            footer = footer_txt
        )
    }

    # ==================================================================
    # SLIDES: Top DE genes table (per contrast)
    # ==================================================================
    for (ctr in contrasts) {
        ctr_tbl <- de_res$tables[[ctr]]
        if (is.null(ctr_tbl) || nrow(ctr_tbl) == 0) next

        # Filter to significant genes
        sig <- ctr_tbl[!is.na(ctr_tbl$padj) & ctr_tbl$padj < p_cut &
                        abs(ctr_tbl$log2FoldChange) >= lfc_cut, , drop = FALSE]
        sig <- sig[order(sig$padj), ]
        top_n <- min(15, nrow(sig))
        if (top_n == 0) next

        top_tbl <- sig[seq_len(top_n), , drop = FALSE]

        # Build display data frame
        display_df <- data.frame(
            Gene       = resolve(as.character(top_tbl$FeatureID)),
            `log2FC`   = round(top_tbl$log2FoldChange, 2),
            `padj`     = signif(top_tbl$padj, 3),
            check.names = FALSE, stringsAsFactors = FALSE
        )

        tbl_context <- sprintf(
            "Top %d DE genes for %s. Top: %s (log2FC=%.2f, padj=%.2e). Total DE: %d.",
            top_n, gsub("_vs_", " vs ", ctr),
            display_df$Gene[1], display_df$`log2FC`[1],
            display_df$padj[1], nrow(sig))
        tbl_bl <- generate_slide_bottom_line(tbl_context, config = config, mode = "rna")

        subtitle <- sprintf("%d significant genes (padj < %s, |log2FC| >= %.2f)",
                             nrow(sig), p_cut, lfc_cut)

        pptx <- pptx_add_table_slide(
            pptx       = pptx,
            title      = paste("Top DE Genes:", gsub("_vs_", " vs ", ctr)),
            display_df = display_df,
            subtitle   = subtitle,
            bl         = tbl_bl,
            footer     = footer_txt
        )
    }

    # ==================================================================
    # SLIDES: DE Heatmap (per contrast from QC_post)
    # ==================================================================
    for (ctr in contrasts) {
        hm_png <- file.path(qc_post_dir, paste0("heatmap_top25_", ctr, ".png"))
        if (!file.exists(hm_png)) next

        hm_context <- sprintf(
            "Heatmap of top 25 DE genes for %s with hierarchical clustering.",
            gsub("_vs_", " vs ", ctr))
        hm_bl <- generate_slide_bottom_line(hm_context, image_path = hm_png,
                                             config = config, mode = "rna")
        pptx <- pptx_add_figure_slide(
            pptx     = pptx,
            title    = paste("DE Heatmap:", gsub("_vs_", " vs ", ctr)),
            img      = hm_png,
            subtitle = "Z-scored normalized expression | hierarchical clustering",
            bl       = hm_bl,
            footer   = footer_txt
        )
    }

    # Fallback: global DE heatmap if no per-contrast heatmaps were found
    global_hm <- file.path(diag_dir, "heatmap_top_DE.png")
    if (!any(file.exists(file.path(qc_post_dir,
                                    paste0("heatmap_top25_", contrasts, ".png")))) &&
        file.exists(global_hm)) {
        hm_context <- "Heatmap of top DE genes across all contrasts with hierarchical clustering."
        hm_bl <- generate_slide_bottom_line(hm_context, image_path = global_hm,
                                             config = config, mode = "rna")
        pptx <- pptx_add_figure_slide(
            pptx     = pptx,
            title    = "Top DE Genes Heatmap",
            img      = global_hm,
            subtitle = "Z-scored normalized expression | hierarchical clustering",
            bl       = hm_bl,
            footer   = footer_txt
        )
    }

    # ==================================================================
    # SLIDES: Enrichment plots
    # ==================================================================
    enrich_plot_dir <- file.path(out_dir, "Enrichment", "plots")
    if (dir.exists(enrich_plot_dir)) {
        enrich_pngs <- list.files(enrich_plot_dir, pattern = "\\.png$",
                                   full.names = TRUE, recursive = TRUE)
        for (ep in enrich_pngs) {
            ep_label <- gsub("_", " ", tools::file_path_sans_ext(basename(ep)))
            ep_label <- tools::toTitleCase(ep_label)
            ep_context <- sprintf("Pathway enrichment: %s.", ep_label)
            ep_bl <- generate_slide_bottom_line(ep_context, image_path = ep,
                                                 config = config, mode = "rna")
            pptx <- pptx_add_figure_slide(
                pptx   = pptx,
                title  = paste("Enrichment:", ep_label),
                img    = ep,
                bl     = ep_bl,
                footer = footer_txt
            )
        }
    }

    # Also check the parent Enrichment directory for any additional PNGs
    enrich_dir <- file.path(out_dir, "Enrichment")
    if (dir.exists(enrich_dir)) {
        enrich_root_pngs <- list.files(enrich_dir, pattern = "\\.png$",
                                        full.names = TRUE, recursive = FALSE)
        for (ep in enrich_root_pngs) {
            ep_label <- gsub("_", " ", tools::file_path_sans_ext(basename(ep)))
            ep_label <- tools::toTitleCase(ep_label)
            ep_context <- sprintf("Pathway enrichment: %s.", ep_label)
            ep_bl <- generate_slide_bottom_line(ep_context, image_path = ep,
                                                 config = config, mode = "rna")
            pptx <- pptx_add_figure_slide(
                pptx   = pptx,
                title  = paste("Enrichment:", ep_label),
                img    = ep,
                bl     = ep_bl,
                footer = footer_txt
            )
        }
    }

    # ---- Save to Results root (parent of mode dir) ----
    results_root <- dirname(out_dir)
    out_file <- file.path(results_root, "rnaseq_summary.pptx")
    print(pptx, target = out_file)
    message("PowerPoint saved: ", out_file)
    out_file
}
