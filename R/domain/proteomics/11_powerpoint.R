# R/domain/proteomics/11_powerpoint.R
#
# Generate a PowerPoint presentation summarizing proteomics pipeline results.
# Each slide includes an AI-generated bottom line via Claude Code CLI,
# following the same pattern as the metabolomics 11_powerpoint.R.
#
# Uses shared helpers from R/core/16_pptx_helpers.R:
#   generate_slide_bottom_line(), pptx_add_title_slide(),
#   pptx_add_design_slide(), pptx_add_figure_slide(),
#   pptx_add_table_slide(), pptx_enabled(),
#   pptx_fp_*() font functions, PPTX_COL_* constants
#
# Requires: officer, flextable (optional, for tables)


# ==============================================================================
# POWERPOINT GENERATION
# ==============================================================================

#' Generate a proteomics summary PowerPoint
#'
#' @param pre         Preprocessed data list (expr_raw, expr_imp_single, meta,
#'                    row_data).
#' @param qc_res      QC results list (with $plots and $files).
#' @param de_res      DE results list (with $summary_df, $de_tables,
#'                    $runs_de_tables).
#' @param pathway_res Pathway / enrichment results list or NULL.
#' @param ppi_res     PPI network analysis results or NULL.
#' @param config      Full pipeline config.
#' @param out_dir     Output directory for the proteomics mode.
#' @return Character path to the generated .pptx file, or character(0) if
#'         generation is disabled or officer is unavailable.
#' @export
generate_proteomics_pptx <- function(pre, qc_res, de_res, pathway_res,
                                      ppi_res, config, out_dir) {

    # ------------------------------------------------------------------
    # Guard: check if PPTX generation is enabled
    # ------------------------------------------------------------------
    if (!pptx_enabled(config, "proteomics")) {
        message("Proteomics PowerPoint generation is disabled in config.")
        return(character(0))
    }

    if (!requireNamespace("officer", quietly = TRUE)) {
        warning("officer package not available -- skipping PowerPoint generation")
        return(character(0))
    }

    message("Generating proteomics PowerPoint presentation...")

    # ------------------------------------------------------------------
    # Config extraction
    # ------------------------------------------------------------------
    prot_cfg  <- config$modes$proteomics
    de_cfg    <- prot_cfg$de %||% list()
    imp_cfg   <- prot_cfg$imputation %||% list()
    norm_cfg  <- prot_cfg$normalization %||% list()
    diag_dir  <- file.path(out_dir, "Diagnostic_plots")
    qc_post   <- file.path(diag_dir, "QC_post")

    # Project info
    project_name <- config$project$name %||% "Proteomics Analysis"
    analyst      <- gsub("_", " ", config$project$analyst %||% "")
    date_str     <- format(Sys.Date(), "%B %d, %Y")
    cond_col     <- de_cfg$condition_column %||% prot_cfg$effects$color %||% "Condition"

    # Group info
    groups <- if (!is.null(pre$meta) && cond_col %in% colnames(pre$meta)) {
        table(pre$meta[[cond_col]])
    } else NULL

    groups_str <- if (!is.null(groups)) {
        paste(sprintf("%s (n=%d)", names(groups), as.integer(groups)),
              collapse = "  vs  ")
    } else "N/A"

    # Protein counts
    n_proteins_raw <- if (!is.null(pre$expr_raw)) nrow(pre$expr_raw) else NA
    n_proteins     <- if (!is.null(pre$expr_imp_single)) {
        nrow(pre$expr_imp_single)
    } else NA
    n_samples      <- if (!is.null(pre$expr_imp_single)) {
        ncol(pre$expr_imp_single)
    } else NA

    # DE stats
    n_de <- if (!is.null(de_res$summary_df) &&
                "pass_any_contrast" %in% colnames(de_res$summary_df)) {
        sum(de_res$summary_df$pass_any_contrast == 1, na.rm = TRUE)
    } else 0
    de_method <- de_res$method %||% de_cfg$method %||% "limma"

    # Thresholds
    p_cut      <- de_cfg$p_cutoff %||% 0.05
    fc_lin_cut <- de_cfg$linear_fc_cutoff %||% 1.5
    lfc_cut    <- log2(fc_lin_cut)

    # Imputation info
    imp_method <- imp_cfg$method %||% "MinProb"

    # Missingness
    pct_missing <- if (!is.null(pre$expr_raw)) {
        total_vals <- prod(dim(pre$expr_raw))
        if (total_vals > 0) round(100 * sum(is.na(pre$expr_raw)) / total_vals, 1) else NA
    } else NA

    # ------------------------------------------------------------------
    # Contrast extraction from summary_df
    # ------------------------------------------------------------------
    contrasts <- character(0)
    if (!is.null(de_res$summary_df)) {
        padj_cols  <- grep("^padj\\.imputs\\.", colnames(de_res$summary_df),
                           value = TRUE)
        contrasts  <- sub("^padj\\.imputs\\.", "", padj_cols)
    }

    # ------------------------------------------------------------------
    # Protein name resolver
    # ------------------------------------------------------------------
    resolve_name <- function(summary_df, idx) {
        # Try Genes, then Protein.Names, then FeatureID
        if ("Genes" %in% colnames(summary_df)) {
            nm <- as.character(summary_df$Genes[idx])
            if (!is.na(nm) && nzchar(nm) && nm != "NA") return(nm)
        }
        if ("Protein.Names" %in% colnames(summary_df)) {
            nm <- as.character(summary_df$Protein.Names[idx])
            if (!is.na(nm) && nzchar(nm) && nm != "NA") return(nm)
        }
        id_col <- prot_cfg$de_table$id_col %||% "FeatureID"
        if (id_col %in% colnames(summary_df)) {
            return(as.character(summary_df[[id_col]][idx]))
        }
        return(as.character(idx))
    }

    resolve_names <- function(summary_df, indices) {
        vapply(indices, function(i) resolve_name(summary_df, i),
               FUN.VALUE = character(1))
    }

    # ------------------------------------------------------------------
    # Footer
    # ------------------------------------------------------------------
    footer_txt <- paste0("multiomics-core | ", analyst, " | ", date_str,
                          "   |   Generated by AI-assisted pipeline")

    # ---- Create presentation ----
    pptx <- officer::read_pptx()

    # ==================================================================
    # SLIDE 1: Title
    # ==================================================================
    title_context <- sprintf(
        "Title slide for proteomics analysis. Project: %s. Analyst: %s. ",
        project_name, analyst)
    if (!is.null(groups)) {
        title_context <- paste0(title_context,
            sprintf("Groups: %s. ", groups_str))
    }
    title_context <- paste0(title_context,
        sprintf("%d proteins (from %s raw), %d samples, %d significant DE proteins. ",
                n_proteins, format(n_proteins_raw, big.mark = ","),
                n_samples, n_de),
        sprintf("DE method: %s with multiple imputation. Imputation: %s.",
                de_method, imp_method))

    title_bl <- generate_slide_bottom_line(title_context, config = config, mode = "proteomics")

    stats_text <- sprintf(
        "%d proteins  |  %d samples  |  %s  |  %d DE proteins",
        n_proteins, n_samples, groups_str, n_de)

    pptx <- pptx_add_title_slide(pptx,
        main_title   = "Proteomics Analysis",
        project_name = project_name,
        analyst      = analyst,
        date_str     = date_str,
        stats_text   = stats_text,
        bl           = title_bl,
        footer       = footer_txt)

    # ==================================================================
    # SLIDE 2: Study Design
    # ==================================================================
    design_context <- sprintf(
        "Study design. %d samples, %d proteins (from %s raw). Groups: %s. ",
        n_samples, n_proteins, format(n_proteins_raw, big.mark = ","),
        groups_str)
    design_context <- paste0(design_context,
        sprintf("DE method: %s. Imputation: %s. ", de_method, imp_method),
        sprintf("Missing data: %.1f%%. p-cutoff: %s, linear FC cutoff: %s.",
                pct_missing %||% NA, p_cut, fc_lin_cut))
    design_bl <- generate_slide_bottom_line(design_context, config = config, mode = "proteomics")

    design_items <- list(
        list(label = "Samples",           value = as.character(n_samples)),
        list(label = "Proteins (filtered)", value = format(n_proteins, big.mark = ",")),
        list(label = "Proteins (raw)",     value = format(n_proteins_raw, big.mark = ",")),
        list(label = "Groups",             value = groups_str),
        list(label = "DE Method",          value = paste0(de_method, " + multiple imputation")),
        list(label = "Imputation",         value = imp_method),
        list(label = "Missing data",       value = sprintf("%.1f%%", pct_missing %||% NA)),
        list(label = "P-value cutoff",     value = as.character(p_cut)),
        list(label = "Linear FC cutoff",   value = as.character(fc_lin_cut)),
        list(label = "DE proteins",        value = as.character(n_de))
    )

    pptx <- pptx_add_design_slide(pptx, items = design_items,
                                    bl = design_bl, footer = footer_txt)

    # ==================================================================
    # SLIDE 3: PCA PC1 vs PC2
    # ==================================================================
    pca_png <- file.path(diag_dir, "PCA_PC1.vs.PC2.png")
    if (file.exists(pca_png)) {
        pca_context <- paste0(
            "PCA plot (PC1 vs PC2) of proteomics expression data. ",
            sprintf("%d samples, groups: %s. %d proteins.",
                    n_samples, groups_str, n_proteins))
        pca_bl <- generate_slide_bottom_line(pca_context,
                                              image_path = pca_png,
                                              config = config, mode = "proteomics")
        pptx <- pptx_add_figure_slide(pptx,
            title    = "PCA: PC1 vs PC2",
            img      = pca_png,
            subtitle = sprintf("%d samples | %s | %d proteins",
                               n_samples, groups_str, n_proteins),
            bl       = pca_bl,
            footer   = footer_txt)
    }

    # ==================================================================
    # SLIDE 4: UMAP (if exists)
    # ==================================================================
    umap_png <- file.path(diag_dir, "UMAP.png")
    if (file.exists(umap_png)) {
        umap_context <- paste0(
            "UMAP dimensionality reduction of proteomics samples. ",
            sprintf("%d samples, groups: %s.", n_samples, groups_str))
        umap_bl <- generate_slide_bottom_line(umap_context,
                                               image_path = umap_png,
                                               config = config, mode = "proteomics")
        pptx <- pptx_add_figure_slide(pptx,
            title    = "UMAP: Protein Expression",
            img      = umap_png,
            subtitle = sprintf("%d samples | %s", n_samples, groups_str),
            bl       = umap_bl,
            footer   = footer_txt)
    }

    # ==================================================================
    # SLIDE 5: Sample Correlation Heatmap
    # ==================================================================
    cor_png <- file.path(diag_dir, "sample_correlation_heatmap.png")
    if (file.exists(cor_png)) {
        cor_context <- paste0(
            "Sample-to-sample Pearson correlation heatmap after imputation. ",
            sprintf("%d samples, groups: %s.", n_samples, groups_str))
        cor_bl <- generate_slide_bottom_line(cor_context,
                                              image_path = cor_png,
                                              config = config, mode = "proteomics")
        pptx <- pptx_add_figure_slide(pptx,
            title    = "Sample Correlation Heatmap",
            img      = cor_png,
            subtitle = sprintf("%d samples | Pearson correlation", n_samples),
            bl       = cor_bl,
            footer   = footer_txt)
    }

    # ==================================================================
    # SLIDE 6+: Volcano Plots (per contrast)
    # ==================================================================
    for (ctr in contrasts) {
        # Try both summary and table1 naming conventions
        volcano_png <- file.path(qc_post,
                                  sprintf("volcano_%s_summary.png", ctr))
        if (!file.exists(volcano_png)) {
            volcano_png <- file.path(qc_post,
                                      sprintf("volcano_%s_table1.png", ctr))
        }
        if (!file.exists(volcano_png)) {
            # Fallback: scan for any volcano file matching this contrast
            candidates <- list.files(qc_post,
                pattern = sprintf("^volcano_%s.*\\.png$", ctr),
                full.names = TRUE)
            if (length(candidates) > 0) volcano_png <- candidates[1]
        }
        if (!file.exists(volcano_png)) next

        # Compute per-contrast stats from summary_df
        padj_vals <- as.numeric(de_res$summary_df[[paste0("padj.imputs.", ctr)]])
        lfc_vals  <- as.numeric(de_res$summary_df[[paste0("linearFC.imputs.", ctr)]])
        log2fc    <- signed_fc_to_log2(lfc_vals)

        is_sig <- !is.na(padj_vals) & padj_vals <= p_cut &
                  !is.na(log2fc) & abs(log2fc) >= lfc_cut
        n_up   <- sum(is_sig & log2fc > 0, na.rm = TRUE)
        n_down <- sum(is_sig & log2fc < 0, na.rm = TRUE)

        # Top hits by absolute fold change
        sig_idx <- which(is_sig)
        top_names <- ""
        if (length(sig_idx) > 0) {
            ordered_idx <- sig_idx[order(abs(log2fc[sig_idx]),
                                          decreasing = TRUE)]
            top_idx <- utils::head(ordered_idx, 8)
            top_names <- paste(resolve_names(de_res$summary_df, top_idx),
                               collapse = ", ")
        }

        vol_stats <- sprintf(
            '{"contrast":"%s","up":%d,"down":%d,"total_tested":%d,"top":"%s"}',
            gsub("_vs_", " vs ", ctr), n_up, n_down,
            sum(!is.na(padj_vals)), top_names)
        vol_context <- sprintf(
            "Volcano plot for proteomics: %s. %d up, %d down ",
            gsub("_vs_", " vs ", ctr), n_up, n_down)
        vol_context <- paste0(vol_context,
            sprintf("(adj.p < %.2g, |log2FC| >= %.2f). Top: %s.",
                    p_cut, lfc_cut, top_names))
        vol_bl <- generate_slide_bottom_line(vol_context,
                                              stats_json = vol_stats,
                                              image_path = volcano_png,
                                              config = config, mode = "proteomics")

        subtitle <- sprintf(
            "Up: %d  |  Down: %d  |  Total DE: %d  (adj.p < %s, linear FC >= %s)",
            n_up, n_down, n_up + n_down, p_cut, fc_lin_cut)

        pptx <- pptx_add_figure_slide(pptx,
            title    = paste("Volcano Plot:", gsub("_vs_", " vs ", ctr)),
            img      = volcano_png,
            subtitle = subtitle,
            bl       = vol_bl,
            footer   = footer_txt)
    }

    # ==================================================================
    # SLIDE: MA Plots (per contrast)
    # ==================================================================
    for (ctr in contrasts) {
        ma_png <- file.path(qc_post,
                             sprintf("ma_%s_summary.png", ctr))
        if (!file.exists(ma_png)) {
            ma_png <- file.path(qc_post,
                                 sprintf("ma_%s_table1.png", ctr))
        }
        if (!file.exists(ma_png)) {
            candidates <- list.files(qc_post,
                pattern = sprintf("^ma_%s.*\\.png$", ctr),
                full.names = TRUE)
            if (length(candidates) > 0) ma_png <- candidates[1]
        }
        if (!file.exists(ma_png)) next

        ma_context <- sprintf(
            "MA plot for proteomics: %s. log2(FC) vs average expression.",
            gsub("_vs_", " vs ", ctr))
        ma_bl <- generate_slide_bottom_line(ma_context,
                                             image_path = ma_png,
                                             config = config, mode = "proteomics")
        pptx <- pptx_add_figure_slide(pptx,
            title    = paste("MA Plot:", gsub("_vs_", " vs ", ctr)),
            img      = ma_png,
            subtitle = sprintf("log2(FC) vs Mean Expression | %s",
                               gsub("_vs_", " vs ", ctr)),
            bl       = ma_bl,
            footer   = footer_txt)
    }

    # ==================================================================
    # SLIDE: Top DE Proteins Table (per contrast)
    # ==================================================================
    for (ctr in contrasts) {
        sdf <- de_res$summary_df
        if (is.null(sdf) || nrow(sdf) == 0) next

        padj_vals <- as.numeric(sdf[[paste0("padj.imputs.", ctr)]])
        lfc_vals  <- as.numeric(sdf[[paste0("linearFC.imputs.", ctr)]])
        log2fc    <- signed_fc_to_log2(lfc_vals)

        is_sig <- !is.na(padj_vals) & padj_vals <= p_cut &
                  !is.na(log2fc) & abs(log2fc) >= lfc_cut
        n_sig  <- sum(is_sig, na.rm = TRUE)
        if (n_sig == 0) next

        # Order by adjusted p-value
        sig_idx <- which(is_sig)
        sig_idx <- sig_idx[order(padj_vals[sig_idx])]
        top_n   <- min(15, length(sig_idx))
        top_idx <- sig_idx[seq_len(top_n)]

        # Build display table
        display_df <- data.frame(
            Protein      = resolve_names(sdf, top_idx),
            `log2(FC)`   = round(log2fc[top_idx], 2),
            `adj.P.Val`  = signif(padj_vals[top_idx], 3),
            `Linear FC`  = round(lfc_vals[top_idx], 2),
            check.names  = FALSE,
            stringsAsFactors = FALSE
        )

        tbl_context <- sprintf(
            "Top %d DE proteins for %s. ",
            top_n, gsub("_vs_", " vs ", ctr))
        tbl_context <- paste0(tbl_context,
            sprintf("Top: %s (log2FC=%.2f, adj.p=%.2e). Total DE: %d.",
                    display_df$Protein[1], display_df$`log2(FC)`[1],
                    display_df$adj.P.Val[1], n_sig))
        tbl_bl <- generate_slide_bottom_line(tbl_context, config = config, mode = "proteomics")

        pptx <- pptx_add_table_slide(pptx,
            title      = paste("Top DE Proteins:",
                               gsub("_vs_", " vs ", ctr)),
            display_df = display_df,
            subtitle   = sprintf(
                "%d significant proteins (adj.p < %s, linear FC >= %s)",
                n_sig, p_cut, fc_lin_cut),
            bl         = tbl_bl,
            footer     = footer_txt)
    }

    # ==================================================================
    # SLIDE: DE Heatmaps (per contrast)
    # ==================================================================
    if (dir.exists(qc_post)) {
        heatmap_pngs <- list.files(qc_post,
                                    pattern = "^heatmap_top.*\\.png$",
                                    full.names = TRUE)
        for (hf in heatmap_pngs) {
            # Extract contrast name: heatmap_top25_ContrastName.png
            cn <- gsub("^heatmap_top[0-9]+_|\\.png$", "", basename(hf))

            hm_context <- sprintf(
                "Heatmap of top DE proteins for %s with hierarchical clustering.",
                gsub("_vs_", " vs ", cn))
            hm_bl <- generate_slide_bottom_line(hm_context,
                                                 image_path = hf,
                                                 config = config, mode = "proteomics")
            pptx <- pptx_add_figure_slide(pptx,
                title    = paste("DE Heatmap:",
                                 gsub("_vs_", " vs ", cn)),
                img      = hf,
                subtitle = "Imputed log2 intensities | hierarchical clustering",
                bl       = hm_bl,
                footer   = footer_txt)
        }
    }

    # ==================================================================
    # SLIDE: PPI Network (if available)
    # ==================================================================
    if (!is.null(ppi_res)) {
        ppi_plots_dir <- file.path(out_dir, "ppi_networks", "plots")

        # Main PPI network plot
        ppi_net_png <- file.path(ppi_plots_dir, "ppi_network.png")
        if (file.exists(ppi_net_png)) {
            ppi_summary <- ppi_res$summary %||% list()
            ppi_context <- paste0(
                "Protein-protein interaction network from STRING. ",
                sprintf("%d nodes, %d edges.",
                        ppi_summary$n_nodes %||% 0,
                        ppi_summary$n_edges %||% 0))
            if (!is.null(ppi_summary$top_hubs) &&
                length(ppi_summary$top_hubs) > 0) {
                ppi_context <- paste0(ppi_context,
                    sprintf(" Top hub proteins: %s.",
                            paste(head(ppi_summary$top_hubs, 5),
                                  collapse = ", ")))
            }
            ppi_bl <- generate_slide_bottom_line(ppi_context,
                                                  image_path = ppi_net_png,
                                                  config = config, mode = "proteomics")

            subtitle_ppi <- sprintf(
                "%d nodes | %d edges | density: %.3f",
                ppi_summary$n_nodes %||% 0,
                ppi_summary$n_edges %||% 0,
                ppi_summary$density %||% 0)

            pptx <- pptx_add_figure_slide(pptx,
                title    = "PPI Network",
                img      = ppi_net_png,
                subtitle = subtitle_ppi,
                bl       = ppi_bl,
                footer   = footer_txt)
        }

        # Community detection plot
        ppi_comm_png <- file.path(ppi_plots_dir, "network_communities.png")
        if (file.exists(ppi_comm_png)) {
            comm_context <- sprintf(
                "PPI network community structure. %d communities, modularity: %.3f.",
                ppi_res$summary$n_communities %||% 0,
                ppi_res$summary$modularity %||% 0)
            comm_bl <- generate_slide_bottom_line(comm_context,
                                                   image_path = ppi_comm_png,
                                                   config = config, mode = "proteomics")
            pptx <- pptx_add_figure_slide(pptx,
                title    = "PPI Network Communities",
                img      = ppi_comm_png,
                subtitle = sprintf("%d communities | modularity: %.3f",
                                   ppi_res$summary$n_communities %||% 0,
                                   ppi_res$summary$modularity %||% 0),
                bl       = comm_bl,
                footer   = footer_txt)
        }

        # Hub proteins plot
        ppi_hub_png <- file.path(ppi_plots_dir, "hub_proteins.png")
        if (file.exists(ppi_hub_png)) {
            hub_context <- "Top hub proteins by degree in PPI network."
            if (!is.null(ppi_res$summary$top_hubs)) {
                hub_context <- paste0(hub_context,
                    sprintf(" Top: %s.",
                            paste(head(ppi_res$summary$top_hubs, 5),
                                  collapse = ", ")))
            }
            hub_bl <- generate_slide_bottom_line(hub_context,
                                                  image_path = ppi_hub_png,
                                                  config = config, mode = "proteomics")
            pptx <- pptx_add_figure_slide(pptx,
                title    = "PPI Hub Proteins",
                img      = ppi_hub_png,
                subtitle = sprintf("%d hub proteins (top 10%% by degree)",
                                   ppi_res$summary$n_hub_proteins %||% 0),
                bl       = hub_bl,
                footer   = footer_txt)
        }

        # Network topology plot
        ppi_topo_png <- file.path(ppi_plots_dir, "network_topology.png")
        if (file.exists(ppi_topo_png)) {
            topo_context <- "PPI network topology: degree distribution and centrality metrics."
            topo_bl <- generate_slide_bottom_line(topo_context,
                                                   image_path = ppi_topo_png,
                                                   config = config, mode = "proteomics")
            pptx <- pptx_add_figure_slide(pptx,
                title    = "PPI Network Topology",
                img      = ppi_topo_png,
                subtitle = "Degree distribution and betweenness centrality",
                bl       = topo_bl,
                footer   = footer_txt)
        }

        # Subcellular localization plot
        ppi_subcell_png <- file.path(ppi_plots_dir,
                                      "subcellular_enrichment.png")
        if (file.exists(ppi_subcell_png)) {
            subcell_context <- paste0(
                "Subcellular localization of significant proteins ",
                "in the PPI network.")
            subcell_bl <- generate_slide_bottom_line(subcell_context,
                image_path = ppi_subcell_png, config = config, mode = "proteomics")
            pptx <- pptx_add_figure_slide(pptx,
                title    = "Subcellular Localization",
                img      = ppi_subcell_png,
                subtitle = "GO Cellular Component annotations",
                bl       = subcell_bl,
                footer   = footer_txt)
        }
    }

    # ==================================================================
    # SLIDE: Enrichment Plots
    # ==================================================================
    enrich_plot_dir <- file.path(out_dir, "Enrichment", "plots")
    if (dir.exists(enrich_plot_dir)) {
        enrich_pngs <- list.files(enrich_plot_dir, pattern = "\\.png$",
                                   full.names = TRUE, recursive = TRUE)
        for (ep in enrich_pngs) {
            ep_label <- gsub("_", " ",
                             tools::file_path_sans_ext(basename(ep)))
            ep_label <- tools::toTitleCase(ep_label)

            ep_context <- sprintf("Pathway enrichment: %s.", ep_label)
            ep_bl <- generate_slide_bottom_line(ep_context,
                                                 image_path = ep,
                                                 config = config, mode = "proteomics")
            pptx <- pptx_add_figure_slide(pptx,
                title    = paste("Enrichment:", ep_label),
                img      = ep,
                bl       = ep_bl,
                footer   = footer_txt)
        }
    }

    # ------------------------------------------------------------------
    # Save to Results root (parent of mode dir)
    # ------------------------------------------------------------------
    results_root <- dirname(out_dir)
    out_file     <- file.path(results_root, "proteomics_summary.pptx")
    print(pptx, target = out_file)
    message("PowerPoint saved: ", out_file)
    out_file
}
