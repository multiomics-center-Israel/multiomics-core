#' Proteomics AI Commentary
#'
#' Extends R/core/12_commentary.R with proteomics-specific figure discovery
#' and interpretation guidance. Reuses core functions:
#'   generate_all_commentary(), run_claude_code_commentary(),
#'   format_commentary_html(), build_figure_context()

# ==============================================================================
# FIGURE DISCOVERY
# ==============================================================================

#' Discover proteomics figures for AI commentary
#'
#' Scans proteomics output directories for PNG files and builds
#' a metadata table compatible with generate_all_commentary().
#'
#' @param out_dir  Proteomics output directory (e.g. .../proteomics)
#' @param config   Full pipeline config
#' @return Data frame with figure_id, filepath, plot_type, section, title, etc.
discover_proteomics_figures <- function(out_dir, config = NULL) {

    diag_dir    <- file.path(out_dir, "Diagnostic_plots")
    qc_post     <- file.path(diag_dir, "QC_post")
    enrich_plt  <- file.path(out_dir, "Enrichment", "plots")
    ppi_plt     <- file.path(out_dir, "ppi_networks", "plots")
    adv_plt     <- file.path(out_dir, "advanced_stats", "plots")

    figures <- list()

    add_if_exists <- function(id, filename, dir, plot_type, section, title, desc,
                               x_axis = NA, y_axis = NA) {
        fp <- file.path(dir, filename)
        if (file.exists(fp)) {
            figures[[id]] <<- data.frame(
                figure_id = id, filepath = fp, plot_type = plot_type,
                section = section, title = title, description = desc,
                x_axis = x_axis, y_axis = y_axis, contrast = NA_character_,
                stringsAsFactors = FALSE
            )
        }
    }

    # --- QC Pre-DE ---
    add_if_exists("pca_plot", "PCA_PC1.vs.PC2.png", diag_dir,
                  "scatter", "QC", "PCA Plot — PC1 vs PC2",
                  "Principal component analysis of protein expression — first two components",
                  "PC1", "PC2")

    add_if_exists("pca_plot_13", "PCA_PC1.vs.PC3.png", diag_dir,
                  "scatter", "QC", "PCA Plot — PC1 vs PC3",
                  "PCA — first and third components", "PC1", "PC3")

    add_if_exists("density_plot", "protein_histograms_summary.png", diag_dir,
                  "density", "QC", "Protein Expression Density",
                  "Per-sample protein expression distribution", "Expression", "Density")

    add_if_exists("umap_plot", "UMAP.png", diag_dir,
                  "scatter", "QC", "UMAP — Protein Expression",
                  "UMAP dimensionality reduction of samples", "UMAP1", "UMAP2")

    # Imputation figures
    imp_w <- NULL; imp_s <- NULL
    if (!is.null(config)) {
        imp_w <- config$modes$proteomics$imputation$width
        imp_s <- config$modes$proteomics$imputation$downshift
    }
    if (!is.null(imp_w) && !is.null(imp_s)) {
        add_if_exists("imputation_histogram",
                      sprintf("imputed_histograms_samples_summary_w%s_s%s.png", imp_w, imp_s),
                      diag_dir, "histogram", "QC", "Imputation QC — Histograms",
                      "Original vs imputed value distributions per sample",
                      "Expression", "Count")
    }
    add_if_exists("imputation_boxplot", "imputed_boxplot.png", diag_dir,
                  "boxplot", "QC", "Post-imputation Boxplots",
                  "Per-sample expression boxplots after imputation", "Sample", "Expression")
    add_if_exists("pre_imputation_boxplot", "pre_imputation_boxplot.png", diag_dir,
                  "boxplot", "QC", "Pre-imputation Boxplots",
                  "Per-sample expression boxplots before imputation", "Sample", "Expression")

    add_if_exists("sample_correlation", "sample_correlation_heatmap.png", diag_dir,
                  "heatmap", "QC", "Sample-to-Sample Correlation",
                  "Pearson correlation between samples", "Samples", "Samples")

    add_if_exists("sample_distance", "sample_distance_heatmap.png", diag_dir,
                  "heatmap", "QC", "Sample-to-Sample Distance",
                  "Euclidean distance between samples", "Samples", "Samples")

    add_if_exists("expr_heatmap", "samples_protein_heatmap.png", diag_dir,
                  "heatmap", "QC", "Protein Expression Heatmap",
                  "Top variable proteins across samples", "Samples", "Proteins")

    # --- Volcano, MA, heatmap plots per contrast ---
    if (dir.exists(qc_post)) {
        for (vf in list.files(qc_post, pattern = "^volcano_.*\\.png$", full.names = TRUE)) {
            cn <- gsub("^volcano_|\\.png$", "", basename(vf))
            fid <- paste0("volcano_", cn)
            figures[[fid]] <- data.frame(
                figure_id = fid, filepath = vf, plot_type = "volcano",
                section = "Differential Expression",
                title = paste("Volcano:", gsub("_", " ", cn)),
                description = "FC vs significance for all proteins",
                x_axis = "log2 Fold Change", y_axis = "-log10(adj p-value)",
                contrast = cn, stringsAsFactors = FALSE
            )
        }

        for (hf in list.files(qc_post, pattern = "^heatmap_top.*\\.png$", full.names = TRUE)) {
            cn <- gsub("^heatmap_top[0-9]+_|\\.png$", "", basename(hf))
            fid <- paste0("heatmap_", gsub("\\.png$", "", basename(hf)))
            figures[[fid]] <- data.frame(
                figure_id = fid, filepath = hf, plot_type = "heatmap",
                section = "Differential Expression",
                title = paste("Top DE Heatmap:", gsub("_", " ", cn)),
                description = "Expression of top DE proteins",
                x_axis = "Samples", y_axis = "Proteins",
                contrast = cn, stringsAsFactors = FALSE
            )
        }
    }

    # --- Pathway plots ---
    if (dir.exists(enrich_plt)) {
        for (pf in list.files(enrich_plt, pattern = "\\.png$", full.names = TRUE)) {
            label <- gsub("^pathway_|\\.png$", "", basename(pf))
            fid <- paste0("pathway_", label)
            figures[[fid]] <- data.frame(
                figure_id = fid, filepath = pf, plot_type = "dotplot",
                section = "Pathway Analysis",
                title = paste("Pathway:", gsub("_", " ", label)),
                description = "Enriched pathways from gene set analysis",
                x_axis = "Enrichment Score", y_axis = "Pathway",
                contrast = label, stringsAsFactors = FALSE
            )
        }
    }

    # --- PPI plots ---
    if (dir.exists(ppi_plt)) {
        add_if_exists("ppi_network", "ppi_network.png", ppi_plt,
                      "network", "PPI Networks", "PPI Network",
                      "STRING protein-protein interaction network", NA, NA)
        add_if_exists("ppi_topology", "network_topology.png", ppi_plt,
                      "bar", "PPI Networks", "Network Topology",
                      "Degree distribution and centrality metrics", "Degree", "Count")
        add_if_exists("ppi_communities", "network_communities.png", ppi_plt,
                      "network", "PPI Networks", "Network Communities",
                      "Community structure of PPI network", NA, NA)
        add_if_exists("ppi_hubs", "hub_proteins.png", ppi_plt,
                      "bar", "PPI Networks", "Hub Proteins",
                      "Top hub proteins by degree", "Protein", "Degree")
        add_if_exists("ppi_subcellular", "subcellular_enrichment.png", ppi_plt,
                      "bar", "PPI Networks", "Subcellular Localization",
                      "Subcellular distribution of significant proteins",
                      "Compartment", "Count")
    }

    # --- Advanced stats plots ---
    if (dir.exists(adv_plt)) {
        add_if_exists("adv_stats_effect_size_forest_plot", "effect_size_forest_plot.png", adv_plt,
                      "forest", "Advanced Statistics", "Effect Size Forest Plot",
                      "Effect sizes with bootstrap confidence intervals",
                      "Effect Size", "Protein")
        add_if_exists("adv_stats_power_analysis", "power_analysis.png", adv_plt,
                      "line", "Advanced Statistics", "Power Analysis",
                      "Statistical power by sample size", "Sample Size", "Power")
        add_if_exists("adv_stats_mtc_comparison", "mtc_comparison.png", adv_plt,
                      "bar", "Advanced Statistics", "MTC Comparison",
                      "Multiple testing correction method comparison",
                      "Method", "Significant Proteins")
        add_if_exists("adv_stats_icc_distribution", "icc_distribution.png", adv_plt,
                      "histogram", "Advanced Statistics", "ICC Distribution",
                      "Intraclass correlation from mixed effects models",
                      "ICC", "Count")
    }

    if (length(figures) == 0) {
        message("No proteomics figures discovered for commentary")
        return(data.frame(
            figure_id = character(), filepath = character(),
            plot_type = character(), section = character(),
            title = character(), description = character(),
            x_axis = character(), y_axis = character(),
            contrast = character(), stringsAsFactors = FALSE
        ))
    }

    figures_tbl <- do.call(rbind, figures)
    rownames(figures_tbl) <- NULL
    message("Discovered ", nrow(figures_tbl), " proteomics figures for commentary")
    figures_tbl
}

# ==============================================================================
# PROTEOMICS-SPECIFIC GUIDANCE
# ==============================================================================

#' Build proteomics-specific figure guidance
#' @return Character string
build_proteomics_figure_guidance <- function() {
    paste0(
        "PROTEOMICS-SPECIFIC INTERPRETATION NOTES:\n",
        "- Proteomics data typically has more missing values than RNA-seq; ",
        "imputation quality should be assessed.\n",
        "- Protein abundance is NOT directly comparable to RNA expression; ",
        "post-translational regulation means mRNA and protein levels often diverge.\n",
        "- Linear fold changes (not log2) are the native scale for proteomics intensity ratios.\n",
        "- Multiple imputation reduces false positives from imputed-vs-real comparisons.\n",
        "- PCA separation in proteomics may be less pronounced due to narrower dynamic range.\n",
        "- In PPI networks, hub proteins with high connectivity may be promiscuous ",
        "interactors rather than biologically central.\n"
    )
}

# ==============================================================================
# EXTRA DATA FOR CONTEXT
# ==============================================================================

#' Build proteomics extra data for commentary context
#'
#' @param de_res     DE results
#' @param qc_pre_obj QC pre-DE results
#' @param config     Full config
#' @return Named list with groups, n_samples, methods, PCA variance, top DE proteins
build_proteomics_extra_data <- function(de_res, qc_pre_obj, config) {
    cfg <- config$modes$proteomics
    extra <- list()

    # Domain label for prompt customisation
    extra$domain_label <- "proteomics"
    extra$domain_guidance <- build_proteomics_figure_guidance()

    # --- Experimental groups and sample count ---
    tryCatch({
        group_col <- cfg$effects$color %||% "Condition"
        pca_scores <- qc_pre_obj$objects$pca_scores
        if (!is.null(pca_scores) && group_col %in% names(pca_scores)) {
            extra$groups <- as.character(unique(pca_scores[[group_col]]))
            extra$n_samples <- nrow(pca_scores)
        }
    }, error = function(e) {
        message("  Could not extract group info: ", e$message)
    })

    # --- PCA variance explained ---
    tryCatch({
        var_expl <- qc_pre_obj$objects$var_expl
        if (!is.null(var_expl) && length(var_expl) >= 2) {
            extra$pca_variance <- as.numeric(var_expl[1:min(3, length(var_expl))])
            message("  PCA variance: PC1=", round(extra$pca_variance[1], 1),
                    "%, PC2=", round(extra$pca_variance[2], 1), "%")
        }
    }, error = function(e) {
        message("  Could not extract PCA variance: ", e$message)
    })

    # --- Methods ---
    extra$norm_method <- cfg$normalization$method %||% "none"
    extra$de_method <- cfg$de$method %||% "limma"
    extra$de_cutoffs <- list(
        p_cutoff  = cfg$de$p_cutoff %||% 0.05,
        fc_cutoff = cfg$de$linear_fc_cutoff %||% 1.5
    )

    # --- Top DE proteins per contrast (up to 5 up + 5 down) ---
    extra$top_genes <- list()
    tryCatch({
        sdf <- de_res$summary_df
        if (!is.null(sdf)) {
            # Detect gene/protein name column
            gene_col <- NULL
            for (candidate in c("Genes", "Protein.Names", "GeneName", "FeatureID")) {
                if (candidate %in% names(sdf)) {
                    gene_col <- candidate
                    break
                }
            }
            if (is.null(gene_col)) gene_col <- names(sdf)[1]

            padj_cols <- grep("^padj\\.imputs\\.", names(sdf), value = TRUE)
            p_cut <- extra$de_cutoffs$p_cutoff %||% 0.05

            for (pc in padj_cols) {
                cn <- sub("^padj\\.imputs\\.", "", pc)
                fc_col <- paste0("linearFC.imputs.", cn)
                if (!(fc_col %in% names(sdf))) next

                padj_vals <- as.numeric(sdf[[pc]])
                lfc_vals  <- as.numeric(sdf[[fc_col]])
                log2fc    <- sign(lfc_vals) * log2(pmax(abs(lfc_vals), 1e-10))
                genes     <- as.character(sdf[[gene_col]])

                sig_mask <- !is.na(padj_vals) & padj_vals <= p_cut
                if (sum(sig_mask) == 0) next

                sig_df <- data.frame(
                    gene   = genes[sig_mask],
                    log2FC = log2fc[sig_mask],
                    padj   = padj_vals[sig_mask],
                    stringsAsFactors = FALSE
                )

                up_df <- sig_df[sig_df$log2FC > 0, ]
                dn_df <- sig_df[sig_df$log2FC < 0, ]
                up_df <- up_df[order(-abs(up_df$log2FC)), ]
                dn_df <- dn_df[order(-abs(dn_df$log2FC)), ]
                top_combined <- rbind(
                    if (nrow(up_df) > 0) cbind(head(up_df, 5), direction = "up") else NULL,
                    if (nrow(dn_df) > 0) cbind(head(dn_df, 5), direction = "down") else NULL
                )

                if (!is.null(top_combined) && nrow(top_combined) > 0) {
                    extra$top_genes[[cn]] <- top_combined
                    message("  Top proteins for ", cn, ": ",
                            min(nrow(up_df), 5), " up + ", min(nrow(dn_df), 5), " down")
                }
            }
        }
    }, error = function(e) {
        message("  Could not extract top DE proteins: ", e$message)
    })

    extra
}

# ==============================================================================
# MAIN COMMENTARY ORCHESTRATOR
# ==============================================================================

#' Run proteomics commentary generation
#'
#' @param de_res     DE results
#' @param qc_pre_obj QC pre-DE results
#' @param config     Full config
#' @param out_dir    Proteomics output directory
#' @return Path to commentary_all.json (or NULL)
run_proteomics_commentary <- function(de_res, qc_pre_obj, config, out_dir) {
    cfg <- config$modes$proteomics
    comm_cfg <- cfg$commentary %||% list()

    if (!isTRUE(comm_cfg$enabled)) {
        message("Proteomics commentary disabled.")
        return(NULL)
    }

    # Discover figures
    figures_tbl <- discover_proteomics_figures(out_dir, config)
    if (nrow(figures_tbl) == 0) return(NULL)

    # Build extra data
    extra_data <- build_proteomics_extra_data(de_res, qc_pre_obj, config)

    # Build DE summary counts for context
    de_summary <- tryCatch(
        build_de_summary_counts_proteomics(de_res$summary_df),
        error = function(e) { message("  Could not build DE summary: ", e$message); NULL }
    )

    # Override config to point at proteomics commentary settings

    # generate_all_commentary reads config$modes$rna$commentary,
    # so we temporarily inject our proteomics commentary config there
    # Deep copy to avoid mutating the original config
    config_for_commentary <- config
    config_for_commentary$modes <- config$modes
    config_for_commentary$modes$rna <- list()
    config_for_commentary$modes$rna$commentary <- comm_cfg
    config_for_commentary$modes$rna$annotation <- cfg$annotation
    config_for_commentary$modes$rna$effects <- cfg$effects
    config_for_commentary$modes$rna$filtering <- cfg$filtering

    # Commentary output dir
    commentary_dir <- file.path(out_dir, "commentary")

    commentary_list <- generate_all_commentary(
        figures_tbl = figures_tbl,
        config      = config_for_commentary,
        de_summary  = de_summary,
        extra_data  = extra_data,
        output_dir  = commentary_dir
    )

    # Return path to combined JSON
    combined_path <- file.path(commentary_dir, "commentary_all.json")
    if (file.exists(combined_path)) combined_path else NULL
}
