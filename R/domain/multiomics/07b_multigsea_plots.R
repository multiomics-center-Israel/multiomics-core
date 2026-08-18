#' MultiGSEA Correlation Plots
#'
#' Functions for visualizing cross-omics enrichment correlations
#' and generating pathway overlay plots using pathview.
#'
#' @name multigsea_plots
NULL

#' Run MultiGSEA cross-omics enrichment summary
#'
#' Builds the combined MultiGSEA summary panel from the per-omics enrichment
#' tables.
#'
#' Pairwise "omics vs omics" scatter plots are deliberately not produced. They
#' were built over the *union* of the two layers' terms, so a term enriched in
#' one layer and simply absent from the other entered as -log10(FDR) = 0 and sat
#' on an axis; the reported correlation then tracked how much term coverage the
#' two GMTs happened to share rather than agreement between the layers.
#'
#' @param enrichment_results List containing enrichment results from run_multiomics_enrichment.
#' @param config Pipeline configuration list.
#' @param out_dir Output directory for plots and tables.
#' @return A named list holding the combined summary plot, or NULL when there is
#'   nothing to plot.
#' @export
run_multigsea_plots <- function(enrichment_results, config, out_dir = NULL) {
    message("=== Running MultiGSEA Cross-Omics Summary ===")

    if (is.null(enrichment_results) || is.null(enrichment_results$per_omics)) {
        message("No enrichment results available for MultiGSEA.")
        return(NULL)
    }

    mg_config <- config$modes$multiomics$enrichment$multigsea %||% list()
    if (!(mg_config$run_multigsea %||% TRUE)) {
        message("MultiGSEA analysis disabled in config.")
        return(NULL)
    }

    if (!is.null(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }

    per_omics <- enrichment_results$per_omics

    if (length(per_omics) < 2) {
        message("Need at least 2 omics with enrichment results for MultiGSEA.")
        return(NULL)
    }

    plots <- list()
    if (!is.null(out_dir)) {
        combined <- plot_multigsea_combined(per_omics, out_dir)
        if (!is.null(combined)) plots[["combined"]] <- combined
    }

    message("MultiGSEA plots generated: ", length(plots))
    plots
}


#' Combined MultiGSEA Enrichment Plot
#'
#' Dot plot of the pathways enriched in more than one omics layer (or, when
#' there are none, the strongest terms overall), one column per layer.
#'
#' @param per_omics Per-omics enrichment results list.
#' @param out_dir Output directory.
#' @return The summary ggplot object, or NULL when no layer carries usable
#'   term / p-value columns.
plot_multigsea_combined <- function(per_omics, out_dir) {
    omics_names <- names(per_omics)
    summary_rows <- list()

    for (om in omics_names) {
        res <- per_omics[[om]]
        if (is.null(res) || !is.data.frame(res)) next

        # "pathway" belongs in this list: the fgsea / ORA tables the gene layers
        # arrive in key their term on `pathway` (`pathway_name` holds the
        # readable label) and carry none of the other three, so without it both
        # gene layers were dropped here and the panel drew metabolomics alone.
        term_col <- NULL
        for (tc in c("term", "ID", "Description", "pathway")) {
            if (tc %in% colnames(res)) { term_col <- tc; break }
        }
        if (is.null(term_col)) next

        # Find p-value column
        padj_col <- NULL
        for (pc in c("padj", "p.adjust", "adj.P.Val", "FDR", "qvalue", "pvalue")) {
            if (pc %in% colnames(res)) { padj_col <- pc; break }
        }
        if (is.null(padj_col)) next

        # Resolve a readable label, falling back to the raw id
        term_ids <- res[[term_col]]
        term_names <- if ("pathway_name" %in% colnames(res)) {
            res$pathway_name
        } else if ("pathway" %in% colnames(res)) {
            res$pathway
        } else {
            term_ids
        }

        df_tmp <- data.frame(
            term = term_ids,
            pathway_name = term_names,
            padj = res[[padj_col]],
            omic = om,
            stringsAsFactors = FALSE
        )
        df_tmp <- df_tmp[!is.na(df_tmp$padj), ]
        summary_rows[[length(summary_rows) + 1]] <- df_tmp
    }

    if (length(summary_rows) == 0) {
        message("No data for summary panel.")
        return(NULL)
    }

    summary_df <- do.call(rbind, summary_rows)

    # Collapse to one row per (omic, term), keeping that layer's best padj. The
    # gene layers arrive with contrasts, databases and methods stacked in one
    # frame, so without this a term significant twice inside a single layer
    # would satisfy the "supported by >= 2 omics" rule on its own.
    summary_df <- summary_df[order(summary_df$padj), ]
    summary_df <- summary_df[!duplicated(summary_df[, c("omic", "term")]), ]

    # Find terms significant in at least 2 omics
    sig_terms <- summary_df[summary_df$padj < 0.05, ]
    term_counts <- table(sig_terms$term)
    co_sig <- names(term_counts[term_counts >= 2])

    if (length(co_sig) == 0) {
        # Fall back to top terms by lowest p-value across any omic
        top_terms <- unique(summary_df$term[order(summary_df$padj)])[1:min(15, length(unique(summary_df$term)))]
        plot_data <- summary_df[summary_df$term %in% top_terms, ]
        panel_title <- "Top Enriched Pathways"
    } else {
        plot_data <- summary_df[summary_df$term %in% co_sig, ]
        # Keep top 15 by mean -log10(padj) across omics
        mean_scores <- tapply(-log10(pmax(plot_data$padj, 1e-300)), plot_data$term, mean, na.rm = TRUE)
        top_terms <- names(sort(mean_scores, decreasing = TRUE))[1:min(15, length(mean_scores))]
        plot_data <- plot_data[plot_data$term %in% top_terms, ]
        panel_title <- "Co-Significant Pathways"
    }

    plot_data$neg_log10_padj <- -log10(pmax(plot_data$padj, 1e-300))
    # Use pathway name, truncate long names
    plot_data$term_short <- ifelse(
        !is.na(plot_data$pathway_name) & nzchar(plot_data$pathway_name),
        plot_data$pathway_name,
        plot_data$term
    )
    plot_data$term_short <- ifelse(
        nchar(plot_data$term_short) > 40,
        paste0(substr(plot_data$term_short, 1, 37), "..."),
        plot_data$term_short
    )
    plot_data$omic_label <- gsub("_", " ", tools::toTitleCase(plot_data$omic))

    combined <- ggplot2::ggplot(plot_data,
        ggplot2::aes(x = omic_label, y = stats::reorder(term_short, neg_log10_padj),
                     size = neg_log10_padj, color = neg_log10_padj)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_viridis_c(name = "-log10(padj)") +
        ggplot2::scale_size_continuous(name = "-log10(padj)", range = c(2, 8)) +
        ggplot2::labs(title = panel_title, x = NULL, y = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 7),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
            axis.title = ggplot2::element_text(face = "bold")
        )

    # Filename kept as-is: the deck builder in scripts/ and the report both
    # look for multigsea_combined_enrichment.png.
    ggplot2::ggsave(
        file.path(out_dir, "multigsea_combined_enrichment.png"),
        plot = combined, width = 10, height = 8, dpi = 300
    )
    message("Saved combined MultiGSEA plot.")

    combined
}


#' Run Pathview Visualization for Agreed Pathways
#'
#' Generates KEGG pathway overlay plots for pathways significant in multiple omics.
#'
#' @param enrichment_results List containing enrichment results.
#' @param mae_data MultiAssayExperiment object (or list with harmonized_omics).
#' @param config Pipeline configuration list.
#' @param out_dir Output directory for pathview plots.
#' @return List of generated plot paths.
#' @export
run_multigsea_pathview <- function(enrichment_results, mae_data, config, out_dir = NULL) {
    message("=== Running MultiGSEA Pathview Visualization ===")

    if (!requireNamespace("pathview", quietly = TRUE)) {
        message("Package 'pathview' not installed. Skipping.")
        return(NULL)
    }

    mg_config <- config$modes$multiomics$enrichment$multigsea %||% list()
    if (!(mg_config$run_pathview %||% TRUE)) {
        message("Pathview analysis disabled in config.")
        return(NULL)
    }

    if (is.null(enrichment_results) || is.null(enrichment_results$per_omics)) {
        message("No enrichment results available.")
        return(NULL)
    }

    # Create output directory
    if (is.null(out_dir)) {
        out_dir <- tempdir()
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Identify Agreed KEGG Pathways
    per_omics <- enrichment_results$per_omics
    omics_names <- names(per_omics)

    kegg_pathways <- list()

    for (omic in omics_names) {
        res <- per_omics[[omic]]
        if (!is.null(res) && nrow(res) > 0) {
            # Check for KEGG results
            if ("type" %in% colnames(res)) {
                kegg_res <- res[res$type == "KEGG", ]
            } else {
                # Heuristic: IDs start with "hsa" or "map" or numeric
                kegg_res <- res[grep("^hsa|^map|^[0-9]+$", res$ID), ]
            }

            if (nrow(kegg_res) > 0) {
                kegg_pathways[[omic]] <- kegg_res$ID
            }
        }
    }

    if (length(kegg_pathways) < 2) {
        message("Less than 2 omics have KEGG results. Skipping Pathview.")
        return(NULL)
    }

    # Find common pathways
    all_kegg <- unlist(kegg_pathways)
    if (length(all_kegg) == 0) {
        message("No KEGG pathways found.")
        return(NULL)
    }

    pathway_counts <- table(all_kegg)
    common_pathways <- names(pathway_counts)[pathway_counts >= 2]

    if (length(common_pathways) == 0) {
        message("No agreed KEGG pathways found between omics.")
        return(NULL)
    }

    message("Found ", length(common_pathways), " agreed KEGG pathways.")

    # Prepare Data for Pathview
    # Helper to get fold changes
    get_logfc <- function(omic_name, id_col = "entrez_id") {
        if (!omic_name %in% names(mae_data$harmonized_omics)) {
            return(NULL)
        }

        dat <- mae_data$harmonized_omics[[omic_name]]
        de <- dat$de_table %||% dat$da_table
        anno <- dat$feature_annotation

        if (is.null(de) || is.null(anno)) {
            return(NULL)
        }

        merged <- merge(de, anno, by = "feature_id")

        if (!id_col %in% colnames(merged)) {
            return(NULL)
        }

        # Get LogFC column
        fc_col <- grep("logFC|log2FoldChange", colnames(merged), ignore.case = TRUE, value = TRUE)[1]
        if (is.na(fc_col)) {
            return(NULL)
        }

        # Create named vector (take mean for duplicates)
        vec <- tapply(merged[[fc_col]], merged[[id_col]], mean, na.rm = TRUE)
        return(vec)
    }

    # Gene data (Transcriptomics + Proteomics)
    gene_data <- NULL
    rna_fc <- get_logfc("transcriptomics", "entrez_id")
    prot_fc <- get_logfc("proteomics", "entrez_id")

    if (!is.null(rna_fc) && !is.null(prot_fc)) {
        all_genes <- unique(c(names(rna_fc), names(prot_fc)))
        gene_data <- matrix(NA,
            nrow = length(all_genes), ncol = 2,
            dimnames = list(all_genes, c("Transcriptomics", "Proteomics"))
        )
        idx_rna <- match(names(rna_fc), all_genes)
        gene_data[idx_rna, 1] <- rna_fc
        idx_prot <- match(names(prot_fc), all_genes)
        gene_data[idx_prot, 2] <- prot_fc
    } else if (!is.null(rna_fc)) {
        gene_data <- rna_fc
    } else if (!is.null(prot_fc)) {
        gene_data <- prot_fc
    }

    # Metabolomics Data (CPD Data)
    cpd_data <- NULL
    met_fc <- get_logfc("metabolomics", "kegg_id")
    if (!is.null(met_fc)) {
        cpd_data <- met_fc
    }

    if (is.null(gene_data) && is.null(cpd_data)) {
        message("No valid Entrez/KEGG IDs found in data for Pathview.")
        return(NULL)
    }

    # Run Pathview
    generated_plots <- list()

    cwd <- getwd()
    setwd(out_dir)
    on.exit(setwd(cwd))

    # Determine organism code
    organism <- config$global$organism %||% "human"
    if (organism == "c_elegans") {
        species_code <- "cel"
    } else {
        species_code <- mg_config$organism_code %||% "hsa"
    }

    for (pid in common_pathways) {
        # Clean ID
        clean_pid <- sub("^[a-z]+", "", pid)

        tryCatch({
            pv_out <- pathview::pathview(
                gene.data = gene_data,
                cpd.data = cpd_data,
                pathway.id = clean_pid,
                species = species_code,
                out.suffix = "multiomics",
                temp.file = TRUE,
                kegg.dir = out_dir,
                keys.align = "y",
                match.data = TRUE,
                multi.state = !is.null(dim(gene_data)) && ncol(gene_data) > 1,
                same.layer = FALSE
            )

            outfile <- paste0(species_code, clean_pid, ".multiomics.png")
            if (file.exists(outfile)) {
                generated_plots[[pid]] <- file.path(out_dir, outfile)
                message("Generated Pathview: ", outfile)
            }
        }, error = function(e) {
            message("Pathview failed for ", pid, ": ", e$message)
        })
    }

    return(generated_plots)
}


# =============================================================================
# Multi-ORA: Combined Over-Representation Analysis across omics
# =============================================================================

#' Run multi-omics ORA (multi-ORA)
#'
#' Pools significant features from all omics layers, maps them to KEGG
#' pathways, and runs a single combined hypergeometric test per pathway.
#' Gene-based omics (RNA, proteomics) are pooled via ENTREZID; metabolomics
#' uses compound-level enrichment separately. Results are combined using
#' Fisher's method.
#'
#' @param de_results Named list of DE results per omics
#' @param harmonization_res Harmonization result with MAE and pre-processing data
#' @param config Full config object
#' @param out_dir Output directory for results and plots
#' @return List with: results (data.frame), plots (list of paths)
run_multi_ora <- function(de_results, harmonization_res, config, out_dir) {
    # KEGG classes this organism cannot have; shared with the cross-omics tables
    # and the pathway maps so every figure filters the same way.
    excl_classes <- config$modes$multiomics$enrichment$exclude_pathway_classes

    message("=== Running Multi-ORA (combined cross-omics ORA) ===")

    if (is.null(de_results) || length(de_results) < 2) {
        message("Multi-ORA requires DE results from at least 2 omics layers")
        return(NULL)
    }

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Track errors/warnings for reporting in HTML
    .ora_issues <- character(0)

    organism <- config$global$organism
    kegg_org <- get_kegg_organism(organism)
    org_db <- get_organism_db(organism)

    if (is.null(kegg_org) || is.null(org_db)) {
        # No KEGG organism / OrgDb (e.g. non-model organisms): fall back to
        # GMT-based ORA when the omics blocks supply custom gene sets.
        gmt_res <- run_multi_ora_gmt(de_results, harmonization_res, config, out_dir)
        if (is.null(gmt_res)) {
            message("Multi-ORA: organism annotation not available for ", organism,
                    " and no usable per-omic GMT gene sets")
        }
        # Without an OrgDb, pathway maps are still reachable: organisms with a
        # KEGG code render natively, and organisms with none render the KEGG
        # reference maps in KO space when a ko_map is configured.
        tryCatch(
            generate_per_omic_union_pathview(de_results, harmonization_res, config, out_dir),
            error = function(e) message("  Union pathview failed: ", conditionMessage(e))
        )
        return(gmt_res)
    }

    # --- Collect per-omics significant KEGG gene IDs ---
    gene_omics <- c("transcriptomics", "proteomics")
    per_omics_sig <- list()
    per_omics_universe <- list()

    for (om in intersect(gene_omics, names(de_results))) {
        de_data <- de_results[[om]]
        de_tables <- extract_de_tables(de_data, om, harmonization_res)
        if (is.null(de_tables) || length(de_tables) == 0) next

        # Map to ENTREZ IDs
        id_map <- tryCatch(
            map_feature_ids_to_entrez(de_tables, om, harmonization_res, org_db),
            error = function(e) {
                msg <- paste0(om, " ENTREZ ID mapping failed: ", e$message)
                message("  Multi-ORA: ", msg)
                .ora_issues <<- c(.ora_issues, msg)
                NULL
            }
        )
        if (is.null(id_map) || nrow(id_map) == 0) {
            msg <- paste0(om, " produced no ENTREZ mappings")
            message("  Multi-ORA: ", msg, ", skipping")
            .ora_issues <<- c(.ora_issues, msg)
            next
        }

        # Collect significant features (padj < 0.05 in any contrast)
        sig_ids <- character(0)
        all_ids <- character(0)
        for (nm in names(de_tables)) {
            df <- de_tables[[nm]]
            df_mapped <- merge(df, id_map, by = "feature_id")
            all_ids <- c(all_ids, df_mapped$ENTREZID)

            sig <- df_mapped$ENTREZID[!is.na(df_mapped$padj) & df_mapped$padj < 0.05]
            if (length(sig) < 5) {
                sig <- df_mapped$ENTREZID[!is.na(df_mapped$pvalue) & df_mapped$pvalue < 0.05]
            }
            sig_ids <- c(sig_ids, sig)
        }

        sig_ids <- unique(sig_ids[!is.na(sig_ids)])
        all_ids <- unique(all_ids[!is.na(all_ids)])

        if (length(sig_ids) > 0) {
            per_omics_sig[[om]] <- sig_ids
            per_omics_universe[[om]] <- all_ids
            message("  ", om, ": ", length(sig_ids), " sig / ", length(all_ids), " total ENTREZ IDs")
        }
    }

    if (length(per_omics_sig) == 0) {
        message("Multi-ORA: no gene-based omics produced significant features")
        .ora_issues <- c(.ora_issues, "No gene-based omics produced significant ENTREZ features")
        writeLines(.ora_issues, file.path(out_dir, "multi_ora_issues.txt"))
        return(NULL)
    }

    # Pool significant genes across gene-based omics
    pooled_sig <- unique(unlist(per_omics_sig))
    pooled_universe <- unique(unlist(per_omics_universe))
    message("  Pooled: ", length(pooled_sig), " sig / ", length(pooled_universe),
            " universe genes across ", length(per_omics_sig), " omics")

    # Convert ENTREZID to KEGG gene IDs
    kegg_conv <- tryCatch(
        convert_entrez_to_kegg(pooled_universe, kegg_org),
        error = function(e) {
            msg <- paste0("KEGG API error during ENTREZ-to-KEGG conversion: ", e$message)
            message("  Multi-ORA: ", msg)
            .ora_issues <<- c(.ora_issues, msg)
            NULL
        }
    )
    if (is.null(kegg_conv)) {
        message("Multi-ORA: KEGG ID conversion failed (possible KEGG API timeout)")
        .ora_issues <- c(.ora_issues, "KEGG ID conversion failed — KEGG REST API may be unavailable")
        # Save issues and return
        writeLines(.ora_issues, file.path(out_dir, "multi_ora_issues.txt"))
        return(NULL)
    }

    pooled_sig_kegg <- unique(kegg_conv[pooled_sig])
    pooled_sig_kegg <- pooled_sig_kegg[!is.na(pooled_sig_kegg)]
    pooled_univ_kegg <- unique(kegg_conv[pooled_universe])
    pooled_univ_kegg <- pooled_univ_kegg[!is.na(pooled_univ_kegg)]

    # Also convert per-omics sig lists to KEGG IDs (for per-omics ORA)
    per_omics_sig_kegg <- lapply(per_omics_sig, function(ids) {
        k <- kegg_conv[ids]
        unique(k[!is.na(k)])
    })

    # --- Run pooled ORA ---
    message("  Running pooled gene ORA...")
    pooled_ora <- run_multi_ora_kegg(
        sig_genes = pooled_sig_kegg,
        universe = pooled_univ_kegg,
        kegg_org = kegg_org,
        label = "pooled"
    )

    # --- Run per-omics ORA (with same universe) ---
    per_omics_ora <- list()
    for (om in names(per_omics_sig_kegg)) {
        message("  Running per-omics ORA for ", om, "...")
        per_omics_ora[[om]] <- run_multi_ora_kegg(
            sig_genes = per_omics_sig_kegg[[om]],
            universe = pooled_univ_kegg,
            kegg_org = kegg_org,
            label = om
        )
    }

    # --- Metabolomics compound ORA (separate) ---
    metab_ora <- NULL
    if ("metabolomics" %in% names(de_results)) {
        metab_ora <- tryCatch({
            de_tables <- extract_de_tables(de_results$metabolomics, "metabolomics", harmonization_res)
            id_map <- map_metabolite_ids_to_kegg(de_tables, harmonization_res)
            if (!is.null(id_map) && nrow(id_map) > 0) {
                full_universe <- unique(id_map$KEGG_CPD[!is.na(id_map$KEGG_CPD)])
                # Merge DE stats with KEGG IDs
                de_df <- do.call(rbind, de_tables)
                de_mapped <- merge(de_df, id_map, by = "feature_id")
                de_mapped$KEGG_ID <- de_mapped$KEGG_CPD
                run_compound_ora(de_mapped, out_dir, 2, 500, 0.1, universe = full_universe)
            } else NULL
        }, error = function(e) {
            message("  Metabolomics compound ORA failed: ", e$message)
            NULL
        })
        if (!is.null(metab_ora) && nrow(metab_ora) > 0) {
            message("  metabolomics: ", nrow(metab_ora), " enriched compound pathways")
            # Save compound ORA results for dedicated report section
            write.csv(metab_ora, file.path(out_dir, "compound_ora_results.csv"),
                      row.names = FALSE)
            # Generate compound ORA barplot
            compound_bp_path <- file.path(out_dir, "compound_ora_barplot.png")
            png(compound_bp_path, width = 1000, height = 700, res = 120)
            tryCatch({
                plot_multi_ora_barplot(metab_ora, "Metabolomics Compound ORA (KEGG Pathways)",
                                       exclude_classes = excl_classes, cache_dir = out_dir)
            }, error = function(e) {
                plot.new()
                text(0.5, 0.5, paste("Plot failed:", e$message), cex = 1.2)
            })
            dev.off()
            message("  Saved compound ORA barplot and results table")
        }
    }

    # --- Combine results into summary table ---
    combined <- build_multi_ora_summary(pooled_ora, per_omics_ora, metab_ora)

    if (is.null(combined) || nrow(combined) == 0) {
        message("Multi-ORA: no enriched pathways found")
        .ora_issues <- c(.ora_issues, "No enriched pathways found after ORA")
        writeLines(.ora_issues, file.path(out_dir, "multi_ora_issues.txt"))
        return(NULL)
    }

    # Save any accumulated issues (even on success, partial issues are useful)
    if (length(.ora_issues) > 0) {
        writeLines(.ora_issues, file.path(out_dir, "multi_ora_issues.txt"))
    }

    # Write results
    write.csv(combined, file.path(out_dir, "multi_ora_results.csv"), row.names = FALSE)
    message("  Multi-ORA found ", nrow(combined), " enriched pathways (pooled)")

    # --- Generate plots ---
    plots <- list()

    # 1. Pooled ORA barplot
    if (!is.null(pooled_ora) && nrow(pooled_ora) > 0) {
        plots$pooled_barplot <- file.path(out_dir, "multi_ora_pooled_barplot.png")
        png(plots$pooled_barplot, width = 1000, height = 700, res = 120)
        tryCatch({
            plot_multi_ora_barplot(pooled_ora, "Pooled Multi-ORA (All Gene-Based Omics)",
                                   exclude_classes = excl_classes, cache_dir = out_dir)
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("Plot failed:", e$message), cex = 1.2)
        })
        dev.off()
    }

    # 2. Multi-ORA dot plot comparing per-omics + pooled
    plots$dotplot <- file.path(out_dir, "multi_ora_dotplot.png")
    tryCatch({
        plot_multi_ora_dotplot(combined, per_omics_ora, metab_ora, out_dir)
    }, error = function(e) {
        message("  Multi-ORA dot plot failed: ", e$message)
    })

    # 3. Upset-style bar showing number of omics supporting each pathway
    plots$support_barplot <- file.path(out_dir, "multi_ora_support_barplot.png")
    tryCatch({
        plot_multi_ora_support(combined, out_dir)
    }, error = function(e) {
        message("  Multi-ORA support plot failed: ", e$message)
    })

    # 4. Pathview maps for pathways supported by >= 2 omics
    plots$pathview_pdf <- tryCatch({
        generate_multi_ora_pathview(
            combined = combined,
            de_results = de_results,
            harmonization_res = harmonization_res,
            config = config,
            out_dir = out_dir,
            min_support = 2
        )
    }, error = function(e) {
        message("  Multi-ORA pathview failed: ", e$message)
        NULL
    })

    # 5. Per-omics pathview: top metabolomics pathways + proteomics overlay,
    #    and top proteomics pathways + metabolomics overlay
    per_omics_pv <- tryCatch({
        generate_per_omics_pathview(
            per_omics_ora = per_omics_ora,
            metab_ora = metab_ora,
            de_results = de_results,
            harmonization_res = harmonization_res,
            config = config,
            out_dir = out_dir,
            top_n = 5
        )
    }, error = function(e) {
        message("  Per-omics pathview failed: ", e$message)
        NULL
    })
    if (!is.null(per_omics_pv)) {
        plots$pathview_metabolomics_pdf <- per_omics_pv$metabolomics_pdf
        plots$pathview_proteomics_pdf <- per_omics_pv$proteomics_pdf
    }

    # --- Per-contrast Multi-ORA ---
    # Re-extract DE tables to get per-contrast names, then run ORA per contrast
    all_de_tables <- list()
    for (om in intersect(gene_omics, names(de_results))) {
        de_tables_tmp <- extract_de_tables(de_results[[om]], om, harmonization_res)
        if (!is.null(de_tables_tmp)) all_de_tables[[om]] <- de_tables_tmp
    }
    # Also extract metabolomics DE tables for per-contrast compound ORA
    metab_de_tables <- NULL
    if ("metabolomics" %in% names(de_results)) {
        metab_de_tables <- extract_de_tables(de_results$metabolomics, "metabolomics",
                                             harmonization_res)
    }
    ora_contrast_names <- unique(unlist(lapply(all_de_tables, names)))

    if (length(ora_contrast_names) > 1) {
        message("  Generating per-contrast Multi-ORA for ",
                length(ora_contrast_names), " contrasts")
        per_contrast_dir <- file.path(out_dir, "per_contrast")

        for (cname in ora_contrast_names) {
            safe_dir <- gsub("[^a-zA-Z0-9._-]", "_", cname)
            contrast_out <- file.path(per_contrast_dir, safe_dir)
            dir.create(contrast_out, recursive = TRUE, showWarnings = FALSE)

            tryCatch({
                .run_multi_ora_contrast_group(
                    all_de_tables = all_de_tables,
                    contrast_name = cname,
                    harmonization_res = harmonization_res,
                    kegg_org = kegg_org,
                    org_db = org_db,
                    out_dir = contrast_out,
                    metab_de_tables = metab_de_tables
                )
            }, error = function(e) {
                message("    Per-contrast Multi-ORA failed for ", cname, ": ", e$message)
            })
        }
        message("  Per-contrast Multi-ORA output saved to: ", per_contrast_dir)
    }

    message("Multi-ORA complete: ", nrow(combined), " pathways")

    list(
        pooled = pooled_ora,
        per_omics = per_omics_ora,
        metabolomics = metab_ora,
        combined = combined,
        plots = plots
    )
}


#' Run KEGG ORA with Fisher's exact test for multi-ORA
#'
#' Run Multi-ORA for a single contrast
#'
#' @param all_de_tables Named list (per omics) of named lists (per contrast) of DE data frames
#' @param contrast_name Contrast name to process
#' @param harmonization_res Harmonization result
#' @param kegg_org KEGG organism code
#' @param org_db Organism annotation database
#' @param out_dir Output directory for this contrast
#' @param metab_de_tables Metabolomics DE tables (named list per contrast), or NULL
#' @return Invisible NULL
.run_multi_ora_contrast_group <- function(all_de_tables, contrast_name,
                                           harmonization_res, kegg_org, org_db,
                                           out_dir, metab_de_tables = NULL) {

    per_omics_sig <- list()
    per_omics_universe <- list()

    for (om in names(all_de_tables)) {
        de_tables <- all_de_tables[[om]]
        if (!contrast_name %in% names(de_tables)) next

        df <- de_tables[[contrast_name]]
        id_map <- tryCatch(
            map_feature_ids_to_entrez(de_tables[contrast_name], om, harmonization_res, org_db),
            error = function(e) NULL
        )
        if (is.null(id_map) || nrow(id_map) == 0) next

        df_mapped <- merge(df, id_map, by = "feature_id")
        all_ids <- unique(df_mapped$ENTREZID[!is.na(df_mapped$ENTREZID)])
        sig <- df_mapped$ENTREZID[!is.na(df_mapped$padj) & df_mapped$padj < 0.05]
        if (length(sig) < 5) {
            sig <- df_mapped$ENTREZID[!is.na(df_mapped$pvalue) & df_mapped$pvalue < 0.05]
        }
        sig <- unique(sig[!is.na(sig)])

        if (length(sig) > 0) {
            per_omics_sig[[om]] <- sig
            per_omics_universe[[om]] <- all_ids
        }
    }

    if (length(per_omics_sig) == 0) return(invisible(NULL))

    pooled_sig <- unique(unlist(per_omics_sig))
    pooled_universe <- unique(unlist(per_omics_universe))

    kegg_conv <- convert_entrez_to_kegg(pooled_universe, kegg_org)
    if (is.null(kegg_conv)) return(invisible(NULL))

    pooled_sig_kegg <- unique(kegg_conv[pooled_sig])
    pooled_sig_kegg <- pooled_sig_kegg[!is.na(pooled_sig_kegg)]
    pooled_univ_kegg <- unique(kegg_conv[pooled_universe])
    pooled_univ_kegg <- pooled_univ_kegg[!is.na(pooled_univ_kegg)]

    pooled_ora <- run_multi_ora_kegg(pooled_sig_kegg, pooled_univ_kegg, kegg_org, "pooled")

    per_omics_ora <- list()
    for (om in names(per_omics_sig)) {
        k <- kegg_conv[per_omics_sig[[om]]]
        k <- unique(k[!is.na(k)])
        per_omics_ora[[om]] <- run_multi_ora_kegg(k, pooled_univ_kegg, kegg_org, om)
    }

    # Run per-contrast metabolomics compound ORA if data is available
    contrast_metab_ora <- NULL
    if (!is.null(metab_de_tables)) {
        # Find matching contrast name (try exact, then fuzzy match)
        metab_cname <- if (contrast_name %in% names(metab_de_tables)) {
            contrast_name
        } else {
            # Fuzzy: try normalized contrast names
            norm_cn <- tolower(gsub("[^a-z0-9]", "", contrast_name))
            metab_norms <- tolower(gsub("[^a-z0-9]", "", names(metab_de_tables)))
            idx <- match(norm_cn, metab_norms)
            if (!is.na(idx)) names(metab_de_tables)[idx] else NULL
        }
        if (!is.null(metab_cname)) {
            contrast_metab_ora <- tryCatch({
                de_list <- metab_de_tables[metab_cname]
                id_map <- map_metabolite_ids_to_kegg(de_list, harmonization_res)
                if (!is.null(id_map) && nrow(id_map) > 0) {
                    full_universe <- unique(id_map$KEGG_CPD[!is.na(id_map$KEGG_CPD)])
                    de_df <- de_list[[1]]
                    de_mapped <- merge(de_df, id_map, by = "feature_id")
                    de_mapped$KEGG_ID <- de_mapped$KEGG_CPD
                    run_compound_ora(de_mapped, out_dir, 2, 500, 0.1,
                                     universe = full_universe)
                } else NULL
            }, error = function(e) {
                message("    Per-contrast compound ORA failed for ", contrast_name,
                        ": ", e$message)
                NULL
            })
        }
    }

    combined <- build_multi_ora_summary(pooled_ora, per_omics_ora, contrast_metab_ora)
    if (!is.null(combined) && nrow(combined) > 0) {
        write.csv(combined, file.path(out_dir, "multi_ora_results.csv"), row.names = FALSE)

        if (!is.null(pooled_ora) && nrow(pooled_ora) > 0) {
            bp_path <- file.path(out_dir, "multi_ora_pooled_barplot.png")
            png(bp_path, width = 1000, height = 700, res = 120)
            tryCatch({
                plot_multi_ora_barplot(pooled_ora,
                    paste0("Multi-ORA (", gsub("_", " ", contrast_name), ")"))
            }, error = function(e) {
                plot.new()
                text(0.5, 0.5, paste("Plot failed:", e$message), cex = 1.2)
            })
            dev.off()
        }
    }

    invisible(NULL)
}


#' @param sig_genes Significant KEGG gene IDs
#' @param universe All KEGG gene IDs (shared universe)
#' @param kegg_org KEGG organism code
#' @param label Label for messages
#' @return data.frame with ORA results
run_multi_ora_kegg <- function(sig_genes, universe, kegg_org,
                                label = "pooled", pval_cutoff = 0.1) {

    if (length(sig_genes) < 3) {
        message("    ", label, ": too few significant genes (", length(sig_genes), ")")
        return(NULL)
    }

    # Try clusterProfiler first — use lenient cutoff, filter manually after
    ora_res <- tryCatch({
        res <- clusterProfiler::enrichKEGG(
            gene = sig_genes,
            universe = universe,
            organism = kegg_org,
            keyType = "kegg",
            minGSSize = 5,
            maxGSSize = 500,
            # Open both cutoffs and let this wrapper threshold below; the default
            # qvalueCutoff = 0.2 would otherwise drop rows with raw p < 0.05 but
            # q >= 0.2 before the raw-p fallback sees them (same as the enricher path).
            pvalueCutoff = 1.0,
            qvalueCutoff = 1.0
        )
        if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
            df <- as.data.frame(res)
            out <- data.frame(
                pathway = df$Description,
                ID = df$ID,
                pvalue = df$pvalue,
                padj = df$p.adjust,
                GeneRatio = df$GeneRatio,
                Count = df$Count,
                geneID = df$geneID,
                stringsAsFactors = FALSE
            )
            # Filter: prefer padj, fall back to pvalue < 0.05
            padj_hits <- out[!is.na(out$padj) & out$padj < pval_cutoff, ]
            if (nrow(padj_hits) > 0) return(padj_hits)
            pval_hits <- out[!is.na(out$pvalue) & out$pvalue < 0.05, ]
            if (nrow(pval_hits) > 0) {
                message("    ", label, ": padj too strict, using pvalue < 0.05 (",
                        nrow(pval_hits), " pathways)")
                return(pval_hits)
            }
        }
        NULL
    }, error = function(e) {
        message("    ", label, " clusterProfiler ORA failed: ", e$message)
        NULL
    })

    if (!is.null(ora_res)) {
        message("    ", label, ": ", nrow(ora_res), " enriched pathways")
        return(ora_res)
    }

    # Fallback: Fisher's exact test
    run_ora_kegg_fisher(sig_genes, universe, kegg_org, 5, 500, pval_cutoff)
}


#' Convert a GMT file to clusterProfiler TERM2GENE / TERM2NAME frames
#'
#' Reuses \code{read_gmt()} (a named list of gene vectors carrying a
#' \code{descriptions} attribute) and reshapes it for
#' \code{clusterProfiler::enricher}.
#'
#' @param gmt_file Path to a GMT file.
#' @return A list with \code{t2g} (data.frame term, gene) and \code{t2n}
#'   (data.frame term, name), or NULL if the file has no usable gene sets.
gmt_to_term2gene <- function(gmt_file) {
    gs <- tryCatch(read_gmt(gmt_file), error = function(e) NULL)
    if (is.null(gs) || length(gs) == 0) return(NULL)
    descr <- attr(gs, "descriptions")
    t2g <- data.frame(
        term = rep(names(gs), lengths(gs)),
        gene = unlist(gs, use.names = FALSE),
        stringsAsFactors = FALSE
    )
    t2n <- if (!is.null(descr)) {
        data.frame(term = names(descr), name = unname(descr), stringsAsFactors = FALSE)
    } else NULL
    list(t2g = t2g, t2n = t2n)
}


#' Run ORA with clusterProfiler::enricher against a custom TERM2GENE
#'
#' GMT-based counterpart to \code{run_multi_ora_kegg()}; returns the same
#' data.frame shape so the shared multi-ORA summary and plots consume it
#' unchanged.
#'
#' @param sig_genes Significant feature IDs (native namespace, e.g. EHI_ / XP_).
#' @param universe All tested feature IDs (shared background).
#' @param term2gene data.frame(term, gene) gene-set membership.
#' @param term2name Optional data.frame(term, name) for readable pathway labels.
#' @param label Label used in progress messages.
#' @param pval_cutoff Adjusted-p cutoff (falls back to raw p < 0.05).
#' @return data.frame(pathway, ID, pvalue, padj, GeneRatio, Count, geneID), or
#'   NULL. Every tested set (not only the enriched ones) travels along in an
#'   \code{all_tested} attribute, the same way \code{run_compound_ora()} does
#'   it: the per-collection audit needs the sets that were tested and came out
#'   flat, and they are invisible in the filtered table.
run_multi_ora_enricher <- function(sig_genes, universe, term2gene, term2name = NULL,
                                   label = "pooled", pval_cutoff = 0.1) {

    if (length(sig_genes) < 3) {
        message("    ", label, ": too few significant genes (", length(sig_genes), ")")
        return(NULL)
    }

    ora_res <- tryCatch({
        res <- clusterProfiler::enricher(
            gene = sig_genes,
            universe = universe,
            TERM2GENE = term2gene,
            TERM2NAME = term2name,
            minGSSize = 5,
            maxGSSize = 500,
            # Ask enricher for every tested set (p and q cutoffs both open); this
            # wrapper does its own padj / raw-p thresholding below. Leaving the
            # default qvalueCutoff = 0.2 would silently drop rows with raw
            # p < 0.05 but q >= 0.2 before the raw-p fallback ever sees them.
            pvalueCutoff = 1.0,
            qvalueCutoff = 1.0
        )
        if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
            df <- as.data.frame(res)
            out <- data.frame(
                pathway   = if (!is.null(df$Description)) df$Description else df$ID,
                ID        = df$ID,
                pvalue    = df$pvalue,
                padj      = df$p.adjust,
                GeneRatio = df$GeneRatio,
                Count     = df$Count,
                geneID    = df$geneID,
                stringsAsFactors = FALSE
            )
            padj_hits <- out[!is.na(out$padj) & out$padj < pval_cutoff, ]
            if (nrow(padj_hits) > 0) {
                attr(padj_hits, "all_tested") <- out
                return(padj_hits)
            }
            pval_hits <- out[!is.na(out$pvalue) & out$pvalue < 0.05, ]
            if (nrow(pval_hits) > 0) {
                message("    ", label, ": padj too strict, using pvalue < 0.05 (",
                        nrow(pval_hits), " pathways)")
                attr(pval_hits, "all_tested") <- out
                return(pval_hits)
            }
        }
        NULL
    }, error = function(e) {
        message("    ", label, " enricher ORA failed: ", e$message)
        NULL
    })

    if (!is.null(ora_res)) {
        message("    ", label, ": ", nrow(ora_res), " enriched pathways")
    }
    ora_res
}


#' Summarise a pooled ORA by gene-set collection
#'
#' Answers the question a mixed-collection bar plot cannot: how many sets each
#' GMT contributed, how many of them the ORA could actually test (a set is only
#' tested when at least one hit gene falls in it and it survives the size
#' filter), and how many came out enriched. A collection with zero enriched
#' terms is then visibly absent-because-flat rather than absent-because-dropped.
#'
#' @param term2gene data.frame(term, gene) of the pooled gene-set membership.
#' @param pooled_ora Pooled ORA table from \code{run_multi_ora_enricher()},
#'   carrying its \code{all_tested} attribute; may be NULL.
#' @return data.frame(collection, n_sets, n_tested, n_p05, n_reported), ordered
#'   by descending n_reported.
summarize_ora_collections <- function(term2gene, pooled_ora) {
    sets <- unique(as.character(term2gene$term))
    set_coll <- classify_pathway_collection(sets)
    colls <- sort(unique(set_coll))

    tested <- attr(pooled_ora, "all_tested")
    tested_coll <- if (is.null(tested)) character(0) else classify_pathway_collection(tested$ID)
    p05_coll <- if (is.null(tested)) character(0) else {
        tested_coll[!is.na(tested$pvalue) & tested$pvalue < 0.05]
    }
    rep_coll <- if (is.null(pooled_ora)) character(0) else classify_pathway_collection(pooled_ora$ID)

    n_in <- function(x, coll) vapply(coll, function(cl) sum(x == cl), integer(1))
    out <- data.frame(
        collection  = colls,
        n_sets      = n_in(set_coll, colls),
        n_tested    = n_in(tested_coll, colls),
        n_p05       = n_in(p05_coll, colls),
        n_reported  = n_in(rep_coll, colls),
        stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    out[order(-out$n_reported, -out$n_tested), , drop = FALSE]
}


#' GMT-based multi-omics ORA (fallback when no KEGG/OrgDb organism)
#'
#' Runs over-representation analysis on the pooled and per-omic significant
#' gene-based DE features using each omic's custom GMT
#' (\code{modes.<omic>.pathway.gmt_file}), via \code{clusterProfiler::enricher}.
#' The per-omic GMTs share pathway IDs across namespaces (e.g. EHI_ for RNA and
#' XP_ for proteomics), so their TERM2GENE tables are row-bound into one
#' collection and a pooled mixed-namespace hit list enriches against it — each
#' omic's hits match that omic's members under the same term. A gene measured in
#' two omics contributes once per omic (pooled multi-omic evidence).
#'
#' Per-contrast ORA is intentionally not run here: the multi-omics runs are
#' single-contrast, so pooled + per-omic already covers the whole comparison.
#'
#' @param de_results Named list of DE results per omics.
#' @param harmonization_res Harmonization result (for \code{extract_de_tables}).
#' @param config Full config object.
#' @param out_dir Output directory for results and plots.
#' @return list(pooled, per_omics, metabolomics, combined, collections, plots),
#'   or NULL. \code{collections} is the per-collection audit also written to
#'   multi_ora_collection_summary.csv.
run_multi_ora_gmt <- function(de_results, harmonization_res, config, out_dir) {
    # KEGG classes this organism cannot have; shared with the cross-omics tables
    # and the pathway maps so every figure filters the same way.
    excl_classes <- config$modes$multiomics$enrichment$exclude_pathway_classes

    gene_omics   <- c("transcriptomics", "proteomics")
    omic_cfg_key <- c(transcriptomics = "rna", proteomics = "proteomics")

    per_omics_sig  <- list()
    per_omics_univ <- list()
    per_omics_t2g  <- list()
    per_omics_t2n  <- list()

    for (om in intersect(gene_omics, names(de_results))) {
        # gmt_file may be a single path or a YAML list of paths (GO + KEGG);
        # read_gmt() already merges several files, so only these guards needed
        # to vectorise — with a list they were comparing length-2 vectors and
        # aborting the whole Multi-ORA step.
        gmt_path <- unlist(config$modes[[omic_cfg_key[[om]]]]$pathway$gmt_file,
                           use.names = FALSE)
        if (length(gmt_path) == 0 || !any(nzchar(gmt_path))) next
        # Resolve like every other user-supplied input (metabolomics enrichment,
        # data files): absolute paths pass through, relative ones resolve under
        # the raw/ data dir. resolve_raw_path() would mangle an absolute path.
        gmt_abs <- resolve_input_path(config, gmt_path)
        if (any(!file.exists(gmt_abs))) {
            message("  Multi-ORA (GMT): ", om, " gmt_file not found: ",
                    paste(gmt_abs[!file.exists(gmt_abs)], collapse = ", "))
            next
        }
        gs <- gmt_to_term2gene(gmt_abs)
        if (is.null(gs) || nrow(gs$t2g) == 0) next

        de_tables <- extract_de_tables(de_results[[om]], om, harmonization_res)
        if (is.null(de_tables) || length(de_tables) == 0) next

        sig <- character(0); univ <- character(0)
        for (nm in names(de_tables)) {
            df  <- de_tables[[nm]]
            # Match GMT members by bare ID: RNA feature_ids may carry a "Gene:"
            # prefix (e.g. "Gene:EHI_012345") while custom GMTs list the bare ID.
            # The single-omics RNA pathway module strips it the same way.
            fid <- sub("^Gene:", "", as.character(df$feature_id))
            univ <- c(univ, fid)
            s <- fid[!is.na(df$padj) & df$padj < 0.05]
            if (length(s) < 5) s <- fid[!is.na(df$pvalue) & df$pvalue < 0.05]
            sig <- c(sig, s)
        }
        sig  <- unique(sig[!is.na(sig)])
        univ <- unique(univ[!is.na(univ)])
        if (length(sig) == 0) next

        per_omics_sig[[om]]  <- sig
        per_omics_univ[[om]] <- univ
        per_omics_t2g[[om]]  <- gs$t2g
        per_omics_t2n[[om]]  <- gs$t2n
        message("  Multi-ORA (GMT) ", om, ": ", length(sig), " sig / ",
                length(univ), " tested features")
    }

    # A cross-omics ORA needs at least two gene-based GMT omics. With only one
    # (e.g. a non-model RNA+metabolomics run, where metabolomics is not covered
    # by this gene-based fallback) the pooled result would be single-omic while
    # presenting as multi-omics, so decline instead.
    if (length(per_omics_sig) < 2) {
        message("Multi-ORA (GMT): need >= 2 gene-based omics with a usable GMT + ",
                "sig features; found ", length(per_omics_sig), " -- skipping")
        return(NULL)
    }

    comb_t2g    <- unique(do.call(rbind, per_omics_t2g))
    comb_t2n    <- unique(do.call(rbind, Filter(Negate(is.null), per_omics_t2n)))
    pooled_sig  <- unique(unlist(per_omics_sig))
    pooled_univ <- unique(unlist(per_omics_univ))

    message("  Running pooled GMT ORA...")
    pooled_ora <- run_multi_ora_enricher(pooled_sig, pooled_univ, comb_t2g, comb_t2n, "pooled")

    per_omics_ora <- list()
    for (om in names(per_omics_sig)) {
        message("  Running per-omics GMT ORA for ", om, "...")
        per_omics_ora[[om]] <- run_multi_ora_enricher(
            per_omics_sig[[om]], per_omics_univ[[om]],
            per_omics_t2g[[om]], per_omics_t2n[[om]], om)
    }

    combined <- build_multi_ora_summary(pooled_ora, per_omics_ora, NULL)
    if (is.null(combined) || nrow(combined) == 0) {
        message("Multi-ORA (GMT): no enriched pathways found")
        return(NULL)
    }
    write.csv(combined, file.path(out_dir, "multi_ora_results.csv"), row.names = FALSE)
    message("  Multi-ORA (GMT) found ", nrow(combined), " enriched pathways (pooled)")

    coll_summary <- summarize_ora_collections(comb_t2g, pooled_ora)
    write.csv(coll_summary, file.path(out_dir, "multi_ora_collection_summary.csv"),
              row.names = FALSE)
    message("  Multi-ORA (GMT) per collection (enriched/tested/sets): ",
            paste(sprintf("%s %d/%d/%d", coll_summary$collection, coll_summary$n_reported,
                          coll_summary$n_tested, coll_summary$n_sets), collapse = "; "))

    plots <- list()
    if (!is.null(pooled_ora) && nrow(pooled_ora) > 0) {
        plots$pooled_barplot <- file.path(out_dir, "multi_ora_pooled_barplot.png")
        png(plots$pooled_barplot, width = 1000, height = 700, res = 120)
        tryCatch(
            plot_multi_ora_barplot(pooled_ora, "Pooled Multi-ORA (GMT gene sets)",
                                   exclude_classes = excl_classes, cache_dir = out_dir),
            error = function(e) { plot.new(); text(0.5, 0.5, paste("Plot failed:", e$message), cex = 1.2) }
        )
        dev.off()
    }
    plots$dotplot <- file.path(out_dir, "multi_ora_dotplot.png")
    tryCatch(plot_multi_ora_dotplot(combined, per_omics_ora, NULL, out_dir),
             error = function(e) message("  Multi-ORA dot plot failed: ", e$message))
    plots$support_barplot <- file.path(out_dir, "multi_ora_support_barplot.png")
    tryCatch(plot_multi_ora_support(combined, out_dir),
             error = function(e) message("  Multi-ORA support plot failed: ", e$message))

    list(pooled = pooled_ora, per_omics = per_omics_ora, metabolomics = NULL,
         combined = combined, collections = coll_summary, plots = plots)
}


#' Build combined multi-ORA summary table
#'
#' Merges pooled ORA with per-omics results into a single table showing
#' which omics support each pathway. Uses full outer join so compound-only
#' pathways from metabolomics are also included.
build_multi_ora_summary <- function(pooled_ora, per_omics_ora, metab_ora) {

    if ((is.null(pooled_ora) || nrow(pooled_ora) == 0) &&
        (is.null(metab_ora) || nrow(metab_ora) == 0)) return(NULL)

    # Start from gene-based pooled ORA if available
    if (!is.null(pooled_ora) && nrow(pooled_ora) > 0) {
        summary <- pooled_ora[, c("pathway", "ID", "pvalue", "padj", "Count"), drop = FALSE]
        colnames(summary)[colnames(summary) == "pvalue"] <- "pooled_pvalue"
        colnames(summary)[colnames(summary) == "padj"] <- "pooled_padj"
        colnames(summary)[colnames(summary) == "Count"] <- "pooled_count"
    } else {
        summary <- data.frame(
            pathway = character(0), ID = character(0),
            pooled_pvalue = numeric(0), pooled_padj = numeric(0),
            pooled_count = integer(0), stringsAsFactors = FALSE
        )
    }

    # Add normalized ID column for cross-prefix matching
    summary$norm_id <- normalize_kegg_pathway_id(summary$ID)

    # Add per-omics support columns (padj and pvalue)
    omics_names <- names(per_omics_ora)
    for (om in omics_names) {
        om_res <- per_omics_ora[[om]]
        if (is.null(om_res) || nrow(om_res) == 0) {
            summary[[paste0(om, "_padj")]] <- NA_real_
            summary[[paste0(om, "_pvalue")]] <- NA_real_
            next
        }
        om_padj <- setNames(om_res$padj, om_res$ID)
        summary[[paste0(om, "_padj")]] <- om_padj[summary$ID]
        if ("pvalue" %in% colnames(om_res)) {
            om_pval <- setNames(om_res$pvalue, om_res$ID)
            summary[[paste0(om, "_pvalue")]] <- om_pval[summary$ID]
        }
    }

    # Add metabolomics via full outer join on normalized pathway IDs
    if (!is.null(metab_ora) && nrow(metab_ora) > 0 && "ID" %in% colnames(metab_ora)) {
        met_norm <- normalize_kegg_pathway_id(metab_ora$ID)
        met_padj <- setNames(metab_ora$padj, met_norm)
        met_pval <- setNames(metab_ora$pvalue, met_norm)

        # Match existing rows
        summary$metabolomics_padj <- met_padj[summary$norm_id]
        summary$metabolomics_pvalue <- met_pval[summary$norm_id]

        # Find compound-only pathways (not in gene ORA results)
        new_ids <- setdiff(met_norm, summary$norm_id)
        if (length(new_ids) > 0) {
            new_rows <- data.frame(
                pathway = metab_ora$pathway[match(new_ids, met_norm)],
                ID = metab_ora$ID[match(new_ids, met_norm)],
                pooled_pvalue = NA_real_,
                pooled_padj = NA_real_,
                pooled_count = NA_integer_,
                norm_id = new_ids,
                stringsAsFactors = FALSE
            )
            # Add per-omics columns as NA
            for (om in omics_names) {
                new_rows[[paste0(om, "_padj")]] <- NA_real_
                new_rows[[paste0(om, "_pvalue")]] <- NA_real_
            }
            new_rows$metabolomics_padj <- met_padj[new_ids]
            new_rows$metabolomics_pvalue <- met_pval[new_ids]

            summary <- rbind(summary, new_rows)
        }
    }

    # Count supporting omics (padj < 0.05)
    padj_cols <- grep("_padj$", colnames(summary), value = TRUE)
    padj_cols <- setdiff(padj_cols, "pooled_padj")
    padj_mat <- sapply(padj_cols, function(col) !is.na(summary[[col]]) & summary[[col]] < 0.05)
    if (!is.matrix(padj_mat)) padj_mat <- matrix(padj_mat, ncol = length(padj_cols))
    summary$n_omics_support <- rowSums(padj_mat)

    # Fallback: count support using pvalue < 0.05 when padj is too strict
    pval_cols <- grep("_pvalue$", colnames(summary), value = TRUE)
    pval_cols <- setdiff(pval_cols, "pooled_pvalue")
    if (length(pval_cols) > 0) {
        pval_mat <- sapply(pval_cols, function(col) !is.na(summary[[col]]) & summary[[col]] < 0.05)
        if (!is.matrix(pval_mat)) pval_mat <- matrix(pval_mat, ncol = length(pval_cols))
        summary$n_omics_support_pval <- rowSums(pval_mat)
    } else {
        summary$n_omics_support_pval <- summary$n_omics_support
    }

    # Sort: gene ORA pathways first (by pooled_pvalue), then compound-only.
    # The GMT fallback passes metab_ora = NULL, so this column may not exist —
    # order() aborts on a NULL sort key.
    metab_sort <- if ("metabolomics_padj" %in% colnames(summary)) {
        summary$metabolomics_padj
    } else {
        rep(NA_real_, nrow(summary))
    }
    summary <- summary[order(is.na(summary$pooled_pvalue), summary$pooled_pvalue,
                             metab_sort), ]

    # Drop internal helper column
    summary$norm_id <- NULL
    summary
}


#' Classify pathway IDs by the gene-set collection they come from
#'
#' Multi-ORA pools gene sets from several GMTs whose IDs live in different
#' namespaces (GO:0006955, map00010, PF00001). Once the TERM2GENE tables are
#' row-bound nothing downstream can tell them apart, so the collection is
#' recovered from the shape of the ID itself.
#'
#' @param ids Character vector of pathway / gene-set IDs.
#' @return Character vector the same length as \code{ids}, each element one of
#'   "GO", "KEGG", "Pfam", "InterPro" or "Other".
classify_pathway_collection <- function(ids) {
    ids <- as.character(ids)
    out <- rep("Other", length(ids))
    out[grepl("^GO:?[0-9]+", ids)] <- "GO"
    # KEGG pathway ids are a 2-4 letter prefix plus 5 digits: the reference maps
    # ("map00010", "ko00010") and every organism code ("hsa04110", "cel04110").
    out[grepl("^[a-z]{2,4}[0-9]{5}$", ids)] <- "KEGG"
    out[grepl("^PF[0-9]{5}", ids)] <- "Pfam"
    out[grepl("^IPR[0-9]{6}", ids)] <- "InterPro"
    out[is.na(ids) | !nzchar(ids)] <- "Other"
    out
}


#' Pick top ORA terms while keeping every collection represented
#'
#' Draws terms round-robin from each collection in ascending p-value. A pooled
#' ORA over GO + KEGG + Pfam is dominated by GO purely by set count (thousands
#' of GO terms against a few hundred KEGG maps), so a plain top-n by p-value
#' leaves the plot looking like a GO-only analysis. Collections enter the
#' rotation best-first, so with fewer slots than collections the strongest ones
#' are still the ones shown. With a single collection this is exactly the plain
#' top-n.
#'
#' @param ora_df ORA table carrying \code{pvalue} and \code{collection} columns.
#' @param top_n Maximum number of rows to keep.
#' @return A subset of \code{ora_df}, ordered by ascending p-value.
select_top_ora_per_collection <- function(ora_df, top_n = 20) {
    # Rank on adjusted p where it exists. Ranking on raw p put terms at the top
    # that multiple testing does not support, and the rest of the report moved to
    # FDR; a figure ordered on a different statistic invites the wrong comparison.
    rank_by <- if ("padj" %in% names(ora_df) && any(is.finite(ora_df$padj))) {
        ora_df$padj
    } else ora_df$pvalue

    ord <- order(rank_by)
    df <- ora_df[ord, , drop = FALSE]
    rank_by <- rank_by[ord]
    if (nrow(df) <= top_n) return(df)

    by_coll <- split(seq_len(nrow(df)), df$collection)
    by_coll <- by_coll[order(vapply(by_coll, function(i) rank_by[i[1]], numeric(1)))]

    keep <- integer(0)
    rank <- 1L
    while (length(keep) < top_n && any(lengths(by_coll) >= rank)) {
        for (idx in by_coll) {
            if (length(keep) >= top_n) break
            if (length(idx) >= rank) keep <- c(keep, idx[rank])
        }
        rank <- rank + 1L
    }
    df[sort(keep), , drop = FALSE]
}


#' Plot multi-ORA pooled barplot
#'
#' Horizontal bar plot of the strongest enriched terms. When the table mixes
#' gene-set collections each bar is tagged and coloured by its collection, so
#' GO terms, KEGG maps and Pfam domains are never silently pooled into one
#' anonymous ranking.
#'
#' A collection that was searched but produced no enriched term cannot have a
#' bar, and an absent bar reads the same as an unconfigured GMT. Where the table
#' carries the \code{all_tested} attribute \code{run_multi_ora_enricher()}
#' attaches, those collections are named in the caption instead of vanishing.
#'
#' @param ora_df ORA table with \code{pathway}, \code{ID} and \code{pvalue},
#'   optionally carrying an \code{all_tested} attribute of every tested set.
#' @param title Plot title.
#' @param top_n Maximum number of terms to draw (default 20).
#' @return Invisible character vector of the collections drawn.
plot_multi_ora_barplot <- function(ora_df, title, top_n = 20,
                                   exclude_classes = NULL, cache_dir = NULL) {
    ids <- if ("ID" %in% colnames(ora_df)) ora_df$ID else ora_df$pathway
    tested <- attr(ora_df, "all_tested")

    # KEGG reference maps are pan-species, so an organism with no KEGG code
    # scores vertebrate-organ and human-disease maps purely on generic orthologs.
    # The cross-omics tables and the pathway maps already drop those classes;
    # without this the pooled barplot was the one figure still showing
    # "Viral myocarditis" and "Bile secretion" for a parasitoid wasp.
    if (!is.null(exclude_classes) && length(exclude_classes) > 0) {
        keep <- keep_kegg_pathways(ids, exclude = unlist(exclude_classes),
                                   cache_dir = cache_dir,
                                   label = "pooled Multi-ORA pathways")
        ora_df <- ora_df[keep, , drop = FALSE]
        ids <- ids[keep]
        if (nrow(ora_df) == 0) {
            plot.new()
            text(0.5, 0.5, "No pathways left after KEGG class filtering", cex = 1.1)
            return(invisible(NULL))
        }
    }

    ora_df$collection <- classify_pathway_collection(ids)
    df <- select_top_ora_per_collection(ora_df, top_n)

    collections <- unique(df$collection)
    mixed <- length(collections) > 1

    tested_coll <- if (is.null(tested) || !"ID" %in% colnames(tested)) {
        character(0)
    } else {
        classify_pathway_collection(tested$ID)
    }
    flat <- setdiff(unique(tested_coll), collections)

    # Tag width has to come out of the name, not be added on top of it: the
    # left margin below is fixed and a longer string just runs off the device.
    name_width <- if (mixed) 42L else 50L
    label <- ifelse(nchar(df$pathway) > name_width,
                    paste0(substr(df$pathway, 1, name_width - 3), "..."),
                    df$pathway)
    if (mixed) label <- paste0("[", df$collection, "] ", label)

    palette <- rep_len(c("#7B2D8E", "#2C7FB8", "#D95F02", "#1B9E77", "#666666"),
                       length(collections))
    bar_col <- if (mixed) palette[match(df$collection, collections)] else "#7B2D8E"

    # Plot whatever the ranking used, so bar length and bar order agree.
    use_padj <- "padj" %in% names(df) && any(is.finite(df$padj))
    stat_vals <- if (use_padj) df$padj else df$pvalue
    stat_label <- if (use_padj) "-log10(FDR)" else "-log10(p-value)"
    neg_log_p <- -log10(stat_vals + 1e-300)
    neg_log_p <- pmin(neg_log_p, 15)

    caption <- character(0)
    if (mixed) {
        caption <- c(caption, paste0("terms drawn in turn from each collection (",
                                     paste(collections, collapse = ", "), ")"))
    }
    if (length(flat) > 0) {
        n_flat <- vapply(flat, function(cl) sum(tested_coll == cl), integer(1))
        caption <- c(caption, paste0(
            "searched, no enriched term: ",
            paste(sprintf("%s (%d sets tested)", flat, n_flat), collapse = ", ")))
    }

    # One caption per line, and one margin line per caption: joined into a single
    # `sub =` string the tail runs off the right edge of the device at this width.
    par(mar = c(4 + max(length(caption), 1L), 17, 3, 2))
    barplot(rev(neg_log_p), horiz = TRUE, names.arg = rev(label),
            las = 1, cex.names = 0.65, col = rev(bar_col),
            xlab = stat_label,
            main = title)
    for (i in seq_along(caption)) {
        mtext(caption[i], side = 1, line = 3 + i, cex = 0.7)
    }
    abline(v = -log10(0.05), col = "red", lty = 2)
    if (mixed) {
        legend("bottomright", legend = collections, fill = palette,
               bty = "n", cex = 0.75)
    }
    invisible(collections)
}


#' Plot multi-ORA dot plot comparing per-omics and pooled results
plot_multi_ora_dotplot <- function(combined, per_omics_ora, metab_ora, out_dir, top_n = 20) {

    if (is.null(combined) || nrow(combined) == 0) return(invisible(NULL))

    top <- combined[seq_len(min(top_n, nrow(combined))), ]

    # Collect per-omics data for these pathways
    plot_data <- list()

    # Pooled
    plot_data[["Pooled"]] <- data.frame(
        pathway = top$pathway,
        source = "Pooled",
        neg_log10_padj = -log10(top$pooled_padj + 1e-300),
        count = top$pooled_count,
        stringsAsFactors = FALSE
    )

    # Per-omics
    for (om in names(per_omics_ora)) {
        om_res <- per_omics_ora[[om]]
        if (is.null(om_res) || nrow(om_res) == 0) next
        padj_map <- setNames(om_res$padj, om_res$pathway)
        count_map <- setNames(om_res$Count, om_res$pathway)

        matched_padj <- padj_map[top$pathway]
        matched_count <- count_map[top$pathway]

        plot_data[[om]] <- data.frame(
            pathway = top$pathway,
            source = gsub("_", " ", tools::toTitleCase(om)),
            neg_log10_padj = -log10(ifelse(is.na(matched_padj), 1, matched_padj) + 1e-300),
            count = ifelse(is.na(matched_count), 0, matched_count),
            stringsAsFactors = FALSE
        )
    }

    # Metabolomics compound ORA (join by normalized KEGG pathway ID)
    if (!is.null(metab_ora) && nrow(metab_ora) > 0) {
        met_padj_map <- setNames(metab_ora$padj, normalize_kegg_pathway_id(metab_ora$ID))
        met_count_map <- setNames(metab_ora$setSize, normalize_kegg_pathway_id(metab_ora$ID))
        norm_top_ids <- normalize_kegg_pathway_id(top$ID)

        matched_padj <- met_padj_map[norm_top_ids]
        matched_count <- met_count_map[norm_top_ids]

        plot_data[["Metabolomics"]] <- data.frame(
            pathway = top$pathway,
            source = "Metabolomics",
            neg_log10_padj = -log10(ifelse(is.na(matched_padj), 1, matched_padj) + 1e-300),
            count = ifelse(is.na(matched_count), 0, matched_count),
            stringsAsFactors = FALSE
        )
    }

    plot_df <- do.call(rbind, plot_data)
    plot_df$neg_log10_padj <- pmin(plot_df$neg_log10_padj, 15)
    # Remove entries where pathway was not found in the per-omics results
    plot_df <- plot_df[plot_df$neg_log10_padj > 0, ]

    if (nrow(plot_df) == 0) return(invisible(NULL))

    # Truncate long pathway names
    plot_df$pathway_short <- ifelse(
        nchar(plot_df$pathway) > 45,
        paste0(substr(plot_df$pathway, 1, 42), "..."),
        plot_df$pathway
    )

    # ggplot dot plot
    source_colors <- c(
        "Pooled" = "#7B2D8E",
        "Transcriptomics" = "#E41A1C",
        "Proteomics" = "#377EB8",
        "Metabolomics" = "#4DAF4A"
    )

    p <- ggplot2::ggplot(plot_df,
        ggplot2::aes(
            x = neg_log10_padj,
            y = stats::reorder(pathway_short, neg_log10_padj),
            color = source,
            size = count
        )) +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::scale_color_manual(values = source_colors, name = "Source") +
        ggplot2::scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
        ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.5) +
        ggplot2::labs(
            title = "Multi-ORA: Pooled vs Per-Omics Enrichment",
            x = "-log10(padj)",
            y = NULL
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 7),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
        )

    ggplot2::ggsave(
        file.path(out_dir, "multi_ora_dotplot.png"),
        plot = p, width = 12, height = max(6, 2 + top_n * 0.3), dpi = 300
    )
    message("  Saved multi-ORA dot plot")
}


#' Plot multi-ORA omics support barplot
#'
#' Shows how many omics layers support each enriched pathway.
plot_multi_ora_support <- function(combined, out_dir, top_n = 25) {
    if (is.null(combined) || nrow(combined) == 0) return(invisible(NULL))

    top <- combined[seq_len(min(top_n, nrow(combined))), ]
    top$label <- ifelse(nchar(top$pathway) > 45,
                         paste0(substr(top$pathway, 1, 42), "..."),
                         top$pathway)

    support_colors <- c("0" = "grey70", "1" = "#FDB863", "2" = "#E66101", "3" = "#B35806")

    p <- ggplot2::ggplot(top,
        ggplot2::aes(
            x = -log10(pooled_pvalue + 1e-300),
            y = stats::reorder(label, -pooled_pvalue),
            fill = factor(n_omics_support)
        )) +
        ggplot2::geom_col(alpha = 0.85) +
        ggplot2::scale_fill_manual(values = support_colors, name = "# Omics\nSupporting") +
        ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.5) +
        ggplot2::labs(
            title = "Multi-ORA: Pathway Enrichment with Omics Support",
            x = "-log10(p-value)",
            y = NULL
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 7),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
        )

    ggplot2::ggsave(
        file.path(out_dir, "multi_ora_support_barplot.png"),
        plot = p, width = 12, height = max(6, 2 + top_n * 0.25), dpi = 300
    )
    message("  Saved multi-ORA support barplot")
}


#' Generate pathview overlays for multi-omics-supported KEGG pathways
#'
#' For pathways enriched in >= min_support omics, renders KEGG pathway maps
#' with per-omics logFC coloring and compiles them into a single PDF.
#'
#' @param combined Multi-ORA summary table (from build_multi_ora_summary)
#' @param de_results Named list of DE results per omics
#' @param harmonization_res Harmonization result
#' @param config Full config
#' @param out_dir Output directory for pathview files
#' @param min_support Minimum n_omics_support to include (default 2)
#' @return Character path to compiled PDF, or NULL
generate_multi_ora_pathview <- function(combined, de_results, harmonization_res,
                                        config, out_dir, min_support = 2,
                                        top_n = 5) {

    if (!requireNamespace("pathview", quietly = TRUE)) {
        message("  Package 'pathview' not installed. Skipping pathway maps.")
        return(NULL)
    }

    if (is.null(combined) || !"n_omics_support" %in% colnames(combined)) return(NULL)

    # Read top_n from config if available
    pv_cfg <- config$modes$multiomics$enrichment$pathview %||% list()
    top_n <- pv_cfg$top_n %||% top_n

    # --- Prioritize pathways supported by both metabolomics AND proteomics ---
    has_met_padj <- "metabolomics_padj" %in% colnames(combined)
    has_prot_padj <- "proteomics_padj" %in% colnames(combined)
    has_met_pval <- "metabolomics_pvalue" %in% colnames(combined)
    has_prot_pval <- "proteomics_pvalue" %in% colnames(combined)

    supported <- NULL

    # Tier 1: both metabolomics + proteomics padj < 0.05
    if (has_met_padj && has_prot_padj) {
        tier1 <- combined[
            !is.na(combined$metabolomics_padj) & combined$metabolomics_padj < 0.05 &
            !is.na(combined$proteomics_padj) & combined$proteomics_padj < 0.05, ]
        if (nrow(tier1) > 0) {
            message("  Pathview: ", nrow(tier1),
                    " pathways supported by both metabolomics & proteomics (padj < 0.05)")
            supported <- tier1
        }
    }

    # Tier 2: both metabolomics + proteomics pvalue < 0.05
    if (is.null(supported) || nrow(supported) == 0) {
        if (has_met_pval && has_prot_pval) {
            tier2 <- combined[
                !is.na(combined$metabolomics_pvalue) & combined$metabolomics_pvalue < 0.05 &
                !is.na(combined$proteomics_pvalue) & combined$proteomics_pvalue < 0.05, ]
            if (nrow(tier2) > 0) {
                message("  Pathview: padj too strict; ", nrow(tier2),
                        " pathways with metabolomics & proteomics pvalue < 0.05")
                supported <- tier2
            }
        }
    }

    # Tier 3: general n_omics_support >= min_support (padj-based)
    if (is.null(supported) || nrow(supported) == 0) {
        supported <- combined[combined$n_omics_support >= min_support, ]
    }

    # Tier 4: general pvalue-based support
    if (nrow(supported) == 0 && "n_omics_support_pval" %in% colnames(combined)) {
        message("  No pathways with padj-based support >= ", min_support,
                ", falling back to pvalue < 0.05")
        supported <- combined[combined$n_omics_support_pval >= min_support, ]
    }

    # Tier 5: relax to single-omics support with pvalue
    if (nrow(supported) == 0 && min_support > 1) {
        message("  No pathways with pvalue-based support >= ", min_support,
                ", relaxing to single-omics support")
        if ("n_omics_support_pval" %in% colnames(combined)) {
            supported <- combined[combined$n_omics_support_pval >= 1, ]
        } else {
            supported <- combined[combined$n_omics_support >= 1, ]
        }
    }

    if (nrow(supported) == 0) {
        message("  No pathways with any omics support for pathview")
        return(NULL)
    }

    # Limit to top_n pathways sorted by pooled_pvalue
    supported <- supported[order(supported$pooled_pvalue, na.last = TRUE), ]
    supported <- head(supported, top_n)

    message("  Generating pathview for ", nrow(supported),
            " pathways (>= ", min_support, " omics support)")

    organism <- config$global$organism %||% "human"
    kegg_org <- get_kegg_organism(organism)
    org_db <- get_organism_db(organism)

    if (is.null(kegg_org)) return(NULL)

    pv_dir <- file.path(out_dir, "pathview")
    dir.create(pv_dir, recursive = TRUE, showWarnings = FALSE)

    # --- Build per-contrast gene logFC and compound logFC ---
    # Extract DE tables and ID maps once per omics
    gene_de_tables <- list()
    gene_id_maps <- list()
    for (om in intersect(c("transcriptomics", "proteomics"), names(de_results))) {
        de_tables <- extract_de_tables(de_results[[om]], om, harmonization_res)
        if (is.null(de_tables) || length(de_tables) == 0) next
        id_map <- tryCatch(
            map_feature_ids_to_entrez(de_tables, om, harmonization_res, org_db),
            error = function(e) NULL
        )
        if (is.null(id_map) || nrow(id_map) == 0) next
        gene_de_tables[[om]] <- de_tables
        gene_id_maps[[om]] <- id_map
    }

    metab_de_tables <- NULL
    metab_id_map <- NULL
    if ("metabolomics" %in% names(de_results)) {
        metab_de_tables <- extract_de_tables(de_results$metabolomics, "metabolomics",
                                              harmonization_res)
        metab_id_map <- tryCatch(
            map_metabolite_ids_to_kegg(metab_de_tables, harmonization_res),
            error = function(e) NULL
        )
    }

    # Determine contrast names (from first omics with multiple contrasts)
    contrast_names <- NULL
    for (om in names(gene_de_tables)) {
        nms <- names(gene_de_tables[[om]])
        if (!is.null(nms) && length(nms) > 0) {
            contrast_names <- nms
            break
        }
    }
    if (is.null(contrast_names)) {
        # Fallback to metabolomics contrast names
        if (!is.null(metab_de_tables)) contrast_names <- names(metab_de_tables)
    }
    if (is.null(contrast_names)) contrast_names <- "contrast_1"

    # pathview requires its 'bods' dataset in the global environment
    if (!exists("bods", envir = globalenv())) {
        utils::data("bods", package = "pathview", envir = globalenv())
    }

    cwd <- getwd()
    setwd(pv_dir)
    on.exit(setwd(cwd), add = TRUE)

    all_generated_pngs <- list()

    for (ci in seq_along(contrast_names)) {
        contrast <- contrast_names[ci]
        safe_contrast <- make.names(contrast)
        message("  Pathview for contrast: ", contrast)

        # Build gene logFC for this contrast
        gene_fc_list <- list()
        for (om in names(gene_de_tables)) {
            de_tbl <- gene_de_tables[[om]]
            idx <- min(ci, length(de_tbl))
            df <- de_tbl[[idx]]
            df_mapped <- merge(df, gene_id_maps[[om]], by = "feature_id")
            fc_arr <- tapply(df_mapped$log2fc, df_mapped$ENTREZID, mean, na.rm = TRUE)
            fc_vec <- as.numeric(fc_arr)
            names(fc_vec) <- names(fc_arr)
            gene_fc_list[[om]] <- fc_vec
        }

        gene_data <- NULL
        if (length(gene_fc_list) == 2) {
            all_genes <- unique(c(names(gene_fc_list[[1]]), names(gene_fc_list[[2]])))
            gene_data <- matrix(NA, nrow = length(all_genes), ncol = 2,
                                dimnames = list(all_genes, names(gene_fc_list)))
            for (i in seq_along(gene_fc_list)) {
                fc <- gene_fc_list[[i]]
                gene_data[names(fc), i] <- fc
            }
        } else if (length(gene_fc_list) == 1) {
            gene_data <- gene_fc_list[[1]]
        }

        # Build compound logFC for this contrast
        cpd_data <- NULL
        if (!is.null(metab_id_map) && nrow(metab_id_map) > 0 &&
            !is.null(metab_de_tables) && length(metab_de_tables) > 0) {
            idx <- min(ci, length(metab_de_tables))
            df <- metab_de_tables[[idx]]
            df_mapped <- merge(df, metab_id_map, by = "feature_id")
            cpd_fc <- tapply(df_mapped$log2fc, df_mapped$KEGG_CPD, mean, na.rm = TRUE)
            cpd_data <- as.numeric(cpd_fc)
            names(cpd_data) <- names(cpd_fc)
        }

        if (is.null(gene_data) && is.null(cpd_data)) next

        # Run pathview per pathway for this contrast
        out_suffix <- paste0("multi_ora_", safe_contrast)
        contrast_pngs <- character(0)

        for (i in seq_len(nrow(supported))) {
            pid <- supported$ID[i]
            pw_name <- supported$pathway[i]
            clean_pid <- normalize_kegg_pathway_id(pid)

            tryCatch({
                pathview::pathview(
                    gene.data  = gene_data,
                    cpd.data   = cpd_data,
                    pathway.id = clean_pid,
                    species    = kegg_org,
                    out.suffix = out_suffix,
                    kegg.dir   = pv_dir,
                    keys.align = "y",
                    match.data = TRUE,
                    multi.state = !is.null(dim(gene_data)) && ncol(gene_data) > 1,
                    same.layer = FALSE
                )

                candidates <- c(
                    paste0(kegg_org, clean_pid, ".", out_suffix, ".multi.png"),
                    paste0(kegg_org, clean_pid, ".", out_suffix, ".png")
                )
                found <- candidates[file.exists(candidates)]
                if (length(found) > 0) {
                    contrast_pngs <- c(contrast_pngs, file.path(pv_dir, found[1]))
                    message("    Pathview: ", pw_name, " (", pid, ")")
                }
            }, error = function(e) {
                message("    Pathview failed for ", pid, ": ", e$message)
            })
        }

        if (length(contrast_pngs) > 0) {
            all_generated_pngs[[contrast]] <- contrast_pngs
        }
    }

    if (length(all_generated_pngs) == 0) {
        message("  No pathview plots generated")
        return(NULL)
    }

    # --- Compile into a single PDF with contrast labels ---
    pdf_path <- file.path(out_dir, "multi_ora_pathview_supported.pdf")
    tryCatch({
        grDevices::pdf(pdf_path, width = 12, height = 8)
        for (contrast in names(all_generated_pngs)) {
            pngs <- all_generated_pngs[[contrast]]
            # Title page for this contrast
            grid::grid.newpage()
            grid::grid.text(
                paste0("Contrast: ", contrast),
                gp = grid::gpar(fontsize = 24, fontface = "bold")
            )
            for (png_file in pngs) {
                img <- png::readPNG(png_file)
                grid::grid.newpage()
                # Add contrast label at top
                grid::grid.text(
                    contrast, x = 0.5, y = 0.98,
                    gp = grid::gpar(fontsize = 10, col = "grey40")
                )
                grid::grid.raster(img, y = 0.48, height = 0.92)
            }
        }
        grDevices::dev.off()
        total <- sum(lengths(all_generated_pngs))
        message("  Compiled ", total, " pathview maps (",
                length(all_generated_pngs), " contrasts) into: ", basename(pdf_path))
        pdf_path
    }, error = function(e) {
        message("  PDF compilation failed: ", e$message)
        tryCatch(grDevices::dev.off(), error = function(e2) NULL)
        NULL
    })
}


#' Load the feature -> KEGG Orthology map used for KO-space pathview
#'
#' Non-model organisms have no KEGG organism code, so their features can only be
#' placed onto KEGG reference maps through KO ids. The map is a user-supplied TSV
#' (see utils/build_spalangia_ko_map.R for how one is built) in long format, with
#' one row per (omics, feature, KO) triple.
#'
#' @param config Full config; reads `modes$multiomics$enrichment$pathview$ko_map`.
#' @return Data frame with columns `omics`, `feature_id`, `KO`, or NULL when the
#'   key is absent, the file is missing, or the columns do not line up.
load_feature_ko_map <- function(config) {
    pv_cfg <- config$modes$multiomics$enrichment$pathview %||% list()
    ko_path <- pv_cfg$ko_map %||% NULL
    if (is.null(ko_path) || !nzchar(ko_path[1])) return(NULL)

    # Resolved like every other user-supplied input (gmt_file, data files):
    # absolute paths pass through, relative ones land under the raw/ data dir.
    ko_abs <- resolve_input_path(config, ko_path)[1]
    if (!file.exists(ko_abs)) {
        message("  Union pathview: ko_map not found: ", ko_abs)
        return(NULL)
    }

    ko_map <- tryCatch(
        utils::read.delim(ko_abs, sep = "\t", header = TRUE,
                          colClasses = "character", stringsAsFactors = FALSE),
        error = function(e) NULL
    )
    if (is.null(ko_map) || !all(c("omics", "feature_id", "KO") %in% names(ko_map))) {
        message("  Union pathview: ko_map needs columns omics, feature_id, KO: ", ko_abs)
        return(NULL)
    }
    ko_map <- ko_map[!is.na(ko_map$KO) & nzchar(ko_map$KO), c("omics", "feature_id", "KO")]
    message("  Union pathview: KO map with ", nrow(ko_map), " (feature, KO) rows")
    ko_map
}


#' Aggregate feature-level log2FC onto KEGG Orthology nodes
#'
#' A KO box on a KEGG reference map stands for an ortholog group, so several
#' features of one layer (paralogues, or a protein group that collapses onto the
#' same ortholog) can land on the same box. The node then carries their mean.
#'
#' @param de_table Standardized DE table with `feature_id` and `log2fc` columns.
#' @param ko_map Feature -> KO map from \code{load_feature_ko_map}.
#' @param omics Omics label selecting the rows of `ko_map` to join against.
#' @return Named numeric vector of log2FC per KO id, or NULL when nothing joins.
aggregate_log2fc_by_ko <- function(de_table, ko_map, omics) {
    if (is.null(de_table) || is.null(ko_map)) return(NULL)
    m <- ko_map[ko_map$omics == omics, c("feature_id", "KO"), drop = FALSE]
    if (nrow(m) == 0) return(NULL)
    joined <- merge(de_table, m, by = "feature_id")
    ok <- !is.na(joined$KO) & nzchar(joined$KO) & is.finite(joined$log2fc)
    if (!any(ok)) return(NULL)
    fc <- tapply(joined$log2fc[ok], joined$KO[ok], mean, na.rm = TRUE)
    stats::setNames(as.numeric(fc), names(fc))
}


#' Choose the emptiest corner of a KEGG map for pathview's colour key
#'
#' pathview draws its colour key at a fixed corner (default "topright"), which on
#' some maps lands squarely on the artwork -- on ko01040 it covered the compound
#' legend panel and made both unreadable. The blank KEGG template is cached on
#' disk, so the ink in each corner can be measured and the quietest one chosen.
#'
#' @param template_png Path to the downloaded blank KEGG map PNG.
#' @param frac Fraction of width and height treated as a corner.
#' @return One of "topright", "bottomleft", "bottomright"; defaults to
#'   "topright" when the template cannot be read.
pick_key_position <- function(template_png, frac = 0.28) {
    if (!file.exists(template_png) || !requireNamespace("png", quietly = TRUE)) {
        return("topright")
    }
    img <- tryCatch(png::readPNG(template_png), error = function(e) NULL)
    if (is.null(img)) return("topright")

    grey <- if (length(dim(img)) == 3) {
        apply(img[, , seq_len(min(3, dim(img)[3])), drop = FALSE], c(1, 2), mean)
    } else img
    nr <- nrow(grey); nc <- ncol(grey)
    rr <- max(1L, floor(nr * frac)); cc <- max(1L, floor(nc * frac))
    ink <- function(rows, cols) mean(grey[rows, cols] < 0.86)

    # topleft is deliberately not a candidate: every KEGG map puts its title box
    # there, and it is thin enough to score as empty -- on ko04510 the key landed
    # squarely on "FOCAL ADHESION".
    scores <- c(
        topright    = ink(seq_len(rr),               seq.int(nc - cc + 1L, nc)),
        bottomleft  = ink(seq.int(nr - rr + 1L, nr), seq_len(cc)),
        bottomright = ink(seq.int(nr - rr + 1L, nr), seq.int(nc - cc + 1L, nc))
    )
    names(scores)[which.min(scores)]
}


#' Thresholds behind the multi-omics KEGG (pathview) maps
#'
#' One place to keep the three numbers, so the pathway filter, the node filter
#' and the figure legend cannot drift apart.
#'
#' \itemize{
#'   \item `fdr_alpha` -- adjusted ORA p below which a pathway is drawn.
#'   \item `node_fc` -- linear fold change (either direction) above which a
#'     feature may colour its node.
#'   \item `node_p` -- raw DE p below which a feature may colour its node,
#'     regardless of effect size.
#' }
PATHVIEW_SIGNIFICANCE <- list(fdr_alpha = 0.05, node_fc = 1.5, node_p = 0.05)


#' Describe, for a figure legend, how the pathview maps were thresholded
#'
#' The report template picks this up so the maps state their own selection
#' rules instead of leaving the reader to guess why a node is white.
#'
#' @param thresholds List as in \code{PATHVIEW_SIGNIFICANCE}, with elements
#'   `fdr_alpha`, `node_fc` and `node_p`.
#' @return A single character string, ready to pass to `figure_legend()`.
pathview_significance_caption <- function(thresholds = PATHVIEW_SIGNIFICANCE) {
    paste0(
        "Pathways were selected at FDR < ", thresholds$fdr_alpha,
        " (adjusted over-representation p) in at least one omics layer; where a ",
        "layer's ORA table carries no adjusted p, raw p < ", thresholds$fdr_alpha,
        " was used for that layer and the substitution is recorded in the run log. ",
        "A node is coloured only where the feature mapping to it both changed by more than ",
        thresholds$node_fc, "-fold in either direction (|log2 fold change| > ",
        round(log2(thresholds$node_fc), 2), ") and reached raw p < ", thresholds$node_p,
        ". A node left at the map's own colour either carries no measured feature ",
        "in that layer or carries only features that failed one or both cut-offs -- ",
        "the map does not separate the two. Where several features map to one node, the ",
        "node carries the mean log2 fold change of those that passed."
    )
}


#' Select KEGG pathways from per-omic ORA tables on adjusted p
#'
#' A pathway qualifies if any single layer calls it enriched, so the union is
#' taken across files. A layer with no usable adjusted p -- the column is absent,
#' or empty over this ID space -- is judged on the raw p instead. That is a laxer
#' scale, so the substitution is announced rather than folded silently into the
#' union.
#'
#' @param ora_files Character vector of per-omic ORA CSV paths.
#' @param id_pattern Regex the `pathway` column must match (organism or KO space).
#' @param alpha Adjusted-p cut-off.
#' @return Character vector of unique pathway ids, possibly empty.
select_ora_pathways_by_fdr <- function(ora_files, id_pattern,
                                       alpha = PATHVIEW_SIGNIFICANCE$fdr_alpha) {
    selected <- character(0)
    for (f in ora_files) {
        d <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
        if (is.null(d) || !"pathway" %in% names(d)) next
        in_space <- grepl(id_pattern, d$pathway)
        if (!any(in_space)) next

        pcol <- intersect(c("padj", "adj.P.Val"), names(d))
        # An adjusted-p column that is empty over this ID space is no adjusted p
        # at all: judging on it would drop the layer without a word.
        if (length(pcol) > 0 && all(is.na(d[[pcol[1]]][in_space]))) {
            pcol <- character(0)
        }
        if (length(pcol) == 0) {
            pcol <- intersect(c("pvalue", "P.Value"), names(d))
            if (length(pcol) == 0) {
                message("    Union pathview: ", basename(f),
                        " has no p-value column; skipped")
                next
            }
            message("    Union pathview: ", basename(f), " has no usable adjusted p; ",
                    "selecting on raw ", pcol[1], " < ", alpha, " for this layer")
        }
        pcol <- pcol[1]
        keep <- in_space & !is.na(d[[pcol]]) & d[[pcol]] < alpha
        selected <- c(selected, d$pathway[keep])
    }
    unique(selected)
}


#' Keep only the features allowed to colour a pathview node
#'
#' Colouring every measured feature makes a KEGG map look like a heat map of
#' noise: most nodes end up tinted by fold changes indistinguishable from zero.
#' Features that neither moved by `node_fc` nor reached `node_p` are dropped, so
#' pathview leaves their nodes at the template colour rather than painting them
#' the mid-point of the scale.
#'
#' @param de_table Standardized DE table with `log2fc` and a p-value column
#'   (`pvalue`, `P.Value`, `padj` or `adj.P.Val`), as returned by
#'   \code{extract_de_tables}. `log2fc` must already be on the log2 scale.
#' @param thresholds List as in \code{PATHVIEW_SIGNIFICANCE}.
#' @param label Layer name used in the console message.
#' @return The subset of rows eligible to colour a node; zero rows if the table
#'   carries no usable log2 fold change.
filter_changed_features <- function(de_table, thresholds = PATHVIEW_SIGNIFICANCE,
                                    label = "features") {
    if (is.null(de_table) || !is.data.frame(de_table) || nrow(de_table) == 0) {
        return(de_table)
    }
    if (!"log2fc" %in% names(de_table)) return(de_table[0, , drop = FALSE])

    pcol <- intersect(c("pvalue", "P.Value", "padj", "adj.P.Val"), names(de_table))
    pcol <- if (length(pcol) > 0) pcol[1] else NA_character_

    usable <- is.finite(de_table$log2fc)
    moved <- usable & abs(de_table$log2fc) > log2(thresholds$node_fc)
    sig <- if (is.na(pcol)) {
        rep(FALSE, nrow(de_table))
    } else {
        !is.na(de_table[[pcol]]) & de_table[[pcol]] < thresholds$node_p
    }
    # Both conditions, not either: a fold change that large with a weak p is
    # usually a noisy low-abundance feature, and a small but confident change is
    # not what a pathway map is meant to highlight. Requiring both keeps the
    # coloured nodes to features that moved AND are supported.
    keep <- usable & moved & sig

    message("    Union pathview: ", label, " ", sum(keep), "/", nrow(de_table),
            " features eligible to colour a node (|log2FC| > ",
            round(log2(thresholds$node_fc), 2), " and ",
            if (is.na(pcol)) "no p column" else paste0(pcol, " < ", thresholds$node_p),
            ")")
    de_table[keep, , drop = FALSE]
}


#' Render KEGG maps for the union of pathways enriched in either gene omic
#'
#' Reads the per-omic ORA tables written by the RNA / proteomics enrichment
#' steps, takes the union of the KEGG pathways they call enriched at
#' FDR < 0.05, and renders one map per pathway with the RNA and protein log2FC
#' overlaid as two states. Only features that actually changed colour a node
#' (see \code{filter_changed_features}); the rest are left at the template
#' colour so a map shows its signal rather than its coverage.
#'
#' Runs in one of two ID spaces:
#' \itemize{
#'   \item organism space -- the species has a KEGG code, so native
#'     `<org>#####` pathways are drawn and feature ids double as KEGG gene ids
#'     (proteins are bridged through the gene-protein mapping).
#'   \item KO space -- the species has no KEGG code at all (non-model organism).
#'     The KEGG reference `map#####` maps are drawn instead, and features are
#'     translated to KEGG Orthology ids through `enrichment.pathview.ko_map`.
#'     Metabolite log2FC is passed as `cpd.data` here, so one map carries all
#'     three layers -- which is the only reason the KO detour is worth taking.
#' }
#'
#' @param de_results Named list of DE results per omics.
#' @param harmonization_res Harmonization result (supplies the gene-protein map
#'   and the metabolomics row data).
#' @param config Full config.
#' @param out_dir Multi-ORA output directory (maps go under out_dir/pathview).
#' @param top_n Max pathways to render (default 25; overridden by
#'   `modes$multiomics$enrichment$pathview$top_n`).
#' @return Path to the compiled PDF, or NULL.
generate_per_omic_union_pathview <- function(de_results, harmonization_res,
                                             config, out_dir, top_n = 25) {
    if (!requireNamespace("pathview", quietly = TRUE)) return(NULL)
    organism <- config$global$organism %||% ""
    kegg_org <- get_kegg_organism(organism)

    pv_cfg <- config$modes$multiomics$enrichment$pathview %||% list()
    top_n <- pv_cfg$top_n %||% top_n
    ko_mode <- identical(pv_cfg$species %||% "", "ko") ||
        is.null(kegg_org) || !nzchar(kegg_org)

    ko_map <- if (ko_mode) load_feature_ko_map(config) else NULL
    if (ko_mode && (is.null(ko_map) || nrow(ko_map) == 0)) {
        # No KEGG code for this organism and no KO map: nothing can be drawn.
        return(NULL)
    }
    pv_species <- if (ko_mode) "ko" else kegg_org
    kegg_re <- if (ko_mode) "^map[0-9]+$" else paste0("^", kegg_org, "[0-9]+$")

    # 1. Union of KEGG pathways enriched at FDR < 0.05 in RNA OR protein, read
    #    from the per-omic ORA tables written by the RNA/proteomics enrichment
    #    steps. Locate them by walking up from out_dir to the run root.
    run_root <- out_dir
    for (i in seq_len(8)) {
        if (dir.exists(file.path(run_root, "rna", "Enrichment")) ||
            dir.exists(file.path(run_root, "proteomics", "Enrichment"))) break
        parent <- dirname(run_root); if (identical(parent, run_root)) break; run_root <- parent
    }
    # The ORA filename prefix carries the contrast and gene-set database, both
    # project-specific, so only the "_ora_up/down.csv" tail is fixed.
    ora_files <- c(list.files(file.path(run_root, "rna", "Enrichment"),
                              "_ora_(up|down)\\.csv$", full.names = TRUE),
                   list.files(file.path(run_root, "proteomics", "Enrichment"),
                              "_ora_(up|down)\\.csv$", full.names = TRUE))
    pathways <- select_ora_pathways_by_fdr(ora_files, kegg_re,
                                           alpha = PATHVIEW_SIGNIFICANCE$fdr_alpha)
    if (length(pathways) == 0) {
        message("  Union pathview: no ", pv_species, " KEGG pathways at FDR < ",
                PATHVIEW_SIGNIFICANCE$fdr_alpha, ".")
        return(NULL)
    }
    message("  Union pathview: ", length(pathways), " pathway(s) at FDR < ",
            PATHVIEW_SIGNIFICANCE$fdr_alpha, " in at least one layer")
    # Same reference-map caveat as the cross-omics tables: drop the KEGG classes
    # this organism cannot have before spending a render on them.
    excl <- config$modes$multiomics$enrichment$exclude_pathway_classes
    if (!is.null(excl) && length(excl) > 0) {
        pathways <- pathways[keep_kegg_pathways(pathways, exclude = unlist(excl),
                                                cache_dir = out_dir,
                                                label = "pathview maps")]
    }
    if (length(pathways) == 0) {
        message("  Union pathview: no pathways left after KEGG class filtering.")
        return(NULL)
    }
    pathways <- head(pathways, top_n)

    # 2. Gene data in the ID space pathview expects: native KEGG gene ids in
    #    organism mode, KO ids in KO mode.
    gpm <- harmonization_res$gene_protein_mapping
    layer_fc <- function(om) {
        tbls <- tryCatch(extract_de_tables(de_results[[om]], om, harmonization_res),
                         error = function(e) NULL)
        if (is.null(tbls) || length(tbls) == 0) return(NULL)
        df <- tbls[[1]]
        if (!all(c("feature_id", "log2fc") %in% names(df))) return(NULL)
        # Filter before aggregating: an unchanged paralogue must not drag the
        # node's mean towards zero, it must simply not be part of it.
        df <- filter_changed_features(df, label = om)
        if (nrow(df) == 0) return(NULL)
        if (ko_mode) return(aggregate_log2fc_by_ko(df, ko_map, om))
        ids <- df$feature_id
        if (identical(om, "proteomics") && !is.null(gpm)) {
            ids <- gpm$gene_id[match(df$feature_id, gpm$protein_id)]
        }
        ok <- !is.na(ids) & is.finite(df$log2fc)
        if (!any(ok)) return(NULL)
        # Many features can collapse onto one node (paralogues, protein groups),
        # so average rather than let an arbitrary one win.
        tapply(df$log2fc[ok], ids[ok], mean, na.rm = TRUE)
    }
    rna_fc  <- layer_fc("transcriptomics")
    prot_fc <- layer_fc("proteomics")
    if (is.null(rna_fc) && is.null(prot_fc)) return(NULL)
    genes <- unique(c(names(rna_fc), names(prot_fc)))
    gene_data <- matrix(NA_real_, length(genes), 2,
                        dimnames = list(genes, c("RNA", "Protein")))
    if (!is.null(rna_fc))  gene_data[names(rna_fc), 1]  <- rna_fc
    if (!is.null(prot_fc)) gene_data[names(prot_fc), 2] <- prot_fc
    message("  Union pathview: ", length(genes), " ",
            if (ko_mode) "KO" else "gene", " nodes coloured (changed features only)")

    # 3. Compound data, so metabolites colour the compound nodes. Only in KO
    #    mode: organism mode keeps the two-layer view the caller already had.
    cpd_data <- NULL
    if (ko_mode && "metabolomics" %in% names(de_results)) {
        cpd_data <- tryCatch({
            metab_tbls <- extract_de_tables(de_results$metabolomics, "metabolomics",
                                            harmonization_res)
            metab_map <- map_metabolite_ids_to_kegg(metab_tbls, harmonization_res)
            if (is.null(metab_map) || nrow(metab_map) == 0 || length(metab_tbls) == 0) {
                NULL
            } else {
                # Same rule as the gene layers, so a compound node is never
                # coloured on a difference the DA analysis did not support.
                mt <- filter_changed_features(metab_tbls[[1]], label = "metabolomics")
                if (nrow(mt) == 0) {
                    NULL
                } else {
                    md <- merge(mt, metab_map, by = "feature_id")
                    fc <- tapply(md$log2fc, md$KEGG_CPD, mean, na.rm = TRUE)
                    stats::setNames(as.numeric(fc), names(fc))
                }
            }
        }, error = function(e) NULL)
        if (!is.null(cpd_data)) {
            message("  Union pathview: ", length(cpd_data),
                    " compound nodes coloured (changed features only)")
        }
    }

    # 4. Render one map per pathway (both gene layers overlaid).
    # pathview needs its 'bods' dataset in the global env when called via ::.
    if (!exists("bods", envir = globalenv())) {
        utils::data("bods", package = "pathview", envir = globalenv())
    }
    pv_dir <- file.path(out_dir, "pathview")
    dir.create(pv_dir, recursive = TRUE, showWarnings = FALSE)

    # Clear overlays from a previous run before writing this one. The report
    # globs this directory, so maps that a narrower pathway selection no longer
    # produces would otherwise linger and be reported alongside the current set.
    # The downloaded KEGG templates (ko#####.png / .xml) are the pathview cache
    # and are deliberately kept.
    stale <- list.files(pv_dir, pattern = "\\.multi_ora[^/]*\\.png$", full.names = TRUE)
    if (length(stale) > 0) {
        unlink(stale)
        message("  Union pathview: cleared ", length(stale), " map(s) from a previous run")
    }

    cwd <- getwd(); setwd(pv_dir); on.exit(setwd(cwd), add = TRUE)
    generated <- character(0)
    for (pid in pathways) {
        clean_pid <- normalize_kegg_pathway_id(pid)
        # Place the colour key in whichever corner the map itself leaves empty.
        # The template is cached by an earlier pathview call, so this is free on
        # a re-run and falls back to the default corner on a first run.
        key_pos <- pick_key_position(
            file.path(pv_dir, paste0(pv_species, clean_pid, ".png")))
        tryCatch({
            pathview::pathview(gene.data = gene_data, cpd.data = cpd_data,
                               pathway.id = clean_pid,
                               species = pv_species, gene.idtype = "KEGG",
                               out.suffix = "multi_ora", kegg.dir = pv_dir,
                               key.pos = key_pos,
                               multi.state = TRUE, same.layer = FALSE)
            f <- c(paste0(pv_species, clean_pid, ".multi_ora.multi.png"),
                   paste0(pv_species, clean_pid, ".multi_ora.png"))
            f <- f[file.exists(f)]
            if (length(f) > 0) {
                generated <- c(generated, file.path(pv_dir, f[1]))
                message("    Union pathview: ", pid)
            }
        }, error = function(e) message("    Union pathview failed for ", pid, ": ", e$message))
    }
    if (length(generated) == 0) return(NULL)

    # 5. Compile into the PDF the report uses as the has_pathview gate.
    pdf_path <- file.path(out_dir, "multi_ora_pathview_supported.pdf")
    tryCatch({
        grDevices::pdf(pdf_path, width = 12, height = 8)
        for (png_file in generated) {
            img <- png::readPNG(png_file)
            grid::grid.newpage()
            grid::grid.raster(img, y = 0.48, height = 0.92)
        }
        grDevices::dev.off()
    }, error = function(e) { tryCatch(grDevices::dev.off(), error = function(e2) NULL) })
    message("  Union pathview: ", length(generated), " maps -> ", basename(pdf_path))
    pdf_path
}


#' Generate pathview plots for top pathways from each omics, overlaying both layers
#'
#' Two sets of plots:
#' 1. Top N metabolomics-enriched pathways with proteomics enzyme logFC overlaid
#' 2. Top N proteomics-enriched pathways with metabolomics compound logFC overlaid
#'
#' @param per_omics_ora Named list of per-omics ORA results
#' @param metab_ora Metabolomics compound ORA results
#' @param de_results Named list of DE results per omics
#' @param harmonization_res Harmonization result
#' @param config Full config
#' @param out_dir Output directory
#' @param top_n Number of top pathways per omics (default 5)
#' @return List with paths to compiled PDFs
generate_per_omics_pathview <- function(per_omics_ora, metab_ora, de_results,
                                        harmonization_res, config, out_dir,
                                        top_n = 5) {

    if (!requireNamespace("pathview", quietly = TRUE)) {
        message("  Package 'pathview' not installed. Skipping per-omics pathview.")
        return(NULL)
    }

    organism <- config$global$organism %||% "human"
    kegg_org <- get_kegg_organism(organism)
    org_db <- get_organism_db(organism)
    if (is.null(kegg_org)) return(NULL)

    pv_cfg <- config$modes$multiomics$enrichment$pathview %||% list()
    top_n <- pv_cfg$top_n %||% top_n

    pv_dir <- file.path(out_dir, "pathview")
    dir.create(pv_dir, recursive = TRUE, showWarnings = FALSE)

    # --- Build proteomics gene logFC vector ---
    gene_data <- NULL
    if ("proteomics" %in% names(de_results)) {
        de_tables <- extract_de_tables(de_results$proteomics, "proteomics",
                                        harmonization_res)
        if (!is.null(de_tables) && length(de_tables) > 0) {
            id_map <- tryCatch(
                map_feature_ids_to_entrez(de_tables, "proteomics",
                                           harmonization_res, org_db),
                error = function(e) NULL
            )
            if (!is.null(id_map) && nrow(id_map) > 0) {
                df <- de_tables[[1]]
                df_mapped <- merge(df, id_map, by = "feature_id")
                fc_arr <- tapply(df_mapped$log2fc, df_mapped$ENTREZID,
                                  mean, na.rm = TRUE)
                gene_data <- as.numeric(fc_arr)
                names(gene_data) <- names(fc_arr)
                message("  Pathview: ", length(gene_data),
                        " proteomics enzymes with logFC")
            }
        }
    }

    # Also try transcriptomics if available
    if (is.null(gene_data) && "transcriptomics" %in% names(de_results)) {
        de_tables <- extract_de_tables(de_results$transcriptomics, "transcriptomics",
                                        harmonization_res)
        if (!is.null(de_tables) && length(de_tables) > 0) {
            id_map <- tryCatch(
                map_feature_ids_to_entrez(de_tables, "transcriptomics",
                                           harmonization_res, org_db),
                error = function(e) NULL
            )
            if (!is.null(id_map) && nrow(id_map) > 0) {
                df <- de_tables[[1]]
                df_mapped <- merge(df, id_map, by = "feature_id")
                fc_arr <- tapply(df_mapped$log2fc, df_mapped$ENTREZID,
                                  mean, na.rm = TRUE)
                gene_data <- as.numeric(fc_arr)
                names(gene_data) <- names(fc_arr)
                message("  Pathview: ", length(gene_data),
                        " transcriptomics genes with logFC")
            }
        }
    }

    # --- Build metabolomics compound logFC vector ---
    cpd_data <- NULL
    if ("metabolomics" %in% names(de_results)) {
        de_tables <- extract_de_tables(de_results$metabolomics, "metabolomics",
                                        harmonization_res)
        id_map <- tryCatch(
            map_metabolite_ids_to_kegg(de_tables, harmonization_res),
            error = function(e) NULL
        )
        if (!is.null(id_map) && nrow(id_map) > 0 && length(de_tables) > 0) {
            df <- de_tables[[1]]
            df_mapped <- merge(df, id_map, by = "feature_id")
            cpd_fc <- tapply(df_mapped$log2fc, df_mapped$KEGG_CPD,
                              mean, na.rm = TRUE)
            cpd_data <- as.numeric(cpd_fc)
            names(cpd_data) <- names(cpd_fc)
            message("  Pathview: ", length(cpd_data),
                    " metabolites with logFC")
        }
    }

    if (is.null(gene_data) && is.null(cpd_data)) {
        message("  No logFC data available for per-omics pathview")
        return(NULL)
    }

    # Ensure pathview bods dataset is loaded
    if (!exists("bods", envir = globalenv())) {
        utils::data("bods", package = "pathview", envir = globalenv())
    }

    cwd <- getwd()
    setwd(pv_dir)
    on.exit(setwd(cwd), add = TRUE)

    results <- list()

    # ------------------------------------------------------------------
    # Set 1: Top metabolomics pathways with proteomics enzyme overlay
    # ------------------------------------------------------------------
    if (!is.null(metab_ora) && nrow(metab_ora) > 0) {
        # Sort by pvalue, take top_n
        met_top <- metab_ora[order(metab_ora$pvalue), ]
        met_top <- head(met_top, top_n)

        message("  Pathview set 1: top ", nrow(met_top),
                " metabolomics pathways + proteomics overlay")

        met_pngs <- character(0)
        for (i in seq_len(nrow(met_top))) {
            pid <- met_top$ID[i]
            pw_name <- met_top$pathway[i]
            clean_pid <- normalize_kegg_pathway_id(pid)

            tryCatch({
                pathview::pathview(
                    gene.data  = gene_data,
                    cpd.data   = cpd_data,
                    pathway.id = clean_pid,
                    species    = kegg_org,
                    out.suffix = "metab_top",
                    kegg.dir   = pv_dir,
                    keys.align = "y",
                    match.data = TRUE,
                    multi.state = FALSE,
                    same.layer  = TRUE
                )

                candidates <- c(
                    paste0(kegg_org, clean_pid, ".metab_top.png"),
                    paste0(kegg_org, clean_pid, ".metab_top.multi.png")
                )
                found <- candidates[file.exists(candidates)]
                if (length(found) > 0) {
                    met_pngs <- c(met_pngs, file.path(pv_dir, found[1]))
                    message("    ", pw_name, " (", pid, ")")
                }
            }, error = function(e) {
                message("    Pathview failed for ", pid, ": ", e$message)
            })
        }

        if (length(met_pngs) > 0) {
            pdf_path <- file.path(out_dir,
                                   "pathview_top_metabolomics_pathways.pdf")
            tryCatch({
                grDevices::pdf(pdf_path, width = 12, height = 8)
                for (png_file in met_pngs) {
                    img <- png::readPNG(png_file)
                    grid::grid.newpage()
                    grid::grid.raster(img)
                }
                grDevices::dev.off()
                message("  Compiled ", length(met_pngs),
                        " metabolomics pathview maps into: ",
                        basename(pdf_path))
                results$metabolomics_pdf <- pdf_path
            }, error = function(e) {
                message("  PDF compilation failed: ", e$message)
                tryCatch(grDevices::dev.off(), error = function(e2) NULL)
            })
        }
    }

    # ------------------------------------------------------------------
    # Set 2: Top proteomics pathways with metabolomics compound overlay
    # ------------------------------------------------------------------
    prot_ora <- per_omics_ora[["proteomics"]]
    if (is.null(prot_ora) || nrow(prot_ora) == 0) {
        prot_ora <- per_omics_ora[["transcriptomics"]]
    }

    if (!is.null(prot_ora) && nrow(prot_ora) > 0) {
        prot_top <- prot_ora[order(prot_ora$pvalue), ]
        prot_top <- head(prot_top, top_n)

        message("  Pathview set 2: top ", nrow(prot_top),
                " proteomics pathways + metabolomics overlay")

        prot_pngs <- character(0)
        for (i in seq_len(nrow(prot_top))) {
            pid <- prot_top$ID[i]
            pw_name <- prot_top$pathway[i]
            clean_pid <- normalize_kegg_pathway_id(pid)

            tryCatch({
                pathview::pathview(
                    gene.data  = gene_data,
                    cpd.data   = cpd_data,
                    pathway.id = clean_pid,
                    species    = kegg_org,
                    out.suffix = "prot_top",
                    kegg.dir   = pv_dir,
                    keys.align = "y",
                    match.data = TRUE,
                    multi.state = FALSE,
                    same.layer  = TRUE
                )

                candidates <- c(
                    paste0(kegg_org, clean_pid, ".prot_top.png"),
                    paste0(kegg_org, clean_pid, ".prot_top.multi.png")
                )
                found <- candidates[file.exists(candidates)]
                if (length(found) > 0) {
                    prot_pngs <- c(prot_pngs, file.path(pv_dir, found[1]))
                    message("    ", pw_name, " (", pid, ")")
                }
            }, error = function(e) {
                message("    Pathview failed for ", pid, ": ", e$message)
            })
        }

        if (length(prot_pngs) > 0) {
            pdf_path <- file.path(out_dir,
                                   "pathview_top_proteomics_pathways.pdf")
            tryCatch({
                grDevices::pdf(pdf_path, width = 12, height = 8)
                for (png_file in prot_pngs) {
                    img <- png::readPNG(png_file)
                    grid::grid.newpage()
                    grid::grid.raster(img)
                }
                grDevices::dev.off()
                message("  Compiled ", length(prot_pngs),
                        " proteomics pathview maps into: ",
                        basename(pdf_path))
                results$proteomics_pdf <- pdf_path
            }, error = function(e) {
                message("  PDF compilation failed: ", e$message)
                tryCatch(grDevices::dev.off(), error = function(e2) NULL)
            })
        }
    }

    if (length(results) == 0) {
        message("  No per-omics pathview plots generated")
        return(NULL)
    }

    results
}
