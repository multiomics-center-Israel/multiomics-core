#' Run DIABLO (mixOmics) integration analysis
#'
#' DIABLO (Data Integration Analysis for Biomarker discovery using Latent cOmponents)
#' is a supervised multi-omics integration method from the mixOmics package.
#'
#' @param mae MultiAssayExperiment object with aligned samples
#' @param config Full config object
#' @param out_dir Output directory for plots
#' @return List with: model, loadings, sample_scores, plots, performance
run_diablo_integration <- function(mae, config, out_dir = NULL) {

    if (!requireNamespace("mixOmics", quietly = TRUE)) {
        stop(
            "Package 'mixOmics' is required for DIABLO integration. ",
            "Install with: BiocManager::install('mixOmics')",
            call. = FALSE
        )
    }

    cfg <- config$modes$multiomics$integration$diablo %||% list()
    ncomp <- cfg$ncomp %||% 2
    design_type <- cfg$design_matrix %||% "full"  # "full", "null", or numeric matrix
    cv_folds <- cfg$cv_folds %||% 5

    message("Running DIABLO integration...")

    # Extract data from MAE
    omics <- names(mae@ExperimentList)
    data_list <- lapply(omics, function(om) {
        t(SummarizedExperiment::assay(mae[[om]], "expr"))  # Samples x Features
    })
    names(data_list) <- omics

    # Extract outcome variable (Y)
    condition_col <- config$modes$multiomics$condition_column %||%
                     config$design$condition_column %||%
                     "condition"

    coldata <- SummarizedExperiment::colData(mae)
    if (!condition_col %in% colnames(coldata)) {
        stop("Condition column '", condition_col, "' not found in MAE colData")
    }

    Y <- as.factor(coldata[[condition_col]])

    # Build design matrix (defines correlation structure between omics)
    if (is.character(design_type)) {
        if (design_type == "full") {
            # Full correlation between all omics
            design <- matrix(1, length(omics), length(omics),
                             dimnames = list(omics, omics))
            diag(design) <- 0
        } else if (design_type == "null") {
            # No correlation assumed
            design <- matrix(0, length(omics), length(omics),
                             dimnames = list(omics, omics))
        } else {
            stop("Unknown design_matrix type: ", design_type)
        }
    } else {
        design <- design_type  # User-provided matrix
    }

    # Initial DIABLO model (without tuning)
    diablo_model <- tryCatch({
        mixOmics::block.plsda(
            X = data_list,
            Y = Y,
            ncomp = ncomp,
            design = design
        )
    }, error = function(e) {
        stop("DIABLO model fitting failed: ", e$message, call. = FALSE)
    })

    # Extract results
    sample_scores <- diablo_model$variates
    loadings <- diablo_model$loadings

    # Compute performance (cross-validation)
    perf <- NULL
    if (cv_folds > 0 && nrow(data_list[[1]]) >= cv_folds * 2) {
        message(sprintf("  Running %d-fold cross-validation...", cv_folds))
        perf <- tryCatch({
            mixOmics::perf(diablo_model, validation = "Mfold", folds = cv_folds,
                           progressBar = FALSE, cpus = 1)
        }, error = function(e) {
            warning("Cross-validation failed: ", e$message)
            NULL
        })
    }

    # Generate plots
    plots <- list()
    if (!is.null(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

        # Sample plot (scores on first 2 components)
        plots$sample_plot <- file.path(out_dir, "diablo_sample_plot.png")
        png(plots$sample_plot, width = 800, height = 600, res = 120)
        tryCatch({
            mixOmics::plotIndiv(diablo_model, comp = c(1, 2), group = Y,
                                ind.names = FALSE, legend = TRUE,
                                title = "DIABLO Sample Scores")
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("Plot failed:", e$message), cex = 1.2)
        })
        dev.off()

        # Variable plot (loadings) — suppress feature labels to avoid clutter
        plots$variable_plot <- file.path(out_dir, "diablo_variable_plot.png")
        png(plots$variable_plot, width = 1000, height = 800, res = 120)
        tryCatch({
            mixOmics::plotVar(diablo_model, comp = c(1, 2), style = "graphics",
                              legend = TRUE, title = "DIABLO Variable Loadings",
                              var.names = FALSE)
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("Plot failed:", e$message), cex = 1.2)
        })
        dev.off()

        # Circos plot (correlations between omics) — suppress feature labels
        plots$circos_plot <- file.path(out_dir, "diablo_circos_plot.png")
        png(plots$circos_plot, width = 1000, height = 1000, res = 120)
        tryCatch({
            mixOmics::circosPlot(diablo_model, cutoff = 0.5, size.variables = 0.5,
                                 var.names = FALSE)
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("Plot failed:", e$message), cex = 1.2)
        })
        dev.off()

        message("  DIABLO plots saved to: ", out_dir)
    }

    # Build feature name map from MAE rowData
    feature_name_map <- build_feature_name_map(mae)

    # Interactive plotly loadings plot (pathway-colored)
    if (!is.null(out_dir)) {
        tryCatch({
            diablo_res_tmp <- list(
                loadings = loadings,
                feature_name_map = feature_name_map
            )
            plots$loadings_interactive <- plot_diablo_loadings_interactive(
                diablo_results = diablo_res_tmp,
                mae = mae,
                config = config,
                out_dir = out_dir,
                top_n = 50,
                color_by = "pathway"
            )
        }, error = function(e) {
            warning("Interactive DIABLO loadings plot failed: ", e$message)
        })
    }

    # Extract top features per component
    top_features <- extract_diablo_top_features(diablo_model, top_n = 50,
                                                 feature_name_map = feature_name_map)

    message(sprintf(
        "DIABLO integration complete: %d components, %d omics layers",
        ncomp, length(omics)
    ))

    list(
        model = diablo_model,
        sample_scores = sample_scores,
        loadings = loadings,
        top_features = top_features,
        feature_name_map = feature_name_map,
        performance = perf,
        plots = plots,
        design = design,
        config = cfg
    )
}


#' Extract top contributing features from DIABLO model
#'
#' @param diablo_model Fitted block.plsda object
#' @param top_n Number of top features per component per omics
#' @param feature_name_map Named vector mapping harmonized IDs to original names
#' @return Named list of data frames (one per omics)
extract_diablo_top_features <- function(diablo_model, top_n = 50,
                                         feature_name_map = NULL) {

    loadings <- diablo_model$loadings
    omics <- names(loadings)

    top_features <- list()

    for (om in omics) {
        load_mat <- loadings[[om]]

        # Get top features per component
        top_list <- list()
        for (comp_idx in seq_len(ncol(load_mat))) {
            comp_name <- colnames(load_mat)[comp_idx]
            comp_loadings <- load_mat[, comp_idx]

            # Rank by absolute loading
            ranked_idx <- order(abs(comp_loadings), decreasing = TRUE)
            n_select <- min(top_n, length(ranked_idx))

            feat_ids <- rownames(load_mat)[ranked_idx[seq_len(n_select)]]
            df <- data.frame(
                feature = feat_ids,
                loading = comp_loadings[ranked_idx[seq_len(n_select)]],
                abs_loading = abs(comp_loadings[ranked_idx[seq_len(n_select)]]),
                component = comp_name,
                stringsAsFactors = FALSE
            )

            # Add original names if available
            if (!is.null(feature_name_map)) {
                df$original_name <- unname(feature_name_map[feat_ids])
            }

            top_list[[comp_name]] <- df
        }

        top_features[[om]] <- do.call(rbind, top_list)
    }

    top_features
}


#' Write DIABLO results to CSV
#'
#' @param diablo_results Output from run_diablo_integration()
#' @param out_dir Output directory
write_diablo_results <- function(diablo_results, out_dir) {

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Sample scores
    for (om in names(diablo_results$sample_scores)) {
        scores <- diablo_results$sample_scores[[om]]
        write.csv(scores, file.path(out_dir, paste0("diablo_scores_", om, ".csv")),
                  row.names = TRUE)
    }

    # Top features
    for (om in names(diablo_results$top_features)) {
        feat_df <- diablo_results$top_features[[om]]
        write.csv(feat_df, file.path(out_dir, paste0("diablo_top_features_", om, ".csv")),
                  row.names = FALSE)
    }

    # Design matrix
    write.csv(diablo_results$design, file.path(out_dir, "diablo_design_matrix.csv"),
              row.names = TRUE)

    message("DIABLO results written to: ", out_dir)

    invisible(NULL)
}


#' Generate multi-omics feature heatmap from DIABLO results
#'
#' Reads DIABLO top_features CSVs and extracts expression data from the MAE
#' to create a cross-omics heatmap of the most important features.
#'
#' @param mae MultiAssayExperiment object (feature-selected subset)
#' @param diablo_dir Directory containing diablo_top_features_*.csv files
#' @param config Full pipeline config
#' @param top_n Number of top features per omics to include
#' @return Path to the saved PNG, or NULL on failure
plot_diablo_feature_heatmap <- function(mae, diablo_dir, config, top_n = 15) {

    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        warning("ComplexHeatmap not available — skipping DIABLO feature heatmap")
        return(NULL)
    }

    # Read top_features CSVs (skip _Y.csv — that's the outcome)
    csv_files <- list.files(diablo_dir, "^diablo_top_features_.*\\.csv$",
                            full.names = TRUE)
    csv_files <- csv_files[!grepl("_Y\\.csv$", csv_files)]

    if (length(csv_files) == 0) {
        message("No DIABLO top_features CSVs found — skipping heatmap")
        return(NULL)
    }

    # Extract condition column name
    condition_col <- config$modes$multiomics$condition_column %||%
                     config$design$condition_column %||%
                     "condition"

    # Collect top features per omics and build expression matrix
    all_expr <- list()
    omics_labels <- character(0)

    for (csv_f in csv_files) {
        omics_name <- sub("^diablo_top_features_(.*)\\.csv$", "\\1", basename(csv_f))

        feat_df <- read.csv(csv_f, stringsAsFactors = FALSE)

        comp1 <- feat_df[feat_df$component == "comp1", ]
        comp1 <- comp1[order(comp1$abs_loading, decreasing = TRUE), ]
        n_take <- min(top_n, nrow(comp1))
        if (n_take == 0) next
        top_feat <- comp1$feature[seq_len(n_take)]

        if (!omics_name %in% names(mae)) {
            message("  Omics '", omics_name, "' not in MAE — skipping")
            next
        }

        expr_mat <- SummarizedExperiment::assay(mae[[omics_name]], "expr")
        avail <- intersect(top_feat, rownames(expr_mat))
        if (length(avail) == 0) next

        sub_mat <- expr_mat[avail, , drop = FALSE]
        # Resolve to original names via rowData
        rd <- SummarizedExperiment::rowData(mae[[omics_name]])
        display <- if ("gene_symbol" %in% colnames(rd)) {
            syms <- as.character(rd[avail, "gene_symbol"])
            ifelse(is.na(syms) | syms == "", avail, syms)
        } else if ("original_id" %in% colnames(rd)) {
            as.character(rd[avail, "original_id"])
        } else if ("Name" %in% colnames(rd)) {
            nms <- as.character(rd[avail, "Name"])
            ifelse(is.na(nms) | nms == "", avail, nms)
        } else if ("Metabolite" %in% colnames(rd)) {
            nms <- as.character(rd[avail, "Metabolite"])
            ifelse(is.na(nms) | nms == "", avail, nms)
        } else if ("HMDB" %in% colnames(rd)) {
            nms <- as.character(rd[avail, "HMDB"])
            ifelse(is.na(nms) | nms == "", avail, nms)
        } else {
            avail
        }
        rownames(sub_mat) <- paste0("[", omics_name, "] ", display)
        all_expr[[omics_name]] <- sub_mat
        omics_labels <- c(omics_labels, rep(omics_name, nrow(sub_mat)))
    }

    if (length(all_expr) == 0) {
        message("No features matched MAE — skipping DIABLO heatmap")
        return(NULL)
    }

    combined <- do.call(rbind, all_expr)

    # Z-score per feature (row-wise)
    combined_z <- t(apply(combined, 1, function(x) {
        s <- sd(x, na.rm = TRUE)
        if (is.na(s) || s == 0) return(rep(0, length(x)))
        (x - mean(x, na.rm = TRUE)) / s
    }))
    colnames(combined_z) <- colnames(combined)

    # Annotations
    coldata <- as.data.frame(SummarizedExperiment::colData(mae))
    conditions <- coldata[[condition_col]]

    # Shorten long feature names for display
    display_names <- rownames(combined_z)
    display_names <- ifelse(nchar(display_names) > 50,
                            paste0(substr(display_names, 1, 47), "..."),
                            display_names)
    rownames(combined_z) <- display_names

    # Color palette
    col_fun <- circlize::colorRamp2(
        c(-2, 0, 2), c("navy", "white", "firebrick3")
    )

    # Omics color annotation (left side)
    omics_colors <- c(
        transcriptomics = "#E41A1C", proteomics = "#377EB8",
        metabolomics = "#4DAF4A"
    )
    omics_colors <- omics_colors[names(omics_colors) %in% unique(omics_labels)]
    row_ha <- ComplexHeatmap::rowAnnotation(
        Omics = omics_labels,
        col = list(Omics = omics_colors),
        show_annotation_name = FALSE
    )

    # Condition annotation (top)
    cond_levels <- unique(conditions)
    cond_colors <- setNames(
        scales::hue_pal()(length(cond_levels)),
        cond_levels
    )
    top_ha <- ComplexHeatmap::HeatmapAnnotation(
        Condition = conditions,
        col = list(Condition = cond_colors),
        show_annotation_name = FALSE
    )

    # Build heatmap — cluster all features together by distance,
    # omics type shown as row annotation color
    ht <- ComplexHeatmap::Heatmap(
        combined_z,
        name = "Z-score",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_column_names = FALSE,
        row_names_side = "right",
        row_names_gp = grid::gpar(fontsize = 7),
        left_annotation = row_ha,
        top_annotation = top_ha,
        column_title = "Top DIABLO Features Across Omics (Z-scored)",
        column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        heatmap_legend_param = list(title = "Z-score")
    )

    # Save heatmap
    n_features <- nrow(combined_z)
    plot_height <- max(8, 3 + n_features * 0.18)
    out_png <- file.path(diablo_dir, "diablo_feature_heatmap.png")
    png(out_png, width = 12, height = plot_height, units = "in", res = 150)
    tryCatch({
        ComplexHeatmap::draw(ht, merge_legend = TRUE)
    }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Heatmap failed:", e$message), cex = 1.2)
    })
    dev.off()

    message("DIABLO feature heatmap saved: ", out_png)
    out_png
}


#' Interactive plotly DIABLO loadings plot
#'
#' Creates a plotly scatter plot of DIABLO loadings (comp1 vs comp2) across
#' all omics layers. Points are shaped by omics type and can be colored by
#' KEGG pathway or GO term membership.
#'
#' @param diablo_results Output from run_diablo_integration()
#' @param mae MultiAssayExperiment used for ID resolution
#' @param config Full pipeline config
#' @param out_dir Output directory for the HTML widget
#' @param top_n Number of top features per omics to include (by abs loading)
#' @param color_by "pathway" (KEGG) or "GO" for coloring; NULL for omics-only
#' @param pathway_top_n Number of top pathways to show (rest grouped as "Other")
#' @return Path to saved HTML file, or NULL on failure
plot_diablo_loadings_interactive <- function(diablo_results, mae, config,
                                             out_dir, top_n = 50,
                                             color_by = "pathway",
                                             pathway_top_n = 8) {

    if (!requireNamespace("plotly", quietly = TRUE)) {
        warning("plotly not available — skipping interactive DIABLO loadings plot")
        return(NULL)
    }
    if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
        warning("htmlwidgets not available — skipping interactive DIABLO loadings plot")
        return(NULL)
    }

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    loadings <- diablo_results$loadings
    feature_name_map <- diablo_results$feature_name_map
    omics <- setdiff(names(loadings), "Y")

    # Build combined data frame: feature, comp1, comp2, omics, display_name
    df_list <- list()
    for (om in omics) {
        load_mat <- loadings[[om]]
        if (ncol(load_mat) < 2) next

        comp1 <- load_mat[, 1]
        comp2 <- load_mat[, 2]
        abs_max <- pmax(abs(comp1), abs(comp2))
        keep_idx <- order(abs_max, decreasing = TRUE)[seq_len(min(top_n, length(abs_max)))]

        feat_ids <- rownames(load_mat)[keep_idx]
        display <- if (!is.null(feature_name_map)) {
            nms <- feature_name_map[feat_ids]
            ifelse(is.na(nms) | nms == "", feat_ids, nms)
        } else {
            feat_ids
        }

        df_list[[om]] <- data.frame(
            feature_id = feat_ids,
            display_name = display,
            comp1 = comp1[keep_idx],
            comp2 = comp2[keep_idx],
            omics = om,
            stringsAsFactors = FALSE
        )
    }

    if (length(df_list) == 0) return(NULL)
    df <- do.call(rbind, df_list)
    rownames(df) <- NULL

    # --- Map features to pathways/GO for coloring ---
    df$pathway_label <- "Unmapped"

    if (!is.null(color_by) && color_by %in% c("pathway", "GO")) {
        organism <- config$global$organism %||% "human"
        org_db <- tryCatch(get_organism_db(organism), error = function(e) NULL)

        if (!is.null(org_db) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
            # Collect all gene-based feature IDs (skip metabolomics)
            gene_omics <- setdiff(omics, "metabolomics")
            gene_features <- df[df$omics %in% gene_omics, ]

            if (nrow(gene_features) > 0) {
                # Resolve harmonized IDs to ENTREZ via rowData
                entrez_map <- character(0)
                for (om in gene_omics) {
                    if (!om %in% names(mae)) next
                    rd <- SummarizedExperiment::rowData(mae[[om]])
                    if ("gene_id" %in% colnames(rd)) {
                        om_feats <- df$feature_id[df$omics == om]
                        avail <- intersect(om_feats, rownames(rd))
                        ids <- as.character(rd[avail, "gene_id"])
                        names(ids) <- avail
                        entrez_map <- c(entrez_map, ids)
                    }
                }

                if (length(entrez_map) > 0) {
                    valid <- entrez_map[!is.na(entrez_map) & entrez_map != ""]

                    if (color_by == "pathway") {
                        # KEGG pathway mapping
                        kegg_map <- tryCatch({
                            AnnotationDbi::select(
                                org_db,
                                keys = unique(valid),
                                keytype = "ENTREZID",
                                columns = "PATH"
                            )
                        }, error = function(e) NULL)

                        if (!is.null(kegg_map) && nrow(kegg_map) > 0) {
                            kegg_map <- kegg_map[!is.na(kegg_map$PATH), ]
                            # Map KEGG IDs to pathway names
                            kegg_org <- tryCatch(get_kegg_organism(organism),
                                                  error = function(e) NULL)
                            if (!is.null(kegg_org)) {
                                pathway_names <- tryCatch({
                                    plist <- clusterProfiler::download_KEGG(kegg_org, keggType = "KEGG")
                                    setNames(plist$KEGGPATHNAME2ID$to,
                                             plist$KEGGPATHNAME2ID$from)
                                }, error = function(e) NULL)

                                # Simpler: use KEGG REST
                                if (is.null(pathway_names)) {
                                    pathway_names <- tryCatch({
                                        kp <- KEGGREST::keggList("pathway", kegg_org)
                                        nm <- sub(paste0(" - .*$"), "", kp)
                                        setNames(nm, sub("path:", "", names(kp)))
                                    }, error = function(e) NULL)
                                }
                            }

                            # Build feature -> pathway name mapping
                            # Reverse: entrez_map is feature_id -> entrezid
                            entrez_to_feat <- split(names(valid), valid)
                            feat_pathways <- list()
                            for (i in seq_len(nrow(kegg_map))) {
                                eid <- kegg_map$ENTREZID[i]
                                pid <- kegg_map$PATH[i]
                                pname <- if (!is.null(pathway_names) && pid %in% names(pathway_names)) {
                                    pathway_names[pid]
                                } else {
                                    pid
                                }
                                fids <- entrez_to_feat[[eid]]
                                if (!is.null(fids)) {
                                    for (fid in fids) {
                                        feat_pathways[[fid]] <- c(feat_pathways[[fid]], pname)
                                    }
                                }
                            }

                            # Count pathway frequency to pick top pathways
                            all_pw <- unlist(feat_pathways)
                            pw_counts <- sort(table(all_pw), decreasing = TRUE)
                            top_pw <- names(pw_counts)[seq_len(min(pathway_top_n, length(pw_counts)))]

                            # Assign: first matching top pathway, else "Other pathway"
                            for (fid in names(feat_pathways)) {
                                pws <- feat_pathways[[fid]]
                                matched <- intersect(pws, top_pw)
                                if (length(matched) > 0) {
                                    df$pathway_label[df$feature_id == fid] <- matched[1]
                                } else {
                                    df$pathway_label[df$feature_id == fid] <- "Other pathway"
                                }
                            }
                        }
                    } else if (color_by == "GO") {
                        # GO biological process mapping
                        go_map <- tryCatch({
                            AnnotationDbi::select(
                                org_db,
                                keys = unique(valid),
                                keytype = "ENTREZID",
                                columns = "GOALL"
                            )
                        }, error = function(e) NULL)

                        if (!is.null(go_map) && nrow(go_map) > 0) {
                            go_bp <- go_map[!is.na(go_map$GOALL) &
                                            go_map$ONTOLOGYALL == "BP", ]
                            # Get GO term names
                            go_terms <- tryCatch({
                                AnnotationDbi::select(
                                    GO.db::GO.db,
                                    keys = unique(go_bp$GOALL),
                                    columns = "TERM"
                                )
                            }, error = function(e) NULL)

                            if (!is.null(go_terms)) {
                                go_name_map <- setNames(go_terms$TERM, go_terms$GOID)
                            } else {
                                go_name_map <- setNames(go_bp$GOALL, go_bp$GOALL)
                            }

                            entrez_to_feat <- split(names(valid), valid)
                            feat_go <- list()
                            for (i in seq_len(nrow(go_bp))) {
                                eid <- go_bp$ENTREZID[i]
                                gid <- go_bp$GOALL[i]
                                gname <- go_name_map[gid] %||% gid
                                fids <- entrez_to_feat[[eid]]
                                if (!is.null(fids)) {
                                    for (fid in fids) {
                                        feat_go[[fid]] <- c(feat_go[[fid]], gname)
                                    }
                                }
                            }

                            all_go <- unlist(feat_go)
                            go_counts <- sort(table(all_go), decreasing = TRUE)
                            top_go <- names(go_counts)[seq_len(min(pathway_top_n, length(go_counts)))]

                            for (fid in names(feat_go)) {
                                gos <- feat_go[[fid]]
                                matched <- intersect(gos, top_go)
                                if (length(matched) > 0) {
                                    df$pathway_label[df$feature_id == fid] <- matched[1]
                                } else {
                                    df$pathway_label[df$feature_id == fid] <- "Other GO term"
                                }
                            }
                        }
                    }
                }
            }

            # Metabolomics features stay as "Unmapped"
        }
    }

    # --- Build plotly figure ---
    omics_shapes <- c(
        transcriptomics = "circle",
        proteomics = "square",
        metabolomics = "diamond"
    )

    # Hover text
    df$hover <- paste0(
        "<b>", df$display_name, "</b><br>",
        "Omics: ", df$omics, "<br>",
        "Comp1: ", round(df$comp1, 4), "<br>",
        "Comp2: ", round(df$comp2, 4), "<br>",
        "Pathway/GO: ", df$pathway_label
    )

    # Assign shape per omics
    df$shape <- omics_shapes[df$omics]
    df$shape[is.na(df$shape)] <- "circle"

    fig <- plotly::plot_ly(
        data = df,
        x = ~comp1,
        y = ~comp2,
        color = ~pathway_label,
        symbol = ~omics,
        symbols = unname(omics_shapes[intersect(names(omics_shapes), unique(df$omics))]),
        text = ~hover,
        hoverinfo = "text",
        type = "scatter",
        mode = "markers",
        marker = list(size = 8, opacity = 0.8, line = list(width = 0.5, color = "grey30"))
    ) |>
        plotly::layout(
            title = list(
                text = "DIABLO Loadings (Interactive)",
                font = list(size = 16)
            ),
            xaxis = list(title = "Component 1 Loading", zeroline = TRUE,
                         zerolinecolor = "grey80"),
            yaxis = list(title = "Component 2 Loading", zeroline = TRUE,
                         zerolinecolor = "grey80"),
            legend = list(title = list(text = "Pathway / Omics")),
            hovermode = "closest"
        )

    # Save as HTML widget
    out_html <- file.path(out_dir, "diablo_loadings_interactive.html")
    htmlwidgets::saveWidget(
        plotly::as_widget(fig),
        file = normalizePath(out_html, mustWork = FALSE),
        selfcontained = TRUE,
        title = "DIABLO Loadings - Interactive"
    )

    message("  Interactive DIABLO loadings plot saved: ", out_html)
    out_html
}
