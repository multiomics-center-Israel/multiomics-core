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


#' Resolve feature IDs to human-readable display names
#'
#' Tries multiple columns in rowData to find display names. For metabolomics,
#' falls back to reading the original abundance file to recover HMDB IDs
#' when rowData annotation is incomplete.
#'
#' @param rd rowData data.frame for the omics layer
#' @param feat_ids Character vector of feature IDs to resolve
#' @param omics_name Name of the omics layer
#' @param config Full config object (for file paths)
#' @return Character vector of display names (same length as feat_ids)
resolve_feature_display_names <- function(rd, feat_ids, omics_name, config) {
    # Try standard annotation columns in priority order
    for (col in c("gene_symbol", "original_id", "Name", "Metabolite",
                   "HMDB", "HMDB_ID")) {
        if (col %in% colnames(rd)) {
            nms <- as.character(rd[feat_ids, col])
            n_valid <- sum(!is.na(nms) & nms != "")
            if (n_valid > length(feat_ids) * 0.3) {
                # At least 30% of features have names — use this column
                return(ifelse(is.na(nms) | nms == "", feat_ids, nms))
            }
        }
    }

    # Fallback for metabolomics: read HMDB_ID from the original abundance file
    if (grepl("metabol", omics_name, ignore.case = TRUE)) {
        hmdb_map <- recover_metabolomics_hmdb_ids(config)
        if (!is.null(hmdb_map)) {
            nms <- hmdb_map[feat_ids]
            if (sum(!is.na(nms)) > 0) {
                return(ifelse(is.na(nms) | nms == "", feat_ids, nms))
            }
        }
    }

    # Last resort: use feature IDs as-is
    feat_ids
}


#' Recover HMDB IDs from metabolomics abundance file
#'
#' When the MAE rowData lacks metabolite names, reads the original abundance
#' CSV to build a mapping from generic feature_N IDs to HMDB IDs.
#'
#' @param config Full config object
#' @return Named character vector: feature_N -> HMDB_ID, or NULL
recover_metabolomics_hmdb_ids <- function(config) {
    # Find the abundance file path
    metab_cfg <- config$modes$metabolomics
    if (is.null(metab_cfg)) return(NULL)

    data_file <- metab_cfg$files$data %||%
                 metab_cfg$files$preprocessed_counts %||%
                 NULL
    if (is.null(data_file) || data_file == "") return(NULL)

    # Resolve relative to project dir
    proj_dir <- config$project$dir %||% "."
    data_path <- if (file.exists(data_file)) {
        data_file
    } else {
        fp <- file.path(proj_dir, config$paths$raw %||% "data", data_file)
        if (file.exists(fp)) fp else NULL
    }
    if (is.null(data_path)) return(NULL)

    tryCatch({
        # Read just the first column to get HMDB IDs
        df <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)
        first_col <- colnames(df)[1]

        # Check if the first column contains HMDB-like IDs
        first_vals <- as.character(df[[first_col]])
        if (!any(grepl("^HMDB", first_vals, ignore.case = TRUE))) return(NULL)

        # Build mapping: feature_N -> HMDB_ID
        # The N corresponds to the row index in the original file
        feature_ids <- paste0("feature_", seq_len(nrow(df)))
        setNames(first_vals, feature_ids)
    }, error = function(e) NULL)
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
        rd <- as.data.frame(SummarizedExperiment::rowData(mae[[omics_name]]))
        display <- resolve_feature_display_names(rd, avail, omics_name, config)
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


#' Plot DIABLO loadings colored by log2FC from DE results
#'
#' Creates a scatter plot of comp1 vs comp2 loadings with points colored by
#' the log2 fold-change from differential expression results. Faceted by omics.
#'
#' @param diablo_results Output from run_diablo_integration()
#' @param de_results Named list of DE results per omics
#' @param config Full pipeline config
#' @param out_dir Output directory
#' @param top_n Number of top features per omics (by abs loading)
#' @param label_n Number of top features to label per omics facet
#' @return Vector of saved PNG paths (one per contrast)
plot_diablo_loadings_log2fc <- function(diablo_results, de_results, config,
                                        out_dir, top_n = 50, label_n = 5,
                                        mae = NULL) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)

    loadings <- diablo_results$loadings
    feature_name_map <- diablo_results$feature_name_map
    omics <- setdiff(names(loadings), "Y")

    if (length(omics) == 0 || is.null(de_results)) return(NULL)

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Build harmonized -> original ID mapping for DE matching and display names
    id_map <- build_harmonized_to_original_map(mae, de_results, config)

    # Build loading data frame across omics
    df_list <- list()
    for (om in omics) {
        load_mat <- loadings[[om]]
        if (ncol(load_mat) < 2) next

        comp1 <- load_mat[, 1]
        comp2 <- load_mat[, 2]
        abs_max <- pmax(abs(comp1), abs(comp2))
        keep_idx <- order(abs_max, decreasing = TRUE)[seq_len(min(top_n, length(abs_max)))]

        feat_ids <- rownames(load_mat)[keep_idx]

        # Resolve display names: feature_name_map first, then id_map fallback
        display <- if (!is.null(feature_name_map)) {
            nms <- feature_name_map[feat_ids]
            ifelse(is.na(nms) | nms == "" | nms == feat_ids, NA_character_, nms)
        } else {
            rep(NA_character_, length(feat_ids))
        }
        # Fill remaining NAs from id_map (covers metabolomics HMDB + proteomics IDs)
        still_na <- is.na(display)
        if (any(still_na) && !is.null(id_map)) {
            display[still_na] <- id_map[feat_ids[still_na]]
        }
        display[is.na(display)] <- feat_ids[is.na(display)]

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

    contrast_names <- get_contrast_names_from_de(de_results)
    saved_paths <- character(0)

    for (ci in seq_along(contrast_names)) {
        contrast <- contrast_names[ci]
        df_plot <- df
        df_plot$log2FC <- NA_real_

        for (om in omics) {
            if (is.null(de_results[[om]])) next
            de_norm <- tryCatch(
                extract_de_contrast(de_results[[om]], ci, om),
                error = function(e) NULL
            )
            if (is.null(de_norm)) next

            idx <- df_plot$omics == om
            feat_ids <- df_plot$feature_id[idx]
            # Try direct match first
            matched <- match(feat_ids, de_norm$feature_id)
            # If no matches, resolve harmonized IDs to original IDs via MAE
            if (all(is.na(matched)) && !is.null(id_map)) {
                orig_ids <- id_map[feat_ids]
                orig_ids[is.na(orig_ids)] <- feat_ids[is.na(orig_ids)]
                matched <- match(orig_ids, de_norm$feature_id)
            }
            df_plot$log2FC[idx] <- de_norm$logFC[matched]
        }

        df_sub <- df_plot[!is.na(df_plot$log2FC), ]
        if (nrow(df_sub) < 3) next

        # Clamp extreme values for symmetric color scale
        lfc_cap <- quantile(abs(df_sub$log2FC), 0.95, na.rm = TRUE)
        lfc_cap <- max(lfc_cap, 0.5)
        df_sub$log2FC_capped <- pmax(pmin(df_sub$log2FC, lfc_cap), -lfc_cap)

        # Select top features per facet for labeling
        df_sub$label <- ""
        for (om in unique(df_sub$omics)) {
            om_idx <- which(df_sub$omics == om)
            top_idx <- om_idx[order(pmax(abs(df_sub$comp1[om_idx]),
                                         abs(df_sub$comp2[om_idx])),
                                    decreasing = TRUE)]
            top_idx <- head(top_idx, label_n)
            df_sub$label[top_idx] <- df_sub$display_name[top_idx]
        }

        p <- ggplot2::ggplot(df_sub,
                ggplot2::aes(x = comp1, y = comp2, color = log2FC_capped)) +
            ggplot2::geom_point(size = 2.5, alpha = 0.7) +
            ggplot2::scale_color_gradient2(
                low = "#4DBBD5", mid = "white", high = "#E64B35",
                midpoint = 0, name = expression(log[2]~FC),
                limits = c(-lfc_cap, lfc_cap)
            ) +
            ggplot2::facet_wrap(~omics, scales = "free") +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
            ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.grid = ggplot2::element_blank(),
                           strip.text = ggplot2::element_text(face = "bold")) +
            ggplot2::labs(
                title = "DIABLO Loadings Colored by log2FC",
                subtitle = contrast,
                x = "Component 1 Loading",
                y = "Component 2 Loading"
            )

        if (requireNamespace("ggrepel", quietly = TRUE)) {
            p <- p + ggrepel::geom_text_repel(
                ggplot2::aes(label = label),
                size = 2.8, color = "grey20", max.overlaps = 15,
                segment.color = "grey50", segment.size = 0.3
            )
        }

        fname <- paste0("diablo_loadings_log2fc_", make.names(contrast), ".png")
        out_path <- file.path(out_dir, fname)
        n_facets <- length(unique(df_sub$omics))
        plot_width <- max(6, 5 * n_facets)
        ggplot2::ggsave(out_path, plot = p, width = plot_width, height = 6, dpi = 300)
        saved_paths <- c(saved_paths, out_path)
    }

    if (length(saved_paths) > 0) {
        message("  DIABLO log2FC loading plots saved: ", length(saved_paths), " contrasts")
    }
    saved_paths
}


#' Plot MOFA weights colored by log2FC from DE results
#'
#' Creates scatter plots of Factor1 vs Factor2 weights with points colored by
#' log2 fold-change. Faceted by omics view.
#'
#' @param mofa_results Output from run_mofa2_integration()
#' @param de_results Named list of DE results per omics
#' @param config Full pipeline config
#' @param out_dir Output directory
#' @param top_n Number of top features per view (by abs weight)
#' @param label_n Number of top features to label per view facet
#' @return Vector of saved PNG paths (one per contrast)
plot_mofa_loadings_log2fc <- function(mofa_results, de_results, config,
                                      out_dir, top_n = 50, label_n = 5,
                                      mae = NULL) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)

    weights <- mofa_results$weights
    feature_name_map <- mofa_results$feature_name_map

    if (is.null(weights) || length(weights) == 0 || is.null(de_results)) return(NULL)

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Build harmonized -> original ID mapping for display names and DE matching
    id_map <- build_harmonized_to_original_map(mae, de_results, config)

    # Build weight data frame across views
    df_list <- list()
    for (view in names(weights)) {
        w <- weights[[view]]
        if (ncol(w) < 2) next

        w1 <- w[, 1]
        w2 <- w[, 2]
        abs_max <- pmax(abs(w1), abs(w2))
        keep_idx <- order(abs_max, decreasing = TRUE)[seq_len(min(top_n, length(abs_max)))]

        feat_ids <- rownames(w)[keep_idx]

        # Resolve display names: feature_name_map first, then id_map fallback
        display <- if (!is.null(feature_name_map)) {
            nms <- feature_name_map[feat_ids]
            still_na <- is.na(nms)
            if (any(still_na)) {
                stripped <- sub("_(transcriptomics|proteomics|metabolomics)$", "",
                                feat_ids[still_na])
                nms[still_na] <- feature_name_map[stripped]
            }
            # Mark entries where name == id as unresolved
            ifelse(is.na(nms) | nms == "" | nms == feat_ids, NA_character_, nms)
        } else {
            rep(NA_character_, length(feat_ids))
        }
        # Fill remaining NAs from id_map (covers metabolomics HMDB + proteomics IDs)
        still_na <- is.na(display)
        if (any(still_na) && !is.null(id_map)) {
            # Try direct, then stripped
            display[still_na] <- id_map[feat_ids[still_na]]
            still_na2 <- is.na(display)
            if (any(still_na2)) {
                stripped <- sub("_(transcriptomics|proteomics|metabolomics)$", "",
                                feat_ids[still_na2])
                display[still_na2] <- id_map[stripped]
            }
        }
        display[is.na(display)] <- feat_ids[is.na(display)]

        df_list[[view]] <- data.frame(
            feature_id = feat_ids,
            display_name = display,
            factor1 = w1[keep_idx],
            factor2 = w2[keep_idx],
            omics = view,
            stringsAsFactors = FALSE
        )
    }
    if (length(df_list) == 0) return(NULL)
    df <- do.call(rbind, df_list)
    rownames(df) <- NULL

    contrast_names <- get_contrast_names_from_de(de_results)
    saved_paths <- character(0)

    for (ci in seq_along(contrast_names)) {
        contrast <- contrast_names[ci]
        df_plot <- df
        df_plot$log2FC <- NA_real_

        for (view in names(weights)) {
            if (is.null(de_results[[view]])) next
            de_norm <- tryCatch(
                extract_de_contrast(de_results[[view]], ci, view),
                error = function(e) NULL
            )
            if (is.null(de_norm)) next

            idx <- df_plot$omics == view
            feat_ids_view <- df_plot$feature_id[idx]
            matched <- match(feat_ids_view, de_norm$feature_id)
            # Try stripping MOFA view suffix
            still_na <- is.na(matched)
            if (any(still_na)) {
                stripped <- sub("_(transcriptomics|proteomics|metabolomics)$", "",
                                feat_ids_view[still_na])
                matched[still_na] <- match(stripped, de_norm$feature_id)
            }
            # Resolve harmonized IDs to original IDs via MAE rowData
            still_na <- is.na(matched)
            if (any(still_na) && !is.null(id_map)) {
                resolve_ids <- feat_ids_view[still_na]
                # Also try stripped version
                stripped2 <- sub("_(transcriptomics|proteomics|metabolomics)$", "",
                                 resolve_ids)
                orig_ids <- id_map[stripped2]
                orig_ids[is.na(orig_ids)] <- id_map[resolve_ids[is.na(orig_ids)]]
                valid <- !is.na(orig_ids)
                if (any(valid)) {
                    matched[still_na][valid] <- match(orig_ids[valid], de_norm$feature_id)
                }
            }
            df_plot$log2FC[idx] <- de_norm$logFC[matched]
        }

        df_sub <- df_plot[!is.na(df_plot$log2FC), ]
        if (nrow(df_sub) < 3) next

        lfc_cap <- quantile(abs(df_sub$log2FC), 0.95, na.rm = TRUE)
        lfc_cap <- max(lfc_cap, 0.5)
        df_sub$log2FC_capped <- pmax(pmin(df_sub$log2FC, lfc_cap), -lfc_cap)

        # Select top features per facet for labeling
        df_sub$label <- ""
        for (view in unique(df_sub$omics)) {
            v_idx <- which(df_sub$omics == view)
            top_idx <- v_idx[order(pmax(abs(df_sub$factor1[v_idx]),
                                        abs(df_sub$factor2[v_idx])),
                                   decreasing = TRUE)]
            top_idx <- head(top_idx, label_n)
            df_sub$label[top_idx] <- df_sub$display_name[top_idx]
        }

        p <- ggplot2::ggplot(df_sub,
                ggplot2::aes(x = factor1, y = factor2, color = log2FC_capped)) +
            ggplot2::geom_point(size = 2.5, alpha = 0.7) +
            ggplot2::scale_color_gradient2(
                low = "#4DBBD5", mid = "white", high = "#E64B35",
                midpoint = 0, name = expression(log[2]~FC),
                limits = c(-lfc_cap, lfc_cap)
            ) +
            ggplot2::facet_wrap(~omics, scales = "free") +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
            ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.grid = ggplot2::element_blank(),
                           strip.text = ggplot2::element_text(face = "bold")) +
            ggplot2::labs(
                title = "MOFA2 Weights Colored by log2FC",
                subtitle = contrast,
                x = "Factor 1 Weight",
                y = "Factor 2 Weight"
            )

        if (requireNamespace("ggrepel", quietly = TRUE)) {
            p <- p + ggrepel::geom_text_repel(
                ggplot2::aes(label = label),
                size = 2.8, color = "grey20", max.overlaps = 15,
                segment.color = "grey50", segment.size = 0.3
            )
        }

        fname <- paste0("mofa_loadings_log2fc_", make.names(contrast), ".png")
        out_path <- file.path(out_dir, fname)
        n_facets <- length(unique(df_sub$omics))
        plot_width <- max(6, 5 * n_facets)
        ggplot2::ggsave(out_path, plot = p, width = plot_width, height = 6, dpi = 300)
        saved_paths <- c(saved_paths, out_path)
    }

    if (length(saved_paths) > 0) {
        message("  MOFA2 log2FC weight plots saved: ", length(saved_paths), " contrasts")
    }
    saved_paths
}


#' Build a mapping from harmonized feature IDs to original IDs using MAE rowData
#'
#' Integration methods (MOFA, DIABLO) use harmonized IDs (e.g. GENE_5113)
#' while DE results use original IDs (e.g. WBGene00011800). This function
#' bridges the two via the MAE's rowData.
#'
#' @param mae MultiAssayExperiment (feature-selected subset)
#' @param de_results Named list of DE results per omics (for cross-referencing)
#' @param config Full config (for reading original data files)
#' @return Named character vector: harmonized_id -> original_id
build_harmonized_to_original_map <- function(mae, de_results = NULL, config = NULL) {
    if (is.null(mae)) return(NULL)

    all_mappings <- character(0)
    for (om in names(mae@ExperimentList)) {
        rd <- as.data.frame(SummarizedExperiment::rowData(mae[[om]]))
        harmonized_ids <- rownames(rd)

        # Try standard rowData columns
        orig_col <- intersect(c("original_id", "gene_id", "feature_id",
                                 "FeatureID", "protein_id"), colnames(rd))
        mapped <- FALSE
        if (length(orig_col) > 0) {
            orig_ids <- as.character(rd[[orig_col[1]]])
            valid <- !is.na(orig_ids) & orig_ids != ""
            if (sum(valid) > length(harmonized_ids) * 0.3) {
                mapping <- setNames(orig_ids[valid], harmonized_ids[valid])
                all_mappings <- c(all_mappings, mapping)
                mapped <- TRUE
            }
        }

        if (mapped) next

        # Fallback for metabolomics: recover HMDB IDs from abundance file
        if (grepl("metabol", om, ignore.case = TRUE) && !is.null(config)) {
            hmdb_map <- recover_metabolomics_hmdb_ids(config)
            if (!is.null(hmdb_map)) {
                matched <- hmdb_map[harmonized_ids]
                valid <- !is.na(matched)
                if (any(valid)) {
                    all_mappings <- c(all_mappings,
                                      setNames(matched[valid], harmonized_ids[valid]))
                }
            }
        }

        # Fallback for proteomics: try reading original DE file for ID column
        if (grepl("proteom", om, ignore.case = TRUE) && !is.null(de_results[[om]])) {
            prot_map <- recover_proteomics_ids(rd, de_results[[om]], config)
            if (!is.null(prot_map)) {
                all_mappings <- c(all_mappings, prot_map)
            }
        }
    }
    all_mappings
}


#' Recover proteomics feature IDs from original DE files
#'
#' When the summary_df uses numeric row indices as FeatureID, falls back
#' to reading the original DE CSV file to get the actual protein IDs,
#' then maps harmonized GENE_ IDs to protein IDs via Wormbase_id or gene_symbol.
#'
#' @param rd rowData data.frame for the proteomics MAE layer
#' @param de_res Proteomics DE result object
#' @param config Full config
#' @return Named character vector: harmonized_id -> protein_id, or NULL
recover_proteomics_ids <- function(rd, de_res, config = NULL) {
    harmonized_ids <- rownames(rd)

    # Check if we can match via a shared column (e.g. Wormbase_id, gene_symbol)
    # rd has columns like: ID, Wormbase_id, original_id, gene_symbol
    # DE file has: ID, Wormbase_id, ...

    # Try to read original DE file
    de_files <- NULL
    if (!is.null(config)) {
        de_files <- config$modes$proteomics$files$de_table
    }
    if (is.null(de_files) || length(de_files) == 0) return(NULL)

    proj_dir <- config$project$dir %||% "."
    data_dir <- file.path(proj_dir, config$paths$raw %||% "data")

    for (de_file in de_files) {
        fp <- if (file.exists(de_file)) de_file else file.path(data_dir, de_file)
        if (!file.exists(fp)) next

        de_df <- tryCatch(read.csv(fp, stringsAsFactors = FALSE, nrows = 5000),
                          error = function(e) NULL)
        if (is.null(de_df)) next

        # The DE file should have protein ID and possibly Wormbase_id
        id_col <- intersect(c("ID", "protein_id", "FeatureID"), colnames(de_df))[1]
        if (is.na(id_col)) next

        # Try matching via Wormbase_id (shared between rd and de_df)
        if ("Wormbase_id" %in% colnames(rd) && "Wormbase_id" %in% colnames(de_df)) {
            rd_wb <- as.character(rd[["Wormbase_id"]])
            de_wb <- as.character(de_df[["Wormbase_id"]])
            de_ids <- as.character(de_df[[id_col]])

            matched_idx <- match(rd_wb, de_wb)
            valid <- !is.na(matched_idx) & !is.na(rd_wb) & rd_wb != ""
            if (sum(valid) > 0) {
                return(setNames(de_ids[matched_idx[valid]], harmonized_ids[valid]))
            }
        }

        # Try matching via original_id (UniProt IDs in rd vs DE file)
        if ("original_id" %in% colnames(rd)) {
            rd_orig <- as.character(rd[["original_id"]])
            de_ids <- as.character(de_df[[id_col]])

            matched_idx <- match(rd_orig, de_ids)
            valid <- !is.na(matched_idx) & !is.na(rd_orig) & rd_orig != ""
            if (sum(valid) > 0) {
                return(setNames(de_ids[matched_idx[valid]], harmonized_ids[valid]))
            }
        }

        break  # Only need the first valid file
    }
    NULL
}


#' Extract contrast names from DE results list
#'
#' @param de_results Named list of DE results per omics
#' @return Character vector of contrast names
get_contrast_names_from_de <- function(de_results) {
    for (om in names(de_results)) {
        de <- de_results[[om]]
        if (is.null(de)) next

        # RNA-seq: list(tables = named list of per-contrast DFs)
        if (!is.null(de$tables) && is.list(de$tables)) {
            nms <- names(de$tables)
            if (length(nms) > 0) return(nms)
        }

        # Proteomics/metabolomics: summary_df with wide logFC columns
        if (!is.null(de$summary_df) && is.data.frame(de$summary_df)) {
            lfc_cols <- grep("^(linearFC\\.imputs\\.|logFC\\.|log2FoldChange\\.|logFC_)",
                             names(de$summary_df), value = TRUE)
            if (length(lfc_cols) > 0) {
                contrasts <- sub("^(linearFC\\.imputs\\.|logFC\\.|log2FoldChange\\.|logFC_)",
                                 "", lfc_cols)
                return(contrasts)
            }
        }
    }
    "contrast_1"
}


#' Extract DE data for a specific contrast index from a single-omics DE result
#'
#' @param de_res Single omics DE result object
#' @param contrast_idx Contrast index (integer)
#' @param omics_type Label for messages
#' @return data.frame with feature_id and logFC columns
extract_de_contrast <- function(de_res, contrast_idx, omics_type = "unknown") {

    # RNA-seq: list(tables = list of DFs)
    if (!is.null(de_res$tables) && is.list(de_res$tables)) {
        idx <- min(contrast_idx, length(de_res$tables))
        tbl <- de_res$tables[[idx]]
        lfc_col <- intersect(c("log2FoldChange", "logFC"), names(tbl))[1]
        id_col <- intersect(c("FeatureID", "gene_id", "feature_id"), names(tbl))[1]
        id_vals <- if (is.na(id_col)) rownames(tbl) else tbl[[id_col]]
        if (is.na(lfc_col)) return(NULL)
        return(data.frame(feature_id = id_vals, logFC = tbl[[lfc_col]],
                          stringsAsFactors = FALSE))
    }

    # Proteomics/metabolomics: summary_df with wide columns
    if (!is.null(de_res$summary_df) && is.data.frame(de_res$summary_df)) {
        df <- de_res$summary_df
        id_col <- intersect(c("FeatureID", "feature_id", "protein_id"), names(df))[1]
        id_vals <- if (is.na(id_col)) rownames(df) else df[[id_col]]

        lfc_cols <- grep("^(linearFC\\.imputs\\.|logFC\\.|log2FoldChange\\.|logFC_)",
                         names(df), value = TRUE)
        if (length(lfc_cols) == 0) return(NULL)

        col_idx <- min(contrast_idx, length(lfc_cols))
        lfc_vals <- df[[lfc_cols[col_idx]]]

        if (grepl("^linearFC", lfc_cols[col_idx])) {
            lfc_vals <- ifelse(lfc_vals >= 0, log2(pmax(lfc_vals, 1e-10)),
                               -log2(pmax(abs(lfc_vals), 1e-10)))
        }
        return(data.frame(feature_id = id_vals, logFC = lfc_vals,
                          stringsAsFactors = FALSE))
    }

    # Flat data frame
    if (is.data.frame(de_res)) {
        if (!"feature_id" %in% names(de_res)) de_res$feature_id <- rownames(de_res)
        lfc_col <- intersect(c("logFC", "log2FoldChange"), names(de_res))[1]
        if (is.na(lfc_col)) return(NULL)
        return(data.frame(feature_id = de_res$feature_id, logFC = de_res[[lfc_col]],
                          stringsAsFactors = FALSE))
    }

    NULL
}
