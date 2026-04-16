#' MOFA2 (Multi-Omics Factor Analysis) integration
#'
#' MOFA2 is an unsupervised multi-omics integration method that identifies
#' latent factors capturing major axes of variation across omics layers.
#'
#' @name mofa2_integration
NULL

#' Run MOFA2 integration analysis
#'
#' @param mae MultiAssayExperiment object with aligned samples
#' @param config Full config object
#' @param out_dir Output directory for plots and results
#' @return List with: model, factors, weights, variance_explained, top_features, plots, config
run_mofa2_integration <- function(mae, config, out_dir = NULL) {

    if (!requireNamespace("MOFA2", quietly = TRUE)) {
        stop(
            "Package 'MOFA2' is required for MOFA2 integration. ",
            "Install with: BiocManager::install('MOFA2')",
            call. = FALSE
        )
    }

    # Extract MOFA2 configuration
    mofa_config <- config$modes$multiomics$integration$mofa2 %||% list()
    num_factors <- mofa_config$num_factors %||% 10
    convergence_mode <- mofa_config$convergence_mode %||% "fast"
    seed <- mofa_config$seed %||% 42

    message("Running MOFA2 integration...")

    # Extract matrices from MAE (features x samples)
    omics <- names(MultiAssayExperiment::experiments(mae))
    matrices <- lapply(omics, function(nm) {
        as.matrix(SummarizedExperiment::assay(mae[[nm]]))
    })
    names(matrices) <- omics

    # Get sample metadata
    coldata <- SummarizedExperiment::colData(mae)
    metadata <- as.data.frame(coldata)

    # Check we have at least 2 omics
    if (length(matrices) < 2) {
        message("MOFA2 requires at least 2 omics layers. Found: ", length(matrices))
        return(NULL)
    }

    # Ensure all matrices share the same samples in the same order
    common_samples <- Reduce(intersect, lapply(matrices, colnames))
    if (length(common_samples) == 0) {
        message("No common samples across omics layers for MOFA2.")
        return(NULL)
    }
    matrices <- lapply(matrices, function(m) m[, common_samples, drop = FALSE])

    message("Preparing data for MOFA2...")
    message("  Number of omics layers: ", length(matrices))
    message("  Number of samples: ", length(common_samples))

    # Create MOFA object from matrices
    mofa_data <- list()
    for (omic in names(matrices)) {
        mat <- matrices[[omic]]

        # Sanitize feature names: replace non-ASCII characters
        rn <- rownames(mat)
        rn_clean <- iconv(rn, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
        if (any(rn != rn_clean)) {
            n_fixed <- sum(rn != rn_clean)
            message("  Sanitized ", n_fixed, " non-ASCII feature names in ", omic)
            rownames(mat) <- rn_clean
        }

        message("  ", omic, ": ", nrow(mat), " features x ", ncol(mat), " samples")

        # Center features (MOFA2 expects centered data for Gaussian likelihood)
        mat_centered <- t(scale(t(mat), center = TRUE, scale = FALSE))
        # Handle remaining NAs
        mat_centered[is.na(mat_centered)] <- 0
        mofa_data[[omic]] <- mat_centered
    }

    # Create MOFA object
    mofa_obj <- tryCatch({
        MOFA2::create_mofa(mofa_data)
    }, error = function(e) {
        message("Error creating MOFA object: ", e$message)
        return(NULL)
    })

    if (is.null(mofa_obj)) return(NULL)

    # Prepare options using MOFA2 API
    data_opts <- MOFA2::get_default_data_options(mofa_obj)

    model_opts <- MOFA2::get_default_model_options(mofa_obj)
    model_opts$num_factors <- min(num_factors, ncol(matrices[[1]]) - 1)

    train_opts <- MOFA2::get_default_training_options(mofa_obj)
    train_opts$convergence_mode <- convergence_mode
    train_opts$seed <- seed
    train_opts$verbose <- FALSE

    # Use prepare_mofa to set all options at once
    mofa_obj <- tryCatch({
        MOFA2::prepare_mofa(
            mofa_obj,
            data_options = data_opts,
            model_options = model_opts,
            training_options = train_opts
        )
    }, error = function(e) {
        message("Error preparing MOFA object: ", e$message)
        return(NULL)
    })

    if (is.null(mofa_obj)) return(NULL)

    # Run MOFA
    message("Training MOFA2 model with ", model_opts$num_factors, " factors...")
    message("  Convergence mode: ", convergence_mode)
    message("  Max iterations: ", train_opts$iter)
    message("  Random seed: ", seed)

    mofa_trained <- tryCatch({
        MOFA2::run_mofa(mofa_obj, use_basilisk = FALSE)
    }, error = function(e) {
        message("Error training MOFA model without basilisk: ", e$message)
        message("  Retrying with basilisk environment...")
        tryCatch({
            MOFA2::run_mofa(mofa_obj, use_basilisk = TRUE)
        }, error = function(e2) {
            message("MOFA training failed with basilisk: ", e2$message)
            return(NULL)
        })
    })

    if (is.null(mofa_trained)) return(NULL)

    message("MOFA2 training complete successfully")

    # Extract results
    mofa_results <- extract_mofa_results(mofa_trained, metadata, config)

    # Build feature name map from MAE rowData for display labels
    feature_name_map <- build_feature_name_map(mae)

    # Create visualizations
    plots <- list()
    if (!is.null(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        plots <- create_mofa_plots(mofa_trained, mofa_results, metadata, config,
                                    out_dir, feature_name_map = feature_name_map)
    }

    list(
        model = mofa_trained,
        factors = mofa_results$factor_df,
        weights = mofa_results$weights,
        variance_explained = mofa_results$var_df,
        top_features = mofa_results$top_features,
        feature_name_map = feature_name_map,
        r2_per_factor = mofa_results$r2_per_factor,
        r2_total = mofa_results$r2_total,
        plots = plots,
        config = mofa_config
    )
}


#' Extract results from trained MOFA model
#'
#' @param mofa_model Trained MOFA model object
#' @param metadata Sample metadata data.frame
#' @param config Full config object
#' @return List with extracted results
extract_mofa_results <- function(mofa_model, metadata, config) {
    message("Extracting MOFA2 results...")

    # Factor values (samples x factors)
    factors <- MOFA2::get_factors(mofa_model)[[1]]
    message("  Extracted factors: ", nrow(factors), " samples x ", ncol(factors), " factors")

    # Weights per view (features x factors per view)
    weights <- MOFA2::get_weights(mofa_model)
    message("  Extracted weights for ", length(weights), " views")

    # Variance explained
    var_explained <- MOFA2::get_variance_explained(mofa_model)

    # Per factor per view
    r2_per_factor <- var_explained$r2_per_factor[[1]]

    # Total per view
    r2_total <- var_explained$r2_total[[1]]
    message("  Total variance explained across views:")
    for (view_name in names(r2_total)) {
        message("    ", view_name, ": ", round(r2_total[[view_name]], 2), "%")
    }

    # Create factor data.frame
    factor_df <- as.data.frame(factors)
    factor_df$sample_id <- rownames(factors)
    factor_df <- factor_df[, c("sample_id", setdiff(colnames(factor_df), "sample_id"))]

    # Add metadata
    condition_col <- config$modes$multiomics$condition_column %||%
                     config$design$condition_column %||%
                     "condition"
    if (condition_col %in% colnames(metadata)) {
        factor_df[[condition_col]] <- metadata[factor_df$sample_id, condition_col]
    }

    # Create variance explained data.frame
    var_df <- as.data.frame(r2_per_factor)
    var_df$view <- rownames(var_df)
    var_df <- var_df[, c("view", setdiff(colnames(var_df), "view"))]

    # Identify top features per factor per view
    top_features <- extract_top_mofa_features(weights, n_top = 50)

    list(
        factors = factors,
        factor_df = factor_df,
        weights = weights,
        variance_explained = var_explained,
        r2_per_factor = r2_per_factor,
        r2_total = r2_total,
        var_df = var_df,
        top_features = top_features
    )
}


#' Extract top features per factor from MOFA weights
#'
#' @param weights Named list of weight matrices per view
#' @param n_top Number of top features per factor
#' @return data.frame of top features
extract_top_mofa_features <- function(weights, n_top = 50) {
    all_top <- list()

    for (view in names(weights)) {
        w <- weights[[view]]
        n_factors <- ncol(w)

        view_top <- list()
        for (k in seq_len(n_factors)) {
            factor_name <- colnames(w)[k]
            loadings <- w[, k]

            # Get top by absolute value
            ord <- order(abs(loadings), decreasing = TRUE)
            top_idx <- head(ord, n_top)

            view_top[[factor_name]] <- data.frame(
                feature_id = rownames(w)[top_idx],
                weight = loadings[top_idx],
                abs_weight = abs(loadings[top_idx]),
                factor = factor_name,
                view = view,
                stringsAsFactors = FALSE
            )
        }

        all_top[[view]] <- do.call(rbind, view_top)
    }

    do.call(rbind, all_top)
}


#' Create MOFA visualization plots
#'
#' @param mofa_model Trained MOFA model
#' @param mofa_results Extracted results from extract_mofa_results()
#' @param metadata Sample metadata
#' @param config Full config object
#' @param out_dir Output directory
#' @return List of plot file paths
create_mofa_plots <- function(mofa_model, mofa_results, metadata, config, out_dir,
                              feature_name_map = NULL) {
    message("Creating MOFA2 plots...")

    condition_col <- config$modes$multiomics$condition_column %||%
                     config$design$condition_column %||%
                     "condition"

    plots <- list()

    # 1. Variance explained heatmap
    tryCatch({
        p <- plot_mofa_variance_heatmap(mofa_results$r2_per_factor)
        plots$variance_heatmap <- p
        ggplot2::ggsave(
            file.path(out_dir, "mofa_variance_heatmap.png"),
            plot = p, width = 8, height = 6, dpi = 300
        )
    }, error = function(e) message("Failed to create variance heatmap: ", e$message))

    # 2. Factor scatter plots
    tryCatch({
        factors <- mofa_results$factors
        if (ncol(factors) >= 2) {
            p <- plot_mofa_factors(factors, metadata, condition_col, c(1, 2))
            plots$factors_1_2 <- p
            ggplot2::ggsave(
                file.path(out_dir, "mofa_factors_1_2.png"),
                plot = p, width = 8, height = 6, dpi = 300
            )
        }
        if (ncol(factors) >= 3) {
            p <- plot_mofa_factors(factors, metadata, condition_col, c(1, 3))
            plots$factors_1_3 <- p
            ggplot2::ggsave(
                file.path(out_dir, "mofa_factors_1_3.png"),
                plot = p, width = 8, height = 6, dpi = 300
            )
        }
    }, error = function(e) message("Failed to create factor plots: ", e$message))

    # 3. Top weights per factor
    tryCatch({
        for (view in names(mofa_results$weights)) {
            p <- plot_mofa_top_weights(mofa_results$weights[[view]], view, n_top = 20,
                                      feature_name_map = feature_name_map)
            plots[[paste0("top_weights_", view)]] <- p
            ggplot2::ggsave(
                file.path(out_dir, paste0("mofa_top_weights_", view, ".png")),
                plot = p, width = 10, height = 8, dpi = 300
            )
        }
    }, error = function(e) message("Failed to create weight plots: ", e$message))

    message("MOFA2 plots complete")
    plots
}


#' Plot variance explained heatmap
#'
#' @param r2_matrix Variance explained matrix (views x factors)
#' @return ggplot object
plot_mofa_variance_heatmap <- function(r2_matrix) {
    # Convert to long format
    df <- as.data.frame(r2_matrix)
    df$view <- rownames(df)
    df_long <- tidyr::pivot_longer(df, -view, names_to = "factor", values_to = "r2")

    ggplot2::ggplot(df_long, ggplot2::aes(x = factor, y = view, fill = r2)) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", r2)), size = 3) +
        ggplot2::scale_fill_gradient(low = "white", high = "steelblue", name = "R\u00B2 (%)") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "MOFA2: Variance Explained per Factor",
            x = "Factor",
            y = "View (Omics)"
        ) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid = ggplot2::element_blank()
        )
}


#' Plot MOFA factors scatter
#'
#' @param factors Factor matrix (samples x factors)
#' @param metadata Sample metadata
#' @param condition_col Condition column name
#' @param factor_idx Vector of two factor indices to plot
#' @return ggplot object
plot_mofa_factors <- function(factors, metadata, condition_col, factor_idx = c(1, 2)) {
    f1 <- factor_idx[1]
    f2 <- factor_idx[2]

    df <- data.frame(
        Factor1 = factors[, f1],
        Factor2 = factors[, f2],
        sample_id = rownames(factors),
        stringsAsFactors = FALSE
    )

    if (!is.null(condition_col) && condition_col %in% colnames(metadata)) {
        df$condition <- metadata[df$sample_id, condition_col]
    } else {
        df$condition <- "All"
    }

    ggplot2::ggplot(df, ggplot2::aes(x = Factor1, y = Factor2, color = condition)) +
        ggplot2::geom_point(size = 3, alpha = 0.8) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = paste0("MOFA2: Factor ", f1, " vs Factor ", f2),
            x = paste0("Factor ", f1),
            y = paste0("Factor ", f2),
            color = condition_col
        )
}


#' Plot top weights for a view
#'
#' @param weights Weight matrix for a single view
#' @param view_name Name of the view
#' @param n_top Number of top features per factor
#' @param n_factors Number of factors to plot
#' @return ggplot object (or patchwork if available)
plot_mofa_top_weights <- function(weights, view_name, n_top = 20, n_factors = 3,
                                  feature_name_map = NULL) {
    n_factors <- min(n_factors, ncol(weights))
    plots <- list()

    for (k in seq_len(n_factors)) {
        w <- weights[, k]
        ord <- order(abs(w), decreasing = TRUE)
        top_idx <- head(ord, n_top)

        feat_ids <- rownames(weights)[top_idx]
        # Drop entries with missing rownames
        valid <- !is.na(feat_ids) & feat_ids != ""
        feat_ids <- feat_ids[valid]
        top_idx <- top_idx[valid]
        if (length(feat_ids) == 0) next

        # Resolve to original names if map is available
        # Strip view suffix first (MOFA appends _viewname for shared features)
        stripped_ids <- sub("_(transcriptomics|proteomics|metabolomics)$", "", feat_ids)
        display_names <- if (!is.null(feature_name_map)) {
            # Try both full and stripped IDs
            mapped <- feature_name_map[feat_ids]
            still_na <- is.na(mapped)
            if (any(still_na)) {
                mapped[still_na] <- feature_name_map[stripped_ids[still_na]]
            }
            ifelse(is.na(mapped), stripped_ids, mapped)
        } else {
            stripped_ids
        }

        # Guard against NA/empty display names
        display_names[is.na(display_names) | display_names == ""] <- feat_ids[is.na(display_names) | display_names == ""]

        # For transcriptomics: ensure GL IDs are shown (not KAE protein IDs)
        if (grepl("transcriptomics", view_name, ignore.case = TRUE)) {
            # Strip suffix and look up feature_id from the name map's original keys
            stripped_ids <- sub("_(transcriptomics|proteomics|metabolomics)$", "", feat_ids)
            # If display shows KAE (protein IDs), try to resolve to GL IDs
            is_protein_id <- grepl("^KAE", display_names)
            if (any(is_protein_id) && !is.null(feature_name_map)) {
                # The feature_name_map maps GENE_xxx to display names
                # For transcriptomics, we want feature_id (GL IDs) from rowData
                # Just use stripped_ids as fallback (GENE_xxx without suffix)
                display_names[is_protein_id] <- stripped_ids[is_protein_id]
            }
            # Also fix any remaining GENE_xxx
            still_generic <- grepl("^GENE_\\d+$", display_names)
            if (any(still_generic)) {
                display_names[still_generic] <- stripped_ids[still_generic]
            }
        }

        # For any view: strip GENE_xxx_viewname suffix
        still_generic <- grepl("^GENE_\\d+_", display_names)
        if (any(still_generic)) {
            display_names[still_generic] <- sub("_(transcriptomics|proteomics|metabolomics)$", "",
                                                 display_names[still_generic])
        }

        # Truncate long metabolite names to first 25 chars
        long_names <- nchar(display_names) > 25
        if (any(long_names)) {
            display_names[long_names] <- paste0(substr(display_names[long_names], 1, 22), "...")
        }

        # Make unique to avoid factor() issues with duplicates
        display_names <- make.unique(display_names, sep = " ")

        df <- data.frame(
            feature = display_names,
            weight = w[top_idx],
            stringsAsFactors = FALSE,
            row.names = NULL
        )
        df$feature <- factor(df$feature, levels = df$feature[order(abs(df$weight))])
        df$direction <- ifelse(df$weight > 0, "positive", "negative")

        plots[[k]] <- ggplot2::ggplot(df, ggplot2::aes(x = weight, y = feature, fill = direction)) +
            ggplot2::geom_col() +
            ggplot2::scale_fill_manual(values = c("positive" = "firebrick", "negative" = "steelblue")) +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = paste0("Factor ", k),
                x = "Weight",
                y = NULL
            ) +
            ggplot2::theme(legend.position = "none")
    }

    # Combine plots
    if (length(plots) > 1 && requireNamespace("patchwork", quietly = TRUE)) {
        combined <- patchwork::wrap_plots(plots, ncol = n_factors) +
            patchwork::plot_annotation(title = paste0("MOFA2 Top Weights: ", view_name))
        return(combined)
    } else {
        return(plots[[1]])
    }
}


#' Test factor associations with phenotypes
#'
#' @param mofa_results Results from run_mofa2_integration()
#' @param metadata Sample metadata
#' @param config Full config object
#' @param out_dir Output directory
#' @return data.frame of association test results
test_mofa_factor_associations <- function(mofa_results, metadata, config, out_dir = NULL) {
    message("Testing MOFA factor associations with phenotypes...")

    factors <- mofa_results$factors
    condition_col <- config$modes$multiomics$condition_column %||%
                     config$design$condition_column %||%
                     "condition"

    results <- list()

    # Test condition association
    if (condition_col %in% colnames(metadata)) {
        condition <- metadata[rownames(factors), condition_col]

        for (k in seq_len(ncol(factors))) {
            factor_name <- colnames(factors)[k]
            factor_values <- factors[, k]

            # t-test or ANOVA depending on number of groups
            n_groups <- length(unique(condition))

            if (n_groups == 2) {
                test_result <- t.test(factor_values ~ condition)
                results[[factor_name]] <- data.frame(
                    factor = factor_name,
                    test = "t-test",
                    statistic = test_result$statistic,
                    pvalue = test_result$p.value,
                    variable = condition_col,
                    stringsAsFactors = FALSE
                )
            } else if (n_groups > 2) {
                test_result <- summary(aov(factor_values ~ condition))[[1]]
                results[[factor_name]] <- data.frame(
                    factor = factor_name,
                    test = "ANOVA",
                    statistic = test_result$`F value`[1],
                    pvalue = test_result$`Pr(>F)`[1],
                    variable = condition_col,
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    if (length(results) > 0) {
        result_df <- do.call(rbind, results)
        result_df$padj <- p.adjust(result_df$pvalue, method = "BH")

        if (!is.null(out_dir)) {
            write.csv(result_df, file.path(out_dir, "mofa_factor_associations.csv"),
                      row.names = FALSE)
        }
        return(result_df)
    }

    return(NULL)
}


#' Write MOFA results to CSV files
#'
#' @param mofa_results Output from run_mofa2_integration()
#' @param out_dir Output directory
write_mofa_results <- function(mofa_results, out_dir) {

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Sample factor values
    if (!is.null(mofa_results$factors)) {
        write.csv(mofa_results$factors,
                  file.path(out_dir, "mofa_factors.csv"),
                  row.names = FALSE)
    }

    # Weights per view
    if (!is.null(mofa_results$weights)) {
        name_map <- mofa_results$feature_name_map
        for (view in names(mofa_results$weights)) {
            w <- mofa_results$weights[[view]]
            w_df <- as.data.frame(w)
            w_df$feature_id <- rownames(w)
            if (!is.null(name_map)) {
                w_df$original_name <- unname(name_map[rownames(w)])
            }
            lead_cols <- intersect(c("feature_id", "original_name"), colnames(w_df))
            w_df <- w_df[, c(lead_cols, setdiff(colnames(w_df), lead_cols))]
            write.csv(w_df,
                      file.path(out_dir, paste0("mofa_weights_", view, ".csv")),
                      row.names = FALSE)
        }
    }

    # Variance explained
    if (!is.null(mofa_results$variance_explained)) {
        write.csv(mofa_results$variance_explained,
                  file.path(out_dir, "mofa_variance_explained.csv"),
                  row.names = FALSE)
    }

    # Top features
    if (!is.null(mofa_results$top_features)) {
        top_df <- mofa_results$top_features
        if (!is.null(name_map)) {
            top_df$original_name <- unname(name_map[top_df$feature_id])
        }
        write.csv(top_df,
                  file.path(out_dir, "mofa_top_features.csv"),
                  row.names = FALSE)
    }

    message("MOFA results written to: ", out_dir)

    invisible(NULL)
}


# =============================================================================
# Python Wrapper for MOFA+ (alternative implementation using mofapy2)
# =============================================================================

#' Run MOFA+ analysis using Python mofapy2 via subprocess
#'
#' @param views Named list of file paths to CSV files containing omics data.
#'   Each CSV should have features as rows and samples as columns.
#' @param outfile Path to save the trained MOFA model (.hdf5 file)
#' @param outdir Directory to save CSV outputs (factors, weights, R2)
#' @param n_factors Number of factors to train (default: 10)
#' @param seed Random seed for reproducibility (default: 42)
#' @param max_iter Maximum number of training iterations (default: 1000)
#' @param convergence_mode Convergence mode: "fast", "medium", or "slow" (default: "fast")
#' @param scale_views Logical, whether to scale views to unit variance (default: FALSE)
#' @param verbose Logical, whether to print verbose output (default: TRUE)
#' @param python_exec Path to Python executable (default: from RETICULATE_PYTHON env var)
#'
#' @return List containing:
#'   \item{success}{Logical indicating if MOFA ran successfully}
#'   \item{model_file}{Path to saved HDF5 model file}
#'   \item{factors}{Data frame of latent factors (samples x factors)}
#'   \item{weights}{Named list of data frames with weights for each view}
#'   \item{variance_explained}{Data frame of variance explained (R2) per view}
#'
#' @export
run_mofa <- function(views,
                     outfile,
                     outdir,
                     n_factors = 10,
                     seed = 42,
                     max_iter = 1000,
                     convergence_mode = c("fast", "medium", "slow"),
                     scale_views = FALSE,
                     verbose = TRUE,
                     python_exec = Sys.getenv("RETICULATE_PYTHON", "python3")) {

    # Input validation
    if (!is.list(views) || length(views) == 0) {
        stop("views must be a non-empty named list")
    }

    if (is.null(names(views)) || any(names(views) == "")) {
        stop("views must be a named list (e.g., list(rna = 'file.csv', protein = 'file2.csv'))")
    }

    # Check that all view files exist
    for (view_name in names(views)) {
        if (!file.exists(views[[view_name]])) {
            stop("View file not found: ", views[[view_name]])
        }
    }

    # Check Python script exists
    script_path <- "scripts/run_mofa.py"
    if (!file.exists(script_path)) {
        stop("MOFA script not found: ", script_path,
             "\nPlease ensure run_mofa.py is in the scripts directory")
    }

    convergence_mode <- match.arg(convergence_mode)

    # Create output directory
    if (!dir.exists(dirname(outfile))) {
        dir.create(dirname(outfile), recursive = TRUE)
    }
    if (!dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
    }

    # Build view arguments
    view_args <- sapply(names(views), function(name) {
        sprintf("%s=%s", name, views[[name]])
    })
    view_args_str <- paste(view_args, collapse = " ")

    # Build command
    cmd <- sprintf(
        "%s '%s' --views %s --outfile '%s' --outdir '%s' --factors %d --seed %d --max_iter %d --convergence_mode %s",
        python_exec,
        script_path,
        view_args_str,
        outfile,
        outdir,
        n_factors,
        seed,
        max_iter,
        convergence_mode
    )

    if (scale_views) {
        cmd <- paste(cmd, "--scale_views")
    }

    if (verbose) {
        cmd <- paste(cmd, "--verbose")
    }

    if (verbose) {
        message("Running MOFA+ analysis...")
        message("  Views: ", paste(names(views), collapse = ", "))
        message("  Factors: ", n_factors)
        message("  Output: ", outdir)
        message()
    }

    # Execute command
    exit_code <- system(cmd, intern = FALSE)

    # Check if successful
    if (exit_code != 0) {
        stop("MOFA script failed with exit code ", exit_code)
    }

    # Load outputs
    if (verbose) {
        message("\nLoading MOFA results...")
    }

    result <- list(
        success = TRUE,
        model_file = outfile
    )

    # Load factors
    factors_file <- file.path(outdir, "factors.csv")
    if (file.exists(factors_file)) {
        result$factors <- read.csv(factors_file, row.names = 1)
        if (verbose) {
            message("  Loaded factors: ", nrow(result$factors), " samples x ",
                    ncol(result$factors), " factors")
        }
    } else {
        warning("Factors file not found: ", factors_file)
    }

    # Load weights for each view
    result$weights <- list()
    for (view_name in names(views)) {
        weights_file <- file.path(outdir, paste0("weights_", view_name, ".csv"))
        if (file.exists(weights_file)) {
            result$weights[[view_name]] <- read.csv(weights_file, row.names = 1)
            if (verbose) {
                message("  Loaded weights for ", view_name, ": ",
                        nrow(result$weights[[view_name]]), " features x ",
                        ncol(result$weights[[view_name]]), " factors")
            }
        } else {
            warning("Weights file not found for view ", view_name, ": ", weights_file)
        }
    }

    # Load variance explained
    r2_total_file <- file.path(outdir, "variance_explained_total.csv")
    r2_per_factor_file <- file.path(outdir, "variance_explained_per_factor.csv")

    if (file.exists(r2_total_file)) {
        result$variance_explained_total <- read.csv(r2_total_file, row.names = 1)
        if (verbose) {
            message("\n  Variance explained (R2):")
            for (view_name in names(views)) {
                if (view_name %in% colnames(result$variance_explained_total)) {
                    r2 <- result$variance_explained_total["Total_R2", view_name]
                    message("    ", view_name, ": ", round(r2, 2), "%")
                }
            }
        }
    }

    if (file.exists(r2_per_factor_file)) {
        result$variance_explained_per_factor <- read.csv(r2_per_factor_file, row.names = 1)
    }

    if (verbose) {
        message("\nMOFA analysis completed successfully!")
    }

    return(result)
}


#' Prepare abundance matrix for MOFA
#'
#' Helper function to prepare an abundance matrix for MOFA analysis by saving
#' it as a CSV file.
#'
#' @param mat Matrix or data.frame with features as rows and samples as columns
#' @param outfile Path to save CSV file
#' @param feature_col Optional column name for feature IDs (if mat is a data.frame)
#' @return Path to saved file (invisibly)
#' @export
prepare_mofa_view <- function(mat, outfile, feature_col = NULL) {

    # Convert to data.frame if needed
    if (is.matrix(mat)) {
        df <- as.data.frame(mat)
    } else if (is.data.frame(mat)) {
        df <- mat
    } else {
        stop("mat must be a matrix or data.frame")
    }

    # Set feature IDs as rownames if specified
    if (!is.null(feature_col)) {
        if (!(feature_col %in% colnames(df))) {
            stop("feature_col '", feature_col, "' not found in data.frame")
        }
        rownames(df) <- df[[feature_col]]
        df[[feature_col]] <- NULL
    }

    # Ensure we have row names
    if (is.null(rownames(df)) || any(rownames(df) == "")) {
        stop("mat must have row names (feature IDs)")
    }

    # Create output directory if needed
    if (!dir.exists(dirname(outfile))) {
        dir.create(dirname(outfile), recursive = TRUE)
    }

    # Write CSV
    write.csv(df, outfile, row.names = TRUE)

    message("Saved MOFA view: ", outfile)
    message("  ", nrow(df), " features x ", ncol(df), " samples")

    invisible(outfile)
}


#' Plot MOFA factor variance explained (for Python wrapper results)
#'
#' Create a barplot showing variance explained by each factor per view
#'
#' @param mofa_result Result object from run_mofa()
#' @param top_n Number of top factors to plot (default: all)
#' @return ggplot object
#' @export
plot_mofa_variance <- function(mofa_result, top_n = NULL) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package required for plotting")
    }

    if (!requireNamespace("tidyr", quietly = TRUE)) {
        stop("tidyr package required for plotting")
    }

    if (is.null(mofa_result$variance_explained_per_factor)) {
        stop("Variance explained data not found in MOFA result")
    }

    df <- mofa_result$variance_explained_per_factor

    # Subset to top N factors if requested
    if (!is.null(top_n) && top_n < nrow(df)) {
        # Sum R2 across views for each factor
        factor_totals <- rowSums(df, na.rm = TRUE)
        top_factors <- names(sort(factor_totals, decreasing = TRUE)[1:top_n])
        df <- df[top_factors, , drop = FALSE]
    }

    # Convert to long format
    df$Factor <- rownames(df)
    df_long <- tidyr::pivot_longer(
        df,
        cols = -Factor,
        names_to = "View",
        values_to = "R2"
    )

    # Create plot
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Factor, y = R2, fill = View)) +
        ggplot2::geom_bar(stat = "identity", position = "dodge") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "Variance Explained by MOFA Factors",
            x = "Factor",
            y = "Variance Explained (R\u00B2%)",
            fill = "View"
        ) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )

    return(p)
}
