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

        # Variable plot (loadings)
        plots$variable_plot <- file.path(out_dir, "diablo_variable_plot.png")
        png(plots$variable_plot, width = 1000, height = 800, res = 120)
        tryCatch({
            mixOmics::plotVar(diablo_model, comp = c(1, 2), style = "graphics",
                              legend = TRUE, title = "DIABLO Variable Loadings")
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("Plot failed:", e$message), cex = 1.2)
        })
        dev.off()

        # Circos plot (correlations between omics)
        plots$circos_plot <- file.path(out_dir, "diablo_circos_plot.png")
        png(plots$circos_plot, width = 1000, height = 1000, res = 120)
        tryCatch({
            mixOmics::circosPlot(diablo_model, cutoff = 0.5, size.variables = 0.5)
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("Plot failed:", e$message), cex = 1.2)
        })
        dev.off()

        message("  DIABLO plots saved to: ", out_dir)
    }

    # Extract top features per component
    top_features <- extract_diablo_top_features(diablo_model, top_n = 50)

    message(sprintf(
        "DIABLO integration complete: %d components, %d omics layers",
        ncomp, length(omics)
    ))

    list(
        model = diablo_model,
        sample_scores = sample_scores,
        loadings = loadings,
        top_features = top_features,
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
#' @return Named list of data frames (one per omics)
extract_diablo_top_features <- function(diablo_model, top_n = 50) {

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

            top_list[[comp_name]] <- data.frame(
                feature = rownames(load_mat)[ranked_idx[seq_len(n_select)]],
                loading = comp_loadings[ranked_idx[seq_len(n_select)]],
                abs_loading = abs(comp_loadings[ranked_idx[seq_len(n_select)]]),
                component = comp_name,
                stringsAsFactors = FALSE
            )
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
