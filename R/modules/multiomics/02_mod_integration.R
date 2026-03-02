#' Module: Multi-omics integration
#'
#' Runs integration methods (DIABLO, SNF, MOFA2) based on config.
#' Includes feature selection, integration analysis, and result exports.
#'
#' @param harmonization_res Output from mod_multiomics_harmonization()
#' @param de_results Named list of DE results per omics (for feature selection)
#' @param config Full config object
#' @param out_dir Output directory
#' @return List with: diablo_results, snf_results, mofa_results, feature_selection, mae_subset
mod_multiomics_integration <- function(harmonization_res, de_results = NULL,
                                        config, out_dir) {

    message("\n=== Multi-Omics Integration ===\n")

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    mae <- harmonization_res$mae
    cfg <- config$modes$multiomics$integration

    # Feature selection for integration
    feature_selection <- select_features_for_integration(
        mae = mae,
        de_results = de_results,
        config = config
    )

    # Subset MAE to selected features
    mae_subset <- subset_mae_by_features(mae, feature_selection$selected_features)

    # Run integration methods
    methods <- cfg$methods %||% c("DIABLO", "SNF")
    results <- list()

    # DIABLO
    if ("DIABLO" %in% methods) {
        message("\nRunning DIABLO integration...")
        diablo_dir <- file.path(out_dir, "diablo")
        dir.create(diablo_dir, showWarnings = FALSE)

        results$diablo <- tryCatch({
            run_diablo_integration(
                mae = mae_subset,
                config = config,
                out_dir = diablo_dir
            )
        }, error = function(e) {
            warning("DIABLO integration failed: ", e$message)
            NULL
        })

        if (!is.null(results$diablo)) {
            write_diablo_results(results$diablo, diablo_dir)
        }
    }

    # SNF
    if ("SNF" %in% methods) {
        message("\nRunning SNF integration...")
        snf_dir <- file.path(out_dir, "snf")
        dir.create(snf_dir, showWarnings = FALSE)

        results$snf <- tryCatch({
            run_snf_integration(
                mae = mae_subset,
                config = config,
                out_dir = snf_dir
            )
        }, error = function(e) {
            warning("SNF integration failed: ", e$message)
            NULL
        })

        if (!is.null(results$snf)) {
            write_snf_results(results$snf, snf_dir)
        }
    }

    # MOFA2
    if ("MOFA2" %in% methods) {
        message("\nRunning MOFA2 integration...")
        mofa_dir <- file.path(out_dir, "mofa")
        dir.create(mofa_dir, showWarnings = FALSE)

        results$mofa_results <- tryCatch({
            run_mofa2_integration(
                mae = mae_subset,
                config = config,
                out_dir = mofa_dir
            )
        }, error = function(e) {
            warning("MOFA2 integration failed: ", e$message)
            NULL
        })

        if (!is.null(results$mofa_results)) {
            write_mofa_results(results$mofa_results, mofa_dir)
        }
    }

    message("\nIntegration analysis complete")

    list(
        diablo_results = results$diablo,
        snf_results = results$snf,
        mofa_results = results$mofa_results,
        feature_selection = feature_selection,
        mae_subset = mae_subset,
        methods_run = methods
    )
}
