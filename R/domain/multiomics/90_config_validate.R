#' Validate multi-omics configuration
#'
#' Checks that multi-omics config contains all required fields and that
#' at least 2 omics modes are configured.
#'
#' @param multiomics_cfg Config object at config$modes$multiomics
#' @return Invisibly returns TRUE if valid, stops with error otherwise
validate_multiomics_config <- function(multiomics_cfg) {

    if (is.null(multiomics_cfg)) {
        stop("config$modes$multiomics is NULL", call. = FALSE)
    }

    # Check for required integration config
    if (is.null(multiomics_cfg$integration)) {
        stop(
            "config$modes$multiomics$integration is required. ",
            "Must specify integration methods (e.g., DIABLO, SNF)",
            call. = FALSE
        )
    }

    # Check that at least one integration method is enabled
    methods <- multiomics_cfg$integration$methods
    if (is.null(methods) || length(methods) == 0) {
        warning(
            "No integration methods specified in config$modes$multiomics$integration$methods. ",
            "Defaulting to: DIABLO, SNF"
        )
        multiomics_cfg$integration$methods <- c("DIABLO", "SNF")
    }

    # Validate method names
    valid_methods <- c("DIABLO", "SNF", "MOFA2")  # MOFA2 deferred for now
    invalid_methods <- setdiff(methods, valid_methods)
    if (length(invalid_methods) > 0) {
        stop(
            "Invalid integration method(s): ", paste(invalid_methods, collapse = ", "), "\n",
            "Valid methods: ", paste(valid_methods, collapse = ", "),
            call. = FALSE
        )
    }

    # Validate DIABLO config if specified
    if ("DIABLO" %in% methods) {
        diablo_cfg <- multiomics_cfg$integration$diablo %||% list()

        if (is.null(diablo_cfg$ncomp)) {
            message("  DIABLO ncomp not specified, defaulting to 2")
            multiomics_cfg$integration$diablo$ncomp <- 2
            diablo_cfg$ncomp <- 2
        }

        if (diablo_cfg$ncomp < 1 || diablo_cfg$ncomp > 10) {
            warning("DIABLO ncomp should be between 1 and 10; got ", diablo_cfg$ncomp)
        }
    }

    # Validate SNF config if specified
    if ("SNF" %in% methods) {
        snf_cfg <- multiomics_cfg$integration$snf %||% list()

        if (is.null(snf_cfg$K)) {
            message("  SNF K (neighbors) not specified, defaulting to 20")
            multiomics_cfg$integration$snf$K <- 20
        }

        if (is.null(snf_cfg$n_clusters)) {
            message("  SNF n_clusters not specified, defaulting to 2")
            multiomics_cfg$integration$snf$n_clusters <- 2
        }
    }

    # Validate MOFA2 config if specified
    if ("MOFA2" %in% methods) {
        mofa2_cfg <- multiomics_cfg$integration$mofa2 %||% list()

        if (is.null(mofa2_cfg$num_factors)) {
            message("  MOFA2 num_factors not specified, defaulting to 10")
            multiomics_cfg$integration$mofa2$num_factors <- 10
        } else if (mofa2_cfg$num_factors < 1 || mofa2_cfg$num_factors > 50) {
            warning("MOFA2 num_factors should be between 1 and 50; got ", mofa2_cfg$num_factors)
        }

        valid_conv_modes <- c("fast", "medium", "slow")
        if (!is.null(mofa2_cfg$convergence_mode) &&
            !(mofa2_cfg$convergence_mode %in% valid_conv_modes)) {
            stop("MOFA2 convergence_mode must be one of: ", paste(valid_conv_modes, collapse = ", "))
        }

        if (is.null(mofa2_cfg$seed)) {
            multiomics_cfg$integration$mofa2$seed <- 42
        }
    }

    # Validate foundational config
    found_cfg <- multiomics_cfg$foundational %||% list()
    if (is.null(found_cfg$run_foundational)) {
        multiomics_cfg$foundational$run_foundational <- TRUE
    }
    if (is.null(found_cfg$correlation_method)) {
        multiomics_cfg$foundational$correlation_method <- "spearman"
    }
    if (is.null(found_cfg$top_variable_features)) {
        multiomics_cfg$foundational$top_variable_features <- 500
    }

    # Validate consensus config
    cons_cfg <- multiomics_cfg$consensus %||% list()
    if (is.null(cons_cfg$compare_methods)) {
        multiomics_cfg$consensus$compare_methods <- TRUE
    }

    # Validate stability config
    stab_cfg <- multiomics_cfg$stability %||% list()
    if (is.null(stab_cfg$run_stability)) {
        multiomics_cfg$stability$run_stability <- FALSE
    }
    if (is.null(stab_cfg$n_bootstrap)) {
        multiomics_cfg$stability$n_bootstrap <- 100
    }

    # Validate enrichment multigsea config
    enrich_cfg <- multiomics_cfg$enrichment %||% list()
    mg_cfg <- enrich_cfg$multigsea %||% list()
    if (is.null(mg_cfg$run_multigsea)) {
        multiomics_cfg$enrichment$multigsea$run_multigsea <- TRUE
    }
    if (is.null(mg_cfg$pvalue_threshold)) {
        multiomics_cfg$enrichment$multigsea$pvalue_threshold <- 0.05
    }

    # Validate pathview config
    pv_cfg <- enrich_cfg$pathview %||% list()
    if (is.null(pv_cfg$run_pathview)) {
        multiomics_cfg$enrichment$pathview$run_pathview <- TRUE
    }
    if (is.null(pv_cfg$top_n)) {
        multiomics_cfg$enrichment$pathview$top_n <- 5
    }

    # Validate feature selection config
    feat_cfg <- multiomics_cfg$feature_selection %||% list()
    if (is.null(feat_cfg$method)) {
        message("  Feature selection method not specified, defaulting to 'variance'")
        multiomics_cfg$feature_selection$method <- "variance"
    }

    if (is.null(feat_cfg$top_n)) {
        message("  Feature selection top_n not specified, defaulting to 500")
        multiomics_cfg$feature_selection$top_n <- 500
    }

    # Validate per-omics feature counts (if provided)
    if (!is.null(feat_cfg$per_omics)) {
        if (!is.list(feat_cfg$per_omics)) {
            stop("feature_selection$per_omics must be a named list (omics_name: count)", call. = FALSE)
        }
        for (nm in names(feat_cfg$per_omics)) {
            val <- feat_cfg$per_omics[[nm]]
            if (!is.numeric(val) || val < 1) {
                stop(sprintf("feature_selection$per_omics$%s must be a positive number, got: %s", nm, val),
                     call. = FALSE)
            }
        }
        message(sprintf("  Per-omics feature counts: %s",
                paste(sprintf("%s=%d", names(feat_cfg$per_omics),
                              unlist(feat_cfg$per_omics)), collapse = ", ")))
    }

    # Check for condition column (needed for supervised methods like DIABLO)
    if ("DIABLO" %in% methods && is.null(multiomics_cfg$condition_column)) {
        warning(
            "DIABLO requires config$modes$multiomics$condition_column. ",
            "Will attempt to use config$design$condition_column as fallback."
        )
    }

    message("Multi-omics config validation passed")
    invisible(multiomics_cfg)
}
