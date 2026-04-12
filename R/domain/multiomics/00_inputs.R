#' Load preprocessed data from individual omics pipelines
#'
#' Expects that rna_pre, prot_pre, and/or metab_pre targets have already run.
#' Returns a named list with available omics data.
#'
#' @param config Full config object
#' @param rna_pre Preprocessed RNA-seq data (from pipe_rnaseq)
#' @param prot_pre Preprocessed proteomics data (from pipe_proteomics)
#' @param metab_pre Preprocessed metabolomics data (from pipe_metabolomics)
#' @return Named list with: omics_available, rna, proteomics, metabolomics
load_multiomics_inputs <- function(config, rna_pre = NULL, prot_pre = NULL, metab_pre = NULL) {

    omics_available <- character(0)
    result <- list()

    # Collect available omics
    if (!is.null(rna_pre)) {
        omics_available <- c(omics_available, "transcriptomics")
        result$transcriptomics <- rna_pre
        message(sprintf("RNA-seq data loaded: %d features x %d samples",
                        nrow(rna_pre$expr_work), ncol(rna_pre$expr_work)))
    }

    if (!is.null(prot_pre)) {
        omics_available <- c(omics_available, "proteomics")
        result$proteomics <- prot_pre
        message(sprintf("Proteomics data loaded: %d features x %d samples",
                        nrow(prot_pre$expr_work), ncol(prot_pre$expr_work)))
    }

    if (!is.null(metab_pre)) {
        omics_available <- c(omics_available, "metabolomics")
        result$metabolomics <- metab_pre
        message(sprintf("Metabolomics data loaded: %d features x %d samples",
                        nrow(metab_pre$expr_work), ncol(metab_pre$expr_work)))
    }

    if (length(omics_available) < 2) {
        stop(
            "Multi-omics integration requires at least 2 omics layers. ",
            "Found: ", paste(omics_available, collapse = ", "),
            call. = FALSE
        )
    }

    result$omics_available <- omics_available
    result$n_omics <- length(omics_available)

    message(sprintf("Multi-omics integration ready: %s",
                    paste(omics_available, collapse = " + ")))

    result
}


#' Validate multi-omics inputs before harmonization
#'
#' Checks that all omics have compatible sample metadata and sufficient overlap.
#'
#' @param inputs Output from load_multiomics_inputs()
#' @param config Full config object
#' @return Invisibly returns TRUE if valid, stops otherwise
validate_multiomics_inputs <- function(inputs, config) {

    omics <- inputs$omics_available

    # Extract sample IDs from each omics
    sample_lists <- list()
    for (om in omics) {
        meta <- inputs[[om]]$meta

        # Map internal omics names to config keys
        om_config_key <- switch(om,
                                "transcriptomics" = "rna",
                                "proteomics" = "proteomics",
                                "metabolomics" = "metabolomics",
                                om)  # fallback

        sample_col <- config$modes[[om_config_key]]$effects$samples %||%
                      config$modes[[om_config_key]]$id_columns$sample_col %||%
                      "SampleID"

        if (!sample_col %in% colnames(meta)) {
            stop(sprintf("Sample column '%s' not found in %s metadata",
                         sample_col, om))
        }

        sample_lists[[om]] <- as.character(meta[[sample_col]])
    }

    # Check sample overlap
    common_samples <- Reduce(intersect, sample_lists)

    if (length(common_samples) < 3) {
        stop(
            "Insufficient sample overlap across omics layers. ",
            "Common samples: ", length(common_samples), " (need ≥3)\n",
            paste(
                sprintf("  %s: %d samples", names(sample_lists),
                        sapply(sample_lists, length)),
                collapse = "\n"
            ),
            call. = FALSE
        )
    }

    message(sprintf(
        "Sample validation passed: %d common samples across %d omics layers",
        length(common_samples), length(omics)
    ))

    invisible(TRUE)
}
