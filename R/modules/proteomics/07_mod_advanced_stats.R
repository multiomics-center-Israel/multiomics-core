#' Proteomics advanced statistics module
#'
#' Bridges DE results to run_advanced_stats() by extracting
#' per-contrast tables.
#'
#' @param de_res  DE results list (with summary_df)
#' @param pre     Preprocessed proteomics object
#' @param inputs  Loaded input data
#' @param config  Full config list
#' @param out_dir Output root directory
#' @return Advanced stats results list or NULL
mod_proteomics_advanced_stats <- function(de_res, pre, inputs, config, out_dir) {
    cfg <- config$modes$proteomics
    if (!isTRUE(cfg$advanced_stats$enabled)) {
        message("Advanced stats disabled.")
        return(NULL)
    }

    summary_df <- de_res$summary_df
    if (is.null(summary_df)) {
        message("No summary_df available. Skipping advanced stats.")
        return(NULL)
    }

    # Build per-contrast DE tables
    contrasts <- sub("^padj\\.imputs\\.", "",
                     grep("^padj\\.imputs\\.", colnames(summary_df), value = TRUE))

    de_tables <- lapply(setNames(contrasts, contrasts), function(cn) {
        extract_de_table_for_pathway(summary_df, cn, config)
    })

    # These analyses fit their own models on the QC matrix: one imputation draw,
    # batch-corrected when a correction ran, with no batch term. That is not the
    # batch-adjusted inference DE performs, and until they model batch
    # themselves (tracked separately) the run says so instead of presenting
    # them as equivalent. The report carries the same caveat.
    if (isTRUE(pre$batch_corrected)) {
        warning("Advanced stats are computed on the batch-corrected single-imputation ",
                "matrix without a batch term. Their intervals and p-values are not ",
                "batch-adjusted inference and can be overconfident; use the DE results ",
                "for inference.", call. = FALSE)
    }

    run_advanced_stats(
        da_results      = de_tables,
        normalized_data = pre$expr_imp_single,
        metadata        = pre$meta,
        config          = config,
        output_dir      = file.path(out_dir, "advanced_stats")
    )
}
