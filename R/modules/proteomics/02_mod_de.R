#' Proteomics DE module (orchestration only)
#'
#' Runs DE according to cfg$de$method (currently supports "limma").
#' No file I/O here.
#'
#' @param pre     output of preprocess_proteomics()
#' @param inputs  output of load_proteomics_inputs()
#' @param config  full config
#' @param verbose logical
#' @return list(method, runs, runs_de_tables, summary_df, de_model)
mod_proteomics_de <- function(pre, inputs, config, verbose = FALSE) {
    assert_pre_contract(pre, stage = "proteomics")

    cfg <- config$modes$proteomics

    # Check for pre-computed DE tables (skip limma if provided)
    de_table_files <- cfg$files$de_table
    has_precomputed <- !is.null(de_table_files) && length(de_table_files) > 0 &&
                       all(nzchar(unlist(de_table_files)))

    if (has_precomputed) {
        message("Proteomics DE: loading pre-computed DE tables (skipping limma)")
        return(load_precomputed_proteomics_de(config, contrasts_df = inputs$contrasts))
    }

    # ---- choose DE method ----
    method <- cfg$de$method %||% "limma"

    if (identical(method, "limma")) {
        # 1) imputations
        imputations <- make_imputations_proteomics(
            expr_mat = pre$expr_filt,
            cfg      = config,
            verbose  = verbose
        )

        # 2) run limma on imputations
        # Use pre$row_data (with custom annotation + contaminant filtering applied)
        # instead of raw inputs$protein so gene names propagate to DE tables
        runs <- run_limma_multimp(
            imputations  = imputations,
            meta         = pre$meta,
            contrasts_df = inputs$contrasts,
            prot_tbl     = pre$row_data,
            cfg          = config,
            verbose      = verbose
        )

        # 3) collect DE tables + summarize
        runs_de_tables <- lapply(runs, function(x) x$de_tables)

        summary_df <- summarize_limma_mult_imputation(
            runs_de_tables = runs_de_tables,
            config         = config
        )

        return(list(
            method = "limma",
            imputations = imputations,
            runs = runs,
            runs_de_tables = runs_de_tables,
            summary_df = summary_df,
            de_model = runs[[1]]$fit2
        ))
    }

    stop("Unsupported proteomics DE method: ", method)
}
