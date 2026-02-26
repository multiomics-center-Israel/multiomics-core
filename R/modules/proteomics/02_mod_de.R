#' Auto-generate all pairwise contrasts from metadata
#' @param meta    data.frame with sample metadata
#' @param group_col character, column name defining biological groups
#' @return data.frame with Contrast_name, Factor, Numerator, Denominator
auto_generate_contrasts <- function(meta, group_col) {
    lvls <- sort(unique(as.character(meta[[group_col]])))
    lvls <- lvls[!is.na(lvls) & nzchar(lvls)]
    if (length(lvls) < 2) stop("Cannot auto-generate contrasts: fewer than 2 groups in '", group_col, "'.")
    pairs <- combn(lvls, 2)
    df <- data.frame(
        Contrast_name = apply(pairs, 2, function(p) paste0(p[1], "_vs_", p[2])),
        Factor        = group_col,
        Numerator     = pairs[1, ],
        Denominator   = pairs[2, ],
        stringsAsFactors = FALSE
    )
    message(sprintf("Auto-generated %d pairwise contrast(s) from '%s': %s",
                    nrow(df), group_col, paste(df$Contrast_name, collapse = ", ")))
    df
}

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

    # Auto-generate contrasts if not provided
    contrasts_df <- inputs$contrasts
    if (is.null(contrasts_df)) {
        group_col <- cfg$effects$color %||% "Condition"
        contrasts_df <- auto_generate_contrasts(pre$meta, group_col)
    }

    # ---- choose DE method ----
    method <- cfg$de$method %||% "limma"

    # --- Shared step: generate multiple imputations ---
    imputations <- make_imputations_proteomics(
        expr_mat = pre$expr_filt,
        cfg      = config,
        verbose  = verbose
    )

    if (identical(method, "limma")) {
        # Run limma on imputations
        runs <- run_limma_multimp(
            imputations  = imputations,
            meta         = pre$meta,
            contrasts_df = contrasts_df,
            prot_tbl     = pre$row_data,
            cfg          = config,
            verbose      = verbose
        )

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

    if (identical(method, "ttest") || identical(method, "welch")) {
        var_equal <- identical(method, "ttest")
        runs <- lapply(seq_along(imputations), function(i) {
            if (isTRUE(verbose)) message(sprintf("%s on imputation: %d / %d", method, i, length(imputations)))
            run_ttest_de(
                expr_imp     = imputations[[i]],
                meta         = pre$meta,
                contrasts_df = contrasts_df,
                prot_tbl     = pre$row_data,
                cfg          = config,
                var_equal    = var_equal
            )
        })

        runs_de_tables <- lapply(runs, function(x) x$de_tables)

        summary_df <- summarize_limma_mult_imputation(
            runs_de_tables = runs_de_tables,
            config         = config
        )

        return(list(
            method = method,
            imputations = imputations,
            runs = runs,
            runs_de_tables = runs_de_tables,
            summary_df = summary_df,
            de_model = NULL
        ))
    }

    if (identical(method, "anova")) {
        runs <- lapply(seq_along(imputations), function(i) {
            if (isTRUE(verbose)) message(sprintf("ANOVA on imputation: %d / %d", i, length(imputations)))
            run_anova_posthoc(
                expr_imp     = imputations[[i]],
                meta         = pre$meta,
                contrasts_df = contrasts_df,
                prot_tbl     = pre$row_data,
                cfg          = config
            )
        })

        runs_de_tables <- lapply(runs, function(x) x$de_tables)

        summary_df <- summarize_limma_mult_imputation(
            runs_de_tables = runs_de_tables,
            config         = config
        )

        return(list(
            method = "anova",
            imputations = imputations,
            runs = runs,
            runs_de_tables = runs_de_tables,
            summary_df = summary_df,
            de_model = NULL
        ))
    }

    stop("Unsupported proteomics DE method: ", method)
}
