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

#' Generate one contrast per non-control level, each against the control
#'
#' Used by the "limma_percontrast" DE method: every other level of `group_col` is
#' contrasted against `control` (the denominator), so ratios read as test/control.
#'
#' @param meta      data.frame with sample metadata
#' @param group_col character, column name defining biological groups
#' @param control   character, the control/reference level (used as denominator)
#' @return data.frame with Contrast_name, Factor, Numerator, Denominator
auto_generate_control_contrasts <- function(meta, group_col, control) {
    lvls <- sort(unique(as.character(meta[[group_col]])))
    lvls <- lvls[!is.na(lvls) & nzchar(lvls)]
    if (!control %in% lvls) {
        stop("control_condition '", control, "' is not a level of '", group_col,
             "'. Available levels: ", paste(lvls, collapse = ", "))
    }
    test_lvls <- setdiff(lvls, control)
    if (length(test_lvls) < 1) {
        stop("No non-control levels in '", group_col,
             "' to contrast against control '", control, "'.")
    }
    df <- data.frame(
        Contrast_name = paste0(test_lvls, "_vs_", control),
        Factor        = group_col,
        Numerator     = test_lvls,
        Denominator   = control,
        stringsAsFactors = FALSE
    )
    message(sprintf("Generated %d control-referenced contrast(s) vs '%s': %s",
                    nrow(df), control, paste(df$Contrast_name, collapse = ", ")))
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

    # Check for pre-computed DE tables (skip limma if provided)
    de_table_files <- cfg$files$de_table
    has_precomputed <- !is.null(de_table_files) && length(de_table_files) > 0 &&
                       all(nzchar(unlist(de_table_files)))

    if (has_precomputed) {
        message("Proteomics DE: loading pre-computed DE tables (skipping limma)")
        return(load_precomputed_proteomics_de(config, contrasts_df = inputs$contrasts))
    }

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

    if (identical(method, "limma_percontrast")) {
        # Each contrast is tested on its own two-group matrix (control + one test
        # condition). When contrasts are not given explicitly, build them as
        # each-non-control vs control from cfg$de$control_condition.
        if (is.null(inputs$contrasts)) {
            control_cond <- cfg$de$control_condition
            if (is.null(control_cond) || !nzchar(control_cond)) {
                stop("de.method 'limma_percontrast' requires de.control_condition ",
                     "(the control/reference level) when contrasts are not provided.")
            }
            group_col <- cfg$de_table$group_col %||% cfg$effects$color %||% "Condition"
            contrasts_df <- auto_generate_control_contrasts(pre$meta, group_col, control_cond)
        }

        # Observed (pre-imputation) mask drives the per-contrast protein filter.
        observed <- !is.na(pre$expr_filt)

        runs <- lapply(seq_along(imputations), function(i) {
            if (isTRUE(verbose)) message(sprintf("limma_percontrast on imputation: %d / %d", i, length(imputations)))
            run_limma_percontrast_proteomics(
                expr_imp     = imputations[[i]],
                observed     = observed,
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
            method = "limma_percontrast",
            imputations = imputations,
            runs = runs,
            runs_de_tables = runs_de_tables,
            summary_df = summary_df,
            de_model = NULL
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
