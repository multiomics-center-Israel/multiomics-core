#' Orchestrator: write all proteomics multi-imputation outputs (legacy-compatible)
write_proteomics_multimpute_outputs <- function(pre, de_res, inputs, config, out_dir, excel_order = NULL, write_runs = FALSE) {
    files <- character(0)
    dirs <- create_legacy_output_dirs(out_dir)

    # 1) datasets
    runs_for_datasets <- if (isTRUE(write_runs)) de_res$runs else NULL
    files <- c(files, write_proteomics_datasets_legacy(pre, runs_for_datasets, config, dirs))

    # 2) summary
    if (!is.null(de_res$summary_df)) {
        files <- c(files, write_limma_multimp_summary_legacy(de_res$summary_df, config, dirs))
    }

    # 3) wide limma per contrast
    if (!is.null(de_res$runs_de_tables) && length(de_res$runs_de_tables) > 0) {
        contrast_names <- names(de_res$runs_de_tables[[1]])
        for (cn in contrast_names) {
            files <- c(files, write_limma_results_multimp_legacy(de_res = de_res, contrast_name = cn, config = config, dirs = dirs))
        }
    }

    # 4) final results TSV
    if (!is.null(inputs$contrasts) && !is.null(de_res$summary_df)) {
        final_results <- build_final_results_proteomics(
            pre = pre,
            summary_df = de_res$summary_df,
            contrasts_df = inputs$contrasts,
            row_data = pre$row_data,
            feature_id_col = config$modes$proteomics$de_table$id_col %||% "FeatureID",
            config = config
        )
        files <- c(files, save_tsv(final_results, dirs$datasets, "final_results.tsv"))

        # 5) Excel outputs (delegated to dedicated file or assume available)
        if (exists("write_final_results_excels_proteomics")) {
            files <- c(files, write_final_results_excels_proteomics(
                final_results = final_results,
                pre = pre,
                config = config,
                out_dir = out_dir,
                excel_order = excel_order
            ))
        }
    }
    unique(files)
}

write_proteomics_datasets_legacy <- function(pre, runs = NULL, config, dirs) {
    cfg <- config$modes$proteomics
    id_col <- cfg$id_columns$protein_id %||% "Protein.Group"

    files <- character(0)

    # Helper: add ID column from rownames
    add_id <- function(x) {
        df <- as.data.frame(x, check.names = FALSE)
        df <- cbind(TEMP_ID = rownames(df), df)
        names(df)[1] <- id_col
        df
    }

    files <- c(files, save_tsv(add_id(pre$expr_filt), dirs$datasets, "protein_log2_filtered_unimputed.tsv"))

    fname_imp <- sprintf("protein_log2_filtered_imputed_once_width_%s_shift_%s.tsv", cfg$imputation$width, cfg$imputation$downshift)
    files <- c(files, save_tsv(add_id(pre$expr_imp_single), dirs$datasets, fname_imp))

    if (!is.null(runs) && length(runs) > 0) {
        rep_dir <- file.path(dirs$datasets, "imputed_repetitions")
        ensure_dir(rep_dir)
        for (i in seq_along(runs)) {
            expr_i <- runs[[i]]$expr_imp
            if (is.null(expr_i)) next
            files <- c(files, save_tsv(add_id(expr_i), rep_dir, sprintf("protein_log2_filtered_imputed_%02d.tsv", i)))
        }
    }
    unique(files)
}

write_limma_multimp_summary_legacy <- function(summary_df, config, dirs) {
    save_tsv(summary_df, dirs$datasets, sprintf("limma_multimp_summary_p%s.tsv",
                                                p_tag_generic(config, "proteomics")))
}

write_limma_results_multimp_legacy <- function(de_res, contrast_name, config, dirs) {
    feature_id_col <- config$modes$proteomics$de_table$id_col %||% "FeatureID"
    wide_df <- build_limma_results_multimp_wide(
        runs_de_tables = de_res$runs_de_tables,
        contrast_name = contrast_name,
        feature_id_col = feature_id_col
    )
    fname <- sprintf("limma_results_multimp_p%s.tsv", p_tag_generic(config, "proteomics"))
    save_tsv(wide_df, dirs$datasets, fname)
}

build_limma_results_multimp_wide <- function(runs_de_tables, contrast_name, stats_cols = c("logFC", "P.Value", "adj.P.Val"), feature_id_col = "FeatureID") {
    stopifnot(length(runs_de_tables) >= 1)
    base <- runs_de_tables[[1]][[contrast_name]]
    if (is.null(base)) stop(sprintf("Contrast '%s' missing in imputation run 1", contrast_name))

    if (!feature_id_col %in% colnames(base)) {
        stop(sprintf("Feature ID column '%s' missing in base table", feature_id_col))
    }

    id_cols <- intersect(c(feature_id_col, "Protein.Names", "Genes", "First.Protein.Description", "Contrast"), colnames(base))
    out <- base[, id_cols, drop = FALSE]

    for (i in seq_along(runs_de_tables)) {
        tab <- runs_de_tables[[i]][[contrast_name]]
        if (is.null(tab)) stop(sprintf("Critical: Contrast '%s' is missing in imputation run %d", contrast_name, i))

        tab <- align_de_table_by_feature_id(tab = tab, ref_ids = out[[feature_id_col]], run_i = i, contrast_name = contrast_name, id_col = feature_id_col)
        stat_block <- tab[, intersect(stats_cols, colnames(tab)), drop = FALSE]
        colnames(stat_block) <- paste0(colnames(stat_block), ".", i)
        out <- cbind(out, stat_block)
    }
    out
}

#' Build Final Results Table for Proteomics
#' 
#' This function aggregates differential expression results with expression data
#' and annotations. It ensures Z-score calculation by passing the correct mode.
#'
#' @param pre List containing preprocessed data (including expr_filt and expr_imp_single).
#' @param summary_df Dataframe containing DE summary statistics.
#' @param contrasts_df Dataframe defining the experimental contrasts.
#' @param row_data Optional annotation data (defaults to pre$row_data).
#' @param feature_id_col The column name for unique identifiers (default "FeatureID").
#'
#' @return A consolidated dataframe with statistics, expression values, and Z-scores.
build_final_results_proteomics <- function(pre, summary_df, contrasts_df, row_data = NULL,
                                            feature_id_col = "FeatureID", config = NULL) {
    cv_cols <- build_group_cv_proteomics(pre, contrasts_df, config)
    mean_cols <- build_group_mean_proteomics(pre, contrasts_df, config)

    build_final_results_generic(
        summary_df = summary_df,
        expr_df = pre$expr_filt,
        contrasts_df = contrasts_df,
        feature_id_col = feature_id_col,
        annot_cols = c(
            "Protein.Names" = "Protein.Names",
            "Genes" = "Genes",
            "First.Protein.Description" = "First.Protein.Description"
        ),
        row_data = row_data %||% pre$row_data,
        fc_is_signed = TRUE, # linearFC is signed
        cv_cols = cv_cols,
        # expr_filt above is what was measured (NAs for unobserved); this is the
        # imputed matrix limma was fitted on. Both are needed to walk from the
        # per-sample values to the reported logFC.
        norm_expr = pre$expr_imp_single,
        mean_cols = mean_cols
    )
}

#' Build per-group mean columns for proteomics final results
#'
#' Means are taken on the imputed log2 matrix (\code{expr_imp_single}) — the
#' scale limma was fitted on — so that
#' \code{Mean.<numerator> - Mean.<denominator>} lands close to the reported
#' \code{log2FC.imputs}. It will not match exactly: the reported statistic
#' averages over all imputation runs while this block shows a single one, and
#' limma reports a moderated coefficient rather than a difference of means.
#'
#' @param pre Proteomics preprocessing results (uses \code{expr_imp_single},
#'   \code{meta}).
#' @param contrasts_df Contrasts table (Factor, Numerator, Denominator).
#' @param config Full pipeline config (feature flag + sample-ID column).
#' @return Feature-indexed data.frame of \code{Mean.<group>} columns, or NULL.
build_group_mean_proteomics <- function(pre, contrasts_df, config = NULL) {
    if (is.null(config)) return(NULL)
    # Same switch as the CV block: one flag governs the whole summary section.
    if (!isTRUE(config$modes$proteomics$excel$group_cv %||% TRUE)) return(NULL)
    if (is.null(pre$expr_imp_single) || is.null(pre$meta)) return(NULL)

    prot_cfg <- config$modes$proteomics %||% list()
    sample_id_col <- prot_cfg$effects$samples %||%
        prot_cfg$id_columns$sample_col %||% "SampleID"

    compute_group_mean_columns(
        expr          = as.matrix(pre$expr_imp_single),
        sample_meta   = pre$meta,
        sample_id_col = sample_id_col,
        contrasts_df  = contrasts_df
    )
}

#' Resolve the log2 pseudocount offset applied to the proteomics assay
#'
#' The pipeline applies \code{log2(x + 1)} only for preprocessed input declared
#' on a linear scale (\code{R/domain/proteomics/01_expression.R}); the DIA-NN
#' path and any already-log2 input use plain \code{log2(x)}. The exact inverse
#' is therefore \code{2^x - offset}.
#'
#' @param config Full pipeline config.
#' @return Numeric offset: 1 for preprocessed+linear input, else 0.
proteomics_log_offset <- function(config) {
    cfg <- config$modes$proteomics %||% list()
    is_preprocessed <- identical(cfg$input$format, "preprocessed")
    scale_in <- cfg$scale_in %||%
        (if (isTRUE(cfg$files$is_logtransformed)) "log2" else "linear")
    if (is_preprocessed && identical(scale_in, "linear")) 1 else 0
}

#' Build per-group CV columns for proteomics final results
#'
#' CV is computed on linear intensities, back-transformed from the unimputed
#' log2 matrix (\code{2^expr_filt - offset}). Imputed values are intentionally
#' excluded — \code{expr_filt} carries NAs for unobserved measurements, so CV
#' uses observed values only and a group with <2 observations yields NA.
#'
#' @param pre Proteomics preprocessing results (uses \code{expr_filt}, \code{meta}).
#' @param contrasts_df Contrasts table (Factor, Numerator, Denominator).
#' @param config Full pipeline config (feature flag, sample-ID column, log scale).
#' @return Feature-indexed data.frame of \code{CV.<group>} columns, or NULL.
build_group_cv_proteomics <- function(pre, contrasts_df, config = NULL) {
    if (is.null(config)) return(NULL)
    enabled <- config$modes$proteomics$excel$group_cv %||% TRUE
    if (!isTRUE(enabled)) return(NULL)
    if (is.null(pre$expr_filt) || is.null(pre$meta)) return(NULL)

    prot_cfg <- config$modes$proteomics %||% list()
    sample_id_col <- prot_cfg$effects$samples %||%
        prot_cfg$id_columns$sample_col %||% "SampleID"

    offset <- proteomics_log_offset(config)
    expr_linear <- 2^as.matrix(pre$expr_filt) - offset  # NAs (unobserved) preserved

    compute_group_cv_columns(
        expr_linear   = expr_linear,
        sample_meta   = pre$meta,
        sample_id_col = sample_id_col,
        contrasts_df  = contrasts_df
    )
}
