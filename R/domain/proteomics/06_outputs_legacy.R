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
        # Guarded on length(), not `%||%`: NULL[[1]] raises "subscript out of
        # bounds" in R, so a config with no stochastic imputation would abort
        # here rather than fall back to the preprocessing matrix.
        expr_model <- if (length(de_res$imputations) > 0) {
            de_res$imputations[[1]]
        } else {
            pre$expr_imp_single
        }

        final_results <- build_final_results_proteomics(
            pre = pre,
            summary_df = de_res$summary_df,
            contrasts_df = inputs$contrasts,
            row_data = pre$row_data,
            feature_id_col = config$modes$proteomics$de_table$id_col %||% "FeatureID",
            config = config,
            expr_model = expr_model
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
                                            feature_id_col = "FeatureID", config = NULL,
                                            expr_model = NULL) {
    # The exported .norm and Mean. columns must describe the SAME imputation
    # limma was fitted on. perseus_like is stochastic and the pipeline draws
    # twice -- once in preprocessing (expr_imp_single) and once per model run
    # via make_imputations_proteomics(), under different seeds. Exporting the
    # preprocessing draw made Mean.<num> - Mean.<den> disagree with the
    # reported log2FC on exactly the features that were imputed, by up to
    # 2.8 in log2. Pass the model's matrix; fall back only if unavailable.
    expr_for_model <- expr_model %||% pre$expr_imp_single
    pre_model <- pre
    pre_model$expr_imp_single <- expr_for_model

    cv_cols <- build_group_cv_proteomics(pre, contrasts_df, config)
    mean_cols <- build_group_mean_proteomics(pre_model, contrasts_df, config)

    # Pre-imputation estimate, computed on expr_filt (NAs still present). This
    # is the only genuinely model-free fold change in the table, and the one the
    # shrinkage check compares against.
    raw_stats <- build_group_raw_stats_proteomics(pre, contrasts_df, config)
    raw_log2fc <- compute_naive_log2fc_columns(
        raw_stats$means, contrasts_df, scale = "log2", prefix = "Mean.raw.")

    # Whether log2FC_from_means is worth emitting depends on how many times the
    # pipeline imputed, so it is decided from imputation$multi_imputation rather
    # than dropped outright:
    #
    #   multi_imputation: false -> one draw. log2FC.imputs is that draw's
    #     coefficient, and on a two-group design a difference of the model
    #     matrix's own group means IS that coefficient. The column would carry
    #     nothing log2FC.imputs does not, so it is omitted.
    #
    #   multi_imputation: true (default) -> log2FC.imputs is the consensus,
    #     log2(mean(2^logFC)) across no_repetitions runs (05_de_summary.R), while
    #     the Mean. columns come from imputations[[1]] alone. On a partially
    #     measured feature those are NOT the same number, and log2FC_from_means
    #     is what makes the gap visible instead of leaving it unstated.
    multi_imp <- config$modes$proteomics$imputation$multi_imputation %||% TRUE
    naive_log2fc <- if (isFALSE(multi_imp)) {
        NULL
    } else {
        compute_naive_log2fc_columns(mean_cols, contrasts_df, scale = "log2")
    }

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
        norm_expr = expr_for_model,
        mean_cols = mean_cols,
        naive_log2fc = naive_log2fc,
        raw_stat_cols = raw_stats$combined,
        raw_log2fc = raw_log2fc,
        # Keep log2FC.imputs adjacent to linearFC.imputs (a pinned contract) and
        # group the two model-free estimates after it.
        naive_after_fc = TRUE
    )
}

#' Pre-imputation per-group summaries for proteomics final results
#'
#' Computes \code{Mean.raw.<group>} and \code{N.observed.<group>} on
#' \code{expr_filt} — the filtered matrix with its NAs intact — rather than on
#' the imputed matrix the model was fitted on.
#'
#' \code{Mean.<group>} answers "what did limma see"; \code{Mean.raw.<group>}
#' answers "what was actually measured". They diverge exactly where imputation
#' did work, and \code{N.observed} is what tells a reader which rows those are.
#'
#' No CV counterpart is emitted: \code{build_group_cv_proteomics()} already
#' computes \code{CV.<group>} from \code{expr_filt}, so a \code{CV.raw.}
#' column would duplicate it exactly.
#'
#' A group with nothing measured yields \code{NA}, not \code{NaN}: see the
#' comment on \code{mean_measured_only} below.
#'
#' @param pre Proteomics preprocessing results (uses \code{expr_filt},
#'   \code{meta}).
#' @param contrasts_df Contrasts table (Factor, Numerator, Denominator).
#' @param config Full pipeline config (feature flag + sample-ID column).
#' @return list(means, combined) where \code{means} holds only the
#'   \code{Mean.raw.} columns (for the naive fold change) and \code{combined}
#'   holds all three blocks, or a list of NULLs when unavailable.
build_group_raw_stats_proteomics <- function(pre, contrasts_df, config = NULL) {
    empty <- list(means = NULL, combined = NULL)
    if (is.null(config)) return(empty)
    if (!isTRUE(config$modes$proteomics$excel$group_cv %||% TRUE)) return(empty)
    if (is.null(pre$expr_filt) || is.null(pre$meta)) return(empty)

    prot_cfg <- config$modes$proteomics %||% list()
    sample_id_col <- prot_cfg$effects$samples %||%
        prot_cfg$id_columns$sample_col %||% "SampleID"

    raw_log2 <- as.matrix(pre$expr_filt)

    # compute_group_stat_columns() directly, rather than
    # compute_group_mean_columns(), for the empty-group case only.
    # rowMeans(na.rm = TRUE) over a group with nothing measured returns NaN,
    # which Excel shows as an error value and which reads in a table as a
    # failed calculation rather than as "this group was never measured". This
    # block is computed on the matrix that still carries its NAs, so the
    # conversion belongs here rather than in the shared helper, whose
    # behaviour is deliberately left alone.
    mean_measured_only <- function(m) {
        mu <- rowMeans(m, na.rm = TRUE)
        mu[rowSums(!is.na(m)) == 0L] <- NA_real_
        mu
    }

    means <- compute_group_stat_columns(
        expr = raw_log2, sample_meta = pre$meta, sample_id_col = sample_id_col,
        contrasts_df = contrasts_df, stat_fn = mean_measured_only,
        prefix = "Mean.raw.")

    n_obs <- compute_group_observed_columns(
        expr = raw_log2, sample_meta = pre$meta, sample_id_col = sample_id_col,
        contrasts_df = contrasts_df)

    # No CV.raw. here: build_group_cv_proteomics() already computes CV. on
    # expr_filt, i.e. on observed values only. A CV.raw. column would be an
    # exact duplicate of CV. for this omic.

    blocks <- Filter(function(x) !is.null(x) && is.data.frame(x) && ncol(x) > 0,
                     list(means, n_obs))
    if (length(blocks) == 0L) return(empty)

    combined <- do.call(cbind, blocks)
    rownames(combined) <- rownames(blocks[[1]])
    list(means = means, combined = combined)
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

#' Per-run reconciliation of the multi-imputation fold change
#'
#' Multi-imputation proteomics reports \code{log2FC.imputs}, which is pooled
#' across independently imputed DE fits as
#' \code{log2( mean( 2^logFC ) )} — see \code{summarize_limma_mult_imputation()}.
#' That rule is invisible in the results table, so a reader who recomputes a
#' fold change from the exported per-sample values lands on a different number
#' and has no way to tell which is wrong. This table shows the pooling itself:
#' every per-run coefficient, its linear ratio, the arithmetic mean of those
#' ratios, and the reported value beside it.
#'
#' It reconciles the POOLING step only, and it does so exactly. Whether a
#' single run's coefficient equals a difference of that run's group means is a
#' separate question, answered by \code{log2FC_from_means} and
#' \code{log2FC_from_raw} in the results table.
#'
#' \code{jensen_gap} is the reason the two cannot be collapsed into one column.
#' The mean of a set of ratios is not the ratio implied by the mean of their
#' logs, so \code{log2FC.imputs} sits at or above \code{mean.log2FC.runs}, with
#' equality exactly when every run agrees — which is the case for a feature
#' that was measured in every sample and so had nothing imputed.
#'
#' Nothing here is a matrix: these are aggregates over fits that already
#' happened, and no model was fitted to any of these values.
#'
#' Built only when there are at least two runs AND those runs actually produced
#' different estimates. Several supported imputation methods return identical
#' matrices for every repetition, and for those the pooling is the identity: a
#' sheet of identical run columns would read as agreement between independent
#' draws when there were no separate draws to agree.
#'
#' @param de_res Proteomics DE result. Uses \code{runs_de_tables} (a list over
#'   imputation runs of per-contrast DE tables, whichever model
#'   \code{modes.proteomics.de.method} selected) and \code{summary_df}.
#' @param config Full pipeline config; only the feature ID column is read.
#' @return A data.frame with one row per feature and contrast, or \code{NULL}
#'   when there is nothing to reconcile. Rows are ordered by feature (in the
#'   order of the first run's table) and then by contrast.
build_de_reconciliation_proteomics <- function(de_res, config = NULL) {
    runs <- de_res$runs_de_tables
    summary_df <- de_res$summary_df

    # Two real runs is the whole precondition. One run means the pooling is the
    # identity and every delta below would be zero by construction, which reads
    # as a verification but demonstrates nothing. That covers both
    # imputation$multi_imputation: false (make_imputations_proteomics() draws
    # once) and precomputed DE input (load_precomputed_proteomics_de() wraps the
    # loaded tables as a single pseudo-run and has no imputations at all).
    if (is.null(runs) || length(runs) < 2L) return(NULL)
    if (is.null(summary_df) || nrow(summary_df) == 0L) return(NULL)

    id_col <- config$modes$proteomics$de_table$id_col %||% "FeatureID"
    contrasts <- names(runs[[1]])
    if (is.null(contrasts) || length(contrasts) == 0L) return(NULL)

    ref_tbl <- runs[[1]][[contrasts[1]]]
    if (is.null(ref_tbl) || !id_col %in% colnames(ref_tbl)) return(NULL)
    ref_ids <- as.character(ref_tbl[[id_col]])
    n_feat <- length(ref_ids)
    n_runs <- length(runs)
    if (n_feat == 0L) return(NULL)

    if (!id_col %in% colnames(summary_df)) return(NULL)
    srow <- match(ref_ids, as.character(summary_df[[id_col]]))

    blocks <- list()
    runs_vary <- FALSE
    for (cn in contrasts) {
        contrast_print <- normalize_contrast_name(cn)
        lr_col  <- paste0("linearRatio.imputs.", contrast_print)
        lfc_col <- paste0("log2FC.imputs.", contrast_print)
        if (!all(c(lr_col, lfc_col) %in% names(summary_df))) {
            message("    DE_reconciliation: no reported columns for contrast '",
                    contrast_print, "'; skipping it.")
            next
        }

        # Matched on the ID rather than on row position: the summary step has
        # already validated that the runs align, and matching keeps this honest
        # if that ever stops being true instead of silently pairing the wrong
        # features.
        per_run <- vapply(runs, function(one_run) {
            tbl <- one_run[[cn]]
            if (is.null(tbl) || !all(c(id_col, "logFC") %in% colnames(tbl))) {
                return(rep(NA_real_, n_feat))
            }
            as.numeric(tbl[["logFC"]])[match(ref_ids, as.character(tbl[[id_col]]))]
        }, numeric(n_feat))
        # vapply drops the dim when FUN.VALUE has length 1, so a single-feature
        # run would come back as a vector and rowMeans() would fail on it.
        per_run <- matrix(per_run, nrow = n_feat, ncol = n_runs)

        # Whether the runs are genuinely separate draws is asked of the numbers,
        # not of the config. Only perseus_like actually varies per run today:
        # impute_proteomics_qrilc() and impute_proteomics_dep2() call set.seed()
        # with a fixed configured seed INSIDE each call, overwriting the per-run
        # seed make_imputations_proteomics() sets, and none/minval/MinDet are
        # deterministic by design. All of those produce N identical matrices, so
        # the pooling is the identity and every column of this sheet would agree
        # by construction -- which reads as agreement between independent draws
        # when there were none. Asking the coefficients keeps this correct if a
        # method's seeding is ever fixed or a new one is added.
        #
        # identical(), not `!=`: any comparison against NA yields NA, and
        # dropping those with na.rm = TRUE reported genuinely different runs as
        # the same whenever run 1 was NA for a feature and another run was not.
        # The pooling above uses na.rm = TRUE, so a feature really can
        # contribute in some runs and not others, which makes that the live
        # case rather than a hypothetical one. Exact and deliberately without a
        # tolerance: the question is whether the estimates are identical,
        # including their NA pattern, and a tolerance would collapse small but
        # real differences between draws.
        if (!runs_vary) {
            runs_vary <- !all(vapply(
                seq_len(n_runs),
                function(r) identical(per_run[, r], per_run[, 1]),
                logical(1)
            ))
        }

        ratios <- 2^per_run
        # na.rm matches summarize_limma_mult_imputation(), so a run that failed
        # on a feature is dropped from both the reported value and this check.
        mean_ratio <- rowMeans(ratios, na.rm = TRUE)
        mean_log2fc_runs <- rowMeans(per_run, na.rm = TRUE)
        log2fc_from_mean_ratio <- log2(mean_ratio)

        reported_ratio <- as.numeric(summary_df[[lr_col]])[srow]
        reported_log2fc <- as.numeric(summary_df[[lfc_col]])[srow]

        block <- data.frame(ref_ids, contrast_print,
                            stringsAsFactors = FALSE, check.names = FALSE)
        names(block) <- c(id_col, "Contrast")

        for (r in seq_len(n_runs)) {
            block[[paste0("run", r, ".log2FC")]] <- per_run[, r]
            block[[paste0("run", r, ".ratio")]]  <- ratios[, r]
        }

        block[["mean.ratio"]]              <- mean_ratio
        block[["linearRatio.imputs"]]      <- reported_ratio
        block[["delta.linearRatio"]]       <- mean_ratio - reported_ratio
        block[["log2FC.from_mean_ratio"]]  <- log2fc_from_mean_ratio
        block[["log2FC.imputs"]]           <- reported_log2fc
        block[["delta.log2FC"]]            <- log2fc_from_mean_ratio - reported_log2fc
        block[["mean.log2FC.runs"]]        <- mean_log2fc_runs
        block[["jensen_gap"]]              <- reported_log2fc - mean_log2fc_runs

        block[[".feature_order"]] <- seq_len(n_feat)
        block[[".contrast_order"]] <- length(blocks) + 1L
        blocks[[length(blocks) + 1L]] <- block
    }

    if (length(blocks) == 0L) return(NULL)
    if (!runs_vary) {
        message("    DE_reconciliation: every imputation run produced the same ",
                "estimates, so there is no pooling to reconcile; skipping the sheet.")
        return(NULL)
    }

    out <- do.call(rbind, blocks)
    # Feature-major: every contrast for one feature sits together, which is how
    # the table gets read -- someone checking one protein by hand.
    out <- out[order(out[[".feature_order"]], out[[".contrast_order"]]), , drop = FALSE]
    out[[".feature_order"]] <- NULL
    out[[".contrast_order"]] <- NULL
    rownames(out) <- NULL
    out
}
