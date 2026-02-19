# R/domain/metabolomics/03_differential.R
#
# Differential expression analysis for metabolomics data.
# Supports limma, Welch t-test, and Wilcoxon rank-sum methods.
#
# Operates on the standard pre-processing contract (expr_work, meta).
# Returns a wide-format summary_df compatible with the DE contract.
#
# Reuses: assert_numeric_matrix, assert_one_of, %||%


# ---- public entry point -----------------------------------------------------

#' Run metabolomics differential analysis
#'
#' @param pre    List from preprocess_metabolomics() (pre contract).
#' @param config Full pipeline config.
#' @return list conforming to the DE contract:
#'   summary_df, method, de_tables (per-contrast list), de_model
run_metabolomics_de <- function(pre, config) {
    cfg <- config$modes$metabolomics
    de_cfg <- cfg$de %||% list()

    method <- tolower(de_cfg$method %||% "limma")
    assert_one_of(method, "de$method", c("limma", "t_test", "wilcoxon"))

    condition_col <- de_cfg$condition_column %||% cfg$effects$color %||% "sample_type"
    sample_col <- cfg$effects$samples %||% "sample_id"

    contrasts <- de_cfg$contrasts
    if (is.null(contrasts) || length(contrasts) == 0) {
        stop("metabolomics DE: config$modes$metabolomics$de$contrasts is required.")
    }
    if (is.list(contrasts)) contrasts <- unlist(contrasts)

    mat  <- pre$expr_work
    # Use pre-scaling (log-transformed) matrix for DE statistical tests.
    # Autoscaling (mean-center + divide by SD) standardises every feature to
    # mean=0, SD=1 which distorts within-group variance and produces
    # uniform p-values.  Autoscaled data is for multivariate methods only.
    mat_for_test <- pre$expr_log %||% mat
    meta <- pre$meta
    assert_numeric_matrix(mat_for_test, "metab_expr_for_test")

    # Align metadata to matrix columns
    meta <- meta[match(colnames(mat_for_test), meta[[sample_col]]), , drop = FALSE]
    condition <- factor(meta[[condition_col]])

    # Thresholds for significance flags
    padj_cutoff <- de_cfg$p_cutoff %||% 0.05
    linear_fc   <- de_cfg$linear_fc_cutoff %||% 1.5
    log2fc_cut  <- log2(linear_fc)

    # Run DE for each contrast
    de_tables <- list()
    de_model  <- NULL

    for (ctr in contrasts) {
        ctr_label <- make_contrast_label(ctr)
        message("metabolomics DE [", method, "]: ", ctr)

        tbl <- switch(method,
            limma    = de_limma(mat_for_test, condition, ctr),
            t_test   = de_t_test(mat_for_test, condition, ctr),
            wilcoxon = de_wilcoxon(mat_for_test, condition, ctr)
        )

        # Capture limma model from first contrast
        if (method == "limma" && is.null(de_model)) {
            de_model <- attr(tbl, "fit")
        }

        de_tables[[ctr_label]] <- tbl
    }

    # Build wide summary_df
    summary_df <- build_de_summary(de_tables, padj_cutoff, log2fc_cut)

    message("metabolomics DE complete: ", nrow(summary_df), " features, ",
            sum(summary_df$pass_any_contrast == 1, na.rm = TRUE), " significant")

    list(
        summary_df = summary_df,
        method     = method,
        de_tables  = de_tables,
        de_model   = de_model
    )
}


# ---- limma -----------------------------------------------------------------

#' Run limma on a single contrast
#'
#' @param mat       Numeric matrix (features x samples).
#' @param condition Factor of conditions.
#' @param contrast_str  Character, e.g. "B - A".
#' @return data.frame with feature_id, logFC, AveExpr, statistic, P.Value, adj.P.Val
de_limma <- function(mat, condition, contrast_str, mat_for_fc = NULL) {
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("Package 'limma' is required for limma DE.")
    }

    # Impute NAs with row means (limma cannot handle NA)
    mat_imp <- mat
    for (i in seq_len(nrow(mat_imp))) {
        nas <- is.na(mat_imp[i, ])
        if (any(nas)) mat_imp[i, nas] <- mean(mat_imp[i, !nas], na.rm = TRUE)
    }

    design <- stats::model.matrix(~ 0 + condition)
    colnames(design) <- levels(condition)

    fit <- limma::lmFit(mat_imp, design)
    contrast_matrix <- limma::makeContrasts(contrasts = contrast_str,
                                             levels = design)
    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
    fit2 <- limma::eBayes(fit2)

    tt <- limma::topTable(fit2, number = Inf, sort.by = "none")

    res <- data.frame(
        feature_id = rownames(tt),
        logFC      = tt$logFC,
        AveExpr    = tt$AveExpr,
        statistic  = tt$t,
        P.Value    = tt$P.Value,
        adj.P.Val  = tt$adj.P.Val,
        stringsAsFactors = FALSE
    )

    # Override logFC from pre-scaling matrix if provided
    if (!is.null(mat_for_fc) && !identical(mat, mat_for_fc)) {
        groups <- parse_metab_contrast(contrast_str)
        idx_A <- which(condition == groups$denominator)
        idx_B <- which(condition == groups$numerator)
        for (i in seq_len(nrow(mat_for_fc))) {
            res$logFC[i] <- mean(mat_for_fc[i, idx_B], na.rm = TRUE) -
                            mean(mat_for_fc[i, idx_A], na.rm = TRUE)
        }
    }

    attr(res, "fit") <- fit2
    res
}


# ---- nonparametric / parametric two-group test (shared) --------------------

#' Per-feature two-group test (shared logic for t-test and Wilcoxon)
#'
#' @param mat          Numeric matrix (features x samples) for statistical test.
#' @param condition    Factor of conditions.
#' @param contrast_str Character, e.g. "B - A".
#' @param mat_for_fc   Optional pre-scaling matrix for logFC computation.
#' @param test_fn      Function(vals_B, vals_A) returning list(statistic, p.value).
#' @return data.frame with feature_id, logFC, AveExpr, statistic, P.Value, adj.P.Val.
de_two_group <- function(mat, condition, contrast_str, mat_for_fc = NULL, test_fn) {
    groups <- parse_metab_contrast(contrast_str)
    idx_A <- which(condition == groups$denominator)
    idx_B <- which(condition == groups$numerator)

    if (length(idx_A) == 0 || length(idx_B) == 0) {
        stop("No samples for one of the groups in contrast: ", contrast_str)
    }

    fc_mat <- if (!is.null(mat_for_fc)) mat_for_fc else mat

    res <- data.frame(
        feature_id = rownames(mat),
        logFC      = NA_real_,
        AveExpr    = NA_real_,
        statistic  = NA_real_,
        P.Value    = NA_real_,
        stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(mat))) {
        fc_A <- fc_mat[i, idx_A]
        fc_B <- fc_mat[i, idx_B]
        res$logFC[i]   <- mean(fc_B, na.rm = TRUE) - mean(fc_A, na.rm = TRUE)
        res$AveExpr[i] <- mean(c(fc_A, fc_B), na.rm = TRUE)

        tt <- tryCatch(
            test_fn(mat[i, idx_B], mat[i, idx_A]),
            error = function(e) NULL
        )
        if (!is.null(tt)) {
            res$statistic[i] <- tt$statistic
            res$P.Value[i]   <- tt$p.value
        }
    }

    res$adj.P.Val <- stats::p.adjust(res$P.Value, method = "BH")
    res
}


# ---- t-test ----------------------------------------------------------------

#' Run Welch t-tests per feature on a single contrast
de_t_test <- function(mat, condition, contrast_str, mat_for_fc = NULL) {
    de_two_group(mat, condition, contrast_str, mat_for_fc,
                 test_fn = function(b, a) stats::t.test(b, a, var.equal = FALSE))
}


# ---- wilcoxon --------------------------------------------------------------

#' Run Wilcoxon rank-sum tests per feature on a single contrast
de_wilcoxon <- function(mat, condition, contrast_str, mat_for_fc = NULL) {
    de_two_group(mat, condition, contrast_str, mat_for_fc,
                 test_fn = function(b, a) stats::wilcox.test(b, a, exact = FALSE))
}


# ---- helpers ----------------------------------------------------------------

#' Parse a metabolomics contrast string "B - A" into numerator/denominator
#'
#' Splits on " - " (space-dash-space) so that group names containing hyphens
#' (e.g. "pre-treatment - post-treatment") are handled correctly.
#'
#' @param contrast_str Character, e.g. "post-treatment - pre-treatment".
#' @return list(numerator, denominator)
parse_metab_contrast <- function(contrast_str) {
    parts <- strsplit(contrast_str, " - ", fixed = TRUE)[[1]]
    if (length(parts) != 2) {
        stop("Cannot parse contrast: '", contrast_str,
             "'. Expected format: 'groupB - groupA'")
    }
    list(numerator = trimws(parts[1]), denominator = trimws(parts[2]))
}


#' Create a clean label from a contrast string
#'
#' Converts "post-treatment - pre-treatment" to "post-treatment_vs_pre-treatment".
make_contrast_label <- function(contrast_str) {
    parts <- strsplit(contrast_str, " - ", fixed = TRUE)[[1]]
    paste(trimws(parts), collapse = "_vs_")
}


#' Build wide summary_df from per-contrast DE tables
#'
#' @param de_tables Named list of per-contrast data.frames.
#' @param padj_cutoff Numeric, adjusted p-value threshold.
#' @param log2fc_cut  Numeric, absolute log2 FC threshold.
#' @return Wide data.frame with feature_id + per-contrast columns + pass flags.
build_de_summary <- function(de_tables, padj_cutoff, log2fc_cut) {
    contrast_names <- names(de_tables)
    first <- de_tables[[1]]

    summary_df <- data.frame(
        feature_id = first$feature_id,
        stringsAsFactors = FALSE
    )

    for (ctr in contrast_names) {
        tbl <- de_tables[[ctr]]

        summary_df[[paste0("logFC_", ctr)]]     <- tbl$logFC
        summary_df[[paste0("AveExpr_", ctr)]]   <- tbl$AveExpr
        summary_df[[paste0("P.Value_", ctr)]]   <- tbl$P.Value
        summary_df[[paste0("adj.P.Val_", ctr)]] <- tbl$adj.P.Val

        # Significance flag
        pass <- as.integer(
            !is.na(tbl$adj.P.Val) &
            tbl$adj.P.Val < padj_cutoff &
            abs(tbl$logFC) >= log2fc_cut
        )
        summary_df[[paste0("pass_", ctr)]] <- pass
    }

    # Aggregate pass flag across contrasts
    pass_cols <- grep("^pass_", colnames(summary_df), value = TRUE)
    summary_df$pass_any_contrast <- as.integer(
        rowSums(summary_df[, pass_cols, drop = FALSE], na.rm = TRUE) > 0
    )

    summary_df
}


#' Extract a per-contrast DE table from summary_df for plotting
#'
#' Returns a data.frame with columns logFC, P.Value, adj.P.Val, AveExpr
#' matching the core plot_volcano / plot_ma signature.
#'
#' @param summary_df Wide DE summary.
#' @param contrast   Contrast label (e.g. "2024_vs_2013").
#' @return data.frame suitable for plot_volcano / plot_ma.
extract_contrast_table <- function(summary_df, contrast) {
    data.frame(
        feature_id = summary_df$feature_id,
        logFC      = summary_df[[paste0("logFC_", contrast)]],
        AveExpr    = summary_df[[paste0("AveExpr_", contrast)]],
        P.Value    = summary_df[[paste0("P.Value_", contrast)]],
        adj.P.Val  = summary_df[[paste0("adj.P.Val_", contrast)]],
        stringsAsFactors = FALSE
    )
}
