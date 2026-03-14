# R/domain/metabolomics/03_differential.R
#
# Differential expression analysis for metabolomics data.
# Supports limma, Welch t-test, and Wilcoxon rank-sum methods.
#
# Operates on the standard pre-processing contract (expr_work, meta).
# Returns a wide-format summary_df compatible with the DE contract.
#
# Contrast loading follows the RNA/proteomics standard: structured table from
# files$contrasts (CSV/TSV with Contrast_name, Factor, Numerator, Denominator)
# loaded via load_omics_inputs(config, mode = "metabolomics").
#
# Reuses: assert_numeric_matrix, assert_one_of, normalize_contrast_name, %||%


# ---- pre-computed DE loader --------------------------------------------------

#' Load pre-computed metabolomics DE tables from config$files$de_table
#'
#' Reads CSV files with columns: FC, log2(FC), raw.pval, -log10(p).
#' Builds a summary_df conforming to the DE contract.
#'
#' @param config Full pipeline config.
#' @return list conforming to the DE contract: summary_df, method, de_tables
load_precomputed_metabolomics_de <- function(config) {
    cfg <- config$modes$metabolomics
    de_cfg <- cfg$de %||% list()

    de_files <- cfg$files$de_table
    if (is.list(de_files)) de_files <- unlist(de_files)

    padj_cutoff <- de_cfg$p_cutoff %||% 0.05
    linear_fc   <- de_cfg$linear_fc_cutoff %||% 1.5
    log2fc_cut  <- log2(linear_fc)

    # Derive contrast labels from file names
    contrast_labels <- vapply(de_files, function(f) {
        bn <- tools::file_path_sans_ext(basename(f))
        sub("^de_", "", bn)
    }, character(1), USE.NAMES = FALSE)

    de_tables <- list()
    for (i in seq_along(de_files)) {
        abs_path <- resolve_raw_path(config, de_files[i])
        if (!file.exists(abs_path)) {
            stop("Pre-computed DE table not found: ", abs_path)
        }

        raw <- read_table_auto(abs_path)

        # Map columns to standard names
        cn <- colnames(raw)

        # Feature IDs: unnamed first column (readr: "...1", base R: "X", or "")
        id_col_idx <- match(TRUE, cn %in% c("...1", "", "X", "V1"))
        feat_ids <- if (!is.na(id_col_idx)) {
            as.character(raw[[id_col_idx]])
        } else {
            rownames(raw)
        }

        # logFC: try common column name variants
        logfc_col <- cn[cn %in% c("log2(FC)", "log2.FC.", "logFC", "log2FC")][1]
        logfc_vals <- if (!is.na(logfc_col)) as.numeric(raw[[logfc_col]]) else NA_real_

        # P-value: try common column name variants
        pval_col <- cn[cn %in% c("raw.pval", "P.Value", "pvalue", "PValue", "p.value")][1]
        pval_vals <- if (!is.na(pval_col)) as.numeric(raw[[pval_col]]) else NA_real_

        tbl <- data.frame(
            feature_id = feat_ids,
            logFC      = logfc_vals,
            P.Value    = pval_vals,
            stringsAsFactors = FALSE
        )
        tbl$AveExpr <- NA_real_
        tbl$adj.P.Val <- stats::p.adjust(tbl$P.Value, method = "BH")

        de_tables[[contrast_labels[i]]] <- tbl
        message("  Loaded ", nrow(tbl), " features from ", basename(de_files[i]),
                " (label: ", contrast_labels[i], ")")
    }

    # Build summary_df with full outer join (tables may have different features)
    all_features <- unique(unlist(lapply(de_tables, function(t) t$feature_id)))
    summary_df <- data.frame(feature_id = all_features, stringsAsFactors = FALSE)

    for (ctr in names(de_tables)) {
        tbl <- de_tables[[ctr]]
        idx <- match(summary_df$feature_id, tbl$feature_id)

        # Signed linear FC from logFC (same transform as build_de_summary)
        lfc <- tbl$logFC[idx]
        linear_fc_signed <- ifelse(lfc >= 0, 2^lfc, -(2^abs(lfc)))

        summary_df[[paste0("linearFC.", ctr)]] <- signif(linear_fc_signed, 3)
        summary_df[[paste0("AveExpr.", ctr)]]  <- tbl$AveExpr[idx]
        summary_df[[paste0("pvalue.", ctr)]]   <- tbl$P.Value[idx]
        summary_df[[paste0("padj.", ctr)]]     <- tbl$adj.P.Val[idx]

        pass <- as.integer(
            !is.na(tbl$adj.P.Val[idx]) &
            tbl$adj.P.Val[idx] < padj_cutoff &
            abs(tbl$logFC[idx]) >= log2fc_cut
        )
        summary_df[[paste0("pass.", ctr)]] <- pass
    }

    summary_df <- add_pass_any_contrast(summary_df, pass_prefix = "^pass\\.")

    message("metabolomics precomputed DE: ", nrow(summary_df), " features, ",
            sum(summary_df$pass_any_contrast == 1, na.rm = TRUE), " significant")

    list(
        summary_df = summary_df,
        method     = "precomputed",
        de_tables  = de_tables,
        de_model   = NULL
    )
}


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

    # ---- Load structured contrast table (RNA/proteomics standard) ----
    inputs <- load_omics_inputs(config, mode = "metabolomics")
    contrast_table <- inputs$contrasts

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

    # Run DE for each row of the contrast table
    de_tables <- list()
    de_model  <- NULL

    for (i in seq_len(nrow(contrast_table))) {
        ctr_name    <- normalize_contrast_name(contrast_table$Contrast_name[i])
        numerator   <- as.character(contrast_table$Numerator[i])
        denominator <- as.character(contrast_table$Denominator[i])

        message("metabolomics DE [", method, "]: ", ctr_name,
                " (", numerator, " vs ", denominator, ")")

        tbl <- switch(method,
            limma    = de_limma(mat_for_test, condition, numerator, denominator),
            t_test   = de_t_test(mat_for_test, condition, numerator, denominator),
            wilcoxon = de_wilcoxon(mat_for_test, condition, numerator, denominator)
        )

        # Capture limma model from first contrast
        if (method == "limma" && is.null(de_model)) {
            de_model <- attr(tbl, "fit")
        }

        de_tables[[ctr_name]] <- tbl
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
#' @param mat         Numeric matrix (features x samples).
#' @param condition   Factor of conditions.
#' @param numerator   Character: numerator group level (e.g. "TreatmentA").
#' @param denominator Character: denominator group level (e.g. "Control").
#' @param mat_for_fc  Optional pre-scaling matrix for logFC override.
#' @return data.frame with feature_id, logFC, AveExpr, statistic, P.Value, adj.P.Val
de_limma <- function(mat, condition, numerator, denominator, mat_for_fc = NULL) {
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
    raw_levels <- levels(condition)
    safe_levels <- make.names(raw_levels)
    colnames(design) <- safe_levels

    # Build contrast formula from structured numerator/denominator
    safe_num <- make.names(numerator)
    safe_den <- make.names(denominator)
    safe_contrast <- paste(safe_num, "-", safe_den)

    fit <- limma::lmFit(mat_imp, design)
    contrast_matrix <- limma::makeContrasts(contrasts = safe_contrast,
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
        idx_A <- which(condition == denominator)
        idx_B <- which(condition == numerator)
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
#' @param numerator    Character: numerator group level.
#' @param denominator  Character: denominator group level.
#' @param mat_for_fc   Optional pre-scaling matrix for logFC computation.
#' @param test_fn      Function(vals_B, vals_A) returning list(statistic, p.value).
#' @return data.frame with feature_id, logFC, AveExpr, statistic, P.Value, adj.P.Val.
de_two_group <- function(mat, condition, numerator, denominator, mat_for_fc = NULL, test_fn) {
    idx_A <- which(condition == denominator)
    idx_B <- which(condition == numerator)

    if (length(idx_A) == 0 || length(idx_B) == 0) {
        stop("No samples for one of the groups in contrast: ",
             numerator, " vs ", denominator)
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
#'
#' @param mat         Numeric matrix (features x samples).
#' @param condition   Factor of conditions.
#' @param numerator   Character: numerator group level.
#' @param denominator Character: denominator group level.
#' @param mat_for_fc  Optional pre-scaling matrix for logFC override.
de_t_test <- function(mat, condition, numerator, denominator, mat_for_fc = NULL) {
    de_two_group(mat, condition, numerator, denominator, mat_for_fc,
                 test_fn = function(b, a) stats::t.test(b, a, var.equal = FALSE))
}


# ---- wilcoxon --------------------------------------------------------------

#' Run Wilcoxon rank-sum tests per feature on a single contrast
#'
#' @param mat         Numeric matrix (features x samples).
#' @param condition   Factor of conditions.
#' @param numerator   Character: numerator group level.
#' @param denominator Character: denominator group level.
#' @param mat_for_fc  Optional pre-scaling matrix for logFC override.
de_wilcoxon <- function(mat, condition, numerator, denominator, mat_for_fc = NULL) {
    de_two_group(mat, condition, numerator, denominator, mat_for_fc,
                 test_fn = function(b, a) stats::wilcox.test(b, a, exact = FALSE))
}


# ---- summary builder -------------------------------------------------------

#' Build wide summary_df from per-contrast DE tables
#'
#' Column naming follows the RNA-style contract (no `.imputs.` infix):
#'   linearFC.<cn>, pvalue.<cn>, padj.<cn>, pass.<cn>, AveExpr.<cn>
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

        # Signed linear FC: preserves directionality from limma logFC
        lfc <- tbl$logFC
        linear_fc_signed <- ifelse(lfc >= 0, 2^lfc, -(2^abs(lfc)))

        summary_df[[paste0("linearFC.", ctr)]] <- signif(linear_fc_signed, 3)
        summary_df[[paste0("AveExpr.", ctr)]]  <- tbl$AveExpr
        summary_df[[paste0("pvalue.", ctr)]]   <- tbl$P.Value
        summary_df[[paste0("padj.", ctr)]]     <- tbl$adj.P.Val

        # Significance flag (logic unchanged — same thresholds, same test)
        pass <- as.integer(
            !is.na(tbl$adj.P.Val) &
            tbl$adj.P.Val < padj_cutoff &
            abs(tbl$logFC) >= log2fc_cut
        )
        summary_df[[paste0("pass.", ctr)]] <- pass
    }

    # Aggregate pass flag across contrasts (reuse shared helper)
    summary_df <- add_pass_any_contrast(summary_df, pass_prefix = "^pass\\.")

    summary_df
}


#' Extract a per-contrast DE table from summary_df for plotting
#'
#' Reads aligned column names (linearFC., pvalue., padj.) from summary_df
#' and returns logFC (back-computed) for plot_volcano / plot_ma compatibility.
#'
#' @param summary_df Wide DE summary.
#' @param contrast   Contrast label (e.g. "2024_vs_2013").
#' @return data.frame suitable for plot_volcano / plot_ma.
extract_contrast_table <- function(summary_df, contrast) {
    linear_fc <- summary_df[[paste0("linearFC.", contrast)]]
    # Back-compute logFC from signed linearFC for plotting
    logfc <- ifelse(linear_fc >= 0, log2(linear_fc), -log2(abs(linear_fc)))

    data.frame(
        feature_id = summary_df$feature_id,
        logFC      = logfc,
        AveExpr    = summary_df[[paste0("AveExpr.", contrast)]],
        P.Value    = summary_df[[paste0("pvalue.", contrast)]],
        adj.P.Val  = summary_df[[paste0("padj.", contrast)]],
        stringsAsFactors = FALSE
    )
}
