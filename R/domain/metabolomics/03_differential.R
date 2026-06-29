# R/domain/metabolomics/03_differential.R
#
# Differential expression analysis for metabolomics data.
# Supports limma, Welch t-test, equal-variance t-test, and Wilcoxon rank-sum methods.
#
# Operates on the standard pre-processing contract (expr_work, meta).
# Returns a wide-format summary_df compatible with the DE contract.
#
# Contrast loading follows the RNA/proteomics standard: structured table from
# files$contrasts (CSV/TSV with Contrast_name, Factor, Numerator, Denominator)
# loaded via load_omics_inputs(config, mode = "metabolomics").
#
# Reuses: assert_numeric_matrix, assert_one_of, normalize_contrast_name, %||%


# ---- local helpers -----------------------------------------------------------

#' Vectorised logical coercion (NA -> FALSE)
#' @keywords internal
isTRUE_vec <- function(x) {
    out <- as.logical(x)
    out[is.na(out)] <- FALSE
    out
}


#' Filter a matrix + metadata to biological samples only (exclude QC/blanks)
#'
#' Removes samples whose \code{condition_col} value matches "qc" or "blank"
#' (case-insensitive), plus any rows flagged by \code{is_QC} or \code{is_blank}
#' metadata columns.
#'
#' @param mat           Numeric matrix (features x samples).
#' @param meta          data.frame with at least \code{sample_col} and
#'                      \code{condition_col}.
#' @param condition_col Column in meta identifying experimental groups.
#' @param sample_col    Column in meta identifying sample IDs.
#' @param label         Character label for log messages.
#' @return list(mat, meta, condition) — filtered matrix, metadata, and factor.
#' @keywords internal
filter_to_biological <- function(mat, meta, condition_col, sample_col,
                                 label = "metabolomics") {
  condition_vals <- trimws(as.character(meta[[condition_col]]))
  sample_vals    <- trimws(as.character(meta[[sample_col]]))
  
  is_qc_or_blank_condition <- grepl("(^|_)(qc|blank|blanks)($|_)",
                                    condition_vals, ignore.case = TRUE)
  is_qc_sample <- grepl("^QC", sample_vals, ignore.case = TRUE)
  
  is_bio <- !(is_qc_or_blank_condition | is_qc_sample)
  
  if ("is_QC" %in% colnames(meta))
    is_bio <- is_bio & !isTRUE_vec(meta[["is_QC"]])
  if ("is_blank" %in% colnames(meta))
    is_bio <- is_bio & !isTRUE_vec(meta[["is_blank"]])
  
  n_excluded <- sum(!is_bio)
  if (n_excluded > 0L) {
    message(sprintf(
      "%s: excluding %d non-biological sample(s) (QC/blank); retaining %d",
      label, n_excluded, sum(is_bio)
    ))
    keep_ids <- meta[[sample_col]][is_bio]
    mat  <- mat[, keep_ids, drop = FALSE]
    meta <- meta[is_bio, , drop = FALSE]
  }
  
  list(mat = mat, meta = meta, condition = factor(meta[[condition_col]]))
}


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
run_metabolomics_de <- function(pre, config, contrast_table) {
    cfg <- config$modes$metabolomics
    de_cfg <- cfg$de %||% list()

    method <- tolower(de_cfg$method %||% "limma")
    assert_one_of(method, "de$method", c("limma", "t_test", "t_test_equal", "wilcoxon"))

    condition_col <- de_cfg$condition_column %||% cfg$effects$color %||% "sample_type"
    sample_col <- cfg$effects$samples %||% "sample_id"

    # Contrast loading: prefer explicit contrast_table argument (from pipeline),
    # fall back to config$de$contrasts for backwards compatibility
    if (!missing(contrast_table) && !is.null(contrast_table)) {
        contrasts <- contrast_table
    } else {
        contrasts <- de_cfg$contrasts
    }
    if (is.null(contrasts) || length(contrasts) == 0) {
        stop("metabolomics DE: contrasts are required (via contrast_table argument or config$de$contrasts).")
    }
    if (is.data.frame(contrasts)) {
        req <- c("Numerator", "Denominator")
        miss <- setdiff(req, colnames(contrasts))
        if (length(miss) > 0) {
            stop("metabolomics DE: contrast_table is missing required column(s): ",
                 paste(miss, collapse = ", "))
        }
        contrast_labels_override <- if ("Contrast_name" %in% colnames(contrasts))
            as.character(contrasts$Contrast_name) else NULL
        contrasts <- paste(as.character(contrasts$Numerator),
                           as.character(contrasts$Denominator),
                           sep = " - ")
    } else {
        contrast_labels_override <- NULL
        if (is.list(contrasts)) contrasts <- unlist(contrasts)
    }

    mat  <- pre$expr_work
    # When variance-scaling (auto/pareto/range) is applied, DE tests must use
    # the pre-scaling matrix (expr_log) because scaling distorts within-group
    # variance and inflates p-values.  But when scaling is "none" or "center",
    # expr_work should be used — it contains the chosen normalization (e.g.
    # EigenMS) which is essential for the test.
    scaling_used <- pre$info$normalization$scaling %||% "none"
    if (scaling_used %in% c("auto", "pareto", "range")) {
        mat_for_test <- pre$expr_log %||% mat
        message("metabolomics DE: using pre-scaling matrix (expr_log) because scaling = '", scaling_used, "'")
    } else {
        mat_for_test <- mat
    }
    meta <- pre$meta
    assert_numeric_matrix(mat_for_test, "metab_expr_for_test")

    # Align metadata to matrix columns
    meta <- meta[match(colnames(mat_for_test), meta[[sample_col]]), , drop = FALSE]
    condition <- factor(meta[[condition_col]])

    # Thresholds for significance flags
    padj_cutoff <- de_cfg$p_cutoff %||% 0.05

    # logfc_cutoff: applied directly to |logFC| in whatever scale the DE matrix uses.
    # linear_fc_cutoff: convenience alias — converted via log2() for backwards compat.
    if (!is.null(de_cfg$logfc_cutoff)) {
        log2fc_cut <- de_cfg$logfc_cutoff
    } else {
        linear_fc  <- de_cfg$linear_fc_cutoff %||% 1.5
        log2fc_cut <- log2(linear_fc)
    }

    # Raw filtered matrix for computing log2(FC) from intensity ratios
    # (MetaboAnalyst-compatible: FC = mean_B / mean_A on raw scale)
    mat_raw <- pre$expr_filt

    # Run DE for each contrast
    de_tables <- list()
    de_model  <- NULL

    for (i in seq_along(contrasts)) {
        ctr <- contrasts[i]
        ctr_label <- if (!is.null(contrast_labels_override))
            contrast_labels_override[i] else make_contrast_label(ctr)
        message("metabolomics DE [", method, "]: ", ctr)

        tbl <- switch(method,
            limma        = de_limma(mat_for_test, condition, ctr, mat_for_fc = mat_raw),
            t_test       = de_t_test(mat_for_test, condition, ctr, mat_for_fc = mat_raw),
            t_test_equal = de_t_test_equal(mat_for_test, condition, ctr, mat_for_fc = mat_raw),
            wilcoxon     = de_wilcoxon(mat_for_test, condition, ctr, mat_for_fc = mat_raw)
        )

        # Capture limma model from first contrast
        if (method == "limma" && is.null(de_model)) {
            de_model <- attr(tbl, "fit")
        }

        de_tables[[ctr_label]] <- tbl
    }

    # Build wide summary_df
    summary_df <- build_de_summary(de_tables, padj_cutoff, log2fc_cut)

    # Add Name column from row_data$original_id (HMDB → compound name)
    if (!is.null(pre$row_data) && "original_id" %in% colnames(pre$row_data)) {
        idx <- match(summary_df$feature_id, rownames(pre$row_data))
        names_resolved <- pre$row_data$original_id[idx]
        # Only add if names differ from feature_id (i.e., actually resolved)
        if (!all(is.na(names_resolved)) && !all(names_resolved == summary_df$feature_id, na.rm = TRUE)) {
            summary_df$Name <- names_resolved
            # Move Name to second column
            cols <- colnames(summary_df)
            summary_df <- summary_df[, c("feature_id", "Name", setdiff(cols, c("feature_id", "Name")))]
        }
    }

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

    # Override logFC with log2(FC) from raw intensities if provided
    if (!is.null(mat_for_fc) && !identical(mat, mat_for_fc)) {
        groups <- parse_metab_contrast(contrast_str)
        idx_A <- which(condition == groups$denominator)
        idx_B <- which(condition == groups$numerator)
        for (i in seq_len(nrow(mat_for_fc))) {
            mean_A <- mean(mat_for_fc[i, idx_A], na.rm = TRUE)
            mean_B <- mean(mat_for_fc[i, idx_B], na.rm = TRUE)
            if (is.na(mean_A) || mean_A <= 0) mean_A <- 1e-10
            if (is.na(mean_B) || mean_B <= 0) mean_B <- 1e-10
            res$logFC[i] <- log2(mean_B / mean_A)
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

    res <- data.frame(
        feature_id = rownames(mat),
        logFC      = NA_real_,
        AveExpr    = NA_real_,
        statistic  = NA_real_,
        P.Value    = NA_real_,
        stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(mat))) {
        # Compute log2(FC) from raw intensities if available
        # (MetaboAnalyst-compatible: FC = mean_B / mean_A, then log2)
        if (!is.null(mat_for_fc)) {
            raw_A <- mat_for_fc[i, idx_A]
            raw_B <- mat_for_fc[i, idx_B]
            mean_A <- mean(raw_A, na.rm = TRUE)
            mean_B <- mean(raw_B, na.rm = TRUE)
            # Avoid division by zero; use small pseudocount
            if (is.na(mean_A) || mean_A <= 0) mean_A <- 1e-10
            if (is.na(mean_B) || mean_B <= 0) mean_B <- 1e-10
            res$logFC[i]   <- log2(mean_B / mean_A)
            res$AveExpr[i] <- log2((mean_A + mean_B) / 2)
        } else {
            # Fallback: difference on the transformed scale
            fc_A <- mat[i, idx_A]
            fc_B <- mat[i, idx_B]
            res$logFC[i]   <- mean(fc_B, na.rm = TRUE) - mean(fc_A, na.rm = TRUE)
            res$AveExpr[i] <- mean(c(fc_A, fc_B), na.rm = TRUE)
        }

        # P-value always from transformed (glog10) data
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


# ---- equal-variance t-test -------------------------------------------------

#' Run equal-variance (Student's) t-tests per feature on a single contrast
#' (MetaboAnalyst default)
de_t_test_equal <- function(mat, condition, contrast_str, mat_for_fc = NULL) {
    de_two_group(mat, condition, contrast_str, mat_for_fc,
                 test_fn = function(b, a) stats::t.test(b, a, var.equal = TRUE))
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

    # Name column will be added from row_data after build_de_summary returns

    for (ctr in contrast_names) {
        tbl <- de_tables[[ctr]]

        # Signed linear FC from logFC (consistent with precomputed loader)
        lfc <- tbl$logFC
        linear_fc_signed <- ifelse(lfc >= 0, 2^lfc, -(2^abs(lfc)))

        summary_df[[paste0("linearFC.", ctr)]]   <- signif(linear_fc_signed, 3)
        summary_df[[paste0("AveExpr.", ctr)]]    <- tbl$AveExpr
        summary_df[[paste0("P.Value.", ctr)]]    <- tbl$P.Value
        summary_df[[paste0("adj.P.Val.", ctr)]]  <- tbl$adj.P.Val

        # Significance flag (uses raw P.Value, not adjusted)
        pass <- as.integer(
            !is.na(tbl$P.Value) &
            tbl$P.Value < padj_cutoff &
            abs(tbl$logFC) >= log2fc_cut
        )
        summary_df[[paste0("pass.", ctr)]] <- pass
    }

    # Aggregate pass flag across contrasts
    pass_cols <- grep("^pass\\.", colnames(summary_df), value = TRUE)
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
    # Convert linearFC back to logFC for plotting compatibility
    linear_fc <- summary_df[[paste0("linearFC.", contrast)]]
    logfc <- ifelse(linear_fc >= 0, log2(linear_fc), -log2(abs(linear_fc)))

    data.frame(
        feature_id = summary_df$feature_id,
        logFC      = logfc,
        AveExpr    = summary_df[[paste0("AveExpr.", contrast)]],
        P.Value    = summary_df[[paste0("P.Value.", contrast)]],
        adj.P.Val  = summary_df[[paste0("adj.P.Val.", contrast)]],
        stringsAsFactors = FALSE
    )
}
