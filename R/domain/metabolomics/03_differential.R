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
  assert_one_of(method, "de$method", c("limma", "t_test", "wilcoxon"))
  
  condition_col <- de_cfg$condition_column %||% cfg$effects$color %||% "sample_type"
  sample_col <- cfg$effects$samples %||% "sample_id"
  
  # Use the current DE-ready metabolomics matrix.
  # In the current pipeline this is expr_work
  # (normalized + transformed + optional drift correction; no feature scaling).
  mat_for_test <- pre$expr_work
  meta <- pre$meta
  assert_numeric_matrix(mat_for_test, "metab_expr_for_test")
  
  # Align metadata to matrix columns
  meta <- meta[match(colnames(mat_for_test), meta[[sample_col]]), , drop = FALSE]
  
  # ---- Filter to biological samples only (exclude QC/blanks) ----
  bio <- filter_to_biological(mat_for_test, meta, condition_col, sample_col,
                              label = "metabolomics DE")
  mat_for_test <- bio$mat
  meta         <- bio$meta
  condition    <- bio$condition
  
  # Thresholds for significance flags
  padj_cutoff <- de_cfg$p_cutoff %||% 0.05
  linear_fc   <- de_cfg$linear_fc_cutoff %||% 1.5
  log2fc_cut  <- log2(linear_fc)
  
  # Run DE — limma fits all contrasts at once (same pattern as proteomics);
  # t_test / wilcoxon use a per-contrast loop.
  de_tables <- list()
  de_model  <- NULL
  design          <- NULL
  contrast_matrix <- NULL
  contrast_formulas <- NULL
  
  if (method == "limma") {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("Package 'limma' is required for limma DE.")
    }
    
    # -- Validate Factor column (same pattern as proteomics) --
    factor_col <- unique(contrast_table$Factor)
    if (length(factor_col) != 1 || factor_col != condition_col) {
      stop(sprintf(
        "All contrasts must use the same Factor column matching condition_col. Got Factor='%s', expected='%s'",
        paste(factor_col, collapse = ", "), condition_col
      ))
    }
    
    # -- Safe factor levels (same as proteomics) --
    raw_levels  <- levels(condition)
    safe_levels <- make.names(raw_levels)
    levels(condition) <- safe_levels
    
    # -- Design matrix (once) --
    design <- stats::model.matrix(~ 0 + condition)
    colnames(design) <- safe_levels
    
    # -- All contrast formulas at once (same API as proteomics) --
    safe_num <- make.names(contrast_table$Numerator)
    safe_den <- make.names(contrast_table$Denominator)
    contrast_names <- vapply(contrast_table$Contrast_name,
                             normalize_contrast_name, character(1))
    contrast_formulas <- setNames(
      paste(safe_num, safe_den, sep = " - "),
      contrast_names
    )
    
    contrast_matrix <- limma::makeContrasts(contrasts = contrast_formulas,
                                            levels = design)
    colnames(contrast_matrix) <- names(contrast_formulas)
    
    # -- Single fit for all contrasts --
    fit  <- limma::lmFit(mat_for_test, design)
    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
    # eBayes - moderated t-statistics with shrinkage of residual variance
    # (matches proteomics + RNA-seq DE for consistency across omics)
    fit2 <- limma::eBayes(fit2)
    de_model <- fit2
    
    # -- Extract per-contrast tables --
    for (cn in colnames(contrast_matrix)) {
      message("metabolomics DE [limma]: ", cn)
      tt <- limma::topTable(fit2, coef = cn, number = Inf, sort.by = "none")
      de_tables[[cn]] <- data.frame(
        feature_id = rownames(tt),
        Contrast   = cn,
        logFC      = tt$logFC,
        AveExpr    = tt$AveExpr,
        t          = tt$t,
        P.Value    = tt$P.Value,
        adj.P.Val  = tt$adj.P.Val,
        # B statistic from eBayes (log-odds of differential expression)
        B          = tt$B,
        stringsAsFactors = FALSE

      )
    }
  } else {
    # t_test / wilcoxon: per-contrast loop
    for (i in seq_len(nrow(contrast_table))) {
      ctr_name    <- normalize_contrast_name(contrast_table$Contrast_name[i])
      numerator   <- as.character(contrast_table$Numerator[i])
      denominator <- as.character(contrast_table$Denominator[i])
      
      message("metabolomics DE [", method, "]: ", ctr_name,
              " (", numerator, " vs ", denominator, ")")
      
      tbl <- switch(method,
                    t_test   = de_t_test(mat_for_test, condition, numerator, denominator),
                    wilcoxon = de_wilcoxon(mat_for_test, condition, numerator, denominator)
      )
      tbl$Contrast <- ctr_name
      de_tables[[ctr_name]] <- tbl
    }
  }
  
  # Build wide summary_df
  summary_df <- build_de_summary(de_tables, padj_cutoff, log2fc_cut)
  
  message("metabolomics DE complete: ", nrow(summary_df), " features, ",
          sum(summary_df$pass_any_contrast == 1, na.rm = TRUE), " significant")
  
  list(
    summary_df        = summary_df,
    method            = method,
    de_tables         = de_tables,
    de_model          = de_model,
    design            = design,
    contrast_formulas = contrast_formulas,
    contrast_matrix   = contrast_matrix
  )
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


# ---- equal-variance t-test -------------------------------------------------

#' Run equal-variance (Student's) t-tests per feature on a single contrast
#' (MetaboAnalyst default)
de_t_test_equal <- function(mat, condition, contrast_str, mat_for_fc = NULL) {
    de_two_group(mat, condition, contrast_str, mat_for_fc,
                 test_fn = function(b, a) stats::t.test(b, a, var.equal = TRUE))
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