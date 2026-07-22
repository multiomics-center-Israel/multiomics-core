#' Batch correction for proteomics intensity matrices
#'
#' Applies a batch-correction transform to the log2 proteomics intensity matrix
#' produced by `preprocess_proteomics()`. The correction is driven by the
#' `modes.proteomics.batch_correction` config block and supports two backends:
#'
#'   * `"combat"`   — empirical-Bayes location/scale adjustment via `sva::ComBat`.
#'   * `"probatch"` — per-feature, per-batch median centering via
#'                    `proBatch::center_feature_batch_medians_dm` (proBatch's
#'                    discrete-batch correction).
#'
#' Both backends estimate and apply the correction on the imputed complete
#' matrix (`pre$expr_imp_single`). The corrected values feed `expr_work` and
#' `expr_imp_single` in full; `expr_filt` receives the corrected values with its
#' original missing-value pattern restored, so it keeps its "filtered,
#' unimputed" contract (the DE step re-imputes it, and CV / Shiny exports still
#' exclude imputed measurements). All downstream DE, QC and clustering consume
#' the corrected data.
#'
#' These functions are side-effect free: they take a matrix (or `pre` list) and
#' return a corrected matrix (or modified `pre`). File I/O, if any, lives in the
#' module wrapper / pipeline target.

#' Apply batch correction to a preprocessed proteomics object
#'
#' Dispatches on `config$modes$proteomics$batch_correction$method`. When the
#' method is `"none"` (or the block is absent / disabled) the input is returned
#' unchanged. Otherwise `expr_work` and `expr_imp_single` become the complete
#' corrected matrix, while `expr_filt` (and `expr_filt_pre_imp`) get the
#' corrected values with the original NA pattern restored.
#'
#' @param pre Output of `preprocess_proteomics()` (a list with at least
#'   `expr_filt` (log2 proteins x samples matrix) and `meta`).
#' @param config Full pipeline config.
#' @param seed Integer seed for the (small) stochastic component of ComBat's
#'   empirical-Bayes estimation; passed through `withr::with_seed()`.
#' @return The `pre` list, with the corrected matrix in `expr_work`/
#'   `expr_imp_single` (complete) and `expr_filt`/`expr_filt_pre_imp` (corrected
#'   observed values, NA where originally missing), plus `batch_corrected`/
#'   `batch_method` metadata fields when a correction was applied.
correct_batch_proteomics <- function(pre, config, seed = 42) {
  cfg <- config$modes$proteomics
  bc  <- get_proteomics_batch_config(cfg)

  if (!isTRUE(bc$enabled) || identical(bc$method, "none")) {
    message("Proteomics batch correction: disabled. Passing through unchanged.")
    return(pre)
  }

  # ComBat/probatch need a complete matrix. preprocess_proteomics() already
  # produces a single-imputation complete matrix (expr_imp_single); estimate and
  # apply the correction on that. The original missing-value pattern is restored
  # to expr_filt afterwards (see below) so the DE step still re-imputes it.
  expr <- pre$expr_imp_single %||% pre$expr_filt
  assert_numeric_matrix(expr, "expr_imp_single")

  # ComBat (and probatch median-centering) misbehave on non-finite input. NAs
  # should already be imputed here; -Inf can still arrive from log2 of zero or
  # negative upstream intensities, which assert_numeric_matrix() does not catch.
  if (!all(is.finite(expr))) {
    cli::cli_abort(c(
      "Proteomics batch correction requires a finite intensity matrix.",
      "i" = "Found {sum(is.na(expr))} NA and {sum(is.infinite(expr))} non-finite value(s) in the imputed matrix.",
      "i" = "Non-finite values usually come from log2 of zero/negative intensities upstream; filter or floor them before this step."
    ))
  }
  meta <- pre$meta

  sample_col <- cfg$effects$samples %||% cfg$id_columns$sample_col
  batch_var  <- resolve_batch_variable(meta, expr, bc$batch_col, sample_col)
  group_var  <- resolve_group_variable(meta, expr, bc$group_col %||% cfg$effects$color,
                                       sample_col, preserve = bc$preserve_condition)
  covariates <- resolve_batch_covariates(meta, expr, bc$covariates, sample_col)

  corrected <- switch(
    bc$method,
    combat   = correct_proteomics_combat(expr, batch_var, group_var, covariates,
                                         preserve_group = bc$preserve_condition, seed = seed),
    probatch = correct_proteomics_probatch(expr, batch_var,
                                           sample_ids = colnames(expr)),
    cli::cli_abort(c(
      "Unknown proteomics batch_correction method: {.val {bc$method}}.",
      "i" = "Supported methods: {.val none}, {.val combat}, {.val probatch}."
    ))
  )

  assert_numeric_matrix(corrected, "batch_corrected")
  message(sprintf("Proteomics batch correction: applied '%s' on batch column '%s'.",
                  bc$method, bc$batch_col))

  # Keep expr_filt's "filtered, unimputed" contract: corrected values at
  # observed positions, NA where the measurement was originally missing. This
  # fixes two things at once —
  #   * the DE step re-imputes expr_filt, so its multiple imputations stay
  #     genuinely different (a complete matrix collapses them to a single draw
  #     and the imputation-consensus stops measuring robustness), and
  #   * build_group_cv_proteomics() and the Shiny expr_raw export still exclude
  #     imputed measurements, which they detect via NAs in expr_filt.
  # expr_work / expr_imp_single stay complete + corrected for QC and plots.
  na_mask <- is.na(pre$expr_filt)[rownames(corrected), colnames(corrected), drop = FALSE]
  corrected_unimputed <- corrected
  corrected_unimputed[na_mask] <- NA_real_

  pre$expr_filt         <- corrected_unimputed
  pre$expr_filt_pre_imp <- corrected_unimputed
  pre$expr_work         <- corrected
  pre$expr_imp_single   <- corrected
  pre$batch_corrected   <- TRUE
  pre$batch_method      <- bc$method
  pre
}

#' Resolve proteomics batch-correction config with defaults
#'
#' @param cfg The `config$modes$proteomics` sub-config.
#' @return A list with `enabled`, `method`, `batch_col`, `group_col`,
#'   `covariates`, `preserve_condition`.
#' @noRd
get_proteomics_batch_config <- function(cfg) {
  bc <- cfg$batch_correction %||% list()
  method <- bc$method %||% "none"
  list(
    enabled            = bc$enabled %||% (!identical(method, "none")),
    method             = method,
    batch_col          = bc$batch_col %||% cfg$qc$batch_col,
    group_col          = bc$group_col,
    covariates         = bc$covariates,
    preserve_condition = bc$preserve_condition %||% TRUE
  )
}

#' Resolve a batch grouping vector aligned to matrix columns
#'
#' @param meta Sample metadata data.frame.
#' @param expr Intensity matrix (proteins x samples).
#' @param batch_col Name of the metadata column holding the batch label.
#' @param sample_col Name of the metadata column holding sample IDs.
#' @return A factor of batch labels, one per matrix column (sample).
#' @noRd
resolve_batch_variable <- function(meta, expr, batch_col, sample_col) {
  if (is.null(batch_col) || !nzchar(batch_col)) {
    cli::cli_abort(c(
      "Proteomics batch correction is enabled but no batch column is set.",
      "i" = "Set {.code modes.proteomics.batch_correction.batch_col} (or {.code modes.proteomics.qc.batch_col})."
    ))
  }
  if (!batch_col %in% colnames(meta)) {
    cli::cli_abort(c(
      "Batch column {.val {batch_col}} not found in proteomics metadata.",
      "i" = "Available columns: {.val {colnames(meta)}}."
    ))
  }
  vals <- align_meta_column(meta, expr, batch_col, sample_col)
  bv <- as.factor(as.character(vals))
  if (nlevels(bv) < 2) {
    cli::cli_abort(c(
      "Batch column {.val {batch_col}} has fewer than 2 levels; nothing to correct.",
      "i" = "Disable batch correction or pick a column with >= 2 batches."
    ))
  }
  bv
}

#' Resolve the biological group vector (for preserving condition), aligned to columns
#'
#' @param meta Sample metadata data.frame.
#' @param expr Intensity matrix (proteins x samples).
#' @param group_col Name of the metadata column holding the biological condition.
#' @param sample_col Name of the metadata column holding sample IDs.
#' @param preserve Logical; whether condition preservation was requested. When
#'   `TRUE`, a `group_col` that is set but absent from `meta` is treated as a
#'   config error rather than silently ignored.
#' @return A factor of group labels per column, or `NULL` when no group column
#'   is configured (nothing to preserve).
#' @noRd
resolve_group_variable <- function(meta, expr, group_col, sample_col, preserve = TRUE) {
  # No condition column configured anywhere -> genuinely nothing to preserve.
  if (is.null(group_col) || !nzchar(group_col)) {
    return(NULL)
  }
  # A column was configured but is absent. Silently dropping it would let ComBat
  # remove the biological condition along with the batch effect (worst when
  # condition and batch are correlated), so fail loudly on the likely typo when
  # preservation was requested.
  if (!group_col %in% colnames(meta)) {
    if (isTRUE(preserve)) {
      cli::cli_abort(c(
        "Condition column {.val {group_col}} for batch-correction preservation not found in metadata.",
        "i" = "Fix {.code modes.proteomics.batch_correction.group_col} (or {.code effects.color}), or set {.code preserve_condition: false} to correct without preserving condition.",
        "i" = "Available columns: {.val {colnames(meta)}}."
      ))
    }
    return(NULL)
  }
  as.factor(as.character(align_meta_column(meta, expr, group_col, sample_col)))
}

#' Resolve optional covariate columns into a model-matrix-ready data.frame
#'
#' @param meta Sample metadata data.frame.
#' @param expr Intensity matrix (proteins x samples).
#' @param covariates Character vector of metadata column names, or `NULL`.
#' @param sample_col Name of the metadata column holding sample IDs.
#' @return A data.frame (one row per sample, aligned to columns) or `NULL`.
#' @noRd
resolve_batch_covariates <- function(meta, expr, covariates, sample_col) {
  if (is.null(covariates) || length(covariates) == 0) return(NULL)
  covariates <- intersect(covariates, colnames(meta))
  if (length(covariates) == 0) return(NULL)
  df <- as.data.frame(lapply(covariates, function(cc) {
    align_meta_column(meta, expr, cc, sample_col)
  }), stringsAsFactors = FALSE)
  colnames(df) <- covariates
  df
}

#' Align a metadata column to matrix columns by sample ID
#'
#' Falls back to assuming row order matches column order when a sample-ID column
#' is unavailable or does not cover all matrix columns.
#'
#' @param meta Sample metadata data.frame.
#' @param expr Intensity matrix (proteins x samples).
#' @param col Metadata column to extract.
#' @param sample_col Name of the metadata column holding sample IDs.
#' @return A vector of length `ncol(expr)`.
#' @noRd
align_meta_column <- function(meta, expr, col, sample_col) {
  samples <- colnames(expr)
  if (!is.null(sample_col) && sample_col %in% colnames(meta) && !is.null(samples) &&
      all(samples %in% as.character(meta[[sample_col]]))) {
    idx <- match(samples, as.character(meta[[sample_col]]))
    return(meta[[col]][idx])
  }
  if (nrow(meta) != ncol(expr)) {
    cli::cli_abort(c(
      "Cannot align metadata column {.val {col}} to the intensity matrix.",
      "i" = "Metadata has {nrow(meta)} rows but matrix has {ncol(expr)} columns, and sample IDs do not match column names."
    ))
  }
  meta[[col]]
}

#' ComBat batch correction for a log2 proteomics matrix
#'
#' Wraps `sva::ComBat`, which expects log-scale data and corrects additive and
#' multiplicative batch effects via empirical Bayes. When `preserve_group` is
#' `TRUE` (the default) and a group vector is supplied, the biological condition
#' (and any covariates) are encoded in the ComBat model matrix so that the
#' design effect is preserved rather than removed.
#'
#' @param expr Numeric log2 matrix (proteins x samples). Must be complete
#'   (no NA) — ComBat does not accept missing values.
#' @param batch_var Factor of batch labels (length `ncol(expr)`).
#' @param group_var Factor of biological group labels, or `NULL`.
#' @param covariates Optional data.frame of additional covariates, or `NULL`.
#' @param preserve_group Logical; encode `group_var`/`covariates` in the ComBat
#'   model so the condition effect is retained.
#' @param seed Integer seed passed through `withr::with_seed()`.
#' @return Corrected numeric matrix with the same dim/dimnames as `expr`.
correct_proteomics_combat <- function(expr, batch_var, group_var = NULL,
                                      covariates = NULL, preserve_group = TRUE,
                                      seed = 42) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg sva} is required for {.val combat} batch correction but is not installed.",
      "i" = "Install it with {.code BiocManager::install('sva')} (a dependency change must be approved by a human)."
    ))
  }

  if (!all(is.finite(expr))) {
    cli::cli_abort(c(
      "ComBat requires a matrix of finite values.",
      "i" = "Found {sum(is.na(expr))} NA and {sum(is.infinite(expr))} non-finite value(s) (e.g. from log2 of zero/negative intensities).",
      "i" = "Impute NAs and floor/remove non-finite values before batch correction."
    ))
  }

  mod <- NULL
  if (isTRUE(preserve_group) && (!is.null(group_var) || !is.null(covariates))) {
    md <- data.frame(.row = seq_len(ncol(expr)))
    if (!is.null(group_var)) md[[".group"]] <- group_var
    if (!is.null(covariates)) md <- cbind(md, covariates)
    md$.row <- NULL
    if (ncol(md) > 0) {
      mod <- stats::model.matrix(~., data = md)
    }
  }

  withr::with_seed(seed, {
    sva::ComBat(
      dat       = as.matrix(expr),
      batch     = batch_var,
      mod       = mod,
      par.prior = TRUE
    )
  })
}

#' proBatch median-centering batch correction for a log2 proteomics matrix
#'
#' Uses `proBatch::center_feature_batch_medians_dm`, proBatch's discrete-batch
#' correction: for each feature, each batch is shifted so its per-batch median
#' matches the global per-feature median. Operates directly on the data matrix
#' (`_dm` interface), keeping this step matrix-in / matrix-out.
#'
#' @param expr Numeric log2 matrix (proteins x samples).
#' @param batch_var Factor of batch labels (length `ncol(expr)`).
#' @param sample_ids Character vector of sample IDs (matrix column names).
#' @return Corrected numeric matrix with the same dim/dimnames as `expr`.
correct_proteomics_probatch <- function(expr, batch_var, sample_ids = colnames(expr)) {
  if (!requireNamespace("proBatch", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg proBatch} is required for {.val probatch} batch correction but is not installed.",
      "i" = "Install it with {.code BiocManager::install('proBatch')} (a dependency change must be approved by a human)."
    ))
  }

  expr <- as.matrix(expr)
  if (is.null(sample_ids)) sample_ids <- paste0("S", seq_len(ncol(expr)))
  colnames(expr) <- sample_ids

  sample_annotation <- data.frame(
    FullRunName = sample_ids,
    batch       = as.character(batch_var),
    stringsAsFactors = FALSE
  )

  corrected <- proBatch::center_feature_batch_medians_dm(
    data_matrix       = expr,
    sample_annotation = sample_annotation,
    sample_id_col     = "FullRunName",
    batch_col         = "batch"
  )

  # proBatch returns the matrix in feature x sample shape; realign columns
  # defensively in case ordering changed.
  corrected <- corrected[rownames(expr), sample_ids, drop = FALSE]
  corrected
}
