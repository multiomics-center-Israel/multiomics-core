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
#' Both backends operate on the imputed log2 matrix (`pre$expr_filt`), and the
#' corrected matrix is propagated to `expr_filt`, `expr_work`, `expr_imp_single`,
#' so that downstream DE, QC and clustering all consume the corrected data.
#'
#' These functions are side-effect free: they take a matrix (or `pre` list) and
#' return a corrected matrix (or modified `pre`). File I/O, if any, lives in the
#' module wrapper / pipeline target.

#' Apply batch correction to a preprocessed proteomics object
#'
#' Dispatches on `config$modes$proteomics$batch_correction$method`. When the
#' method is `"none"` (or the block is absent / disabled) the input is returned
#' unchanged. Otherwise the corrected matrix replaces `expr_filt`, `expr_work`
#' and `expr_imp_single` in the returned `pre` list.
#'
#' @param pre Output of `preprocess_proteomics()` (a list with at least
#'   `expr_filt` (log2 proteins x samples matrix) and `meta`).
#' @param config Full pipeline config.
#' @param seed Integer seed for the (small) stochastic component of ComBat's
#'   empirical-Bayes estimation; passed through `withr::with_seed()`.
#' @return The `pre` list, with `expr_filt`/`expr_work`/`expr_imp_single`
#'   replaced by the corrected matrix and `batch_corrected`/`batch_method`
#'   metadata fields added when a correction was applied.
correct_batch_proteomics <- function(pre, config, seed = 42) {
  cfg <- config$modes$proteomics
  bc  <- get_proteomics_batch_config(cfg)

  if (!isTRUE(bc$enabled) || identical(bc$method, "none")) {
    message("Proteomics batch correction: disabled. Passing through unchanged.")
    return(pre)
  }

  # ComBat/probatch need a complete matrix. preprocess_proteomics() already
  # produces a single-imputation complete matrix (expr_imp_single); correct that
  # rather than the NA-containing expr_filt. Downstream DE re-runs imputation on
  # the (now complete) corrected matrix, where it is a no-op, so the corrected
  # values flow straight into the limma fits.
  expr <- pre$expr_imp_single %||% pre$expr_filt
  assert_numeric_matrix(expr, "expr_imp_single")
  meta <- pre$meta

  sample_col <- cfg$effects$samples %||% cfg$id_columns$sample_col
  batch_var  <- resolve_batch_variable(meta, expr, bc$batch_col, sample_col)
  group_var  <- resolve_group_variable(meta, expr, bc$group_col %||% cfg$effects$color, sample_col)
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

  pre$expr_filt        <- corrected
  pre$expr_work        <- corrected
  pre$expr_imp_single  <- corrected
  pre$batch_corrected  <- TRUE
  pre$batch_method     <- bc$method
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
#' @return A factor of group labels per column, or `NULL` if `group_col` is unset
#'   or absent.
#' @noRd
resolve_group_variable <- function(meta, expr, group_col, sample_col) {
  if (is.null(group_col) || !nzchar(group_col) || !group_col %in% colnames(meta)) {
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

  if (anyNA(expr)) {
    cli::cli_abort(c(
      "ComBat cannot run on a matrix containing NA.",
      "i" = "Run batch correction after imputation, or impute before this step."
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
