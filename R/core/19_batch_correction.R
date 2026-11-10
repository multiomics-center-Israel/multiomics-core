#' Batch correction of a log-scale expression matrix
#'
#' Shared by every omics mode. Each mode resolves its own config into the
#' vectors below (batch labels, group labels, covariates) and calls one of the
#' two backends:
#'
#'   * [correct_batch_combat()] -- `sva::ComBat`, empirical Bayes. Adjusts each
#'     batch's mean and, unless `mean_only = TRUE`, its spread.
#'   * [correct_batch_limma()]  -- `limma::removeBatchEffect`. Adjusts each
#'     batch's mean only; accepts a second batch factor.
#'
#' Both take the biological group (and any covariates) into account, so the
#' correction is estimated after the group effect rather than blind to it.
#' Whether batch and biology are actually separable depends on the design:
#' both backends stop when they are fully confounded
#' ([assert_batch_not_confounded()]), and a partly confounded design still needs
#' to be checked by eye on the before/after plots.
#'
#' The input must be a continuous, log-scale matrix with no missing or
#' non-finite values (never raw counts). The corrected matrix is for display
#' and clustering; callers keep it apart from the DE input.
#'
#' These functions are pure: matrix in, matrix out, no file I/O.

#' Build the model matrix that carries the biology to keep
#'
#' @param group_var Factor of group labels per sample, or `NULL`.
#' @param covariates Data.frame of extra columns to keep (one row per sample),
#'   or `NULL`.
#' @return A model matrix with an intercept, or `NULL` when there is nothing to
#'   keep.
#' @noRd
batch_keep_design <- function(group_var = NULL, covariates = NULL) {
  if (is.null(group_var) && is.null(covariates)) return(NULL)
  n <- if (!is.null(group_var)) length(group_var) else nrow(covariates)
  md <- data.frame(.row = seq_len(n))
  if (!is.null(group_var)) md[[".group"]] <- group_var
  if (!is.null(covariates)) md <- cbind(md, covariates)
  md$.row <- NULL
  stats::model.matrix(~., data = md)
}

#' Stop when batch cannot be separated from the biology to keep
#'
#' Builds the design of the kept biology (group and covariates) plus the batch
#' factor(s) and requires it to be of full rank. A rank-deficient design means
#' at least one batch effect is indistinguishable from a group or covariate
#' effect -- typically every batch holds a single group -- and no correction
#' can remove one without the other.
#'
#' Only full confounding is caught. A partly confounded design (batches that
#' hold uneven shares of the groups) passes and must be judged on the
#' before/after plots.
#'
#' @param batch_var Factor of batch labels per sample.
#' @param group_var Factor of group labels per sample, or `NULL`.
#' @param batch2 Optional second batch factor, or `NULL`.
#' @param covariates Optional data.frame of extra kept columns, or `NULL`.
#' @return `invisible(TRUE)` when the design is estimable; otherwise an error.
assert_batch_not_confounded <- function(batch_var, group_var = NULL,
                                        batch2 = NULL, covariates = NULL) {
  keep <- batch_keep_design(group_var, covariates)
  if (is.null(keep)) keep <- matrix(1, nrow = length(batch_var), ncol = 1)

  batch_cols <- function(b) {
    b <- as.factor(b)
    if (nlevels(b) < 2) return(NULL)
    stats::model.matrix(~b)[, -1, drop = FALSE]
  }
  x <- cbind(keep, batch_cols(batch_var), batch_cols(batch2))

  if (qr(x)$rank < ncol(x)) {
    per_batch <- if (!is.null(group_var)) {
      tab <- table(as.factor(batch_var), as.factor(group_var))
      paste(vapply(rownames(tab), function(b) {
        sprintf("%s: %s", b, paste(colnames(tab)[tab[b, ] > 0], collapse = ", "))
      }, character(1)), collapse = "; ")
    }
    cli::cli_abort(c(
      "Batch and the biology to keep are confounded, so batch correction cannot separate them.",
      "i" = if (!is.null(per_batch)) "Groups in each batch: {per_batch}.",
      "i" = "Any correction would remove the group effect together with the batch effect.",
      "i" = "Turn batch correction off, or check that the batch and group columns are the intended ones."
    ))
  }
  invisible(TRUE)
}

#' Stop unless the matrix is complete and finite
#' @noRd
assert_batch_input_finite <- function(expr, method) {
  if (!all(is.finite(expr))) {
    cli::cli_abort(c(
      "{method} batch correction requires a matrix of finite values.",
      "i" = "Found {sum(is.na(expr))} NA and {sum(is.infinite(expr))} non-finite value(s) (e.g. from log2 of zero).",
      "i" = "Impute NAs and floor/remove non-finite values before batch correction."
    ))
  }
}

#' ComBat batch correction for a log-scale matrix
#'
#' Wraps `sva::ComBat` (parametric empirical Bayes). When `preserve_group` is
#' `TRUE` and a group vector or covariates are given, they enter ComBat's model
#' matrix, so the correction is estimated after their effect.
#'
#' @param expr Numeric log-scale matrix (features x samples), complete and
#'   finite.
#' @param batch_var Factor of batch labels (length `ncol(expr)`).
#' @param group_var Factor of biological group labels, or `NULL`.
#' @param covariates Optional data.frame of additional covariates, or `NULL`.
#' @param preserve_group Logical; put `group_var`/`covariates` in the ComBat
#'   model.
#' @param mean_only Logical, passed to ComBat as `mean.only`. `TRUE` adjusts
#'   only the batch means, not their spread; use it when batches are small
#'   (e.g. two samples per batch), where the per-batch spread is unreliable.
#' @param seed Integer seed passed through `withr::with_seed()`.
#' @return Corrected numeric matrix with the same dim/dimnames as `expr`.
#' @examples
#' \dontrun{
#' expr <- matrix(rnorm(40, 10), nrow = 5,
#'                dimnames = list(paste0("F", 1:5), paste0("S", 1:8)))
#' correct_batch_combat(expr, batch_var = factor(rep(c("B1", "B2"), each = 4)),
#'                      group_var = factor(rep(c("A", "B"), times = 4)))
#' }
correct_batch_combat <- function(expr, batch_var, group_var = NULL,
                                 covariates = NULL, preserve_group = TRUE,
                                 mean_only = FALSE, seed = 42) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg sva} is required for {.val combat} batch correction but is not installed.",
      "i" = "Install it with {.code BiocManager::install('sva')} (a dependency change must be approved by a human)."
    ))
  }
  assert_batch_input_finite(expr, "ComBat")

  mod <- if (isTRUE(preserve_group)) batch_keep_design(group_var, covariates)
  if (!is.null(mod)) {
    assert_batch_not_confounded(batch_var, group_var = group_var,
                                covariates = covariates)
  }

  withr::with_seed(seed, {
    sva::ComBat(
      dat       = as.matrix(expr),
      batch     = batch_var,
      mod       = mod,
      par.prior = TRUE,
      mean.only = isTRUE(mean_only)
    )
  })
}

#' limma batch correction for a log-scale matrix
#'
#' Wraps `limma::removeBatchEffect`, which fits each feature on the kept
#' biology plus the batch factor(s) and subtracts the batch terms. The group
#' and covariates go in as `design`, so the batch shift is estimated after the
#' group effect rather than from raw batch means; without that, a batch that
#' holds more of one group would also lose part of that group's effect.
#' Deterministic: no seed.
#'
#' @param expr Numeric log-scale matrix (features x samples), complete and
#'   finite.
#' @param batch_var Factor of batch labels (length `ncol(expr)`).
#' @param batch2 Optional second batch factor, or `NULL`.
#' @param group_var Factor of biological group labels, or `NULL`.
#' @param covariates Optional data.frame of additional covariates, or `NULL`.
#' @return Corrected numeric matrix with the same dim/dimnames as `expr`.
#' @examples
#' \dontrun{
#' expr <- matrix(rnorm(40, 10), nrow = 5,
#'                dimnames = list(paste0("F", 1:5), paste0("S", 1:8)))
#' correct_batch_limma(expr, batch_var = factor(rep(c("B1", "B2"), each = 4)),
#'                     group_var = factor(rep(c("A", "B"), times = 4)))
#' }
correct_batch_limma <- function(expr, batch_var, batch2 = NULL,
                                group_var = NULL, covariates = NULL) {
  assert_batch_input_finite(expr, "limma")
  assert_batch_not_confounded(batch_var, group_var = group_var,
                              batch2 = batch2, covariates = covariates)

  design <- batch_keep_design(group_var, covariates)
  if (is.null(design)) design <- matrix(1, nrow = ncol(expr), ncol = 1)

  limma::removeBatchEffect(
    as.matrix(expr),
    batch  = batch_var,
    batch2 = batch2,
    design = design
  )
}
