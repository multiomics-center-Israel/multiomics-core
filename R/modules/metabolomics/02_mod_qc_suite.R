# R/modules/metabolomics/02_mod_qc_suite.R
#
# QC suite module for the metabolomics preprocessing pipeline.
#
# Provides:
#   compute_linear_rsd()  — authoritative RSD computation with scale-aware
#                           back-transformation from Log2 to linear.
#
# Future additions (Phase 2+):
#   qc_full_metabolomics_suite()
#   mod_met_qc_suite()


# ==============================================================================
# compute_linear_rsd — scale-aware RSD on back-transformed linear intensities
# ==============================================================================

#' Compute median RSD after back-transforming Log2 values to linear scale
#'
#' This function is the \strong{authoritative source} for Relative Standard
#' Deviation (RSD) computation in the metabolomics QC pipeline.  All RSD
#' values reported in QC metrics TSVs and the normalization comparison table
#' must be produced by this function.
#'
#' \strong{Why back-transformation is required:}
#' RSD is defined as \eqn{\mathrm{SD}(x) / |\mathrm{mean}(x)|} for linear
#' intensities.  Computing this ratio directly on log2 values yields the
#' coefficient of variation of the log-transformed data, which has no standard
#' interpretive meaning in metabolomics QC and must not be used.
#'
#' \strong{Scale Contract alignment:}
#' Every Log2-output target in the pipeline (\code{met_log}, \code{met_norm_tss},
#' \code{met_norm_pqn}, \code{met_norm_median}) used a transformation of the
#' form \eqn{\log_2(x + p)} where \eqn{p} is the pseudocount.  The exact
#' inverse is \eqn{x = 2^{\text{val}} - p}.  For \code{met_norm_median} this
#' inversion is approximate (see \code{backtransform_exact} in the return
#' value); for all other stages it is exact.
#'
#' @param mat         Numeric matrix (features \eqn{\times} samples) on the
#'   \strong{Log2} scale.  Must have at least 1 row and 2 columns for
#'   meaningful RSD computation.
#' @param stage       Character scalar identifying the preprocessing stage.
#'   Must be one of \code{"log"}, \code{"norm_tss"}, \code{"norm_pqn"}, or
#'   \code{"norm_median"}.  Used in warning messages and to set the
#'   \code{backtransform_exact} flag.
#' @param pseudocount Numeric scalar \eqn{\geq 0}.  The pseudocount added
#'   before the log2 transform in \code{transform_metab()} (default 1).
#'   \strong{Must be the same value used when the matrix was created.}
#'   Do not infer this from the data.
#'
#' @return Named list with:
#'   \describe{
#'     \item{\code{median_rsd}}{Numeric scalar.  Median RSD across all
#'       features (NA if no finite RSD values could be computed).}
#'     \item{\code{rsd_per_feature}}{Named numeric vector of per-feature RSD
#'       values (NA for features with zero or non-finite mean after
#'       back-transformation).}
#'     \item{\code{n_features}}{Integer.  Number of rows in \code{mat}.}
#'     \item{\code{n_samples}}{Integer.  Number of columns in \code{mat}.}
#'     \item{\code{backtransform_exact}}{Logical.  \code{TRUE} for stages
#'       \code{"log"}, \code{"norm_tss"}, \code{"norm_pqn"} where
#'       \eqn{2^{\text{val}} - p} exactly reverses the log2 transform.
#'       \code{FALSE} for \code{"norm_median"} where the log-shift
#'       operation does not preserve the clean \eqn{\log_2(x + p)}
#'       relationship, making the inversion approximate.}
#'     \item{\code{pseudocount_used}}{Numeric.  The pseudocount value
#'       supplied; echoed back for audit purposes.}
#'   }
#'
#' @details
#' \strong{Back-transformation:}
#' \deqn{\text{linear\_val}_{ij} = 2^{\text{mat}_{ij}} - p}
#' where \eqn{p} is the pseudocount.
#'
#' \strong{Negative-value guard:}
#' Theoretically, \eqn{2^{\text{val}} \geq p} for any value produced by
#' \eqn{\log_2(x + p)} with \eqn{x \geq 0}.  In practice, numerical
#' imprecision or approximate back-transformation (median stage) can yield
#' small negative values.  Any cell where \eqn{2^{\text{val}} - p < 0} is
#' floored to \eqn{\varepsilon = 2.2 \times 10^{-16}}
#' (\code{.Machine$double.eps}) and a warning is issued.
#'
#' \strong{Per-feature RSD:}
#' \deqn{\mathrm{RSD}_i = \frac{\mathrm{SD}_j(\text{linear\_val}_{ij})}{|\mathrm{mean}_j(\text{linear\_val}_{ij})|}  }
#' Features where the absolute mean is 0 or non-finite receive \code{NA}.
#'
#' \strong{Median RSD:}
#' \deqn{\mathrm{median\_rsd} = \mathrm{median}_i(\mathrm{RSD}_i)}
#' computed with \code{na.rm = TRUE}.
compute_linear_rsd <- function(mat, stage, pseudocount = 1) {
  # ---- Input validation -------------------------------------------------------
  valid_stages <- c("log", "norm_tss", "norm_pqn", "norm_median")
  if (!stage %in% valid_stages) {
    warning(sprintf(
      "compute_linear_rsd: unknown stage '%s'; expected one of: %s. ",
      stage, paste(valid_stages, collapse = ", ")
    ))
  }
  if (!is.numeric(pseudocount) || length(pseudocount) != 1 || pseudocount < 0) {
    stop("compute_linear_rsd: 'pseudocount' must be a single non-negative number.")
  }
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("compute_linear_rsd: 'mat' must be a numeric matrix.")
  }
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    warning(sprintf(
      "compute_linear_rsd [stage='%s']: matrix has 0 rows or 0 columns; returning NA.", stage
    ))
    return(list(
      median_rsd          = NA_real_,
      rsd_per_feature     = setNames(numeric(0), character(0)),
      n_features          = nrow(mat),
      n_samples           = ncol(mat),
      backtransform_exact = stage %in% c("log", "norm_tss", "norm_pqn"),
      pseudocount_used    = pseudocount
    ))
  }

  # ---- Back-transformation: linear_val = 2^(log2_val) - pseudocount ----------
  #
  # This exactly inverts transform_metab(mat, method = "log2", pseudocount = p),
  # which computes log2(x + p).  The inverse is 2^(val) - p = x (linear).
  # For met_norm_median, the log-shift does not preserve the log2(x + p) form
  # exactly, so this inversion is approximate for that stage.
  mat_linear <- 2^mat - pseudocount

  # ---- Negative-value guard --------------------------------------------------
  #
  # Negative linear values are physically impossible for peak intensities and
  # indicate numerical imprecision or approximation error.  Floor to epsilon.
  epsilon <- .Machine$double.eps
  n_neg   <- sum(mat_linear < 0, na.rm = TRUE)
  if (n_neg > 0) {
    warning(sprintf(
      "compute_linear_rsd [stage='%s']: %d cell(s) produced negative values after ",
      "back-transform (2^val - %g); floored to %g.  If stage is 'norm_median', ",
      "this is expected (approximate back-transform).",
      stage, n_neg, pseudocount, epsilon
    ))
    mat_linear[mat_linear < 0] <- epsilon
  }

  # ---- Per-feature RSD: SD / |mean| across samples ---------------------------
  feat_sds      <- apply(mat_linear, 1, stats::sd, na.rm = TRUE)
  feat_means    <- rowMeans(mat_linear, na.rm = TRUE)
  feat_means_ab <- abs(feat_means)

  # Guard: zero or non-finite absolute mean → RSD undefined → NA
  feat_means_ab[feat_means_ab == 0 | !is.finite(feat_means_ab)] <- NA_real_

  rsd <- feat_sds / feat_means_ab
  rsd[!is.finite(rsd)] <- NA_real_

  # Preserve feature names if present
  if (!is.null(rownames(mat))) {
    names(rsd) <- rownames(mat)
  }

  # ---- Aggregate metric ------------------------------------------------------
  median_rsd <- stats::median(rsd, na.rm = TRUE)

  list(
    median_rsd          = median_rsd,
    rsd_per_feature     = rsd,
    n_features          = nrow(mat_linear),
    n_samples           = ncol(mat_linear),
    backtransform_exact = stage %in% c("log", "norm_tss", "norm_pqn"),
    pseudocount_used    = pseudocount
  )
}
