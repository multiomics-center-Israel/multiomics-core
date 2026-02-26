# R/domain/metabolomics/09_imputation_met.R
#
# Targeted imputation for metabolomics:
#   MNAR features → min/2 per feature
#   MAR  features → KNN (impute::impute.knn)
#
# Functions:
#   impute_mnar_min_half()        — min/2 imputation for MNAR rows
#   impute_mar_knn()              — KNN imputation for MAR rows
#   impute_metabolomics()         — dispatcher


#' Impute MNAR features using half the per-feature minimum observed value
#'
#' For each feature classified as "MNAR", replaces NA values with
#' \code{min(observed) / 2}.  If the entire row is NA, falls back to
#' \code{global_min / 2} with a warning.
#'
#' @param mat Numeric matrix (features x samples), possibly with NAs.
#' @param mnar_class Character vector (length == nrow(mat)) with values
#'   "MNAR", "MAR", "all_missing", or "all_observed".  Names must match
#'   rownames(mat).
#'
#' @return list with:
#'   \describe{
#'     \item{mat}{Matrix with MNAR NAs replaced.}
#'     \item{n_imputed}{Integer: number of NA cells imputed.}
#'   }
#'
impute_mnar_min_half <- function(mat, mnar_class) {
  mat  <- as.matrix(mat)

  # Align mnar_class to row order of mat
  if (!is.null(names(mnar_class))) {
    mnar_class <- mnar_class[rownames(mat)]
  }

  is_mnar <- mnar_class == "MNAR"
  if (!any(is_mnar, na.rm = TRUE)) {
    return(list(mat = mat, n_imputed = 0L))
  }

  global_min <- min(mat, na.rm = TRUE)
  if (!is.finite(global_min)) global_min <- 0

  n_imputed <- 0L

  for (i in which(is_mnar)) {
    row_vals  <- mat[i, ]
    na_mask   <- is.na(row_vals)
    obs_vals  <- row_vals[!na_mask]

    fill_val <- if (length(obs_vals) > 0 && any(is.finite(obs_vals))) {
      min(obs_vals, na.rm = TRUE) / 2
    } else {
      warning(sprintf(
        "impute_mnar_min_half: feature '%s' is entirely NA; using global min/2.",
        rownames(mat)[i]
      ))
      global_min / 2
    }

    mat[i, na_mask] <- fill_val
    n_imputed       <- n_imputed + sum(na_mask)
  }

  list(mat = mat, n_imputed = n_imputed)
}


#' Impute MAR features using KNN
#'
#' Subsets the matrix to MAR rows, runs \code{impute::impute.knn()}, and
#' writes imputed values back.
#'
#' Guards:
#' \itemize{
#'   \item All-NA MAR rows are skipped (logged; they should not appear after
#'         \code{filter_by_missingness()} with a 20\% feature threshold).
#'   \item If \code{k >= ncol(mat)}, \code{k} is clamped to
#'         \code{ncol(mat) - 1}.  If that yields \code{k < 1}, column-mean
#'         fallback is used with a warning.
#' }
#'
#' @param mat Numeric matrix (features x samples), possibly with NAs.
#' @param mar_mask Logical vector (length == nrow(mat)); TRUE = MAR feature.
#' @param k Integer KNN neighbour count (default 10).
#'
#' @return list with:
#'   \describe{
#'     \item{mat}{Matrix with MAR NAs imputed.}
#'     \item{n_imputed}{Integer: cells imputed.}
#'     \item{n_skipped}{Integer: all-NA MAR rows skipped.}
#'   }
#'
impute_mar_knn <- function(mat, mar_mask, k = 10L) {
  if (!requireNamespace("impute", quietly = TRUE)) {
    stop("Package 'impute' (Bioconductor) is required for KNN imputation.")
  }

  mat <- as.matrix(mat)
  k   <- as.integer(k)

  # Align mar_mask to mat rows
  if (!is.null(names(mar_mask))) {
    mar_mask <- mar_mask[rownames(mat)]
  }

  mar_idx <- which(mar_mask)
  if (length(mar_idx) == 0L) {
    return(list(mat = mat, n_imputed = 0L, n_skipped = 0L))
  }

  # Guard: k must be < ncol(mat)
  k_eff <- min(k, ncol(mat) - 1L)

  if (k_eff < 1L) {
    warning(
      "impute_mar_knn: k < 1 after clamping (ncol(mat) = ", ncol(mat), "). ",
      "Falling back to column-mean imputation for MAR features."
    )
    n_imputed <- 0L
    for (i in mar_idx) {
      na_mask <- is.na(mat[i, ])
      if (any(na_mask)) {
        col_means        <- colMeans(mat, na.rm = TRUE)
        mat[i, na_mask]  <- col_means[na_mask]
        n_imputed        <- n_imputed + sum(na_mask)
      }
    }
    return(list(mat = mat, n_imputed = n_imputed, n_skipped = 0L))
  }

  if (k_eff < k) {
    message(sprintf(
      "impute_mar_knn: k clamped from %d to %d (ncol = %d).",
      k, k_eff, ncol(mat)
    ))
  }

  # Identify all-NA MAR rows (should not exist post-filter, but guard anyway)
  all_na_mar <- mar_idx[apply(mat[mar_idx, , drop = FALSE], 1,
                              function(r) all(is.na(r)))]
  run_idx    <- setdiff(mar_idx, all_na_mar)

  n_skipped <- length(all_na_mar)
  if (n_skipped > 0L) {
    warning(sprintf(
      "impute_mar_knn: %d all-NA MAR feature(s) skipped (still NA after imputation).",
      n_skipped
    ))
  }

  if (length(run_idx) == 0L) {
    return(list(mat = mat, n_imputed = 0L, n_skipped = n_skipped))
  }

  sub_mat <- mat[run_idx, , drop = FALSE]
  n_before <- sum(is.na(sub_mat))

  imp_result <- impute::impute.knn(sub_mat, k = k_eff)
  mat[run_idx, ] <- imp_result$data

  list(
    mat       = mat,
    n_imputed = n_before - sum(is.na(mat[run_idx, , drop = FALSE])),
    n_skipped = n_skipped
  )
}


#' Targeted metabolomics imputation dispatcher (MNAR → min/2, MAR → KNN)
#'
#' Classifies features using \code{miss_stats} (data.frame from
#' \code{classify_missingness()}) or a logical MNAR mask, then applies
#' MNAR min/2 first, followed by KNN for MAR features.
#'
#' @param mat Numeric matrix (features x samples).
#' @param miss_stats data.frame with columns \code{feature_id} and
#'   \code{mnar_class} ("MNAR", "MAR", "all_missing", "all_observed"), as
#'   returned by \code{classify_missingness()}.  Alternatively, a named
#'   logical vector (TRUE = MNAR).
#' @param knn_k Integer KNN k (default 10).
#' @param epsilon Numeric guard; imputed values below this are floored to
#'   \code{epsilon} (default 1e-8).
#'
#' @return list with:
#'   \describe{
#'     \item{mat_imputed}{Fully imputed matrix.}
#'     \item{imputed_flag}{Logical matrix: TRUE where a value was imputed.}
#'     \item{info}{list: n_mnar, n_mar, n_skipped, n_all_missing, knn_k}
#'   }
#'
impute_metabolomics <- function(mat, miss_stats, knn_k = 10L, epsilon = 1e-8) {
  mat_orig <- as.matrix(mat)
  na_orig  <- is.na(mat_orig)

  # Build mnar_class character vector aligned to mat rows
  if (is.data.frame(miss_stats)) {
    mnar_class <- stats::setNames(miss_stats$mnar_class, miss_stats$feature_id)
    mnar_class <- mnar_class[rownames(mat_orig)]
  } else {
    # Named logical vector (TRUE = MNAR)
    lv         <- miss_stats[rownames(mat_orig)]
    mnar_class <- ifelse(lv, "MNAR", "MAR")
    names(mnar_class) <- rownames(mat_orig)
  }

  n_mnar       <- sum(mnar_class == "MNAR", na.rm = TRUE)
  n_mar        <- sum(mnar_class == "MAR",  na.rm = TRUE)
  n_all_missing <- sum(mnar_class == "all_missing", na.rm = TRUE)

  # Step 1: MNAR → min/2
  mnar_result <- impute_mnar_min_half(mat_orig, mnar_class)
  mat_work    <- mnar_result$mat

  # Step 2: MAR → KNN
  mar_mask    <- mnar_class == "MAR"
  mar_result  <- impute_mar_knn(mat_work, mar_mask, k = knn_k)
  mat_work    <- mar_result$mat

  # Step 3: epsilon floor for imputed values
  imputed_flag <- na_orig & !is.na(mat_work)
  if (epsilon > 0 && any(imputed_flag, na.rm = TRUE)) {
    below_eps <- imputed_flag & (mat_work < epsilon)
    if (any(below_eps, na.rm = TRUE)) {
      mat_work[below_eps] <- epsilon
    }
  }

  list(
    mat_imputed  = mat_work,
    imputed_flag = imputed_flag,
    info = list(
      n_mnar       = n_mnar,
      n_mar        = n_mar,
      n_skipped    = mar_result$n_skipped,
      n_all_missing = n_all_missing,
      knn_k        = knn_k
    )
  )
}
