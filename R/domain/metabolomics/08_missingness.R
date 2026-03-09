# R/domain/metabolomics/08_missingness.R
#
# Missingness classification and filtering for metabolomics preprocessing.
#
# Functions:
#   classify_missingness()       — per-feature MNAR/MAR classification
#   compute_missingness_stats()  — full stats table + pass/fail lists
#   filter_by_missingness()      — threshold-based feature + sample removal


#' Classify missingness mechanism per feature (MNAR vs MAR)
#'
#' Uses Spearman correlation between a feature's missingness indicator and
#' sample intensity rank.  Negative correlation indicates that values are
#' missing in low-intensity samples (Not At Random, MNAR).
#'
#' @param mat Numeric matrix (features x samples).
#' @param mnar_cor_threshold Numeric scalar (default -0.3). Features with
#'   Spearman correlation < \code{mnar_cor_threshold} are classified as MNAR.
#'   Note the sign: the threshold is negative.
#'
#' @return data.frame with one row per feature and columns:
#'   \describe{
#'     \item{feature_id}{rowname of \code{mat}}
#'     \item{n_missing}{integer count of missing values}
#'     \item{pct_missing}{fraction of samples missing (0–1)}
#'     \item{spearman_cor}{Spearman correlation of missingness indicator with
#'       sample intensity rank; NA for fully observed/missing features}
#'     \item{mnar_class}{character: "MNAR", "MAR", "all_missing", or
#'       "all_observed"}
#'   }
#'
classify_missingness <- function(mat, mnar_cor_threshold = -0.3) {
  mat <- as.matrix(mat)
  n_samp <- ncol(mat)

  # Sample intensity rank (rank of average observed intensity across samples)
  sample_intensity_rank <- rank(colMeans(mat, na.rm = TRUE),
                                ties.method = "average")

  result <- lapply(seq_len(nrow(mat)), function(i) {
    row_vals <- mat[i, ]
    miss_ind <- as.integer(is.na(row_vals))  # 1 = missing, 0 = observed
    n_miss   <- sum(miss_ind)

    mnar_class <- if (n_miss == 0L) {
      "all_observed"
    } else if (n_miss == n_samp) {
      "all_missing"
    } else {
      cor_val <- stats::cor(miss_ind, sample_intensity_rank,
                            method = "spearman", use = "complete.obs")
      if (!is.na(cor_val) && cor_val < mnar_cor_threshold) "MNAR" else "MAR"
    }

    cor_val_out <- if (n_miss == 0L || n_miss == n_samp) NA_real_ else {
      stats::cor(miss_ind, sample_intensity_rank,
                 method = "spearman", use = "complete.obs")
    }

    data.frame(
      feature_id   = rownames(mat)[i],
      n_missing    = n_miss,
      pct_missing  = n_miss / n_samp,
      spearman_cor = cor_val_out,
      mnar_class   = mnar_class,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, result)
}


#' Compute full missingness statistics for a metabolomics matrix
#'
#' Calls \code{classify_missingness()} and also computes per-sample missingness,
#' identifying features and samples that fail the given thresholds.
#'
#' @param mat Numeric matrix (features x samples).
#' @param meta data.frame of sample metadata (optional; used only for alignment
#'   in the returned object).
#' @param sample_col Column in \code{meta} matching \code{colnames(mat)}.
#' @param feat_miss_threshold Numeric (default 0.20).  Features with
#'   \code{pct_missing > feat_miss_threshold} are flagged.
#' @param samp_miss_threshold Numeric (default 0.30).  Samples with fraction
#'   missing > \code{samp_miss_threshold} are flagged.
#' @param mnar_cor_threshold Numeric (default -0.3).  Passed to
#'   \code{classify_missingness()}.
#'
#' @return list with:
#'   \describe{
#'     \item{stats_df}{data.frame from \code{classify_missingness()} with an
#'       extra column \code{fails_threshold} (logical)}
#'     \item{samp_miss_df}{data.frame: sample_id, n_missing, pct_missing,
#'       fails_threshold}
#'     \item{feat_pass}{character vector of feature IDs passing the threshold}
#'     \item{samp_pass}{character vector of sample IDs passing the threshold}
#'     \item{n_feat_removed}{integer}
#'     \item{n_samp_removed}{integer}
#'   }
#'
compute_missingness_stats <- function(mat,
                                      meta = NULL,
                                      sample_col = "sample_id",
                                      feat_miss_threshold = 0.20,
                                      samp_miss_threshold = 0.30,
                                      mnar_cor_threshold  = -0.3) {
  mat <- as.matrix(mat)

  stats_df <- classify_missingness(mat, mnar_cor_threshold)
  stats_df$fails_threshold <- stats_df$pct_missing > feat_miss_threshold

  # Per-sample missingness
  samp_miss_pct <- colMeans(is.na(mat))
  samp_miss_df  <- data.frame(
    sample_id       = colnames(mat),
    n_missing       = colSums(is.na(mat)),
    pct_missing     = samp_miss_pct,
    fails_threshold = samp_miss_pct > samp_miss_threshold,
    stringsAsFactors = FALSE
  )
  colnames(samp_miss_df)[1] <- sample_col

  feat_pass <- stats_df$feature_id[!stats_df$fails_threshold]
  samp_pass <- samp_miss_df[[sample_col]][!samp_miss_df$fails_threshold]

  list(
    stats_df       = stats_df,
    samp_miss_df   = samp_miss_df,
    feat_pass      = feat_pass,
    samp_pass      = samp_pass,
    n_feat_removed = sum(stats_df$fails_threshold),
    n_samp_removed = sum(samp_miss_df$fails_threshold)
  )
}


#' Filter features and samples by missingness thresholds
#'
#' Applies two sequential filters:
#' 1. Drop samples where the fraction of missing values > \code{sample_threshold}.
#' 2. Drop features where the fraction of missing values (on remaining samples)
#'    > \code{feat_threshold}.
#'
#' Samples are filtered first so that feature missingness is evaluated on the
#' retained sample set.
#'
#' @param mat Numeric matrix (features x samples).
#' @param meta data.frame with a column matching \code{sample_col}.
#' @param sample_col Column in \code{meta} containing sample identifiers that
#'   match \code{colnames(mat)}.
#' @param feat_threshold Numeric (default 0.20).
#' @param sample_threshold Numeric (default 0.30).
#'
#' @return list with:
#'   \describe{
#'     \item{mat}{Filtered numeric matrix}
#'     \item{meta}{Metadata filtered to retained samples}
#'     \item{dropped_features}{character vector of removed feature IDs}
#'     \item{dropped_samples}{character vector of removed sample IDs}
#'   }
#'
filter_by_missingness <- function(mat, meta, sample_col,
                                  feat_threshold   = 0.20,
                                  sample_threshold = 0.30) {
  mat  <- as.matrix(mat)
  meta <- as.data.frame(meta)

  # 1. Filter samples first
  samp_miss_pct    <- colMeans(is.na(mat))
  keep_samps       <- colnames(mat)[samp_miss_pct <= sample_threshold]
  dropped_samples  <- setdiff(colnames(mat), keep_samps)

  if (length(dropped_samples) > 0) {
    message(sprintf(
      "filter_by_missingness: dropping %d sample(s) with > %.0f%% missing: %s",
      length(dropped_samples),
      sample_threshold * 100,
      paste(head(dropped_samples, 5), collapse = ", "),
      if (length(dropped_samples) > 5) "..." else ""
    ))
  }

  mat  <- mat[, keep_samps, drop = FALSE]
  meta <- meta[meta[[sample_col]] %in% keep_samps, , drop = FALSE]

  # 2. Filter features on the retained sample set
  feat_miss_pct   <- rowMeans(is.na(mat))
  keep_feats      <- rownames(mat)[feat_miss_pct <= feat_threshold]
  dropped_feats   <- setdiff(rownames(mat), keep_feats)

  if (length(dropped_feats) > 0) {
    message(sprintf(
      "filter_by_missingness: dropping %d feature(s) with > %.0f%% missing.",
      length(dropped_feats),
      feat_threshold * 100
    ))
  }

  mat <- mat[keep_feats, , drop = FALSE]

  message(sprintf(
    "filter_by_missingness: retained %d features x %d samples.",
    nrow(mat), ncol(mat)
  ))

  list(
    mat             = mat,
    meta            = meta,
    dropped_features = dropped_feats,
    dropped_samples  = dropped_samples
  )
}
