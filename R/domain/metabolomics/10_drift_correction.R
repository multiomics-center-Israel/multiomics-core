# R/domain/metabolomics/10_drift_correction.R
#
# LOESS-based instrument drift correction for metabolomics.
#
# Functions:
#   check_drift_preconditions()    — validate metadata requirements
#   correct_drift_loess()          — per-feature LOESS correction
#   apply_drift_correction()       — config-aware dispatcher
#   plot_drift_pca_comparison()    — 2-panel PCA before/after (pure)


#' Check whether LOESS drift correction can be applied
#'
#' @param meta data.frame of sample metadata.
#' @param injection_order_col Column name containing numeric injection order.
#' @param qc_flag_col Column name containing a logical QC flag (TRUE = QC sample).
#'
#' @return list with:
#'   \describe{
#'     \item{can_apply}{logical}
#'     \item{reason}{character explanation (empty string if can_apply = TRUE)}
#'   }
#'
check_drift_preconditions <- function(meta, injection_order_col, qc_flag_col) {
  meta <- as.data.frame(meta)

  if (!injection_order_col %in% colnames(meta)) {
    return(list(can_apply = FALSE,
                reason = sprintf("Column '%s' not found in metadata.", injection_order_col)))
  }
  if (!qc_flag_col %in% colnames(meta)) {
    return(list(can_apply = FALSE,
                reason = sprintf("Column '%s' not found in metadata.", qc_flag_col)))
  }

  order_vals <- suppressWarnings(as.numeric(meta[[injection_order_col]]))
  if (any(is.na(order_vals))) {
    return(list(can_apply = FALSE,
                reason = sprintf(
                  "Column '%s' contains NA or non-numeric values.",
                  injection_order_col
                )))
  }

  is_qc   <- as.logical(meta[[qc_flag_col]])
  n_qc    <- sum(is_qc, na.rm = TRUE)
  min_qc  <- 4L

  if (n_qc < min_qc) {
    return(list(can_apply = FALSE,
                reason = sprintf(
                  "Only %d QC sample(s) found (minimum %d required for stable LOESS).",
                  n_qc, min_qc
                )))
  }

  list(can_apply = TRUE, reason = "")
}


#' Per-feature LOESS drift correction
#'
#' Fits a LOESS curve through QC sample intensities along injection order,
#' then divides each sample's intensity by the predicted trend normalised to
#' the QC median.
#'
#' @param mat Numeric matrix (features x samples), log-scale intensities.
#' @param meta data.frame aligned to \code{colnames(mat)}.
#' @param sample_col Column in \code{meta} with sample identifiers.
#' @param injection_order_col Column with numeric injection order.
#' @param qc_flag_col Column with logical QC indicator.
#' @param loess_span Numeric LOESS span (default 0.75). Clamped to
#'   \code{max(loess_span, 4 / n_qc)} to prevent unstable fits.
#' @param epsilon Numeric guard; features where any predicted value <
#'   \code{epsilon} are skipped with a warning.
#'
#' @return list with:
#'   \describe{
#'     \item{mat_corrected}{Drift-corrected matrix.}
#'     \item{n_corrected}{Integer: features successfully corrected.}
#'     \item{n_skipped}{Integer: features skipped.}
#'     \item{skipped_ids}{Character vector of skipped feature IDs.}
#'   }
#'
correct_drift_loess <- function(mat, meta, sample_col,
                                injection_order_col,
                                qc_flag_col,
                                loess_span = 0.75,
                                epsilon    = 1e-6) {
  mat       <- as.matrix(mat)
  meta      <- as.data.frame(meta)
  order_all <- as.numeric(meta[[injection_order_col]])
  is_qc     <- as.logical(meta[[qc_flag_col]])

  n_qc       <- sum(is_qc, na.rm = TRUE)
  span_use   <- max(loess_span, 4 / n_qc)
  order_qc   <- order_all[is_qc]

  mat_corr   <- mat
  n_corrected <- 0L
  skipped_ids <- character(0)

  for (i in seq_len(nrow(mat))) {
    feat_id <- rownames(mat)[i]
    y_qc    <- mat[i, is_qc]

    if (any(is.na(y_qc))) {
      skipped_ids <- c(skipped_ids, feat_id)
      next
    }

    fit <- tryCatch(
      stats::loess(y_qc ~ order_qc, span = span_use),
      error = function(e) {
        warning(sprintf("correct_drift_loess: LOESS failed for '%s': %s", feat_id, e$message))
        NULL
      }
    )
    if (is.null(fit)) {
      skipped_ids <- c(skipped_ids, feat_id)
      next
    }

    predicted_all <- tryCatch(
      stats::predict(fit, newdata = data.frame(order_qc = order_all)),
      error = function(e) NULL
    )
    if (is.null(predicted_all) || any(is.na(predicted_all))) {
      skipped_ids <- c(skipped_ids, feat_id)
      next
    }

    if (any(predicted_all < epsilon)) {
      warning(sprintf(
        "correct_drift_loess: predicted values < epsilon for '%s'; skipping.",
        feat_id
      ))
      skipped_ids <- c(skipped_ids, feat_id)
      next
    }

    qc_median <- stats::median(predicted_all[is_qc], na.rm = TRUE)
    if (!is.finite(qc_median) || qc_median < epsilon) {
      skipped_ids <- c(skipped_ids, feat_id)
      next
    }

    correction_factor  <- predicted_all / qc_median
    mat_corr[i, ]      <- mat[i, ] / correction_factor
    n_corrected        <- n_corrected + 1L
  }

  list(
    mat_corrected = mat_corr,
    n_corrected   = n_corrected,
    n_skipped     = length(skipped_ids),
    skipped_ids   = skipped_ids
  )
}


#' Apply drift correction using pipeline configuration
#'
#' Reads drift correction settings from \code{cfg$preprocessing$drift_correction},
#' calls \code{check_drift_preconditions()}, then either applies
#' \code{correct_drift_loess()} or returns \code{mat} unchanged with a message.
#'
#' @param mat Numeric matrix (features x samples).
#' @param meta data.frame of sample metadata aligned to \code{colnames(mat)}.
#' @param cfg Mode config: \code{config$modes$metabolomics}.
#'
#' @return list with:
#'   \describe{
#'     \item{mat}{(Possibly corrected) numeric matrix.}
#'     \item{applied}{logical: TRUE if drift correction was applied.}
#'     \item{info}{list with correction details.}
#'   }
#'
apply_drift_correction <- function(mat, meta, cfg) {
  mat  <- as.matrix(mat)
  dc   <- cfg$preprocessing$drift_correction %||% list()
  enabled <- isTRUE(dc$enabled %||% TRUE)

  if (!enabled) {
    message("apply_drift_correction: disabled in config; skipping.")
    return(list(mat = mat, applied = FALSE,
                info = list(reason = "disabled in config")))
  }

  inj_col <- dc$injection_order_col %||% "injection_order"
  qc_col  <- dc$qc_flag_col         %||% "is_QC"
  span    <- dc$loess_span           %||% 0.75
  eps     <- dc$epsilon              %||% 1e-8

  pre <- check_drift_preconditions(meta, inj_col, qc_col)
  if (!pre$can_apply) {
    message(sprintf(
      "apply_drift_correction: skipped — %s", pre$reason
    ))
    return(list(mat = mat, applied = FALSE,
                info = list(reason = pre$reason)))
  }

  result <- correct_drift_loess(
    mat                 = mat,
    meta                = meta,
    sample_col          = cfg$effects$samples %||% "sample_id",
    injection_order_col = inj_col,
    qc_flag_col         = qc_col,
    loess_span          = span,
    epsilon             = eps
  )

  message(sprintf(
    "apply_drift_correction: corrected %d features, skipped %d.",
    result$n_corrected, result$n_skipped
  ))

  list(
    mat     = result$mat_corrected,
    applied = TRUE,
    info    = list(
      n_corrected  = result$n_corrected,
      n_skipped    = result$n_skipped,
      skipped_ids  = result$skipped_ids,
      loess_span   = span,
      inj_col      = inj_col,
      qc_col       = qc_col
    )
  )
}


#' 2-panel PCA comparison: before vs after drift correction
#'
#' PURE function: no config access, no file I/O.
#' Returns a patchwork/ggplot grid coloured by injection order (continuous).
#'
#' @param mat_before Numeric matrix before correction (features x samples).
#' @param mat_after  Numeric matrix after correction.
#' @param meta data.frame aligned to \code{colnames(mat_before)}.
#' @param injection_order_col Column in \code{meta} with numeric injection order.
#'
#' @return patchwork ggplot object (2 panels side by side).
#'
plot_drift_pca_comparison <- function(mat_before, mat_after, meta,
                                      injection_order_col = "injection_order") {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for plot_drift_pca_comparison().")
  }

  meta <- as.data.frame(meta)
  inj_order <- as.numeric(meta[[injection_order_col]])

  make_pca_df <- function(mat, label) {
    res    <- compute_pca_scores(mat, pcs = c(1, 2), center = TRUE, scale = FALSE)
    scores <- res$scores
    var_pct <- round(100 * res$var_expl[1:2], 1)
    scores$injection_order <- inj_order
    scores$.panel <- label
    scores$.pc1_label <- sprintf("PC1 (%.1f%%)", var_pct[1])
    scores$.pc2_label <- sprintf("PC2 (%.1f%%)", var_pct[2])
    scores
  }

  df_before <- make_pca_df(mat_before, "Before correction")
  df_after  <- make_pca_df(mat_after,  "After correction")

  make_panel <- function(df) {
    ggplot2::ggplot(df, ggplot2::aes(
      x = PC1, y = PC2, colour = injection_order
    )) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::scale_colour_viridis_c(name = "Injection\norder") +
      ggplot2::labs(
        title = unique(df$.panel),
        x     = unique(df$.pc1_label),
        y     = unique(df$.pc2_label)
      ) +
      ggplot2::theme_minimal(base_size = 12)
  }

  patchwork::wrap_plots(make_panel(df_before), make_panel(df_after), ncol = 2)
}
