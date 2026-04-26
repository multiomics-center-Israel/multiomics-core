# R/domain/metabolomics/01_normalization.R
#
# Three-step normalization for metabolomics data.
# Three independent steps (applied in order):
#   1. Sample normalization  (sum, median, PQN, internal standard)
#   2. Transformation        (log2, log10)
#   3. Scaling               (auto, pareto, range)
#
# Each function operates on a numeric matrix [features x samples].
# Reuses: assert_numeric_matrix


# ==== 1. SAMPLE NORMALIZATION ================================================

#' Apply sample normalization
#'
#' @param mat     Numeric matrix (features x samples).
#' @param method  One of "none", "sum", "median", "pqn", "is".
#' @param ref_col Character: column in row_data that flags the internal standard
#'                (only used when method = "is").
#' @param row_data data.frame with at least a column matching ref_col (for IS).
#' @return Normalized numeric matrix (same dimensions).
normalize_samples <- function(mat, method = "none", ref_col = NULL, row_data = NULL) {
    method <- tolower(method)
    assert_one_of(method, "sample_norm",
                  c("none", "sum", "median", "pqn", "is"))

    switch(method,
        none   = mat,
        sum    = norm_total_sum(mat),
        median = norm_median(mat),
        pqn    = norm_pqn(mat),
        is     = norm_internal_standard(mat, ref_col, row_data),
        stop("Unknown sample normalization method: ", method)
    )
}


#' Total sum normalization
#' Scale each sample so that column sums are equal (to the median column sum).
norm_total_sum <- function(mat) {
    col_sums <- colSums(mat, na.rm = TRUE)
    if (any(col_sums == 0)) warning("Some samples have zero total intensity.")
    target <- stats::median(col_sums[col_sums > 0])
    scale_factors <- target / col_sums
    scale_factors[!is.finite(scale_factors)] <- 1
    sweep(mat, 2, scale_factors, FUN = "*")
}


#' Median normalization
#' Scale each sample by its median intensity (relative to the global median of medians).
norm_median <- function(mat) {
    col_medians <- apply(mat, 2, stats::median, na.rm = TRUE)
    if (any(col_medians == 0 | is.na(col_medians))) {
        warning("Some samples have zero/NA median; they will not be scaled.")
    }
    target <- stats::median(col_medians[col_medians > 0 & is.finite(col_medians)])
    scale_factors <- target / col_medians
    scale_factors[!is.finite(scale_factors)] <- 1
    sweep(mat, 2, scale_factors, FUN = "*")
}

#' Probabilistic Quotient Normalization (PQN)
#'
#' Applies PQN after an initial median normalization step.
#' The reference spectrum is defined as the feature-wise median across samples.
#'
#' Procedure:
#' 1. Apply median normalization.
#' 2. Compute the reference spectrum as the feature-wise median across samples.
#' 3. For each sample, calculate quotients relative to the reference.
#' 4. Use the sample's median quotient as a scaling factor.
#' 5. Divide each sample by this factor.
#'
#' Features with zero, missing, or non-finite reference values are excluded
#' from quotient calculation.
norm_pqn <- function(mat, qc_idx = NULL) {
  # Step 1: initial median normalization
  mat_mn <- norm_median(mat)
  
  # Step 2: reference spectrum
  if (!is.null(qc_idx)) {
    # Use QC samples as reference
    ref <- apply(mat_mn[, qc_idx, drop = FALSE], 1, stats::median, na.rm = TRUE)
  } else {
    # Fallback: use all samples
    ref <- apply(mat_mn, 1, stats::median, na.rm = TRUE)
  }
  
  # Exclude features with invalid reference values
  ref_ok <- is.finite(ref) & abs(ref) > .Machine$double.eps
  if (!any(ref_ok)) {
    warning("norm_pqn: all reference values are zero, missing, or non-finite; returning median-normalized data.")
    return(mat_mn)
  }
  
  # Step 3: compute sample/reference quotients
  quotients <- sweep(mat_mn[ref_ok, , drop = FALSE], 1, ref[ref_ok], FUN = "/")
  
  # Step 4: estimate sample-specific scaling factors
  med_quotients <- apply(quotients, 2, stats::median, na.rm = TRUE)
  bad <- !is.finite(med_quotients) | med_quotients == 0
  if (any(bad)) {
    warning(sprintf(
      "norm_pqn: %d sample(s) had invalid median quotient; using scale factor 1.",
      sum(bad)
    ))
    med_quotients[bad] <- 1
  }
  
  # Step 5: apply PQN scaling
  sweep(mat_mn, 2, med_quotients, FUN = "/")
}


#' Internal standard normalization
#'
#' Divides each sample's values by the intensity of the internal standard(s)
#' in that sample (using the median if multiple IS features).
#'
#' @param mat      Numeric matrix (features x samples).
#' @param ref_col  Column name in row_data that flags IS rows (logical or "IS" string).
#' @param row_data data.frame (same number of rows as mat).
norm_internal_standard <- function(mat, ref_col = NULL, row_data = NULL) {
    if (is.null(ref_col) || is.null(row_data)) {
        stop("Internal standard normalization requires ref_col and row_data.")
    }
    if (!ref_col %in% colnames(row_data)) {
        stop("ref_col '", ref_col, "' not found in row_data.")
    }

    is_flag <- row_data[[ref_col]]
    if (is.logical(is_flag)) {
        is_rows <- which(is_flag)
    } else {
        is_rows <- which(tolower(as.character(is_flag)) %in% c("true", "yes", "is", "1"))
    }

    if (length(is_rows) == 0) {
        stop("No internal standard features found in row_data$", ref_col)
    }

    # IS intensity per sample (median if multiple)
    is_mat <- mat[is_rows, , drop = FALSE]
    is_intensity <- apply(is_mat, 2, stats::median, na.rm = TRUE)
    is_intensity[!is.finite(is_intensity) | is_intensity == 0] <- 1

    sweep(mat, 2, is_intensity, FUN = "/")
}


# ==== 2. TRANSFORMATION ======================================================

#' Apply transformation
#'
#' @param mat         Numeric matrix.
#' @param method      One of "none", "log2", "log10".
#' @param pseudocount Numeric value added before log (default 1).
#' @return Transformed matrix.
transform_metab <- function(mat, method = "none", pseudocount = 1) {
    method <- tolower(method)
    assert_one_of(method, "transform", c("none", "log2", "log10"))

    if (method == "none") return(mat)

    if (pseudocount > 0) {
        mat <- mat + pseudocount
    }

    # Warn about non-positive values
    if (any(mat <= 0, na.rm = TRUE)) {
        n_neg <- sum(mat <= 0, na.rm = TRUE)
        warning(sprintf(
            "transform_metab: %d non-positive values remain after pseudocount; will produce -Inf/NaN.",
            n_neg
        ))
    }

    switch(method,
        log2  = log2(mat),
        log10 = log10(mat)
    )
}


# ==== 3. SCALING ==============================================================

#' Apply feature scaling (row-wise)
#'
#' @param mat    Numeric matrix (features x samples).
#' @param method One of "none", "auto", "pareto", "range".
#' @return Scaled matrix.
scale_metab <- function(mat, method = "none") {
    method <- tolower(method)
    assert_one_of(method, "scaling", c("none", "center", "auto", "pareto", "range"))

    if (method == "none") return(mat)

    # Compute row statistics
    row_means <- rowMeans(mat, na.rm = TRUE)

    # Center first (all scaling methods center)
    mat_centered <- sweep(mat, 1, row_means, FUN = "-")

    if (method == "center") return(mat_centered)

    row_sds <- apply(mat, 1, stats::sd, na.rm = TRUE)

    switch(method,
        auto = {
            # Unit variance scaling (autoscaling)
            row_sds[row_sds == 0 | !is.finite(row_sds)] <- 1
            sweep(mat_centered, 1, row_sds, FUN = "/")
        },
        pareto = {
            # Divide by sqrt(sd) instead of sd
            pareto_factor <- sqrt(row_sds)
            pareto_factor[pareto_factor == 0 | !is.finite(pareto_factor)] <- 1
            sweep(mat_centered, 1, pareto_factor, FUN = "/")
        },
        range = {
            # Range scaling: divide by (max - min)
            row_ranges <- apply(mat, 1, function(x) {
                r <- diff(range(x, na.rm = TRUE))
                if (r == 0 || !is.finite(r)) 1 else r
            })
            sweep(mat_centered, 1, row_ranges, FUN = "/")
        }
    )
}


# ==== PIPELINE WRAPPER ========================================================
#' 
#' #' Apply full normalization pipeline: sample_norm → transform → scaling
#' #'
#' #' @param mat         Numeric matrix (features x samples).
#' #' @param norm_cfg    normalization config section.
#' #' @param row_data    data.frame (for IS normalization only).
#' #' @return list(expr_norm, applied) where applied records what was done.
#' apply_normalization_pipeline <- function(mat, norm_cfg, row_data = NULL) {
#'     sample_norm <- norm_cfg$sample_norm %||% "none"
#'     transform   <- norm_cfg$transform   %||% "none"
#'     scaling     <- norm_cfg$scaling     %||% "none"
#'     pseudocount <- norm_cfg$pseudocount %||% 1
#' 
#'     # IS-specific config
#'     is_ref_col <- norm_cfg$is_ref_col %||% NULL
#' 
#'     # Step 1: sample normalization
#'     mat <- normalize_samples(mat, method = sample_norm,
#'                              ref_col = is_ref_col, row_data = row_data)
#' 
#'     # Step 2: transformation
#'     mat <- transform_metab(mat, method = transform, pseudocount = pseudocount)
#' 
#'     # Step 3: scaling
#'     mat <- scale_metab(mat, method = scaling)
#' 
#'     list(
#'         expr_norm = mat,
#'         applied = list(
#'             sample_norm = sample_norm,
#'             transform   = transform,
#'             scaling     = scaling,
#'             pseudocount = if (transform != "none") pseudocount else NA
#'         )
#'     )
#' }
