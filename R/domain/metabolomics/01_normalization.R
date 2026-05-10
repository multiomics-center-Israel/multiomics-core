# R/domain/metabolomics/01_normalization.R
#
# Three-step normalization for metabolomics data.
# Three independent steps (applied in order):
#   1. Sample normalization  (sum, median, PQN, internal standard)
#   2. Transformation        (log2, log10, glog10)
#   3. Scaling               (auto, pareto, range)
#
# Each function operates on a numeric matrix [features x samples].
# Reuses: assert_numeric_matrix


# ==== 1. SAMPLE NORMALIZATION ================================================

#' Apply sample normalization
#'
#' @param mat     Numeric matrix (features x samples).
#' @param method  One of "none", "sum", "median", "pqn", "eigenms", "is".
#' @param ref_col Character: column in row_data that flags the internal standard
#'                (only used when method = "is").
#' @param row_data data.frame with at least a column matching ref_col (for IS).
#' @param groups  Character/factor vector of group labels per sample
#'                (required for eigenms; ignored by other methods).
#' @return Normalized numeric matrix (same dimensions).
normalize_samples <- function(mat, method = "none", ref_col = NULL, row_data = NULL,
                              groups = NULL, ref_samples = NULL,
                              meta = NULL, bio_factor_col = NULL) {
    method <- tolower(method)
    assert_one_of(method, "sample_norm",
                  c("none", "sum", "median", "pqn", "eigenms", "eigenms_forced",
                    "is", "bio_factor"))

    switch(method,
        none           = mat,
        sum            = norm_total_sum(mat),
        median         = norm_median(mat),
        pqn            = norm_pqn(mat, ref_samples = ref_samples),
        eigenms        = norm_eigenms(mat, groups),
        eigenms_forced = norm_eigenms_forced(mat, groups),
        is             = norm_internal_standard(mat, ref_col, row_data),
        bio_factor     = norm_biological_factor(mat, meta, bio_factor_col),
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
#' 1. Perform median normalization first.
#' 2. Compute reference spectrum:
#'    - If \code{ref_samples} is provided, use the median of those samples
#'      (e.g. QC pools) as the reference spectrum.
#'    - Otherwise, use the median across all samples.
#' 3. For each sample, compute quotients (sample / reference) and take the median.
#' 4. Divide each sample by its median quotient.
#'
#' Features where the reference is 0 or NA are excluded from the quotient
#' calculation (they cannot inform sample-level scaling).
#'
#' @param mat          Numeric matrix (features x samples).
#' @param ref_samples  Optional character vector of sample (column) names to use
#'                     as the reference spectrum. If a single sample is given, that
#'                     sample's values are used directly. If multiple, their median
#'                     per feature is used. If NULL (default), the median across
#'                     all samples is used.
#' @return Normalized numeric matrix (same dimensions).
#'
#' @references Dieterle F et al. (2006) Anal Chem 78:4281-4290.
norm_pqn <- function(mat, ref_samples = NULL) {
    # Step 1: initial median normalization
    mat_mn <- norm_median(mat)

    # Step 2: reference spectrum
    if (!is.null(ref_samples)) {
        ref_cols <- intersect(ref_samples, colnames(mat_mn))
        if (length(ref_cols) == 0) {
            warning("norm_pqn: none of ref_samples found in data columns; ",
                    "falling back to median of all samples.")
            ref <- apply(mat_mn, 1, stats::median, na.rm = TRUE)
        } else if (length(ref_cols) == 1) {
            message(sprintf("norm_pqn: using single reference sample '%s'", ref_cols))
            ref <- mat_mn[, ref_cols]
        } else {
            message(sprintf("norm_pqn: using median of %d reference samples (%s)",
                            length(ref_cols), paste(ref_cols, collapse = ", ")))
            ref <- apply(mat_mn[, ref_cols, drop = FALSE], 1, stats::median, na.rm = TRUE)
        }
    } else {
        ref <- apply(mat_mn, 1, stats::median, na.rm = TRUE)
    }

    # Exclude features with zero/NA reference — they produce Inf quotients
    ref_ok <- is.finite(ref) & ref != 0
    if (!any(ref_ok)) {
        warning("norm_pqn: all reference values are 0 or NA; returning median-normalized data.")
        return(mat_mn)
    }

    # Step 3 & 4: per-sample quotient normalization (using valid features only)
    quotients <- sweep(mat_mn[ref_ok, , drop = FALSE], 1,
                       ref[ref_ok], FUN = "/")
    med_quotients <- apply(quotients, 2, stats::median, na.rm = TRUE)
    med_quotients[!is.finite(med_quotients) | med_quotients == 0] <- 1

    sweep(mat_mn, 2, med_quotients, FUN = "/")
}


#' EigenMS normalization (Eigenvector-based normalization)
#'
#' Uses the ProteoMM Bioconductor package (Karpievitch et al., 2009, 2014) to
#' perform EigenMS normalization via SVD on residuals from a group-aware linear
#' model, identifying and removing systematic technical bias while preserving
#' biological signal.
#'
#' References:
#'   Karpievitch YV et al. (2009) Bioinformatics — original EigenMS method.
#'   Karpievitch YV et al. (2014) PLoS ONE — EigenMS for metabolomics.
#'
#' @param mat    Numeric matrix (features x samples), log-scale or linear.
#' @param groups Character/factor vector of group labels per sample.
#' @return Normalized numeric matrix (same dimensions).
norm_eigenms <- function(mat, groups = NULL) {
    if (is.null(groups)) {
        warning("norm_eigenms: groups not provided — falling back to PQN normalization.")
        return(norm_pqn(mat))
    }

    if (!requireNamespace("ProteoMM", quietly = TRUE)) {
        stop("norm_eigenms requires the ProteoMM Bioconductor package.\n",
             "Install with: BiocManager::install('ProteoMM')")
    }

    groups <- as.factor(groups)
    if (length(groups) != ncol(mat)) {
        stop("norm_eigenms: length(groups) must equal ncol(mat)")
    }

    # ProteoMM requires rownames
    if (is.null(rownames(mat))) rownames(mat) <- paste0("feat_", seq_len(nrow(mat)))
    if (is.null(colnames(mat))) colnames(mat) <- paste0("S", seq_len(ncol(mat)))

    # prot.info: data.frame with at least one column of feature IDs
    prot_info <- data.frame(id = rownames(mat), stringsAsFactors = FALSE)

    # Suppress ProteoMM's built-in plotting (it calls par/plot internally)
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    # Step 1: Identify eigentrends (SVD + parallel analysis)
    set.seed(12345)  # reproducibility as recommended by ProteoMM
    rv <- suppressWarnings(suppressMessages(
        ProteoMM::eig_norm1(m = mat, treatment = groups, prot.info = prot_info)
    ))

    n_eigentrends <- rv$h.c
    message(sprintf("norm_eigenms (ProteoMM): found %d significant eigentrend(s)",
                    n_eigentrends))
    # Store eigentrend count as attribute for downstream reporting
    attr_info <- list(n_eigentrends = n_eigentrends, method = "ProteoMM_bootstrap")

    if (n_eigentrends == 0) {
        message("norm_eigenms: no significant systematic bias detected — returning original data.")
        attr(mat, "eigenms_info") <- attr_info
        return(mat)
    }

    # Step 2: Normalize by removing identified eigentrends
    res <- suppressWarnings(suppressMessages(
        ProteoMM::eig_norm2(rv)
    ))

    # res$norm_m is the normalized matrix (features x samples) with row/colnames
    mat_norm <- res$norm_m

    # Ensure output has same dimensions and names as input
    # ProteoMM may drop features that are entirely NA in any group
    if (nrow(mat_norm) == nrow(mat)) {
        rownames(mat_norm) <- rownames(mat)
        colnames(mat_norm) <- colnames(mat)
        attr(mat_norm, "eigenms_info") <- attr_info
        return(mat_norm)
    }

    # If ProteoMM dropped some features, map back to original dimensions
    mat_out <- mat
    shared <- intersect(rownames(mat_norm), rownames(mat))
    if (length(shared) > 0) {
        mat_out[shared, ] <- mat_norm[shared, ]
    }
    # Features not returned by ProteoMM keep their original values
    n_kept <- length(shared)
    if (n_kept < nrow(mat)) {
        message(sprintf("norm_eigenms: ProteoMM normalized %d/%d features; %d kept original values.",
                        n_kept, nrow(mat), nrow(mat) - n_kept))
    }
    attr(mat_out, "eigenms_info") <- attr_info
    mat_out
}


#' EigenMS normalization — forced eigentrend removal (NOREVA-style)
#'
#' Unlike the ProteoMM version which uses bootstrap parallel analysis to
#' determine the number of significant eigentrends, this implementation forces
#' removal of a fixed number of eigentrends (default: 20% of samples, minimum 1).
#' This is the approach used in NOREVA (Liang et al., 2020) and by the Technion
#' metabolomics unit (Ifat Abramovich).
#'
#' WARNING: This method does not test whether the eigentrends are statistically
#' significant. It may remove real biological variation in addition to technical
#' bias. Use with caution and compare results to the ProteoMM version.
#'
#' References:
#'   Karpievitch YV et al. (2014) PLoS ONE — EigenMS for metabolomics.
#'   Liang et al. (2020) Nature Protocols — NOREVA.
#'
#' @param mat    Numeric matrix (features x samples), log-scale.
#' @param groups Character/factor vector of group labels per sample.
#' @param pct_eigentrends Fraction of samples to use as eigentrend count (default 0.2).
#' @return Normalized numeric matrix (same dimensions).
norm_eigenms_forced <- function(mat, groups = NULL, pct_eigentrends = 0.2) {
    if (is.null(groups)) {
        warning("norm_eigenms_forced: groups not provided — falling back to PQN.")
        return(norm_pqn(mat))
    }

    groups <- as.factor(groups)
    if (length(groups) != ncol(mat)) {
        stop("norm_eigenms_forced: length(groups) must equal ncol(mat)")
    }

    n_features <- nrow(mat)
    n_samples  <- ncol(mat)

    # Transpose to samples x metabolites for lm fitting
    m <- t(mat)

    residuals_m <- matrix(NA_real_, n_samples, n_features)
    fitted_m    <- matrix(NA_real_, n_samples, n_features)

    for (j in seq_len(n_features)) {
        y <- m[, j]
        if (any(is.na(y))) next
        fit <- stats::lm(y ~ groups)
        residuals_m[, j] <- stats::resid(fit)
        fitted_m[, j]    <- stats::fitted(fit)
    }

    ok_cols <- which(colSums(is.na(residuals_m)) == 0)
    if (length(ok_cols) < 2) {
        warning("norm_eigenms_forced: too few complete features; returning original.")
        return(mat)
    }

    R <- residuals_m[, ok_cols]
    svd_res <- svd(R)

    # Force removal of n_bias eigentrends (20% of samples, min 1)
    n_bias <- max(1L, round(pct_eigentrends * n_samples))
    n_bias <- min(n_bias, length(svd_res$d) - 1L)

    message(sprintf("norm_eigenms_forced: forcing removal of %d eigentrend(s) from %d samples",
                    n_bias, n_samples))

    U_bias <- svd_res$u[, seq_len(n_bias), drop = FALSE]
    R_corrected <- R - U_bias %*% (t(U_bias) %*% R)

    result <- m
    result[, ok_cols] <- fitted_m[, ok_cols] + R_corrected

    # Transpose back to features x samples
    out <- t(result)
    attr(out, "eigenms_info") <- list(n_eigentrends = n_bias, method = "forced_SVD",
                                      pct_eigentrends = pct_eigentrends)
    out
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


#' Biological-factor normalization
#'
#' Divide each sample column by a per-sample numeric value (e.g., mg of total
#' protein) read from a metadata column. Used when intensities should be
#' normalized to a measured biological covariate rather than to a
#' sample-level statistic.
#'
#' @param mat        Numeric matrix (features x samples).
#' @param meta       data.frame with per-sample rows; must contain
#'                   `factor_col` and a sample-id column matching colnames(mat)
#'                   (or be in the same row order as colnames(mat)).
#' @param factor_col Name of the column in meta with per-sample numeric values.
#' @return Normalized matrix.
norm_biological_factor <- function(mat, meta, factor_col) {
    if (is.null(meta) || is.null(factor_col) || !nzchar(factor_col)) {
        stop("biological_factor normalization requires both meta and ",
             "biological_factor_col to be set in config.")
    }
    if (!factor_col %in% colnames(meta)) {
        stop("biological_factor: column '", factor_col,
             "' not found in metadata.")
    }

    # Match meta rows to sample column order via a sample-id column when present
    sample_cols <- intersect(c("sample_id", "SampleID", "Sample.ID", "sample"),
                             colnames(meta))
    if (length(sample_cols) > 0) {
        idx <- match(colnames(mat), meta[[sample_cols[1]]])
        if (any(is.na(idx))) {
            missing <- colnames(mat)[is.na(idx)]
            stop("biological_factor: samples missing from meta: ",
                 paste(utils::head(missing, 5), collapse = ", "))
        }
        factors <- as.numeric(meta[[factor_col]][idx])
    } else {
        if (nrow(meta) != ncol(mat)) {
            stop("biological_factor: cannot align meta to samples; ",
                 "no sample-id column found and row count mismatch.")
        }
        factors <- as.numeric(meta[[factor_col]])
    }

    if (any(is.na(factors)) || any(factors <= 0)) {
        bad <- colnames(mat)[is.na(factors) | factors <= 0]
        stop("biological_factor: invalid (NA or non-positive) values for ",
             "samples: ", paste(utils::head(bad, 5), collapse = ", "))
    }

    sweep(mat, 2, factors, FUN = "/")
}


# ==== 2. TRANSFORMATION ======================================================

#' Apply transformation
#'
#' @param mat         Numeric matrix.
#' @param method      One of "none", "log2", "log10", "glog10".
#' @param pseudocount Numeric value added before log (default 1).
#'                    Ignored when method = "glog10" (uses adaptive min.val).
#' @return Transformed matrix.
transform_metab <- function(mat, method = "none", pseudocount = 1) {
    method <- tolower(method)
    assert_one_of(method, "transform", c("none", "log2", "log10", "glog10"))

    if (method == "none") return(mat)

    # Generalized log10 (MetaboAnalyst convention):
    #   glog10(x) = log10( (x + sqrt(x^2 + min.val^2)) / 2 )
    # where min.val = 1/10 of the smallest positive value.
    # Smoothly handles zeros and near-zero values without a pseudocount.
    if (method == "glog10") {
        pos_vals <- mat[mat > 0 & is.finite(mat)]
        min_val <- if (length(pos_vals) > 0) min(pos_vals) / 10 else 1
        return(log10((mat + sqrt(mat^2 + min_val^2)) / 2))
    }

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

#' Apply full normalization pipeline: sample_norm → transform → scaling
#'
#' @param mat         Numeric matrix (features x samples).
#' @param norm_cfg    normalization config section.
#' @param row_data    data.frame (for IS normalization only).
#' @param groups      Character/factor vector of group labels per sample
#'                    (required for eigenms; ignored by other methods).
#' @return list(expr_norm, applied) where applied records what was done.
apply_normalization_pipeline <- function(mat, norm_cfg, row_data = NULL, groups = NULL,
                                         meta = NULL) {
    sample_norm <- norm_cfg$sample_norm %||% "none"
    transform   <- norm_cfg$transform   %||% "none"
    scaling     <- norm_cfg$scaling     %||% "none"
    pseudocount <- norm_cfg$pseudocount %||% 1

    # IS-specific config
    is_ref_col <- norm_cfg$is_ref_col %||% NULL

    # bio_factor-specific config (per-sample numeric column in meta)
    bio_factor_col <- norm_cfg$biological_factor_col %||% NULL

    # Step 1: sample normalization
    mat <- normalize_samples(mat, method = sample_norm,
                             ref_col = is_ref_col, row_data = row_data,
                             groups = groups,
                             meta = meta, bio_factor_col = bio_factor_col)

    # Step 2: transformation
    mat <- transform_metab(mat, method = transform, pseudocount = pseudocount)

    # Step 3: scaling
    mat <- scale_metab(mat, method = scaling)

    list(
        expr_norm = mat,
        applied = list(
            sample_norm = sample_norm,
            transform   = transform,
            scaling     = scaling,
            pseudocount = if (transform != "none") pseudocount else NA
        )
    )
}
