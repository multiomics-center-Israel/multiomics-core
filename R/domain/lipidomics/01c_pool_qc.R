# R/domain/lipidomics/01c_pool_qc.R
#
# QC-pool Coefficient-of-Variation computation.
# Pool sub-matrix is captured in preprocess_lipidomics() (section 2.5) before
# the sample filter strips pools, so this stays a pure downstream function.
#
# Plotting helpers live alongside in the same file (added in step 3).


#' Compute per-feature CV across QC pool samples
#'
#' For each feature, CV (%) = 100 * SD / |mean|, computed over non-missing
#' pool values (NA and zero treated as missing). Features with fewer than
#' `min_valid` valid pool measurements get CV = NA.
#'
#' @param pool_matrix  Numeric matrix, features x pool samples. If NULL or
#'                     has fewer than `min_pools` columns, returns NULL.
#' @param min_pools    Minimum number of pool columns required to compute CV
#'                     (default 2 — CV is undefined with one observation).
#' @param min_valid    Minimum non-missing pool values per feature for that
#'                     feature's CV to be computed (default 2).
#' @return list(per_feature = data.frame, median_cv, mean_cv,
#'              n_features, n_pools) or NULL if pools are absent/insufficient.
compute_pool_cv <- function(pool_matrix, min_pools = 2L, min_valid = 2L) {
    if (is.null(pool_matrix) || ncol(pool_matrix) < min_pools) {
        n <- if (is.null(pool_matrix)) 0L else ncol(pool_matrix)
        message(sprintf(
            "compute_pool_cv: %d pool sample(s) available (need >= %d) -- skipping.",
            n, min_pools
        ))
        return(NULL)
    }

    # Treat NA and zero as missing -- pool intensities of 0 indicate non-detect.
    pm <- pool_matrix
    pm[!is.na(pm) & pm == 0] <- NA_real_

    n_valid <- rowSums(!is.na(pm))
    row_mean <- rowMeans(pm, na.rm = TRUE)
    row_sd   <- apply(pm, 1L, stats::sd, na.rm = TRUE)

    cv <- ifelse(
        n_valid >= min_valid & is.finite(row_mean) & abs(row_mean) > 0,
        100 * row_sd / abs(row_mean),
        NA_real_
    )

    feature_id <- rownames(pool_matrix)
    if (is.null(feature_id)) feature_id <- as.character(seq_len(nrow(pool_matrix)))

    per_feature <- data.frame(
        feature_id = feature_id,
        n_pools    = n_valid,
        mean       = row_mean,
        sd         = row_sd,
        cv         = cv,
        stringsAsFactors = FALSE
    )

    list(
        per_feature = per_feature,
        median_cv   = stats::median(cv, na.rm = TRUE),
        mean_cv     = mean(cv, na.rm = TRUE),
        n_features  = sum(!is.na(cv)),
        n_pools     = ncol(pool_matrix)
    )
}


#' Plot density of per-feature pool CVs
#'
#' Dashed line = median CV; dotted greys = `cv_thresholds` reference lines
#' (default 20 and 30%). X-axis clipped so reference lines stay visible without
#' being dominated by long-tail outliers.
#'
#' @param cv_result      List from `compute_pool_cv()`. Returns NULL if NULL.
#' @param cv_thresholds  Numeric vector of percent thresholds to mark.
#' @param title          Optional plot title.
#' @return ggplot object, or NULL when there are no valid CVs to plot.
plot_pool_cv_density <- function(cv_result,
                                  cv_thresholds = c(20, 30),
                                  title = NULL) {
    if (is.null(cv_result) || is.null(cv_result$per_feature)) return(NULL)

    df <- cv_result$per_feature
    df <- df[!is.na(df$cv), , drop = FALSE]
    if (nrow(df) == 0L) return(NULL)

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("plot_pool_cv_density: ggplot2 not available.")
        return(NULL)
    }

    med <- cv_result$median_cv
    upper <- max(c(cv_thresholds * 1.5,
                   stats::quantile(df$cv, 0.99, na.rm = TRUE) * 1.05),
                 na.rm = TRUE)

    title <- title %||% sprintf(
        "Pool CV distribution (n=%d features, %d pools)",
        nrow(df), cv_result$n_pools
    )
    subtitle <- sprintf("Median CV = %.1f%%", med)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$cv)) +
        ggplot2::geom_density(fill = "#4575b4", alpha = 0.35,
                              colour = "#4575b4") +
        ggplot2::geom_vline(xintercept = med, linetype = "dashed",
                            colour = "black") +
        ggplot2::labs(title = title, subtitle = subtitle,
                      x = "Coefficient of variation (%)",
                      y = "Density") +
        ggplot2::coord_cartesian(xlim = c(0, upper)) +
        ggplot2::theme_bw()

    for (thr in cv_thresholds) {
        p <- p + ggplot2::geom_vline(xintercept = thr,
                                      linetype = "dotted",
                                      colour = "grey45")
    }
    p
}
