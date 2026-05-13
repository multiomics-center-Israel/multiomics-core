#' Proteomics Advanced QC Functions
#'
#' UMAP dimensionality reduction, outlier detection, top DE heatmaps,
#' and filtering statistics for proteomics data.

# ==============================================================================
# UMAP
# ==============================================================================

#' Run UMAP on proteomics expression matrix
#'
#' @param expr_mat Numeric matrix (proteins x samples)
#' @param meta     Sample metadata data.frame
#' @param cfg      Proteomics config section (config$modes$proteomics)
#' @return list(coords, plot) or NULL
run_proteomics_umap <- function(expr_mat, meta, cfg) {
    if (!requireNamespace("umap", quietly = TRUE)) {
        message("Package 'umap' not available. Skipping UMAP.")
        return(NULL)
    }

    n_neighbors <- cfg$qc$umap_n_neighbors %||% 15
    min_dist    <- cfg$qc$umap_min_dist %||% 0.1

    # UMAP operates on samples (rows), so transpose
    mat_t <- t(expr_mat)

    # Remove features with zero variance
    keep <- apply(expr_mat, 1, function(x) var(x, na.rm = TRUE) > 0)
    mat_t <- mat_t[, keep, drop = FALSE]

    # Handle NAs by replacing with column median
    for (j in seq_len(ncol(mat_t))) {
        nas <- is.na(mat_t[, j])
        if (any(nas)) mat_t[nas, j] <- median(mat_t[, j], na.rm = TRUE)
    }

    # Adjust n_neighbors if fewer samples
    n_neighbors <- min(n_neighbors, nrow(mat_t) - 1)
    if (n_neighbors < 2) {
        message("Too few samples for UMAP.")
        return(NULL)
    }

    umap_cfg <- umap::umap.defaults
    umap_cfg$n_neighbors <- n_neighbors
    umap_cfg$min_dist    <- min_dist

    umap_res <- tryCatch(
        umap::umap(mat_t, config = umap_cfg),
        error = function(e) { message("UMAP failed: ", e$message); NULL }
    )
    if (is.null(umap_res)) return(NULL)

    coords <- as.data.frame(umap_res$layout)
    colnames(coords) <- c("UMAP1", "UMAP2")

    color_col  <- cfg$effects$color %||% "Condition"
    sample_col <- cfg$effects$samples %||% "SampleID"

    coords$Sample <- if (sample_col %in% colnames(meta)) meta[[sample_col]] else rownames(mat_t)
    if (color_col %in% colnames(meta)) {
        coords$Group <- meta[[color_col]]
    } else {
        coords$Group <- "all"
    }

    p <- ggplot2::ggplot(coords, ggplot2::aes(x = UMAP1, y = UMAP2, color = Group)) +
        ggplot2::geom_point(size = 3) +
        ggplot2::labs(title = "UMAP — Protein Expression",
                      x = "UMAP1", y = "UMAP2", color = color_col) +
        ggplot2::theme_minimal()

    list(coords = coords, plot = p)
}

# ==============================================================================
# OUTLIER DETECTION
# ==============================================================================

#' Detect outlier samples via Mahalanobis distance + low correlation
#'
#' @param expr_mat   Numeric matrix (proteins x samples)
#' @param pca_result prcomp result (or NULL; extracted from QC pre objects)
#' @param cor_mat    Sample correlation matrix
#' @param cfg        Proteomics config section
#' @return list(flagged_samples_df, n_outliers)
detect_proteomics_outliers <- function(expr_mat, pca_result, cor_mat, cfg) {
    sd_thresh  <- cfg$qc$outlier_sd_threshold %||% 3
    min_corr   <- cfg$qc$min_sample_correlation %||% 0.8
    sample_col <- cfg$effects$samples %||% "SampleID"

    sample_names <- colnames(expr_mat)
    n <- length(sample_names)

    pca_outlier   <- rep(FALSE, n)
    corr_outlier  <- rep(FALSE, n)
    reasons       <- rep("", n)

    # --- Mahalanobis on PCA scores (PC1:PC3) ---
    if (!is.null(pca_result) && inherits(pca_result, "prcomp")) {
        n_pcs <- min(3, ncol(pca_result$x))
        scores <- pca_result$x[, seq_len(n_pcs), drop = FALSE]
        if (nrow(scores) == n && n_pcs >= 2) {
            tryCatch({
                center <- colMeans(scores)
                cov_mat <- cov(scores)
                md <- stats::mahalanobis(scores, center, cov_mat)
                # Flag samples beyond sd_thresh standard deviations
                md_threshold <- median(md) + sd_thresh * mad(md)
                pca_outlier <- md > md_threshold
                reasons[pca_outlier] <- paste0(reasons[pca_outlier],
                    "Mahalanobis distance above threshold; ")
            }, error = function(e) {
                message("Mahalanobis outlier detection failed: ", e$message)
            })
        }
    }

    # --- Low correlation flagging ---
    if (!is.null(cor_mat) && nrow(cor_mat) == n) {
        # Median correlation of each sample with all others
        diag(cor_mat) <- NA
        med_corr <- apply(cor_mat, 1, median, na.rm = TRUE)
        corr_outlier <- med_corr < min_corr
        reasons[corr_outlier] <- paste0(reasons[corr_outlier],
            sprintf("Median correlation %.3f < %.2f; ", med_corr[corr_outlier], min_corr))
    }

    is_outlier <- pca_outlier | corr_outlier

    flagged_df <- data.frame(
        Sample         = sample_names,
        pca_outlier    = pca_outlier,
        corr_outlier   = corr_outlier,
        is_outlier     = is_outlier,
        reason         = trimws(reasons),
        stringsAsFactors = FALSE
    )

    # Only return flagged samples
    flagged_df <- flagged_df[flagged_df$is_outlier, , drop = FALSE]

    list(
        flagged_samples_df = flagged_df,
        n_outliers         = sum(is_outlier)
    )
}

# ==============================================================================
# TOP DE HEATMAP
# ==============================================================================

#' Plot heatmap of top DE proteins
#'
#' @param de_tbl        DE table with FeatureID, adj.P.Val, logFC columns
#' @param expr_mat      Numeric matrix (proteins x samples)
#' @param meta          Sample metadata
#' @param cfg           Proteomics config section
#' @param n_top         Number of top DE proteins
#' @param contrast_name Name of the contrast (for title)
#' @return pheatmap object or NULL
plot_top_de_heatmap <- function(de_tbl, expr_mat, meta, cfg, n_top = 50,
                                 contrast_name = "contrast") {
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        message("Package 'pheatmap' not available. Skipping heatmap.")
        return(NULL)
    }

    if (is.null(de_tbl) || nrow(de_tbl) == 0) return(NULL)

    # Sort by adjusted p-value, take top N
    de_tbl <- de_tbl[order(de_tbl$adj.P.Val, na.last = TRUE), ]
    top_ids <- head(de_tbl$FeatureID[!is.na(de_tbl$adj.P.Val)], n_top)
    top_ids <- top_ids[top_ids %in% rownames(expr_mat)]

    if (length(top_ids) < 2) return(NULL)

    mat_sub <- expr_mat[top_ids, , drop = FALSE]

    # Z-score rows
    mat_z <- t(scale(t(mat_sub)))
    mat_z[is.na(mat_z)] <- 0

    # Annotation
    color_col <- cfg$effects$color %||% "Condition"
    annot_col <- NULL
    if (color_col %in% colnames(meta)) {
        annot_col <- data.frame(
            Condition = meta[[color_col]],
            row.names = colnames(mat_z)
        )
    }

    tryCatch({
        pheatmap::pheatmap(
            mat_z,
            annotation_col   = annot_col,
            cluster_rows     = TRUE,
            cluster_cols     = TRUE,
            show_rownames    = n_top <= 50,
            show_colnames    = TRUE,
            fontsize_row     = max(4, 10 - n_top / 10),
            main             = sprintf("Top %d DE Proteins — %s", length(top_ids), contrast_name),
            color            = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100),
            silent           = TRUE,
          border_color = NA
        )
    }, error = function(e) {
        message("Heatmap failed for ", contrast_name, ": ", e$message)
        NULL
    })
}

# ==============================================================================
# FILTERING STATISTICS
# ==============================================================================

#' Build filtering statistics table
#'
#' @param pre    Preprocessed proteomics object
#' @param config Full config list
#' @return data.frame with filtering metrics
build_filtering_stats <- function(pre, config) {
    n_raw <- if (!is.null(pre$expr_raw)) nrow(pre$expr_raw) else NA
    n_filtered <- nrow(pre$expr_imp_single)
    n_removed <- if (!is.na(n_raw)) n_raw - n_filtered else NA
    pct_retained <- if (!is.na(n_raw) && n_raw > 0) round(100 * n_filtered / n_raw, 1) else NA

    # Missingness
    pct_missing_before <- NA
    if (!is.null(pre$expr_raw)) {
        total_vals <- prod(dim(pre$expr_raw))
        if (total_vals > 0) {
            pct_missing_before <- round(100 * sum(is.na(pre$expr_raw)) / total_vals, 1)
        }
    }

    n_imputed <- NA
    if (!is.null(pre$imputation_qc) && !is.null(pre$imputation_qc$imputed_flag)) {
        n_imputed <- sum(pre$imputation_qc$imputed_flag, na.rm = TRUE)
    }

    data.frame(
        Metric = c("Proteins (raw input)", "Proteins (after filtering)",
                    "Proteins removed", "Retained (%)",
                    "Missing values before filtering (%)",
                    "Imputed values (single imp.)"),
        Value = c(
            format_metric(n_raw), format_metric(n_filtered),
            format_metric(n_removed), format_metric(pct_retained),
            format_metric(pct_missing_before), format_metric(n_imputed)
        ),
        stringsAsFactors = FALSE
    )
}

#' Format a metric for display (handles NA)
#' @noRd
format_metric <- function(x) {
    if (is.na(x)) return("N/A")
    if (is.numeric(x)) return(format(x, big.mark = ","))
    as.character(x)
}
