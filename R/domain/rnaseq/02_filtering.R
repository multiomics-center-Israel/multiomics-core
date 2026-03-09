#' Determine which features pass the filter based on min_count per condition
pass_filter <- function(expr_mat, group, min_per_group) {
    expr_mat <- as.matrix(expr_mat)
    group <- as.character(group)
    groups <- unique(group)

    if (length(min_per_group) == 1 && is.null(names(min_per_group))) {
        min_per_group <- setNames(rep(min_per_group, length(groups)), groups)
    }

    passes_per_group <- sapply(groups, function(g) {
        cols <- which(group == g)
        if (length(cols) == 0) {
            return(rep(FALSE, nrow(expr_mat)))
        }
        sums <- rowSums(!is.na(expr_mat[, cols, drop = FALSE]))
        sums >= min_per_group[[g]]
    })

    apply(passes_per_group, 1, any)
}

# ==============================================================================
# 1. Automatic Threshold Detection (KDE-based)
# ==============================================================================
find_optimal_threshold_safe <- function(cpm_mat, min_limit = 0.5, max_limit = 2.0, fallback = 1.0) {
    # Log2 transformation for density estimation
    log_vals <- c(as.matrix(log2(cpm_mat + 0.1)))
    log_vals <- log_vals[!is.na(log_vals)]
    d <- density(log_vals)

    # Identify peaks and valleys (local maxima and minima)
    peaks <- which(diff(sign(diff(d$y))) == -2) + 1
    valleys <- which(diff(sign(diff(d$y))) == 2) + 1

    calculated_threshold <- NA

    # Look for a valley between the noise peak and the signal peak
    if (length(peaks) >= 2 && length(valleys) >= 1) {
        # Assume the first peak is noise
        noise_peak_idx <- peaks[which.min(d$x[peaks])]

        # Find valleys to the right of the noise peak
        relevant_valleys <- valleys[d$x[valleys] > d$x[noise_peak_idx]]

        if (length(relevant_valleys) > 0) {
            optimal_log_val <- d$x[relevant_valleys[1]]
            calculated_threshold <- 2^optimal_log_val
        }
    }

    # --- Safety Checks (Guardrails) ---

    if (is.na(calculated_threshold)) {
        warning("Auto-threshold failed to find a valley. Using fallback.")
        return(fallback)
    }

    if (calculated_threshold < min_limit) {
        warning(sprintf(
            "Auto-threshold (%.2f) too low. Clamping to min_limit (%.2f).",
            calculated_threshold, min_limit
        ))
        return(min_limit)
    }

    if (calculated_threshold > max_limit) {
        warning(sprintf(
            "Auto-threshold (%.2f) too high. Clamping to max_limit (%.2f).",
            calculated_threshold, max_limit
        ))
        return(max_limit)
    }

    return(calculated_threshold)
}

# ==============================================================================
# 2. Optimized Filtering Engine (Vectorized)
# ==============================================================================
filter_features_optimized <- function(norm_mat, meta, sample_col, group_col, threshold) {
    # Validity checks
    stopifnot(sample_col %in% colnames(meta))
    stopifnot(group_col %in% colnames(meta))

    norm_mat <- as.matrix(norm_mat)

    # Ensure samples match metadata
    if (!all(colnames(norm_mat) %in% meta[[sample_col]])) {
        stop("Samples missing in metadata")
    }

    # Align metadata
    meta2 <- meta[match(colnames(norm_mat), meta[[sample_col]]), , drop = FALSE]
    group_list <- split(meta2[[sample_col]], meta2[[group_col]])

    # Initialize keep vector (FALSE for all genes)
    keep_vec <- rep(FALSE, nrow(norm_mat))

    # Iterate over groups
    for (grp in names(group_list)) {
        samp <- group_list[[grp]]
        sub_mat <- norm_mat[, samp, drop = FALSE]

        # Calculate replicates present per gene
        n_reps <- rowSums(!is.na(sub_mat))

        # Determine minimum required replicates (at least half, rounded up)
        min_pass <- ceiling(n_reps / 2)

        # Count how many samples pass the threshold
        n_above_threshold <- rowSums(sub_mat >= threshold, na.rm = TRUE)

        # Update keep vector using OR logic
        keep_vec <- keep_vec | (n_above_threshold >= min_pass)
    }

    list(filtered = norm_mat[keep_vec, , drop = FALSE], keep_vec = keep_vec)
}

# ==============================================================================
# 3. Main Pipeline Wrapper
# ==============================================================================
run_auto_filter_pipeline <- function(cpm_mat, meta, sample_col, group_col,
                                     output_plot = "threshold_qc.png") {
    message("--- Starting Auto-Filtering Pipeline ---")

    # Step 1: Calculate Optimal Threshold
    auto_thresh <- find_optimal_threshold_safe(cpm_mat)
    message(sprintf("Optimal Threshold: %.3f CPM", auto_thresh))

    # Step 2: Generate QC Plot
    if (!is.null(output_plot)) {
        # Ensure directory exists for plot
        ensure_dir(dirname(output_plot))

        png(output_plot, width = 800, height = 600)

        log_vals <- c(as.matrix(log2(cpm_mat + 0.1)))
        log_vals <- log_vals[!is.na(log_vals)]
        d <- density(log_vals)

        plot(d, main = "Automated Threshold Selection", xlab = "log2(CPM)", lwd = 2)
        polygon(d, col = "gray90", border = "gray")
        abline(v = log2(auto_thresh), col = "red", lwd = 2)
        text(
            x = log2(auto_thresh), y = max(d$y) * 0.9,
            labels = paste("Cutoff:", round(auto_thresh, 2)),
            pos = 4, col = "red"
        )

        dev.off()
        message(paste("QC Plot saved to:", output_plot))
    }

    # Step 3: Execute Filtering
    filter_res <- filter_features_optimized(
        norm_mat = cpm_mat,
        meta = meta,
        sample_col = sample_col,
        group_col = group_col,
        threshold = auto_thresh
    )

    n_original <- nrow(cpm_mat)
    n_kept <- sum(filter_res$keep_vec)
    message(sprintf(
        "Done: Kept %d / %d genes (%.1f%%)",
        n_kept, n_original, 100 * n_kept / n_original
    ))

    return(list(
        filtered_data = filter_res$filtered,
        keep_vec = filter_res$keep_vec,
        used_threshold = auto_thresh
    ))
}

