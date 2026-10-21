#' Determine which features pass the filter based on min_count per condition
#'
#' @param expr_mat Features x samples matrix; a cell counts as observed when it
#'   is not `NA`.
#' @param group Group label per column of `expr_mat`, in column order.
#' @param min_per_group Minimum observations a feature needs within a group,
#'   either one number for every group or a vector named by group.
#' @param min_groups How many groups a feature must reach that minimum in.
#' @return One logical per feature, in row order, named by `rownames(expr_mat)`.
pass_filter <- function(expr_mat, group, min_per_group, min_groups = 1) {
    expr_mat <- as.matrix(expr_mat)
    group <- as.character(group)
    groups <- unique(group)

    if (length(min_per_group) == 1 && is.null(names(min_per_group))) {
        min_per_group <- setNames(rep(min_per_group, length(groups)), groups)
    }

    # sapply() simplifies on the common length of the per-group results, so with
    # a single feature it hands back a plain vector and the rowSums() below
    # aborts. State the features x groups shape instead of inferring it.
    passes_per_group <- matrix(
        sapply(groups, function(g) {
            cols <- which(group == g)
            if (length(cols) == 0) {
                return(rep(FALSE, nrow(expr_mat)))
            }
            sums <- rowSums(!is.na(expr_mat[, cols, drop = FALSE]))
            sums >= min_per_group[[g]]
        }),
        nrow = nrow(expr_mat), ncol = length(groups),
        dimnames = list(rownames(expr_mat), groups)
    )

    # Count how many groups each feature passes in, compare to min_groups
    n_groups_passed <- rowSums(passes_per_group)
    n_groups_passed >= min_groups
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
                                     output_plot = "threshold_qc.png",
                                     min_limit = 0.5,
                                     max_limit = 2.0,
                                     fallback  = 1.0) {
    message("--- Starting Auto-Filtering Pipeline ---")

    # Step 1: Calculate Optimal Threshold
    auto_thresh <- find_optimal_threshold_safe(cpm_mat,
                                               min_limit = min_limit,
                                               max_limit = max_limit,
                                               fallback  = fallback)
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

# ==============================================================================
# 4. Reader-facing description of the filter that ran
# ==============================================================================

#' Describe the RNA-seq gene filter that actually ran
#'
#' The one reading of what \code{preprocess_rna()} recorded about gene
#' filtering, shared by the HTML report and the pipeline summary so that both
#' describe the same run. The mode, the matrix filtered on, the threshold and
#' the grouping column come from \code{pre$info}, not from the config: the
#' runtime treats an unrecognised mode as adaptive and skips filtering for
#' preprocessed input whatever the mode says. The config is read only for the
#' adaptive search range and fallback, with the defaults
#' \code{preprocess_rna()} applies.
#'
#' @param info The \code{info} list of the object returned by
#'   \code{preprocess_rna()}, or NULL when it is not available (e.g. a report
#'   re-knit by hand); the description then says so instead of guessing.
#' @param filtering_cfg The \code{modes$rna$filtering} config list (may be NULL).
#' @return A list with \code{mode} (\code{"adaptive"}, \code{"fixed"},
#'   \code{"deseq2_only"}, \code{"none"}, \code{"preprocessed"} or
#'   \code{"unknown"}), \code{applied} (logical: a gene filter removed genes by
#'   rule), \code{matrix} (matrix the threshold was applied to, e.g. \code{"CPM"},
#'   or NA when no threshold was applied), \code{group_col}
#'   (grouping column of the group rule, or NA), \code{threshold}
#'   (numeric or NA), \code{threshold_label}, \code{rule} (or NULL),
#'   \code{short} (one-line label), \code{methods} (reader-facing paragraph)
#'   and \code{de_note} (DESeq2 independent filtering, or NULL).
describe_rna_filtering <- function(info, filtering_cfg = NULL) {
    filter_norm <- info$filter_norm %||% NA_character_
    mode <- info$filter_mode %||% (
        if (filter_norm %in% c("preprocessed", "none")) filter_norm else "unknown"
    )

    unit <- if (grepl("^TPM", filter_norm)) "TPM" else if (identical(filter_norm, "CPM")) "CPM" else NA_character_
    unit_label <- if (identical(filter_norm, "TPM (tximport)")) "TPM (tximport abundance)" else unit
    num <- function(x) format(signif(as.numeric(x), 3))

    thr <- suppressWarnings(as.numeric(info$filter_threshold %||% NA))[1]
    has_thr <- mode %in% c("adaptive", "fixed") && !is.na(thr)
    threshold_label <- if (has_thr) {
        sprintf("%s %s", num(thr), unit)
    } else if (mode == "unknown") {
        "N/A (not recorded for this report)"
    } else {
        "N/A (no threshold-based filtering)"
    }

    group <- info$filter_group_col %||% NA_character_
    group_txt <- if (is.na(group)) "the grouping column" else sprintf("'%s'", group)
    rule <- if (mode %in% c("adaptive", "fixed")) {
        sprintf(paste0(
            "A gene was kept if, in at least one group of %s, at least half ",
            "(rounded up) of the samples with a non-missing value had %s \u2265 the threshold."
        ), group_txt, unit)
    } else {
        NULL
    }

    methods <- switch(mode,
        adaptive = paste(
            sprintf(paste0(
                "Lowly expressed genes were filtered on %s values. The threshold was selected ",
                "automatically from a kernel density estimate of log2(%s + 0.1) values pooled over ",
                "all genes and samples: with v the position of the first local minimum to the ",
                "right of the lowest-expression peak, the threshold is 2^v. A threshold outside ",
                "%s\u2013%s %s is set to the nearest limit; when the density has fewer than two ",
                "peaks or no such minimum, the fallback of %s %s is used as configured. The ",
                "threshold applied in this run was %s."
            ), unit_label, unit,
            num(filtering_cfg$auto_filter$min %||% 0.5), num(filtering_cfg$auto_filter$max %||% 2.0),
            unit, num(filtering_cfg$auto_filter$fallback %||% 1.0), unit, threshold_label),
            rule),
        fixed = paste(
            sprintf("Lowly expressed genes were filtered on %s values with a fixed threshold of %s.",
                    unit_label, threshold_label),
            rule),
        deseq2_only = paste(
            "Only genes with a total count of zero across all samples were removed;",
            "no expression threshold was applied."),
        none = paste0(
            "No gene filtering was applied; all genes were kept.",
            if (identical(info$expr_work_method, "VST") &&
                isTRUE(info$n_genes_filt < info$n_genes_raw)) {
                " Genes with a zero count in every sample were dropped before the variance-stabilizing transformation."
            } else ""),
        preprocessed = paste(
            "Gene filtering was skipped because the input was a preprocessed expression",
            "matrix; all features were kept."),
        "The gene filter applied in this run is not recorded for this report."
    )

    short <- switch(mode,
        adaptive = sprintf("adaptive (KDE valley), threshold %s", threshold_label),
        fixed = sprintf("fixed threshold %s", threshold_label),
        deseq2_only = "zero-count genes removed (no expression threshold)",
        none = "none",
        preprocessed = "skipped (preprocessed input)",
        "not recorded"
    )

    # DESeq2 runs on count input only; it is a separate DE-stage step that the
    # deseq2_only mode name invites readers to confuse with the gene filter.
    de_note <- if (mode %in% c("adaptive", "fixed", "deseq2_only", "none")) {
        paste(
            "Separately, when DESeq2 performs the differential-expression testing, its",
            "independent filtering does not change which genes were kept above; it sets",
            "adjusted p-values to NA for genes with low mean normalized counts."
        )
    } else {
        NULL
    }

    list(
        mode = mode,
        applied = mode %in% c("adaptive", "fixed", "deseq2_only"),
        matrix = if (mode %in% c("adaptive", "fixed")) unit_label else NA_character_,
        group_col = if (mode %in% c("adaptive", "fixed")) group else NA_character_,
        threshold = if (has_thr) thr else NA_real_,
        threshold_label = threshold_label,
        rule = rule,
        short = short,
        methods = methods,
        de_note = de_note
    )
}

