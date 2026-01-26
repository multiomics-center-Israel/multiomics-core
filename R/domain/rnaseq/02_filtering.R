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

filter_features_dynamic <- function(norm_mat, meta, sample_col, group_col, threshold = 3) {
    stopifnot(sample_col %in% colnames(meta))
    stopifnot(group_col %in% colnames(meta))

    if (!all(colnames(norm_mat) %in% meta[[sample_col]])) {
        miss <- setdiff(colnames(norm_mat), meta[[sample_col]])
        stop("Samples missing in metadata: ", paste(miss, collapse = ", "))
    }
    meta2 <- meta[match(colnames(norm_mat), meta[[sample_col]]), , drop = FALSE]

    group_list <- split(meta2[[sample_col]], meta2[[group_col]])

    keep_feature <- function(expr_values) {
        for (grp in names(group_list)) {
            samp <- group_list[[grp]]
            sub_vals <- expr_values[samp]
            n_reps <- sum(!is.na(sub_vals))
            min_pass <- ceiling(n_reps / 2)
            if (sum(sub_vals >= threshold, na.rm = TRUE) >= min_pass) {
                return(TRUE)
            }
        }
        FALSE
    }
    keep_vec <- apply(norm_mat, 1, keep_feature)
    list(filtered = norm_mat[keep_vec, , drop = FALSE], keep_vec = keep_vec)
}
