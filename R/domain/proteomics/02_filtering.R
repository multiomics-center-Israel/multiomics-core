#' Filter proteomics data using min_count per condition
filter_proteomics_by_min_count <- function(expr_mat, row_data, meta, cfg, group_col = NULL) {
    eff <- cfg$effects
    sample_id_col <- eff$samples
    if (is.null(group_col)) group_col <- eff$color

    check_has_cols(meta, sample_id_col, df_name = "metadata")
    check_has_cols(meta, group_col, df_name = "metadata")

    sample_ids <- meta[[sample_id_col]]
    expr_mat <- expr_mat[, sample_ids, drop = FALSE]

    group <- meta[[group_col]]
    min_cfg <- cfg$filtering$min_count
    min_per_group <- extract_min_count(min_cfg, group)

    keep <- pass_filter(expr_mat = expr_mat, group = group, min_per_group = min_per_group)

    expr_mat_filt <- expr_mat[keep, , drop = FALSE]
    row_data_filt <- row_data[keep, , drop = FALSE]

    ids <- NULL
    if (!is.null(rownames(expr_mat))) {
        ids <- rownames(expr_mat)[keep]
    } else if (!is.null(cfg$id_columns$protein_id) && cfg$id_columns$protein_id %in% colnames(row_data)) {
        ids <- row_data[[cfg$id_columns$protein_id]][keep]
    }

    message("Filtering step: kept ", sum(keep), " of ", length(keep), " proteins.")

    list(expr_mat = expr_mat_filt, row_data = row_data_filt, ids = ids, keep = keep)
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

# --- Generic helpers (duplicated from 01_preprocessing.R to ensure domain independence or assume loaded) ---
# Since these are pure logic, keeping them here or in core is fine. I'll include them here for safety.

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

extract_min_count <- function(min_cfg, groups) {
    groups <- unique(as.character(groups))
    if (is.numeric(min_cfg) && length(min_cfg) == 1 && is.null(names(min_cfg))) {
        return(setNames(rep(min_cfg, length(groups)), groups))
    }
    if (!is.null(min_cfg$default)) {
        out <- setNames(rep(min_cfg$default, length(groups)), groups)
        overrides <- min_cfg[setdiff(names(min_cfg), "default")]
        for (nm in names(overrides)) {
            if (nm %in% names(out)) out[nm] <- overrides[[nm]]
        }
        return(out)
    }
    out <- unlist(min_cfg)
    names(out) <- names(min_cfg)
    out
}
