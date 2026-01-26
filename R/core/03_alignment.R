#' Align a DE table to a reference FeatureID order
align_de_table_by_feature_id <- function(tab, ref_ids, run_i = NA_integer_, contrast_name = NA_character_, id_col = "FeatureID") {
    if (is.null(tab)) stop(sprintf("DE table is NULL (run %s, contrast '%s').", run_i, contrast_name))
    if (!is.data.frame(tab)) tab <- as.data.frame(tab)

    if (!id_col %in% colnames(tab)) {
        stop(sprintf("DE table missing '%s' column (run %s, contrast '%s').", id_col, run_i, contrast_name))
    }

    ids <- tab[[id_col]]
    if (anyDuplicated(ids)) {
        dup <- unique(ids[duplicated(ids)])
        stop(sprintf("Duplicate %s values in DE table (run %s, contrast '%s'). Examples: %s", id_col, run_i, contrast_name, paste(head(dup, 3), collapse = ", ")))
    }

    idx <- match(ref_ids, ids)
    if (anyNA(idx)) {
        missing_ids <- ref_ids[is.na(idx)]
        stop(sprintf("Run %s: %d %s values missing in DE table for contrast '%s'. Examples: %s", run_i, length(missing_ids), id_col, contrast_name, paste(head(missing_ids, 3), collapse = ", ")))
    }

    tab[idx, , drop = FALSE]
}

#' Align metadata rows to match matrix column order
align_meta_to_matrix <- function(sample_ids, meta, sample_col) {
    stopifnot(sample_col %in% colnames(meta))
    meta <- as.data.frame(meta)

    idx <- match(sample_ids, meta[[sample_col]])
    if (any(is.na(idx))) {
        missing <- sample_ids[is.na(idx)]
        stop("Some samples are missing from metadata[[", sample_col, "]]: ", paste(head(missing), collapse = ", "))
    }

    meta_sub <- meta[idx, , drop = FALSE]
    rownames(meta_sub) <- sample_ids
    meta_sub
}

#' Apply sample map to rename matrix columns
apply_sample_map_to_colnames <- function(mat, sample_map, map_from, map_to) {
    stopifnot(!is.null(colnames(mat)))
    check_has_cols(sample_map, c(map_from, map_to), df_name = "sample_map")

    raw_names <- colnames(mat)
    new_names <- sample_map[[map_to]][match(raw_names, sample_map[[map_from]])]

    unmatched <- is.na(new_names)
    if (any(unmatched)) {
        warning("These matrix columns did not match any row in sample_map$", map_from, ": ", paste(head(raw_names[unmatched], 20), collapse = ", "))
    }

    colnames(mat) <- ifelse(unmatched, raw_names, new_names)
    mat
}

#' Align matrix columns to match meta sample order
align_matrix_to_meta <- function(mat, meta, sample_col, strict = TRUE) {
    check_has_cols(meta, sample_col, df_name = "meta")
    samp <- as.character(meta[[sample_col]])

    missing_in_mat <- setdiff(samp, colnames(mat))
    if (length(missing_in_mat) > 0) {
        msg <- paste0("meta contains samples missing in matrix columns: ", paste(head(missing_in_mat, 10), collapse = ", "))
        if (isTRUE(strict)) stop(msg) else warning(msg)
    }

    keep <- intersect(samp, colnames(mat))
    mat[, keep, drop = FALSE]
}

#' Align sample metadata to an expression matrix by sample ID
align_meta_to_expr <- function(expr_mat, meta, cfg) {
    sample_col <- get_sample_col(cfg)
    assert_expr_meta_alignment(expr_mat, meta, cfg, strict = FALSE)

    mi <- match(colnames(expr_mat), meta[[sample_col]])
    if (anyNA(mi)) {
        missing <- colnames(expr_mat)[is.na(mi)]
        stop("meta missing samples present in expr_mat: ", paste(head(missing, 5), collapse = ", "))
    }

    meta[mi, , drop = FALSE]
}

#' Align differential expression results to an expression matrix
align_de_to_expr <- function(de, expr_mat, contrast_name = NULL) {
    m <- match(rownames(expr_mat), rownames(de))
    if (anyNA(m)) {
        missing <- rownames(expr_mat)[is.na(m)]
        stop(paste0("topTable missing features", if (!is.null(contrast_name)) paste0(" for contrast '", contrast_name, "'"), ": ", paste(head(missing, 5), collapse = ", ")))
    }
    de[m, , drop = FALSE]
}

#' Align feature annotations to an expression matrix
align_feature_tbl_to_mat <- function(mat, feature_tbl, feature_id_col, annot_cols) {
    if (!is.matrix(mat)) stop("'mat' must be a matrix.")
    if (!is.data.frame(feature_tbl)) stop("'feature_tbl' must be a data.frame.")
    check_has_cols(feature_tbl, c(feature_id_col, annot_cols), df_name = "feature_tbl")

    idx <- match(rownames(mat), feature_tbl[[feature_id_col]])
    if (anyNA(idx)) {
        missing_ids <- rownames(mat)[is.na(idx)]
        stop("feature_tbl is missing metadata for ", sum(is.na(idx)), " features. First 5: ", paste(head(missing_ids, 5), collapse = ", "))
    }

    aligned_df <- feature_tbl[idx, annot_cols, drop = FALSE]
    rownames(aligned_df) <- rownames(mat)
    aligned_df
}

#' Align feature annotations to an expression matrix (Proteomics specific alias)
align_annotations_to_expr <- function(expr_mat, prot_tbl, protein_id_col, annot_cols) {
    align_feature_tbl_to_mat(expr_mat, prot_tbl, protein_id_col, annot_cols)
}

#' Apply sample filtering rules to (meta, expr)
apply_sample_filter <- function(sample_col, meta, expr, rules, mode = "omics", strict_cols = FALSE) {
    stopifnot(is.data.frame(meta))
    expr <- as.matrix(expr)

    if (is.null(rownames(meta)) || anyNA(rownames(meta))) stop("[", mode, "] meta must have non-NA rownames.")
    if (!identical(colnames(expr), rownames(meta))) stop("[", mode, "] colnames(expr) != rownames(meta).")

    if (is.null(rules) || length(rules) == 0) {
        return(list(meta = meta, expr = expr, info = list(n_before = nrow(meta), n_after = nrow(meta))))
    }

    keep <- rep(TRUE, nrow(meta))
    for (col in names(rules)) {
        if (!col %in% names(meta)) {
            if (isTRUE(strict_cols)) stop("sample_filter col not found: ", col) else next
        }
        vals <- unlist(rules[[col]], use.names = FALSE)
        keep <- keep & (meta[[col]] %in% vals)
    }

    meta2 <- meta[keep, , drop = FALSE]
    if (nrow(meta2) == 0) stop("sample_filter removed all samples.")
    expr2 <- expr[, rownames(meta2), drop = FALSE]

    list(meta = meta2, expr = expr2, info = list(n_before = nrow(meta), n_after = nrow(meta2)))
}
