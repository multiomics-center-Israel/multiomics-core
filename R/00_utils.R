# TODO [multiomics]:
# - Generalize build_final_results_* functions into a single generic helper:
#     * input: pre, summary_df, contrasts_df, id_col, expr_matrix
#     * output: final_results table with identical column structure across omics
# - Remove remaining proteomics-specific assumptions (e.g. Protein vs Gene)
# - Allow per-mode annotation adapters (e.g. protein annotations, gene annotations)
# - Keep Excel writers unchanged (they already support generic inputs)
# - Goal: fully unified Final_results API for RNA / proteomics / metabolomics

# TODO:
# - allow RNASeq clustering to export excel_order
# - optionally write order + zscore also to final_results.tsv (configurable)


# simple helper for "x %||% default"
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

get_sample_filter_rules <- function(config, mode) {
  cfg <- config$modes[[mode]]
  if (is.null(cfg)) return(NULL)
  sf <- cfg$sample_filter
  if (is.null(sf) || !isTRUE(sf$enabled)) return(NULL)
  sf$rules %||% NULL
}

ensure_dir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  invisible(dir)
}

save_tsv <- function(x, dir, filename) {
  ensure_dir(dir)
  path <- file.path(dir, filename)
  readr::write_tsv(x, path)
  path
}

save_tsv_path <- function(x, path) {
  ensure_dir(dirname(path))
  readr::write_tsv(x, path)
  path
}

resolve_project_path <- function(config, ...) {
  root <- config$project$dir %||% stop("config$project$dir is missing.")
  file.path(root, ...)
}

resolve_raw_path <- function(config, rel_path) {
  base <- config$paths$raw %||% "data"
  resolve_project_path(config, base, rel_path)
}

resolve_out_path <- function(config, rel_path = "") {
  base <- config$paths$out %||% "outputs"
  resolve_project_path(config, base, rel_path)
}

#' Align a DE table to a reference FeatureID order (fail-fast)
#'
#' @param tab A data.frame containing a FeatureID column.
#' @param ref_ids Character vector of FeatureID in the desired order.
#' @param run_i Integer; imputation run index (for error messages).
#' @param contrast_name Character; contrast name (for error messages).
#' @param id_col Character; name of the ID column (default "FeatureID").
#' @return tab reordered to match ref_ids (same number of rows, same order).
align_de_table_by_feature_id <- function(tab,
                                         ref_ids,
                                         run_i = NA_integer_,
                                         contrast_name = NA_character_,
                                         id_col = "FeatureID") {
  if (is.null(tab)) {
    stop(sprintf("DE table is NULL (run %s, contrast '%s').", run_i, contrast_name))
  }
  if (!is.data.frame(tab)) tab <- as.data.frame(tab)
  
  if (!id_col %in% colnames(tab)) {
    stop(sprintf("DE table missing '%s' column (run %s, contrast '%s').",
                 id_col, run_i, contrast_name))
  }
  
  ids <- tab[[id_col]]
  
  # Duplicates in tab are ambiguous for match()
  if (anyDuplicated(ids)) {
    dup <- unique(ids[duplicated(ids)])
    stop(sprintf(
      "Duplicate %s values in DE table (run %s, contrast '%s'). Examples: %s",
      id_col, run_i, contrast_name, paste(head(dup, 3), collapse = ", ")
    ))
  }
  
  # Reorder by reference ids
  idx <- match(ref_ids, ids)
  
  if (anyNA(idx)) {
    missing_ids <- ref_ids[is.na(idx)]
    stop(sprintf(
      "Run %s: %d %s values missing in DE table for contrast '%s'. Examples: %s",
      run_i, length(missing_ids), id_col, contrast_name,
      paste(head(missing_ids, 3), collapse = ", ")
    ))
  }
  
  tab[idx, , drop = FALSE]
}

# Row-wise z-score of a matrix
zscore_rows <- function(mat) {
  mat <- as.matrix(mat)
  row_means <- rowMeans(mat, na.rm = TRUE)
  row_sds <- apply(mat, 1, stats::sd, na.rm = TRUE)
  row_sds[row_sds == 0 | is.na(row_sds)] <- 1
  
  scaled <- sweep(mat, 1, row_means, FUN = "-")
  scaled <- sweep(scaled, 1, row_sds, FUN = "/")
  scaled
}

# compute hierarchical order for rows 
get_hclust_row_order <- function(mat, row_distance = "correlation", hclust_method = "complete") {
  mat <- as.matrix(mat)
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  if (nrow(mat) < 2) return(rownames(mat))
  
  d <- switch(
    row_distance,
    correlation = stats::as.dist(1 - stats::cor(t(mat), use = "pairwise.complete.obs")),
    euclidean   = stats::dist(mat),
    stop("Unsupported row_distance: ", row_distance)
  )
  
  hc <- stats::hclust(d, method = hclust_method)
  rownames(mat)[hc$order]
}

# Align metadata rows to match matrix column order
align_meta_to_matrix <- function(sample_ids, meta, sample_col) {
  stopifnot(sample_col %in% colnames(meta))
  meta <- as.data.frame(meta)
  
  idx <- match(sample_ids, meta[[sample_col]])
  if (any(is.na(idx))) {
    missing <- sample_ids[is.na(idx)]
    stop("Some samples are missing from metadata[[", sample_col, "]]: ",
         paste(head(missing), collapse = ", "))
  }
  
  meta_sub <- meta[idx, , drop = FALSE]
  rownames(meta_sub) <- sample_ids
  meta_sub
}
#' Apply sample filtering rules to (meta, expr)
#'
#' @return list(meta = meta_filtered, expr = expr_filtered, info = list(n_before, n_after))
apply_sample_filter <- function(sample_col,
                                meta,
                                expr,
                                rules,
                                mode = "omics",
                                strict_cols = FALSE) {
  stopifnot(is.data.frame(meta))
  stopifnot(is.matrix(expr) || is.data.frame(expr))
  expr <- as.matrix(expr)
  
  if (is.null(rownames(meta)) || anyNA(rownames(meta))) {
    stop("[", mode, "] meta must have non-NA rownames (SampleIDs).")
  }
  if (is.null(colnames(expr))) {
    stop("[", mode, "] expr must have colnames (SampleIDs).")
  }
  
  # enforce alignment assumption
  if (!identical(colnames(expr), rownames(meta))) {
    stop("[", mode, "] apply_sample_filter requires identical colnames(expr) == rownames(meta).")
  }
  
  if (is.null(rules) || length(rules) == 0) {
    return(list(meta = meta, expr = expr, info = list(n_before = nrow(meta), n_after = nrow(meta))))
  }
  
  keep <- rep(TRUE, nrow(meta))
  
  for (col in names(rules)) {
    if (!col %in% names(meta)) {
      msg <- paste0("[", mode, "] sample_filter column not found in metadata: ", col, " (skipping)")
      if (isTRUE(strict_cols)) stop(msg) else warning(msg)
      next
    }
    
    vals <- rules[[col]]
    if (is.null(vals)) next
    vals <- unlist(vals, use.names = FALSE)
    
    keep <- keep & (meta[[col]] %in% vals)
  }
  
  n_before <- nrow(meta)
  meta2 <- meta[keep, , drop = FALSE]
  
  if (nrow(meta2) == 0) {
    stop("[", mode, "] sample_filter removed all samples. Check your rules.")
  }
  
  expr2 <- expr[, meta2[[sample_col]], drop = FALSE]
  
  message(sprintf("[%s] sample_filter kept %d/%d samples.", mode, nrow(meta2), n_before))
  
  list(meta = meta2, expr = expr2, info = list(n_before = n_before, n_after = nrow(meta2)))
}
# Minimal contracts for a standard omics object.
assert_pre_contract <- function(pre, stage = c("proteomics", "rna", "metabolomics")) {
  stage <- match.arg(stage)
  
  stopifnot(is.list(pre))
  
  base_fields <- c(
    "expr_raw",
    "expr_filt",
    "expr_work",
    "meta",
    "row_data"
  )
  
  stage_fields <- switch(
    stage,
    proteomics = c("expr_imp_single"),
    rna        = character(0),
    metabolomics = character(0)
  )
  
  required <- c(base_fields, stage_fields)
  missing <- setdiff(required, names(pre))
  
  if (length(missing) > 0) {
    stop(
      sprintf(
        "Preprocess contract failed for %s. Missing fields: %s",
        stage,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

assert_meta_contract <- function(meta, sample_col) {
  if (!is.data.frame(meta)) stop("meta must be a data.frame.")
  check_has_cols(meta, sample_col, df_name = "meta")
  samp <- as.character(meta[[sample_col]])
  if (any(!nzchar(samp))) stop("meta$", sample_col, " contains empty sample IDs.")
  if (anyDuplicated(samp) > 0) stop("meta$", sample_col, " contains duplicated sample IDs.")
  invisible(TRUE)
}

assert_expr_meta_alignment <- function(expr_mat, meta, cfg, strict = TRUE) {
  sample_col <- get_sample_col(cfg)
  
  assert_numeric_matrix(expr_mat, "expr_mat")
  assert_meta_contract(meta, sample_col)
  
  expr_ids <- colnames(expr_mat)
  meta_ids <- as.character(meta[[sample_col]])
  
  missing_in_expr <- setdiff(meta_ids, expr_ids)
  if (length(missing_in_expr) > 0) {
    stop(
      "expr_mat is missing samples from meta (sample_col='", sample_col, "'): ",
      paste(head(missing_in_expr, 10), collapse = ", "),
      if (length(missing_in_expr) > 10) sprintf(" ... (+%d more)", length(missing_in_expr) - 10) else ""
    )
  }
  
  # expr -> meta (עודפים ב-matrix)
  if (isTRUE(strict)) {
    extra_in_expr <- setdiff(expr_ids, meta_ids)
    if (length(extra_in_expr) > 0) {
      stop(
        "expr_mat has samples not present in meta (sample_col='", sample_col, "'): ",
        paste(head(extra_in_expr, 10), collapse = ", "),
        if (length(extra_in_expr) > 10) sprintf(" ... (+%d more)", length(extra_in_expr) - 10) else ""
      )
    }
  }
  
  invisible(TRUE)
}

assert_omics_obj <- function(obj, stage = c("raw", "work"), sample_col) {
  stage <- match.arg(stage)
  
  if (!is.list(obj)) stop("omics obj must be a list.")
  if (is.null(obj$assay) || !is.list(obj$assay)) stop("omics obj missing $assay list.")
  if (is.null(obj$col_data) || !is.data.frame(obj$col_data)) stop("omics obj missing $col_data (data.frame).")
  
  if (missing(sample_col) || is.null(sample_col) || !nzchar(sample_col)) {
    stop("assert_omics_obj: 'sample_col' must be provided (no default).")
  }
  
  mat <- obj$assay[[stage]]
  if (is.null(mat)) stop("omics obj missing assay$", stage)
  
  assert_numeric_matrix(mat, paste0("assay$", stage))
  assert_meta_contract(obj$col_data, sample_col)
  
  meta_ids <- as.character(obj$col_data[[sample_col]])
  if (!identical(meta_ids, colnames(mat))) {
    # give a helpful diff
    missing_in_mat <- setdiff(meta_ids, colnames(mat))
    extra_in_mat   <- setdiff(colnames(mat), meta_ids)
    
    stop(
      "assay$", stage, " and col_data are not aligned.\n",
      "- Expected colnames(mat) to match col_data[['", sample_col, "']] exactly.\n",
      if (length(missing_in_mat)) paste0("  Missing in mat: ", paste(head(missing_in_mat, 10), collapse = ", "), "\n") else "",
      if (length(extra_in_mat))   paste0("  Extra in mat: ", paste(head(extra_in_mat, 10), collapse = ", "), "\n") else "",
      "Tip: call align_meta_to_expr(mat, col_data, sample_col) upstream."
    )
  }
  
  invisible(TRUE)
}


#' Align feature annotations to an expression matrix by feature ID
#'
#' Fail-fast helper that subsets and reorders `prot_tbl` rows to match the
#' row order of an expression matrix (`expr_mat`) using a feature ID column.
#'
#' The function:
#' - checks `protein_id_col` exists in `prot_tbl`
#' - checks all requested `annot_cols` exist in `prot_tbl`
#' - ensures every `rownames(expr_mat)` feature is present in `prot_tbl`
#' - returns an annotation data.frame aligned to `rownames(expr_mat)`
#'
#' @param expr_mat A matrix (features × samples) with feature IDs in `rownames(expr_mat)`.
#' @param prot_tbl A data.frame containing feature annotations (e.g. DIA-NN protein table).
#' @param protein_id_col Character. Column in `prot_tbl` containing the feature IDs.
#' @param annot_cols Character vector. Columns to return (should include `protein_id_col`).
#'
#' @return A data.frame of annotations with `nrow == nrow(expr_mat)` and rows aligned to `rownames(expr_mat)`.
#' @export
align_annotations_to_expr <- function(expr_mat, prot_tbl, protein_id_col, annot_cols) {
  if (!is.matrix(expr_mat)) stop("expr_mat must be a matrix.")
  if (!is.data.frame(prot_tbl)) stop("prot_tbl must be a data.frame.")
  if (!(protein_id_col %in% colnames(prot_tbl))) stop("prot_tbl missing protein_id_col: ", protein_id_col)
  
  missing_cols <- setdiff(annot_cols, colnames(prot_tbl))
  if (length(missing_cols) > 0) {
    stop("prot_tbl is missing annotation columns: ", paste(missing_cols, collapse = ", "))
  }
  
  idx <- match(rownames(expr_mat), prot_tbl[[protein_id_col]])
  if (anyNA(idx)) {
    missing_ids <- rownames(expr_mat)[is.na(idx)]
    stop(
      "prot_tbl missing ", sum(is.na(idx)), " proteins (protein_id_col='", protein_id_col, "'). ",
      "Example: ", paste(head(missing_ids, 5), collapse = ", ")
    )
  }
  
  prot_tbl[idx, annot_cols, drop = FALSE]
}


get_mode_cfg <- function(config, mode) {
  cfg <- config$modes[[mode]]
  if (is.null(cfg)) stop("Missing config$modes$", mode)
  cfg
}

get_sample_col <- function(cfg) {
  # source of truth: cfg$effects$samples
  sample_col <- cfg$effects$samples %||% cfg$id_columns$sample_col
  if (is.null(sample_col) || !nzchar(sample_col)) stop("Missing sample column in cfg$effects$samples / cfg$id_columns$sample_col")
  sample_col
}

# Rename matrix columns using sample_map (from -> to).
# Unmatched columns keep original name (warn).
apply_sample_map_to_colnames <- function(mat, sample_map, map_from, map_to) {
  stopifnot(!is.null(colnames(mat)))
  check_has_cols(sample_map, c(map_from, map_to), df_name = "sample_map")
  
  raw_names <- colnames(mat)
  new_names <- sample_map[[map_to]][match(raw_names, sample_map[[map_from]])]
  
  unmatched <- is.na(new_names)
  if (any(unmatched)) {
    warning(
      "These matrix columns did not match any row in sample_map$", map_from, ": ",
      paste(head(raw_names[unmatched], 20), collapse = ", "),
      if (sum(unmatched) > 20) sprintf(" ... (+%d more)", sum(unmatched) - 20) else ""
    )
  }
  
  colnames(mat) <- ifelse(unmatched, raw_names, new_names)
  mat
}

# Reorder columns of mat to match meta sample order.
# If strict=TRUE -> fail if any meta samples are missing in mat.
align_matrix_to_meta <- function(mat, meta, sample_col, strict = TRUE) {
  check_has_cols(meta, sample_col, df_name = "meta")
  samp <- as.character(meta[[sample_col]])
  
  missing_in_mat <- setdiff(samp, colnames(mat))
  if (length(missing_in_mat) > 0) {
    msg <- paste0(
      "meta contains samples missing in matrix columns: ",
      paste(head(missing_in_mat, 10), collapse = ", "),
      if (length(missing_in_mat) > 10) sprintf(" ... (+%d more)", length(missing_in_mat) - 10) else ""
    )
    if (isTRUE(strict)) stop(msg) else warning(msg)
  }
  
  keep <- intersect(samp, colnames(mat))
  mat[, keep, drop = FALSE]
}

#' Align sample metadata to an expression matrix by sample ID
#'
#' Fail-fast helper that reorders `meta` rows to match the column order of
#' an expression matrix (`expr_mat`). This prevents silent misalignment bugs
#' in downstream modeling (e.g., limma).
#'
#' @param expr_mat A numeric matrix (features × samples) with sample IDs in `colnames(expr_mat)`.
#' @param meta A data.frame with one row per sample.
#' @param sample_col Character. Column name in `meta` holding sample IDs.
#'
#' @return A reordered `meta` data.frame with rows aligned to `colnames(expr_mat)`.
#' @export

align_meta_to_expr <- function(expr_mat, meta, cfg) {
  sample_col <- get_sample_col(cfg)
  assert_expr_meta_alignment(expr_mat, meta, cfg, strict = FALSE)
  
  mi <- match(colnames(expr_mat), meta[[sample_col]])
  if (anyNA(mi)) {
    missing <- colnames(expr_mat)[is.na(mi)]
    stop(
      "meta missing ", sum(is.na(mi)), " samples present in expr_mat (sample_col='", sample_col, "'). ",
      "Example: ", paste(head(missing, 5), collapse = ", ")
    )
  }
  
  meta[mi, , drop = FALSE]
}


#' Align differential expression results to an expression matrix
#'
#' Fail-fast helper that reorders a `topTable()` result to match the feature
#' order of the input expression matrix. This prevents silent mis-annotation
#' when binding DE statistics with feature annotations.
#' 
#' @param de A data.frame returned by `limma::topTable()` (rownames are feature IDs).
#' @param expr_mat A matrix (features × samples) with feature IDs in `rownames(expr_mat)`.
#' @param contrast_name Optional character. Used only to improve error messages.
#'
#' @return `de` reordered to match `rownames(expr_mat)`.
#' @export
align_de_to_expr <- function(de, expr_mat, contrast_name = NULL) {
  m <- match(rownames(expr_mat), rownames(de))
  if (anyNA(m)) {
    missing <- rownames(expr_mat)[is.na(m)]
    msg <- paste0(
      "topTable missing ", sum(is.na(m)), " features",
      if (!is.null(contrast_name)) paste0(" for contrast '", contrast_name, "'") else "",
      ". Example: ", paste(head(missing, 5), collapse = ", ")
    )
    stop(msg)
  }
  de[m, , drop = FALSE]
}

align_feature_tbl_to_mat <- function(mat, feature_tbl, feature_id_col, annot_cols) {
  
  # 1. Guards
  if (!is.matrix(mat)) stop("'mat' must be a matrix.")
  if (!is.data.frame(feature_tbl)) stop("'feature_tbl' must be a data.frame.")
  
  if (!(feature_id_col %in% colnames(feature_tbl))) {
    stop("feature_tbl is missing id column: '", feature_id_col, "'")
  }
  
  missing_cols <- setdiff(annot_cols, colnames(feature_tbl))
  if (length(missing_cols) > 0) {
    stop("feature_tbl is missing annotation columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # 2. Alignment Logic
  # Align feature_tbl rows to match rownames(mat)
  idx <- match(rownames(mat), feature_tbl[[feature_id_col]])
  
  if (anyNA(idx)) {
    missing_ids <- rownames(mat)[is.na(idx)]
    stop(
      "feature_tbl is missing metadata for ", sum(is.na(idx)), " features found in the matrix.\n",
      "ID Column: '", feature_id_col, "'\n",
      "First 5 missing IDs: ", paste(head(missing_ids, 5), collapse = ", ")
    )
  }
  
  # 3. Output
  aligned_df <- feature_tbl[idx, annot_cols, drop = FALSE]
  rownames(aligned_df) <- rownames(mat)
  
  return(aligned_df)
}

# Create a standard omics object (generic container)
init_omics_obj <- function(mode, engine = NULL, assay_raw, col_data, row_data = NULL, assay_work = NULL, params = list()) {
  list(
    mode   = mode,
    engine = engine,
    assay  = list(
      raw  = assay_raw,
      work = assay_work
    ),
    row_data = row_data,
    col_data = col_data,
    params   = params
  )
}

#' Assert a strict numeric matrix contract (features × samples)
#'
#' Fail-fast validator for expression matrices used throughout the pipeline.
#' Ensures the object is a numeric matrix with **non-empty, unique** row and
#' column names. This catches silent alignment bugs early (e.g. lost rownames,
#' duplicated sample IDs, character coercions).
#'
#' @param x An object expected to be a numeric matrix.
#' @param name Character. A label used in error messages (helps debugging).
#'
#' @return Invisibly returns TRUE if all checks pass.
#'
#' @examples
#' \dontrun{
#' assert_numeric_matrix(expr_mat, "expr_filt")
#' }
assert_numeric_matrix <- function(x, name = "x") {
  if (!is.matrix(x)) stop(sprintf("`%s` must be a matrix.", name))
  if (!is.numeric(x)) stop(sprintf("`%s` must be numeric.", name))
  if (is.null(rownames(x)) || anyNA(rownames(x)) || any(rownames(x) == "")) {
    stop(sprintf("`%s` must have non-empty rownames (feature IDs).", name))
  }
  if (is.null(colnames(x)) || anyNA(colnames(x)) || any(colnames(x) == "")) {
    stop(sprintf("`%s` must have non-empty colnames (sample IDs).", name))
  }
  if (anyDuplicated(rownames(x)) > 0) stop(sprintf("`%s` has duplicated rownames.", name))
  if (anyDuplicated(colnames(x)) > 0) stop(sprintf("`%s` has duplicated colnames.", name))
  invisible(TRUE)
}
#' Coerce a data.frame to a numeric matrix (optionally set rownames)
#'
#' Converts a data.frame (or already-matrix) to a matrix with storage mode
#' forced to numeric. This is useful for I/O sources (e.g. DIA-NN exports)
#' where columns may be read as character due to blanks or formatting.
#'
#' If `rownames_vec` is provided, it is used to set matrix rownames (feature IDs),
#' and its length must exactly match the number of rows.
#'
#' @param df A data.frame or matrix.
#' @param rownames_vec Optional vector of feature IDs to use as rownames.
#' @param name Character. A label used in error messages (helps debugging).
#'
#' @return A numeric matrix.
#'
#' @examples
#' \dontrun{
#' m <- coerce_df_to_numeric_matrix(df, rownames_vec = row_data$Protein.Group)
#' }
coerce_df_to_numeric_matrix <- function(df, rownames_vec = NULL, name = "df") {
  if (is.matrix(df)) return(df)
  if (!is.data.frame(df)) stop(sprintf("`%s` must be a data.frame or matrix.", name))
  
  m <- as.matrix(df)
  suppressWarnings(storage.mode(m) <- "numeric")
  
  if (!is.null(rownames_vec)) {
    if (length(rownames_vec) != nrow(m)) {
      stop(sprintf("`rownames_vec` length (%d) != nrow(%s) (%d).",
                   length(rownames_vec), name, nrow(m)))
    }
    rownames(m) <- as.character(rownames_vec)
  }
  m
}
#' Build output root directory for a specific project run
#'
#' Creates a stable "run folder" name to separate outputs between analysis rounds.
#' Example: outputs/Results_E_Pick_Analysis_02
#'
#' @param config Full config list.
#' @return Character path to the run output directory.
get_run_out_dir <- function(config) {
  proj  <- config$project$name %||% "Project"
  round <- config$project$analysis_round %||% "Analysis"
  resolve_out_path(config, sprintf("Results_%s_%s", proj, round))
}

#' Get output directory for a specific omics mode
#'
#' @param out_dir Base run directory (from get_run_out_dir)
#' @param mode    Omics mode name (e.g. "proteomics", "rnaseq")
#'
#' @return Path to mode-specific output directory
get_mode_out_dir <- function(out_dir, mode) {
  stopifnot(is.character(out_dir), length(out_dir) == 1)
  stopifnot(is.character(mode), length(mode) == 1)
  file.path(out_dir, mode)
}


#' Create legacy-style output folder structure for a run
#'
#' Mirrors the original Neat proteomics output tree inside the run folder:
#'   Datasets/, Diagnostic_plots/, Clustering/, Enrichment/, GSEA_enrichment/
#'
#' @param out_dir Run output root directory (e.g., outputs/Results_E_Pick_Analysis_02).
#' @return Named list of important subdirectories.
create_legacy_output_dirs <- function(out_dir, create = TRUE) {
  stopifnot(is.character(out_dir), length(out_dir) == 1)
  
  dirs <- list(
    datasets         = file.path(out_dir, "Datasets"),
    diagnostic_plots = file.path(out_dir, "Diagnostic_plots"),
    clustering       = file.path(out_dir, "Clustering"),
    enrichment       = file.path(out_dir, "Enrichment"),
    gsea_enrichment  = file.path(out_dir, "GSEA_enrichment")
  )
  
  if (isTRUE(create)) {
    for (d in unique(unname(dirs))) ensure_dir(d)
  }
  
  dirs
}
# multiomics-core — run execution info utilities
#
# Goal:
# Snapshot the exact inputs needed to reproduce a run:
# - config snapshot
# - git commit + timestamp + session info (existing code)
# - optionally copy _targets.R
# - optionally save project workspace (project.Rdata) using save.image()
#
# Notes:
# - safe both interactively and via {targets}
# - writes under <run_dir>/execution_info/
# - returns a character vector of file paths (targets-friendly)

write_execution_info <- function(config,
                                 run_dir,
                                 config_path = NULL,
                                 manifest_df = NULL) {
  
  exec_dir <- file.path(run_dir, "execution_info")
  dir.create(exec_dir, recursive = TRUE, showWarnings = FALSE)
  
  written_files <- character(0)
  
  # --- config snapshot (full resolved list) ---
  cfg_path <- file.path(exec_dir, "config_used.yaml")
  yaml::write_yaml(config, cfg_path)
  written_files <- c(written_files, cfg_path)
  
  # --- record which config file was used (path) ---
  if (!is.null(config_path)) {
    cfgp_path <- file.path(exec_dir, "config_path.txt")
    writeLines(normalizePath(config_path, winslash = "/", mustWork = FALSE), cfgp_path)
    written_files <- c(written_files, cfgp_path)
  }
  
  # --- git commit ---
  git_path <- file.path(exec_dir, "git_commit.txt")
  git_commit <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) "UNKNOWN")
  writeLines(git_commit, git_path)
  written_files <- c(written_files, git_path)
  
  # --- session info ---
  sess_path <- file.path(exec_dir, "sessionInfo.txt")
  writeLines(capture.output(sessionInfo()), sess_path)
  written_files <- c(written_files, sess_path)
  
  # --- timestamp ---
  time_path <- file.path(exec_dir, "timestamp.txt")
  writeLines(as.character(Sys.time()), time_path)
  written_files <- c(written_files, time_path)
  
  # --- (A) copy renv.lock ---
  renv_lock_src <- "renv.lock"
  if (file.exists(renv_lock_src)) {
    renv_lock_dst <- file.path(exec_dir, "renv.lock")
    ok <- file.copy(renv_lock_src, renv_lock_dst, overwrite = TRUE)
    if (isTRUE(ok)) written_files <- c(written_files, renv_lock_dst)
  } else {
    warning("renv.lock not found in project root; run may be harder to reproduce.")
  }
  
  # --- (B) copy _targets.R ---
  save_targets <- isTRUE(config$execution_info$save_targets_file %||% TRUE)
  if (save_targets) {
    src_targets <- "_targets.R"
    if (!file.exists(src_targets)) stop("Cannot find _targets.R in project root.")
    dst_targets <- file.path(exec_dir, "_targets.R")
    ok <- file.copy(src_targets, dst_targets, overwrite = TRUE)
    if (!isTRUE(ok)) stop("Failed to copy _targets.R into execution_info.")
    written_files <- c(written_files, dst_targets)
  }
  
  
  # --- (C) save targets manifest (recommended to pass in from a target) ---
  if (!is.null(manifest_df)) {
    man_path <- file.path(exec_dir, "targets_manifest.csv")
    utils::write.csv(manifest_df, man_path, row.names = FALSE)
    written_files <- c(written_files, man_path)
  }
  
  # --- (D) optional workspace snapshot (project.Rdata) ---
  save_ws <- isTRUE(config$execution_info$save_workspace %||% FALSE)
  if (save_ws) {
    ws_name <- config$execution_info$workspace_filename %||% "project.Rdata"
    ws_path <- file.path(exec_dir, ws_name)
    save.image(file = ws_path)
    if (!file.exists(ws_path)) stop("save.image() did not create project workspace file.")
    written_files <- c(written_files, ws_path)
  }
  
  written_files
}

# Optional helper if you want an explicit function name to create project.Rdata
# (still uses save.image(), but keeps callsite clearer).
write_project_rdata <- function(run_dir, ..., filename = "project.Rdata") {
  exec_dir <- file.path(run_dir, "execution_info")
  dir.create(exec_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(exec_dir, filename)
  
  objects_list <- list(...)
  save(list = names(objects_list), envir = list2env(objects_list), file = out)
  
  if (!file.exists(out)) stop("save() did not create project.Rdata.")
  out
}


#' Validate global pipeline configuration
#'
#' Runs fail-fast checks on the YAML-derived config object.
#' Ensures critical parameters exist and have valid types/values before execution.
#'
#' @param config Full config list (from config.yaml).
#' @return Invisibly TRUE if all checks pass; otherwise stops with an error.
validate_config <- function(config) {
  assert_named_list(config, "config")
  
  # Validate basic required fields structure
  assert_named_list(config$paths,  "config$paths")
  assert_named_list(config$modes,  "config$modes")
  assert_named_list(config$params, "config$params")
  
  # Validate specific scalar values
  assert_scalar_chr(config$paths$out, "config$paths$out")
  assert_scalar_num(config$params$seed, "config$params$seed")
  
  # Mode-specific validations (only if the mode is enabled/present)
  if (!is.null(config$modes$proteomics))   validate_proteomics_config(config$modes$proteomics)
  if (!is.null(config$modes$rna))          validate_rna_config(config$modes$rna)
  if (!is.null(config$modes$metabolomics)) validate_metabolomics_config(config$modes$metabolomics)
  if (!is.null(config$modes$lipidomics))   validate_lipidomics_config(config$modes$lipidomics)
  
  invisible(TRUE)
}










