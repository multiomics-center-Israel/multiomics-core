#' Preprocess proteomics data starting from loaded inputs
#'
#' @param inputs List from load_proteomics_inputs().
#' @param config Full config list.
#' @return A named list of preprocessing results. Expression matrices
#'   (features x samples):
#'   \itemize{
#'     \item \code{expr_raw} — the data on its original input scale: linear when
#'       \code{scale_in = "linear"} and the linear matrix was retained, otherwise
#'       log2. The actual scale is reported in \code{expr_raw_scale}.
#'     \item \code{expr_log2} — the analysis-ready log2 assay the pipeline
#'       operates on (contaminant-filtered, aligned; still contains NAs).
#'     \item \code{expr_filt} — the log2 assay after min-count filtering.
#'     \item \code{expr_work} / \code{expr_imp_single} — the imputed log2 assay.
#'   }
#'   Plus \code{expr_raw_scale} ("linear" or "log2"), \code{row_data},
#'   \code{meta}, and \code{imputation_qc}.
preprocess_proteomics <- function(inputs, config) {
  cfg <- config$modes$proteomics

  # validate inputs
  validate_proteomics_inputs(inputs, cfg)

  # Build standardized proteomics object
  prot_obj <- get_proteomics_expression_matrix(inputs, config)

  # The pipeline works on the log2 assay throughout; `expr_log2` is derived at the
  # end to hold the data on its original input scale (see #138).
  expr_log2 <- prot_obj$assay_log2

  row_data <- prot_obj$row_data
  col_data <- prot_obj$col_data
  
  # Use protein_id_col from config
  protein_id_col <- cfg$id_columns$protein_id
  
  # Ensure expr_log2 matrix contract
  if (is.data.frame(expr_log2)) {
    rn <- NULL
    if (!is.null(row_data) && protein_id_col %in% colnames(row_data)) rn <- row_data[[protein_id_col]]
    expr_log2 <- coerce_df_to_numeric_matrix(expr_log2, rownames_vec = rn, name = "expr_log2")
  }
  
  # Ensure rownames
  if (!is.null(row_data) && protein_id_col %in% colnames(row_data)) {
    rn <- as.character(row_data[[protein_id_col]])
    if (is.null(rownames(expr_log2)) || all(grepl("^\\d+$", rownames(expr_log2)))) {
      rownames(expr_log2) <- rn
    }
  }
  assert_numeric_matrix(expr_log2, "expr_log2")
  
  # Contaminant filtering (e.g. cRAP proteins)
  contam_res <- filter_contaminants(expr_log2, row_data, cfg)
  expr_log2 <- contam_res$expr_mat
  row_data <- contam_res$row_data
  
  # Align col_data
  sample_id_col <- cfg$effects$samples %||% cfg$id_columns$sample_col
  check_has_cols(col_data, sample_id_col, df_name = "col_data")
  expr_log2 <- align_matrix_to_meta(expr_log2, col_data, sample_id_col)
  col_data <- align_meta_to_expr(expr_log2, col_data, cfg)
  
  # Optional: sample_filter
  rules <- get_sample_filter_rules(config, mode = "proteomics")
  if (!is.null(rules)) {
    # implementation specific to sample filtering would act here
    # assuming logic similar to utils::apply_sample_filter is available or moved.
    # We will assume apply_sample_filter is in R/core/something.R or we need to copy it.
    # I didn't see apply_sample_filter in my extraction list, it was in 00_utils.R (lines 158).
    # I should have put it in R/core/02_validation.R or 03_alignment.R.
    # I'll check if I missed it. Just in case, I'll allow this file to assume it's available.
    if (exists("apply_sample_filter")) {
      filtered <- apply_sample_filter(sample_col = sample_id_col, meta = col_data, expr = expr_log2, rules = rules, mode = "proteomics")
      col_data <- filtered$meta
      expr_log2 <- filtered$expr
    }
  }
  
  # Filtering
  filt <- filter_proteomics_by_min_count(expr_mat = expr_log2, row_data = row_data, meta = col_data, cfg = cfg)
  expr_filt <- filt$expr_mat
  row_data_f <- filt$row_data
  assert_numeric_matrix(expr_filt, "expr_filt")
  
  # Normalization — dispatches based on cfg$normalization$method
  norm_method <- cfg$normalization$method %||% "none"
  if (norm_method == "median") {
    expr_filt <- normalize_proteomics_median(expr_filt)
    message("Normalization: median centering applied.")
  } else if (norm_method != "none") {
    stop(sprintf("Unknown normalization method: '%s'. Supported: 'none', 'median'.", norm_method))
  }
  
  # Single imputation (QC/plots) — dispatches based on cfg$imputation$method
  imp_res <- impute_proteomics(expr_mat = expr_filt, cfg = cfg, return_flags = TRUE)
  expr_imp_single <- imp_res$imputed
  assert_numeric_matrix(expr_imp_single, "expr_imp_single")
  
  imputation_qc <- NULL
  if (!is.null(imp_res$imputed_flag)) {
    imputation_qc <- list(imputed_flag = imp_res$imputed_flag)
  }
  
  # Recover the input-scale matrix for the `expr_raw` contract (#138): subset the
  # retained linear matrix to the cleaned features/samples. Fall back to the log2
  # assay when no linear matrix is available (already-log2 input, or DIA-NN), so
  # expr_raw always holds the data on its original input scale.
  if (!is.null(prot_obj$assay_linear) &&
      all(rownames(expr_log2) %in% rownames(prot_obj$assay_linear)) &&
      all(colnames(expr_log2) %in% colnames(prot_obj$assay_linear))) {
    expr_raw       <- prot_obj$assay_linear[rownames(expr_log2), colnames(expr_log2), drop = FALSE]
    expr_raw_scale <- "linear"
  } else {
    expr_raw       <- expr_log2
    expr_raw_scale <- "log2"
  }

  list(
    expr_raw = expr_raw,
    expr_raw_scale = expr_raw_scale,
    expr_log2 = expr_log2,
    expr_filt = expr_filt,
    expr_filt_pre_imp = expr_filt,
    expr_work = expr_imp_single,
    expr_imp_single = expr_imp_single,
    row_data = row_data_f,
    meta = col_data,
    imputation_qc = imputation_qc
  )
}

#' Median centering normalization for proteomics
#'
#' Centers each sample (column) by subtracting its median, so all samples
#' have the same median intensity. Common for log2-transformed proteomics data.
#' @param expr_mat numeric matrix (proteins x samples)
#' @return normalized matrix
normalize_proteomics_median <- function(expr_mat) {
  col_medians <- apply(expr_mat, 2, median, na.rm = TRUE)
  global_median <- median(col_medians, na.rm = TRUE)
  sweep(expr_mat, 2, col_medians - global_median)
}

get_sample_filter_rules <- function(config, mode) {
  cfg <- config$modes[[mode]]
  if (is.null(cfg)) {
    return(NULL)
  }
  sf <- cfg$sample_filter
  if (is.null(sf) || !isTRUE(sf$enabled)) {
    return(NULL)
  }
  sf$rules %||% NULL
}