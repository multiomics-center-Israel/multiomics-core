# R/modules/metabolomics/00_mod_preprocessing.R
#
# Target-ready wrappers for the new met_* preprocessing pipeline.
#
# Each function is called directly from the {targets} target bodies in
# pipe_metabolomics().  They orchestrate domain functions from:
#   R/domain/metabolomics/00_inputs.R      (load / parse)
#   R/domain/metabolomics/08_missingness.R (classify / filter)
#   R/domain/metabolomics/09_imputation_met.R (impute)
#   R/domain/metabolomics/10_drift_correction.R (LOESS)
#   R/domain/metabolomics/01_normalization.R   (norm_*, transform_metab)
#   R/core/08_qc.R                         (qc_pca_scatter, norm_boxplot)
#   R/core/06_plots.R                      (qc_missingness_heatmap via 08_qc.R)


# ==============================================================================
# mod_met_raw — load, parse, apply sample filter, align metadata
# ==============================================================================

#' Load and parse raw metabolomics inputs
#'
#' Dispatches to the appropriate format parser, applies optional sample filter
#' (QC/blank exclusion), and aligns metadata to the expression matrix.
#'
#' @param inp    Pre-loaded inputs list from \code{load_metabolomics_inputs()}.
#' @param config Full pipeline config list.
#' @return list with: \code{expr_raw}, \code{meta}, \code{row_data},
#'   \code{sample_col}, \code{format}.
#'
mod_met_raw <- function(inp, config) {
  cfg        <- config$modes$metabolomics
  fmt        <- inp$format %||% cfg$input$format %||% "cd_raw"
  sample_col <- cfg$effects$samples %||% "sample_id"
  map_from   <- cfg$id_columns$map_from
  map_to     <- cfg$id_columns$map_to

  # Effective per-file format: for multi_level each level is parsed with
  # level_format; for single-file modes level_format == fmt.
  level_format <- if (fmt == "multi_level")
    cfg$input[["level_format"]] %||% "cd_raw"
  else
    fmt

  # processed_wide: columns in the raw data ARE the sample identifiers;
  # parse_processed_wide matches them against metadata IDs -> rename BEFORE
  # parsing.  For multi_level we rename each level's data_df individually.
  inp_data <- inp$data
  if (!is.null(inp$sample_map) && level_format == "processed_wide") {
    if (is.null(map_from) || is.null(map_to))
      stop("mod_met_raw: id_columns$map_from and id_columns$map_to are required ",
           "when files$sample_map is set")
    if (fmt == "multi_level") {
      inp_data <- lapply(inp_data, function(item) {
        item$data_df <- apply_sample_map_to_colnames(
          item$data_df, inp$sample_map, map_from, map_to)
        item
      })
    } else {
      inp_data <- apply_sample_map_to_colnames(inp_data, inp$sample_map, map_from, map_to)
    }
  }

  parsed <- switch(fmt,
    cd_raw         = parse_cd_raw(inp_data, cfg),
    processed_wide = parse_processed_wide(inp_data, cfg, inp$metadata),
    multi_level    = parse_multi_level(inp_data, cfg, inp$metadata),
    stop("mod_met_raw: unsupported format: '", fmt, "'")
  )

  expr_raw <- parsed$expr_raw
  row_data <- parsed$row_data

  # cd_raw: parse_cd_raw extracts sample IDs from Area: column names via regex;
  # the external sample_map remaps those extracted IDs -> rename AFTER parsing.
  # Applies to both single-file cd_raw and multi_level with level_format cd_raw.
  if (!is.null(inp$sample_map) && level_format == "cd_raw") {
    if (is.null(map_from) || is.null(map_to))
      stop("mod_met_raw: id_columns$map_from and id_columns$map_to are required ",
           "when files$sample_map is set")
    message("mod_met_raw: applying sample map to expression matrix column names.")
    expr_raw <- apply_sample_map_to_colnames(expr_raw, inp$sample_map, map_from, map_to)
  }

  # Align metadata; build minimal stub if none provided
  meta <- inp$metadata %||% build_minimal_meta(colnames(expr_raw))
  meta <- align_meta_to_matrix(colnames(expr_raw), meta, sample_col)

  # Apply optional sample filter (QC/blank exclusion)
  rules <- get_sample_filter_rules_metab(cfg)
  if (!is.null(rules)) {
    keep_ids <- apply_sample_filter_metab(colnames(expr_raw), meta, rules, sample_col)
    if (length(keep_ids) < ncol(expr_raw)) {
      message(sprintf(
        "mod_met_raw: sample filter removed %d sample(s); retaining %d.",
        ncol(expr_raw) - length(keep_ids), length(keep_ids)
      ))
      expr_raw <- expr_raw[, keep_ids, drop = FALSE]
      meta     <- meta[meta[[sample_col]] %in% keep_ids, , drop = FALSE]
    }
  }

  # Normalise row_data rownames to feature_id so all downstream subsetting
  # (mod_met_filtered, etc.) can safely use rownames(row_data).
  if (!is.null(row_data) && !is.null(row_data$feature_id)) {
    rownames(row_data) <- row_data$feature_id
  }

  list(
    expr_raw      = expr_raw,
    meta          = meta,
    row_data      = row_data,
    sample_col    = sample_col,
    format        = fmt,
    duplicate_log = parsed$duplicate_log
  )
}


# ==============================================================================
# mod_met_missingness_stats — classify MNAR / MAR per feature
# ==============================================================================

#' Classify features as MNAR or MAR based on missingness vs. intensity correlation
#'
#' @param raw    List returned by \code{mod_met_raw()}.
#' @param config Full pipeline config list.
#' @return List from \code{compute_missingness_stats()} plus
#'   \code{mnar_mask} (named logical: TRUE = MNAR),
#'   \code{mar_mask}  (named logical: TRUE = MAR).
#'
mod_met_missingness_stats <- function(raw, config) {
  pre_cfg  <- config$modes$metabolomics$preprocessing %||% list()
  mnar_thr <- pre_cfg$mnar_threshold %||% 0.3

  result <- compute_missingness_stats(
    mat                = raw$expr_raw,
    meta               = raw$meta,
    sample_col         = raw$sample_col,
    feat_miss_threshold = pre_cfg$feat_missing_threshold   %||% 0.20,
    samp_miss_threshold = pre_cfg$sample_missing_threshold %||% 0.30,
    mnar_cor_threshold  = -abs(mnar_thr)
  )

  stats_df <- result$stats_df
  mnar_mask <- stats::setNames(stats_df$mnar_class == "MNAR", stats_df$feature_id)
  mar_mask  <- stats::setNames(
    stats_df$mnar_class == "MAR",
    stats_df$feature_id
  )

  c(result, list(mnar_mask = mnar_mask, mar_mask = mar_mask))
}


# ==============================================================================
# mod_met_missingness_plot — save binary heatmap as PNG, return file path
# ==============================================================================

#' Save the missingness heatmap to disk
#'
#' @param miss_stats List returned by \code{mod_met_missingness_stats()}.
#' @param raw        List returned by \code{mod_met_raw()}.
#' @param out_dir    Mode output directory (met_out_dir).
#' @param config     Full pipeline config list.
#' @return Character scalar: path to saved PNG (for \code{format = "file"} target).
#'
mod_met_missingness_plot <- function(miss_stats, raw, out_dir, config) {
  diag_dir <- file.path(out_dir, "diagnostic_plots")
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(diag_dir, "missingness_heatmap.png")

  qc_missingness_heatmap(raw$expr_raw, miss_stats$stats_df, out_file)

  out_file
}


# ==============================================================================
# mod_met_filtered — apply missingness thresholds
# ==============================================================================

#' Filter features and samples by missingness thresholds
#'
#' @param raw    List returned by \code{mod_met_raw()}.
#' @param config Full pipeline config list.
#' @return list with: \code{mat} (expr_filt), \code{meta}, \code{row_data},
#'   \code{dropped_features}, \code{dropped_samples}.
#'
mod_met_filtered <- function(raw, config) {
  pre_cfg <- config$modes$metabolomics$preprocessing %||% list()

  filt <- filter_by_missingness(
    mat              = raw$expr_raw,
    meta             = raw$meta,
    sample_col       = raw$sample_col,
    feat_threshold   = pre_cfg$feat_missing_threshold   %||% 0.20,
    sample_threshold = pre_cfg$sample_missing_threshold %||% 0.30
  )

  row_data <- raw$row_data
  if (!is.null(row_data)) {
    keep     <- intersect(rownames(filt$mat), rownames(row_data))
    row_data <- row_data[keep, , drop = FALSE]
  }

  list(
    mat              = filt$mat,
    meta             = filt$meta,
    row_data         = row_data,
    dropped_features = filt$dropped_features,
    dropped_samples  = filt$dropped_samples
  )
}


# ==============================================================================
# mod_met_imputed — MNAR (min/2) + MAR (KNN) imputation
# ==============================================================================

#' Impute missing values (MNAR -> min/2, MAR -> KNN)
#'
#' @param filtered   List returned by \code{mod_met_filtered()}.
#' @param config     Full pipeline config list.
#' @return list with: \code{mat}, \code{imputed_flag}, \code{meta},
#'   \code{row_data}, \code{info}.
#'
mod_met_imputed <- function(filtered, config) {
  pre_cfg  <- config$modes$metabolomics$preprocessing %||% list()
  mnar_thr <- pre_cfg$mnar_threshold %||% 0.3
  knn_k    <- as.integer(pre_cfg$knn_k  %||% 10L)
  epsilon  <- pre_cfg$epsilon %||% 1e-8

  miss_stats_filt <- classify_missingness(
    filtered$mat,
    mnar_cor_threshold = -abs(mnar_thr)
  )

  result <- impute_metabolomics(
    mat        = filtered$mat,
    miss_stats = miss_stats_filt,
    knn_k      = knn_k,
    epsilon    = epsilon
  )

  message(sprintf(
    "mod_met_imputed: MNAR=%d, MAR=%d, skipped=%d, all-missing=%d.",
    result$info$n_mnar, result$info$n_mar,
    result$info$n_skipped, result$info$n_all_missing
  ))

  list(
    mat          = result$mat_imputed,
    imputed_flag = result$imputed_flag,
    meta         = filtered$meta,
    row_data     = filtered$row_data,
    info         = result$info
  )
}


# ==============================================================================
# mod_met_log — log2 transformation
# ==============================================================================

#' Apply log2 transformation to the imputed matrix
#'
#' @param imputed List returned by \code{mod_met_imputed()}.
#' @param config  Full pipeline config list.
#' @return list with: \code{mat}, \code{meta}, \code{row_data}.
#'
mod_met_log <- function(imputed, config) {
  norm_cfg    <- config$modes$metabolomics$normalization %||% list()
  method      <- norm_cfg$transform   %||% "log2"
  pseudocount <- norm_cfg$pseudocount %||% 1

  mat_log <- transform_metab(imputed$mat, method = method,
                             pseudocount = pseudocount)

  list(
    mat      = mat_log,
    meta     = imputed$meta,
    row_data = imputed$row_data
  )
}


# ==============================================================================
# mod_met_norm_comparison — compute RSD + PCA metrics, save TSV
# ==============================================================================

#' Compute normalization QC metrics for all three methods and save TSV
#'
#' @param norm_tss    List returned by \code{mod_met_normalize_linear(method="tss")}.
#' @param norm_median List returned by \code{mod_met_normalize_log()}.
#' @param norm_pqn    List returned by \code{mod_met_normalize_linear(method="pqn")}.
#' @param logged      List returned by \code{mod_met_log()}.
#' @param out_dir     Mode output directory.
#' @param config      Full pipeline config list.
#' @return Character scalar: path to saved TSV.
#'
mod_met_norm_comparison <- function(norm_tss, norm_median, norm_pqn,
                                    logged, out_dir, config) {
  methods <- list(
    list(label = "tss",      obj = norm_tss),
    list(label = "median",   obj = norm_median),
    list(label = "pqn",      obj = norm_pqn),
    list(label = "log_only", obj = list(mat = logged$mat))
  )

  rows <- lapply(methods, function(m) {
    mat <- m$obj$mat
    if (is.null(mat) || ncol(mat) < 2) {
      return(data.frame(method = m$label, median_rsd = NA_real_,
                        pc1_var = NA_real_, stringsAsFactors = FALSE))
    }

    feat_means <- rowMeans(mat, na.rm = TRUE)
    feat_sds   <- apply(mat, 1, stats::sd, na.rm = TRUE)
    rsd        <- feat_sds / abs(feat_means)
    median_rsd <- stats::median(rsd[is.finite(rsd)], na.rm = TRUE)

    pc1_var <- tryCatch({
      res <- compute_pca_scores(mat, pcs = 1L, center = TRUE, scale = FALSE)
      res$var_expl[1]
    }, error = function(e) NA_real_)

    data.frame(method = m$label, median_rsd = round(median_rsd, 4),
               pc1_var = round(pc1_var, 4), stringsAsFactors = FALSE)
  })

  comp_df  <- do.call(rbind, rows)
  ds_dir   <- file.path(out_dir, "datasets")
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(ds_dir, "normalization_comparison.tsv")
  save_tsv(comp_df, ds_dir, "normalization_comparison.tsv")

  out_file
}


# ==============================================================================
# mod_met_corrected — select norm, apply drift correction, return final matrix
# ==============================================================================

#' Select normalization method, apply optional LOESS drift correction
#'
#' @param norm_tss    List returned by \code{mod_met_normalize_linear(method="tss")}.
#' @param norm_median List returned by \code{mod_met_normalize_log()}.
#' @param norm_pqn    List returned by \code{mod_met_normalize_linear(method="pqn")}.
#' @param logged      List returned by \code{mod_met_log()}.
#' @param meta        data.frame of sample metadata.
#' @param out_dir     Mode output directory.
#' @param config      Full pipeline config list.
#' @return list with: \code{mat}, \code{meta}, \code{row_data}, \code{info}.
#'
mod_met_corrected <- function(norm_tss, norm_median, norm_pqn,
                              logged, meta, out_dir, config) {
  cfg_mode <- config$modes$metabolomics
  pre_cfg  <- cfg_mode$preprocessing %||% list()
  norm_cfg <- cfg_mode$normalization  %||% list()

  chosen_norm <- tolower(pre_cfg$chosen_norm)

  chosen_mat <- switch(chosen_norm,
    tss    = norm_tss$mat,
    median = norm_median$mat,
    pqn    = norm_pqn$mat,
    stop(sprintf("mod_met_corrected: unknown chosen_norm '%s'. ",
                 "Valid options: tss, median, pqn.", chosen_norm))
  )

  drift_result <- apply_drift_correction(chosen_mat, meta, cfg_mode)
  final_mat    <- drift_result$mat

  if (drift_result$applied) {
    diag_dir <- file.path(out_dir, "diagnostic_plots")
    dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

    dc_cfg  <- pre_cfg$drift_correction %||% list()
    inj_col <- dc_cfg$injection_order_col %||% "injection_order"

    if (inj_col %in% colnames(meta)) {
      p_drift <- tryCatch(
        plot_drift_pca_comparison(chosen_mat, final_mat, meta,
                                  injection_order_col = inj_col),
        error = function(e) {
          warning("mod_met_corrected: drift PCA plot failed: ", e$message)
          NULL
        }
      )
      if (!is.null(p_drift)) {
        ggplot2::ggsave(file.path(diag_dir, "drift_correction_pca.png"),
                        plot = p_drift, width = 12, height = 6, dpi = 150)
      }
    }
  }

  norm_info <- list(
    sample_norm   = chosen_norm,
    transform     = norm_cfg$transform %||% "log2",
    scaling       = "none",
    chosen_norm   = chosen_norm,
    drift_applied = drift_result$applied
  )

  list(
    mat      = final_mat,
    meta     = meta,
    row_data = logged$row_data,
    info     = list(
      chosen_norm   = chosen_norm,
      drift_applied = drift_result$applied,
      drift_info    = drift_result$info,
      normalization = norm_info
    )
  )
}


# ==============================================================================
# mod_met_normalize_linear — TSS / PQN on linear scale, then log2 transform
# ==============================================================================

#' Normalize a linear-scale matrix then apply log2 transformation
#'
#' @param data   List returned by \code{mod_met_imputed()} (Linear scale).
#' @param method Character: \code{"tss"} or \code{"pqn"}.
#' @param config Full pipeline config list.
#' @return list with: \code{mat} (Log2 scale), \code{meta}, \code{row_data}.
#'
mod_met_normalize_linear <- function(data, method, config) {
  method      <- tolower(method)
  norm_cfg    <- config$modes$metabolomics$normalization %||% list()
  pseudocount <- norm_cfg$pseudocount %||% 1

  mat_norm <- switch(method,
    tss = norm_total_sum(data$mat),
    pqn = norm_pqn(data$mat),
    stop(
      "mod_met_normalize_linear: 'method' must be \"tss\" or \"pqn\"; got \"",
      method, "\".  Use mod_met_normalize_log() for median normalization."
    )
  )

  mat_log <- transform_metab(mat_norm, method = "log2", pseudocount = pseudocount)

  list(
    mat      = mat_log,
    meta     = data$meta,
    row_data = data$row_data
  )
}


# ==============================================================================
# mod_met_normalize_log — Median normalization as a log-shift (subtraction)
# ==============================================================================

#' Apply median normalization as a log-shift on a Log2-scale matrix
#'
#' @param data   List returned by \code{mod_met_log()} (Log2 scale).
#' @param config Full pipeline config list.
#' @return list with: \code{mat} (Log2 scale after shift), \code{meta},
#'   \code{row_data}.
#'
mod_met_normalize_log <- function(data, config) {
  mat_log     <- data$mat
  col_medians <- apply(mat_log, 2, stats::median, na.rm = TRUE)

  valid <- is.finite(col_medians) & col_medians != 0
  if (!any(valid)) {
    warning(
      "mod_met_normalize_log: no samples have finite/non-zero log2-medians; ",
      "returning matrix unchanged."
    )
    return(list(mat = mat_log, meta = data$meta, row_data = data$row_data))
  }
  if (!all(valid)) {
    warning(sprintf(
      "mod_met_normalize_log: %d sample(s) have non-finite or zero log2-median; ",
      "they will not be shifted.", sum(!valid)
    ))
  }

  target_median      <- stats::median(col_medians[valid])
  shifts             <- col_medians - target_median
  shifts[!valid]     <- 0
  mat_shifted        <- sweep(mat_log, 2, shifts, FUN = "-")

  list(
    mat      = mat_shifted,
    meta     = data$meta,
    row_data = data$row_data
  )
}
