# R/modules/metabolomics/00_mod_preprocessing.R
#
# Target-ready wrappers for the new met_* preprocessing pipeline.
#
# Each function is called directly from the {targets} target bodies in
# pipe_metabolomics().  They orchestrate domain functions from:
#   R/domain/metabolomics/00_inputs.R         (load / parse)
#   R/domain/metabolomics/10_drift_correction.R (LOESS)
#   R/domain/metabolomics/01_normalization.R  (norm_*, transform_metab)
#   R/core/08_qc.R                            (qc_pca_scatter, norm_boxplot)
#
# Inlined helpers (previously in deleted files):
#   filter_by_missingness()  — threshold-based feature + sample filtering

#' Detect QC pool sample names from metadata
#'
#' Looks for samples where is_QC == TRUE/1/"yes", or Group/sample_type contains
#' "Pool"/"QC", or sample name starts with "Pool".
#'
#' @param meta     Metadata data.frame.
#' @param cfg_mode metabolomics mode config.
#' @return Character vector of pool sample names (may be empty).
.detect_pool_samples <- function(meta, cfg_mode) {
  sample_col <- cfg_mode$effects$samples %||% "sample_id"
  ids <- as.character(meta[[sample_col]])
  
  # Check is_QC column
  if ("is_QC" %in% colnames(meta)) {
    qc_flag <- meta$is_QC
    is_pool <- !is.na(qc_flag) & (qc_flag %in% c(TRUE, 1, "1", "yes", "Yes", "TRUE"))
    if (any(is_pool)) return(ids[is_pool])
  }
  
  # Check group/condition column for "Pool" or "QC"
  group_col <- cfg_mode$effects$color %||% cfg_mode$de$condition_column
  if (!is.null(group_col) && group_col %in% colnames(meta)) {
    grp <- as.character(meta[[group_col]])
    is_pool <- grepl("(?i)^(pool|qc)", grp)
    if (any(is_pool)) return(ids[is_pool])
  }
  
  # Fallback: sample names starting with "Pool"
  is_pool <- grepl("(?i)^pool", ids)
  if (any(is_pool)) return(ids[is_pool])
  
  character(0)
}


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
  expr_raw <- align_matrix_to_meta(expr_raw, meta, sample_col)
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

  # Populate KEGG column from HMDB → KEGG mapping so downstream modules
  # (network, etc.) can use a uniform row_data$KEGG without re-reading the file.
  row_data <- add_kegg_from_hmdb(row_data,
                                 cfg$enrichment$mapping_file %||% NULL)

  # Sanitize the KEGG column: keep only real KEGG compound IDs and route
  # ChemSpider (CSID…) ids into their own column. Runs AFTER the HMDB fill so it
  # also cleans anything the mapping introduced, and so a naive non-empty count
  # of KEGG reflects real coverage (e.g. the 438-non-empty → 315-KEGG case).
  row_data <- clean_kegg_chemspider(row_data)

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
# filter_by_missingness — low-level helper used by mod_met_filtered()
# ==============================================================================

#' Filter features and samples by missingness thresholds
#'
#' Applies two sequential filters:
#' 1. Drop samples where the fraction of missing values > \code{sample_threshold}.
#' 2. Drop features where the fraction of missing values (on remaining samples)
#'    > \code{feat_threshold}.
#'
#' Samples are filtered first so that feature missingness is evaluated on the
#' retained sample set.
#'
#' @param mat Numeric matrix (features x samples).
#' @param meta data.frame with a column matching \code{sample_col}.
#' @param sample_col Column in \code{meta} containing sample identifiers that
#'   match \code{colnames(mat)}.
#' @param feat_threshold Numeric (default 0.50 — soft filter).  Drop features
#'   missing in more than this fraction of samples.
#' @param sample_threshold Numeric (default 1.0 — sample filtering disabled).
#'   Drop samples missing more than this fraction of features.  Set to 1.0 to
#'   keep all samples.
#'
#' @return list with: \code{mat}, \code{meta}, \code{dropped_features},
#'   \code{dropped_samples}.
#'
filter_by_missingness <- function(mat, meta, sample_col,
                                  feat_threshold   = 0.50,
                                  sample_threshold = 1.0) {
  mat  <- as.matrix(mat)
  meta <- as.data.frame(meta)
  
  # 1. Filter samples first
  samp_miss_pct    <- colMeans(is.na(mat))
  keep_samps       <- colnames(mat)[samp_miss_pct <= sample_threshold]
  dropped_samples  <- setdiff(colnames(mat), keep_samps)
  
  if (length(dropped_samples) > 0) {
    message(sprintf(
      "filter_by_missingness: dropping %d sample(s) with > %.0f%% missing: %s",
      length(dropped_samples),
      sample_threshold * 100,
      paste(head(dropped_samples, 5), collapse = ", "),
      if (length(dropped_samples) > 5) "..." else ""
    ))
  }
  
  mat  <- mat[, keep_samps, drop = FALSE]
  meta <- meta[meta[[sample_col]] %in% keep_samps, , drop = FALSE]
  
  # 2. Filter features on the retained sample set
  feat_miss_pct   <- rowMeans(is.na(mat))
  keep_feats      <- rownames(mat)[feat_miss_pct <= feat_threshold]
  dropped_feats   <- setdiff(rownames(mat), keep_feats)
  
  if (length(dropped_feats) > 0) {
    message(sprintf(
      "filter_by_missingness: dropping %d feature(s) with > %.0f%% missing.",
      length(dropped_feats),
      feat_threshold * 100
    ))
  }
  
  mat <- mat[keep_feats, , drop = FALSE]
  
  message(sprintf(
    "filter_by_missingness: retained %d features x %d samples.",
    nrow(mat), ncol(mat)
  ))
  
  list(
    mat              = mat,
    meta             = meta,
    dropped_features = dropped_feats,
    dropped_samples  = dropped_samples
  )
}


# ==============================================================================
# mod_met_filtered — apply missingness thresholds
# ==============================================================================

#' Filter features and samples by missingness thresholds (target-ready wrapper)
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
    feat_threshold   = pre_cfg$feat_missing_threshold   %||% 0.50,
    sample_threshold = pre_cfg$sample_missing_threshold %||% 1.0
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
# mod_met_log — log2 transformation
# ==============================================================================

#' Apply log2 transformation to the filtered matrix
#'
#' @param filtered List returned by \code{mod_met_filtered()}.
#' @param config  Full pipeline config list.
#' @return list with: \code{mat}, \code{meta}, \code{row_data}.
#'
mod_met_log <- function(filtered, config) {
  norm_cfg    <- config$modes$metabolomics$preprocessing %||% list()
  method      <- norm_cfg$transform   %||% "log2"
  pseudocount <- norm_cfg$pseudocount %||% 1
  
  mat_log <- transform_metab(filtered$mat, method = method,
                             pseudocount = pseudocount)
  
  list(
    mat      = mat_log,
    meta     = filtered$meta,
    row_data = filtered$row_data
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
                              logged, meta, out_dir, config,
                              norm_eigenms = NULL,
                              norm_eigenms_forced = NULL,
                              norm_bio_factor = NULL) {
  cfg_mode <- config$modes$metabolomics
  pre_cfg  <- cfg_mode$preprocessing %||% list()
  norm_cfg <- cfg_mode$preprocessing  %||% list()
  
  chosen_norm <- tolower(pre_cfg$chosen_norm)
  
  chosen_mat <- switch(chosen_norm,
                       # "none": table is already normalized upstream — skip
                       # sample normalization and use the transform-only matrix.
                       # transform/scaling still apply via preprocessing.transform
                       # and preprocessing.scaling (set transform: "none" too if
                       # the table also arrives log-scaled).
                       none           = logged$mat,
                       tss            = norm_tss$mat,
                       median         = norm_median$mat,
                       pqn            = norm_pqn$mat,
                       eigenms        = if (!is.null(norm_eigenms)) norm_eigenms$mat else stop("EigenMS target not available"),
                       eigenms_forced = if (!is.null(norm_eigenms_forced)) norm_eigenms_forced$mat else stop("EigenMS_forced target not available"),
                       bio_factor     = if (!is.null(norm_bio_factor)) norm_bio_factor$mat else stop("mod_met_corrected: chosen_norm = 'bio_factor' but the normalization returned NULL. Set preprocessing.biological_factor_col to a per-sample metadata column (e.g. total protein)."),
                       stop(sprintf("mod_met_corrected: unknown chosen_norm '%s'. ",
                                    "Valid options: none, tss, median, pqn, eigenms, eigenms_forced, bio_factor.", chosen_norm))
  )
  
  # Apply scaling if configured
  scaling_method <- tolower(norm_cfg$scaling %||% "none")
  if (scaling_method != "none") {
    message(sprintf("mod_met_corrected: applying '%s' scaling", scaling_method))
    chosen_mat <- scale_metab(chosen_mat, method = scaling_method)
  }
  
  # Phase 2 of sample filter: remove QC/pool samples AFTER normalization
  # (only when exclude_after_norm = true)
  if (isTRUE(cfg_mode$sample_filter$exclude_after_norm)) {
    rules <- get_sample_filter_rules_metab(cfg_mode)
    if (!is.null(rules)) {
      sample_col <- cfg_mode$effects$samples %||% "sample_id"
      all_ids <- colnames(chosen_mat)
      keep_ids <- apply_sample_filter_metab(all_ids, meta, rules, sample_col)
      if (length(keep_ids) < length(all_ids)) {
        removed <- setdiff(all_ids, keep_ids)
        message(sprintf(
          "mod_met_corrected: post-normalization filter removed %d sample(s): %s",
          length(removed), paste(removed, collapse = ", ")
        ))
        chosen_mat <- chosen_mat[, keep_ids, drop = FALSE]
        meta <- meta[meta[[sample_col]] %in% keep_ids, , drop = FALSE]
      }
    }
  }
  
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
    scaling       = scaling_method,
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
#' @param data   List returned by \code{mod_met_filtered()} (Linear scale).
#' @param method Character: \code{"tss"} or \code{"pqn"}.
#' @param config Full pipeline config list.
#' @return list with: \code{mat} (Log2 scale), \code{meta}, \code{row_data}.
#'
mod_met_normalize_linear <- function(data, method, config) {
  method      <- tolower(method)
  norm_cfg    <- config$modes$metabolomics$preprocessing %||% list()
  pseudocount <- norm_cfg$pseudocount %||% 1
  
  pre_cfg <- config$modes$metabolomics$preprocessing %||% list()
  dc_cfg  <- pre_cfg$drift_correction %||% list()
  qc_col  <- dc_cfg$qc_flag_col %||% "is_QC"
  
  qc_idx <- NULL
  
  if (method == "pqn") {
    meta <- data$meta
    
    if (is.null(meta) || !qc_col %in% colnames(meta)) {
      warning(sprintf(
        "mod_met_normalize_linear: QC flag column '%s' not found; PQN will use all samples as reference.",
        qc_col
      ))
    } else {
      qc_flag <- isTRUE_vec(meta[[qc_col]])
      idx <- which(qc_flag)
      
      if (length(idx) >= 2L) {
        qc_idx <- idx
      } else {
        warning(sprintf(
          "mod_met_normalize_linear: QC flag column '%s' has fewer than 2 QC samples; PQN will use all samples as reference.",
          qc_col
        ))
      }
    }
  }
  
  # PQN reference samples: use preprocessing.pqn_reference to specify
  # which sample(s) to use as the reference spectrum (e.g. a QC pool).
  # Options: a sample name, "pools" (auto-detect from is_QC/Pool columns),
  #          "median_pool" (middle pool by injection order), or null (default: all).
  pqn_ref <- norm_cfg$pqn_reference %||% NULL
  ref_samples <- NULL
  if (method == "pqn" && !is.null(pqn_ref)) {
    if (tolower(pqn_ref) == "pools") {
      # Auto-detect pool samples from metadata
      ref_samples <- .detect_pool_samples(data$meta, config$modes$metabolomics)
      if (length(ref_samples) > 0) {
        message(sprintf("norm_pqn: using %d pool sample(s) as reference: %s",
                        length(ref_samples), paste(ref_samples, collapse = ", ")))
      }
    } else if (tolower(pqn_ref) == "median_pool") {
      pools <- .detect_pool_samples(data$meta, config$modes$metabolomics)
      if (length(pools) > 0) {
        # Pick the middle pool (by position in metadata / injection order)
        mid_idx <- ceiling(length(pools) / 2)
        ref_samples <- pools[mid_idx]
        message(sprintf("norm_pqn: using median pool '%s' as single reference", ref_samples))
      }
    } else {
      # Treat as explicit sample name(s)
      ref_samples <- unlist(strsplit(as.character(pqn_ref), ",\\s*"))
    }
  }
  
  mat_norm <- switch(method,
                     tss = norm_total_sum(data$mat),
                     pqn = norm_pqn(data$mat, qc_idx = qc_idx),
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
# mod_met_normalize_eigenms — EigenMS on linear scale, then log2 transform
# ==============================================================================

#' Apply EigenMS normalization then log2 transformation
#'
#' EigenMS uses SVD on residuals from a group-aware model to identify and
#' remove systematic technical bias while preserving biological signal.
#'
#' @param data   List returned by \code{mod_met_imputed()} (Linear scale).
#' @param config Full pipeline config list.
#' @return list with: \code{mat} (Log2 scale), \code{meta}, \code{row_data}.
mod_met_normalize_eigenms <- function(data, config) {
  norm_cfg    <- config$modes$metabolomics$preprocessing %||% list()
  pseudocount <- norm_cfg$pseudocount %||% 1
  
  # Extract group labels for EigenMS
  cfg_mode  <- config$modes$metabolomics
  group_col <- cfg_mode$effects$color %||% cfg_mode$de$condition_column %||% "sample_type"
  groups <- if (!is.null(data$meta) && group_col %in% colnames(data$meta)) {
    as.character(data$meta[[group_col]])
  } else {
    NULL
  }
  
  seed     <- config$params$seed %||% 1
  mat_norm <- norm_eigenms(data$mat, groups = groups, seed = seed)
  eigenms_info <- attr(mat_norm, "eigenms_info")
  mat_log  <- transform_metab(mat_norm, method = "log2", pseudocount = pseudocount)
  
  list(
    mat          = mat_log,
    meta         = data$meta,
    row_data     = data$row_data,
    eigenms_info = eigenms_info
  )
}


# ==============================================================================
# mod_met_normalize_eigenms_forced — Forced EigenMS (NOREVA-style)
# ==============================================================================

#' Apply forced EigenMS normalization then log2 transformation
#'
#' Forces removal of eigentrends without statistical significance testing.
#' This replicates the NOREVA/Ifat Abramovich approach.
#'
#' @param data   List returned by \code{mod_met_imputed()} (Linear scale).
#' @param config Full pipeline config list.
#' @return list with: \code{mat} (Log2 scale), \code{meta}, \code{row_data}.
mod_met_normalize_eigenms_forced <- function(data, config) {
  norm_cfg    <- config$modes$metabolomics$preprocessing %||% list()
  pseudocount <- norm_cfg$pseudocount %||% 1
  
  cfg_mode  <- config$modes$metabolomics
  group_col <- cfg_mode$effects$color %||% cfg_mode$de$condition_column %||% "sample_type"
  groups <- if (!is.null(data$meta) && group_col %in% colnames(data$meta)) {
    as.character(data$meta[[group_col]])
  } else {
    NULL
  }
  
  # Forced EigenMS operates on log2 scale (like Ifat's pipeline)
  pre_cfg <- config$modes$metabolomics$preprocessing %||% list()
  n_forced <- pre_cfg$n_eigentrends_forced
  # If user specified an explicit count, convert to fraction of samples
  pct <- if (!is.null(n_forced) && is.numeric(n_forced)) {
    n_forced / ncol(data$mat)
  } else {
    0.2  # default: 20% of samples
  }
  mat_log  <- transform_metab(data$mat, method = "log2", pseudocount = pseudocount)
  mat_norm <- norm_eigenms_forced(mat_log, groups = groups, pct_eigentrends = pct)
  eigenms_info <- attr(mat_norm, "eigenms_info")
  
  list(
    mat          = mat_norm,
    meta         = data$meta,
    row_data     = data$row_data,
    eigenms_info = eigenms_info
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


# ==============================================================================
# mod_met_normalize_bio_factor — divide each sample by a measured biological
# covariate (e.g. total protein), on the linear scale, then log2 transform
# ==============================================================================

#' Apply biological-factor normalization then log2 transformation
#'
#' Divides each sample by a per-sample numeric value (e.g. mg of total protein
#' from a NanoDrop / BCA assay) read from a metadata column, then log2s — the
#' same linear-then-log2 shape as \code{mod_met_normalize_linear()}. Use when
#' intensities should be scaled to a measured input amount rather than to a
#' within-sample statistic. The column is named by
#' \code{preprocessing.biological_factor_col}; the values live in the sample
#' metadata. Returns \code{NULL} when that column is unset (the target is built
#' on every run but only used when \code{chosen_norm = "bio_factor"}).
#'
#' @param data   List returned by \code{mod_met_imputed()} (Linear scale); must
#'   carry a \code{meta} table containing \code{biological_factor_col}.
#' @param config Full pipeline config list.
#' @return list with: \code{mat} (Log2 scale), \code{meta}, \code{row_data}, or
#'   \code{NULL} when \code{biological_factor_col} is unconfigured.
mod_met_normalize_bio_factor <- function(data, config) {
  cfg_mode    <- config$modes$metabolomics
  norm_cfg    <- cfg_mode$preprocessing %||% list()
  pseudocount <- norm_cfg$pseudocount %||% 1
  factor_col  <- norm_cfg$biological_factor_col

  # This target is built on every run (targets has no conditional build) but is
  # only consumed when chosen_norm == "bio_factor". Return NULL when the covariate
  # column is unconfigured so non-bio_factor runs don't fail here; the guard in
  # mod_met_corrected() raises a clear error only if bio_factor is actually chosen.
  if (is.null(factor_col) || !nzchar(factor_col)) {
    message("mod_met_normalize_bio_factor: preprocessing.biological_factor_col is not ",
            "set; returning NULL (only an error if chosen_norm = 'bio_factor').")
    return(NULL)
  }

  # Align the covariate to sample columns via the configured sample-id column.
  sample_col <- cfg_mode$effects$samples %||% NULL

  mat_norm <- normalize_samples(data$mat, method = "bio_factor",
                                meta = data$meta, bio_factor_col = factor_col,
                                sample_col = sample_col)
  mat_log  <- transform_metab(mat_norm, method = "log2", pseudocount = pseudocount)

  list(
    mat      = mat_log,
    meta     = data$meta,
    row_data = data$row_data
  )
}
