# R/modules/metabolomics/02_mod_qc_suite.R
#
# QC suite module for the metabolomics preprocessing pipeline.
#
# Provides:
#   compute_linear_rsd()           — authoritative RSD (back-transform + SD/|mean|)
#   qc_full_metabolomics_suite()   — all plots + metrics for one stage × one subset
#   mod_met_qc_suite()             — target-facing wrapper (with_qc + no_qc)
#   mod_met_qc_comparison_table()  — aggregate benchmark TSV across all stages
#
# Scale policy: every matrix passed to qc_full_metabolomics_suite() is on the
# Log2 scale.  No internal transformation is applied.  All QC visualisations
# (PCA, density, boxplot, heatmaps) consume Log2 values directly; Euclidean
# distance and Pearson correlation are therefore computed in log2 space, which
# is correct and standard for omics QC.


# ==============================================================================
# compute_linear_rsd — scale-aware RSD on back-transformed linear intensities
# ==============================================================================

#' Compute median RSD after back-transforming Log2 values to linear scale
#'
#' This function is the \strong{authoritative source} for Relative Standard
#' Deviation (RSD) computation in the metabolomics QC pipeline.  All RSD
#' values reported in QC metrics TSVs and the normalization comparison table
#' must be produced by this function.
#'
#' \strong{Why back-transformation is required:}
#' RSD is defined as \eqn{\mathrm{SD}(x) / |\mathrm{mean}(x)|} for linear
#' intensities.  Computing this ratio directly on log2 values yields the
#' coefficient of variation of the log-transformed data, which has no standard
#' interpretive meaning in metabolomics QC and must not be used.
#'
#' \strong{Scale Contract alignment:}
#' Every Log2-output target in the pipeline (\code{met_log}, \code{met_norm_tss},
#' \code{met_norm_pqn}, \code{met_norm_median}) used a transformation of the
#' form \eqn{\log_2(x + p)} where \eqn{p} is the pseudocount.  The exact
#' inverse is \eqn{x = 2^{\text{val}} - p}.  For \code{met_norm_median} this
#' inversion is approximate (see \code{backtransform_exact} in the return
#' value); for all other stages it is exact.
#'
#' @param mat         Numeric matrix (features \eqn{\times} samples) on the
#'   \strong{Log2} scale.  Must have at least 1 row and 2 columns for
#'   meaningful RSD computation.
#' @param stage       Character scalar identifying the preprocessing stage.
#'   Must be one of \code{"log"}, \code{"norm_tss"}, \code{"norm_pqn"}, or
#'   \code{"norm_median"}.  Used in warning messages and to set the
#'   \code{backtransform_exact} flag.
#' @param pseudocount Numeric scalar \eqn{\geq 0}.  The pseudocount added
#'   before the log2 transform in \code{transform_metab()} (default 1).
#'   \strong{Must be the same value used when the matrix was created.}
#'   Do not infer this from the data.
#'
#' @return Named list with:
#'   \describe{
#'     \item{\code{median_rsd}}{Numeric scalar.  Median RSD across all
#'       features (NA if no finite RSD values could be computed).}
#'     \item{\code{rsd_per_feature}}{Named numeric vector of per-feature RSD
#'       values (NA for features with zero or non-finite mean after
#'       back-transformation).}
#'     \item{\code{n_features}}{Integer.  Number of rows in \code{mat}.}
#'     \item{\code{n_samples}}{Integer.  Number of columns in \code{mat}.}
#'     \item{\code{backtransform_exact}}{Logical.  \code{TRUE} for stages
#'       \code{"log"}, \code{"norm_tss"}, \code{"norm_pqn"} where
#'       \eqn{2^{\text{val}} - p} exactly reverses the log2 transform.
#'       \code{FALSE} for \code{"norm_median"} where the log-shift
#'       operation does not preserve the clean \eqn{\log_2(x + p)}
#'       relationship, making the inversion approximate.}
#'     \item{\code{pseudocount_used}}{Numeric.  The pseudocount value
#'       supplied; echoed back for audit purposes.}
#'   }
#'
#' @details
#' \strong{Back-transformation:}
#' \deqn{\text{linear\_val}_{ij} = 2^{\text{mat}_{ij}} - p}
#' where \eqn{p} is the pseudocount.
#'
#' \strong{Negative-value guard:}
#' Theoretically, \eqn{2^{\text{val}} \geq p} for any value produced by
#' \eqn{\log_2(x + p)} with \eqn{x \geq 0}.  In practice, numerical
#' imprecision or approximate back-transformation (median stage) can yield
#' small negative values.  Any cell where \eqn{2^{\text{val}} - p < 0} is
#' floored to \eqn{\varepsilon = 2.2 \times 10^{-16}}
#' (\code{.Machine$double.eps}) and a warning is issued.
#'
#' \strong{Per-feature RSD:}
#' \deqn{\mathrm{RSD}_i =
#'   \frac{\mathrm{SD}_j(\text{linear\_val}_{ij})}{|\mathrm{mean}_j(\text{linear\_val}_{ij})|}
#' }
#' Features where the absolute mean is 0 or non-finite receive \code{NA}.
#'
#' \strong{Median RSD:}
#' \deqn{\mathrm{median\_rsd} = \mathrm{median}_i(\mathrm{RSD}_i)}
#' computed with \code{na.rm = TRUE}.
compute_linear_rsd <- function(mat, stage, pseudocount = 1) {
  # ---- Input validation -------------------------------------------------------
  valid_stages <- c("log", "norm_tss", "norm_pqn", "norm_median")
  if (!stage %in% valid_stages) {
    warning(sprintf(
      "compute_linear_rsd: unknown stage '%s'; expected one of: %s.",
      stage, paste(valid_stages, collapse = ", ")
    ))
  }
  if (!is.numeric(pseudocount) || length(pseudocount) != 1 || pseudocount < 0) {
    stop("compute_linear_rsd: 'pseudocount' must be a single non-negative number.")
  }
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("compute_linear_rsd: 'mat' must be a numeric matrix.")
  }
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    warning(sprintf(
      "compute_linear_rsd [stage='%s']: matrix has 0 rows or 0 columns; returning NA.", stage
    ))
    return(list(
      median_rsd          = NA_real_,
      rsd_per_feature     = setNames(numeric(0), character(0)),
      n_features          = nrow(mat),
      n_samples           = ncol(mat),
      backtransform_exact = stage %in% c("log", "norm_tss", "norm_pqn"),
      pseudocount_used    = pseudocount
    ))
  }

  # ---- Back-transformation: linear_val = 2^(log2_val) - pseudocount ----------
  #
  # Exactly inverts transform_metab(mat, method = "log2", pseudocount = p),
  # which computes log2(x + p).  The inverse is 2^(val) - p = x (linear).
  # For met_norm_median, the log-shift does not preserve the log2(x + p) form
  # exactly, so this inversion is approximate for that stage.
  mat_linear <- 2^mat - pseudocount

  # ---- Negative-value guard --------------------------------------------------
  epsilon <- .Machine$double.eps
  n_neg   <- sum(mat_linear < 0, na.rm = TRUE)
  if (n_neg > 0) {
    warning(sprintf(
      "compute_linear_rsd [stage='%s']: %d cell(s) produced negative values after back-transform (2^val - %g); floored to %g.  If stage is 'norm_median', this is expected (approximate back-transform).",
      stage, n_neg, pseudocount, epsilon
    ))
    mat_linear[mat_linear < 0] <- epsilon
  }

  # ---- Per-feature RSD: SD / |mean| across samples ---------------------------
  feat_sds      <- apply(mat_linear, 1, stats::sd, na.rm = TRUE)
  feat_means    <- rowMeans(mat_linear, na.rm = TRUE)
  feat_means_ab <- abs(feat_means)

  # Guard: zero or non-finite absolute mean → RSD undefined → NA
  feat_means_ab[feat_means_ab == 0 | !is.finite(feat_means_ab)] <- NA_real_

  rsd <- feat_sds / feat_means_ab
  rsd[!is.finite(rsd)] <- NA_real_

  if (!is.null(rownames(mat))) names(rsd) <- rownames(mat)

  list(
    median_rsd          = stats::median(rsd, na.rm = TRUE),
    rsd_per_feature     = rsd,
    n_features          = nrow(mat_linear),
    n_samples           = ncol(mat_linear),
    backtransform_exact = stage %in% c("log", "norm_tss", "norm_pqn"),
    pseudocount_used    = pseudocount
  )
}


# ==============================================================================
# .qc_write_metrics — private helper: write metrics_summary.tsv
# ==============================================================================

# Not exported.  Called at the end of every qc_full_metabolomics_suite() call,
# including skip and partial cases, so {targets} always has a valid file to hash.
.qc_write_metrics <- function(out_dir, stage, subset_mode, mat, meta,
                               cfg_mode, pseudocount, skipped, status,
                               pca_var) {
  color_col  <- cfg_mode$effects$color   %||% "sample_id"
  sample_col <- cfg_mode$effects$samples %||% "sample_id"

  # ---- Overall RSD -----------------------------------------------------------
  rsd_result <- tryCatch(
    compute_linear_rsd(mat, stage = stage, pseudocount = pseudocount),
    error = function(e) {
      warning(sprintf(".qc_write_metrics: RSD failed for [%s][%s]: %s",
                      stage, subset_mode, e$message))
      NULL
    }
  )

  # ---- Per-group RSD ---------------------------------------------------------
  group_rsd_str <- NA_character_
  if (!is.null(rsd_result) && color_col %in% colnames(meta)) {
    groups <- unique(meta[[color_col]])
    group_vals <- vapply(groups, function(g) {
      gids <- meta[[sample_col]][meta[[color_col]] == g]
      cols <- intersect(as.character(gids), colnames(mat))
      if (length(cols) < 2L) return(NA_real_)
      r <- tryCatch(
        compute_linear_rsd(mat[, cols, drop = FALSE], stage, pseudocount),
        error = function(e) NULL
      )
      if (is.null(r)) NA_real_ else r$median_rsd
    }, numeric(1L))
    valid_pairs <- !is.na(group_vals)
    if (any(valid_pairs)) {
      group_rsd_str <- paste(
        paste0(groups[valid_pairs], "=", round(group_vals[valid_pairs], 6L)),
        collapse = ";"
      )
    }
  }

  # ---- PCA variance ----------------------------------------------------------
  pc1 <- if (!is.null(pca_var) && length(pca_var) >= 1L) pca_var[[1L]] else NA_real_
  pc2 <- if (!is.null(pca_var) && length(pca_var) >= 2L) pca_var[[2L]] else NA_real_
  pc3 <- if (!is.null(pca_var) && length(pca_var) >= 3L) pca_var[[3L]] else NA_real_

  # ---- Assemble TSV row ------------------------------------------------------
  tsv_row <- data.frame(
    stage                   = stage,
    subset_mode             = subset_mode,
    n_samples               = ncol(mat),
    n_features              = nrow(mat),
    pca_pc1_var_expl        = pc1,
    pca_pc2_var_expl        = pc2,
    pca_pc3_var_expl        = pc3,
    median_rsd_overall      = if (!is.null(rsd_result)) rsd_result$median_rsd else NA_real_,
    median_rsd_per_group    = group_rsd_str,
    rsd_backtransform_exact = if (!is.null(rsd_result)) rsd_result$backtransform_exact else NA,
    pseudocount_used        = pseudocount,
    plots_skipped           = if (length(skipped) > 0L) paste(skipped, collapse = "; ") else "",
    status                  = status,
    stringsAsFactors        = FALSE
  )

  tsv_path <- file.path(out_dir, "metrics_summary.tsv")
  utils::write.table(tsv_row, tsv_path, sep = "\t", row.names = FALSE, quote = FALSE)
  tsv_path
}


# ==============================================================================
# qc_full_metabolomics_suite — all plots + metrics for one stage × one subset
# ==============================================================================

#' Generate a full QC visualisation suite for one preprocessing stage and one
#' sample subset
#'
#' This is the core plotting function of the QC layer.  It is called twice by
#' \code{mod_met_qc_suite()}: once with \code{subset_mode = "with_qc"} (all
#' samples) and once with \code{subset_mode = "no_qc"} (biological samples
#' only, QC-flagged samples removed).
#'
#' \strong{Scale policy:} \code{mat} must be on the \strong{Log2} scale.  All
#' visualisations (PCA, density, boxplot, heatmaps) consume Log2 values
#' directly.  Euclidean distance and Pearson correlation are computed in log2
#' space, which is correct and standard for omics QC diagnostics.
#'
#' @param mat         Numeric matrix (features \eqn{\times} samples) on the
#'   \strong{Log2} scale.
#' @param meta        Data frame, one row per sample.  Must be aligned to
#'   \code{colnames(mat)} (same order and set, as guaranteed by upstream
#'   module functions).
#' @param stage       Character scalar matching a row in the Scale Contract
#'   Table.  One of \code{"log"}, \code{"norm_tss"}, \code{"norm_median"},
#'   \code{"norm_pqn"}.  Used for plot titles, file-name prefixes, and the
#'   \code{backtransform_exact} flag in metrics.
#' @param pseudocount Numeric scalar.  The pseudocount used when log2-
#'   transforming the matrix (read from
#'   \code{config$modes$metabolomics$normalization$pseudocount}).  Passed
#'   unchanged to \code{compute_linear_rsd()} for exact back-transformation.
#' @param cfg         The QC config section:
#'   \code{config$modes$metabolomics$qc}.  Controls toggles
#'   (\code{cfg$plots$*}), thresholds (\code{cfg$thresholds$*}), and
#'   \code{cfg$qc_flag_column}.
#' @param cfg_mode    The metabolomics mode config:
#'   \code{config$modes$metabolomics}.  Passed to the underlying
#'   \code{qc_*} functions from \code{R/core/08_qc.R} which require
#'   \code{cfg$effects$samples}, \code{cfg$effects$color}, etc.
#' @param out_dir     Absolute path to the output directory for this
#'   stage × subset combination (e.g.
#'   \code{{metab_out_dir}/qc/norm_tss/with_qc/}).  Must already exist.
#' @param subset_mode Character scalar, \code{"with_qc"} or \code{"no_qc"}.
#'
#' @return Character vector of absolute paths to every file written (PNGs +
#'   \code{metrics_summary.tsv}).  Always includes at least the TSV.
#'   The returned vector is suitable for use as the value of a
#'   \code{format = "file"} target in \code{{targets}}.
#'
#' @details
#' \strong{Subsetting (\code{subset_mode = "no_qc"}):}
#' Reads the column named \code{cfg$qc_flag_column} (default \code{"is_QC"})
#' from \code{meta}.  Samples where this column equals \code{TRUE} (logical),
#' \code{"TRUE"}, \code{"true"}, or \code{"1"} (character) are removed.
#' If the column is absent, or if no QC samples are found, the function writes
#' a metrics TSV with the appropriate \code{status} field and returns without
#' generating any plots.
#'
#' \strong{Graceful skip semantics:}
#' Each plot call is individually wrapped in \code{tryCatch()}.  A failure in
#' one plot never prevents subsequent plots from being attempted.  All skipped
#' plots are logged via \code{message()} and recorded in the
#' \code{plots_skipped} column of \code{metrics_summary.tsv}.
#'
#' \strong{Status values in metrics TSV:}
#' \describe{
#'   \item{\code{"ok"}}{All requested plots were generated.}
#'   \item{\code{"partial"}}{Some plots were skipped (threshold or toggle).}
#'   \item{\code{"skipped_no_flag_column"}}{no_qc requested but flag absent.}
#'   \item{\code{"skipped_no_qc_samples"}}{no_qc requested, flag present, no QC samples found.}
#' }
qc_full_metabolomics_suite <- function(mat, meta, stage, pseudocount,
                                       cfg, cfg_mode, out_dir, subset_mode) {
  # ---- Read toggles and thresholds (with defaults) ---------------------------
  qc_flag_col  <- cfg$qc_flag_column %||% "is_QC"
  max_heatmap  <- cfg$thresholds$max_samples_for_heatmaps %||% 120L
  min_pca      <- cfg$thresholds$min_samples_for_pca      %||% 3L
  top_features <- cfg$thresholds$top_variable_features    %||% 100L

  do_pca   <- isTRUE(cfg$plots$pca                  %||% TRUE)
  do_dens  <- isTRUE(cfg$plots$density               %||% TRUE)
  do_box   <- isTRUE(cfg$plots$boxplot               %||% TRUE)
  do_hist  <- isTRUE(cfg$plots$histogram             %||% TRUE)
  do_rle   <- isTRUE(cfg$plots$rle                    %||% TRUE)
  do_hcor  <- isTRUE(cfg$plots$heatmap_correlation   %||% TRUE)
  do_hdist <- isTRUE(cfg$plots$heatmap_distance      %||% TRUE)
  do_hexp  <- isTRUE(cfg$plots$heatmap_expression    %||% TRUE)

  files   <- character(0L)
  skipped <- character(0L)
  status  <- "ok"
  pca_var <- NULL          # full var_expl vector; populated after first PCA call

  mat_sub  <- as.matrix(mat)
  meta_sub <- as.data.frame(meta)

  sample_col <- cfg_mode$effects$samples %||% "sample_id"
  color_col  <- cfg_mode$effects$color   %||% "sample_id"
  shape_col  <- cfg_mode$effects$shape

  # ---- Subsetting ------------------------------------------------------------
  if (subset_mode == "no_qc") {

    if (!qc_flag_col %in% colnames(meta_sub)) {
      message(sprintf(
        "[QC][%s][no_qc] Flag column '%s' absent from metadata; skipping plots.",
        stage, qc_flag_col
      ))
      tsv <- .qc_write_metrics(out_dir, stage, subset_mode, mat_sub, meta_sub,
                                cfg_mode, pseudocount, skipped,
                                status = "skipped_no_flag_column", pca_var = NULL)
      return(c(files, tsv))
    }

    flag_vals     <- meta_sub[[qc_flag_col]]
    qc_flag_value <- cfg$qc_flag_value %||% "qc"

    if (is.logical(flag_vals)) {
      # Boolean column (e.g. is_QC = TRUE / FALSE).
      qc_mask <- !is.na(flag_vals) & flag_vals
    } else {
      char_vals <- tolower(trimws(as.character(flag_vals)))
      char_vals[is.na(flag_vals)] <- NA_character_

      # Two acceptance criteria, OR-combined:
      #   1. Backward-compat boolean-string columns: "TRUE" / "1"
      #      (e.g. a metadata column exported as character from Excel)
      #   2. Value-match: case-insensitive compare to qc_flag_value
      #      (e.g. treatment == "qc", treatment == "QC", treatment == "Qc")
      bool_match  <- char_vals %in% c("true", "1")
      value_match <- !is.na(char_vals) &
                     char_vals == tolower(trimws(qc_flag_value))
      qc_mask <- bool_match | value_match
    }
    qc_mask[is.na(qc_mask)] <- FALSE

    if (!any(qc_mask)) {
      message(sprintf(
        "[QC][%s][no_qc] No QC samples found in '%s'; skipping plots.",
        stage, qc_flag_col
      ))
      tsv <- .qc_write_metrics(out_dir, stage, subset_mode, mat_sub, meta_sub,
                                cfg_mode, pseudocount, skipped,
                                status = "skipped_no_qc_samples", pca_var = NULL)
      return(c(files, tsv))
    }

    qc_ids   <- as.character(meta_sub[[sample_col]][qc_mask])
    keep_ids <- setdiff(colnames(mat_sub), qc_ids)
    mat_sub  <- mat_sub[, keep_ids, drop = FALSE]
    meta_sub <- meta_sub[meta_sub[[sample_col]] %in% keep_ids, , drop = FALSE]
    rownames(meta_sub) <- NULL
  }

  n_samp <- ncol(mat_sub)
  n_feat <- nrow(mat_sub)

  # ---- Fatal validation (any subset_mode) ------------------------------------
  if (n_samp == 0L) {
    stop(sprintf(
      "qc_full_metabolomics_suite [%s][%s]: 0 samples after subsetting.",
      stage, subset_mode
    ))
  }
  if (n_feat == 0L) {
    stop(sprintf(
      "qc_full_metabolomics_suite [%s][%s]: 0 features in matrix.",
      stage, subset_mode
    ))
  }

  # ---- PCA -------------------------------------------------------------------
  # n_pcs: maximum computable PCs = min(n_samples - 1, n_features).
  # We cap at 3 since we only need PC1/PC2/PC3 for plots and metrics.
  n_pcs_possible <- min(n_samp - 1L, n_feat, 3L)

  if (do_pca && n_samp >= min_pca && n_pcs_possible >= 2L) {

    # Compute PCA once for metrics; var_expl covers all n_pcs_possible PCs.
    pca_ok <- tryCatch({
      pca_res <- compute_pca_scores(mat_sub, pcs = seq_len(n_pcs_possible))
      pca_var <- pca_res$var_expl   # full variance-explained vector
      TRUE
    }, error = function(e) {
      skipped  <<- c(skipped, sprintf("pca (compute): %s", conditionMessage(e)))
      status   <<- "partial"
      message(sprintf("[QC][%s][%s] PCA computation failed: %s",
                      stage, subset_mode, conditionMessage(e)))
      FALSE
    })
    
    if (pca_ok) {
      # Save scores + color column as TSV so the HTML report can render
      # interactive PCA without needing the original matrices.
      tryCatch({
        scores_df <- pca_res$scores   # has PC1..PCn + sample columns
        idx <- match(scores_df$sample, as.character(meta_sub[[sample_col]]))
        if (color_col %in% colnames(meta_sub))
          scores_df[[color_col]] <- meta_sub[[color_col]][idx]
        if (!is.null(shape_col) && shape_col %in% colnames(meta_sub))
          scores_df[[shape_col]] <- meta_sub[[shape_col]][idx]
        # Embed variance explained as extra columns so the report has axis labels
        for (k in seq_along(pca_var)) {
          scores_df[[paste0("var_PC", k)]] <- pca_var[[k]]
        }
        scores_file <- file.path(out_dir, "pca_scores.tsv")
        utils::write.table(scores_df, scores_file,
                           sep = "\t", row.names = FALSE, quote = FALSE)
        files <- c(files, scores_file)
      }, error = function(e) {
        message(sprintf("[QC][%s][%s] pca_scores.tsv skipped: %s",
                        stage, subset_mode, conditionMessage(e)))
      })

      # PC1 vs PC2
      pca12_file <- file.path(out_dir, "pca_pc1_pc2.png")
      pca12_success <- tryCatch({
        qc_pca_scatter(mat_sub, meta_sub, cfg_mode, pcs = c(1L, 2L), out_file = pca12_file)
        TRUE
        }, error = function(e) {
        skipped <<- c(skipped, sprintf("pca_pc1_pc2: %s", conditionMessage(e)))
        status  <<- "partial"
        message(sprintf("[QC][%s][%s] Skipped pca_pc1_pc2: %s",
                        stage, subset_mode, conditionMessage(e)))
        FALSE
      })
      if (isTRUE(pca12_success)) files <- c(files, pca12_file)

      # PC1 vs PC3 (only when at least 3 PCs are computable)
      if (n_pcs_possible >= 3L) {
        pca13_file <- file.path(out_dir, "pca_pc1_pc3.png")
        pca13_success <- tryCatch({
          qc_pca_scatter(mat_sub, meta_sub, cfg_mode, pcs = c(1L, 3L), out_file = pca13_file)
          TRUE
        }, error = function(e) {
          skipped <<- c(skipped, sprintf("pca_pc1_pc3: %s", conditionMessage(e)))
          status  <<- "partial"
          message(sprintf("[QC][%s][%s] Skipped pca_pc1_pc3: %s",
                          stage, subset_mode, conditionMessage(e)))
          FALSE
        })
        
        if (pca13_success) files <- c(files, pca13_file)
      } else {
        skipped <- c(skipped,
          sprintf("pca_pc1_pc3: only %d PC(s) computable (need 3)", n_pcs_possible))
        status <- "partial"
      }
    }

  } else if (do_pca) {
    reason <- if (n_samp < min_pca) {
      sprintf("n_samples (%d) < min_samples_for_pca (%d)", n_samp, min_pca)
    } else {
      sprintf("only %d PC(s) computable (need 2)", n_pcs_possible)
    }
    skipped <- c(skipped, sprintf("pca: %s", reason))
    status  <- "partial"
    message(sprintf("[QC][%s][%s] PCA skipped: %s", stage, subset_mode, reason))
  }

  # ---- Density ---------------------------------------------------------------
  if (do_dens && n_samp >= 2L) {
    dens_file <- file.path(out_dir, "density.png")
    tryCatch({
      qc_omic_density(
        mat_sub, meta_sub, cfg_mode,
        out_file = dens_file,
        title    = sprintf("Intensity density — %s [%s]", stage, subset_mode)
      )
      files <- c(files, dens_file)
    }, error = function(e) {
      skipped <<- c(skipped, sprintf("density: %s", conditionMessage(e)))
      status  <<- "partial"
      message(sprintf("[QC][%s][%s] Skipped density: %s",
                      stage, subset_mode, conditionMessage(e)))
    })
  } else if (do_dens) {
    skipped <- c(skipped, sprintf("density: n_samples (%d) < 2", n_samp))
    status  <- "partial"
  }

  # ---- Boxplot ---------------------------------------------------------------
  if (do_box && n_samp >= 2L) {
    box_file <- file.path(out_dir, "boxplot.png")
    tryCatch({
      norm_boxplot(
        mat_sub, meta_sub, cfg_mode,
        out_file = box_file,
        title    = sprintf("Sample intensities — %s [%s]", stage, subset_mode)
      )
      files <- c(files, box_file)
    }, error = function(e) {
      skipped <<- c(skipped, sprintf("boxplot: %s", conditionMessage(e)))
      status  <<- "partial"
      message(sprintf("[QC][%s][%s] Skipped boxplot: %s",
                      stage, subset_mode, conditionMessage(e)))
    })
  } else if (do_box) {
    skipped <- c(skipped, sprintf("boxplot: n_samples (%d) < 2", n_samp))
    status  <- "partial"
  }

  # ---- RLE (Relative Log Expression) -----------------------------------------
  if (do_rle && n_samp >= 2L) {
    rle_file <- file.path(out_dir, "rle_boxplot.png")
    tryCatch({
      qc_rle_boxplot(
        mat_sub, meta_sub, cfg_mode,
        out_file = rle_file
      )
      files <- c(files, rle_file)
    }, error = function(e) {
      skipped <<- c(skipped, sprintf("rle: %s", conditionMessage(e)))
      status  <<- "partial"
      message(sprintf("[QC][%s][%s] Skipped RLE: %s",
                      stage, subset_mode, conditionMessage(e)))
    })
  } else if (do_rle) {
    skipped <- c(skipped, sprintf("rle: n_samples (%d) < 2", n_samp))
    status  <- "partial"
  }

  # ---- Histogram (2–8 groups) ------------------------------------------------
  if (do_hist) {
    n_groups <- if (color_col %in% colnames(meta_sub)) {
      length(unique(meta_sub[[color_col]]))
    } else { 0L }

    if (n_groups >= 2L && n_groups <= 8L) {
      hist_file <- file.path(out_dir, "histogram_by_group.png")
      tryCatch({
        norm_histogram_summary(mat_sub, meta_sub, cfg_mode, out_file = hist_file)
        files <- c(files, hist_file)
      }, error = function(e) {
        skipped <<- c(skipped, sprintf("histogram: %s", conditionMessage(e)))
        status  <<- "partial"
        message(sprintf("[QC][%s][%s] Skipped histogram: %s",
                        stage, subset_mode, conditionMessage(e)))
      })
    } else {
      reason <- sprintf("n_groups=%d (must be 2–8)", n_groups)
      skipped <- c(skipped, sprintf("histogram: %s", reason))
      if (n_groups != 0L) status <- "partial"
    }
  }

  # ---- Heatmaps (sample-count gated) ----------------------------------------
  # Validate and resolve heatmap annotation columns once, shared across all three.
  annot_cols <- cfg_mode$effects$heatmap_annotations
  if (!is.null(annot_cols)) {
    annot_cols <- intersect(as.character(annot_cols), colnames(meta_sub))
    if (length(annot_cols) == 0L) annot_cols <- NULL
  }

  heatmap_eligible <- (n_samp >= 3L && n_samp <= max_heatmap)
  heatmap_requested <- do_hcor || do_hdist || do_hexp

  if (!heatmap_eligible && heatmap_requested) {
    reason <- if (n_samp > max_heatmap) {
      sprintf("n_samples (%d) > max_samples_for_heatmaps (%d)", n_samp, max_heatmap)
    } else {
      sprintf("n_samples (%d) < 3", n_samp)
    }
    skipped <- c(skipped, sprintf("heatmaps (all): %s", reason))
    status  <- "partial"
    message(sprintf("[QC][%s][%s] All heatmaps skipped: %s", stage, subset_mode, reason))
  }

  if (heatmap_eligible) {

    # Correlation heatmap
    if (do_hcor) {
      corr_file <- file.path(out_dir, "heatmap_correlation.png")
      tryCatch({
        qc_sample_correlation_heatmap(mat_sub, meta_sub, cfg_mode,
                                      out_file   = corr_file,
                                      annot_cols = annot_cols)
        files <- c(files, corr_file)
      }, error = function(e) {
        skipped <<- c(skipped, sprintf("heatmap_correlation: %s", conditionMessage(e)))
        status  <<- "partial"
        message(sprintf("[QC][%s][%s] Skipped heatmap_correlation: %s",
                        stage, subset_mode, conditionMessage(e)))
      })
    }

    # Distance heatmap
    if (do_hdist) {
      dist_file <- file.path(out_dir, "heatmap_distance.png")
      tryCatch({
        qc_sample_distance_heatmap(mat_sub, meta_sub, cfg_mode,
                                   out_file   = dist_file,
                                   annot_cols = annot_cols)
        files <- c(files, dist_file)
      }, error = function(e) {
        skipped <<- c(skipped, sprintf("heatmap_distance: %s", conditionMessage(e)))
        status  <<- "partial"
        message(sprintf("[QC][%s][%s] Skipped heatmap_distance: %s",
                        stage, subset_mode, conditionMessage(e)))
      })
    }

    # Expression heatmap (requires >= 5 features for a meaningful display)
    if (do_hexp) {
      if (n_feat >= 5L) {
        hexp_file <- file.path(out_dir, "heatmap_expression.png")
        tryCatch({
          wrap_qc_heatmap(mat_sub, meta_sub, cfg_mode,
                          stage    = stage,
                          out_file = hexp_file)
          files <- c(files, hexp_file)
        }, error = function(e) {
          skipped <<- c(skipped, sprintf("heatmap_expression: %s", conditionMessage(e)))
          status  <<- "partial"
          message(sprintf("[QC][%s][%s] Skipped heatmap_expression: %s",
                          stage, subset_mode, conditionMessage(e)))
        })
      } else {
        skipped <- c(skipped,
          sprintf("heatmap_expression: n_features (%d) < 5", n_feat))
        status <- "partial"
      }
    }
  }

  # ---- Metrics TSV (always written) ------------------------------------------
  tsv <- .qc_write_metrics(
    out_dir     = out_dir,
    stage       = stage,
    subset_mode = subset_mode,
    mat         = mat_sub,
    meta        = meta_sub,
    cfg_mode    = cfg_mode,
    pseudocount = pseudocount,
    skipped     = skipped,
    status      = status,
    pca_var     = pca_var
  )

  c(files, tsv)
}


# ==============================================================================
# mod_met_qc_suite — target-facing wrapper (with_qc + no_qc)
# ==============================================================================

#' Run the full QC suite for one preprocessing stage (both sample subsets)
#'
#' This is the function called directly from \code{{targets}} target
#' expressions.  It creates output directories, resolves config parameters,
#' calls \code{qc_full_metabolomics_suite()} for both \code{"with_qc"} and
#' \code{"no_qc"} subsets, and returns the combined file-path vector for
#' \code{format = "file"} tracking.
#'
#' @param data    List with elements \code{mat} (Log2 numeric matrix),
#'   \code{meta} (data frame), and \code{row_data}.  Matches the return
#'   contract of \code{mod_met_log()}, \code{mod_met_normalize_linear()}, and
#'   \code{mod_met_normalize_log()}.
#' @param stage   Character scalar.  One of \code{"log"}, \code{"norm_tss"},
#'   \code{"norm_median"}, \code{"norm_pqn"}.  Determines the subdirectory
#'   name and is forwarded to \code{qc_full_metabolomics_suite()}.
#' @param out_dir Mode output directory (\code{metab_out_dir}).  The function
#'   creates \code{out_dir/qc/<stage>/with_qc/} and
#'   \code{out_dir/qc/<stage>/no_qc/} beneath it.
#' @param config  Full pipeline config list.  The function extracts:
#'   \itemize{
#'     \item \code{config$modes$metabolomics$qc} (toggles, thresholds)
#'     \item \code{config$modes$metabolomics} (effects, for plot functions)
#'     \item \code{config$modes$metabolomics$normalization$pseudocount}
#'   }
#'
#' @return Character vector of all file paths written by both subset calls
#'   (PNGs + two \code{metrics_summary.tsv} files).  An empty
#'   \code{character(0)} is returned when \code{qc$enabled} is \code{FALSE}.
mod_met_qc_suite <- function(data, stage, out_dir, config) {
  # Skip QC plots when chosen_norm is already set (pass 2)
  chosen <- config$modes$metabolomics$preprocessing$chosen_norm
  if (!is.null(chosen)) {
    message(sprintf("[mod_met_qc_suite] chosen_norm = '%s'; skipping QC for stage '%s'.",
                    chosen, stage))
    return(character(0))
  }

  cfg_qc      <- config$modes$metabolomics$qc           %||% list()
  cfg_mode    <- config$modes$metabolomics
  norm_cfg    <- config$modes$metabolomics$normalization %||% list()
  pseudocount <- norm_cfg$pseudocount %||% 1

  # Short-circuit when QC is globally disabled
  if (!isTRUE(cfg_qc$enabled %||% TRUE)) {
    message(sprintf(
      "[mod_met_qc_suite] QC disabled (qc.enabled = false); skipping stage '%s'.", stage
    ))
    return(character(0L))
  }

  dir_with <- file.path(out_dir, "qc", stage, "with_qc")
  dir_no   <- file.path(out_dir, "qc", stage, "no_qc")
  dir.create(dir_with, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_no,   recursive = TRUE, showWarnings = FALSE)

  common_args <- list(
    mat         = as.matrix(data$mat),
    meta        = as.data.frame(data$meta),
    stage       = stage,
    pseudocount = pseudocount,
    cfg         = cfg_qc,
    cfg_mode    = cfg_mode
  )

  files_with <- do.call(qc_full_metabolomics_suite,
                        c(common_args, list(out_dir = dir_with, subset_mode = "with_qc")))
  files_no   <- do.call(qc_full_metabolomics_suite,
                        c(common_args, list(out_dir = dir_no,   subset_mode = "no_qc")))

  c(files_with, files_no)
}


# ==============================================================================
# mod_met_qc_comparison_table — aggregate benchmark TSV across all stages
# ==============================================================================

#' Aggregate per-stage QC metrics into a cross-stage benchmark table
#'
#' Reads the \code{with_qc/metrics_summary.tsv} from each QC stage target,
#' row-binds them into a single table, adds a \code{pca_pc1_delta} column
#' (PC1 variance gained relative to the pre-normalisation log stage), and
#' saves the result to \code{qc/comparison/normalization_qc_benchmark.tsv}.
#'
#' @param log_qc_files    Character vector returned by the \code{met_log_qc}
#'   file target.
#' @param tss_qc_files    Character vector returned by the
#'   \code{met_norm_tss_qc} file target.
#' @param median_qc_files Character vector returned by the
#'   \code{met_norm_median_qc} file target.
#' @param pqn_qc_files    Character vector returned by the
#'   \code{met_norm_pqn_qc} file target.
#' @param out_dir         Mode output directory (\code{metab_out_dir}).
#' @param config          Full pipeline config list (unused in computation;
#'   included for interface consistency).
#'
#' @return Character scalar: absolute path to the written TSV file.
#'   Suitable for use as the value of a \code{format = "file"} target.
#'
#' @details
#' Only \code{with_qc/metrics_summary.tsv} files are used for the
#' cross-stage comparison (the full-dataset view).  The \code{no_qc} metrics
#' are retained in their per-stage TSVs for manual inspection.
#'
#' \code{pca_pc1_delta} is defined as:
#' \deqn{\Delta\text{PC1} = \text{PC1}_{\text{stage}} - \text{PC1}_{\text{log}}}
#' A positive delta indicates that normalization increased the proportion of
#' variance captured by the first principal component (i.e., reduced unwanted
#' technical variation).
mod_met_qc_comparison_table <- function(log_qc_files, tss_qc_files,
                                        median_qc_files, pqn_qc_files,
                                        out_dir, config,
                                        imputed_data = NULL) {
  # Skip comparison when chosen_norm is already set (pass 2)
  chosen <- config$modes$metabolomics$preprocessing$chosen_norm
  if (!is.null(chosen)) {
    message(sprintf("[mod_met_qc_comparison_table] chosen_norm = '%s'; skipping QC comparison.", chosen))
    return(character(0))
  }

  all_files <- c(log_qc_files, tss_qc_files, median_qc_files, pqn_qc_files)

  # Locate with_qc/metrics_summary.tsv in each file vector.
  # Use a regex that matches the platform path separator (/ or \).
  metrics_paths <- grep(
    "with_qc[\\/]metrics_summary\\.tsv$",
    all_files, value = TRUE, perl = TRUE
  )

  comp_dir <- file.path(out_dir, "qc", "comparison")
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(comp_dir, "normalization_qc_benchmark.tsv")

  # ---- raw_linear row: RSD on met_imputed (linear scale; no back-transform) --
  # The data is already on linear scale, so SD/|mean| is computed directly.
  # rsd_backtransform_exact = TRUE because no back-transformation is needed.
  # pseudocount_used = 0 because the matrix has not been log-transformed.
  # PCA columns are NA because PCA is not computed on linear-scale matrices.
  raw_row <- NULL
  if (!is.null(imputed_data) && !is.null(imputed_data$mat)) {
    mat_lin       <- as.matrix(imputed_data$mat)
    feat_sds      <- apply(mat_lin, 1, stats::sd, na.rm = TRUE)
    feat_means_ab <- abs(rowMeans(mat_lin, na.rm = TRUE))
    feat_means_ab[feat_means_ab == 0 | !is.finite(feat_means_ab)] <- NA_real_
    rsd_vals      <- feat_sds / feat_means_ab
    rsd_vals[!is.finite(rsd_vals)] <- NA_real_

    raw_row <- data.frame(
      stage                   = "raw_linear",
      subset_mode             = "all",
      n_samples               = ncol(mat_lin),
      n_features              = nrow(mat_lin),
      pca_pc1_var_expl        = NA_real_,
      pca_pc2_var_expl        = NA_real_,
      pca_pc3_var_expl        = NA_real_,
      median_rsd_overall      = stats::median(rsd_vals, na.rm = TRUE),
      median_rsd_per_group    = NA_character_,
      rsd_backtransform_exact = TRUE,
      pseudocount_used        = 0,
      plots_skipped           = "",
      status                  = "ok",
      stringsAsFactors        = FALSE
    )
  }

  if (length(metrics_paths) == 0L) {
    warning("mod_met_qc_comparison_table: no with_qc/metrics_summary.tsv files found; writing empty benchmark.")
    empty <- if (!is.null(raw_row)) raw_row else data.frame(stringsAsFactors = FALSE)
    utils::write.table(empty, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
    return(out_file)
  }

  rows <- lapply(metrics_paths, function(f) {
    tryCatch(
      utils::read.table(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                        check.names = FALSE),
      error = function(e) {
        warning(sprintf(
          "mod_met_qc_comparison_table: failed to read '%s': %s", f, e$message
        ))
        NULL
      }
    )
  })
  rows <- Filter(Negate(is.null), rows)

  if (length(rows) == 0L) {
    warning("mod_met_qc_comparison_table: all metrics files failed to parse; writing empty benchmark.")
    empty <- if (!is.null(raw_row)) raw_row else data.frame(stringsAsFactors = FALSE)
    utils::write.table(empty, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
    return(out_file)
  }

  combined <- do.call(rbind, rows)

  # Prepend raw_linear row so it appears first in the benchmark table.
  if (!is.null(raw_row)) combined <- rbind(raw_row, combined)

  # ---- pca_pc1_delta: PC1 variance relative to log stage --------------------
  # Reference is the "log" stage (first log2-scale checkpoint).
  # raw_linear gets NA automatically (NA - ref = NA).
  if ("stage" %in% names(combined) && "pca_pc1_var_expl" %in% names(combined)) {
    log_idx <- which(combined$stage == "log")
    ref_pc1 <- if (length(log_idx) == 1L &&
                   !is.na(combined$pca_pc1_var_expl[log_idx])) {
      combined$pca_pc1_var_expl[log_idx]
    } else {
      NA_real_
    }
    combined$pca_pc1_delta <- combined$pca_pc1_var_expl - ref_pc1
  } else {
    combined$pca_pc1_delta <- NA_real_
  }

  utils::write.table(combined, out_file,
                     sep = "\t", row.names = FALSE, quote = FALSE)
  out_file
}
