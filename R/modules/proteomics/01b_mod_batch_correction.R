#' Module wrapper: Proteomics batch correction
#'
#' Thin orchestration layer over `correct_batch_proteomics()`
#' (R/domain/proteomics/04b_batch_correction.R). Returns the (possibly
#' corrected) `pre` object so downstream DE / QC / clustering consume the
#' corrected matrices. No file I/O here.
#'
#' @param pre Output of `preprocess_proteomics()`.
#' @param config Full pipeline config.
#' @param verbose Logical; emit a status message when active.
#' @return The `pre` list, unchanged when batch correction is disabled.
mod_proteomics_batch_correction <- function(pre, config, verbose = FALSE) {
  assert_pre_contract(pre, stage = "proteomics")

  bc <- config$modes$proteomics$batch_correction %||% list()
  method <- bc$method %||% "none"
  enabled <- bc$enabled %||% (!identical(method, "none"))

  if (!isTRUE(enabled) || identical(method, "none")) {
    if (isTRUE(verbose)) {
      message("[mod_proteomics_batch_correction] Disabled. Passing pre through unchanged.")
    }
    return(pre)
  }

  if (isTRUE(verbose)) {
    message(sprintf("[mod_proteomics_batch_correction] Running '%s' correction.", method))
  }
  correct_batch_proteomics(pre, config)
}


#' Before/after batch-correction PCA diagnostic
#'
#' Draws one PNG with four PC1-vs-PC2 panels: uncorrected and batch-corrected,
#' each coloured by condition (`effects.color`) and by the batch column. It is a
#' visual check of the two questions a correction raises -- does the batch
#' structure fade, does the condition separation remain -- not a test of either.
#' PCA is fitted separately on each side, so axes are not comparable across
#' panels.
#'
#' The two sides differ only in the correction. `pre_raw` is the input that
#' `correct_batch_proteomics()` was fitted on, so both matrices hold the same
#' samples and the same imputation draw; the feature set is fixed by
#' [select_shared_pca_features()] so neither PCA drops a protein the other keeps.
#'
#' Nothing is drawn unless `pre$batch_corrected` records that a correction ran
#' -- the runtime evidence, not the config that requested it. The PNG is removed
#' first either way, so a panel from an earlier run never outlives its
#' correction.
#'
#' @param pre_raw `pre` before batch correction (the `prot_pre_raw` target).
#' @param pre `pre` after [mod_proteomics_batch_correction()].
#' @param config Full pipeline config.
#' @param out_dir Proteomics output directory.
#' @return Path of the written PNG, or `character(0)` when no correction ran.
mod_proteomics_batch_pca <- function(pre_raw, pre, config, out_dir) {
  out_qc <- create_legacy_output_dirs(out_dir)$diagnostic_plots
  f_out <- file.path(out_qc, "PCA_batch_before_after.png")
  if (file.exists(f_out)) file.remove(f_out)

  if (!isTRUE(pre$batch_corrected)) {
    message("[mod_proteomics_batch_pca] No batch correction ran -- no before/after PCA.")
    return(character(0))
  }

  cfg <- config$modes$proteomics
  batch_col <- get_proteomics_batch_config(cfg)$batch_col

  # The correction returns the matrix it was given, so the samples should
  # already match; check rather than assume, since a mismatch would make the
  # two sides differ by more than the correction.
  samples <- colnames(pre_raw$expr_imp_single)
  if (!setequal(samples, colnames(pre$expr_imp_single))) {
    cli::cli_abort(c(
      "Before/after batch PCA: the uncorrected and corrected matrices hold different samples.",
      "i" = "Only in uncorrected: {.val {setdiff(samples, colnames(pre$expr_imp_single))}}.",
      "i" = "Only in corrected: {.val {setdiff(colnames(pre$expr_imp_single), samples)}}.",
      "i" = "Both are expected to come from the same preprocessing run (prot_pre_raw -> prot_pre)."
    ))
  }
  ids <- select_shared_pca_features(pre_raw$expr_imp_single, pre$expr_imp_single)
  before <- pre_raw$expr_imp_single[ids, samples, drop = FALSE]
  after  <- pre$expr_imp_single[ids, samples, drop = FALSE]

  # A numeric batch label (e.g. pair_no 1, 2, 3) would get a colour gradient;
  # batches are categories.
  meta <- pre$meta
  meta[[batch_col]] <- factor(meta[[batch_col]])
  cfg_by_batch <- cfg
  cfg_by_batch$effects$color <- batch_col

  panel <- function(mat, panel_cfg, title) {
    qc_pca_scatter(mat, meta, panel_cfg, pcs = c(1, 2)) + ggplot2::ggtitle(title)
  }
  cond_label <- cfg$effects$color
  panels <- list(
    panel(before, cfg,          sprintf("Uncorrected — by %s", cond_label)),
    panel(after,  cfg,          sprintf("Batch-corrected — by %s", cond_label)),
    panel(before, cfg_by_batch, sprintf("Uncorrected — by %s", batch_col)),
    panel(after,  cfg_by_batch, sprintf("Batch-corrected — by %s", batch_col))
  )
  combined <- patchwork::wrap_plots(panels, ncol = 2)
  ggplot2::ggsave(f_out, plot = combined, width = 12, height = 10)
  f_out
}


#' Proteins that both PCA matrices keep
#'
#' `compute_pca_scores()` drops, per matrix, features with fewer than two
#' finite values or zero variance. Run separately on the uncorrected and the
#' corrected matrix, that rule could keep different proteins on each side, and
#' the comparison would no longer differ only in the correction. This applies
#' the same rule to both and keeps the proteins that pass in both.
#'
#' @param before,after Numeric matrices, proteins x samples, with row names.
#' @return Character vector of protein IDs, in the row order of `before`.
#' @examples
#' m <- matrix(c(1, 2, 3, 5, 5, 5), nrow = 2, byrow = TRUE,
#'             dimnames = list(c("P1", "P2"), c("S1", "S2", "S3")))
#' select_shared_pca_features(m, m)  # "P1" -- P2 has zero variance
select_shared_pca_features <- function(before, after) {
  usable <- function(m) {
    apply(m, 1, function(x) {
      v <- x[is.finite(x)]
      length(v) >= 2 && stats::var(v) > 0
    })
  }
  ids <- intersect(rownames(before), rownames(after))
  ids[usable(before[ids, , drop = FALSE]) & usable(after[ids, , drop = FALSE])]
}
