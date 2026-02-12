#' Module wrapper: RNA-seq Batch Correction
#'
#' Returns modified pre (or original if disabled).
#' Delegates to R/domain/rnaseq/04_batch_correction.R
#'
#' @param pre Preprocessed data list
#' @param config Full config
#' @param out_dir Output directory
#' @return pre list (potentially with corrected expr_work/expr_filt)
mod_rnaseq_batch_correction <- function(pre, config, out_dir) {
    rna_cfg <- config$modes$rna

    bc_cfg <- rna_cfg$batch_correction %||% list(enabled = FALSE)

    if (!isTRUE(bc_cfg$enabled)) {
        message("[mod_batch_correction] Batch correction disabled. Passing pre through unchanged.")
        return(pre)
    }

    message("[mod_batch_correction] Running batch correction...")
    run_batch_analysis(pre, rna_cfg, out_dir)
}
