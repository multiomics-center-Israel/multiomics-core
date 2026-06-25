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
