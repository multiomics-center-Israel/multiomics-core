# R/modules/rnaseq/12_mod_fact_sheet.R
#
# Target-ready wrapper around build_rnaseq_fact_sheet().

#' Write the RNA-seq results fact sheet
#'
#' Assembles the run's headline numbers, each paired with the artifact it can be
#' checked against, and writes them to \code{results_fact_sheet.tsv} in the
#' mode's results directory.
#'
#' Runs after the legacy outputs and the enrichment module so that the files it
#' reads are already on disk; those targets are taken as arguments purely to
#' order the pipeline, not for their values.
#'
#' @param pre Preprocessing results.
#' @param inputs Loaded inputs (for contrast names).
#' @param config Full pipeline config.
#' @param out_dir The RNA-seq results directory.
#' @param run_dir Run root, for the execution_info rows.
#' @param outputs_legacy Ignored; declares the dependency on the written TSVs.
#' @param pathway_res Ignored; declares the dependency on the enrichment files.
#' @return Path to the written fact sheet, or character(0) if nothing was written.
#' @export
mod_rnaseq_fact_sheet <- function(pre, inputs, config, out_dir, run_dir = NULL,
                                  outputs_legacy = NULL, pathway_res = NULL) {
    sheet <- tryCatch(
        build_rnaseq_fact_sheet(
            out_dir = out_dir, config = config, pre = pre,
            inputs = inputs, run_dir = run_dir
        ),
        error = function(e) {
            # A fact sheet is a convenience artifact. It must never be the
            # reason a completed run reports failure.
            warning(sprintf("RNA-seq fact sheet: could not assemble it (%s); skipping.",
                            conditionMessage(e)))
            NULL
        }
    )
    if (is.null(sheet) || nrow(sheet) == 0) return(character(0))

    message(sprintf("[RNA-seq] Fact sheet: %d numbers recorded with their source files.", nrow(sheet)))
    write_results_fact_sheet(sheet, out_dir)
}
