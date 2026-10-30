# R/modules/metabolomics/05c_mod_mummichog_report.R
#
# Target-ready wrappers that turn the pinned stage's pathway tables into the
# per-contrast report sections (06e) and write those sections' standalone
# exports to disk.
#
# These live outside 06_mod_report.R on purpose. Both steps used to run inside
# mod_metabolomics_report(), after its `rmarkdown` guard, which made two
# analytical outputs — the supporting evidence and the complementary GSEA —
# conditional on HTML rendering: they were skipped whenever `rmarkdown` was
# absent, and skipped entirely in multiomics mode, where skip_outputs drops the
# report target while the mummichog stage itself still runs in analysis_core.
#
# Wired into 00_pipe_metabolomics.R as metab_mummichog_sections (ordinary
# target, consumed by the report) and metab_mummichog_exports (format = "file",
# a leaf). Both sit in the CORE targets beside the stage they describe, so the
# exports are produced in both single-omics and multiomics mode. See issue #306.


#' Build the per-contrast mummichog report sections
#'
#' Thin wrapper around `build_mummichog_report_sections()` so the sections are a
#' target of their own: computed once, then consumed by both the exports target
#' and the HTML report rather than rebuilt by each.
#'
#' Writes nothing — that is `mod_mummichog_exports()`'s job, and keeping the two
#' apart means a deleted export re-runs only the writer, never the GSEA
#' permutations behind it.
#'
#' @param pathways Named list keyed by contrast, each a mummichog pathway tibble
#'   (`metab_mummichog_report_pathways`); NULL when mummichog is disabled.
#' @param files    Flat character vector of every mummichog output file across
#'   all contrasts (`metab_mummichog_report_files`). Activates the supporting-
#'   evidence layer and the complementary GSEA, both of which read the
#'   EmpiricalCompound tables the pathway tables alone do not carry.
#' @param pre      Preprocessing results; `pre$row_data` supplies the feature
#'   annotations the evidence layer compares against.
#' @param de_res   DE results from `mod_metabolomics_de()`. Supplies both the
#'   ranking statistic and `de_res$method`, which decides which column that is.
#' @param config   Full pipeline config.
#' @param enrichment_res Identified-metabolite enrichment results
#'   (`metab_enrichment_res`), so each contrast's Summary can say whether a
#'   pathway is also supported there. NULL makes every comparable pathway read
#'   "Not tested".
#' @return A named list keyed by contrast — `list(title, subtitle, plot, table,
#'   gsea, gsea_plot, evidence, summary, slug)` per section — exactly as
#'   `build_mummichog_report_sections()` returns it. Empty list when there is
#'   nothing to show.
mod_mummichog_sections <- function(pathways, files, pre, de_res, config,
                                   enrichment_res = NULL) {
    build_mummichog_report_sections(
        pathways, config,
        files    = files,
        row_data = pre$row_data,
        de_res   = de_res,
        enrichment_res = enrichment_res)
}


#' Write the standalone mummichog exports for every contrast
#'
#' Saves each section's presentation exports into `mummichog_pinned/` (plots as
#' PNG + PDF, tables as TSV + CSV: the ORA bubble and pathway table, the GSEA
#' scatter and table, the three evidence tables when they were built, and the
#' Summary), named
#' by the section's de-duplicated slug so labels that collapse to the same token
#' cannot overwrite each other. The engine's own mcg_* result files are
#' untouched.
#'
#' @param sections Named list from `mod_mummichog_sections()`; NULL or empty
#'   when mummichog is disabled or no contrast has a result.
#' @param files    Flat character vector of every mummichog output file
#'   (`metab_mummichog_report_files`), used to find contrasts the stage ran that
#'   produced no section — their stale GSEA exports would otherwise survive.
#' @param out_dir  Metabolomics mode output directory (`metab_out_dir`).
#' @return Character vector of every path written, for the `format = "file"`
#'   target to track; `character(0)` when nothing was written.
mod_mummichog_exports <- function(sections, files, out_dir) {
    written <- character(0)
    for (contrast in names(sections)) {
        s <- sections[[contrast]]
        written <- c(written,
                     save_mummichog_exports(s$plot, s$table, out_dir,
                                            contrast_label = s$slug,
                                            gsea      = s$gsea,
                                            gsea_plot = s$gsea_plot,
                                            evidence  = s$evidence,
                                            summary   = s$summary))
    }

    # A contrast the stage ran but that produced no section never reaches the
    # loop above, so the by-name GSEA and Summary clearing inside
    # save_mummichog_exports()
    # never runs for it and a previous run's table would survive. Clearing is
    # safe for those contrasts because there is nothing of theirs to keep;
    # whether they should get a GSEA analysis at all is a section-contract
    # question, tracked separately.
    clear_orphan_mummichog_gsea_exports(files, sections, out_dir)

    # De-duplicated for the same reason the report module did it before the
    # split: a `format = "file"` target should not be handed the same path
    # twice. The slugs are already made unique, so this should be a no-op.
    unique(written)
}
