# R/modules/proteomics/11_mod_powerpoint.R
#
# Proteomics PowerPoint module wrapper for targets integration.

#' Generate proteomics summary PowerPoint presentation
#'
#' @param pre          List from preprocess.
#' @param qc_res       QC results list.
#' @param de_res       DE results list.
#' @param pathway_res  Pathway enrichment results (or NULL).
#' @param ppi_res      PPI network results (or NULL).
#' @param config       Full pipeline config.
#' @param out_dir      Output directory for this mode.
#' @return Character path to the generated .pptx file.
mod_proteomics_powerpoint <- function(pre, qc_res, de_res, pathway_res,
                                       ppi_res, config, out_dir) {
    generate_proteomics_pptx(
        pre         = pre,
        qc_res      = qc_res,
        de_res      = de_res,
        pathway_res = pathway_res,
        ppi_res     = ppi_res,
        config      = config,
        out_dir     = out_dir
    )
}
