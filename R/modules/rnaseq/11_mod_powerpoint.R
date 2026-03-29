# R/modules/rnaseq/11_mod_powerpoint.R
#
# RNA-seq PowerPoint module wrapper for targets integration.

#' Generate RNA-seq summary PowerPoint presentation
#'
#' @param pre            List from preprocess.
#' @param qc_pre_obj     QC pre results.
#' @param de_res         DE results list.
#' @param clustering_obj Clustering results (or NULL).
#' @param pathway_res    Pathway enrichment results (or NULL).
#' @param config         Full pipeline config.
#' @param out_dir        Output directory for this mode.
#' @return Character path to the generated .pptx file.
mod_rnaseq_powerpoint <- function(pre, qc_pre_obj, de_res,
                                    clustering_obj, pathway_res,
                                    config, out_dir) {
    generate_rnaseq_pptx(
        pre            = pre,
        qc_pre_obj     = qc_pre_obj,
        de_res         = de_res,
        clustering_obj = clustering_obj,
        pathway_res    = pathway_res,
        config         = config,
        out_dir        = out_dir
    )
}
