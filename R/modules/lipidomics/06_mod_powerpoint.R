# R/modules/lipidomics/06_mod_powerpoint.R
#
# Lipidomics PowerPoint module: wraps generate_lipidomics_pptx()
# for targets integration.


#' Generate lipidomics summary PowerPoint presentation
#'
#' @param pre             List from preprocess_lipidomics().
#' @param qc_res          List from mod_lipidomics_qc_pre().
#' @param de_res          List from mod_lipidomics_de().
#' @param feature_sel_res List from mod_lipidomics_feature_selection() (or NULL).
#' @param class_res       List from mod_lipidomics_lipid_class() (or NULL).
#' @param config          Full pipeline config.
#' @param out_dir         Output directory for this mode.
#' @param biomarker_res   List from mod_lipidomics_biomarker() (or NULL).
#' @param pathway_res     List from mod_lipidomics_pathway() (or NULL).
#' @param qc_enhanced     List from mod_lipidomics_qc_enhanced() (or NULL).
#' @param clustering_res  List from mod_lipidomics_clustering() (or NULL).
#' @param pool_cv         List from compute_pool_cv() target (or NULL).
#' @return Character path to the generated .pptx file.
mod_lipidomics_powerpoint <- function(pre, qc_res, de_res, feature_sel_res,
                                       class_res, config, out_dir,
                                       biomarker_res  = NULL,
                                       pathway_res    = NULL,
                                       qc_enhanced    = NULL,
                                       clustering_res = NULL,
                                       pool_cv        = NULL) {
    generate_lipidomics_pptx(
        pre             = pre,
        qc_res          = qc_res,
        de_res          = de_res,
        feature_sel_res = feature_sel_res,
        class_res       = class_res,
        config          = config,
        out_dir         = out_dir,
        biomarker_res   = biomarker_res,
        pathway_res     = pathway_res,
        qc_enhanced     = qc_enhanced,
        clustering_res  = clustering_res,
        pool_cv         = pool_cv
    )
}
