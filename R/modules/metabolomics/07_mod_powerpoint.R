# R/modules/metabolomics/07_mod_powerpoint.R
#
# Metabolomics PowerPoint module: wraps generate_metabolomics_pptx()
# for targets integration.


#' Generate metabolomics summary PowerPoint presentation
#'
#' @param pre             List from preprocess_metabolomics().
#' @param qc_res          List from mod_metabolomics_qc_pre().
#' @param de_res          List from mod_metabolomics_de().
#' @param feature_sel_res List from mod_metabolomics_feature_selection() (or NULL).
#' @param enrichment_res  List from mod_metabolomics_enrichment() (or NULL).
#' @param config          Full pipeline config.
#' @param out_dir         Output directory for this mode.
#' @return Character path to the generated .pptx file.
mod_metabolomics_powerpoint <- function(pre, qc_res, de_res, feature_sel_res,
                                         enrichment_res, config, out_dir) {
    generate_metabolomics_pptx(
        pre             = pre,
        qc_res          = qc_res,
        de_res          = de_res,
        feature_sel_res = feature_sel_res,
        enrichment_res  = enrichment_res,
        config          = config,
        out_dir         = out_dir
    )
}
