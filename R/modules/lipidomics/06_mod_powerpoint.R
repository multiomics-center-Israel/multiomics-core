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
#' @return Character path to the generated .pptx file.
mod_lipidomics_powerpoint <- function(pre, qc_res, de_res, feature_sel_res,
                                       class_res, config, out_dir) {
    generate_lipidomics_pptx(
        pre             = pre,
        qc_res          = qc_res,
        de_res          = de_res,
        feature_sel_res = feature_sel_res,
        class_res       = class_res,
        config          = config,
        out_dir         = out_dir
    )
}
