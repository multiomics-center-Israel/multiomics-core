# R/domain/lipidomics/03_differential.R
#
# Differential expression analysis for lipidomics data.
# Reuses the metabolomics DE engine (de_limma, de_t_test, de_wilcoxon,
# build_de_summary, parse_metab_contrast, make_contrast_label,
# extract_contrast_table) but reads from config$modes$lipidomics.


# ---- public entry point -----------------------------------------------------

#' Run lipidomics differential analysis
#'
#' @param pre    List from preprocess_lipidomics() (pre contract).
#' @param config Full pipeline config.
#' @return list conforming to the DE contract:
#'   summary_df, method, de_tables, de_model
run_lipidomics_de <- function(pre, config) {
    cfg <- config$modes$lipidomics
    de_cfg <- cfg$de %||% list()

    method <- tolower(de_cfg$method %||% "limma")
    assert_one_of(method, "de$method", c("limma", "t_test", "t_test_equal", "wilcoxon"))

    condition_col <- de_cfg$condition_column %||% cfg$effects$color %||% "sample_type"
    sample_col <- cfg$effects$samples %||% "sample_id"

    contrasts <- de_cfg$contrasts
    if (is.null(contrasts) || length(contrasts) == 0) {
        stop("lipidomics DE: config$modes$lipidomics$de$contrasts is required.")
    }
    if (is.list(contrasts)) contrasts <- unlist(contrasts)

    mat  <- pre$expr_work
    mat_for_test <- pre$expr_log %||% mat
    meta <- pre$meta
    assert_numeric_matrix(mat_for_test, "lipid_expr_for_test")

    meta <- meta[match(colnames(mat_for_test), meta[[sample_col]]), , drop = FALSE]
    condition <- factor(meta[[condition_col]])

    # Thresholds
    padj_cutoff <- de_cfg$p_cutoff %||% 0.05
    if (!is.null(de_cfg$logfc_cutoff)) {
        log2fc_cut <- de_cfg$logfc_cutoff
    } else {
        linear_fc  <- de_cfg$linear_fc_cutoff %||% 1.5
        log2fc_cut <- log2(linear_fc)
    }

    mat_raw <- pre$expr_filt

    de_tables <- list()
    de_model  <- NULL

    for (ctr in contrasts) {
        ctr_label <- make_contrast_label(ctr)
        message("lipidomics DE [", method, "]: ", ctr)

        tbl <- switch(method,
            limma        = de_limma(mat_for_test, condition, ctr, mat_for_fc = mat_raw),
            t_test       = de_t_test(mat_for_test, condition, ctr, mat_for_fc = mat_raw),
            t_test_equal = de_t_test_equal(mat_for_test, condition, ctr, mat_for_fc = mat_raw),
            wilcoxon     = de_wilcoxon(mat_for_test, condition, ctr, mat_for_fc = mat_raw)
        )

        if (method == "limma" && is.null(de_model)) {
            de_model <- attr(tbl, "fit")
        }

        de_tables[[ctr_label]] <- tbl
    }

    summary_df <- build_de_summary(de_tables, padj_cutoff, log2fc_cut)

    message("lipidomics DE complete: ", nrow(summary_df), " features, ",
            sum(summary_df$pass_any_contrast == 1, na.rm = TRUE), " significant")

    list(
        summary_df = summary_df,
        method     = method,
        de_tables  = de_tables,
        de_model   = de_model
    )
}
