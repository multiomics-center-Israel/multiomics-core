#' RNA-seq QC plots of the batch-corrected matrices
#'
#' The "after" side of the before/after comparison. The "before" side is the
#' established QC output of [mod_rnaseq_qc_pre()] and [mod_rnaseq_qc_post()],
#' which this module does not touch. It draws, for the batch-corrected copies
#' ([compute_rna_batch_correction()]):
#'
#'   * every sample-structure plot of the QC matrix, by the same code as the
#'     originals ([draw_rna_sample_structure_plots()]), under the same names
#'     plus `_batch_corr_<method>`, with titles that say batch-corrected;
#'   * the top-DE heatmap of the analysis matrix
#'     ([draw_rna_top_de_heatmap()]), same suffix;
#'   * `PCA_batch_before_after.png`: PC1 vs PC2 before and after, coloured by
#'     the primary colour column and by each batch column.
#'
#' Nothing is drawn unless `batch_corrected$ran` records that a correction ran
#' -- the runtime evidence, not the config that requested it. Plots of this
#' module left by an earlier run are removed first either way, so they never
#' outlive their correction.
#'
#' @param pre List returned by `preprocess_rna()` (the uncorrected matrices).
#' @param batch_corrected Output of [compute_rna_batch_correction()].
#' @param de_res DE results list from `mod_rnaseq_de()`.
#' @param config Full config list.
#' @param out_dir RNA output directory.
#' @return Paths of the written plots, or `character(0)` when no correction ran.
mod_rnaseq_batch_qc <- function(pre, batch_corrected, de_res, config, out_dir) {
    out_qc <- create_legacy_output_dirs(out_dir, create = TRUE)$diagnostic_plots
    stale <- c(
        list.files(out_qc, pattern = "_batch_corr_(limma|combat)\\.(png|html)$",
                   full.names = TRUE),
        file.path(out_qc, "PCA_batch_before_after.png")
    )
    stale <- stale[file.exists(stale)]
    if (length(stale) > 0) file.remove(stale)

    if (!isTRUE(batch_corrected$ran)) {
        message("[mod_rnaseq_batch_qc] No batch correction ran -- no batch-corrected QC plots.")
        return(character(0))
    }

    assert_pre_contract(pre, stage = "rna")
    cfg <- config$modes$rna
    validate_effects_config(cfg$effects)
    meta <- pre$meta
    rec <- batch_corrected$record
    suffix <- paste0("_batch_corr_", rec$method)
    title_tag <- sprintf("batch-corrected (%s)", rec$method)

    files <- draw_rna_sample_structure_plots(
        batch_corrected$qc_expr, meta, cfg, out_qc,
        stage = "rna", suffix = suffix, title_tag = title_tag
    )$files

    # Top-DE heatmap, when the DE summary can be built (as in mod_rnaseq_qc_post).
    sumdf <- de_res$summary_df %||%
        tryCatch(build_rnaseq_summary_df(de_res$tables, cfg$de), error = function(e) NULL)
    if (!is.null(sumdf)) {
        files <- c(files, draw_rna_top_de_heatmap(
            batch_corrected$expr_work, meta, sumdf, cfg,
            out_file = file.path(out_qc, sprintf("heatmap_top_DE%s.png", suffix)),
            title_tag = title_tag
        ))
    }

    files <- c(files, draw_rna_batch_pca_panel(
        before     = rna_qc_expr(pre, cfg$qc),
        after      = batch_corrected$qc_expr,
        meta       = meta,
        cfg        = cfg,
        batch_cols = strsplit(rec$batch_cols, ";", fixed = TRUE)[[1]],
        method     = rec$method,
        out_file   = file.path(out_qc, "PCA_batch_before_after.png")
    ))

    # Only what is on disk: the target tracks these paths as files, and a
    # plot whose device or widget could not be saved is not one of them.
    files <- unique(files)
    files[file.exists(files)]
}

#' Before/after PCA panel, by condition and by batch
#'
#' One row per colour column (the primary colour column, then each batch
#' column), uncorrected on the left and batch-corrected on the right. It is a
#' visual check of the two questions a correction raises -- does the batch
#' structure fade, does the group separation remain -- not a test of either.
#' PCA is fitted separately on each side, so axes are not comparable across
#' panels.
#'
#' @param before,after The QC matrix before and after correction (same
#'   features and samples).
#' @param meta Sample metadata.
#' @param cfg The `config$modes$rna` list.
#' @param batch_cols Batch column(s) that were corrected for.
#' @param method Correction method that ran, for the panel titles.
#' @param out_file PNG path.
#' @return `out_file`.
draw_rna_batch_pca_panel <- function(before, after, meta, cfg, batch_cols, method, out_file) {
    # A numeric batch label (e.g. pair 1, 2, 3) would get a colour gradient;
    # batches are categories.
    for (b in batch_cols) meta[[b]] <- factor(as.character(meta[[b]]))

    color_cols <- unique(c(rna_qc_color_values(cfg)[1], batch_cols))
    panels <- list()
    for (col in color_cols) {
        panel_cfg <- cfg
        panel_cfg$effects$color <- col
        # Colour and shape may not be the same column.
        if (col %in% as.character(unlist(panel_cfg$effects$shape))) {
            panel_cfg$effects$shape <- NULL
        }
        panels <- c(panels, list(
            qc_pca_scatter(before, meta, panel_cfg, pcs = c(1, 2)) +
                ggplot2::ggtitle(sprintf("Uncorrected — by %s", col)),
            qc_pca_scatter(after, meta, panel_cfg, pcs = c(1, 2)) +
                ggplot2::ggtitle(sprintf("Batch-corrected (%s) — by %s", method, col))
        ))
    }
    combined <- patchwork::wrap_plots(panels, ncol = 2)
    ggplot2::ggsave(out_file, plot = combined, width = 12, height = 5 * length(color_cols))
    out_file
}
