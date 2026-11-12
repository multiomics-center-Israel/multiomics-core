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
#'     the group the correction kept, by the primary colour column when that
#'     is another one, and by each batch column.
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
        group_col  = rec$group_col,
        batch_cols = strsplit(rec$batch_cols, ";", fixed = TRUE)[[1]],
        method     = rec$method,
        out_file   = file.path(out_qc, "PCA_batch_before_after.png")
    ))

    # Only what is on disk: the target tracks these paths as files, and a
    # plot whose device or widget could not be saved is not one of them.
    files <- unique(files)
    files[file.exists(files)]
}

#' Colour columns of the before/after PCA panel
#'
#' The group that the correction kept comes first and is the only one
#' labelled as kept. The report's primary colour column (`effects$color`) is
#' shown as well when it is a different column, under its own name.
#'
#' @param cfg The `config$modes$rna` list.
#' @param group_col The group column recorded by the correction.
#' @param batch_cols Batch column(s) that were corrected for.
#' @return A data.frame with `column` and `label` (the text used in the panel
#'   titles), one row per panel row.
#' @examples
#' rna_batch_pca_panel_columns(list(effects = list(color = "Tissue")),
#'                             "Condition", "Batch")
rna_batch_pca_panel_columns <- function(cfg, group_col, batch_cols) {
    columns <- unique(c(group_col, rna_qc_color_values(cfg)[1], batch_cols))
    data.frame(
        column = columns,
        label  = ifelse(columns == group_col,
                        paste0(columns, " (group kept by the correction)"),
                        ifelse(columns %in% batch_cols,
                               paste0(columns, " (batch)"), columns)),
        stringsAsFactors = FALSE
    )
}

#' Genes that both PCA matrices keep
#'
#' `compute_pca_scores()` drops, per matrix, features with fewer than two
#' finite values or zero variance. Run separately on the uncorrected and the
#' corrected matrix, that rule could keep different genes on each side, and
#' the comparison would no longer differ only in the correction. This applies
#' the same rule to both and keeps the genes that pass in both.
#'
#' @param before,after Numeric matrices, genes x samples, with row names.
#' @return Character vector of gene IDs, in the row order of `before`.
#' @examples
#' m <- matrix(c(1, 2, 3, 5, 5, 5), nrow = 2, byrow = TRUE,
#'             dimnames = list(c("g1", "g2"), c("S1", "S2", "S3")))
#' rna_shared_pca_features(m, m)  # "g1" -- g2 has zero variance
rna_shared_pca_features <- function(before, after) {
    usable <- function(m) {
        apply(m, 1, function(x) {
            v <- x[is.finite(x)]
            length(v) >= 2 && stats::var(v) > 0
        })
    }
    ids <- intersect(rownames(before), rownames(after))
    ids[usable(before[ids, , drop = FALSE]) & usable(after[ids, , drop = FALSE])]
}

#' Before/after PCA panel, by kept group and by batch
#'
#' One row per colour column ([rna_batch_pca_panel_columns()]), uncorrected on
#' the left and batch-corrected on the right. It is a visual check of the two
#' questions a correction raises -- does the batch structure fade, does the
#' separation of the kept group remain -- not a test of either.
#'
#' Both sides use the same genes and samples
#' ([rna_shared_pca_features()]), so they differ only in the correction. The
#' two PCAs are still fitted independently: the same genes do not make the
#' axes the same, and positions are not comparable between left and right.
#'
#' @param before,after The QC matrix before and after correction.
#' @param meta Sample metadata.
#' @param cfg The `config$modes$rna` list.
#' @param group_col The group column the correction kept (from its record).
#' @param batch_cols Batch column(s) that were corrected for.
#' @param method Correction method that ran, for the panel titles.
#' @param out_file PNG path.
#' @return `out_file`.
draw_rna_batch_pca_panel <- function(before, after, meta, cfg, group_col, batch_cols,
                                     method, out_file) {
    if (!setequal(colnames(before), colnames(after))) {
        stop("[batch_correction] Before/after PCA: the uncorrected and the corrected QC ",
             "matrix hold different samples.", call. = FALSE)
    }
    ids <- rna_shared_pca_features(before, after)
    before <- before[ids, , drop = FALSE]
    after  <- after[ids, colnames(before), drop = FALSE]

    panel_cols <- rna_batch_pca_panel_columns(cfg, group_col, batch_cols)
    # A numeric label (e.g. pair 1, 2, 3) would get a colour gradient; the
    # kept group and the batches are categories.
    for (col in c(group_col, batch_cols)) meta[[col]] <- factor(as.character(meta[[col]]))

    panels <- list()
    for (i in seq_len(nrow(panel_cols))) {
        col <- panel_cols$column[i]
        panel_cfg <- cfg
        panel_cfg$effects$color <- col
        # Colour and shape may not be the same column.
        if (col %in% as.character(unlist(panel_cfg$effects$shape))) {
            panel_cfg$effects$shape <- NULL
        }
        panels <- c(panels, list(
            qc_pca_scatter(before, meta, panel_cfg, pcs = c(1, 2)) +
                ggplot2::ggtitle(sprintf("Uncorrected — by %s", panel_cols$label[i])),
            qc_pca_scatter(after, meta, panel_cfg, pcs = c(1, 2)) +
                ggplot2::ggtitle(sprintf("Batch-corrected (%s) — by %s", method,
                                         panel_cols$label[i]))
        ))
    }
    combined <- patchwork::wrap_plots(panels, ncol = 2)
    ggplot2::ggsave(out_file, plot = combined, width = 12, height = 5 * nrow(panel_cols))
    out_file
}
