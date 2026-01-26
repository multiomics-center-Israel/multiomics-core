mod_rnaseq_qc_post <- function(pre, de_res, config, out_dir) {
    assert_pre_contract(pre, stage = "rna")
    dirs <- create_legacy_output_dirs(out_dir, create = TRUE)

    mat <- pre$expr_work
    meta <- pre$meta
    cfg <- config$modes$rna
    eff_color <- cfg$effects$color %||% "Group"

    sumdf <- de_res$summary_df

    if (is.null(sumdf)) {
        # If using output of run_deseq2_de (which returns tables list, not summary_df directly in the list usually?),
        # we need to ensure summary_df is passed correctly.
        # run_deseq2_de returns list(method, tables, contrasts, info).
        # summary_df is built via build_rnaseq_summary_df in the pipeline.
        # The pipeline step 'rna_de_res' returns the DE obj.
        # Wait, in 00_pipe_rnaseq.R, 'rna_de_res' is the result of mod_rnaseq_de.
        # 'mod_rnaseq_de' (legacy) returned list(method, tables, contrasts, info).
        # 'mod_rnaseq_qc_post' (legacy) accepted 'de_res'.
        # But legacy 'mod_rnaseq_qc_post' used 'de_res$summary_df' (line 10).
        # Where did summary_df come from in legacy pipeline?
        # In legacy 00_pipe_rnaseq.R (Step 520 lines 17 & 33), 'rna_de_res' is passed to mod_rnaseq_qc_post.
        # BUT 'rna_de_res' comes from 'mod_rnaseq_de'.
        # 'mod_rnaseq_de' (Step 461) DOES NOT return summary_df.
        # This implies there was a bug or I missed where summary_df is attached.
        # Actually, in Step 520, line 21 calls 'write_rnaseq_outputs_legacy'.
        # 'write_rnaseq_outputs_legacy' calls 'build_rnaseq_summary_df'.
        # But 'mod_rnaseq_qc_post' is a separate target 'rna_qc_post_obj' that takes 'rna_de_res'.
        # If 'rna_de_res' doesn't have 'summary_df', mod_rnaseq_qc_post will fail or 'rna_de_res' was modified?
        # Targets are immutable.
        # Ah, maybe I should construct summary_df inside mod_rnaseq_qc_post if missing?
        # Or maybe the legacy code I read in 461/520 was incomplete/broken or I misread it.
        # In Step 520 line 36: de_res  = rna_de_res.
        # In Step 519 line 10: sumdf <- de_res$summary_df.
        # This confirms the legacy code EXPECTED summary_df.
        # But mod_rnaseq_de (Step 461) returned list(method, tables, contrasts, info).
        # THIS IS A BUG in the legacy code or my read of it.
        # I will FIX it in the refactor.
        # I will make mod_rnaseq_qc_post robust: build summary_df if missing using build_rnaseq_summary_df from domain.

        # Correction: I'll use build_rnaseq_summary_df here if de_res$summary_df is null.
        # Assuming config is available.
        sumdf <- tryCatch(build_rnaseq_summary_df(de_res, config), error = function(e) NULL)
        if (is.null(sumdf)) {
            warning("RNA QC-post: Could not build summary_df. Skipping heatmap.")
            return(list(files = character(0)))
        }
    }

    # Recalculate de_ids logic
    de_ids <- sumdf$FeatureID[sumdf$pass_any_contrast == 1]

    if (length(de_ids) < 2) {
        warning("RNA QC-post: <2 DE features; skipping heatmap.")
        return(list(files = character(0)))
    }

    # choose top N by min padj across contrasts
    padj_cols <- grep("^padj\\.", names(sumdf), value = TRUE)
    min_padj <- apply(sumdf[, padj_cols, drop = FALSE], 1, function(x) {
        if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
    })
    ord <- order(min_padj, na.last = NA)
    top_n <- min(2000L, length(ord))
    top_ids <- sumdf$FeatureID[ord][seq_len(top_n)]

    m <- mat[top_ids, , drop = FALSE]
    # z-score by gene
    m <- t(scale(t(m)))
    m <- m[complete.cases(m), , drop = FALSE]

    ann <- data.frame(Color = meta[[eff_color]])
    rownames(ann) <- rownames(meta)

    hm_png <- file.path(dirs$diagnostic_plots, "heatmap_top_DE.png")
    grDevices::png(hm_png, width = 1400, height = 1200, res = 150)
    pheatmap::pheatmap(
        m,
        show_rownames = FALSE,
        annotation_col = ann,
        main = sprintf("Top DE (%d)", nrow(m)),
        cluster_cols = TRUE
    )
    grDevices::dev.off()

    list(files = c(hm_png))
}
