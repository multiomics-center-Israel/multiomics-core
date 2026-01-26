mod_rnaseq_qc_post <- function(pre, de_res, config, out_dir) {
  assert_pre_contract(pre, stage = "rna")
  dirs <- create_legacy_output_dirs(out_dir, create = TRUE)
  
  mat <- pre$expr_work
  meta <- pre$meta
  cfg <- config$modes$rna
  eff_color <- cfg$effects$color %||% "Group"
  
  sumdf <- de_res$summary_df
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
