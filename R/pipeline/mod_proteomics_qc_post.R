# R/pipeline/mod_proteomics_qc_post.R

#' Proteomics QC post-DE module
#'
#' Writes QC artifacts that depend on DE results.
#' Supports DE source:
#' - "summary": builds per-contrast tables from summary_df (wide)
#' - "table1" : uses first imputation run tables (legacy-ish)
#'
#' @param pre        preprocessed proteomics object
#' @param de_res     DE results list
#' @param config     full YAML config list
#' @param run_dir    output run directory
#' @param de_source  "summary" or "table1"
#' @return list(plots, files)
mod_proteomics_qc_post <- function(pre, de_res, config, run_dir, de_source = c("summary", "table1")) {
  de_source <- match.arg(de_source)
  
  cfg <- config$modes$proteomics
  
  dirs <- create_legacy_output_dirs(run_dir)
  out_qc_post <- file.path(dirs$diagnostic_plots, "QC_post")
  ensure_dir(out_qc_post)
  
  plots <- list()
  files <- character(0)
  
  # always write the summary if available
  if (!is.null(de_res$summary_df)) {
    f_all <- file.path(out_qc_post, "de_summary_all.tsv")
    save_tsv_path(de_res$summary_df, f_all)
    files <- c(files, f_all)
  }
  
  # unified: get tables per contrast (list)
  tables <- get_de_tables_qc_post(de_res, cfg, de_source = de_source)

  
  # write per-contrast tables (and later: volcano)
  for (cn in names(tables)) {
    de_tbl <- tables[[cn]]
    if (is.null(de_tbl)) next
    
    p <- plot_volcano(de_tbl, cfg = cfg, title = paste0("Volcano: ", cn), use_adj = TRUE)
    
    f_png <- file.path(out_qc_post, sprintf("volcano_%s.png", cn))
    ggplot2::ggsave(f_png, plot = p, width = 8, height = 6, dpi = 150)
    
    plots[[paste0("volcano_", cn)]] <- p
    files <- c(files, f_png)
  }
  
  list(plots = plots, files = unique(files))
}


# ---- helpers for QC_post ----

get_de_tables_qc_post <- function(de_res, cfg, de_source = c("summary", "table1"), use_adj = TRUE) {
  de_source <- match.arg(de_source)
  
  if (de_source == "summary") {
    if (is.null(de_res$summary_df)) stop("QC_post: de_res$summary_df is missing.")
    return(qc_post_tables_from_summary(de_res$summary_df, cfg, use_adj = use_adj))
  }
  
  # table1: list(contrast -> limma topTable) already, so return as-is
  if (is.null(de_res$runs_de_tables) || length(de_res$runs_de_tables) < 1) {
    stop("QC_post: de_source='table1' requested but de_res$runs_de_tables is missing.")
  }
  de_res$runs_de_tables[[1]]
}


# returns named list: contrast -> de_tbl (with logFC + P.Value/adj.P.Val)
qc_post_tables_from_summary <- function(summary_df, cfg, use_adj = TRUE) {
  stopifnot(is.data.frame(summary_df))
  id_col <- cfg$de_table$id_col %||% "FeatureID"
  if (!(id_col %in% colnames(summary_df))) stop("summary_df missing id col: ", id_col)
  
  contrasts <- sub("^padj\\.imputs\\.", "", grep("^padj\\.imputs\\.", colnames(prot_de_res$summary_df), value = TRUE))
  if (length(contrasts) == 0) stop("No contrasts found in summary_df (expected columns like padj.imputs.<contrast>).")
  
  out <- setNames(vector("list", length(contrasts)), contrasts)
  
  for (cn in contrasts) {
    p_col   <- if (isTRUE(use_adj)) paste0("padj.imputs.", cn) else paste0("pvalue.imputs.", cn)
    fc_col  <- paste0("linearFC.imputs.", cn)   # this is linear fold-change in your summary
    if (!(p_col %in% colnames(summary_df))) stop("summary_df missing: ", p_col)
    if (!(fc_col %in% colnames(summary_df))) stop("summary_df missing: ", fc_col)
    
    de_tbl <- data.frame(
      FeatureID = summary_df[[id_col]],
      logFC = signed_fc_to_log2(summary_df[[fc_col]]),
      P.Value   = as.numeric(summary_df[[paste0("pvalue.imputs.", cn)]]),
      adj.P.Val = as.numeric(summary_df[[paste0("padj.imputs.", cn)]]),
      stringsAsFactors = FALSE
    )
    
    # optional annotations if exist
    for (a in intersect(c("Protein.Names","Genes","First.Protein.Description"), colnames(summary_df))) {
      de_tbl[[a]] <- summary_df[[a]]
    }
    
    out[[cn]] <- de_tbl
  }
  
  out
}

