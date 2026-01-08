# R/plots/plot_de.R

#' Volcano plot for a single DE table (one contrast)
#'
#' Expects columns:
#'   - logFC
#'   - P.Value and/or adj.P.Val
#'
#' Uses cfg$de$p_cutoff and cfg$de$linear_fc_cutoff (linear FC) to mark "pass".
#' Converts linear_fc_cutoff to log2FC cutoff.
#'
#' @param de_tbl data.frame with DE results for one contrast
#' @param cfg mode config (e.g., config$modes$proteomics)
#' @param title optional plot title
#' @param use_adj logical; if TRUE uses adj.P.Val, else P.Value
#'
#' @return ggplot object
plot_volcano <- function(de_tbl, cfg, title = NULL, use_adj = TRUE) {
  stopifnot(is.data.frame(de_tbl))
  if (!("logFC" %in% colnames(de_tbl))) stop("plot_volcano: de_tbl missing 'logFC' column.")
  
  p_col <- if (isTRUE(use_adj) && "adj.P.Val" %in% colnames(de_tbl)) "adj.P.Val" else "P.Value"
  if (!(p_col %in% colnames(de_tbl))) {
    stop("plot_volcano: de_tbl missing p-value column. Need 'P.Value' or 'adj.P.Val'.")
  }
  
  # thresholds from config (with safe defaults)
  p_cut <- cfg$de$p_cutoff %||% 0.1
  lin_fc_cut <- cfg$de$linear_fc_cutoff %||% 1.5
  log2fc_cut <- log2(lin_fc_cut)
  
  df <- de_tbl
  df$.p <- as.numeric(df[[p_col]])
  df$.logFC <- as.numeric(df[["logFC"]])
  # df$.neglog10p <- -log10(df$.p)
  
  # pass flag (robust to NAs)
  df$.pass <- !is.na(df$.p) & !is.na(df$.logFC) &
    (df$.p <= p_cut) & (abs(df$.logFC) >= log2fc_cut)
  
  df$.neglog10p <- -log10(pmax(df$.p, .Machine$double.xmin))
  
  ggplot2::ggplot(df, ggplot2::aes(x = .logFC, y = .neglog10p)) +
    ggplot2::geom_point(ggplot2::aes(color = .pass), size = 1.5, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "red"), guide = "none") +
    ggplot2::geom_vline(xintercept = c(-log2fc_cut, log2fc_cut), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = "dashed") +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.9), guide = "none") +
    ggplot2::labs(
      title = title %||% "Volcano plot",
      x = "log2 Fold Change",
      y = paste0("-log10(", p_col, ")")
    ) +
    ggplot2::theme_minimal()
}
