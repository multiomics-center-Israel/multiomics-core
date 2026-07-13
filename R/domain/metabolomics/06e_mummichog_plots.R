# R/domain/metabolomics/06e_mummichog_plots.R
#
# Presentation layer for the pinned mummichog v2 stage (06c/05b): a
# MetaboAnalyst-style pathway bubble plot and a clean results table for the
# metabolomics HTML report. These consume the single-run pathway table the stage
# already produces (mcg_pathwayanalysis_*.tsv, read via read_mummichog_pathways);
# they run NO analysis of their own.
#
# All builders are PURE — they return objects (a ggplot / a data.frame / a list)
# and never touch disk. The report module (06_mod_report) calls them and passes
# the results into the Rmd. Columns arrive verbatim from mummichog 2.7.0, so the
# p-value column is the literal `p-value` (hyphenated) — located tolerantly via
# .mmc_find_col() (06c).
#
# Dependencies (R): ggplot2, RColorBrewer (YlOrRd), optionally ggrepel (guarded).
# All already in renv.lock.


#' Build the mummichog pathway bubble plot
#'
#' MetaboAnalyst-style bubble plot of single-run mummichog pathway results:
#' enrichment ratio (x) vs significance (y), bubbles sized by pathway size and
#' coloured by -log10(p). A dashed line marks the significance cutoff and the
#' most significant pathways are labelled.
#'
#' @param pathways  A data.frame of mummichog pathway results (from
#'   `read_mummichog_pathways()`): columns `pathway`, `overlap_size`,
#'   `pathway_size`, and a p-value column (literal `p-value`).
#' @param title     Plot title.
#' @param subtitle  Plot subtitle (e.g. the contrast). Optional.
#' @param p_cutoff  Significance cutoff for the dashed reference line.
#' @param label_top Number of most-significant pathways to label (needs
#'   `ggrepel`; silently unlabelled if it is unavailable).
#' @return A ggplot object, or `NULL` when there is nothing to plot (so the
#'   report section hides gracefully).
plot_mummichog_bubble <- function(pathways, title, subtitle = NULL,
                                  p_cutoff = 0.05, label_top = 8) {
  if (is.null(pathways) || !is.data.frame(pathways) || nrow(pathways) == 0) {
    return(NULL)
  }
  p_col <- .mmc_find_col(pathways,
                         c("p-value", "p.value", "pvalue", "p_value", "P.Value"),
                         "^p[._-]?value$")
  needed <- c("pathway", "overlap_size", "pathway_size")
  if (is.null(p_col) || !all(needed %in% names(pathways))) return(NULL)

  df <- data.frame(
    pathway      = as.character(pathways$pathway),
    overlap_size = suppressWarnings(as.numeric(pathways$overlap_size)),
    pathway_size = suppressWarnings(as.numeric(pathways$pathway_size)),
    p            = suppressWarnings(as.numeric(pathways[[p_col]])),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$p) & is.finite(df$pathway_size) & df$pathway_size > 0, ,
           drop = FALSE]
  if (nrow(df) == 0) return(NULL)

  df$enrichment_ratio <- df$overlap_size / df$pathway_size
  # -log10 with a floor so an exact-zero p (never expected here) can't blow up.
  df$neg_log10_p <- -log10(pmax(df$p, .Machine$double.eps))

  p <- ggplot2::ggplot(
    df, ggplot2::aes(x = .data$enrichment_ratio, y = .data$neg_log10_p)) +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff),
                        linetype = "dashed", colour = "grey50") +
    ggplot2::annotate("text", x = Inf, y = -log10(p_cutoff),
                      label = sprintf("p = %g", p_cutoff),
                      hjust = 1.05, vjust = -0.5, size = 3, colour = "grey40") +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$pathway_size, fill = .data$neg_log10_p),
      shape = 21, colour = "grey20", stroke = 0.3, alpha = 0.9) +
    ggplot2::scale_fill_distiller(palette = "YlOrRd", direction = 1,
                                  name = "-log10(p)") +
    ggplot2::scale_size_continuous(range = c(2, 10), name = "pathway size") +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      x        = "Enrichment ratio  (overlap / pathway size)",
      y        = "-log10(p-value)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5))

  # Label the most significant pathways (smallest p). ggrepel is optional across
  # the repo, so guard it and simply skip labels when it is not installed.
  if (label_top > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    lab_df <- df[order(df$p), , drop = FALSE]
    lab_df <- lab_df[seq_len(min(label_top, nrow(lab_df))), , drop = FALSE]
    p <- p + ggrepel::geom_text_repel(
      data = lab_df,
      ggplot2::aes(label = .data$pathway),
      size = 3, max.overlaps = Inf, box.padding = 0.5,
      min.segment.length = 0, seed = 42)
  }
  p
}

#' Build a clean, sorted mummichog pathway results table
#'
#' Reshapes the raw mummichog pathway table into a compact, report-ready frame:
#' friendly column names, an added enrichment ratio, sorted by ascending
#' p-value. Intended for the report's `make_dt()`.
#'
#' @param pathways A data.frame of mummichog pathway results (see
#'   `plot_mummichog_bubble()`).
#' @return A data.frame with columns `Pathway`, `Overlap`, `Pathway size`,
#'   `Enrichment ratio`, `p.value`, sorted by `p.value`; or `NULL` when empty.
#'   The p-value column is named `p.value` (dot, not hyphen) so the report's
#'   `make_dt()` formatter renders it in scientific notation rather than rounding.
build_mummichog_pathway_table <- function(pathways) {
  if (is.null(pathways) || !is.data.frame(pathways) || nrow(pathways) == 0) {
    return(NULL)
  }
  p_col <- .mmc_find_col(pathways,
                         c("p-value", "p.value", "pvalue", "p_value", "P.Value"),
                         "^p[._-]?value$")
  if (is.null(p_col) || !all(c("pathway", "overlap_size", "pathway_size") %in%
                             names(pathways))) {
    return(NULL)
  }
  overlap      <- suppressWarnings(as.numeric(pathways$overlap_size))
  pathway_size <- suppressWarnings(as.numeric(pathways$pathway_size))
  pval         <- suppressWarnings(as.numeric(pathways[[p_col]]))

  out <- data.frame(
    check.names = FALSE,
    stringsAsFactors = FALSE,
    "Pathway"          = as.character(pathways$pathway),
    "Overlap"          = overlap,
    "Pathway size"     = pathway_size,
    "Enrichment ratio" = round(overlap / pathway_size, 3),
    "p.value"          = pval
  )
  out[order(out[["p.value"]]), , drop = FALSE]
}

#' Compose the mummichog report plot title and subtitle
#'
#' Derives a human-readable title/subtitle from the pipeline config and DE
#' results: `"Mummichog pathway analysis — <organism> (<model> model)"` with a
#' `"<contrast>, all features"` subtitle. The model tag is a heuristic short name
#' taken from the configured model (the leading token of a `model_ref` URL
#' basename, e.g. `cre_kegg_20260711.json` -> `cre`; a `model_json` basename;
#' else the built-in `human_mfn`). The contrast is the first DE table's name
#' (the pinned engine runs on a single contrast), prettified `A_vs_B` -> `A vs B`.
#'
#' @param config Full pipeline config.
#' @param de_res Optional DE results; when `contrast` is not given, the subtitle
#'   uses `names(de_res$de_tables)[1]`.
#' @param contrast Optional explicit contrast label for the subtitle (takes
#'   precedence over `de_res`). Used to title per-contrast mummichog sections.
#' @return A list with character scalars `title` and `subtitle` (subtitle may be
#'   `NULL` when no contrast is available).
mummichog_report_titles <- function(config, de_res = NULL, contrast = NULL) {
  mummi_cfg <- config$modes$metabolomics$enrichment$mummichog %||% list()

  organism <- config$modes$metabolomics$organism
  organism <- trimws(paste(unlist(organism), collapse = " "))

  model_tag <- if (!is.null(mummi_cfg$model_ref$url) &&
                   nzchar(mummi_cfg$model_ref$url)) {
    base <- sub("\\.json$", "", basename(mummi_cfg$model_ref$url),
                ignore.case = TRUE)
    sub("_.*$", "", base)                     # cre_kegg_20260711 -> cre
  } else if (!is.null(mummi_cfg$model_json) && nzchar(mummi_cfg$model_json)) {
    sub("\\.json$", "", basename(mummi_cfg$model_json), ignore.case = TRUE)
  } else {
    "human_mfn"
  }

  title <- if (nzchar(organism)) {
    sprintf("Mummichog pathway analysis — %s (%s model)", organism, model_tag)
  } else {
    sprintf("Mummichog pathway analysis (%s model)", model_tag)
  }

  if (is.null(contrast) && !is.null(de_res$de_tables) &&
      length(de_res$de_tables) > 0) {
    contrast <- names(de_res$de_tables)[1]
  }
  subtitle <- if (!is.null(contrast) && nzchar(contrast)) {
    paste0(gsub("_vs_", " vs ", contrast, fixed = TRUE), ", all features")
  } else {
    NULL
  }

  list(title = title, subtitle = subtitle)
}

#' Build per-contrast mummichog report sections
#'
#' Turns the per-contrast pathway tables (from `metab_mummichog_report_pathways`)
#' into ready-to-render report sections. Pure: builds a bubble plot + a sorted
#' table + a per-contrast title for each contrast that has a usable result, and
#' drops contrasts with no result. No file I/O — the report module saves the
#' standalone exports separately.
#'
#' @param pathways_by_contrast Named list keyed by contrast, each a mummichog
#'   pathway tibble (or NULL). Also tolerates a single tibble (one contrast).
#' @param config Full pipeline config (for titles + the p-value cutoff line).
#' @return A named list keyed by contrast, each `list(title, subtitle, plot,
#'   table)`; empty list when there is nothing to show.
build_mummichog_report_sections <- function(pathways_by_contrast, config) {
  if (is.null(pathways_by_contrast)) return(list())
  # Tolerate a bare single tibble (name it "contrast").
  if (is.data.frame(pathways_by_contrast)) {
    pathways_by_contrast <- list(contrast = pathways_by_contrast)
  }
  p_cut <- config$modes$metabolomics$enrichment$mummichog$p_cutoff %||% 0.05

  nms <- names(pathways_by_contrast)
  if (is.null(nms)) nms <- paste0("contrast_", seq_along(pathways_by_contrast))

  sections <- list()
  for (i in seq_along(pathways_by_contrast)) {
    contrast <- nms[[i]]
    pw       <- pathways_by_contrast[[i]]
    if (is.null(pw) || !is.data.frame(pw) || nrow(pw) == 0) next
    ttl  <- mummichog_report_titles(config, contrast = contrast)
    plot <- plot_mummichog_bubble(pw, title = ttl$title,
                                  subtitle = ttl$subtitle, p_cutoff = p_cut)
    tbl  <- build_mummichog_pathway_table(pw)
    if (is.null(plot)) next            # nothing plottable -> skip this contrast
    sections[[contrast]] <- list(title = ttl$title, subtitle = ttl$subtitle,
                                 plot = plot, table = tbl)
  }
  sections
}

#' Save the mummichog report plot and table as standalone files
#'
#' Writes the presentation exports for the pinned mummichog stage into the
#' `mummichog_pinned/` folder, next to the engine's result tree: the bubble plot
#' as PNG + PDF and the sorted pathway table as TSV + CSV. These are report
#' artefacts derived from the pinned results — the engine's own result files
#' (mcg_*.tsv) are left untouched.
#'
#' Each ggsave is guarded (matching the DE/enrichment modules): a device
#' failure warns and skips that file rather than aborting report generation, so
#' the return vector lists only the files actually written.
#'
#' @param plot   A ggplot from `plot_mummichog_bubble()`, or NULL.
#' @param table  A data.frame from `build_mummichog_pathway_table()`, or NULL.
#' @param out_dir Metabolomics mode output directory (metab_out_dir); files land
#'   in its `mummichog_pinned/` subdirectory.
#' @param contrast_label Optional contrast tag woven into the filenames
#'   (e.g. `"LL_vs_HL"`); sanitised to a safe token.
#' @return Character vector of the files written (may be empty).
save_mummichog_exports <- function(plot, table, out_dir, contrast_label = NULL) {
  if (is.null(plot) && is.null(table)) return(character(0))

  # Write next to the engine's result tree (same dir as mod_mummichog_pinned()),
  # so every mummichog artefact sits together under mummichog_pinned/. The stage
  # runs first and only wipes its v2/ subdir, so these top-level files persist.
  dest_dir <- file.path(out_dir, "mummichog_pinned")
  ensure_dir(dest_dir)
  suffix <- if (!is.null(contrast_label) && nzchar(contrast_label)) {
    paste0("_", gsub("[^A-Za-z0-9]+", "_", contrast_label))
  } else {
    ""
  }
  written <- character(0)

  if (!is.null(plot)) {
    for (f in c(file.path(dest_dir, paste0("mummichog_pathway_bubble", suffix, ".png")),
                file.path(dest_dir, paste0("mummichog_pathway_bubble", suffix, ".pdf")))) {
      ok <- tryCatch({
        ggplot2::ggsave(f, plot, width = 11, height = 8, dpi = 300)
        TRUE
      }, error = function(e) {
        warning("mummichog export failed for ", basename(f), ": ",
                conditionMessage(e))
        FALSE
      })
      if (isTRUE(ok) && file.exists(f)) written <- c(written, f)
    }
  }

  if (!is.null(table)) {
    tsv <- file.path(dest_dir, paste0("mummichog_pathway_table", suffix, ".tsv"))
    csv <- file.path(dest_dir, paste0("mummichog_pathway_table", suffix, ".csv"))
    readr::write_tsv(table, tsv)
    readr::write_csv(table, csv)
    written <- c(written, tsv, csv)
  }

  written
}
