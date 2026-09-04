# R/domain/metabolomics/06e_mummichog_plots.R
#
# Presentation layer for the pinned mummichog v2 stage (06c/05b). It assembles,
# per contrast, the three views the report shows:
#
#   * the MetaboAnalyst-style GSEA summary scatter (06g) — the section's pathway
#     plot, replacing the earlier ORA bubble plot, which is kept here as the
#     fallback for when GSEA cannot run (no readable metabolic model, no fgsea,
#     no signed statistic);
#   * the mummichog ORA pathway table, unchanged;
#   * the pathway supporting-evidence drill-down (06f).
#
# These builders run NO analysis of their own beyond calling the GSEA layer;
# mummichog's ORA statistics are passed through verbatim. Columns arrive as
# mummichog 2.7.0 wrote them, so the ORA p-value column is the literal `p-value`
# (hyphenated) — located tolerantly via .mmc_find_col() (06c).
#
# All builders are PURE — they return objects (a ggplot / a data.frame / a list)
# and never touch disk, except save_mummichog_exports() which exists to write.
# The report module (06_mod_report) calls them and passes the results to the Rmd.
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
#' into ready-to-render report sections: for each contrast with a usable result,
#' the ORA pathway table (unchanged), the complementary GSEA analysis + its
#' MetaboAnalyst-style scatter (06g), and the pathway supporting-evidence
#' drill-down (06f). Contrasts with no result are dropped.
#'
#' The section's `plot` is the GSEA scatter when GSEA could run, and the ORA
#' bubble plot otherwise — so the report always has a pathway plot even without
#' a readable metabolic model, fgsea, or a signed DE statistic. Which one it is
#' shows in `plot_kind`. The ORA table and its evidence are never replaced by
#' the GSEA view; they are complementary answers to different questions.
#'
#' Every added layer is optional and fail-soft: a GSEA or evidence failure warns
#' and leaves that part `NULL` rather than losing the contrast's ORA results.
#'
#' @param pathways_by_contrast Named list keyed by contrast, each a mummichog
#'   pathway tibble (or NULL). Also tolerates a single tibble (one contrast).
#' @param config Full pipeline config (for titles + the p-value cutoff line).
#' @param files  Flat character vector of every mummichog output file across all
#'   contrasts (the `format = "file"` target's value). Needed for the EC-level
#'   tables that GSEA and the evidence layer read; omit it and the section falls
#'   back to ORA-only.
#' @param de_res DE results; the per-contrast DE table supplies GSEA's signed
#'   ranking statistic (moderated `t` when present — see
#'   `mmc_gsea_rank_metric()`). Omit it and GSEA is skipped.
#' @param row_data Feature annotations (`pre$row_data`) for the original
#'   annotations the evidence layer compares against. Omit it and every feature
#'   reads "Not assessed".
#' @param model  Metabolic model index from `mmc_load_model_pathways()`; when
#'   `NULL` it is resolved from `config` once for all contrasts.
#' @return A named list keyed by contrast, each `list(title, subtitle, plot,
#'   plot_kind, ora_plot, table, gsea, gsea_plot, evidence, slug)`. `plot` is the
#'   one the report renders and `plot_kind` says which it is; `ora_plot` and
#'   `gsea_plot` are always the individual plots (the latter `NULL` when GSEA did
#'   not run), so the exports can write both. `slug` is a filesystem-safe,
#'   de-duplicated token for the contrast (mirrors the engine's per-contrast
#'   directory naming) so the standalone exports never collide. Empty list when
#'   there is nothing to show.
build_mummichog_report_sections <- function(pathways_by_contrast, config,
                                            files = NULL, de_res = NULL,
                                            row_data = NULL, model = NULL) {
  if (is.null(pathways_by_contrast)) return(list())
  # Tolerate a bare single tibble (name it "contrast").
  if (is.data.frame(pathways_by_contrast)) {
    pathways_by_contrast <- list(contrast = pathways_by_contrast)
  }
  mummi_cfg <- config$modes$metabolomics$enrichment$mummichog %||% list()
  p_cut <- mummi_cfg$p_cutoff %||% 0.05

  nms <- names(pathways_by_contrast)
  if (is.null(nms)) nms <- paste0("contrast_", seq_along(pathways_by_contrast))

  # Shared, resolved once for every contrast: the model's pathway membership
  # (mummichog exports none, so it comes from the model the stage ran on) and the
  # dataset's own annotations flattened to the generic contract.
  if (is.null(model) && !is.null(files) && length(files) > 0) {
    model <- mmc_load_model_pathways(config)
  }
  # The HMDB -> KEGG mapping is optional; a config without a resolvable paths
  # block must not cost us the whole evidence layer.
  mapping_file <- tryCatch(
    resolve_input_path(config, config$modes$metabolomics$enrichment$mapping_file),
    error = function(e) NULL
  )
  annot <- normalize_metab_annotation(row_data, mapping_file = mapping_file)
  files_by_contrast <- if (is.null(files)) list() else
    group_mummichog_files_by_contrast(files)

  sections <- list()
  for (i in seq_along(pathways_by_contrast)) {
    contrast <- nms[[i]]
    pw       <- pathways_by_contrast[[i]]
    if (is.null(pw) || !is.data.frame(pw) || nrow(pw) == 0) next
    ttl  <- mummichog_report_titles(config, contrast = contrast)
    ora_plot <- plot_mummichog_bubble(pw, title = ttl$title,
                                      subtitle = ttl$subtitle, p_cutoff = p_cut)
    tbl  <- build_mummichog_pathway_table(pw)
    if (is.null(ora_plot)) next        # nothing plottable -> skip this contrast

    cfiles <- files_by_contrast[[contrast]]

    # -- complementary GSEA (additive; never touches the ORA results) --------
    gsea <- NULL
    if (!is.null(cfiles) && !is.null(model)) {
      de_tbl <- .mmc_de_table_for_contrast(de_res, contrast)
      if (!is.null(de_tbl)) {
        gsea <- tryCatch(
          run_mummichog_gsea(
            files    = cfiles,
            de_table = de_tbl,
            model    = model,
            contrast = contrast,
            n_perm   = mummi_cfg$gsea_permutations %||% 1000,
            seed     = mummi_cfg$gsea_seed %||% 42
          ),
          error = function(e) {
            warning("mummichog GSEA: contrast '", contrast, "' failed: ",
                    conditionMessage(e), call. = FALSE)
            NULL
          }
        )
      }
    }
    gsea_plot <- if (is.null(gsea)) NULL else
      plot_mummichog_gsea_scatter(
        gsea$table,
        title = sub("^Mummichog pathway analysis",
                    "Mummichog GSEA (peaks to pathways)", ttl$title),
        subtitle = ttl$subtitle, p_cutoff = p_cut)

    # -- supporting evidence (fail-soft) -------------------------------------
    evidence <- if (is.null(cfiles) || is.null(model)) NULL else tryCatch(
      build_mummichog_pathway_evidence(pw, cfiles, model, annot,
                                       p_cutoff = p_cut),
      error = function(e) {
        warning("mummichog evidence: contrast '", contrast, "' failed: ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )

    sections[[contrast]] <- list(
      title     = ttl$title,
      subtitle  = ttl$subtitle,
      plot      = gsea_plot %||% ora_plot,
      plot_kind = if (is.null(gsea_plot)) "ora_bubble" else "gsea_scatter",
      ora_plot  = ora_plot,
      table     = tbl,
      gsea      = gsea,
      gsea_plot = gsea_plot,
      evidence  = evidence
    )
  }

  # Assign a de-duplicated, filesystem-safe slug per section (same sanitise +
  # make.unique idiom the engine uses for contrast directories), so two labels
  # that collapse to the same token get distinct export filenames instead of
  # overwriting each other. Used by save_mummichog_exports() for the suffix.
  if (length(sections) > 0) {
    slugs <- make.unique(gsub("[^A-Za-z0-9]+", "_", names(sections)), sep = "_")
    for (i in seq_along(sections)) sections[[i]]$slug <- slugs[[i]]
  }
  sections
}

#' Locate one contrast's DE table inside a DE results object
#'
#' Mirrors how the stage itself resolves contrasts (05b): the named
#' `de_tables` entry first, then the single-table fallbacks a non-standard
#' `de_res` can take. Returns `NULL` rather than guessing when the contrast
#' cannot be identified — GSEA is then skipped for it, which is preferable to
#' ranking one contrast's ECs with another contrast's statistics.
#'
#' @param de_res   DE results (or a bare data.frame).
#' @param contrast Contrast label to look up.
#' @return A data.frame, or `NULL`.
#' @noRd
.mmc_de_table_for_contrast <- function(de_res, contrast) {
  if (is.null(de_res)) return(NULL)
  tabs <- de_res$de_tables
  if (!is.null(tabs) && length(tabs) > 0) {
    if (!is.null(contrast) && contrast %in% names(tabs)) return(tabs[[contrast]])
    # A single unnamed contrast is unambiguous; anything else is not.
    if (length(tabs) == 1L) return(tabs[[1]])
    return(NULL)
  }
  if (is.data.frame(de_res)) return(de_res)
  NULL
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
#' @param gsea   Optional GSEA result from `run_mummichog_gsea()`; its table is
#'   written as `mummichog_gsea_table_<contrast>.{tsv,csv}`.
#' @param gsea_plot Optional GSEA scatter from `plot_mummichog_gsea_scatter()`;
#'   written as `mummichog_gsea_scatter_<contrast>.{png,pdf}`.
#' @param evidence Optional evidence object from
#'   `build_mummichog_pathway_evidence()`; its three frames are written as
#'   `mummichog_evidence_{pathways,empirical_compounds,features}_<contrast>.tsv`.
#' @return Character vector of the files written (may be empty).
save_mummichog_exports <- function(plot, table, out_dir, contrast_label = NULL,
                                   gsea = NULL, gsea_plot = NULL,
                                   evidence = NULL) {
  if (is.null(plot) && is.null(table) && is.null(gsea) &&
      is.null(gsea_plot) && is.null(evidence)) {
    return(character(0))
  }

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

  # Each ggsave is guarded: a device failure warns and skips that file rather
  # than aborting the whole report.
  save_plot <- function(p, stem) {
    out <- character(0)
    if (is.null(p)) return(out)
    for (f in c(file.path(dest_dir, paste0(stem, suffix, ".png")),
                file.path(dest_dir, paste0(stem, suffix, ".pdf")))) {
      ok <- tryCatch({
        ggplot2::ggsave(f, p, width = 11, height = 8, dpi = 300)
        TRUE
      }, error = function(e) {
        warning("mummichog export failed for ", basename(f), ": ",
                conditionMessage(e))
        FALSE
      })
      if (isTRUE(ok) && file.exists(f)) out <- c(out, f)
    }
    out
  }
  save_tsv_csv <- function(df, stem) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(character(0))
    tsv <- file.path(dest_dir, paste0(stem, suffix, ".tsv"))
    csv <- file.path(dest_dir, paste0(stem, suffix, ".csv"))
    readr::write_tsv(df, tsv)
    readr::write_csv(df, csv)
    c(tsv, csv)
  }

  written <- c(written, save_plot(plot, "mummichog_pathway_bubble"))
  written <- c(written, save_tsv_csv(table, "mummichog_pathway_table"))

  # GSEA is the complementary analysis, exported alongside — never instead of —
  # the ORA table above.
  written <- c(written, save_plot(gsea_plot, "mummichog_gsea_scatter"))
  if (!is.null(gsea)) {
    written <- c(written, save_tsv_csv(gsea$table, "mummichog_gsea_table"))
  }

  if (!is.null(evidence)) {
    written <- c(
      written,
      save_tsv_csv(evidence$pathway_summary, "mummichog_evidence_pathways"),
      save_tsv_csv(evidence$ec_table, "mummichog_evidence_empirical_compounds"),
      save_tsv_csv(evidence$feature_table, "mummichog_evidence_features")
    )
  }

  written
}
