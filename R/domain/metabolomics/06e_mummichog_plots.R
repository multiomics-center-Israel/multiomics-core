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

#' Find one contrast's DE table
#'
#' The ranking statistic GSEA needs lives in the DE result, keyed by the same
#' contrast label the mummichog stage ran under. Falls back to the single-table
#' shapes the stage itself accepts (05b), in the same order, so a config the
#' stage runs happily does not silently lose GSEA.
#'
#' Returns NULL rather than guessing when several tables are present and none
#' matches: picking one arbitrarily would rank a contrast on another contrast's
#' statistics. The caller reports that by name.
#'
#' @param de_res   DE result list, or NULL.
#' @param contrast Contrast label to look up.
#' @return A data.frame, or NULL when none can be resolved unambiguously.
#' @noRd
.mmc_de_table_for_contrast <- function(de_res, contrast) {
  if (is.null(de_res)) return(NULL)
  tabs <- de_res$de_tables
  if (!is.null(tabs) && length(tabs) > 0) {
    if (!is.null(contrast) && contrast %in% names(tabs)) return(tabs[[contrast]])
    if (length(tabs) == 1L) return(tabs[[1]])
    return(NULL)
  }
  if (is.data.frame(de_res$summary_df)) return(de_res$summary_df)
  if (is.data.frame(de_res)) return(de_res)
  NULL
}

#' Build per-contrast mummichog report sections
#'
#' Turns the per-contrast pathway tables (from `metab_mummichog_report_pathways`)
#' into ready-to-render report sections: a bubble plot + a sorted table + a
#' per-contrast title for each contrast that has a usable result, and, when the
#' mummichog output files and a readable metabolic model are supplied, that
#' contrast's supporting evidence. Reads files but writes none — the report
#' module saves the standalone exports separately.
#'
#' The evidence layer is additive and fail-soft: it recomputes no statistics,
#' and a contrast whose evidence cannot be built still renders its ORA plot and
#' table with `evidence = NULL`.
#'
#' @param pathways_by_contrast Named list keyed by contrast, each a mummichog
#'   pathway tibble (or NULL). Also tolerates a single tibble (one contrast).
#' @param config Full pipeline config (for titles + the p-value cutoff line).
#' @param files  Flat character vector of every mummichog output file across all
#'   contrasts (the `format = "file"` target's value). The evidence layer reads
#'   the EmpiricalCompound tables from it; omit it and every section's
#'   `evidence` is `NULL`.
#' @param row_data Feature annotations (`pre$row_data`), normalised to the
#'   generic contract for the agreement comparison. Omit it and every feature
#'   reads "Not assessed".
#' @param de_res DE result list (`metab_de_res`). Supplies both the per-contrast
#'   ranking statistic and `de_res$method`, which decides WHICH column that is —
#'   see `mmc_gsea_rank_metric()`. Omit it and GSEA is skipped for every
#'   contrast; the ORA and evidence outputs are unaffected.
#' @param model  Metabolic model index from `mmc_load_model_pathways()`; when
#'   `NULL` it is resolved from `config` once for all contrasts. mummichog
#'   exports no pathway -> compound membership, so without a readable model
#'   there is no evidence to build.
#' @return A named list keyed by contrast, each `list(title, subtitle, plot,
#'   table, gsea, gsea_plot, evidence, slug)`. `slug` is a filesystem-safe,
#'   de-duplicated token for the contrast (mirrors the engine's per-contrast
#'   directory naming) so the standalone exports never collide. `plot` remains
#'   the ORA bubble and keeps its meaning; `gsea` is the result list from
#'   `run_mummichog_gsea()` and `gsea_plot` its scatter, both `NULL` when GSEA
#'   did not run. `evidence` is the three-frame list from
#'   `build_mummichog_pathway_evidence()`, or `NULL`. Empty list when there is
#'   nothing to show.
build_mummichog_report_sections <- function(pathways_by_contrast, config,
                                            files = NULL, row_data = NULL,
                                            de_res = NULL, model = NULL) {
  if (is.null(pathways_by_contrast)) return(list())
  # Tolerate a bare single tibble (name it "contrast").
  if (is.data.frame(pathways_by_contrast)) {
    pathways_by_contrast <- list(contrast = pathways_by_contrast)
  }
  mummi_cfg <- config$modes$metabolomics$enrichment$mummichog %||% list()
  p_cut <- mummi_cfg$p_cutoff %||% 0.05
  # Independent switch, defaulting to ON so configs written before GSEA existed
  # keep working and gain it. The stage's own `enabled` still gates everything:
  # without it there are no files to read.
  gsea_enabled <- isTRUE(mummi_cfg$gsea_enabled %||% TRUE)

  nms <- names(pathways_by_contrast)
  if (is.null(nms)) nms <- paste0("contrast_", seq_along(pathways_by_contrast))

  # Resolved once for every contrast: the model's pathway membership (mummichog
  # exports none, so it comes from the model the stage ran on) and the dataset's
  # own annotations flattened to the generic contract.
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
    plot <- plot_mummichog_bubble(pw, title = ttl$title,
                                  subtitle = ttl$subtitle, p_cutoff = p_cut)
    tbl  <- build_mummichog_pathway_table(pw)
    if (is.null(plot)) next            # nothing plottable -> skip this contrast

    # -- supporting evidence (additive, fail-soft) ---------------------------
    # Needs this contrast's own output files AND a readable model; without
    # either there is nothing to trace, and a failure here must not cost the
    # contrast its ORA plot and table.
    cfiles   <- files_by_contrast[[contrast]]
    evidence <- if (is.null(cfiles) || is.null(model)) NULL else tryCatch(
      build_mummichog_pathway_evidence(pw, cfiles, model, annot),
      error = function(e) {
        warning("mummichog evidence: contrast '", contrast, "' failed: ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )

    # -- complementary GSEA (additive, fail-soft) ----------------------------
    # Same preconditions as the evidence layer, plus a DE table for the ranking
    # statistic. Three outcomes are kept apart deliberately: a MISSING input is
    # silent or named (nothing ran, nothing to report), a VALID EMPTY result is
    # the engine's own message (it ran and found nothing), and an EXECUTION
    # ERROR warns. Conflating them would let an environment failure read as a
    # finding about the data.
    gsea <- NULL
    if (gsea_enabled && !is.null(cfiles) && !is.null(model)) {
      de_tbl <- .mmc_de_table_for_contrast(de_res, contrast)
      if (is.null(de_tbl)) {
        # Not silent: the stage ran this contrast, so losing GSEA here is a
        # mapping failure the reader should see, not "GSEA not applicable".
        message("mummichog GSEA: contrast '", contrast, "' has no matching DE ",
                "table — skipping GSEA for it (ORA and evidence are unaffected)")
      } else {
        gsea <- tryCatch(
          run_mummichog_gsea(
            files     = cfiles,
            de_table  = de_tbl,
            model     = model,
            contrast  = contrast,
            # THE critical argument: without it the engine will not trust the
            # generic `statistic` column and every limma run silently ranks on
            # logFC. See mmc_gsea_rank_metric().
            #
            # Guarded on the shape: in the bare-data.frame fallback above,
            # `de_res$method` would be a COLUMN if one happened to be named
            # `method`, and the engine would take its first value as the DE
            # method. A legacy table carries no method, so NULL is both the
            # honest answer and the safe one -- the engine then ranks on logFC
            # and says why.
            de_method = if (is.data.frame(de_res)) NULL else de_res$method,
            # MetaboAnalyst-equivalent defaults unless the config overrides
            # them, in which case run_mummichog_gsea() records the deviation.
            n_perm   = mummi_cfg$gsea_permutations %||% .MMC_GSEA_DEFAULTS$n_perm,
            seed     = mummi_cfg$gsea_seed         %||% .MMC_GSEA_DEFAULTS$seed,
            min_size = mummi_cfg$gsea_min_size     %||% .MMC_GSEA_DEFAULTS$min_size,
            max_size = mummi_cfg$gsea_max_size     %||% .MMC_GSEA_DEFAULTS$max_size
          ),
          error = function(e) {
            warning("mummichog GSEA: contrast '", contrast, "' failed: ",
                    conditionMessage(e), call. = FALSE)
            NULL
          }
        )
      }
    }
    # Same cutoff the ORA plot and the report's figure legends use: the legend
    # states the value from config, so a scatter drawn at the function default
    # would put its dashed line somewhere the caption does not claim.
    gsea_plot <- if (is.null(gsea)) NULL else
      plot_mummichog_gsea_scatter(gsea$table, title = ttl$title,
                                  subtitle = ttl$subtitle, p_cutoff = p_cut)

    sections[[contrast]] <- list(title = ttl$title, subtitle = ttl$subtitle,
                                 plot = plot, table = tbl,
                                 gsea = gsea, gsea_plot = gsea_plot,
                                 evidence = evidence)
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
#' @param evidence Optional three-frame list from
#'   `build_mummichog_pathway_evidence()`. Each frame is written as TSV + CSV
#'   beside the pathway table; the HTML report caps its drill-down, so these are
#'   the complete tables. This contrast's six evidence files are removed first,
#'   by name, so a run that comes back without evidence (or without one frame)
#'   cannot leave the previous run's tables behind. Exports of other kinds, of
#'   other contrasts, and anything else in the directory are left alone.
#' @param gsea Optional GSEA result from `run_mummichog_gsea()`; its table is
#'   written as `mummichog_gsea_table_<contrast>.{tsv,csv}`.
#' @param gsea_plot Optional scatter from `plot_mummichog_gsea_scatter()`;
#'   written as `mummichog_gsea_scatter_<contrast>.{png,pdf}`. This contrast's
#'   four GSEA files are cleared by name first, for the same reason the evidence
#'   ones are: GSEA can legitimately disappear from a rerun.
#' @return Character vector of the files written (may be empty).
save_mummichog_exports <- function(plot, table, out_dir, contrast_label = NULL,
                                   gsea = NULL, gsea_plot = NULL,
                                   evidence = NULL) {
  # Write next to the engine's result tree (same dir as mod_mummichog_pinned()),
  # so every mummichog artefact sits together under mummichog_pinned/. The stage
  # runs first and only wipes its v2/ subdir, so these top-level files persist.
  dest_dir <- file.path(out_dir, "mummichog_pinned")
  suffix <- if (!is.null(contrast_label) && nzchar(contrast_label)) {
    paste0("_", gsub("[^A-Za-z0-9]+", "_", contrast_label))
  } else {
    ""
  }

  # Clear this contrast's evidence exports before writing whatever this run
  # produced. A rerun can legitimately come back without evidence -- no readable
  # model, a partial mummichog run salvaged by .mmc_check_run_outputs(), or the
  # fail-soft branch in build_mummichog_report_sections() -- and a single frame
  # can be missing from an otherwise good set. {targets} stops tracking paths a
  # later run omits but never deletes them, so without this the report correctly
  # omits the tab while last run's tables stay on disk reading as current.
  #
  # Ahead of the nothing-to-save return, so a contrast that comes back empty is
  # covered too, and unconditional, so the per-frame case is. Only these six
  # names, never a pattern: the pathway and plot exports, other contrasts, and
  # anything else in the directory are untouched.
  evidence_exports <- file.path(
    dest_dir,
    paste0(rep(c("mummichog_evidence_pathways",
                 "mummichog_evidence_empirical_compounds",
                 "mummichog_evidence_features"), each = 2),
           suffix, c(".tsv", ".csv")))
  # Same reasoning, same by-name discipline, for this contrast's GSEA exports:
  # GSEA can legitimately vanish from a rerun (switched off, no DE table for the
  # contrast, no pathway through the size filter) and the stale table would read
  # as current. Scoped to these four names only -- the general stale-output
  # problem is #299 and is not touched here.
  gsea_exports <- c(
    file.path(dest_dir, paste0("mummichog_gsea_table", suffix, c(".tsv", ".csv"))),
    file.path(dest_dir, paste0("mummichog_gsea_scatter", suffix, c(".png", ".pdf")))
  )
  unlink(gsea_exports[file.exists(gsea_exports)])

  if (is.null(plot) && is.null(table) && is.null(evidence) &&
      is.null(gsea) && is.null(gsea_plot)) {
    return(character(0))
  }

  ensure_dir(dest_dir)
  written <- character(0)

  # One guarded ggsave for both plots: a device failure warns and skips that
  # file rather than aborting the report, and the GSEA scatter gets the same
  # treatment the bubble already had instead of a second copy of it.
  save_plot <- function(p, stem) {
    if (is.null(p)) return(character(0))
    out <- character(0)
    for (f in file.path(dest_dir, paste0(stem, suffix, c(".png", ".pdf")))) {
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
  written <- c(written, save_plot(plot, "mummichog_pathway_bubble"))
  written <- c(written, save_plot(gsea_plot, "mummichog_gsea_scatter"))

  # NULL only: a zero-row frame is still a result, and the pathway table's
  # existing behaviour (written whenever it is not NULL) must not change.
  save_tsv_csv <- function(df, stem) {
    if (is.null(df) || !is.data.frame(df)) return(character(0))
    tsv <- file.path(dest_dir, paste0(stem, suffix, ".tsv"))
    csv <- file.path(dest_dir, paste0(stem, suffix, ".csv"))
    readr::write_tsv(df, tsv)
    readr::write_csv(df, csv)
    c(tsv, csv)
  }

  if (!is.null(table)) {
    written <- c(written, save_tsv_csv(table, "mummichog_pathway_table"))
  }

  # The report caps its GSEA table too; this is the complete one.
  if (!is.null(gsea)) {
    written <- c(written, save_tsv_csv(gsea$table, "mummichog_gsea_table"))
  }

  # The report's drill-down is capped; these carry every row.
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
