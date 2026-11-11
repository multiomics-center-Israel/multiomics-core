# R/modules/metabolomics/05_mod_enrichment.R
#
# Metabolomics enrichment module: runs QEA, ssGSEA, ORA, and GSEA,
# generates plots, and saves result tables.
#
# Reuses domain: run_metabolomics_qea, run_metabolomics_ssgsea,
#   run_metabolomics_ora, run_metabolomics_gsea,
#   plot_enrichment_barplot, plot_ssgsea_boxplots,
#   plot_gsea_nes_barplot, plot_ora_dotplot


#' Metabolomics enrichment module
#'
#' Runs up to four enrichment methods (each independently optional):
#'   QEA (globaltest), ssGSEA (GSVA), ORA (Fisher), GSEA (fgsea).
#'
#' Each result unit is written to its own directory under
#' \code{<out_dir>/Enrichment/} (see \code{06p_enrichment_layout.R}), and
#' \code{Enrichment/enrichment_manifest.csv} lists the units this run wrote.
#' Once the run has written its results it removes the stage's own files that
#' an earlier run left and this one did not write. A disabled run, a run that
#' fails and a run that produces no result remove nothing.
#' The first-contrast ORA/GSEA copies under the flat legacy names are still
#' written (\code{.save_legacy_enrichment()}).
#'
#' @param pre     List from preprocess_metabolomics().
#' @param de_res  DE results (needed by ORA and GSEA; NULL if unavailable).
#' @param config  Full config.
#' @param out_dir Output directory for this mode.
#' @param inputs  Loaded inputs (uses \code{$contrasts} for the per-contrast
#'   QEA/ssGSEA runs); NULL skips those.
#' @return list(qea, ssgsea, ora, gsea, qea_by_contrast, ssgsea_by_contrast,
#'   ora_by_contrast, gsea_by_contrast, plots, files).
mod_metabolomics_enrichment <- function(pre, de_res = NULL, config, out_dir,
                                        inputs = NULL) {
    assert_pre_contract(pre, stage = "metabolomics")

    enr_cfg <- config$modes$metabolomics$enrichment %||% list()

    if (!isTRUE(enr_cfg$run_enrichment)) {
        message("metabolomics enrichment: disabled — skipping")
        return(NULL)
    }

    dirs    <- create_legacy_output_dirs(out_dir)
    out_enr <- dirs$enrichment
    out_ds  <- dirs$datasets

    # The manifest is this run's record of what it wrote, and it is written
    # last. Removing the previous one first means a run that stops part-way
    # leaves no manifest, instead of an earlier run's manifest beside new files.
    manifest_file <- file.path(out_enr, .METAB_ENRICHMENT_MANIFEST_FILE)
    stale_manifest <- c(manifest_file, paste0(manifest_file, ".tmp"))
    unlink(stale_manifest[file.exists(stale_manifest)])

    contrast_names <- if (!is.null(de_res) && !is.null(de_res$de_tables)) {
        names(de_res$de_tables)
    } else {
        character(0)
    }
    contrast_dirs <- metab_enrichment_contrast_dirs(
        contrast_names, .sanitize_contrast(contrast_names))

    plots <- list()
    files <- character(0)
    units <- list()

    # ---- QEA ----
    qea_res <- tryCatch(
        run_metabolomics_qea(pre, config),
        error = function(e) {
            warning("metabolomics QEA failed: ", e$message)
            NULL
        }
    )

    if (!is.null(qea_res) && !is.null(qea_res$table)) {
        qea_rel <- metab_enrichment_unit_rel_dir("qea")
        qea_dir <- file.path(out_enr, qea_rel)
        f_qea <- save_tsv(qea_res$table, qea_dir, "results.tsv")
        files <- c(files, f_qea)
        units[[length(units) + 1L]] <- list(method = "qea", contrast = NA_character_,
                                            rel_dir = qea_rel, table = qea_res$table)

        # QEA barplot
        f_qea_plot <- file.path(qea_dir, "barplot.png")
        p_qea <- tryCatch({
            pq <- plot_enrichment_barplot(qea_res$table, top_n = 20,
                                           title = "Pathway Enrichment (QEA)")
            if (!is.null(pq)) {
                ggplot2::ggsave(f_qea_plot, pq, width = 12, height = 8, dpi = 300)
            }
            pq
        }, error = function(e) {
            warning("QEA barplot failed: ", e$message)
            NULL
        })
        if (!is.null(p_qea)) {
            files <- c(files, f_qea_plot)
            plots$qea_barplot <- p_qea
        }

        # QEA lollipop plot
        f_qea_loll <- file.path(qea_dir, "lollipop.png")
        p_qea_loll <- tryCatch({
            pl <- plot_qea_lollipop(qea_res$table, top_n = 20)
            if (!is.null(pl)) {
                ggplot2::ggsave(f_qea_loll, pl, width = 12, height = 8, dpi = 300)
            }
            pl
        }, error = function(e) NULL)
        if (!is.null(p_qea_loll)) {
            files <- c(files, f_qea_loll)
            plots$qea_lollipop <- p_qea_loll
        }

        # Per-library barplots
        if (!is.null(qea_res$table$library)) {
            for (lib_name in unique(qea_res$table$library)) {
                lib_df <- qea_res$table[qea_res$table$library == lib_name, ]
                lib_label <- toupper(gsub("_pathway$", "", lib_name))
                f_lib <- file.path(qea_dir, paste0(
                    .METAB_ENRICHMENT_LIBRARY_PLOT_PREFIX, tolower(lib_label), ".png"))
                p_lib <- tryCatch({
                    pl <- plot_enrichment_barplot(lib_df, top_n = 20,
                                                  title = paste0(lib_label, " Enrichment"))
                    if (!is.null(pl)) {
                        ggplot2::ggsave(f_lib, pl, width = 12, height = 8, dpi = 300)
                    }
                    pl
                }, error = function(e) NULL)
                if (!is.null(p_lib)) {
                    files <- c(files, f_lib)
                    plots[[paste0("qea_", tolower(lib_label))]] <- p_lib
                }
            }
        }
    }

    # ---- ssGSEA ----
    ssgsea_res <- tryCatch(
        run_metabolomics_ssgsea(pre, config),
        error = function(e) {
            warning("metabolomics ssGSEA failed: ", e$message)
            NULL
        }
    )

    if (!is.null(ssgsea_res) && !is.null(ssgsea_res$table)) {
        ssgsea_rel <- metab_enrichment_unit_rel_dir("ssgsea")
        ssgsea_dir <- file.path(out_enr, ssgsea_rel)
        f_ssgsea <- save_tsv(ssgsea_res$table, ssgsea_dir, "results.tsv")
        files <- c(files, f_ssgsea)
        units[[length(units) + 1L]] <- list(method = "ssgsea", contrast = NA_character_,
                                            rel_dir = ssgsea_rel,
                                            table = ssgsea_res$table)

        # ssGSEA barplot
        f_ssgsea_plot <- file.path(ssgsea_dir, "barplot.png")
        p_ssgsea <- tryCatch({
            ps <- plot_enrichment_barplot(ssgsea_res$table, top_n = 20,
                                           title = "ssGSEA Pathway Enrichment")
            if (!is.null(ps)) {
                ggplot2::ggsave(f_ssgsea_plot, ps, width = 12, height = 8, dpi = 300)
            }
            ps
        }, error = function(e) {
            warning("ssGSEA barplot failed: ", e$message)
            NULL
        })
        if (!is.null(p_ssgsea)) {
            files <- c(files, f_ssgsea_plot)
            plots$ssgsea_barplot <- p_ssgsea
        }

        # ssGSEA boxplots (if Wilcoxon results available)
        if (!is.null(ssgsea_res$scores) &&
            any(ssgsea_res$table$significant, na.rm = TRUE)) {
            sample_col <- config$modes$metabolomics$effects$samples %||% "sample_id"
            cond_col <- enr_cfg$condition_column %||%
                        config$modes$metabolomics$de$condition_column %||%
                        config$modes$metabolomics$effects$color %||% "sample_type"
            conditions <- factor(pre$meta[[cond_col]][
                match(colnames(ssgsea_res$scores), pre$meta[[sample_col]])
            ])

            f_box <- file.path(ssgsea_dir, "boxplots.png")
            p_box <- tryCatch({
                n_sig <- sum(ssgsea_res$table$significant, na.rm = TRUE)
                plot_h <- max(6, ceiling(min(n_sig, 12) / 3) * 3.5)
                pb <- plot_ssgsea_boxplots(ssgsea_res$scores, conditions,
                                            ssgsea_res$table)
                if (!is.null(pb)) {
                    ggplot2::ggsave(f_box, pb, width = 12, height = plot_h,
                                     dpi = 300)
                }
                pb
            }, error = function(e) {
                warning("ssGSEA boxplots failed: ", e$message)
                NULL
            })
            if (!is.null(p_box)) {
                files <- c(files, f_box)
                plots$ssgsea_boxplots <- p_box
            }
        }

        # ssGSEA lollipop plot
        f_ssgsea_loll <- file.path(ssgsea_dir, "lollipop.png")
        p_ssgsea_loll <- tryCatch({
            pl <- plot_ssgsea_lollipop(ssgsea_res$table, top_n = 20)
            if (!is.null(pl)) {
                ggplot2::ggsave(f_ssgsea_loll, pl, width = 12, height = 8, dpi = 300)
            }
            pl
        }, error = function(e) NULL)
        if (!is.null(p_ssgsea_loll)) {
            files <- c(files, f_ssgsea_loll)
            plots$ssgsea_lollipop <- p_ssgsea_loll
        }
    }

    # ---- ORA (per contrast) ----
    # ORA is DE-driven, so run it once per contrast in de_res$de_tables. The
    # first contrast is also exposed as `ora` + plots$ora_* and written under the
    # legacy (non-suffixed) file names for backward compatibility (Shiny payload,
    # PowerPoint, scripts/build_per_omics_pptx.py); the full set is returned as
    # ora_by_contrast for the report.
    ora_by_contrast <- list()
    if (!is.null(de_res) && !is.null(de_res$de_tables) &&
        length(de_res$de_tables) > 0) {
        for (ct in names(de_res$de_tables)) {
            rel <- metab_enrichment_unit_rel_dir("ora", contrast_dirs[[ct]])
            one <- .mod_metab_ora_contrast(pre, de_res, config, ct,
                                           file.path(out_enr, rel))
            if (!is.null(one)) {
                ora_by_contrast[[ct]] <- one
                units[[length(units) + 1L]] <- list(method = "ora", contrast = ct,
                                                    rel_dir = rel,
                                                    table = one$res$table)
            }
        }
    }
    ora_res <- NULL
    if (length(ora_by_contrast) > 0) {
        first <- ora_by_contrast[[1]]
        ora_res <- first$res
        if (!is.null(first$plots$dotplot))  plots$ora_dotplot  <- first$plots$dotplot
        if (!is.null(first$plots$lollipop)) plots$ora_lollipop <- first$plots$lollipop
        if (!is.null(first$plots$barplot))  plots$ora_barplot  <- first$plots$barplot
        for (nm in names(ora_by_contrast)) files <- c(files, ora_by_contrast[[nm]]$files)
        files <- c(files, .save_legacy_enrichment("ora", ora_res, first$plots,
                                                   out_enr, out_ds))
    }

    # ---- GSEA (per contrast) ----
    # Same per-contrast treatment as ORA (see comment above).
    gsea_by_contrast <- list()
    if (!is.null(de_res) && !is.null(de_res$de_tables) &&
        length(de_res$de_tables) > 0) {
        for (ct in names(de_res$de_tables)) {
            rel <- metab_enrichment_unit_rel_dir("gsea", contrast_dirs[[ct]])
            one <- .mod_metab_gsea_contrast(pre, de_res, config, ct,
                                            file.path(out_enr, rel))
            if (!is.null(one)) {
                gsea_by_contrast[[ct]] <- one
                units[[length(units) + 1L]] <- list(method = "gsea", contrast = ct,
                                                    rel_dir = rel,
                                                    table = one$res$table)
            }
        }
    }
    gsea_res <- NULL
    if (length(gsea_by_contrast) > 0) {
        first <- gsea_by_contrast[[1]]
        gsea_res <- first$res
        if (!is.null(first$plots$nes_barplot)) plots$gsea_nes_barplot <- first$plots$nes_barplot
        if (!is.null(first$plots$lollipop))    plots$gsea_lollipop    <- first$plots$lollipop
        if (!is.null(first$plots$barplot))     plots$gsea_barplot     <- first$plots$barplot
        for (nm in names(gsea_by_contrast)) files <- c(files, gsea_by_contrast[[nm]]$files)
        files <- c(files, .save_legacy_enrichment("gsea", gsea_res, first$plots,
                                                   out_enr, out_ds))
    }

    # ---- Per-contrast QEA & ssGSEA ----
    # QEA/ssGSEA are group-based (not DE-driven), so the global runs above use
    # the whole condition column (and ssGSEA yields no stats when there are more
    # than two groups). To get per-contrast results we subset the samples to each
    # contrast's two groups — read from inputs$contrasts (Factor + Numerator/
    # Denominator) — and re-run each method on that subset. The global qea/ssgsea
    # are kept for backward compatibility; the per-contrast sets feed the report.
    qea_by_contrast    <- list()
    ssgsea_by_contrast <- list()
    cond_col <- enr_cfg$condition_column %||%
                config$modes$metabolomics$de$condition_column %||%
                config$modes$metabolomics$effects$color %||% "sample_type"
    sample_col <- config$modes$metabolomics$effects$samples %||% "sample_id"
    if (!is.null(de_res) && !is.null(de_res$de_tables) &&
        length(de_res$de_tables) > 0 && !is.null(inputs)) {
        for (ct in names(de_res$de_tables)) {
            grp <- .resolve_contrast_groups(inputs, ct, cond_col)
            if (is.null(grp)) {
                warning("metabolomics QEA/ssGSEA: could not resolve the two groups ",
                        "for contrast '", ct, "' from inputs$contrasts — skipping.")
                next
            }
            pre_sub <- .subset_pre_to_groups(pre, grp$group_col, grp$groups, sample_col)
            if (is.null(pre_sub)) {
                warning("metabolomics QEA/ssGSEA: too few samples for contrast '",
                        ct, "' (groups ", paste(grp$groups, collapse = " vs "),
                        ") — skipping.")
                next
            }
            q_rel <- metab_enrichment_unit_rel_dir("qea", contrast_dirs[[ct]])
            q <- .mod_metab_qea_contrast(pre_sub, config, ct,
                                         file.path(out_enr, q_rel))
            if (!is.null(q)) {
                qea_by_contrast[[ct]] <- q; files <- c(files, q$files)
                units[[length(units) + 1L]] <- list(method = "qea", contrast = ct,
                                                    rel_dir = q_rel,
                                                    table = q$res$table)
            }
            s_rel <- metab_enrichment_unit_rel_dir("ssgsea", contrast_dirs[[ct]])
            s <- .mod_metab_ssgsea_contrast(pre_sub, config, ct,
                                            file.path(out_enr, s_rel),
                                            groups = grp$groups,
                                            group_col = grp$group_col)
            if (!is.null(s)) {
                ssgsea_by_contrast[[ct]] <- s; files <- c(files, s$files)
                units[[length(units) + 1L]] <- list(method = "ssgsea", contrast = ct,
                                                    rel_dir = s_rel,
                                                    table = s$res$table)
            }
        }
    }

    # mummichog now runs as its own pinned {targets} stage (06c/05b), wired into
    # the metabolomics DAG when modes.metabolomics.enrichment.mummichog.enabled is
    # true — it is no longer produced here. See mod_mummichog_pinned().

    # The manifest goes first: if it cannot be written the run stops here, with
    # nothing removed and no manifest, exactly like any other failed run.
    .write_metab_enrichment_manifest(build_metab_enrichment_manifest(units),
                                     manifest_file)
    files <- c(files, manifest_file)

    # Every result of this run is now on disk and recorded, so what an earlier
    # run left under this stage's own file names is no longer current. Done at
    # the end, not the start, so a run that fails keeps the earlier results it
    # did not replace. A run that produced no result at all (a wrong GMT path,
    # say) replaces nothing, so it removes nothing either: its empty manifest
    # already says that no unit on disk belongs to it.
    if (length(units) > 0) {
        .remove_stale_metab_enrichment_files(out_enr, out_ds, written = files)
    } else {
        message("metabolomics enrichment: this run produced no result, so ",
                "files left by an earlier run were not removed; the manifest ",
                "lists none of them.")
    }

    list(
        qea                = qea_res,
        ssgsea             = ssgsea_res,
        ora                = ora_res,
        gsea               = gsea_res,
        qea_by_contrast    = qea_by_contrast,
        ssgsea_by_contrast = ssgsea_by_contrast,
        ora_by_contrast    = ora_by_contrast,
        gsea_by_contrast   = gsea_by_contrast,
        plots              = plots,
        files              = unique(files)
    )
}


# ==============================================================================
# Output housekeeping
# ==============================================================================

#' Delete files and say which are really gone
#'
#' \code{unlink()} does not fail when a file cannot be deleted (it is open in
#' another program, or the folder is read-only), so its caller would otherwise
#' report a deletion that did not happen.
#'
#' @param paths Paths of files (not directories).
#' @return Logical vector, one per path: TRUE when the file no longer exists.
.delete_metab_enrichment_files <- function(paths) {
    if (length(paths) == 0) return(logical(0))
    unlink(paths)
    !file.exists(paths)
}

#' Remove this stage's files that the current run did not write
#'
#' Results left by an earlier run would otherwise sit beside the new ones and
#' read as current: a contrast that is gone or no longer yields a result, a
#' plot that was not drawn this time, the flat files from before the unit
#' layout. Only files the stage owns are candidates
#' (\code{list_metab_enrichment_owned_files()}); anything else in these folders
#' is left alone. Files are removed one by one; directories are never removed,
#' so an emptied unit directory stays.
#'
#' A file that could not be deleted is not counted as removed. It is named in a
#' warning instead, because it is still on disk and is not a result of this
#' run. The run itself is complete either way: its results and its manifest
#' were written before this is called.
#'
#' @param out_enr,out_ds The mode's Enrichment and Datasets directories.
#' @param written Paths of the files the current run wrote.
#' @return Invisibly, the paths that were actually removed.
.remove_stale_metab_enrichment_files <- function(out_enr, out_ds, written) {
    key   <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)
    owned <- list_metab_enrichment_owned_files(out_enr, out_ds)
    stale <- owned[!key(owned) %in% key(written)]
    if (length(stale) == 0) return(invisible(character(0)))

    gone    <- .delete_metab_enrichment_files(stale)
    removed <- stale[gone]
    left    <- stale[!gone]
    if (length(removed) > 0) {
        message("metabolomics enrichment: removed ", length(removed),
                " file(s) left by an earlier run under ", out_enr, " and ", out_ds)
    }
    if (length(left) > 0) {
        warning("metabolomics enrichment: ", length(left), " file(s) left by an ",
                "earlier run could not be removed. They are NOT results of this ",
                "run and are not in its manifest; delete them by hand: ",
                paste(left, collapse = ", "), call. = FALSE)
    }
    invisible(removed)
}

#' Write the enrichment manifest so that it is either complete or absent
#'
#' The table is written to a temporary file beside the target and then renamed
#' onto it. A rename within one directory replaces the target in one step, so a
#' failure while writing leaves no partial manifest: the target is untouched and
#' the temporary file is removed.
#'
#' @param manifest The manifest data.frame.
#' @param path Path of the manifest file.
#' @return \code{path}, invisibly.
.write_metab_enrichment_manifest <- function(manifest, path) {
    tmp <- paste0(path, ".tmp")
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    utils::write.csv(manifest, tmp, row.names = FALSE)
    if (!file.rename(tmp, path)) {
        stop("metabolomics enrichment: could not move the manifest into place at ",
             path, ". Check that the folder is writable and the file is not ",
             "open in another program.")
    }
    invisible(path)
}


# ==============================================================================
# Per-contrast helpers (ORA / GSEA)
# ==============================================================================

#' Sanitize a contrast label for use in a file name
#'
#' @param x Contrast name (any string).
#' @return A file-safe string: non-alphanumerics collapsed to "_", trimmed.
.sanitize_contrast <- function(x) {
    s <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
    gsub("^_+|_+$", "", s)
}

#' Render + save a set of enrichment plots for one result unit
#'
#' Each plot is written as \code{<key>.png} in the unit directory, which already
#' carries the method and the contrast.
#'
#' @param table    Result table to plot.
#' @param specs    Named list of \code{key -> function(table) -> ggplot}.
#' @param unit_dir Directory of the result unit (must exist).
#' @param label    Method and contrast, for the warning when a plot fails.
#' @return list(plots, files): ggplot objects (by key) and written PNG paths.
.render_enrichment_plots <- function(table, specs, unit_dir, label) {
    plots <- list()
    files <- character(0)
    for (key in names(specs)) {
        f <- file.path(unit_dir, paste0(key, ".png"))
        p <- tryCatch({
            pp <- specs[[key]](table)
            if (!is.null(pp)) {
                ggplot2::ggsave(f, pp, width = 12, height = 8, dpi = 300)
            }
            pp
        }, error = function(e) {
            warning(label, ": ", key, " plot failed: ", e$message)
            NULL
        })
        if (!is.null(p)) {
            plots[[key]] <- p
            files <- c(files, f)
        }
    }
    list(plots = plots, files = files)
}

#' Run ORA for one contrast and render/save its table + plots
#'
#' @param pre,de_res,config Passed through to \code{run_metabolomics_ora()}.
#' @param contrast Contrast name (a key of \code{de_res$de_tables}).
#' @param unit_dir Directory of this contrast's result unit; the table is
#'   written there as \code{results.tsv} beside its plots.
#' @return list(res, plots, files) or NULL when the contrast yields no table.
.mod_metab_ora_contrast <- function(pre, de_res, config, contrast, unit_dir) {
    res <- tryCatch(
        run_metabolomics_ora(pre, de_res, config, contrast = contrast),
        error = function(e) {
            warning("metabolomics ORA failed for '", contrast, "': ", e$message)
            NULL
        }
    )
    if (is.null(res) || is.null(res$table)) return(NULL)

    files <- save_tsv(res$table, unit_dir, "results.tsv")
    rendered <- .render_enrichment_plots(
        res$table,
        specs = list(
            dotplot  = function(t) plot_ora_dotplot(t, top_n = 20,
                          title = paste0("ORA — ", contrast)),
            lollipop = function(t) plot_ora_lollipop(t, top_n = 20,
                          title = paste0("ORA — ", contrast)),
            barplot  = function(t) plot_enrichment_barplot(t, top_n = 20,
                          title = paste0("ORA — ", contrast))
        ),
        unit_dir = unit_dir,
        label    = paste0("metabolomics ORA '", contrast, "'")
    )
    list(res = res, plots = rendered$plots,
         files = c(files, rendered$files))
}

#' Run GSEA for one contrast and render/save its table + plots
#'
#' @inheritParams .mod_metab_ora_contrast
#' @return list(res, plots, files) or NULL when the contrast yields no table.
.mod_metab_gsea_contrast <- function(pre, de_res, config, contrast, unit_dir) {
    res <- tryCatch(
        run_metabolomics_gsea(pre, de_res, config, contrast = contrast),
        error = function(e) {
            warning("metabolomics GSEA failed for '", contrast, "': ", e$message)
            NULL
        }
    )
    if (is.null(res) || is.null(res$table)) return(NULL)

    files <- save_tsv(res$table, unit_dir, "results.tsv")
    rendered <- .render_enrichment_plots(
        res$table,
        specs = list(
            nes_barplot = function(t) plot_gsea_nes_barplot(t, top_n = 20,
                             title = paste0("GSEA — ", contrast)),
            lollipop    = function(t) plot_gsea_lollipop(t, top_n = 20,
                             title = paste0("GSEA — ", contrast)),
            barplot     = function(t) plot_enrichment_barplot(t, top_n = 20,
                             title = paste0("GSEA — ", contrast))
        ),
        unit_dir = unit_dir,
        label    = paste0("metabolomics GSEA '", contrast, "'")
    )
    list(res = res, plots = rendered$plots,
         files = c(files, rendered$files))
}

#' Write the first-contrast enrichment outputs under legacy (non-suffixed) names
#'
#' External consumers (Shiny payload plumbing, scripts/build_per_omics_pptx.py)
#' expect the non-suffixed \code{enrichment_<method>_results.tsv} and plot PNGs.
#' Re-emit them from the already-built first-contrast objects so those tools keep
#' working alongside the new per-contrast files.
#'
#' @param method  "ora" or "gsea".
#' @param res     First-contrast result (with \code{$table}).
#' @param plots   First-contrast plot list (keyed as produced above).
#' @param out_enr,out_ds Output directories.
#' @return Character vector of written file paths.
.save_legacy_enrichment <- function(method, res, plots, out_enr, out_ds) {
    files <- save_tsv(res$table, out_ds,
                      paste0("enrichment_", method, "_results.tsv"))
    legacy <- if (method == "ora") {
        c(dotplot     = "enrichment_ora_dotplot.png",
          lollipop    = "enrichment_ora_lollipop.png",
          barplot     = "enrichment_ora_barplot.png")
    } else {
        c(nes_barplot = "enrichment_gsea_nes_barplot.png",
          lollipop    = "enrichment_gsea_lollipop.png",
          barplot     = "enrichment_gsea_barplot.png")
    }
    for (key in names(legacy)) {
        if (!is.null(plots[[key]])) {
            f <- file.path(out_enr, legacy[[key]])
            tryCatch(
                ggplot2::ggsave(f, plots[[key]], width = 12, height = 8, dpi = 300),
                error = function(e) NULL
            )
            files <- c(files, f)
        }
    }
    files
}


# ==============================================================================
# Per-contrast QEA / ssGSEA helpers
# ==============================================================================

#' Resolve a contrast label to its two groups (from inputs$contrasts)
#'
#' Looks the contrast up in \code{inputs$contrasts} (a data.frame with
#' Numerator/Denominator, optionally Contrast_name and Factor) and returns the
#' grouping column and the two group values. Matches on Contrast_name first,
#' then on the "Num - Den" / "Num_vs_Den" label styles.
#'
#' @param inputs        Loaded inputs list (uses \code{$contrasts}).
#' @param contrast_label A key of \code{de_res$de_tables}.
#' @param fallback_col  Grouping column to use when the contrast has no Factor.
#' @return list(group_col, groups) or NULL when it cannot be resolved.
.resolve_contrast_groups <- function(inputs, contrast_label, fallback_col) {
    ct_df <- inputs$contrasts
    if (is.null(ct_df) || !is.data.frame(ct_df) ||
        !all(c("Numerator", "Denominator") %in% colnames(ct_df))) {
        return(NULL)
    }
    row <- NULL
    if ("Contrast_name" %in% colnames(ct_df)) {
        hit <- which(as.character(ct_df$Contrast_name) == contrast_label)
        if (length(hit) >= 1) row <- ct_df[hit[1], , drop = FALSE]
    }
    if (is.null(row)) {
        lab_dash <- paste(ct_df$Numerator, ct_df$Denominator, sep = " - ")
        lab_vs   <- paste(ct_df$Numerator, ct_df$Denominator, sep = "_vs_")
        hit <- which(lab_dash == contrast_label | lab_vs == contrast_label)
        if (length(hit) >= 1) row <- ct_df[hit[1], , drop = FALSE]
    }
    if (is.null(row)) return(NULL)

    group_col <- if ("Factor" %in% colnames(ct_df) &&
                     !is.na(row$Factor) && nzchar(as.character(row$Factor))) {
        as.character(row$Factor)
    } else {
        fallback_col
    }
    list(group_col = group_col,
         groups    = c(as.character(row$Numerator), as.character(row$Denominator)))
}

#' Subset a preprocessing object to the samples in a given set of groups
#'
#' Keeps the metadata rows (and the matching expression-matrix columns) whose
#' \code{group_col} value is one of \code{groups}. row_data (features) is left
#' unchanged. Returns NULL when fewer than two samples remain.
#'
#' @param pre        Preprocessing list (expr_*, meta, row_data).
#' @param group_col  Metadata column defining the groups.
#' @param groups     Character vector of group values to keep.
#' @param sample_col Metadata column holding sample ids (matches matrix columns).
#' @return A pre-like list subset to the requested samples, or NULL.
.subset_pre_to_groups <- function(pre, group_col, groups, sample_col) {
    meta <- pre$meta
    if (is.null(meta) || !group_col %in% colnames(meta)) return(NULL)
    keep_rows <- as.character(meta[[group_col]]) %in% groups
    if (sum(keep_rows) < 2) return(NULL)
    keep_ids <- as.character(meta[[sample_col]][keep_rows])

    sub <- pre
    sub$meta <- meta[keep_rows, , drop = FALSE]
    for (m in c("expr_raw", "expr_work", "expr_filt", "expr_log", "expr_unscaled")) {
        if (!is.null(pre[[m]])) {
            cols   <- intersect(colnames(pre[[m]]), keep_ids)
            sub[[m]] <- pre[[m]][, cols, drop = FALSE]
        }
    }
    sub
}

#' Run QEA for one contrast (on a two-group subset) and render/save its outputs
#'
#' @param pre_sub Two-group subset of the preprocessing object.
#' @param config,contrast,unit_dir As in the other per-contrast helpers.
#' @return list(res, plots, files) or NULL when the contrast yields no table.
.mod_metab_qea_contrast <- function(pre_sub, config, contrast, unit_dir) {
    res <- tryCatch(
        run_metabolomics_qea(pre_sub, config),
        error = function(e) {
            warning("metabolomics QEA failed for '", contrast, "': ", e$message)
            NULL
        }
    )
    if (is.null(res) || is.null(res$table)) return(NULL)

    files <- save_tsv(res$table, unit_dir, "results.tsv")
    rendered <- .render_enrichment_plots(
        res$table,
        specs = list(
            lollipop = function(t) plot_qea_lollipop(t, top_n = 20,
                          title = paste0("QEA — ", contrast)),
            barplot  = function(t) plot_enrichment_barplot(t, top_n = 20,
                          title = paste0("QEA — ", contrast))
        ),
        unit_dir = unit_dir,
        label    = paste0("metabolomics QEA '", contrast, "'")
    )
    list(res = res, plots = rendered$plots, files = c(files, rendered$files))
}

#' Run ssGSEA for one contrast (on a two-group subset) and render/save outputs
#'
#' The two-group subset is exactly what ssGSEA's Wilcoxon step needs, so a
#' dataset with more than two groups still yields per-contrast ssGSEA results.
#'
#' @inheritParams .mod_metab_qea_contrast
#' @param groups `c(numerator, denominator)` for this contrast, from
#'   `.resolve_contrast_groups()`. Passed through so `score_diff` is reported in
#'   the contrast's own direction rather than in alphabetical order. Groups the
#'   scored samples do not carry abort the scoring, and the `tryCatch` below
#'   turns that into a warning naming this contrast and skips it — a contrast
#'   whose direction cannot be established yields no result rather than a
#'   result with an arbitrary sign.
#' @param group_col The column those groups live in, from
#'   `.resolve_contrast_groups()$group_col` — the contrast's own `Factor` when
#'   it has one. It travels WITH the groups: the subset above was built from
#'   this column, so leaving ssGSEA to read the configured
#'   `enrichment.condition_column` instead would label the samples by an
#'   unrelated variable, and the groups would not be found at all.
#' @return list(res, plots, files) or NULL when the contrast yields no table.
.mod_metab_ssgsea_contrast <- function(pre_sub, config, contrast, unit_dir,
                                       groups = NULL, group_col = NULL) {
    res <- tryCatch(
        run_metabolomics_ssgsea(pre_sub, config, contrast_groups = groups,
                                contrast_column = group_col),
        error = function(e) {
            warning("metabolomics ssGSEA failed for '", contrast, "': ", e$message)
            NULL
        }
    )
    if (is.null(res) || is.null(res$table)) return(NULL)

    files <- save_tsv(res$table, unit_dir, "results.tsv")
    rendered <- .render_enrichment_plots(
        res$table,
        specs = list(
            lollipop = function(t) plot_ssgsea_lollipop(t, top_n = 20,
                          title = paste0("ssGSEA — ", contrast)),
            barplot  = function(t) plot_enrichment_barplot(t, top_n = 20,
                          title = paste0("ssGSEA — ", contrast))
        ),
        unit_dir = unit_dir,
        label    = paste0("metabolomics ssGSEA '", contrast, "'")
    )
    list(res = res, plots = rendered$plots, files = c(files, rendered$files))
}
