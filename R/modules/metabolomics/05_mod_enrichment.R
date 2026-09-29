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
#' @param pre     List from preprocess_metabolomics().
#' @param de_res  DE results (needed by ORA and GSEA; NULL if unavailable).
#' @param config  Full config.
#' @param out_dir Output directory for this mode.
#' @return list(qea, ssgsea, ora, gsea, plots, files).
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

    plots <- list()
    files <- character(0)

    # ---- QEA ----
    qea_res <- tryCatch(
        run_metabolomics_qea(pre, config),
        error = function(e) {
            warning("metabolomics QEA failed: ", e$message)
            NULL
        }
    )

    if (!is.null(qea_res) && !is.null(qea_res$table)) {
        f_qea <- save_tsv(qea_res$table, out_ds, "enrichment_qea_results.tsv")
        files <- c(files, f_qea)

        # QEA barplot
        f_qea_plot <- file.path(out_enr, "enrichment_qea_barplot.png")
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
        f_qea_loll <- file.path(out_enr, "enrichment_qea_lollipop.png")
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
                f_lib <- file.path(out_enr,
                                    paste0("enrichment_qea_", tolower(lib_label), ".png"))
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
        f_ssgsea <- save_tsv(ssgsea_res$table, out_ds, "enrichment_ssgsea_results.tsv")
        files <- c(files, f_ssgsea)

        # ssGSEA barplot
        f_ssgsea_plot <- file.path(out_enr, "enrichment_ssgsea_barplot.png")
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

            f_box <- file.path(out_enr, "enrichment_ssgsea_boxplots.png")
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
        f_ssgsea_loll <- file.path(out_enr, "enrichment_ssgsea_lollipop.png")
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
            one <- .mod_metab_ora_contrast(pre, de_res, config, ct, out_enr, out_ds)
            if (!is.null(one)) ora_by_contrast[[ct]] <- one
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
            one <- .mod_metab_gsea_contrast(pre, de_res, config, ct, out_enr, out_ds)
            if (!is.null(one)) gsea_by_contrast[[ct]] <- one
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
            q <- .mod_metab_qea_contrast(pre_sub, config, ct, out_enr, out_ds)
            if (!is.null(q)) { qea_by_contrast[[ct]] <- q; files <- c(files, q$files) }
            s <- .mod_metab_ssgsea_contrast(pre_sub, config, ct, out_enr, out_ds)
            if (!is.null(s)) { ssgsea_by_contrast[[ct]] <- s; files <- c(files, s$files) }
        }
    }

    # mummichog now runs as its own pinned {targets} stage (06c/05b), wired into
    # the metabolomics DAG when modes.metabolomics.enrichment.mummichog.enabled is
    # true — it is no longer produced here. See mod_mummichog_pinned().

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

#' Render + save a set of enrichment plots for one contrast
#'
#' @param table   Result table to plot.
#' @param specs   Named list of \code{key -> function(table) -> ggplot}.
#' @param prefix  File-name prefix, e.g. "enrichment_ora".
#' @param sfx     Sanitized contrast suffix.
#' @param out_enr Enrichment output directory.
#' @return list(plots, files): ggplot objects (by key) and written PNG paths.
.render_enrichment_plots <- function(table, specs, prefix, sfx, out_enr) {
    plots <- list()
    files <- character(0)
    for (key in names(specs)) {
        f <- file.path(out_enr, paste0(prefix, "_", key, "_", sfx, ".png"))
        p <- tryCatch({
            pp <- specs[[key]](table)
            if (!is.null(pp)) {
                ggplot2::ggsave(f, pp, width = 12, height = 8, dpi = 300)
            }
            pp
        }, error = function(e) {
            warning(prefix, " ", key, " failed for '", sfx, "': ", e$message)
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
#' @param out_enr,out_ds Enrichment / datasets output directories.
#' @return list(res, plots, files) or NULL when the contrast yields no table.
.mod_metab_ora_contrast <- function(pre, de_res, config, contrast,
                                    out_enr, out_ds) {
    res <- tryCatch(
        run_metabolomics_ora(pre, de_res, config, contrast = contrast),
        error = function(e) {
            warning("metabolomics ORA failed for '", contrast, "': ", e$message)
            NULL
        }
    )
    if (is.null(res) || is.null(res$table)) return(NULL)

    sfx   <- .sanitize_contrast(contrast)
    files <- save_tsv(res$table, out_ds,
                      paste0("enrichment_ora_", sfx, "_results.tsv"))
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
        prefix = "enrichment_ora", sfx = sfx, out_enr = out_enr
    )
    list(res = res, plots = rendered$plots,
         files = c(files, rendered$files))
}

#' Run GSEA for one contrast and render/save its table + plots
#'
#' @inheritParams .mod_metab_ora_contrast
#' @return list(res, plots, files) or NULL when the contrast yields no table.
.mod_metab_gsea_contrast <- function(pre, de_res, config, contrast,
                                     out_enr, out_ds) {
    res <- tryCatch(
        run_metabolomics_gsea(pre, de_res, config, contrast = contrast),
        error = function(e) {
            warning("metabolomics GSEA failed for '", contrast, "': ", e$message)
            NULL
        }
    )
    if (is.null(res) || is.null(res$table)) return(NULL)

    sfx   <- .sanitize_contrast(contrast)
    files <- save_tsv(res$table, out_ds,
                      paste0("enrichment_gsea_", sfx, "_results.tsv"))
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
        prefix = "enrichment_gsea", sfx = sfx, out_enr = out_enr
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
    for (m in c("expr_raw", "expr_work", "expr_filt", "expr_log")) {
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
#' @param config,contrast,out_enr,out_ds As in the other per-contrast helpers.
#' @return list(res, plots, files) or NULL when the contrast yields no table.
.mod_metab_qea_contrast <- function(pre_sub, config, contrast, out_enr, out_ds) {
    res <- tryCatch(
        run_metabolomics_qea(pre_sub, config),
        error = function(e) {
            warning("metabolomics QEA failed for '", contrast, "': ", e$message)
            NULL
        }
    )
    if (is.null(res) || is.null(res$table)) return(NULL)

    sfx   <- .sanitize_contrast(contrast)
    files <- save_tsv(res$table, out_ds,
                      paste0("enrichment_qea_", sfx, "_results.tsv"))
    rendered <- .render_enrichment_plots(
        res$table,
        specs = list(
            lollipop = function(t) plot_qea_lollipop(t, top_n = 20,
                          title = paste0("QEA — ", contrast)),
            barplot  = function(t) plot_enrichment_barplot(t, top_n = 20,
                          title = paste0("QEA — ", contrast))
        ),
        prefix = "enrichment_qea", sfx = sfx, out_enr = out_enr
    )
    list(res = res, plots = rendered$plots, files = c(files, rendered$files))
}

#' Run ssGSEA for one contrast (on a two-group subset) and render/save outputs
#'
#' The two-group subset is exactly what ssGSEA's Wilcoxon step needs, so a
#' dataset with more than two groups still yields per-contrast ssGSEA results.
#'
#' @inheritParams .mod_metab_qea_contrast
#' @return list(res, plots, files) or NULL when the contrast yields no table.
.mod_metab_ssgsea_contrast <- function(pre_sub, config, contrast, out_enr, out_ds) {
    res <- tryCatch(
        run_metabolomics_ssgsea(pre_sub, config),
        error = function(e) {
            warning("metabolomics ssGSEA failed for '", contrast, "': ", e$message)
            NULL
        }
    )
    if (is.null(res) || is.null(res$table)) return(NULL)

    sfx   <- .sanitize_contrast(contrast)
    files <- save_tsv(res$table, out_ds,
                      paste0("enrichment_ssgsea_", sfx, "_results.tsv"))
    rendered <- .render_enrichment_plots(
        res$table,
        specs = list(
            lollipop = function(t) plot_ssgsea_lollipop(t, top_n = 20,
                          title = paste0("ssGSEA — ", contrast)),
            barplot  = function(t) plot_enrichment_barplot(t, top_n = 20,
                          title = paste0("ssGSEA — ", contrast))
        ),
        prefix = "enrichment_ssgsea", sfx = sfx, out_enr = out_enr
    )
    list(res = res, plots = rendered$plots, files = c(files, rendered$files))
}
