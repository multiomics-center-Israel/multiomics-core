#' Multi-Omics Report Rendering
#'
#' Copies the canonical Rmd template into the run directory and renders
#' a self-contained HTML report.

#' Render the multi-omics integration report
#'
#' @param run_dir  The results run directory (e.g. outputs/project/Results_.../multiomics)
#' @param config   Full pipeline config list
#' @param config_file Path to the original YAML config file. Retained for the
#'   pipeline call sites; the config snapshot is serialized from \code{config}
#'   rather than copied from this path, so that it matches what the
#'   \code{execution_info_files} target writes.
#' @return Path to the rendered HTML file (character, format = "file")
#' @export
render_multiomics_report <- function(run_dir, config, config_file = NULL) {

    # Locate the canonical Rmd template shipped with multiomics-core
    template_path <- system.file("report_template_multiomics.Rmd",
                                  package = "multiomics.core",
                                  mustWork = FALSE)

    # Fallback: template lives alongside this source file
    if (!nzchar(template_path) || !file.exists(template_path)) {
        src_dir <- system.file("R", "domain", "multiomics",
                                package = "multiomics.core",
                                mustWork = FALSE)
        if (!nzchar(src_dir)) {
            # When running via targets::tar_source(), try the working directory first,
            # then fall back to the project root from config
            src_dir <- file.path("R", "domain", "multiomics")
            if (!file.exists(file.path(src_dir, "report_template_multiomics.Rmd"))) {
                proj_dir <- config$project$dir %||% "."
                src_dir <- file.path(proj_dir, "R", "domain", "multiomics")
            }
        }
        template_path <- file.path(src_dir, "report_template_multiomics.Rmd")
    }

    if (!file.exists(template_path)) {
        warning("Multi-omics report template not found at: ", template_path,
                ". Skipping report generation.")
        return(NA_character_)
    }

    # Ensure run_dir exists, then copy template into it (the knit directory)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    dest_rmd <- file.path(run_dir, "report_multiomics.Rmd")
    file.copy(template_path, dest_rmd, overwrite = TRUE)

    # Write execution_info/config_used.yaml (needed by the template). Rewritten on
    # every render rather than only when absent: the old guard left the first run's
    # snapshot in place forever, so config edits appeared to do nothing.
    #
    # Always serialized from `config`, never copied from the source YAML. This can
    # be the same file the execution_info_files target owns, and that target writes
    # it with this same yaml::write_yaml(config, ...). Copying the raw YAML here
    # would put different content at a tracked path -- load_config() adds
    # .config_path and .config_mtime, which the source file does not carry --
    # leaving the target outdated the moment the report finished.
    exec_dir <- file.path(run_dir, "execution_info")
    dir.create(exec_dir, recursive = TRUE, showWarnings = FALSE)
    yaml::write_yaml(config, file.path(exec_dir, "config_used.yaml"))

    # Render
    out_html <- file.path(run_dir, "report_multiomics.html")

    # Ensure any stale report from a previous run is gone, as the proteomics
    # renderer already does. The config snapshot above is now rewritten every
    # time, so a render that fails afterwards would otherwise leave last run's
    # HTML sitting beside this run's snapshot, reading as though that report had
    # been produced from these settings.
    if (file.exists(out_html)) try(file.remove(out_html), silent = TRUE)

    message("Rendering multi-omics report to: ", out_html)

    rmarkdown::render(
        input       = dest_rmd,
        output_file = out_html,
        output_dir  = run_dir,
        quiet       = TRUE,
        envir       = new.env(parent = globalenv())
    )

    if (file.exists(out_html)) {
        message("Multi-omics report rendered successfully: ", out_html)
    } else {
        warning("Multi-omics report rendering did not produce output file.")
    }

    out_html
}


#' How many contrasts the multi-omics report shows, and the name of a lone one
#'
#' The report collapses each section to a single tab when a run compared one
#' contrast, because the per-contrast view would then repeat the combined one.
#' Whether that is so is read from what the run produced as well as from what
#' was configured: the design block records the intent, but the pipeline writes
#' one `per_contrast/` directory per contrast it actually found in the data,
#' and a run fed per-omics DE tables can hold more contrasts than its design
#' lists. Each of the three per-contrast report sections has its own producer
#' -- enrichment, MultiGSEA and Multi-ORA, the last of which runs even when
#' MultiGSEA does not -- so all three are counted. The count is the largest of
#' the design and those outputs, so a run is only treated as single-contrast
#' when neither the config nor any output says otherwise.
#'
#' @param design_contrasts The config's \code{design$contrasts}, or NULL.
#' @param enrichment_dir Cross-omics enrichment output directory.
#' @param multigsea_dir MultiGSEA output directory; Multi-ORA writes under its
#'   \code{multi_ora/} subdirectory.
#' @return List with \code{n_contrasts}, \code{single_contrast} (logical) and
#'   \code{label}: the lone contrast's name for a single-contrast run (from its
#'   output directory, else the design), or "Results" when none is known.
#' @examples
#' report_contrast_layout(list("A_vs_B"), tempdir(), tempdir())$single_contrast
report_contrast_layout <- function(design_contrasts, enrichment_dir, multigsea_dir) {
    contrast_dirs <- function(d) {
        d <- file.path(d, "per_contrast")
        if (dir.exists(d)) sort(list.dirs(d, full.names = FALSE, recursive = FALSE))
        else character(0)
    }
    design <- as.character(unlist(design_contrasts, use.names = FALSE))
    enrich <- contrast_dirs(enrichment_dir)
    mg <- contrast_dirs(multigsea_dir)
    mora <- contrast_dirs(file.path(multigsea_dir, "multi_ora"))

    n_contrasts <- max(length(design), length(enrich), length(mg), length(mora))
    label <- if (length(enrich) == 1) {
        gsub("_", " ", enrich)
    } else if (length(mg) == 1) {
        gsub("_", " ", mg)
    } else if (length(mora) == 1) {
        gsub("_", " ", mora)
    } else if (length(design) == 1) {
        design
    } else {
        "Results"
    }
    list(n_contrasts = n_contrasts, single_contrast = n_contrasts <= 1,
         label = label)
}


#' Human-readable title for a rendered pathview map
#'
#' pathview leaves the KGML beside the PNG it renders, and that file carries the
#' pathway's title. Reading it avoids both a network call and a bare numeric id
#' as the section heading.
#'
#' @param png_path Path to the rendered pathview PNG.
#' @return Character title, or NA when no KGML sits beside the image.
pathview_map_title <- function(png_path) {
    stem <- sub("\\..*$", "", basename(png_path))
    kgml <- file.path(dirname(png_path), paste0(stem, ".xml"))
    if (!file.exists(kgml)) return(NA_character_)
    hit <- grep("title=", readLines(kgml, n = 40, warn = FALSE), value = TRUE)
    if (length(hit) == 0) return(NA_character_)
    ttl <- sub('.*title="([^"]+)".*', "\\1", hit[1])
    if (identical(ttl, hit[1])) NA_character_ else ttl
}


#' Contrast a Multi-ORA pathway map was drawn for, read from its filename
#'
#' Multi-ORA names each map `<org><id>.multi_ora_<contrast>[.multi].png`, keyed
#' on the canonical contrast, which keeps apart two contrasts that make.names()
#' would collide. The readable spelling comes from the renderer's own record of
#' that key; a key it did not record falls back to the key with its dots read as
#' spaces. A single-contrast map (`.multi_ora.png`) and the per-layer maps
#' (`.metab_top`, `.prot_top`) name no contrast and get "".
#'
#' @param png_path Character vector of rendered pathview PNG paths.
#' @param contrast_labels Named list, filename key -> readable contrast, as the
#'   union renderer records it in its sidecar YAML; empty for the other one.
#' @return Character vector as long as \code{png_path}.
#' @examples
#' pathview_map_contrast("ko00010.multi_ora_A.vs.B.multi.png")  # "A vs B"
pathview_map_contrast <- function(png_path, contrast_labels = list()) {
    bn <- basename(png_path)
    key <- ifelse(grepl("\\.multi_ora", bn),
                  sub("\\.(multi\\.)?png$", "", sub("^.*\\.multi_ora_?", "", bn)),
                  "")
    vapply(key, function(k) {
        if (!nzchar(k)) return("")
        lab <- contrast_labels[[k]]
        if (is.null(lab)) gsub("\\.", " ", k) else as.character(lab)[1]
    }, character(1), USE.NAMES = FALSE)
}


#' Every pathway map the multi-omics report shows, one row per map
#'
#' The report lists all maps in one table and draws them in one tab per map set,
#' and both read this index, so a map's name and contrast cannot differ between
#' the table and its tab. Each set is included only when the report says it was
#' produced (its PDF exists), as the tabs always were.
#'
#' @param pathview_dir The Multi-ORA \code{pathview/} directory.
#' @param multi_ora,metab_top,prot_top Logical: include the cross-omics maps, the
#'   top metabolomics pathways, the top proteomics pathways.
#' @param multi_ora_label Display name of the cross-omics set, from
#'   \code{pathview_multi_ora_support()}.
#' @param contrast_labels Passed to \code{pathview_map_contrast()}.
#' @return Data frame with columns \code{set} (\code{"multi_ora"},
#'   \code{"metab_top"} or \code{"prot_top"}), \code{map_set} (its display
#'   name), \code{contrast} ("" where the map names none), \code{kegg_id},
#'   \code{pathway} ("" without a KGML title), \code{heading} and \code{png};
#'   zero rows when there is nothing to show.
pathview_map_index <- function(pathview_dir, multi_ora = TRUE, metab_top = TRUE,
                               prot_top = TRUE, multi_ora_label = "Cross-omics pathways",
                               contrast_labels = list()) {
    sets <- list(
        list(on = multi_ora, key = "multi_ora",
             pattern = "\\.multi_ora.*\\.(multi\\.)?png$", strip = "\\.multi_ora.*$",
             label = multi_ora_label),
        # The per-layer renderer keeps whichever of `.png` and `.multi.png`
        # pathview wrote, so both are maps.
        list(on = metab_top, key = "metab_top",
             pattern = "\\.metab_top(\\.multi)?\\.png$",
             strip = "\\.metab_top(\\.multi)?\\.png$", label = "Top metabolomics pathways"),
        list(on = prot_top, key = "prot_top",
             pattern = "\\.prot_top(\\.multi)?\\.png$",
             strip = "\\.prot_top(\\.multi)?\\.png$", label = "Top proteomics pathways"))
    empty <- data.frame(set = character(0), map_set = character(0),
                        contrast = character(0), kegg_id = character(0),
                        pathway = character(0), heading = character(0),
                        png = character(0), stringsAsFactors = FALSE)
    if (!dir.exists(pathview_dir)) return(empty)

    rows <- lapply(sets, function(s) {
        if (!isTRUE(s$on)) return(NULL)
        pngs <- sort(list.files(pathview_dir, pattern = s$pattern, full.names = TRUE))
        if (length(pngs) == 0) return(NULL)
        ids <- sub("^[a-z]+", "", sub(s$strip, "", basename(pngs)))
        titles <- vapply(pngs, pathview_map_title, character(1), USE.NAMES = FALSE)
        contrast <- pathview_map_contrast(pngs, contrast_labels)
        name <- ifelse(is.na(titles), paste0("Pathway ", ids),
                       paste0(titles, " (", ids, ")"))
        data.frame(set = s$key, map_set = s$label, contrast = contrast,
                   kegg_id = ids, pathway = ifelse(is.na(titles), "", titles),
                   heading = ifelse(nzchar(contrast), paste0(name, " - ", contrast), name),
                   png = pngs, stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, c(list(empty), rows))
    rownames(out) <- NULL
    out
}


#' How the report names the cross-omics pathway maps, and what it says of them
#'
#' The supported renderer prefers pathways enriched in two or more layers but
#' falls back to a single layer when none qualify, and records which in its
#' sidecar; the no-OrgDb renderer unions the hits of either gene layer. The
#' overview's set name and the tab's opening sentence both come from here, so
#' neither can claim more layers than the maps rest on. A run from before the
#' record existed gets wording that claims no number of layers.
#'
#' @param is_union Logical: the maps came from the no-OrgDb union renderer.
#' @param support_layers The supported renderer's recorded
#'   \code{support_layers}, or NULL when its sidecar has none.
#' @return List with \code{label} (the overview's map-set name) and
#'   \code{intro} (the tab's opening sentence, Markdown).
#' @examples
#' pathview_multi_ora_support(FALSE, 1L)$label  # "Enriched in one layer or more"
pathview_multi_ora_support <- function(is_union, support_layers = NULL) {
    if (isTRUE(is_union)) {
        return(list(
            label = "Enriched in a gene layer",
            intro = paste("KEGG reference maps, in KO space, for pathways enriched in",
                          "**at least one** gene-based omics layer, shown **per contrast**.",
                          "A map here can rest on transcriptomics or proteomics alone.")))
    }
    n <- suppressWarnings(as.integer(unlist(support_layers))[1])
    if (length(n) == 0 || is.na(n) || n < 1) {
        list(label = "Cross-omics pathways",
             intro = "KEGG pathway maps for the top cross-omics pathways, shown **per contrast**.")
    } else if (n == 1) {
        list(label = "Enriched in one layer or more",
             intro = paste("KEGG pathway maps for pathways enriched in **at least one**",
                           "omics layer, shown **per contrast**. No pathway reached two",
                           "layers in this run, so a map here can rest on one layer alone."))
    } else {
        list(label = sprintf("Enriched in >= %d layers", n),
             intro = sprintf(paste("KEGG pathway maps for pathways enriched in",
                                   "**>= %d omics layers**, shown **per contrast**."), n))
    }
}


#' What the pathway-map section says when it has no maps
#'
#' Maps can be missing because the project switched them off or because the run
#' drew none, and the two need different sentences. The switch is read with the
#' same key, default and coercion as \code{run_multi_ora()}, which owns it for
#' every renderer the section shows.
#'
#' @param multi_cfg The config's \code{modes$multiomics} section.
#' @return A one-line Markdown note.
pathview_missing_note <- function(multi_cfg) {
    run_pathview <- isTRUE((multi_cfg$enrichment$pathview$run_pathview %||% TRUE))
    if (!run_pathview) {
        "*Pathway maps are switched off for this run (`enrichment.pathview.run_pathview: false`).*"
    } else {
        "*No pathway maps are available for this run.*"
    }
}
