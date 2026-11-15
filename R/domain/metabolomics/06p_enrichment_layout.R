# R/domain/metabolomics/06p_enrichment_layout.R
#
# Output layout of the metabolomics enrichment stage (QEA, ssGSEA, ORA, GSEA).
#
# Every path under <mode>/Enrichment/ is built HERE, so the tree is defined once.
# It follows the RNA-seq layout (R/core/09_enrichment.R: ora_unit_dir(),
# gsea_unit_dir()): the analysis context lives in the DIRECTORY, and file names
# inside a unit are short and fixed (results.tsv, dotplot.png, ...).
#
#   Enrichment/
#     enrichment_manifest.csv
#     ORA/<db>/contrasts/<contrast>/de/
#     GSEA/<db>/ranking_by_fc/<contrast>/
#     QEA/<db>/all_samples/            QEA/<db>/contrasts/<contrast>/
#       (all_samples also holds library_<name>.png, one per GMT library)
#     ssGSEA/<db>/all_samples/         ssGSEA/<db>/contrasts/<contrast>/
#
# This file organises files only. It computes no statistic and reads no result
# value other than the FDR column it counts for the manifest.

# The four methods test all configured GMT libraries as ONE family (FDR is
# adjusted across libraries), so there is a single database level. The library
# of each pathway stays a column of the result table where the method records it.
.METAB_ENRICHMENT_DB <- "combined"

# Method key -> directory name.
.METAB_ENRICHMENT_METHOD_DIR <- c(ora = "ORA", gsea = "GSEA",
                                  qea = "QEA", ssgsea = "ssGSEA")

.METAB_ENRICHMENT_MANIFEST_FILE <- "enrichment_manifest.csv"

# The plot files each method writes into a unit, by file stem. These, and
# results.tsv, are the only names the stage treats as its own inside a unit.
.METAB_ENRICHMENT_PLOT_KEYS <- list(
    ora    = c("dotplot", "lollipop", "barplot"),
    gsea   = c("nes_barplot", "lollipop", "barplot"),
    qea    = c("barplot", "lollipop"),
    ssgsea = c("barplot", "boxplots", "lollipop")
)

# The all-samples QEA unit also holds one barplot per GMT library. The library
# name is free text, so the prefix is what marks the file as the stage's.
.METAB_ENRICHMENT_LIBRARY_PLOT_PREFIX <- "library_"

# FDR level the manifest counts pathways at. It is the value the enrichment code
# already uses for its `significant` flag, its plot threshold lines and its log
# messages; there is no config key for it. Display/bookkeeping only.
.METAB_ENRICHMENT_SIG_FDR <- 0.05


#' Directory names for a set of contrasts, guaranteed distinct
#'
#' A contrast label becomes one path component. Two labels can sanitise to the
#' same token (for example "A vs B" and "A/vs/B"), and two tokens that differ
#' only in case are the same directory on a case-insensitive file system. Either
#' would make one contrast overwrite another, so a colliding token gets a numeric
#' suffix and a warning names the contrasts involved.
#'
#' @param contrasts Character vector of contrast labels (the keys of
#'   \code{de_res$de_tables}).
#' @param tokens Character vector, same length: each label already sanitised for
#'   use in a file name.
#' @return Named character vector (names = \code{contrasts}) of directory names
#'   that are non-empty, never "." or "..", never a Windows device name, and
#'   distinct ignoring case.
metab_enrichment_contrast_dirs <- function(contrasts, tokens) {
    contrasts <- as.character(contrasts)
    tokens    <- as.character(tokens)
    if (length(contrasts) != length(tokens)) {
        stop("metab_enrichment_contrast_dirs(): 'contrasts' and 'tokens' must ",
             "have the same length (got ", length(contrasts), " and ",
             length(tokens), ").")
    }
    if (length(contrasts) == 0) return(stats::setNames(character(0), character(0)))

    # A trailing dot is dropped by Windows, and "." / ".." are not directories
    # of their own; an empty token would collapse a path level.
    tokens <- sub("\\.+$", "", tokens)
    tokens[is.na(tokens) | !nzchar(tokens)] <- "contrast"
    # Windows reserves these names for devices, with or without an extension,
    # and refuses to create a directory called by one of them.
    reserved <- grepl("^(con|prn|aux|nul|com[1-9]|lpt[1-9])(\\..*)?$", tokens,
                      ignore.case = TRUE)
    tokens[reserved] <- paste0("contrast_", tokens[reserved])

    out  <- character(length(tokens))
    used <- character(0)
    for (i in seq_along(tokens)) {
        candidate <- tokens[i]
        k <- 1L
        while (tolower(candidate) %in% used) {
            k <- k + 1L
            candidate <- paste0(tokens[i], "_", k)
        }
        out[i] <- candidate
        used   <- c(used, tolower(candidate))
    }

    renamed <- out != tokens
    if (any(renamed)) {
        warning("metabolomics enrichment: contrast names collide once made ",
                "file-safe, so these were given a numeric suffix to keep their ",
                "outputs apart: ",
                paste0("'", contrasts[renamed], "' -> '", out[renamed], "'",
                       collapse = ", "), ".", call. = FALSE)
    }
    stats::setNames(out, contrasts)
}


#' Directory of one enrichment result unit, relative to Enrichment/
#'
#' @param method One of "ora", "gsea", "qea", "ssgsea".
#' @param contrast_dir Directory name of the contrast, from
#'   \code{metab_enrichment_contrast_dirs()}. Required for ORA and GSEA. For QEA
#'   and ssGSEA, NULL means the run over all samples.
#' @return Relative path with "/" separators (not created).
metab_enrichment_unit_rel_dir <- function(method, contrast_dir = NULL) {
    method <- as.character(method)
    if (length(method) != 1L || !method %in% names(.METAB_ENRICHMENT_METHOD_DIR)) {
        stop("metab_enrichment_unit_rel_dir(): unknown method '",
             paste(method, collapse = ", "), "'. Expected one of: ",
             paste(names(.METAB_ENRICHMENT_METHOD_DIR), collapse = ", "), ".")
    }
    if (!is.null(contrast_dir)) {
        contrast_dir <- as.character(contrast_dir)
        if (length(contrast_dir) != 1L || is.na(contrast_dir) ||
            !nzchar(contrast_dir) || contrast_dir %in% c(".", "..") ||
            grepl("[/\\\\]", contrast_dir)) {
            stop("metab_enrichment_unit_rel_dir(): '",
                 paste(contrast_dir, collapse = ", "),
                 "' is not a usable directory name. Pass a name from ",
                 "metab_enrichment_contrast_dirs().")
        }
    } else if (method %in% c("ora", "gsea")) {
        stop("metab_enrichment_unit_rel_dir(): ", toupper(method),
             " runs per contrast, so 'contrast_dir' is required.")
    }

    base <- paste(.METAB_ENRICHMENT_METHOD_DIR[[method]], .METAB_ENRICHMENT_DB,
                  sep = "/")
    switch(method,
        ora  = paste(base, "contrasts", contrast_dir, "de", sep = "/"),
        gsea = paste(base, "ranking_by_fc", contrast_dir, sep = "/"),
        if (is.null(contrast_dir)) paste(base, "all_samples", sep = "/")
        else paste(base, "contrasts", contrast_dir, sep = "/")
    )
}


#' Directory of one enrichment result unit
#'
#' @param enr_root The mode's Enrichment directory.
#' @inheritParams metab_enrichment_unit_rel_dir
#' @return Path of the unit directory (not created).
metab_enrichment_unit_dir <- function(enr_root, method, contrast_dir = NULL) {
    file.path(enr_root, metab_enrichment_unit_rel_dir(method, contrast_dir))
}


#' Flat (pre-layout) file name of an enrichment plot
#'
#' Before the unit layout, every enrichment plot sat directly in Enrichment/ and
#' carried its context in its name, e.g. \code{enrichment_ora_dotplot_<contrast>.png}.
#' Consumers that label or order plots by that name (the summary PowerPoint) use
#' this to keep doing so. A file that is already directly in Enrichment/, or is
#' outside it, keeps its own name.
#'
#' @param path Path(s) of plot files.
#' @param enr_root The mode's Enrichment directory.
#' @return Character vector of file names, same length as \code{path}.
metab_enrichment_legacy_plot_name <- function(path, enr_root) {
    to_slash <- function(x) gsub("\\\\", "/", as.character(x))
    root <- sub("/+$", "", to_slash(enr_root))
    method_of <- stats::setNames(names(.METAB_ENRICHMENT_METHOD_DIR),
                                 .METAB_ENRICHMENT_METHOD_DIR)

    vapply(to_slash(path), function(p) {
        prefix <- paste0(root, "/")
        if (!startsWith(p, prefix)) return(basename(p))
        parts <- strsplit(substring(p, nchar(prefix) + 1L), "/", fixed = TRUE)[[1]]
        # <METHOD>/<db>/<scope>/[<contrast>/[de/]]<file>
        if (length(parts) < 4L || !parts[1] %in% names(method_of)) return(basename(p))

        file <- parts[length(parts)]
        stem <- tools::file_path_sans_ext(file)
        ext  <- tools::file_ext(file)
        sfx  <- if (identical(parts[3], "all_samples")) "" else paste0("_", parts[4])
        # The per-library QEA plots were enrichment_qea_<library>.png.
        if (identical(parts[3], "all_samples") &&
            startsWith(stem, .METAB_ENRICHMENT_LIBRARY_PLOT_PREFIX)) {
            stem <- substring(stem, nchar(.METAB_ENRICHMENT_LIBRARY_PLOT_PREFIX) + 1L)
        }
        paste0("enrichment_", method_of[[parts[1]]], "_", stem, sfx,
               if (nzchar(ext)) paste0(".", ext) else "")
    }, character(1), USE.NAMES = FALSE)
}


#' Enrichment files on disk that this stage owns
#'
#' The files only the metabolomics enrichment stage writes. A file is listed
#' only when its whole path under the mode folder has a shape the stage
#' produces AND its name is one the stage gives, so a file another stage or a
#' person put in these folders is never listed, whatever its extension:
#' \itemize{
#'   \item in a unit directory of the layout (see
#'     \code{metab_enrichment_unit_rel_dir()}): \code{results.tsv} and the
#'     method's own plots (\code{.METAB_ENRICHMENT_PLOT_KEYS}), plus
#'     \code{library_<name>.png} in the all-samples QEA unit;
#'   \item directly in \code{Enrichment/}, the flat names used before the unit
#'     layout and by the first-contrast copies still written:
#'     \code{enrichment_<method>_<plot>[_<contrast>].png} for the method's own
#'     plots, and \code{enrichment_qea_<name>.png} (the old per-library QEA
#'     plots carried a free-text library name);
#'   \item directly in \code{Datasets/}:
#'     \code{enrichment_<method>[_<contrast>]_results.tsv}.
#' }
#' The manifest is not listed: the stage removes and rewrites it itself. The
#' mummichog stage's \code{mummichog_pinned/} tree is outside both folders.
#'
#' @param enr_root The mode's Enrichment directory.
#' @param datasets_dir The mode's Datasets directory.
#' @return Character vector of existing file paths (possibly empty).
list_metab_enrichment_owned_files <- function(enr_root, datasets_dir) {
    keys <- .METAB_ENRICHMENT_PLOT_KEYS
    alt  <- function(x) paste0("(", paste(x, collapse = "|"), ")")
    seg  <- "[^/]+"
    unit <- function(method, scope, plot_keys, extra = NULL) {
        names_ok <- c("results\\.tsv", paste0(alt(plot_keys), "\\.png"), extra)
        paste0("^", .METAB_ENRICHMENT_METHOD_DIR[[method]], "/",
               .METAB_ENRICHMENT_DB, "/", scope, "/", alt(names_ok), "$")
    }
    lib_png <- paste0(.METAB_ENRICHMENT_LIBRARY_PLOT_PREFIX, seg, "\\.png")
    flat <- function(method) {
        paste0("^enrichment_", method, "_", alt(keys[[method]]),
               "(_", seg, ")?\\.png$")
    }

    enr_patterns <- c(
        unit("ora",    paste0("contrasts/", seg, "/de"), keys$ora),
        unit("gsea",   paste0("ranking_by_fc/", seg),    keys$gsea),
        unit("qea",    "all_samples",                    keys$qea, extra = lib_png),
        unit("qea",    paste0("contrasts/", seg),        keys$qea),
        unit("ssgsea", "all_samples",                    keys$ssgsea),
        unit("ssgsea", paste0("contrasts/", seg),        keys$ssgsea),
        flat("ora"), flat("gsea"), flat("ssgsea"),
        paste0("^enrichment_qea_", seg, "\\.png$")
    )
    methods  <- paste(names(.METAB_ENRICHMENT_METHOD_DIR), collapse = "|")
    flat_tsv <- paste0("^enrichment_(", methods, ")(_.*)?_results\\.tsv$")

    owned_enr <- character(0)
    if (dir.exists(enr_root)) {
        # all.files: a contrast whose name starts with "." has a hidden unit
        # directory, which list.files() otherwise does not descend into.
        rel <- list.files(enr_root, recursive = TRUE, all.files = TRUE)
        hit <- Reduce(`|`, lapply(enr_patterns, grepl, x = rel),
                      rep(FALSE, length(rel)))
        owned_enr <- file.path(enr_root, rel[hit])
    }
    owned_ds <- character(0)
    if (dir.exists(datasets_dir)) {
        owned_ds <- list.files(datasets_dir, pattern = flat_tsv, full.names = TRUE)
        # list.files() also returns directories that match the pattern.
        owned_ds <- owned_ds[!dir.exists(owned_ds)]
    }
    c(owned_enr, owned_ds)
}


#' Manifest of the enrichment units a run wrote
#'
#' One row per result unit, in the column schema of the RNA-seq enrichment
#' manifest (\code{.empty_enrichment_manifest()}), so both omics describe their
#' units the same way. Only units that produced a result table are listed: a
#' method that skipped a contrast has no row, because "skipped" and "tested,
#' nothing significant" are different facts and the run does not record the
#' first one.
#'
#' \code{status} is "significant" when at least one pathway has
#' FDR < \code{.METAB_ENRICHMENT_SIG_FDR}, otherwise "empty".
#'
#' The columns are the RNA-seq ones, but \code{storage_key} does NOT mean the
#' same thing. In RNA-seq it is the key of the result table inside the
#' pathway-results object (the Shiny payload builder looks tables up by it).
#' Here it is the unit's directory relative to Enrichment/, i.e. where the files
#' are. It locates a unit on disk; it is not a key into any R object, so code
#' written against the RNA-seq manifest must not use it as one.
#'
#' @param units List of units, each a list with \code{method} ("ora", "gsea",
#'   "qea", "ssgsea"), \code{contrast} (the contrast label, or NA for a run over
#'   all samples), \code{rel_dir} (from \code{metab_enrichment_unit_rel_dir()})
#'   and \code{table} (the unit's result data.frame, with an \code{FDR} column).
#' @return data.frame with the manifest columns, one row per unit.
build_metab_enrichment_manifest <- function(units) {
    empty <- .empty_enrichment_manifest()
    if (length(units) == 0) return(empty)

    rows <- lapply(units, function(u) {
        method       <- as.character(u$method)
        has_contrast <- !is.null(u$contrast) && !is.na(u$contrast)
        item         <- if (has_contrast) as.character(u$contrast) else "all_samples"
        is_ora       <- identical(method, "ora")

        fdr   <- if (is.data.frame(u$table)) u$table$FDR else NULL
        n_sig <- if (is.null(fdr)) 0L
                 else sum(!is.na(fdr) & fdr < .METAB_ENRICHMENT_SIG_FDR)

        data.frame(
            analysis      = .METAB_ENRICHMENT_METHOD_DIR[[method]],
            database      = .METAB_ENRICHMENT_DB,
            # ORA and GSEA use the RNA-seq vocabulary: ORA is the direction-
            # agnostic result of a contrast, GSEA is ranked by fold change.
            group         = switch(method,
                                   ora  = "contrasts",
                                   gsea = "fc",
                                   if (has_contrast) "contrasts" else "all_samples"),
            item          = item,
            subgroup      = if (is_ora) "all" else NA_character_,
            group_type    = if (is_ora) "direction" else NA_character_,
            group_label   = if (is_ora) "All" else NA_character_,
            evaluated     = TRUE,
            status        = if (n_sig > 0L) "significant" else "empty",
            n_significant = as.integer(n_sig),
            has_simplify  = FALSE,
            storage_key   = as.character(u$rel_dir),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, c(list(empty), rows))
    rownames(out) <- NULL
    out
}
