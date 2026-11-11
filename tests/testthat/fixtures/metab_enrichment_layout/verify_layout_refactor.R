# tests/testthat/fixtures/metab_enrichment_layout/verify_layout_refactor.R
#
# Before/after check for the metabolomics enrichment layout refactor.
#
# It proves on a real run what the unit tests prove on stubs: that moving the
# enrichment files changed where they are and nothing about what they contain.
# It is run by hand, in the pipeline's own R environment; testthat does not run
# it. It reads results, so it never changes a run.
#
# How to run
#
# This file exists only on the refactor branch, and the "before" run is made on
# `main`, where it is absent. So do not source it from the working tree: put one
# copy OUTSIDE the repo and use that same copy for both runs. It calls no
# function of the pipeline, so one copy serves both.
#
#   0. Once, from the repo root (any branch):
#
#        git fetch origin refactor/metab-enrichment-layout
#        git show origin/refactor/metab-enrichment-layout:tests/testthat/fixtures/metab_enrichment_layout/verify_layout_refactor.R > <outside>/verify_layout_refactor.R
#
#      and choose an empty folder outside the repo for the snapshots.
#
#   1. Check out `main`. Point the config at a NEW, empty output folder and run
#      the pipeline. Then, in R:
#
#        source("<outside>/verify_layout_refactor.R")
#        snapshot_enrichment_layout("before",
#            mode_dir = "<project>/<out>/Results_<name>_<round>/metabolomics",
#            snap_dir = "<outside>/snapshots")
#
#   2. Check out the refactor branch. Point the SAME config at ANOTHER new,
#      empty output folder and run the pipeline. Then, in R:
#
#        source("<outside>/verify_layout_refactor.R")
#        snapshot_enrichment_layout("after",
#            mode_dir = "<the new .../metabolomics>",
#            snap_dir = "<outside>/snapshots")
#        compare_enrichment_layout(snap_dir = "<outside>/snapshots")
#
# Instead of switching branches in one checkout you can keep both side by side:
#
#        git worktree add ../multiomics-core-main main
#
# and run step 1 in ../multiomics-core-main and step 2 in this checkout. Each
# worktree needs its own renv library restored (renv::restore()), and the
# script is still sourced from the outside copy.
#
# Use a config with enrichment enabled and at least two contrasts; two GMT
# files exercise the per-library QEA plots as well. Use the same {targets}
# store settings in both runs, and pass `store =` when the run used a store
# other than the default (run.R prints the one it uses).
#
# DATA SAFETY: the snapshots hold result tables of the run. `snap_dir` has no
# default so that they are never written into the repo by accident. Keep it out
# of git, and report which checks passed, not the contents.

.el_md5 <- function(paths) unname(tools::md5sum(paths))

# The output folders this script looks into, in the spelling the comparison
# uses. The pipeline does not spell them one way everywhere ("Datasets" and
# "datasets" both occur), and on a file system that ignores case the folder
# keeps whichever spelling created it first. Matching the spelling exactly
# therefore found no enrichment table on such a system.
.EL_TOP_DIRS <- c("Enrichment", "Datasets", "mummichog_pinned")

# Relative paths with "/" separators and the top-level folder in the canonical
# spelling, whatever its case on disk. Only that first component is rewritten.
.el_canon_paths <- function(paths) {
    paths <- gsub("\\\\", "/", as.character(paths))
    for (d in .EL_TOP_DIRS) {
        paths <- sub(paste0("^", d, "(/|$)"), paste0(d, "\\1"), paths,
                     ignore.case = TRUE)
    }
    paths
}

# Which relative paths are enrichment result tables / plots of either layout.
.el_is_enrichment_tsv <- function(paths) {
    grepl("^Enrichment/.*\\.tsv$|^Datasets/enrichment_[^/]*\\.tsv$",
          .el_canon_paths(paths))
}
.el_is_enrichment_png <- function(paths) {
    grepl("^Enrichment/.*\\.png$", .el_canon_paths(paths))
}

# What is on disk under a mode folder: every file (canonical relative path) and
# the content hash of the enrichment tables and of the mummichog files. The
# hashes are read through the path as it is on disk and named by the canonical
# one.
.el_inventory <- function(mode_dir) {
    # all.files: a contrast whose name starts with "." has a hidden directory.
    on_disk <- sort(list.files(mode_dir, recursive = TRUE, all.files = TRUE))
    canon   <- .el_canon_paths(on_disk)
    md5_of  <- function(sel) stats::setNames(
        .el_md5(file.path(mode_dir, on_disk[sel])), canon[sel])
    list(files         = canon,
         tsv_md5       = md5_of(.el_is_enrichment_tsv(canon)),
         mummichog_md5 = md5_of(grepl("^mummichog_pinned/", canon)))
}

# TRUE when a snapshot hashed every enrichment table its file list shows, and
# there is at least one. A snapshot with none proves nothing about the tables,
# so the table checks must fail on it rather than compare two empty sets.
.el_tsv_collected <- function(snap) {
    listed <- .el_canon_paths(snap$files)
    listed <- listed[.el_is_enrichment_tsv(listed)]
    hashed <- .el_canon_paths(names(snap$tsv_md5))
    length(hashed) > 0 && setequal(hashed, listed)
}

# Plot data and labels compared by value, exactly: no numeric tolerance, and
# attributes count. identical() is the wrong test here because it compares
# environments and other references by address, and those are never the same
# object in two R sessions.
.el_same_plots <- function(before, after) {
    isTRUE(all.equal(before, after, tolerance = 0, check.attributes = TRUE))
}

# The enrichment result with every `files` element removed, at any depth: the
# paths are what the refactor changes, and `files` occurs at the top level and
# inside each *_by_contrast entry. `plots` are taken out here and compared
# separately (see .el_plots), because a ggplot object carries environments that
# never compare identical between two R sessions.
.el_results <- function(x) {
    if (inherits(x, "ggplot") || is.data.frame(x) || !is.list(x)) return(x)
    keep <- !(names(x) %in% c("files", "plots"))
    if (is.null(names(x))) keep <- rep(TRUE, length(x))
    lapply(x[keep], .el_results)
}

# What each plot was drawn from: the data frame and the labels given to
# ggplot(). Same recursion, keeping only the `plots` branches.
.el_plots <- function(x) {
    if (inherits(x, "ggplot")) {
        return(tryCatch(list(data = x$data, labels = unclass(x$labels)),
                        error = function(e) list(class = class(x))))
    }
    if (is.data.frame(x) || !is.list(x)) return(NULL)
    out <- lapply(x[names(x) != "files"], .el_plots)
    out[!vapply(out, is.null, logical(1))]
}

#' Record what one run produced
#'
#' @param label "before" or "after".
#' @param mode_dir The run's metabolomics output folder.
#' @param snap_dir Folder to write the snapshot to. No default: choose one
#'   outside the repo.
#' @param store The {targets} store the run used.
#' @return The snapshot, invisibly.
snapshot_enrichment_layout <- function(label, mode_dir, snap_dir,
                                       store = targets::tar_config_get("store")) {
    stopifnot(label %in% c("before", "after"), dir.exists(mode_dir))
    dir.create(snap_dir, showWarnings = FALSE, recursive = TRUE)
    res <- targets::tar_read_raw("metab_enrichment_res", store = store)
    if (is.null(res)) {
        stop("metab_enrichment_res is NULL in this store: enrichment did not ",
             "run. Enable modes.metabolomics.enrichment.run_enrichment.")
    }

    inv <- .el_inventory(mode_dir)
    if (length(inv$tsv_md5) == 0) {
        stop("No enrichment result table was found under ", mode_dir,
             ". Check that this is the run's metabolomics folder and that ",
             "enrichment ran; a snapshot without tables cannot be compared.")
    }

    rds <- list.files(mode_dir, pattern = "^shiny_payload_metabolomics_.*\\.rds$",
                      full.names = TRUE)
    payload <- if (length(rds) == 1) readRDS(rds) else NULL
    payload_tables <- if (is.null(payload)) NULL else c(
        payload[grep("^enrichment", names(payload), value = TRUE)],
        list(mummichog_tables = lapply(payload$mummichog, function(s) s$table)))

    pptx <- list.files(dirname(mode_dir), pattern = "^metabolomics_summary\\.pptx$",
                       full.names = TRUE)
    slide_titles <- if (length(pptx) == 1 && requireNamespace("officer", quietly = TRUE)) {
        s <- officer::pptx_summary(officer::read_pptx(pptx))
        s <- s[order(s$slide_id), ]
        s$text[grepl("^Enrichment:", s$text)]
    } else NULL

    snap <- list(
        results        = .el_results(res),
        plots          = .el_plots(res),
        result_names   = names(res),
        by_contrast_names = lapply(
            res[grep("_by_contrast$", names(res), value = TRUE)],
            function(b) lapply(b, names)),
        files          = inv$files,
        tsv_md5        = inv$tsv_md5,
        mummichog_md5  = inv$mummichog_md5,
        payload_names  = names(payload),
        payload_tables = payload_tables,
        slide_titles   = slide_titles
    )
    out <- file.path(snap_dir, paste0(label, ".rds"))
    saveRDS(snap, out)
    message("Saved ", out)
    invisible(snap)
}

# Where a flat enrichment TSV of the old layout lives in the unit layout. The
# contrast token in the old file name is the contrast's directory name.
.el_tsv_new_path <- function(old) {
    m <- regmatches(basename(old), regexec(
        "^enrichment_(qea|ssgsea|ora|gsea)(?:_(.*))?_results\\.tsv$",
        basename(old), perl = TRUE))[[1]]
    if (length(m) == 0) return(NA_character_)
    method <- m[2]; ct <- m[3]
    top <- c(qea = "QEA", ssgsea = "ssGSEA", ora = "ORA", gsea = "GSEA")[[method]]
    base <- paste("Enrichment", top, "combined", sep = "/")
    if (!nzchar(ct)) {
        # Global QEA / ssGSEA moved; the first-contrast ORA / GSEA copy did not.
        if (method %in% c("ora", "gsea")) return(old)
        return(paste(base, "all_samples", "results.tsv", sep = "/"))
    }
    switch(method,
        ora  = paste(base, "contrasts", ct, "de", "results.tsv", sep = "/"),
        gsea = paste(base, "ranking_by_fc", ct, "results.tsv", sep = "/"),
        paste(base, "contrasts", ct, "results.tsv", sep = "/"))
}

# Where a flat enrichment PNG of the old layout lives in the unit layout, given
# the contrast directory names of the run (longest first, so a contrast whose
# name ends another's is not mistaken for it).
.el_png_new_path <- function(old, contrast_dirs) {
    b <- sub("\\.png$", "", basename(old))
    m <- regmatches(b, regexec("^enrichment_(qea|ssgsea|ora|gsea)_(.*)$", b))[[1]]
    if (length(m) == 0) return(NA_character_)
    method <- m[2]; rest <- m[3]
    top <- c(qea = "QEA", ssgsea = "ssGSEA", ora = "ORA", gsea = "GSEA")[[method]]
    base <- paste("Enrichment", top, "combined", sep = "/")
    for (ct in contrast_dirs[order(-nchar(contrast_dirs))]) {
        sfx <- paste0("_", ct)
        if (endsWith(rest, sfx)) {
            key <- substr(rest, 1, nchar(rest) - nchar(sfx))
            return(switch(method,
                ora  = paste(base, "contrasts", ct, "de", paste0(key, ".png"), sep = "/"),
                gsea = paste(base, "ranking_by_fc", ct, paste0(key, ".png"), sep = "/"),
                paste(base, "contrasts", ct, paste0(key, ".png"), sep = "/")))
        }
    }
    # No contrast suffix: a global QEA / ssGSEA plot moved; a first-contrast
    # ORA / GSEA copy did not.
    if (method %in% c("ora", "gsea")) return(old)
    # Any global QEA plot other than these two was a per-library barplot.
    if (method == "qea" && !rest %in% c("barplot", "lollipop")) {
        rest <- paste0("library_", rest)
    }
    paste(base, "all_samples", paste0(rest, ".png"), sep = "/")
}

#' Compare the two snapshots
#'
#' @param snap_dir The folder both snapshots were written to.
#' @return TRUE if every check passed, invisibly, with the per-check results in
#'   the attribute "checks". Prints one line per check.
compare_enrichment_layout <- function(snap_dir) {
    b <- readRDS(file.path(snap_dir, "before.rds"))
    a <- readRDS(file.path(snap_dir, "after.rds"))
    check <- function(name, cond, hint = "") {
        pass <- isTRUE(cond)
        cat(sprintf("[%s] %s\n", if (pass) " OK " else "FAIL", name))
        if (!pass && nzchar(hint)) cat("       ", hint, "\n", sep = "")
        pass
    }
    # Names of what differs, never values.
    differing <- function(x, y, same = identical) {
        nms <- union(names(x), names(y))
        bad <- nms[!vapply(nms, function(n) same(x[[n]], y[[n]]), logical(1))]
        if (length(bad)) paste("differs in:", paste(bad, collapse = ", ")) else ""
    }

    # Paths in one spelling on both sides, also for snapshots taken before the
    # script canonicalised them.
    b_files <- .el_canon_paths(b$files)
    a_files <- .el_canon_paths(a$files)
    b_tsv <- stats::setNames(unname(b$tsv_md5), .el_canon_paths(names(b$tsv_md5)))
    a_tsv <- stats::setNames(unname(a$tsv_md5), .el_canon_paths(names(a$tsv_md5)))
    b_png <- b_files[.el_is_enrichment_png(b_files)]
    a_png <- a_files[.el_is_enrichment_png(a_files)]

    tsv_ok <- .el_tsv_collected(b) && .el_tsv_collected(a)
    tsv_hint <- if (tsv_ok) "" else paste(
        "a snapshot holds no enrichment table, or not all the tables its file",
        "list shows. Take both snapshots again with this version of the script;",
        "that only reads the existing run folders.")
    tsv_new <- vapply(names(b_tsv), .el_tsv_new_path, character(1))

    # Contrast directory names, read off the unit paths of the new layout.
    contrast_dirs <- unique(basename(dirname(sub(
        "/de/results\\.tsv$", "/results.tsv",
        grep("/(contrasts|ranking_by_fc)/", names(a_tsv), value = TRUE)))))
    png_new <- vapply(b_png, .el_png_new_path, character(1),
                      contrast_dirs = contrast_dirs)

    passed <- c(
        results = check("1. every result in metab_enrichment_res is identical (all fields except `files` and `plots`, at every level)",
              identical(b$results, a$results), differing(b$results, a$results)),
        plots = check("2. every plot was drawn from equal data with equal labels (exact values, attributes included)",
              .el_same_plots(b$plots, a$plots),
              differing(b$plots, a$plots, same = .el_same_plots)),
        fields = check("3. metab_enrichment_res has the same fields, top level and per contrast",
              identical(b$result_names, a$result_names) &&
                  identical(b$by_contrast_names, a$by_contrast_names)),
        tsv_content = check("4. every enrichment table of the old layout has a file with identical content in the new one",
              tsv_ok && !anyNA(tsv_new) && all(tsv_new %in% names(a_tsv)) &&
                  identical(unname(b_tsv), unname(a_tsv[tsv_new])),
              if (tsv_ok) "an old .tsv has no counterpart, or its content changed"
              else tsv_hint),
        tsv_complete = check("5. the new layout holds no enrichment table beyond those",
              tsv_ok && setequal(tsv_new, names(a_tsv)),
              if (tsv_ok) paste("unexpected:",
                                paste(setdiff(names(a_tsv), tsv_new), collapse = ", "))
              else tsv_hint),
        png = check("6. every enrichment plot of the old layout exists in the new one, and no others",
              length(b_png) > 0 && !anyNA(png_new) && setequal(png_new, a_png),
              if (length(b_png) == 0) "no enrichment plot was found in the 'before' run"
              else paste("missing:", paste(setdiff(png_new, a_png), collapse = ", "),
                         "| unexpected:", paste(setdiff(a_png, png_new), collapse = ", "))),
        manifest = check("7. Enrichment/enrichment_manifest.csv was written, and no temporary file is left",
              "Enrichment/enrichment_manifest.csv" %in% a_files &&
                  !("Enrichment/enrichment_manifest.csv.tmp" %in% a_files)),
        mummichog = check("8. mummichog_pinned/ holds the same files",
              identical(names(b$mummichog_md5), names(a$mummichog_md5))),
        payload_keys = check("9. Shiny payload has the same keys",
              identical(b$payload_names, a$payload_names),
              "NULL on both sides means the payload .rds was not found"),
        payload_tables = check("10. Shiny enrichment keys and mummichog tables are identical",
              !is.null(b$payload_tables) &&
                  identical(b$payload_tables, a$payload_tables),
              differing(b$payload_tables, a$payload_tables)),
        slides = check("11. PPTX enrichment slides: same titles in the same order",
              !is.null(b$slide_titles) && identical(b$slide_titles, a$slide_titles),
              "NULL means the deck or the 'officer' package was not found")
    )
    cat(sprintf("\n%d of %d checks passed.\n", sum(passed), length(passed)))
    if (!identical(b$mummichog_md5, a$mummichog_md5)) {
        cat("Note: some mummichog_pinned/ files differ in content between the two runs.\n",
            "The refactor does not touch that stage; compare two runs of `main` to see\n",
            "whether those files are reproducible on your setup.\n", sep = "")
    }
    cat("Not covered here: open both HTML reports and compare the Enrichment tabs.\n")
    invisible(structure(all(passed), checks = passed))
}
