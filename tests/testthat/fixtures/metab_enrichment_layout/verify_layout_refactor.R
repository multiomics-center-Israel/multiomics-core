# tests/testthat/fixtures/metab_enrichment_layout/verify_layout_refactor.R
#
# Before/after check for the metabolomics enrichment layout refactor.
#
# It proves on a real run what the unit tests prove on stubs: that moving the
# enrichment files changed where they are and nothing about what they contain.
# It is run by hand, in the pipeline's own R environment; testthat does not run
# it. It reads results, so it never changes a run.
#
# How to run (from the repo root):
#
#   1. Check out `main`. Point the config at a NEW, empty output folder and run
#      the pipeline. Then:
#
#        source("tests/testthat/fixtures/metab_enrichment_layout/verify_layout_refactor.R")
#        snapshot_enrichment_layout("before",
#            mode_dir = "<project>/<out>/Results_<name>_<round>/metabolomics",
#            snap_dir = "<a folder OUTSIDE the repo>")
#
#   2. Check out the refactor branch. Point the SAME config at ANOTHER new,
#      empty output folder and run the pipeline. Then:
#
#        source("tests/testthat/fixtures/metab_enrichment_layout/verify_layout_refactor.R")
#        snapshot_enrichment_layout("after",
#            mode_dir = "<the new .../metabolomics>",
#            snap_dir = "<the same folder as in step 1>")
#        compare_enrichment_layout(snap_dir = "<the same folder>")
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

    all_files <- sort(list.files(mode_dir, recursive = TRUE))
    is_tsv <- grepl("^Enrichment/.*\\.tsv$|^Datasets/enrichment_[^/]*\\.tsv$", all_files)
    is_png <- grepl("^Enrichment/.*\\.png$", all_files)
    is_mmc <- grepl("^mummichog_pinned/", all_files)
    md5_of <- function(sel) stats::setNames(.el_md5(file.path(mode_dir, all_files[sel])),
                                            all_files[sel])

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
        files          = all_files,
        tsv_md5        = md5_of(is_tsv),
        png_files      = all_files[is_png],
        mummichog_md5  = md5_of(is_mmc),
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
    paste(base, "all_samples", paste0(rest, ".png"), sep = "/")
}

#' Compare the two snapshots
#'
#' @param snap_dir The folder both snapshots were written to.
#' @return TRUE if every check passed, invisibly. Prints one line per check.
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
    differing <- function(x, y) {
        nms <- union(names(x), names(y))
        bad <- nms[!vapply(nms, function(n) identical(x[[n]], y[[n]]), logical(1))]
        if (length(bad)) paste("differs in:", paste(bad, collapse = ", ")) else ""
    }

    tsv_new <- vapply(names(b$tsv_md5), .el_tsv_new_path, character(1))
    # Contrast directory names, read off the unit paths of the new layout.
    contrast_dirs <- unique(basename(dirname(sub(
        "/de/results\\.tsv$", "/results.tsv",
        grep("/(contrasts|ranking_by_fc)/", names(a$tsv_md5), value = TRUE)))))
    png_new <- vapply(b$png_files, .el_png_new_path, character(1),
                      contrast_dirs = contrast_dirs)

    passed <- c(
        check("1. every result in metab_enrichment_res is identical (all fields except `files` and `plots`, at every level)",
              identical(b$results, a$results), differing(b$results, a$results)),
        check("2. every plot was drawn from identical data with identical labels",
              identical(b$plots, a$plots), differing(b$plots, a$plots)),
        check("3. metab_enrichment_res has the same fields, top level and per contrast",
              identical(b$result_names, a$result_names) &&
                  identical(b$by_contrast_names, a$by_contrast_names)),
        check("4. every enrichment table of the old layout has a file with identical content in the new one",
              !anyNA(tsv_new) && all(tsv_new %in% names(a$tsv_md5)) &&
                  identical(unname(b$tsv_md5), unname(a$tsv_md5[tsv_new])),
              "an old .tsv has no counterpart, or its content changed"),
        check("5. the new layout holds no enrichment table beyond those",
              setequal(tsv_new, names(a$tsv_md5)),
              paste("unexpected:", paste(setdiff(names(a$tsv_md5), tsv_new), collapse = ", "))),
        check("6. every enrichment plot of the old layout exists in the new one, and no others",
              !anyNA(png_new) && setequal(png_new, a$png_files),
              paste("missing:", paste(setdiff(png_new, a$png_files), collapse = ", "),
                    "| unexpected:", paste(setdiff(a$png_files, png_new), collapse = ", "))),
        check("7. Enrichment/enrichment_manifest.csv was written, and no temporary file is left",
              "Enrichment/enrichment_manifest.csv" %in% a$files &&
                  !("Enrichment/enrichment_manifest.csv.tmp" %in% a$files)),
        check("8. mummichog_pinned/ holds the same files",
              identical(names(b$mummichog_md5), names(a$mummichog_md5))),
        check("9. Shiny payload has the same keys",
              identical(b$payload_names, a$payload_names),
              "NULL on both sides means the payload .rds was not found"),
        check("10. Shiny enrichment keys and mummichog tables are identical",
              !is.null(b$payload_tables) &&
                  identical(b$payload_tables, a$payload_tables),
              differing(b$payload_tables, a$payload_tables)),
        check("11. PPTX enrichment slides: same titles in the same order",
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
    invisible(all(passed))
}
