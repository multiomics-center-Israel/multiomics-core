# Guard against duplicate top-level function definitions in the proteomics
# domain layer.
#
# R/ files are source()d in filename order (see helper.R / _targets.R), so a
# second top-level definition of the same name silently shadows the first —
# a latent divergence hazard (see issue #137). This test fails loudly if any
# proteomics domain function name is defined in more than one file.
#
# Known pre-existing collision tracked in its own issue (#140):
#   build_de_contrast_summary — 05_de_summary.R vs 07_shiny_export.R, different
#   signatures; the fix is a rename (not a deletion) and touches shiny-export
#   code. Remove it from `known_duplicates` below once #140 is resolved so the
#   guard starts enforcing it.
#
# Scope note: this guard scans only R/domain/proteomics. Because the whole R/
# tree is sourced into one environment, proteomics functions can also be
# shadowed from other directories (e.g. filter_features_optimized in
# domain/rnaseq), and other layers collide among themselves. A codebase-wide
# audit of every such collision is tracked in #142; broaden this scan to the
# full sourced file set once that lands.

test_that("no proteomics domain function is defined more than once", {
    root_dir <- normalizePath(if (dir.exists("R")) "." else "../..")
    prot_dir <- file.path(root_dir, "R", "domain", "proteomics")
    skip_if_not(dir.exists(prot_dir), "proteomics domain directory not found")

    # Pre-existing collisions tracked in their own issues, each pinned to the
    # EXACT set of files it may appear in. A definition of an allowlisted name
    # in any *other* file is still a failure, so a third shadow can't hide.
    known_duplicates <- list(
        build_de_contrast_summary = c("05_de_summary.R", "07_shiny_export.R")  # #140
    )

    # Collect every top-level `name <- function(...)` and the file(s) it's in.
    # Parse-based (not regex) so commented-out or string occurrences don't count.
    # recursive = TRUE to match production: _targets.R / helper.R source R/
    # recursively, so a function added under proteomics/<subdir>/ is live too.
    r_files <- list.files(prot_dir, pattern = "\\.[Rr]$", full.names = TRUE,
                          recursive = TRUE)
    defs <- list()
    for (f in r_files) {
        for (e in parse(f, keep.source = FALSE)) {
            if (is.call(e) && length(e) == 3L &&
                as.character(e[[1L]]) %in% c("<-", "=") &&
                is.name(e[[2L]]) &&
                is.call(e[[3L]]) &&
                identical(as.character(e[[3L]][[1L]]), "function")) {
                nm <- as.character(e[[2L]])
                defs[[nm]] <- c(defs[[nm]], basename(f))
            }
        }
    }

    dup <- defs[vapply(defs, length, integer(1L)) > 1L]
    # Drop only the exact allowlisted collisions; any extra definition keeps it
    # failing. Compare as a multiset (identical() on sorted vectors) rather than
    # setequal(): a second definition inside one allowlisted file changes the
    # count but not the set, and must still be reported as a new shadow.
    for (nm in names(known_duplicates)) {
        if (!is.null(dup[[nm]]) &&
            identical(sort(dup[[nm]]), sort(known_duplicates[[nm]]))) {
            dup[[nm]] <- NULL
        }
    }

    detail <- paste(vapply(names(dup), function(nm) {
        sprintf("  %s: %s", nm, paste(dup[[nm]], collapse = ", "))
    }, character(1L)), collapse = "\n")

    expect(
        length(dup) == 0L,
        sprintf(paste0(
            "Duplicate top-level function definitions in R/domain/proteomics/ ",
            "(a later-sorted file silently shadows the earlier one):\n%s\n",
            "Keep a single definition per name, or add to the tracked allowlist."),
            detail)
    )
})
