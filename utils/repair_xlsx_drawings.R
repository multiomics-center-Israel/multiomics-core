#!/usr/bin/env Rscript
#' Repair .xlsx workbooks that reference drawing parts the archive does not hold
#'
#' openxlsx (through 4.2.8.1) registers a drawing relationship, a vmlDrawing
#' relationship and a `[Content_Types].xml` Override for every worksheet, but
#' only writes those parts when the sheet really carries a drawing or a comment.
#' Workbooks built without `insertImage()` therefore ship with relationships
#' pointing at parts that are not in the zip: Excel offers to repair the file
#' and `openpyxl.load_workbook()` aborts with a missing-part error.
#'
#' This script rewrites such a workbook without those dangling references. Every
#' other part is copied through unchanged, so cell values, styles, named regions
#' and formulas are untouched.
#'
#' New files are pipeline output going forward — `R/core/05_export_excel.R`
#' prunes the references before saving. This script exists for the workbooks
#' that were already handed out.
#'
#' Usage:
#'   Rscript utils/repair_xlsx_drawings.R <file-or-directory> [options]
#'
#'   --in-place      Overwrite the source workbook. WITHOUT this flag the
#'                   repaired copy is written next to it as
#'                   <name>.repaired.xlsx and the original is left alone.
#'   --recursive     When the path is a directory, descend into subdirectories.
#'   --dry-run       Report what would change; write nothing.
#'   --quiet         Only print the final summary.
#'
#' Examples:
#'   Rscript utils/repair_xlsx_drawings.R outputs/run1/Final_results_ALL_P_0.05.xlsx
#'   Rscript utils/repair_xlsx_drawings.R outputs/run1 --recursive --dry-run

if (!requireNamespace("zip", quietly = TRUE)) {
    stop("Package 'zip' is required (it ships with openxlsx). Run renv::restore().")
}

#' Parse this script's command line
#'
#' @param args Character vector of arguments (typically
#'   \code{commandArgs(trailingOnly = TRUE)}).
#' @return List with \code{path} plus the logical flags \code{in_place},
#'   \code{recursive}, \code{dry_run} and \code{quiet}.
parse_repair_args <- function(args) {
    known <- c("--in-place", "--recursive", "--dry-run", "--quiet")
    flags <- args[startsWith(args, "-")]
    unknown <- setdiff(flags, known)
    if (length(unknown) > 0) {
        stop("Unknown option(s): ", paste(unknown, collapse = ", "),
             "\nSupported: ", paste(known, collapse = " "), call. = FALSE)
    }
    positional <- args[!startsWith(args, "-")]
    if (length(positional) != 1L) {
        stop("Expected exactly one path (a .xlsx file or a directory), got ",
             length(positional), ".", call. = FALSE)
    }
    list(
        path      = positional[[1]],
        in_place  = "--in-place" %in% flags,
        recursive = "--recursive" %in% flags,
        dry_run   = "--dry-run" %in% flags,
        quiet     = "--quiet" %in% flags
    )
}

#' Collect the .xlsx files a path refers to
#'
#' @param path A single .xlsx file or a directory containing some.
#' @param recursive Logical; descend into subdirectories when \code{path} is a
#'   directory.
#' @return Character vector of absolute paths. Files already named
#'   \code{*.repaired.xlsx} are skipped so a second run does not repair its own
#'   output.
collect_xlsx_files <- function(path, recursive = FALSE) {
    if (!file.exists(path)) stop("Path does not exist: ", path, call. = FALSE)
    files <- if (dir.exists(path)) {
        list.files(path, pattern = "\\.xlsx$", full.names = TRUE,
                   recursive = recursive)
    } else {
        path
    }
    # Only when scanning a directory: an explicitly named *.repaired.xlsx is a
    # deliberate request and should still be processed.
    if (dir.exists(path)) files <- files[!grepl("\\.repaired\\.xlsx$", files)]
    normalizePath(files, mustWork = FALSE)
}

#' Strip references to drawing parts that are missing from an unpacked workbook
#'
#' Operates on a directory holding an unzipped .xlsx. Only
#' \code{[Content_Types].xml}, the worksheet .rels files and (when a removed
#' relationship was still referenced) the sheet XML are rewritten; a reference
#' whose target part is present in the archive is always kept, so a workbook
#' with a genuine image survives untouched.
#'
#' @param dir Directory containing the unpacked workbook.
#' @return List with \code{rels} (number of relationships removed),
#'   \code{overrides} (Content_Types Overrides removed) and \code{sheet_refs}
#'   (\code{<drawing>}/\code{<legacyDrawing>} elements removed).
strip_missing_drawing_refs <- function(dir) {
    present <- basename(list.files(file.path(dir, "xl", "drawings")))
    n_rels <- 0L
    n_overrides <- 0L
    n_sheet_refs <- 0L

    rels_files <- list.files(file.path(dir, "xl", "worksheets", "_rels"),
                             pattern = "\\.rels$", full.names = TRUE)
    missing_parts <- character(0)

    for (rf in rels_files) {
        txt <- readLines(rf, warn = FALSE)
        one <- paste(txt, collapse = "\n")
        # Relationships are self-closing single elements; split them out so a
        # target can be tested one at a time.
        rels <- regmatches(one, gregexpr("<Relationship\\b[^>]*/>", one, perl = TRUE))[[1]]
        if (length(rels) == 0L) next

        targets <- sub('.*Target="([^"]*)".*', "\\1", rels)
        ids     <- sub('.*Id="([^"]*)".*', "\\1", rels)
        points_at_drawings <- grepl("^\\.\\./drawings/", targets)
        target_files <- basename(targets)
        drop <- points_at_drawings & !(target_files %in% present)
        if (!any(drop)) next

        missing_parts <- c(missing_parts, target_files[drop])
        dropped_ids <- ids[drop]
        for (rel in rels[drop]) one <- sub(rel, "", one, fixed = TRUE)
        n_rels <- n_rels + sum(drop)
        writeLines(one, rf)

        # A rel is only half the reference: the sheet may still point at the id.
        sheet <- file.path(dir, "xl", "worksheets", sub("\\.rels$", "", basename(rf)))
        if (file.exists(sheet) && length(dropped_ids) > 0L) {
            s <- paste(readLines(sheet, warn = FALSE), collapse = "\n")
            before <- s
            for (id in dropped_ids) {
                pat <- sprintf('<(legacyDrawing|drawing)\\b[^>]*r:id="%s"[^>]*/>', id)
                hits <- regmatches(s, gregexpr(pat, s, perl = TRUE))[[1]]
                n_sheet_refs <- n_sheet_refs + length(hits)
                s <- gsub(pat, "", s, perl = TRUE)
            }
            if (!identical(s, before)) writeLines(s, sheet)
        }
    }

    ct_file <- file.path(dir, "[Content_Types].xml")
    if (file.exists(ct_file) && length(missing_parts) > 0L) {
        ct <- paste(readLines(ct_file, warn = FALSE), collapse = "\n")
        overrides <- regmatches(ct, gregexpr("<Override\\b[^>]*/>", ct, perl = TRUE))[[1]]
        part_names <- sub('.*PartName="([^"]*)".*', "\\1", overrides)
        drop <- grepl("^/xl/drawings/", part_names) &
            !(basename(part_names) %in% present)
        if (any(drop)) {
            for (ov in overrides[drop]) ct <- sub(ov, "", ct, fixed = TRUE)
            n_overrides <- sum(drop)
            writeLines(ct, ct_file)
        }
    }

    list(rels = n_rels, overrides = n_overrides, sheet_refs = n_sheet_refs)
}

#' Repair a single .xlsx workbook
#'
#' @param path Path to the workbook to repair.
#' @param in_place Logical; \code{TRUE} overwrites \code{path},
#'   \code{FALSE} (the default) writes \code{<name>.repaired.xlsx} beside it.
#'   Nothing is written when the workbook has no dangling references.
#' @param dry_run Logical; report the counts without writing anything.
#' @return List with \code{path}, \code{out} (the file written, or NA), and the
#'   removal counts from \code{\link{strip_missing_drawing_refs}}.
repair_xlsx_drawings <- function(path, in_place = FALSE, dry_run = FALSE) {
    tmp <- tempfile("xlsx_repair_")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

    utils::unzip(path, exdir = tmp)
    counts <- strip_missing_drawing_refs(tmp)
    changed <- counts$rels > 0L || counts$overrides > 0L || counts$sheet_refs > 0L

    out <- NA_character_
    if (changed && !dry_run) {
        # Write to a sibling temp file first, so an interrupted --in-place run
        # cannot leave a truncated workbook where the deliverable used to be.
        dest <- if (in_place) path else sub("\\.xlsx$", ".repaired.xlsx", path)
        staging <- tempfile("xlsx_repaired_", tmpdir = dirname(dest), fileext = ".xlsx")
        entries <- list.files(tmp, recursive = TRUE, all.files = TRUE, no.. = TRUE)
        # [Content_Types].xml conventionally leads the archive; some readers
        # expect to find it without scanning the whole central directory.
        entries <- c(intersect("[Content_Types].xml", entries),
                     setdiff(entries, "[Content_Types].xml"))
        # "mirror" keeps each entry's path relative to root; "cherry-pick"
        # would flatten xl/worksheets/sheet1.xml to sheet1.xml.
        zip::zip(zipfile = staging, files = entries, root = tmp, mode = "mirror")
        file.rename(staging, dest)
        out <- dest
    }

    c(list(path = path, out = out, changed = changed), counts)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
    opt <- parse_repair_args(args)
    files <- collect_xlsx_files(opt$path, recursive = opt$recursive)
    if (length(files) == 0L) {
        message("No .xlsx files found under: ", opt$path)
        return(invisible(NULL))
    }
    if (opt$in_place && !opt$dry_run) {
        message("--in-place: ", length(files),
                " workbook(s) will be OVERWRITTEN in place.")
    }

    n_changed <- 0L
    for (f in files) {
        res <- tryCatch(
            repair_xlsx_drawings(f, in_place = opt$in_place, dry_run = opt$dry_run),
            error = function(e) {
                message("FAILED  ", f, ": ", conditionMessage(e))
                NULL
            }
        )
        if (is.null(res)) next
        if (res$changed) n_changed <- n_changed + 1L
        if (!opt$quiet) {
            message(sprintf(
                "%-8s %s  (rels: %d, overrides: %d, sheet refs: %d)%s",
                if (!res$changed) "clean" else if (opt$dry_run) "would fix" else "fixed",
                basename(f), res$rels, res$overrides, res$sheet_refs,
                if (!is.na(res$out)) paste0(" -> ", basename(res$out)) else ""
            ))
        }
    }

    message(sprintf("%d/%d workbook(s) %s.", n_changed, length(files),
                    if (opt$dry_run) "would be repaired" else "repaired"))
    invisible(NULL)
}

# Run only when executed as a script; sourcing the file just defines the
# functions, so the repair can also be driven from an R session.
if (sys.nframe() == 0L && !interactive()) {
    main()
}
