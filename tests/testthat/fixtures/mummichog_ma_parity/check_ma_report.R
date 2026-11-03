# tests/testthat/fixtures/mummichog_ma_parity/check_ma_report.R
#
# Check-only validation of a report run built by build_ma_report_run.R: compares
# the new run's report and Shiny payload with the original run's. Reads both,
# writes only <dest>/report_check.txt. Run from the multiomics-core root:
#
#   Rscript tests/testthat/fixtures/mummichog_ma_parity/check_ma_report.R \
#       --orig /path/to/Results_<run> --dest /path/to/Results_<run>_MA
#
# build_ma_report_run.R sources this file and calls check_ma_report() itself.
#
# Report structure is compared on heading TEXT with the section numbers
# stripped (number_sections renumbers everything after a section that changed
# size) and without the title block (its date is the render date):
#   * outside the Mummichog section: every original heading must be present,
#     and no section may hold fewer figures, interactive widgets / tables or
#     AI-commentary blocks than before;
#   * the Mummichog section: the original per-polarity tabs
#     ("<contrast>, positive ions" / ", negative ions") must map onto exactly
#     the new per-contrast tabs, and each new tab must carry Summary, ORA plot,
#     GSEA plot and Advanced / Downloads.
# Shiny payload: same keys as the original; every key the original filled is
# filled; the embedded final Excel files are byte-identical to the original's.

.parse_report <- function(html) {
  m <- gregexpr("(?s)<h([1-4])([^>]*)>(.*?)</h\\1>", html, perl = TRUE)[[1]]
  if (m[1] < 0) return(data.frame(level = integer(), text = character(), start = integer()))
  tags <- regmatches(html, list(m))[[1]]
  level <- as.integer(sub("^<h([1-4]).*", "\\1", tags))
  attrs <- sub("^<h[1-4]([^>]*)>.*", "\\1", tags)
  text <- trimws(gsub("<[^>]+>", "", sub("(?s)^<h[1-4][^>]*>(.*)</h[1-4]>$", "\\1", tags,
                                         perl = TRUE)))
  text <- gsub("\\s+", " ", text)
  text <- trimws(sub("^[0-9]+(\\.[0-9]+)*\\s+", "", text))      # section numbers
  keep <- !grepl("class=\"(title|subtitle|author|date)", attrs)  # title block
  data.frame(level = level, text = text, start = as.integer(m), stringsAsFactors = FALSE)[keep, ]
}

.section_counts <- function(html, h) {
  h1 <- h[h$level == 1, ]
  if (nrow(h1) == 0) return(data.frame())
  ends <- c(h1$start[-1] - 1L, nchar(html))
  cnt <- function(seg, pat) lengths(regmatches(seg, gregexpr(pat, seg, perl = TRUE)))
  segs <- substring(html, h1$start, ends)
  data.frame(section = h1$text, img = cnt(segs, "<img\\b"),
             widget = cnt(segs, "html-widget"),
             commentary = cnt(segs, "commentary-block"), stringsAsFactors = FALSE)
}

.mummichog_tabs <- function(h) {
  top <- which(h$level == 1 & grepl("^Mummichog Pathway Analysis", h$text))
  if (length(top) == 0) return(list())
  nxt <- which(h$level == 1 & seq_len(nrow(h)) > top[1])
  last <- if (length(nxt)) nxt[1] - 1L else nrow(h)
  blk <- h[(top[1] + 1L):last, ]
  tabs <- list(); cur <- NULL
  for (i in seq_len(nrow(blk))) {
    if (blk$level[i] == 2) {
      cur <- blk$text[i]
      if (!identical(cur, "Methodology")) tabs[[cur]] <- character()
    } else if (blk$level[i] == 3 && !is.null(cur) && !identical(cur, "Methodology")) {
      tabs[[cur]] <- c(tabs[[cur]], blk$text[i])
    }
  }
  tabs
}

#' Compare a rebuilt report run with its original
#'
#' @param orig  Original run directory.
#' @param dest  New run directory.
#' @param check Function(name, ok, detail) recording one PASS/FAIL line.
#' @return Invisibly, the lines written to `<dest>/report_check.txt`.
check_ma_report <- function(orig, dest, check) {
  lines <- character()
  read_html <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
  f_old <- file.path(orig, "report_metabolomics.html")
  f_new <- file.path(dest, "report_metabolomics.html")
  if (!file.exists(f_old) || !file.exists(f_new)) {
    check("both reports found", FALSE, paste(f_old, f_new, sep = " | "))
    return(invisible(lines))
  }
  old <- read_html(f_old); new <- read_html(f_new)
  h_old <- .parse_report(old); h_new <- .parse_report(new)

  # -- outside the Mummichog section ------------------------------------------
  in_mm <- function(h) {
    top <- which(h$level == 1 & grepl("^Mummichog Pathway Analysis", h$text))
    out <- rep(FALSE, nrow(h))
    if (length(top)) {
      nxt <- which(h$level == 1 & seq_len(nrow(h)) > top[1])
      out[top[1]:(if (length(nxt)) nxt[1] - 1L else nrow(h))] <- TRUE
    }
    out
  }
  miss <- setdiff(h_old$text[!in_mm(h_old)], h_new$text[!in_mm(h_new)])
  check("every non-Mummichog heading of the original report is present",
        length(miss) == 0, paste(miss, collapse = " | "))
  lines <- c(lines, "# non-Mummichog headings missing", miss)

  c_old <- .section_counts(old, h_old); c_new <- .section_counts(new, h_new)
  cmp <- merge(c_old, c_new, by = "section", suffixes = c(".old", ".new"), all.x = TRUE)
  cmp <- cmp[!grepl("^Mummichog Pathway Analysis", cmp$section), ]
  for (col in c("img.new", "widget.new", "commentary.new")) {
    cmp[[col]][is.na(cmp[[col]])] <- -1L                       # section gone
  }
  short <- cmp$section[cmp$img.new < cmp$img.old |
                         cmp$widget.new < cmp$widget.old |
                         cmp$commentary.new < cmp$commentary.old]
  check("no non-Mummichog section lost figures, tables/widgets or AI commentary",
        length(short) == 0, paste(short, collapse = " | "))
  lines <- c(lines, "# per-section counts (old vs new)",
             utils::capture.output(print(cmp, row.names = FALSE)))

  # -- the Mummichog section --------------------------------------------------
  t_old <- .mummichog_tabs(h_old); t_new <- .mummichog_tabs(h_new)
  base_old <- unique(sub(", (positive|negative) ions$", "", names(t_old)))
  check("Mummichog: the new per-contrast tabs cover exactly the original contrasts",
        setequal(base_old, names(t_new)) && length(t_new) > 0,
        sprintf("original %s | new %s", paste(base_old, collapse = ", "),
                paste(names(t_new), collapse = ", ")))
  check("Mummichog: no per-polarity tab left", !any(grepl("(positive|negative) ions$", names(t_new))))
  want <- c("Summary", "ORA plot", "GSEA plot", "Advanced / Downloads")
  lacking <- names(t_new)[!vapply(t_new, function(s) all(want %in% s), logical(1))]
  check("Mummichog: every contrast tab has Summary, ORA plot, GSEA plot, Advanced / Downloads",
        length(lacking) == 0, paste(lacking, collapse = ", "))
  lines <- c(lines, "# Mummichog tabs (original)", names(t_old),
             "# Mummichog tabs (new)", names(t_new))

  # -- Shiny payload ------------------------------------------------------------
  p_old <- list.files(file.path(orig, "metabolomics"), "^shiny_payload_metabolomics_.*\\.rds$",
                      full.names = TRUE)
  p_new <- list.files(file.path(dest, "metabolomics"), "^shiny_payload_metabolomics_.*\\.rds$",
                      full.names = TRUE)
  if (length(p_old) == 1L && length(p_new) == 1L) {
    a <- readRDS(p_old); b <- readRDS(p_new)
    check("Shiny payload: same keys as the original", identical(sort(names(a)), sort(names(b))))
    lost <- names(a)[!vapply(names(a), function(k) is.null(a[[k]]) || !is.null(b[[k]]), logical(1))]
    check("Shiny payload: every key the original filled is filled", length(lost) == 0,
          paste(lost, collapse = ", "))
    for (k in c("all_final_xlsx", "de_final_xlsx")) {
      check(sprintf("Shiny payload: %s identical to the original", k), identical(a[[k]], b[[k]]),
            if (is.null(b[[k]])) "NULL in the new payload" else "")
    }
    check("Shiny payload: one mummichog entry per new report tab",
          setequal(gsub("_vs_", " vs ", names(b$mummichog), fixed = TRUE), names(t_new)))
  } else {
    check("one Shiny payload in each run", FALSE,
          sprintf("%d original, %d new", length(p_old), length(p_new)))
  }
  writeLines(lines, file.path(dest, "report_check.txt"))
  invisible(lines)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- function(flag) {
    i <- match(flag, args)
    if (is.na(i) || i == length(args)) stop("usage: check_ma_report.R --orig <run> --dest <new run>")
    args[[i + 1L]]
  }
  results <- list()
  check <- function(name, ok, detail = "") {
    results[[name]] <<- isTRUE(ok)
    cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", name,
                if (nzchar(detail)) paste0(" -- ", detail) else ""))
  }
  check_ma_report(normalizePath(opt("--orig"), mustWork = TRUE),
                  normalizePath(opt("--dest"), mustWork = TRUE), check)
  if (!all(unlist(results))) quit(status = 1)
}
