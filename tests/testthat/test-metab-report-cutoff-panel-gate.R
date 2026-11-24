# The metabolomics report lets the cutoff panel recompute DE calls, so they
# change when the user switches between adjusted and raw p or moves a cutoff.
# At the configured defaults those calls must still equal the pipeline's
# pass.<contrast>, which is p < cutoff on the exact logFC (03_differential.R).
# Three things keep them equal: no passFlag gate, the exact logFC rather than
# one rebuilt from the 3-significant-digit linearFC, and a strict p < cutoff.
# Values below are synthetic.

METAB_ROOT <- normalizePath(if (dir.exists("R")) "." else "../..",
                            mustWork = FALSE)
METAB_TEMPLATE <- file.path(METAB_ROOT, "R", "pipeline", "metabolomics",
                            "templates", "report_metabolomics.Rmd")

# Evaluate the first complete top-level `<name> <- function(...)` in the Rmd.
load_template_fn <- function(name) {
  lines <- readLines(METAB_TEMPLATE, warn = FALSE)
  start <- grep(paste0("^", name, " <- function"), lines)[1]
  if (is.na(start)) stop(name, "() not found in ", METAB_TEMPLATE)
  for (end in start:length(lines)) {
    ex <- tryCatch(parse(text = lines[start:end]), error = function(e) NULL)
    if (!is.null(ex) && length(ex) == 1L) return(eval(ex[[1]]))
  }
  stop("could not parse ", name, "() from ", METAB_TEMPLATE)
}

test_that("metabolomics plot registrations carry no passFlag column", {
  code <- readLines(METAB_TEMPLATE, warn = FALSE)
  code <- code[!grepl("^\\s*#", code)]
  expect_false(any(grepl("\\$passFlag\\s*<-", code)))
  expect_false(any(grepl("passFlag\\s*=", code)))
})

test_that("both metabolomics registrations ask for a strict p < cutoff", {
  code <- paste(readLines(METAB_TEMPLATE, warn = FALSE), collapse = "\n")
  calls <- regmatches(code, gregexpr("cutoff_register_plot\\((.|\n)*?\\)\\)", code))[[1]]
  expect_length(calls, 2L)
  expect_true(all(grepl("strict_p = TRUE", calls, fixed = TRUE)))
})

test_that("cutoff_register_plot emits strictP, default false", {
  df <- data.frame(name = "F1", logFC = 1, pval = 0.05, adjPval = 0.05)
  expect_true(grepl("strictP: false", cutoff_register_plot("v1", df), fixed = TRUE))
  expect_true(grepl("strictP: true",
                    cutoff_register_plot("v2", df, strict_p = TRUE), fixed = TRUE))
  panel <- cutoff_panel_html()
  expect_true(grepl("strictP ? pval < pCut : pval <= pCut", panel, fixed = TRUE))
  expect_true(grepl("reg.strictP === true", panel, fixed = TRUE))
})

test_that("metab_panel_logfc returns the exact logFC, not the rebuilt one", {
  metab_panel_logfc <- load_template_fn("metab_panel_logfc")
  cut <- log2(1.5)
  exact <- c(0.5840, -2, 1)                  # F1 sits just below the FC cutoff
  linear <- signif(2^abs(exact), 3) * sign(exact)
  rebuilt <- sign(linear) * log2(abs(linear))  # what extract_contrast_table() rebuilds
  expect_gte(abs(rebuilt[1]), cut)           # rounding pushes F1 over the cutoff
  ctr_tbl <- data.frame(feature_id = c("F1", "F2", "F9"), logFC = rebuilt)
  de_tables <- list(A_vs_B = data.frame(feature_id = c("F2", "F1"),
                                        logFC = exact[c(2, 1)]))
  got <- metab_panel_logfc(ctr_tbl, de_tables, "A_vs_B")
  expect_equal(got[1:2], exact[1:2])
  expect_lt(abs(got[1]), cut)                # so the panel agrees with the pipeline
  expect_equal(got[3], rebuilt[3])           # absent from the DE table: fallback
  expect_equal(metab_panel_logfc(ctr_tbl, NULL, "A_vs_B"), rebuilt)
})
