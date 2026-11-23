# The metabolomics report must not send the pipeline pass flag to the cutoff
# panel. The panel treats passFlag as a hard gate, which froze the volcano and
# DE counts when switching to raw p or loosening cutoffs. Metabolomics
# pass.<contrast> is a plain p / |log2FC| rule the panel recomputes itself.

METAB_ROOT <- normalizePath(if (dir.exists("R")) "." else "../..",
                            mustWork = FALSE)
METAB_TEMPLATE <- file.path(METAB_ROOT, "R", "pipeline", "metabolomics",
                            "templates", "report_metabolomics.Rmd")

test_that("metabolomics plot registrations carry no passFlag column", {
  code <- readLines(METAB_TEMPLATE, warn = FALSE)
  code <- code[!grepl("^\\s*#", code)]
  expect_false(any(grepl("\\$passFlag\\s*<-", code)))
  expect_false(any(grepl("passFlag\\s*=", code)))
})
