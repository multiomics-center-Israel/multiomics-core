# The Metabolite Explorer text must name the transform its values carry.
# It used to say "glog10" whatever preprocessing.transform was.

METAB_ROOT <- normalizePath(if (dir.exists("R")) "." else "../..",
                            mustWork = FALSE)
METAB_TEMPLATE <- file.path(METAB_ROOT, "R", "pipeline", "metabolomics",
                            "templates", "report_metabolomics.Rmd")

# Evaluate the template's `explorer_value_desc <- ...` expression with the
# given `pre` and `metab_cfg`, and return the resulting text.
explorer_desc_for <- function(pre, metab_cfg) {
  lines <- readLines(METAB_TEMPLATE, warn = FALSE)
  start <- grep("^explorer_value_desc <- ", lines)[1]
  if (is.na(start)) stop("explorer_value_desc not found in ", METAB_TEMPLATE)
  for (end in start:length(lines)) {
    ex <- tryCatch(parse(text = lines[start:end]), error = function(e) NULL)
    if (!is.null(ex) && length(ex) == 1L) break
  }
  env <- new.env()
  env$`%||%` <- function(x, y) if (is.null(x)) y else x
  env$pre <- pre
  env$metab_cfg <- metab_cfg
  eval(ex[[1]], envir = env)
  env$explorer_value_desc
}

test_that("explorer text names the configured transform, default log2", {
  pre <- list(expr_log = matrix(1, 1, 1))
  expect_match(explorer_desc_for(pre, list(preprocessing = list(transform = "log10"))),
               "^log10-transformed intensities")
  expect_match(explorer_desc_for(pre, list(preprocessing = list())),
               "^log2-transformed intensities")
  expect_match(explorer_desc_for(pre, list(preprocessing = list(transform = "none"))),
               "^untransformed intensities")
  expect_match(explorer_desc_for(pre, list()), "before sample normalization and scaling")
})

test_that("explorer text falls back when expr_log is absent", {
  desc <- explorer_desc_for(list(expr_log = NULL), list(preprocessing = list(transform = "log2")))
  expect_false(grepl("transformed intensities", desc, fixed = TRUE))
})

test_that("the template no longer hard-codes glog10 in the explorer text", {
  txt <- readLines(METAB_TEMPLATE, warn = FALSE)
  expect_false(any(grepl("transformed expression (glog10", txt, fixed = TRUE)))
})
