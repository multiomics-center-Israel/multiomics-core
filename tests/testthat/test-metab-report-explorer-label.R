# The Metabolite Explorer text must name the transform its values carry.
# It used to say "glog10" whatever preprocessing.transform was. The transform
# is resolved once, by mod_met_log(), and the report reads what it recorded.

METAB_ROOT <- normalizePath(if (dir.exists("R")) "." else "../..",
                            mustWork = FALSE)
METAB_TEMPLATE <- file.path(METAB_ROOT, "R", "pipeline", "metabolomics",
                            "templates", "report_metabolomics.Rmd")

# Evaluate the template's `explorer_value_desc <- ...` expression with the
# given `pre` (and a config that must not be consulted), return the text.
explorer_desc_for <- function(pre) {
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
  # A config that disagrees with what ran: the text must follow `pre`.
  env$metab_cfg <- list(preprocessing = list(transform = "glog10"))
  eval(ex[[1]], envir = env)
  env$explorer_value_desc
}

test_that("mod_met_log records the transform it applied", {
  filtered <- list(mat = matrix(c(1, 3, 7, 15), 2), meta = NULL, row_data = NULL)
  cfg <- function(tf) list(modes = list(metabolomics = list(preprocessing = tf)))
  expect_identical(mod_met_log(filtered, cfg(list(transform = "log10")))$transform, "log10")
  expect_identical(mod_met_log(filtered, cfg(list()))$transform, "log2")
})

test_that("explorer text names the transform mod_met_log recorded", {
  m <- matrix(1, 1, 1)
  expect_match(explorer_desc_for(list(expr_log = m, expr_log_transform = "log10")),
               "^log10-transformed intensities")
  expect_match(explorer_desc_for(list(expr_log = m, expr_log_transform = "none")),
               "^untransformed intensities")
  # Cached pre without the field: no transform is named rather than a guess.
  expect_match(explorer_desc_for(list(expr_log = m)), "^transformed intensities")
  expect_match(explorer_desc_for(list(expr_log = m, expr_log_transform = "log2")),
               "before sample normalization and scaling")
})

test_that("explorer text does not claim the DE matrix when it falls back", {
  desc <- explorer_desc_for(list(expr_log = NULL))
  expect_false(grepl("differential testing", desc, fixed = TRUE))
  expect_match(desc, "feature-scaled when scaling is configured", fixed = TRUE)
})

test_that("the template no longer hard-codes glog10 in the explorer text", {
  txt <- readLines(METAB_TEMPLATE, warn = FALSE)
  expect_false(any(grepl("transformed expression (glog10", txt, fixed = TRUE)))
})
