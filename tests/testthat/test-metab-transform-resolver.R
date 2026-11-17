# tests/testthat/test-metab-transform-resolver.R
#
# resolve_metab_transform(): the one answer to "what scale is the metabolomics
# working matrix on?". preprocessing.transform reaches only median and none
# (both built on met_log); tss, pqn, eigenms, eigenms_forced and bio_factor
# always log2, whatever the setting says, while info$normalization$transform
# records the configured value. Labels and the group-CV back-transform must
# agree with what ran, so they all ask the resolver.
#
# No data: the resolver only reads two settings.

FIXED_LOG2_NORMS <- c("tss", "pqn", "eigenms", "eigenms_forced", "bio_factor")
TRANSFORMS       <- c("none", "log2", "log10", "glog10")

test_that("linear-scale normalizers resolve to log2 whatever the setting says", {
  for (norm in FIXED_LOG2_NORMS) {
    for (tf in TRANSFORMS) {
      expect_identical(resolve_metab_transform(norm, tf), "log2",
                       info = paste(norm, tf))
    }
    expect_identical(resolve_metab_transform(norm, NULL), "log2", info = norm)
  }
})

test_that("median and none resolve to the configured transform", {
  for (norm in c("median", "none")) {
    for (tf in TRANSFORMS) {
      expect_identical(resolve_metab_transform(norm, tf), tf,
                       info = paste(norm, tf))
    }
  }
})

test_that("an unset transform takes mod_met_log()'s default, log2", {
  expect_identical(resolve_metab_transform("median", NULL), "log2")
  expect_identical(resolve_metab_transform("none", NULL), "log2")
})

test_that("method and transform names are matched case-insensitively", {
  expect_identical(resolve_metab_transform("PQN", "glog10"), "log2")
  expect_identical(resolve_metab_transform("Median", "LOG10"), "log10")
})

test_that("an unknown or missing chosen_norm resolves to NA, not a guess", {
  for (norm in list(NULL, NA, NA_character_, "", "quantile",
                    c("pqn", "median"))) {
    expect_identical(resolve_metab_transform(norm, "log2"), NA_character_,
                     info = paste(deparse(norm), collapse = ""))
  }
})

test_that("group CV back-transforms exactly when the resolver says log2", {
  meta <- data.frame(SampleID = paste0("S", 1:6),
                     group    = rep(c("ctrl", "trt"), each = 3),
                     stringsAsFactors = FALSE)
  lin <- matrix(c(100, 110, 90, 200, 220, 180,
                  50,  55,  45, 60,  66,  54),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("f1", "f2"), meta$SampleID))
  contrasts <- data.frame(Contrast_name = "trt - ctrl", Factor = "group",
                          Numerator = "trt", Denominator = "ctrl",
                          stringsAsFactors = FALSE)
  pre <- list(expr_work = log2(lin + 1), meta = meta)
  cfg_for <- function(norm, tf) list(modes = list(metabolomics = list(
    effects = list(samples = "SampleID"),
    preprocessing = list(chosen_norm = norm, transform = tf,
                         scaling = "none", pseudocount = 1))))

  # eigenms with a configured glog10 still yields a log2 workspace.
  cv <- build_group_cv_metabolomics(pre, contrasts, cfg_for("eigenms", "glog10"))
  expect_false(is.null(cv))
  expect_equal(cv$CV.ctrl[1], 100 * stats::sd(lin[1, 1:3]) / mean(lin[1, 1:3]))

  # median with glog10 does not, so CV is skipped rather than computed wrongly.
  expect_warning(
    cv_skip <- build_group_cv_metabolomics(pre, contrasts, cfg_for("median", "glog10")),
    "provably-log2")
  expect_null(cv_skip)
})

test_that("the report template names the scale only through the resolver", {
  root <- normalizePath(if (dir.exists("R")) "." else "../..", mustWork = FALSE)
  tmpl <- readLines(file.path(root, "R", "pipeline", "metabolomics", "templates",
                              "report_metabolomics.Rmd"), warn = FALSE)
  code <- sub("#.*$", "", tmpl)  # comments may mention the field
  reads <- grep("\\$transform\\b", code, value = TRUE)
  # The single read is the argument passed to resolve_metab_transform() for
  # transform_applied; every label uses that value.
  expect_length(reads, 1L)
  expect_match(reads, "pre$info$normalization$transform", fixed = TRUE)
})
