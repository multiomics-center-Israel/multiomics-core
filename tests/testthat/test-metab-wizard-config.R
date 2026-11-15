# tests/testthat/test-metab-wizard-config.R
#
# Both config wizards (the terminal wizard in run.R and the browser wizard in
# wizard.html) must write a metabolomics config the pipeline accepts: settings
# under `preprocessing:`, `chosen_norm` in the pipeline's method names, and no
# `normalization:` block, which validate_metabolomics_config() rejects.
#
# The YAML checked here is produced by the wizards' own generating code:
# metab_preprocessing_yaml() is evaluated out of run.R, and
# metabPreprocessingYaml() out of wizard.html (the latter needs the V8 package
# and is skipped without it).
#
# All data are synthetic.

wizard_repo_root <- normalizePath(if (dir.exists("R")) "." else "../..", mustWork = FALSE)

# Every method the wizards offer, by the name they write to chosen_norm.
wizard_norms     <- c("pqn", "median", "tss", "eigenms", "eigenms_forced", "bio_factor", "none")
# Methods whose pipeline branch always applies log2.
log2_only_norms  <- c("tss", "pqn", "eigenms", "eigenms_forced", "bio_factor")
wizard_transforms <- c("log2", "log10", "glog10", "none")

# Evaluate named top-level definitions out of run.R. run.R cannot be sourced:
# it ends by calling main().
run_r_definitions <- function(names) {
  env <- new.env(parent = globalenv())
  exprs <- parse(file.path(wizard_repo_root, "run.R"), keep.source = FALSE, encoding = "UTF-8")
  for (e in exprs) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% names) {
      eval(e, env)
    }
  }
  env
}

# A whole config around a wizard-written preprocessing block, parsed the way a
# config file reaches the pipeline.
config_with_preprocessing <- function(preprocessing_yaml) {
  yaml::yaml.load(paste(
    "paths:",
    "  raw: \"data/x\"",
    "  out: \"outputs/x\"",
    "params:",
    "  seed: 1",
    "modes:",
    "  metabolomics:",
    preprocessing_yaml,
    "    input:",
    "      format: \"processed_wide\"",
    "    files:",
    "      data: \"x.csv\"",
    "      metadata: \"m.csv\"",
    "    effects:",
    "      samples: \"sample_id\"",
    sep = "\n"
  ))
}

# One stage per normalization target, each with a different constant matrix, so
# the matrix mod_met_corrected() returns identifies the branch it took.
wizard_stages <- function() {
  ids <- c("S1", "S2", "S3", "S4")
  stage <- function(value) {
    list(mat = matrix(value, nrow = 3, ncol = 4,
                      dimnames = list(c("f1", "f2", "f3"), ids)),
         meta = data.frame(sample_id = ids, stringsAsFactors = FALSE),
         row_data = data.frame(feature_id = c("f1", "f2", "f3")))
  }
  stats::setNames(lapply(seq_along(wizard_norms), stage), wizard_norms)
}

expect_selection_applied <- function(cfg, method) {
  st  <- wizard_stages()
  out <- suppressMessages(mod_met_corrected(
    norm_tss = st$tss, norm_median = st$median, norm_pqn = st$pqn,
    logged = st$none, meta = st$none$meta, out_dir = tempdir(), config = cfg,
    norm_eigenms = st$eigenms, norm_eigenms_forced = st$eigenms_forced,
    norm_bio_factor = st$bio_factor
  ))
  expect_identical(out$info$chosen_norm, method, info = method)
  expect_identical(out$mat, st[[method]]$mat, info = method)
}

test_that("run.R wizard configs validate and select the chosen normalization", {
  run_r <- run_r_definitions("metab_preprocessing_yaml")
  expect_true(is.function(run_r$metab_preprocessing_yaml))

  for (method in wizard_norms) {
    transforms <- if (method %in% log2_only_norms) "log2" else wizard_transforms
    for (transform in transforms) {
      label <- paste(method, transform)
      bio   <- if (method == "bio_factor") "mg_protein" else NULL
      cfg <- config_with_preprocessing(run_r$metab_preprocessing_yaml(
        method, transform, "none", bio_factor_col = bio))
      pre <- cfg$modes$metabolomics$preprocessing

      expect_null(cfg$modes$metabolomics$normalization, info = label)
      expect_error(validate_config(cfg), NA, info = label)
      expect_identical(pre$chosen_norm, method, info = label)
      expect_identical(pre$transform, transform, info = label)
      expect_identical(pre$scaling, "none", info = label)
      expect_identical(pre$biological_factor_col, bio, info = label)
    }
    expect_selection_applied(cfg, method)
  }
})

test_that("run.R wizard keeps scaling, QC-only mode and multi-level thresholds", {
  run_r <- run_r_definitions("metab_preprocessing_yaml")

  for (scaling in c("none", "auto", "pareto", "range")) {
    cfg <- config_with_preprocessing(run_r$metab_preprocessing_yaml("median", "log2", scaling))
    expect_error(validate_config(cfg), NA, info = scaling)
    expect_identical(cfg$modes$metabolomics$preprocessing$scaling, scaling)
  }

  # QC-review mode is chosen_norm: null.
  qc <- config_with_preprocessing(run_r$metab_preprocessing_yaml(
    "pqn", "log2", "none", qc_only = TRUE, is_multi_level = TRUE))
  pre <- qc$modes$metabolomics$preprocessing
  expect_error(validate_config(qc), NA)
  expect_null(pre$chosen_norm)
  expect_identical(pre$feat_missing_threshold, 0.80)
  expect_identical(pre$sample_missing_threshold, 0.99)
})

test_that("neither wizard still writes the rejected legacy keys", {
  run_r <- run_r_definitions("wizard_metabolomics")
  body  <- paste(deparse(run_r$wizard_metabolomics), collapse = "\n")
  html  <- readLines(file.path(wizard_repo_root, "wizard.html"), warn = FALSE, encoding = "UTF-8")
  from  <- grep("^function generateMetaboYaml", html)
  to    <- grep("^function generateLipidYaml", html)
  expect_length(from, 1L)
  expect_length(to, 1L)
  js <- paste(html[from:to], collapse = "\n")

  for (legacy in c("input_already_normalized", "is_ref_col", "na_policy")) {
    expect_false(grepl(legacy, body, fixed = TRUE), info = paste("run.R", legacy))
    expect_false(grepl(legacy, js, fixed = TRUE), info = paste("wizard.html", legacy))
  }
})

test_that("both wizards offer the same methods and the same log2-only set", {
  run_r <- run_r_definitions("metab_log2_only_norms")
  html  <- paste(readLines(file.path(wizard_repo_root, "wizard.html"), warn = FALSE,
                           encoding = "UTF-8"), collapse = "\n")

  select  <- regmatches(html, regexpr("(?s)<select id=\"metab_sample_norm\".*?</select>", html, perl = TRUE))
  expect_length(select, 1L)
  offered <- regmatches(select, gregexpr("(?<=value=\")[^\"]+", select, perl = TRUE))[[1]]
  expect_setequal(offered, wizard_norms)

  js_line <- regmatches(html, regexpr("var METAB_LOG2_ONLY_NORMS = \\[[^]]*\\]", html))
  expect_length(js_line, 1L)
  js_set  <- regmatches(js_line, gregexpr("(?<=')[a-z0-9_]+(?=')", js_line, perl = TRUE))[[1]]
  expect_setequal(js_set, log2_only_norms)
  expect_setequal(run_r$metab_log2_only_norms, log2_only_norms)
})

test_that("wizard.html configs validate and select the chosen normalization", {
  skip_if_not_installed("V8")
  html <- readLines(file.path(wizard_repo_root, "wizard.html"), warn = FALSE, encoding = "UTF-8")
  from <- grep("^function metabPreprocessingYaml", html)
  expect_length(from, 1L)
  to   <- from + which(html[-seq_len(from)] == "}")[1]
  ctx  <- V8::v8()
  ctx$eval(paste(c(grep("^function yamlStr", html, value = TRUE),
                   grep("^var METAB_LOG2_ONLY_NORMS", html, value = TRUE),
                   html[from:to]), collapse = "\n"))

  for (method in wizard_norms) {
    for (transform in wizard_transforms) {
      label <- paste(method, transform)
      cfg <- config_with_preprocessing(ctx$call(
        "metabPreprocessingYaml", method, transform, "none", "mg_protein", FALSE))
      pre <- cfg$modes$metabolomics$preprocessing

      expect_null(cfg$modes$metabolomics$normalization, info = label)
      expect_error(validate_config(cfg), NA, info = label)
      expect_identical(pre$chosen_norm, method, info = label)
      # The form locks these methods to log2; the generator must agree with it.
      expect_identical(pre$transform,
                       if (method %in% log2_only_norms) "log2" else transform,
                       info = label)
      expect_identical(pre$biological_factor_col,
                       if (method == "bio_factor") "mg_protein" else NULL,
                       info = label)
    }
    expect_selection_applied(cfg, method)
  }
})
