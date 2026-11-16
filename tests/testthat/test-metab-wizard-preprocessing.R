# tests/testthat/test-metab-wizard-preprocessing.R
#
# The metabolomics block both config wizards write (run.R terminal wizard and
# wizard.html). They used to write a `normalization:` block, which config
# validation rejects, with method names the pipeline does not know ("sum",
# "is") and a transform choice for methods that always use log2. Now:
#   - the settings go under `preprocessing:`, the block the pipeline reads;
#   - the method is a chosen_norm value mod_met_corrected() accepts;
#   - a transform other than log2 is offered only for median and none.
# No data are involved: configs are built in memory and matrices are synthetic.

transforms <- c("log2", "log10", "glog10", "none")

# Mode config as the pipeline would read it from a wizard-written file.
wizard_mode_cfg <- function(...) {
  block <- build_metab_wizard_preprocessing_yaml(...)
  txt <- paste0("modes:\n  metabolomics:\n",
                "    input:\n      format: \"processed_wide\"\n",
                "    files:\n      data: \"metabolomics/data.csv\"\n",
                "    effects:\n      samples: \"sample_id\"\n",
                block, "\n")
  yaml::yaml.load(txt)$modes$metabolomics
}

# Every method / transform pair the wizards may write.
supported_pairs <- function() {
  do.call(rbind, lapply(metab_wizard_norm_methods, function(m) {
    data.frame(method = m,
               transform = if (metab_norm_honours_transform(m)) transforms else "log2",
               stringsAsFactors = FALSE)
  }))
}

test_that("only median and none honour the configurable transform", {
  honours <- vapply(metab_wizard_norm_methods, metab_norm_honours_transform, logical(1))
  expect_setequal(names(honours)[honours], c("median", "none"))
})

test_that("every supported wizard config passes metabolomics config validation", {
  pairs <- supported_pairs()
  for (i in seq_len(nrow(pairs))) {
    for (qc_only in c(FALSE, TRUE)) {
      info <- sprintf("%s / %s / qc_only = %s", pairs$method[i], pairs$transform[i], qc_only)
      cfg <- wizard_mode_cfg(pairs$method[i], transform = pairs$transform[i],
                             scaling = "pareto", qc_only = qc_only,
                             biological_factor_col = "mg_protein")
      expect_null(cfg$normalization, info = info)
      expect_true(validate_metabolomics_config(cfg), info = info)

      pre <- cfg$preprocessing
      if (qc_only) {
        expect_null(pre$chosen_norm, info = info)
      } else {
        expect_identical(pre$chosen_norm, pairs$method[i], info = info)
      }
      expect_identical(pre$transform, pairs$transform[i], info = info)
      expect_identical(pre$scaling, "pareto", info = info)
      expect_equal(pre$pseudocount, 1, info = info)
    }
  }
})

test_that("a transform the method would not apply is refused, not written", {
  fixed_log2 <- metab_wizard_norm_methods[
    !vapply(metab_wizard_norm_methods, metab_norm_honours_transform, logical(1))]
  expect_setequal(fixed_log2, c("pqn", "tss", "eigenms", "eigenms_forced", "bio_factor"))
  for (m in fixed_log2) {
    for (tr in setdiff(transforms, "log2")) {
      expect_error(
        build_metab_wizard_preprocessing_yaml(m, transform = tr,
                                              biological_factor_col = "mg_protein"),
        "always uses log2", info = paste(m, tr))
    }
  }
})

test_that("methods the pipeline does not have are refused", {
  for (m in c("is", "sum", "PQN", NA_character_)) {
    expect_error(build_metab_wizard_preprocessing_yaml(m), "chosen_norm must be one of")
  }
})

test_that("bio_factor needs its column, and other methods do not write one", {
  expect_error(build_metab_wizard_preprocessing_yaml("bio_factor"), "biological_factor_col")
  expect_error(build_metab_wizard_preprocessing_yaml("bio_factor", biological_factor_col = " "),
               "biological_factor_col")

  cfg <- wizard_mode_cfg("bio_factor", biological_factor_col = "protein's \"mg\"")
  expect_identical(cfg$preprocessing$biological_factor_col, "protein's \"mg\"")

  cfg <- wizard_mode_cfg("pqn", biological_factor_col = "mg_protein")
  expect_null(cfg$preprocessing$biological_factor_col)
})

test_that("QC-review mode never writes the bio_factor column", {
  # The bio_factor target is built on every run and normalizes whenever the key
  # is set, so writing it would let a missing column stop a review run.
  cfg <- wizard_mode_cfg("bio_factor", qc_only = TRUE, biological_factor_col = "mg_protein")
  expect_null(cfg$preprocessing$chosen_norm)
  expect_null(cfg$preprocessing$biological_factor_col)
  expect_null(suppressMessages(mod_met_normalize_bio_factor(
    list(mat = matrix(1, 1, 1), meta = NULL, row_data = NULL),
    list(modes = list(metabolomics = cfg)))))

  # Nor is the column required there.
  expect_error(build_metab_wizard_preprocessing_yaml("bio_factor", qc_only = TRUE), NA)
})

test_that("the looser missingness thresholds are written only on request", {
  pre <- wizard_mode_cfg("pqn", relaxed_missingness = TRUE)$preprocessing
  expect_equal(pre$feat_missing_threshold, 0.80)
  expect_equal(pre$sample_missing_threshold, 0.99)

  pre <- wizard_mode_cfg("pqn")$preprocessing
  expect_null(pre$feat_missing_threshold)
  expect_null(pre$sample_missing_threshold)
})

test_that("every wizard method is a chosen_norm that mod_met_corrected() selects", {
  cols <- c("S1", "S2", "S3", "S4")
  meta <- data.frame(sample_id = cols, stringsAsFactors = FALSE)
  base <- matrix(as.numeric(1:8), nrow = 2, dimnames = list(c("f1", "f2"), cols))
  # A different matrix per stage, so the selected branch is identifiable.
  offsets <- c(none = 0, tss = 100, median = 200, pqn = 300, eigenms = 400,
               eigenms_forced = 500, bio_factor = 600)
  expect_setequal(names(offsets), metab_wizard_norm_methods)
  stage <- function(m) list(mat = base + offsets[[m]], meta = meta,
                            row_data = data.frame(feature_id = rownames(base)))

  for (m in metab_wizard_norm_methods) {
    cfg <- list(modes = list(metabolomics = wizard_mode_cfg(
      m, biological_factor_col = "mg_protein")))
    out <- suppressMessages(mod_met_corrected(
      norm_tss = stage("tss"), norm_median = stage("median"), norm_pqn = stage("pqn"),
      logged = stage("none"), meta = meta, out_dir = tempdir(), config = cfg,
      norm_eigenms = stage("eigenms"), norm_eigenms_forced = stage("eigenms_forced"),
      norm_bio_factor = stage("bio_factor")))
    expect_equal(out$mat, base + offsets[[m]], info = m)
    expect_identical(out$info$chosen_norm, m, info = m)
  }
})

# ---- the two wizards themselves ---------------------------------------------
# Neither runs under testthat (one is interactive, one is JavaScript), so these
# read the source: each must build the block through the shared builder and
# offer exactly the methods above.

# Lines of `src` from the line matching `start` up to the next line matching `end`.
source_section <- function(src, start, end) {
  from <- grep(start, src)
  expect_length(from, 1L)
  to <- from + grep(end, src[-seq_len(from)])[1] - 1L
  src[from:to]
}

test_that("the terminal wizard writes the block through the shared builder", {
  src  <- readLines(file.path(root_dir, "run.R"), warn = FALSE)
  body <- source_section(src, "^wizard_metabolomics <- function", "^wizard_lipidomics <- function")

  expect_true(any(grepl("build_metab_wizard_preprocessing_yaml(", body, fixed = TRUE)))
  expect_true(any(grepl("sample_norm <- metab_wizard_norm_methods[norm_idx]", body, fixed = TRUE)))
  expect_true(any(grepl("if (metab_norm_honours_transform(sample_norm))", body, fixed = TRUE)))
  expect_false(any(grepl("^\\s*normalization:", body)))
  expect_false(any(grepl("is_ref_col|Internal standard", body)))

  # The menu and the vector it indexes must have the same length and order.
  menu <- source_section(body, "ask_choice\\(\"Sample normalization method:\"", "default = 1\\)")
  labels <- grep("^\\s*c?\\(?\"", menu[-1], value = TRUE)
  expect_length(labels, length(metab_wizard_norm_methods))
  expect_match(labels[match("tss", metab_wizard_norm_methods)], "Sum normalization")
  expect_match(labels[match("eigenms_forced", metab_wizard_norm_methods)], "EigenMS-forced")
  expect_match(labels[match("none", metab_wizard_norm_methods)], "None")
})

test_that("the browser wizard offers the same methods and the same transform rule", {
  src <- readLines(file.path(root_dir, "wizard.html"), warn = FALSE)

  select <- source_section(src, "<select id=\"metab_sample_norm\"", "</select>")
  values <- sub(".*<option value=\"([^\"]*)\".*", "\\1", grep("<option ", select, value = TRUE))
  expect_setequal(values, metab_wizard_norm_methods)

  list_line <- grep("^var METAB_TRANSFORM_CONFIGURABLE_NORMS = ", src, value = TRUE)
  expect_length(list_line, 1L)
  js_norms <- regmatches(list_line, gregexpr("'[^']*'", list_line))[[1]]
  expect_setequal(gsub("'", "", js_norms),
                  metab_wizard_norm_methods[vapply(metab_wizard_norm_methods,
                                                   metab_norm_honours_transform, logical(1))])

  generator <- source_section(src, "^function generateMetaboYaml\\(", "^function updateMetabNormCondFields\\(")
  expect_true(any(grepl("buildMetabPreprocessingYaml(", generator, fixed = TRUE)))
  expect_true(any(grepl("yaml = '    preprocessing:\\n'", generator, fixed = TRUE)))
  expect_false(any(grepl("normalization:", generator, fixed = TRUE)))
  expect_false(any(grepl("is_ref", generator, fixed = TRUE)))
})
