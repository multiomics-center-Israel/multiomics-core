# tests/testthat/test-metab-qc-samples.R
#
# Which metabolomics samples are QC is decided by qc.qc_flag_column and
# qc.qc_flag_value alone (metab_is_qc_sample()). These tests pin:
#   - the config contract: values, two explicit nulls, or an error;
#   - that a missing key and an explicit null stay distinguishable after a
#     config is read from a YAML file;
#   - exact matching, with no detection from sample or group names;
#   - that DE, QEA/ssGSEA, clustering and the "without QC" plots receive the
#     same biological samples;
#   - that the functions shared with lipidomics behave as before when they are
#     not given an explicit mask.
#
# All data is synthetic.

# ---- fixtures ----------------------------------------------------------------

QC_NONE <- list(qc_flag_column = NULL, qc_flag_value = NULL)
QC_TREATMENT <- list(qc_flag_column = "treatment", qc_flag_value = "QC")

# Three biological groups, one of them named like a pool, one biological sample
# whose id starts with "QC", and two pooled-QC samples whose value differs from
# the configured one in case and surrounding spaces only.
.qc_fixture <- function(shuffle_meta = FALSE) {
  ids       <- c("S1", "S2", "S3", "S4", "S5", "S6", "QC_run1", "Q1", "Q2")
  treatment <- c("A", "A", "B", "B", "pool_A", "pool_A", "A", "QC", " qc ")
  mat <- withr::with_seed(11, matrix(
    stats::rnorm(6 * length(ids), mean = 10, sd = 0.5), nrow = 6,
    dimnames = list(paste0("f", 1:6), ids)
  ))
  meta <- data.frame(sample_id = ids, treatment = treatment,
                     stringsAsFactors = FALSE)
  if (shuffle_meta) {
    meta <- meta[withr::with_seed(5, sample(nrow(meta))), , drop = FALSE]
  }
  pre <- list(expr_raw = mat, expr_filt = mat, expr_work = mat, expr_log = mat,
              meta = meta,
              row_data = data.frame(feature_id = rownames(mat),
                                    stringsAsFactors = FALSE),
              info = list(normalization = list(scaling = "none")))
  list(mat = mat, meta = meta, pre = pre,
       bio_ids = c("S1", "S2", "S3", "S4", "S5", "S6", "QC_run1"),
       qc_ids  = c("Q1", "Q2"))
}

.qc_config <- function(qc) {
  list(modes = list(metabolomics = list(
    effects    = list(samples = "sample_id", color = "treatment"),
    de         = list(method = "t_test", condition_column = "treatment"),
    clustering = list(enabled = FALSE),
    qc         = qc
  )))
}

# A config that passes every other validation rule, so the QC keys are the only
# thing a test varies.
.valid_mode_cfg <- function(qc) {
  cfg <- list(input = list(format = "processed_wide"),
              files = list(data = "data.csv"))
  cfg["qc"] <- list(qc)
  cfg
}

# Stub by assignment into the functions' own environment (these are sourced
# functions, not a package, so with_mocked_bindings() has nothing to rebind).
.local_qc_stubs <- function(stubs, env = parent.frame()) {
  target <- environment(filter_to_biological)
  old <- lapply(names(stubs), get, envir = target, inherits = FALSE)
  names(old) <- names(stubs)
  withr::defer(
    for (nm in names(old)) assign(nm, old[[nm]], envir = target),
    envir = env
  )
  for (nm in names(stubs)) assign(nm, stubs[[nm]], envir = target)
  invisible(NULL)
}


# ---- the config contract -----------------------------------------------------

test_that("a column and a value name the QC samples by exact, case-insensitive match", {
  fx <- .qc_fixture()
  is_qc <- metab_is_qc_sample(QC_TREATMENT, fx$meta)

  expect_type(is_qc, "logical")
  expect_length(is_qc, nrow(fx$meta))
  # "QC" and " qc " match; "pool_A" and the sample called QC_run1 do not.
  expect_identical(fx$meta$sample_id[is_qc], fx$qc_ids)
})

test_that("the value is matched whole, never as a fragment", {
  meta <- data.frame(sample_id = paste0("S", 1:5),
                     treatment = c("QC", "QC2", "pooled_QC", "not qc", NA),
                     stringsAsFactors = FALSE)
  expect_identical(metab_is_qc_sample(QC_TREATMENT, meta),
                   c(TRUE, FALSE, FALSE, FALSE, FALSE))
})

test_that("a logical flag column works with a logical value", {
  meta <- data.frame(sample_id = paste0("S", 1:3),
                     is_pooled = c(FALSE, TRUE, NA))
  expect_identical(
    metab_is_qc_sample(list(qc_flag_column = "is_pooled", qc_flag_value = TRUE),
                       meta),
    c(FALSE, TRUE, FALSE))
})

test_that("two explicit nulls mean no QC samples, without a warning", {
  fx <- .qc_fixture()
  expect_silent(is_qc <- metab_is_qc_sample(QC_NONE, fx$meta))
  expect_identical(is_qc, rep(FALSE, nrow(fx$meta)))
})

test_that("missing QC keys are a config error that shows both ways to write them", {
  fx <- .qc_fixture()
  for (qc in list(NULL, list(), list(enabled = TRUE))) {
    expect_error(metab_is_qc_sample(qc, fx$meta), "are not set")
    expect_error(metab_is_qc_sample(qc), "are not set")
  }
  # Only one of the two keys written.
  expect_error(metab_is_qc_sample(list(qc_flag_column = "treatment"), fx$meta),
               "qc\\.qc_flag_value is not set")
  expect_error(metab_is_qc_sample(list(qc_flag_value = "QC"), fx$meta),
               "qc\\.qc_flag_column is not set")

  # The message carries a paste-ready example of each valid state.
  msg <- tryCatch(metab_is_qc_sample(list()), error = conditionMessage)
  expect_match(msg, "qc_flag_column: \"treatment\"", fixed = TRUE)
  expect_match(msg, "qc_flag_value: \"QC\"", fixed = TRUE)
  expect_match(msg, "qc_flag_column: null", fixed = TRUE)
  expect_match(msg, "qc_flag_value: null", fixed = TRUE)
})

test_that("one key with a value and the other null is a config error", {
  expect_error(
    metab_is_qc_sample(list(qc_flag_column = "treatment", qc_flag_value = NULL)),
    "set together or both be null")
  expect_error(
    metab_is_qc_sample(list(qc_flag_column = NULL, qc_flag_value = "QC")),
    "set together or both be null")
  expect_error(
    metab_is_qc_sample(list(qc_flag_column = "", qc_flag_value = "QC")),
    "name of one metadata column")
  expect_error(
    metab_is_qc_sample(list(qc_flag_column = "treatment",
                            qc_flag_value = c("QC", "pool"))),
    "single value")
})

test_that("a configured column that the metadata lacks stops the run", {
  fx <- .qc_fixture()
  expect_error(
    metab_is_qc_sample(list(qc_flag_column = "sample_kind", qc_flag_value = "QC"),
                       fx$meta),
    "no such column")
})

test_that("a value no sample carries warns and marks nothing", {
  fx <- .qc_fixture()
  qc <- list(qc_flag_column = "treatment", qc_flag_value = "reference")
  expect_warning(is_qc <- metab_is_qc_sample(qc, fx$meta),
                 "no sample has the value 'reference'")
  expect_false(any(is_qc))
  # Callers working on one contrast's subset turn the warning off.
  expect_silent(metab_is_qc_sample(qc, fx$meta, warn_unmatched = FALSE))
})

test_that("validation fails early on missing QC keys and accepts both valid states", {
  expect_error(validate_metabolomics_config(.valid_mode_cfg(NULL)),
               "are not set")
  expect_error(validate_metabolomics_config(.valid_mode_cfg(list(enabled = TRUE))),
               "are not set")
  expect_error(
    validate_metabolomics_config(.valid_mode_cfg(list(qc_flag_column = "treatment"))),
    "qc\\.qc_flag_value is not set")

  expect_silent(validate_metabolomics_config(.valid_mode_cfg(QC_NONE)))
  expect_silent(validate_metabolomics_config(.valid_mode_cfg(QC_TREATMENT)))
})


# ---- missing key vs explicit null, through a YAML file -----------------------

test_that("a missing key and an explicit null stay distinct after loading YAML", {
  skip_if_not_installed("yaml")

  write_cfg <- function(qc_lines) {
    path <- withr::local_tempfile(fileext = ".yaml", .local_envir = parent.frame())
    writeLines(c(
      "modes:",
      "  metabolomics:",
      "    input:",
      "      format: \"processed_wide\"",
      "    files:",
      "      data: \"data.csv\"",
      "    qc:",
      "      enabled: true",
      qc_lines
    ), path)
    load_config(path)$modes$metabolomics
  }
  fx <- .qc_fixture()

  # Explicit nulls, in both YAML spellings: the keys survive loading.
  for (null in c("null", "~")) {
    cfg <- write_cfg(c(paste("      qc_flag_column:", null),
                       paste("      qc_flag_value:", null)))
    expect_true(all(c("qc_flag_column", "qc_flag_value") %in% names(cfg$qc)),
                info = null)
    expect_silent(validate_metabolomics_config(cfg))
    expect_false(any(metab_is_qc_sample(cfg$qc, fx$meta)), info = null)
  }

  # Keys absent: the same file without the two lines is rejected.
  cfg <- write_cfg(character(0))
  expect_false(any(c("qc_flag_column", "qc_flag_value") %in% names(cfg$qc)))
  expect_error(validate_metabolomics_config(cfg), "are not set")

  # One key absent.
  cfg <- write_cfg("      qc_flag_column: \"treatment\"")
  expect_error(validate_metabolomics_config(cfg), "qc\\.qc_flag_value is not set")

  # Values.
  cfg <- write_cfg(c("      qc_flag_column: \"treatment\"",
                     "      qc_flag_value: \"QC\""))
  expect_silent(validate_metabolomics_config(cfg))
  expect_identical(fx$meta$sample_id[metab_is_qc_sample(cfg$qc, fx$meta)],
                   fx$qc_ids)
})

test_that("the shipped metabolomics template states its QC samples", {
  skip_if_not_installed("yaml")
  root <- normalizePath(if (dir.exists("R")) "." else "../..")
  qc <- yaml::read_yaml(
    file.path(root, "config/templates/metabolomics_template.yaml")
  )$modes$metabolomics$qc
  expect_true(all(c("qc_flag_column", "qc_flag_value") %in% names(qc)))
  expect_silent(metab_is_qc_sample(qc))
})


# ---- filter_to_biological(): explicit mask vs name-based detection -----------

test_that("an explicit mask removes exactly the masked samples, by sample id", {
  fx <- .qc_fixture(shuffle_meta = TRUE)
  is_qc <- metab_is_qc_sample(QC_TREATMENT, fx$meta)

  res <- suppressMessages(
    filter_to_biological(fx$mat, fx$meta, "treatment", "sample_id",
                         label = "test", is_qc = is_qc))

  expect_setequal(colnames(res$mat), fx$bio_ids)
  # Matrix columns follow the metadata that produced the mask, value for value.
  expect_identical(colnames(res$mat), res$meta$sample_id)
  expect_identical(res$mat, fx$mat[, res$meta$sample_id, drop = FALSE])
  # The biological group called pool_A is a level of the condition.
  expect_setequal(levels(res$condition), c("A", "B", "pool_A"))
})

test_that("an all-FALSE mask keeps every sample; only a missing mask uses names", {
  fx <- .qc_fixture()

  # Declared QC-free: nothing goes, although ids and groups look like QC/pool.
  none <- filter_to_biological(fx$mat, fx$meta, "treatment", "sample_id",
                               label = "test",
                               is_qc = rep(FALSE, nrow(fx$meta)))
  expect_identical(colnames(none$mat), colnames(fx$mat))
  expect_identical(none$meta, fx$meta)

  # No mask at all is the lipidomics path, and it still detects by name.
  legacy <- suppressMessages(
    filter_to_biological(fx$mat, fx$meta, "treatment", "sample_id",
                         label = "test"))
  expect_setequal(colnames(legacy$mat), c("S1", "S2", "S3", "S4"))
})

test_that("a mask that does not fit the metadata is refused", {
  fx <- .qc_fixture()
  for (bad in list(c(TRUE, FALSE), rep(NA, nrow(fx$meta)),
                   rep("FALSE", nrow(fx$meta)))) {
    expect_error(
      filter_to_biological(fx$mat, fx$meta, "treatment", "sample_id",
                           label = "test", is_qc = bad),
      "is_qc must be a logical vector")
  }
})


# ---- build_qc_subsets(): explicit mask vs the treatment-column rule ----------

test_that("build_qc_subsets uses the explicit mask when given one", {
  fx <- .qc_fixture(shuffle_meta = TRUE)
  # Point the config at a column other than `treatment`, so the old hard-coded
  # rule and the configured one disagree.
  fx$meta$run_type <- ifelse(fx$meta$sample_id %in% c("S1", "Q1"),
                             "reference", "study")
  is_qc <- metab_is_qc_sample(
    list(qc_flag_column = "run_type", qc_flag_value = "reference"), fx$meta)

  subsets <- build_qc_subsets(fx$mat, fx$mat, fx$meta, "sample_id",
                              is_qc = is_qc)
  expect_length(subsets, 2L)
  expect_identical(subsets[[2]]$tag, "_noQC")
  expect_setequal(colnames(subsets[[2]]$expr_work),
                  setdiff(colnames(fx$mat), c("S1", "Q1")))
  expect_identical(colnames(subsets[[2]]$expr_work),
                   subsets[[2]]$meta$sample_id)
})

test_that("an experiment without QC gets no 'without QC' subset", {
  fx <- .qc_fixture()
  # treatment still contains "QC", which the old rule would have picked up.
  subsets <- build_qc_subsets(fx$mat, fx$mat, fx$meta, "sample_id",
                              is_qc = metab_is_qc_sample(QC_NONE, fx$meta))
  expect_length(subsets, 1L)
  expect_identical(subsets[[1]]$tag, "")
})

test_that("build_qc_subsets without a mask keeps the treatment-column rule", {
  fx <- .qc_fixture()
  subsets <- build_qc_subsets(fx$mat, fx$mat, fx$meta, "sample_id")
  expect_length(subsets, 2L)
  expect_setequal(colnames(subsets[[2]]$expr_work), fx$bio_ids)

  # No treatment column and no mask: a single subset, as before.
  meta <- fx$meta[, "sample_id", drop = FALSE]
  expect_length(build_qc_subsets(fx$mat, fx$mat, meta, "sample_id"), 1L)
})


# ---- every consumer receives the same biological samples ---------------------

# Runs each metabolomics consumer on the same fixture and returns the sample ids
# each one is left with after QC exclusion.
.consumer_samples <- function(fx, qc) {
  config <- .qc_config(qc)
  seen <- list()

  real_filter  <- filter_to_biological
  real_metrics <- .qc_write_metrics
  .local_qc_stubs(list(
    filter_to_biological = function(...) {
      res <- real_filter(...)
      seen[[list(...)$label]] <<- colnames(res$mat)
      res
    },
    # Record what the suite would plot and report, without drawing anything.
    .qc_write_metrics = function(out_dir, stage, subset_mode, mat, meta, ...,
                                 status) {
      seen[["suite"]]        <<- colnames(mat)
      seen[["suite_status"]] <<- status
      file.path(out_dir, "metrics_summary.tsv")
    }
  ))

  contrasts <- data.frame(Numerator = "B", Denominator = "A",
                          Contrast_name = "B_vs_A", stringsAsFactors = FALSE)
  # Only the sample set DE starts from is under test, not the statistics, so a
  # failure after the filter does not matter here.
  try(suppressWarnings(suppressMessages(
    run_metabolomics_de(fx$pre, config, contrasts))), silent = TRUE)

  seen[["enrichment"]] <- colnames(suppressMessages(
    .metab_de_matrix_condition(fx$pre, config))$mat)

  de_res <- list(summary_df = data.frame(feature_id = character(0),
                                         pass_any_contrast = integer(0)))
  suppressMessages(
    mod_metabolomics_clustering(fx$pre, de_res, config, out_dir = tempdir()))

  subsets <- build_qc_subsets(
    fx$pre$expr_work, fx$pre$expr_filt, fx$pre$meta, "sample_id",
    is_qc = metab_is_qc_sample(qc, fx$pre$meta))
  seen[["qc_pre"]] <- colnames(subsets[[length(subsets)]]$expr_work)

  plots_off <- list(pca = FALSE, density = FALSE, boxplot = FALSE,
                    histogram = FALSE, rle = FALSE, heatmap_correlation = FALSE,
                    heatmap_distance = FALSE, heatmap_expression = FALSE)
  suite_files <- suppressMessages(qc_full_metabolomics_suite(
    mat = fx$pre$expr_work, meta = fx$pre$meta, stage = "log", pseudocount = 1,
    cfg = c(qc, list(plots = plots_off)), cfg_mode = config$modes$metabolomics,
    out_dir = tempdir(), subset_mode = "no_qc"))
  seen[["suite_files"]] <- suite_files
  seen
}

test_that("DE, enrichment, clustering and the no-QC plots share one sample set", {
  fx   <- .qc_fixture(shuffle_meta = TRUE)
  seen <- .consumer_samples(fx, QC_TREATMENT)

  for (consumer in c("metabolomics DE", "enrichment",
                     "metabolomics clustering", "qc_pre", "suite")) {
    expect_setequal(seen[[consumer]], fx$bio_ids)
  }
  # DE and enrichment both align the metadata to the matrix first, so they
  # agree on column order as well as membership.
  expect_identical(seen[["metabolomics DE"]], seen[["enrichment"]])
})

test_that("an experiment without QC runs every consumer on all samples", {
  fx   <- .qc_fixture(shuffle_meta = TRUE)
  all_ids <- colnames(fx$mat)
  expect_no_error(seen <- .consumer_samples(fx, QC_NONE))

  for (consumer in c("metabolomics DE", "enrichment",
                     "metabolomics clustering", "qc_pre")) {
    expect_setequal(seen[[consumer]], all_ids)
  }
  # The no-QC suite has nothing to remove: it reports that and draws nothing.
  expect_identical(seen[["suite_status"]], "skipped_no_qc_samples")
  expect_false(any(grepl("\\.png$", seen[["suite_files"]])))
})

test_that("the no-QC suite stops on a configured column the metadata lacks", {
  fx <- .qc_fixture()
  config <- .qc_config(list(qc_flag_column = "sample_kind",
                            qc_flag_value = "QC"))
  expect_error(
    qc_full_metabolomics_suite(
      mat = fx$mat, meta = fx$meta, stage = "log", pseudocount = 1,
      cfg = config$modes$metabolomics$qc,
      cfg_mode = config$modes$metabolomics,
      out_dir = tempdir(), subset_mode = "no_qc"),
    "no such column")
})


# ---- lipidomics is not touched -----------------------------------------------

test_that("lipidomics callers pass no explicit mask", {
  # Lipidomics shares filter_to_biological() and build_qc_subsets(). It keeps
  # the name-based behaviour only as long as it does not pass is_qc.
  for (fn in list(mod_lipidomics_clustering, mod_lipidomics_qc_pre)) {
    body_txt <- paste(deparse(body(fn)), collapse = " ")
    expect_false(grepl("is_qc", body_txt, fixed = TRUE))
    expect_false(grepl("metab_is_qc_sample", body_txt, fixed = TRUE))
  }
  expect_null(formals(filter_to_biological)$is_qc)
  expect_null(formals(build_qc_subsets)$is_qc)
})
