# tests/testthat/test-mummichog-seed.R
#
# Reproducibility of the pinned mummichog v2 stage (#144). mummichog 2.7.0
# draws its permutation nulls with Python's global stdlib `random` and never
# seeds it, so 06c seeds it through a launcher prelude and runs the subprocess
# with PYTHONHASHSEED=0.
#
#   * seed validation and argument/env builders: pure R, always run;
#   * the launcher prelude: run against a stand-in module with any python3, so
#     the seeding itself is checked without a mummichog install;
#   * two full mummichog runs: only with the pinned venv ($MUMMICHOG_PYTHON).

# ---- seed validation --------------------------------------------------------

test_that(".mmc_as_seed accepts whole numbers in range and returns an integer", {
  expect_identical(.mmc_as_seed(0), 0L)
  expect_identical(.mmc_as_seed(7), 7L)
  expect_identical(.mmc_as_seed(7L), 7L)
  expect_identical(.mmc_as_seed(.Machine$integer.max), .Machine$integer.max)
})

test_that(".mmc_as_seed refuses values as.integer() would quietly change", {
  bad <- list(1.5, NA, NA_integer_, NA_real_, Inf, -Inf, NaN, -1,
              .Machine$integer.max + 1, "7", TRUE, c(1, 2), numeric(0), NULL)
  for (x in bad) {
    expect_error(.mmc_as_seed(x, key = "params.seed"), "params.seed",
                 info = paste(deparse(x), collapse = ""))
  }
})

test_that("mod_mummichog_pinned rejects a bad params.seed before any work", {
  config <- list(
    params = list(seed = 1.5),
    modes  = list(metabolomics = list(enrichment = list(
      mummichog = list(enabled = TRUE))))
  )
  expect_error(
    mod_mummichog_pinned(pre = NULL, de_res = NULL, config = config,
                         out_dir = withr::local_tempdir()),
    "params.seed"
  )
})

# ---- argument and environment builders -------------------------------------

test_that(".mmc_launcher_args seeds random and runs the module as __main__", {
  args <- .mmc_launcher_args(7)
  expect_length(args, 3L)
  expect_identical(args[[1]], "-c")
  expect_match(args[[2]], "random.seed(int(sys.argv[1]))", fixed = TRUE)
  expect_match(args[[2]], "runpy.run_module('mummichog.main', run_name='__main__'",
               fixed = TRUE)
  expect_identical(args[[3]], "7")
  expect_error(.mmc_launcher_args(1.5), "seed")
})

test_that(".mmc_build_cli_args returns mummichog's flags only", {
  args <- .mmc_build_cli_args(
    infile = "in.tsv", project = "proj", network = "human_mfn",
    mode = "pos_default", instrument_ppm = 10, permutations = 100,
    cutoff = 0.05)
  expect_false(any(args == "mummichog.main"))
  # Python's -c belongs to the launcher, so the only -c here is the cutoff.
  expect_identical(sum(args == "-c"), 1L)
})

test_that("the mummichog subprocess runs with PYTHONHASHSEED=0 and a headless backend", {
  env <- .mmc_subprocess_env()
  expect_identical(unname(env[[1]]), "current")
  expect_identical(unname(env[["PYTHONHASHSEED"]]), "0")
  expect_identical(unname(env[["MPLBACKEND"]]), "Agg")
})

# ---- the launcher prelude, against a stand-in module -----------------------

test_that("the launcher seeds Python's RNG and passes the CLI through unchanged", {
  skip_if_not_installed("processx")
  py <- Sys.which("python3")
  if (!nzchar(py)) skip("python3 not available")

  dir <- withr::local_tempdir()
  # The stand-in prints what mummichog would see: a draw from the global RNG,
  # a string hash (fixed only by PYTHONHASHSEED) and its CLI arguments.
  writeLines(c(
    "import random, sys",
    "print(repr(random.random()))",
    "print(hash('mummichog'))",
    "print(' '.join(sys.argv[1:]))"
  ), file.path(dir, "mmc_probe.py"))

  probe <- function(seed) {
    res <- processx::run(
      py, c(.mmc_launcher_args(seed, module = "mmc_probe"), "-f", "in.tsv", "-p", "5"),
      wd = dir, env = .mmc_subprocess_env(), error_on_status = FALSE)
    expect_identical(res$status, 0L, info = res$stderr)
    strsplit(trimws(res$stdout), "\n", fixed = TRUE)[[1]]
  }

  a1 <- probe(42)
  a2 <- probe(42)
  b  <- probe(43)
  expect_length(a1, 3L)
  # Same seed, separate processes: same draw, same hash, same arguments.
  expect_identical(a1, a2)
  # The seed is what drives the draw.
  expect_false(identical(a1[[1]], b[[1]]))
  # The seed argument is consumed; mummichog's getopt sees only its own flags.
  expect_identical(a1[[3]], "-f in.tsv -p 5")
})

# ---- two full mummichog runs (pinned venv only) ----------------------------

test_that("two seeded mummichog runs in separate processes give the same pathway table", {
  py <- Sys.getenv("MUMMICHOG_PYTHON", "")
  if (!nzchar(py) || !file.exists(py)) {
    skip("MUMMICHOG_PYTHON not set; the seeded-rerun check needs the pinned venv")
  }

  # Synthetic input: the smoke test's metabolite-like M+H masses carry the
  # small p-values, over a seeded background of unrelated features.
  n_bg <- 200L
  bg <- withr::with_seed(11, data.frame(
    feature_id     = paste0("bg_", seq_len(n_bg)),
    mz             = stats::runif(n_bg, 100, 900),
    retention_time = stats::runif(n_bg, 10, 300),
    p_value        = stats::runif(n_bg, 0.2, 1),
    statistic      = stats::rnorm(n_bg),
    stringsAsFactors = FALSE
  ))
  hits <- data.frame(
    feature_id     = paste0("hit_", 1:8),
    mz             = c(176.1263, 181.0703, 205.0972, 167.0932,
                       269.0877, 348.0704, 260.0295, 147.0768),
    retention_time = c(199, 82.5, 120, 210, 145, 260, 175, 55),
    p_value        = c(0.004, 0.005, 0.002, 0.001, 0.003, 0.004, 0.002, 0.006),
    statistic      = c(-2.25, 1.77, 2.40, 3.10, -1.90, 2.20, 1.90, -2.00),
    stringsAsFactors = FALSE
  )
  work  <- withr::local_tempdir()
  input <- write_mummichog_input(rbind(hits, bg), file.path(work, "input.tsv"),
                                 id_col = "feature_id")

  run_once <- function(tag) {
    files <- run_mummichog_v2(
      infile = input, out_dir = file.path(work, tag), project = "seed_test",
      python = py, permutations = 50, cutoff = 0.05, seed = 42
    )
    as.data.frame(read_mummichog_pathways(files))
  }
  t1 <- run_once("run1")
  t2 <- run_once("run2")

  # A meaningful comparison: pathways were reported, with usable p-values.
  expect_gt(nrow(t1), 0L)
  expect_true("p-value" %in% names(t1))
  p <- t1[["p-value"]]
  expect_true(is.numeric(p))
  expect_true(all(is.finite(p) & p >= 0 & p <= 1))

  # mummichog writes rows with tied p-values, and the ids inside an overlap
  # cell, in the iteration order of sets of Python objects, which follows
  # memory addresses and is not something a seed controls. Compare rows keyed
  # by pathway and those id lists as sets; every value, p-values included,
  # must match exactly.
  canonical <- function(t) {
    sort_ids <- function(x, sep) {
      vapply(strsplit(as.character(x), sep, fixed = TRUE),
             function(v) paste(sort(v), collapse = sep), character(1))
    }
    for (col in c("overlap_EmpiricalCompounds (id)", "overlap_features (id)")) {
      if (col %in% names(t)) t[[col]] <- sort_ids(t[[col]], ",")
    }
    if ("overlap_features (name)" %in% names(t)) {
      t[["overlap_features (name)"]] <- sort_ids(t[["overlap_features (name)"]], "$")
    }
    t <- t[order(t$pathway, t[["p-value"]]), , drop = FALSE]
    rownames(t) <- NULL
    t
  }
  expect_identical(canonical(t2), canonical(t1))
})
