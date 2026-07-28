# tests/testthat/helper-e2e.R
#
# Shared scaffolding for the end-to-end per-omic regression tests
# (test-e2e-rnaseq.R, test-e2e-proteomics.R, test-e2e-metabolomics.R).
#
# These tests run a single-omic pipeline for real via targets::tar_make() on the
# repo's shipped synthetic data and assert on the differential-abundance result.
# They are heavy, so every test file guards on RUN_E2E_OMICS == "1" (see
# e2e_should_run()); the fast default suite — including CI — skips them.
#
# Isolation contract (never touch the user's working state):
#   * project.dir is redirected to a throwaway temp dir, so all outputs land
#     there, never in the repo tree.
#   * tar_make() runs against a private store under that temp dir, never the
#     real _targets/ cache.
#   * MULTIOMICS_CONFIG is set only for the duration of the run and restored.

#' Repo root, resolved from either the repo root or tests/testthat.
#'
#' helper.R sources all of R/ using the same heuristic; the pipeline's
#' _targets.R sources R/ with paths relative to the repo root, so tar_make()
#' must run from here.
#'
#' @return Absolute path to the repository root.
e2e_repo_root <- function() {
  normalizePath(if (dir.exists("R")) "." else file.path("..", ".."))
}

#' Whether the heavy end-to-end omic tests should run.
#'
#' @return TRUE when the RUN_E2E_OMICS environment variable is "1".
e2e_should_run <- function() {
  identical(Sys.getenv("RUN_E2E_OMICS"), "1")
}

#' Whether the extra (heavier) file-output assertions should run.
#'
#' Building the export targets pulls in QC/clustering/writers, so it is gated
#' behind its own flag on top of RUN_E2E_OMICS.
#'
#' @return TRUE when RUN_E2E_OMICS_FULL == "1".
e2e_should_run_full <- function() {
  identical(Sys.getenv("RUN_E2E_OMICS_FULL"), "1")
}

#' Build a data-staging function that copies shipped synthetic inputs.
#'
#' @param src_rel Path under the repo's `data/` dir to copy from
#'   (e.g. "example_metabolomics/metabolomics").
#' @param dest_sub Sub-directory name to create under the temp project's raw dir
#'   (e.g. "metabolomics"), matching the config's `files:` paths.
#' @return A function `(raw_dir, repo_root)` that stages the data.
e2e_stage_shipped <- function(src_rel, dest_sub) {
  function(raw_dir, repo_root) {
    src <- file.path(repo_root, "data", src_rel)
    if (!dir.exists(src)) {
      testthat::skip(sprintf("Shipped example data not found: %s", src))
    }
    dest <- file.path(raw_dir, dest_sub)
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    files <- list.files(src, pattern = "\\.csv$", full.names = TRUE)
    file.copy(files, dest, overwrite = TRUE)
    invisible(dest)
  }
}

#' Generate a seeded synthetic RNA-seq count dataset.
#'
#' Writes counts / metadata / contrasts into `dir` matching e2e_rnaseq_config.yaml.
#' Counts are negative-binomial; the first `n_de` genes carry a strong (~4x)
#' up/down signal in Treatment vs Control across a balanced 3-vs-3 design — large
#' enough for a genuine DESeq2 run (the repo's shipped 5-gene toy is not).
#'
#' @param dir Output directory (created if missing).
#' @param seed Integer RNG seed for reproducibility.
#' @param n_genes Number of genes.
#' @param n_de Number of genes carrying injected DE signal (the first `n_de`).
#' @return Invisibly, the character vector of injected (truth) gene IDs.
generate_synthetic_rna <- function(dir, seed = 1L, n_genes = 300L, n_de = 40L) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  samples_ctrl <- c("CTRL_01", "CTRL_02", "CTRL_03")
  samples_trt  <- c("TRT_01", "TRT_02", "TRT_03")
  all_samples  <- c(samples_ctrl, samples_trt)
  gene_ids     <- sprintf("GENE_%04d", seq_len(n_genes))

  withr::with_seed(seed, {
    size      <- 8                                   # NB dispersion (size param)
    base_mean <- exp(stats::rnorm(n_genes, mean = 4.5, sd = 1.2))
    counts <- matrix(
      stats::rnbinom(n_genes * length(all_samples),
                     mu = rep(base_mean, times = length(all_samples)), size = size),
      nrow = n_genes, ncol = length(all_samples),
      dimnames = list(gene_ids, all_samples)
    )
    de_idx   <- seq_len(n_de)
    fold     <- sample(c(0.25, 4), n_de, replace = TRUE)  # strong down / up
    trt_cols <- match(samples_trt, all_samples)
    for (i in seq_along(de_idx)) {
      g <- de_idx[i]
      counts[g, trt_cols] <- stats::rnbinom(length(trt_cols),
                                            mu = base_mean[g] * fold[i], size = size)
    }
  })

  counts_df <- data.frame(gene_id = gene_ids, counts,
                          check.names = FALSE, stringsAsFactors = FALSE)
  utils::write.csv(counts_df, file.path(dir, "counts.csv"), row.names = FALSE)

  meta_df <- data.frame(
    SampleID = all_samples,
    group    = c(rep("Control", length(samples_ctrl)),
                 rep("Treatment", length(samples_trt))),
    stringsAsFactors = FALSE
  )
  utils::write.csv(meta_df, file.path(dir, "metadata.csv"), row.names = FALSE)

  contrasts_df <- data.frame(
    Contrast_name = "Treatment_vs_Control",
    Factor        = "group",
    Numerator     = "Treatment",
    Denominator   = "Control",
    stringsAsFactors = FALSE
  )
  utils::write.csv(contrasts_df, file.path(dir, "contrasts.csv"), row.names = FALSE)

  invisible(gene_ids[seq_len(n_de)])
}

#' Staging function for the synthetic RNA dataset.
#'
#' @param seed RNG seed forwarded to generate_synthetic_rna().
#' @return A function `(raw_dir, repo_root)` suitable for run_omic_e2e().
e2e_stage_rna <- function(seed = 1L) {
  function(raw_dir, repo_root) {
    generate_synthetic_rna(file.path(raw_dir, "rna"), seed = seed)
  }
}

#' Run one omic end-to-end in a fully isolated sandbox.
#'
#' Stages input data into a temp project, redirects `project.dir` there, sets
#' MULTIOMICS_CONFIG, and runs tar_make() against a private store for the
#' requested targets only (their upstream deps are built automatically; nothing
#' downstream is). The temp project is removed when the calling test finishes.
#'
#' @param config_fixture Basename of a YAML fixture under tests/testthat/fixtures/.
#' @param stage_fn Function `(raw_dir, repo_root)` that populates the raw data dir.
#' @param target_names Character vector of target names to build and read back.
#' @return list(values = named list of target values, proj, out_dir, config_path).
run_omic_e2e <- function(config_fixture, stage_fn, target_names) {
  testthat::skip_if_not_installed("targets")
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("withr")
  testthat::skip_if_not_installed("tidyselect")

  repo_root    <- e2e_repo_root()
  fixture_path <- testthat::test_path("fixtures", config_fixture)

  proj <- tempfile("e2e_proj_")
  dir.create(proj)
  # Clean up when the *calling* test_that() block finishes, so file-output
  # assertions can still see the outputs before they are removed.
  withr::defer(unlink(proj, recursive = TRUE, force = TRUE), envir = parent.frame())

  # Symlink R/ and _targets.R into the sandbox and run tar_make() from there.
  # _targets.R sources R/ with paths relative to the working directory, so this
  # keeps the pipeline runnable while ALSO trapping any relative output path a
  # stage writes (e.g. preprocess_rna()'s "outputs/rnaseq/filtering_threshold_qc.png",
  # built from config$paths$out without project$dir) inside proj rather than the
  # repository tree.
  file.symlink(file.path(repo_root, "R"), file.path(proj, "R"))
  file.symlink(file.path(repo_root, "_targets.R"), file.path(proj, "_targets.R"))

  raw_dir <- file.path(proj, "data")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  stage_fn(raw_dir, repo_root)

  cfg <- yaml::read_yaml(fixture_path)
  cfg$project$dir <- proj
  cfg_path <- file.path(proj, "config.yaml")
  yaml::write_yaml(cfg, cfg_path)

  store <- file.path(proj, "_targets_store")

  values <- withr::with_dir(proj, {
    withr::local_envvar(c(MULTIOMICS_CONFIG = cfg_path))
    suppressMessages(
      targets::tar_make(
        names          = tidyselect::all_of(target_names),
        store          = store,
        callr_function = NULL,
        reporter       = "silent"
      )
    )
    stats::setNames(
      lapply(target_names, function(nm) targets::tar_read_raw(nm, store = store)),
      target_names
    )
  })

  out_dir <- file.path(
    proj, cfg$paths$out,
    sprintf("Results_%s_%s", cfg$project$name, cfg$project$analysis_round)
  )

  list(values = values, proj = proj, out_dir = out_dir,
       config_path = cfg_path, config = cfg)
}

#' Extract the feature-id column from a DE summary_df, tolerating naming
#' differences across omics (metabolomics uses `feature_id`; rnaseq/proteomics
#' use `FeatureID`).
#'
#' @param df A summary_df from a DE result.
#' @return Character vector of feature IDs.
e2e_feature_ids <- function(df) {
  candidates <- c("feature_id", "FeatureID", "original_id", "Protein.Group", "gene_id")
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) {
    stop("No recognized feature-id column in summary_df. Columns: ",
         paste(names(df), collapse = ", "), call. = FALSE)
  }
  as.character(df[[hit[1]]])
}

#' Logical vector marking significant rows via the DE contract's
#' `pass_any_contrast` column (numeric 1/NA or logical TRUE/NA).
#'
#' @param df A summary_df from a DE result.
#' @return Logical vector, one per row.
e2e_pass_flag <- function(df) {
  if (!"pass_any_contrast" %in% names(df)) {
    stop("summary_df is missing the 'pass_any_contrast' column.", call. = FALSE)
  }
  x <- df[["pass_any_contrast"]]
  !is.na(x) & (x == 1 | x == TRUE)
}

#' Feature IDs flagged significant in at least one contrast.
#'
#' @param df A summary_df from a DE result.
#' @return Character vector of significant feature IDs.
e2e_significant_ids <- function(df) {
  e2e_feature_ids(df)[e2e_pass_flag(df)]
}
