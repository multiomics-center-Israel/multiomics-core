#!/usr/bin/env Rscript
# utils/benchmark_tabpfn_metab.R
#
# Standalone driver for the metabolomics TabPFN-vs-RF classifier benchmark POC.
# It reuses the preprocessed matrix and feature-selection ranking your pipeline
# already produced, so run the metabolomics pipeline first (tar_make()), then:
#
#   # optional, enables the TabPFN side (RF-only otherwise):
#   bash setup-tabpfn-venv.sh
#   export TABPFN_PYTHON="$(pwd)/envs/tabpfn/bin/python"
#
#   Rscript utils/benchmark_tabpfn_metab.R [path/to/config.yaml]
#
# The config path is optional; it defaults to $MULTIOMICS_CONFIG. Prints a 2-row
# AUC table (RF, TabPFN). Re-running with the same seed gives identical numbers.

suppressWarnings(suppressMessages({
    # Source every function in R/ (late binding makes ordering within R/ moot).
    r_files <- list.files("R", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
    invisible(lapply(sort(r_files), source))
}))

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[[1]] else Sys.getenv("MULTIOMICS_CONFIG", "config.yaml")
if (!file.exists(cfg_path)) {
    stop("Config not found: '", cfg_path, "'. Pass a path or set MULTIOMICS_CONFIG.",
         call. = FALSE)
}

config <- validate_config(load_config(cfg_path))
if (is.null(config$modes$metabolomics)) {
    stop("This config has no metabolomics mode — nothing to benchmark.", call. = FALSE)
}

# Read the preprocessed object and feature-selection result from the {targets}
# store the pipeline built. A clear message beats a cryptic tar_read() error.
read_target <- function(name) {
    tryCatch(
        targets::tar_read_raw(name),
        error = function(e) stop(
            "Could not read target '", name, "' from the _targets store. ",
            "Run the metabolomics pipeline first (targets::tar_make()).",
            call. = FALSE
        )
    )
}

pre <- read_target("metab_pre")
fs  <- tryCatch(read_target("metab_feature_sel_res"), error = function(e) NULL)

bench <- compute_metab_classifier_benchmark(pre, fs, config)
tbl <- format_classifier_benchmark(bench)

if (is.null(tbl)) {
    message("Benchmark could not run (see messages above — e.g. not a 2-class problem).")
} else {
    cat("\nClassifier benchmark — CV AUC (higher is better):\n\n")
    print(tbl, row.names = FALSE)
    if (is.na(tbl$AUC[tbl$model == "TabPFN"])) {
        cat("\nTabPFN was skipped (venv/config). See setup-tabpfn-venv.sh to enable it.\n")
    }
    if (bench$n_samples < 40) {
        cat("\nNote: with", bench$n_samples,
            "samples, CV AUC is high-variance — treat differences cautiously.\n")
    }
}
