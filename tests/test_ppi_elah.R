#!/usr/bin/env Rscript
# Standalone test: run PPI module on Elah Pick protein data
# without re-running the full targets pipeline.

project_root <- "/home/multiomics-core"
setwd(project_root)

# --- renv library ---
renv_lib <- file.path(project_root, "renv/library/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu")
if (dir.exists(renv_lib)) .libPaths(c(renv_lib, .libPaths()))

# --- Source all R files: core -> domain -> modules ---
core_files   <- sort(list.files("R/core",               "\\.R$", full.names = TRUE, recursive = TRUE))
domain_files <- sort(list.files("R/domain/proteomics",  "\\.R$", full.names = TRUE, recursive = TRUE))
module_files <- sort(list.files("R/modules/proteomics", "\\.R$", full.names = TRUE, recursive = TRUE))

for (f in c(core_files, domain_files, module_files)) {
    tryCatch(source(f, local = FALSE), error = function(e)
        message("  [skip] ", f, ": ", conditionMessage(e)))
}

# --- Load config ---
config <- load_config("config/prot_elah_pick_protein.yaml")

# --- Load existing DE summary ---
summary_path <- file.path(
    project_root,
    "outputs/elah_pick_protein/Results_Elah_Pick_Protein_1/proteomics/Datasets",
    "limma_multimp_summary_p0.05.tsv"
)
stopifnot(file.exists(summary_path))
summary_df <- read.delim(summary_path, stringsAsFactors = FALSE)
message("Loaded summary_df: ", nrow(summary_df), " rows, ", ncol(summary_df), " cols")

# --- Load metadata ---
meta_path <- file.path(project_root, "data/elah_pick_protein/proteomics/metadata_EP.csv")
stopifnot(file.exists(meta_path))
meta <- read.csv(meta_path, stringsAsFactors = FALSE)
message("Loaded metadata: ", nrow(meta), " rows")

# --- Reconstruct inputs expected by mod_proteomics_ppi ---
de_res <- list(summary_df = summary_df)
pre    <- list(meta = meta)
out_dir <- file.path(
    project_root,
    "outputs/elah_pick_protein/Results_Elah_Pick_Protein_1/proteomics"
)

# --- Run PPI ---
message("\n=== Running PPI module ===")
t0 <- Sys.time()
result <- tryCatch(
    mod_proteomics_ppi(de_res, pre, config, out_dir),
    error = function(e) { message("ERROR: ", conditionMessage(e)); NULL }
)
elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
message("Elapsed: ", elapsed, " seconds")

# --- Report ---
ppi_dir <- file.path(out_dir, "ppi_networks")
if (!is.null(result)) {
    message("\nSUCCESS — PPI analysis completed.")
    files <- list.files(ppi_dir, recursive = TRUE)
    message("Generated ", length(files), " files:")
    for (f in files) message("  ", f)
} else {
    message("\nFAILED — PPI returned NULL. Check messages above.")
    if (dir.exists(ppi_dir)) {
        files <- list.files(ppi_dir, recursive = TRUE)
        if (length(files)) {
            message("Partial outputs in ", ppi_dir, ":")
            for (f in files) message("  ", f)
        }
    }
}
