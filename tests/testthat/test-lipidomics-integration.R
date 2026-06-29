#!/usr/bin/env Rscript
# Test script for lipidomics pipeline
# Runs each stage manually to verify correctness

cat("=== Lipidomics Pipeline Test ===\n\n")

# ---- 1. Source all R files in dependency order ----
cat("1. Sourcing R files...\n")
core_files <- sort(list.files("R/core", pattern = "\\.R$", full.names = TRUE, recursive = TRUE))
service_files <- sort(list.files("R/services", pattern = "\\.R$", full.names = TRUE, recursive = TRUE))
domain_files <- sort(list.files("R/domain", pattern = "\\.R$", full.names = TRUE, recursive = TRUE))
module_files <- sort(list.files("R/modules", pattern = "\\.R$", full.names = TRUE, recursive = TRUE))

suppressWarnings(suppressMessages({
  for (f in c(core_files, service_files, domain_files, module_files)) {
    tryCatch(source(f, local = FALSE), error = function(e) {
      cat("  WARN: could not source", basename(f), ":", e$message, "\n")
    })
  }
}))
cat("   Done sourcing.\n\n")

# ---- 2. Load config ----
cat("2. Loading config...\n")
config <- load_config(testthat::test_path("fixtures", "test_lipidomics_config.yaml"))
# Anchor project$dir to the repo root derived from the test location so data
# paths resolve regardless of the working directory testthat runs under
# (tests/testthat) or where the repo is checked out.
config$project$dir <- normalizePath(testthat::test_path("..", ".."))
cat("   Config loaded. Modes:", paste(names(config$modes), collapse = ", "), "\n\n")

# ---- 3. Validate config ----
cat("3. Validating lipidomics config...\n")
tryCatch({
  validate_lipidomics_config(config$modes$lipidomics)
  cat("   PASS: config validation\n\n")
}, error = function(e) {
  cat("   FAIL: config validation -", e$message, "\n\n")
})

# ---- 4. Load inputs ----
cat("4. Loading lipidomics inputs...\n")
tryCatch({
  inputs <- load_lipidomics_inputs(config)
  cat("   Format:", inputs$format, "\n")
  cat("   Data dimensions:", nrow(inputs$data), "rows x", ncol(inputs$data), "cols\n")
  cat("   Metadata:", if (is.null(inputs$metadata)) "auto-generated from group row" else "provided", "\n")
  if (!is.null(inputs$metadata)) {
    cat("   Metadata dims:", nrow(inputs$metadata), "rows x", ncol(inputs$metadata), "cols\n")
    cat("   Metadata cols:", paste(colnames(inputs$metadata), collapse = ", "), "\n")
    cat("   Conditions:", paste(unique(inputs$metadata$condition), collapse = ", "), "\n")
  }
  cat("   PASS: load inputs\n\n")
}, error = function(e) {
  cat("   FAIL: load inputs -", e$message, "\n\n")
  stop("Cannot continue without inputs")
})

# ---- 5. Preprocess ----
cat("5. Preprocessing...\n")
tryCatch({
  pre <- preprocess_lipidomics(inputs, config)
  cat("   Raw features:", nrow(pre$expr_raw), "\n")
  cat("   Filtered features:", nrow(pre$expr_filt), "\n")
  cat("   Samples:", ncol(pre$expr_work), "\n")
  cat("   Lipid classes:", pre$info$n_lipid_classes, "\n")
  cat("   Top classes:", paste(head(names(sort(table(pre$row_data$lipid_class), decreasing = TRUE)), 5), collapse = ", "), "\n")
  cat("   Metadata cols:", paste(colnames(pre$meta), collapse = ", "), "\n")
  cat("   Sample IDs:", paste(colnames(pre$expr_work), collapse = ", "), "\n")
  cat("   NA count:", pre$info$missingness$total_na, "\n")
  cat("   Normalization:", paste(pre$info$normalization$sample_norm,
                                 pre$info$normalization$transform,
                                 pre$info$normalization$scaling, sep = " / "), "\n")

  # Check chain info
  has_chains <- sum(!is.na(pre$row_data$total_carbons))
  cat("   Features with chain info:", has_chains, "/", nrow(pre$row_data), "\n")

  cat("   PASS: preprocessing\n\n")
}, error = function(e) {
  cat("   FAIL: preprocessing -", e$message, "\n")
  cat("   Traceback:\n")
  traceback()
  cat("\n")
  stop("Cannot continue without pre")
})

# ---- 6. QC module ----
cat("6. Running QC module...\n")
out_dir <- file.path("test_outputs", "lipidomics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tryCatch({
  qc_res <- mod_lipidomics_qc_pre(pre, config, out_dir)
  cat("   Plots generated:", length(qc_res$plots), "\n")
  cat("   Files generated:", length(qc_res$files), "\n")
  cat("   Plot names:", paste(names(qc_res$plots), collapse = ", "), "\n")
  cat("   PASS: QC module\n\n")
}, error = function(e) {
  cat("   FAIL: QC module -", e$message, "\n\n")
  qc_res <<- NULL
})

# ---- 7. DE module ----
cat("7. Running DE module...\n")
tryCatch({
  de_res <- mod_lipidomics_de(pre, config, out_dir)
  n_sig <- sum(de_res$summary_df$pass_any_contrast == 1, na.rm = TRUE)
  cat("   Method:", de_res$method, "\n")
  cat("   Total features:", nrow(de_res$summary_df), "\n")
  cat("   Significant:", n_sig, "\n")
  cat("   Plots:", length(de_res$plots), "\n")
  cat("   Files:", length(de_res$files), "\n")

  # Show top 5 by p-value
  top5 <- head(de_res$de_tables[[1]][order(de_res$de_tables[[1]]$P.Value), ], 5)
  cat("   Top 5 lipids by p-value:\n")
  for (i in seq_len(nrow(top5))) {
    cat("     ", top5$feature_id[i], ": logFC=", round(top5$logFC[i], 3),
        ", p=", signif(top5$P.Value[i], 3), "\n")
  }
  cat("   PASS: DE module\n\n")
}, error = function(e) {
  cat("   FAIL: DE module -", e$message, "\n\n")
  de_res <<- NULL
})

# ---- 8. Feature selection module ----
cat("8. Running feature selection module...\n")
tryCatch({
  fs_res <- mod_lipidomics_feature_selection(pre, config, out_dir)
  if (!is.null(fs_res)) {
    cat("   RF:", if (!is.null(fs_res$rf)) "completed" else "skipped", "\n")
    cat("   PLS-DA:", if (!is.null(fs_res$plsda)) "completed" else "skipped", "\n")
    if (!is.null(fs_res$rf)) {
      cat("   Top 3 RF features:", paste(head(fs_res$rf$importance_df$feature_id, 3), collapse = ", "), "\n")
    }
    if (!is.null(fs_res$plsda)) {
      cat("   Top 3 VIP features:", paste(head(fs_res$plsda$vip_df$feature_id, 3), collapse = ", "), "\n")
      cat("   Explained variance:", paste(round(fs_res$plsda$explained_variance * 100, 1), "%", collapse = ", "), "\n")
    }
    cat("   PASS: feature selection\n\n")
  } else {
    cat("   SKIP: both methods unavailable\n\n")
  }
  fs_res <<- fs_res
}, error = function(e) {
  cat("   FAIL: feature selection -", e$message, "\n\n")
  fs_res <<- NULL
})

# ---- 9. Lipid class analysis module ----
cat("9. Running lipid class analysis module...\n")
tryCatch({
  class_res <- mod_lipidomics_class_analysis(pre, de_res, config, out_dir)
  cat("   Class composition:", if (!is.null(class_res$class_comp)) "computed" else "failed", "\n")
  cat("   Chain saturation:", if (!is.null(class_res$chain_sat)) "computed" else "failed", "\n")
  cat("   Chain length:", if (!is.null(class_res$chain_len)) "computed" else "failed", "\n")
  cat("   Class ORA:", if (!is.null(class_res$class_ora)) paste(nrow(class_res$class_ora), "classes tested") else "skipped/failed", "\n")
  cat("   Plots:", length(class_res$plots), "\n")
  cat("   Files:", length(class_res$files), "\n")

  if (!is.null(class_res$class_comp)) {
    cat("   Top 5 classes (normalized %):\n")
    norm_means <- rowMeans(class_res$class_comp$class_norm)
    top_cl <- head(sort(norm_means, decreasing = TRUE), 5)
    for (nm in names(top_cl)) {
      cat("     ", nm, ":", round(top_cl[nm], 1), "%\n")
    }
  }

  cat("   PASS: lipid class analysis\n\n")
}, error = function(e) {
  cat("   FAIL: lipid class analysis -", e$message, "\n\n")
  class_res <<- NULL
})

# ---- 10. Report rendering ----
cat("10. Rendering HTML report...\n")
tryCatch({
  report_path <- mod_lipidomics_report(pre, qc_res, de_res, fs_res, class_res, config, out_dir)
  cat("   Report:", report_path, "\n")
  cat("   Exists:", file.exists(report_path), "\n")
  cat("   Size:", round(file.info(report_path)$size / 1024, 1), "KB\n")
  cat("   PASS: report rendering\n\n")
}, error = function(e) {
  cat("   FAIL: report rendering -", e$message, "\n\n")
})

# ---- Summary ----
cat("=== Test Complete ===\n")
cat("Output files:\n")
all_files <- list.files("test_outputs", recursive = TRUE, full.names = TRUE)
for (f in all_files) {
  cat("  ", f, "(", round(file.info(f)$size / 1024, 1), "KB)\n")
}
