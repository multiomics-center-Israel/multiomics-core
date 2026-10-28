# tests/testthat/test-e2e-proteomics-batch.R
#
# End-to-end check of proteomics batch handling on the shipped synthetic data
# (data/example_proteomics), with a synthetic Batch column added to a sandbox
# copy of the metadata. Uses the isolation in helper-e2e.R: temp project, private
# targets store, nothing written to the repo or the real _targets/ cache.
#
# Heavy + opt-in: RUN_E2E_OMICS=1 and RUN_E2E_OMICS_FULL=1. Needs sva (ComBat).
# The report is rendered only when pandoc is available.

# 3 Control + 3 Treatment. Each batch holds both conditions, so
# ~ 0 + Condition + Batch has 3 columns and 3 residual degrees of freedom.
PROT_BATCH_OF <- c(CTRL_01 = "B1", CTRL_02 = "B1", TRT_01 = "B1",
                   CTRL_03 = "B2", TRT_02 = "B2", TRT_03 = "B2")

stage_prot_with_batch <- function(raw_dir, repo_root) {
    dest <- e2e_stage_shipped("example_proteomics/proteomics", "proteomics")(raw_dir, repo_root)
    meta_path <- file.path(dest, "metadata.csv")
    meta <- utils::read.csv(meta_path, stringsAsFactors = FALSE)
    stopifnot(setequal(meta$SampleID, names(PROT_BATCH_OF)))
    meta$Batch <- unname(PROT_BATCH_OF[meta$SampleID])
    utils::write.csv(meta, meta_path, row.names = FALSE)
    invisible(dest)
}

# Numeric equality that ignores names, dims and other attributes.
same_values <- function(a, b, tolerance = 0) {
    isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tolerance,
                     check.attributes = FALSE))
}

find_one <- function(root, pattern) {
    hits <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
    if (length(hits) != 1L) stop(sprintf("expected one '%s' under the sandbox, found %d",
                                         pattern, length(hits)), call. = FALSE)
    hits
}

# Pandoc wraps lines, and the template hides inactive sections inside <!-- -->
# (kept in the HTML), so flatten whitespace and drop comments before matching.
flat_html <- function(path) {
    html <- gsub("\\s+", " ", paste(readLines(path, warn = FALSE), collapse = " "))
    gsub("<!--.*?-->", "", html, perl = TRUE)
}

test_that("proteomics batch handling end to end: ComBat for QC, batch covariate in DE", {
    if (!e2e_should_run_full()) skip("Set RUN_E2E_OMICS=1 and RUN_E2E_OMICS_FULL=1 to run.")
    skip_if_not_installed("limma")
    skip_if_not_installed("sva")
    skip_if_not_installed("openxlsx")

    render <- requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()
    tgts <- c("prot_pre_raw", "prot_pre", "prot_de_res", "prot_qc_pre_obj",
              "prot_clustering_obj", "prot_exports", if (render) "prot_report")

    res <- run_omic_e2e("e2e_proteomics_batch_config.yaml", stage_prot_with_batch, tgts)
    v   <- res$values
    raw <- v$prot_pre_raw
    pre <- v$prot_pre
    de  <- v$prot_de_res

    # --- ComBat ran, and only on the QC matrices ---------------------------------
    expect_true(isTRUE(pre$batch_corrected))
    expect_identical(pre$batch_method, "combat")
    expect_true(same_values(pre$expr_filt, raw$expr_filt))
    expect_identical(is.na(pre$expr_filt), is.na(raw$expr_filt))
    expect_false(same_values(pre$expr_imp_single, raw$expr_imp_single, tolerance = 1e-8))
    expect_true(same_values(pre$expr_work, pre$expr_imp_single))

    # --- DE: uncorrected input, batch in every design, every draw ---------------
    expect_identical(get_proteomics_batch_config(res$config$modes$proteomics)$de_mode,
                     "covariate")
    obs <- !is.na(pre$expr_filt)
    expect_gt(length(de$imputations), 1L)
    for (i in seq_along(de$imputations)) {
        imp <- as.matrix(de$imputations[[i]])[rownames(pre$expr_filt), colnames(pre$expr_filt)]
        expect_true(same_values(imp[obs], pre$expr_filt[obs]), info = sprintf("draw %d", i))
        expect_identical(unname(de$runs[[i]]$design), unname(de$runs[[1]]$design))
    }
    expect_true("batch_B2" %in% colnames(de$runs[[1]]$design))
    expect_false(is.null(de$batch_design))
    expect_equal(de$batch_design$df_residual, 3)

    # Each draw's fit equals a direct limma fit built from the metadata alone.
    meta <- pre$meta[match(colnames(de$imputations[[1]]), pre$meta$SampleID), ]
    X <- cbind(Control   = as.numeric(meta$Condition == "Control"),
               Treatment = as.numeric(meta$Condition == "Treatment"),
               batch_B2  = as.numeric(meta$Batch == "B2"))
    cm <- limma::makeContrasts(Treatment_vs_Control = Treatment - Control, levels = X)
    for (i in seq_along(de$imputations)) {
        ref <- limma::eBayes(limma::contrasts.fit(limma::lmFit(de$imputations[[i]], X), cm))
        got <- de$runs[[i]]$fit2
        for (f in c("coefficients", "t", "p.value")) {
            expect_true(same_values(got[[f]][, "Treatment_vs_Control"], ref[[f]][, 1],
                                    tolerance = 1e-10), info = sprintf("draw %d %s", i, f))
        }
    }

    # --- Rubin pooling over the batch-adjusted runs ------------------------------
    pooled <- pool_limma_imputations(de$runs, "Treatment_vs_Control", observed = obs)
    expect_equal(nrow(pooled), nrow(pre$expr_filt))
    expect_true(all(is.finite(pooled$pooled_logFC)))
    expect_true(all(pooled$n_samples == 6L))

    # --- PCA and clustering read the corrected QC matrix -------------------------
    pca <- v$prot_qc_pre_obj$objects$norm_log_counts_pca
    ctr <- pca$center
    expect_true(same_values(ctr, rowMeans(pre$expr_imp_single)[names(ctr)], tolerance = 1e-10))
    expect_false(same_values(ctr, rowMeans(raw$expr_imp_single)[names(ctr)], tolerance = 1e-8))
    # The hierarchical step returns zscore_rows() of the matrix it clustered.
    # Row z-scores ignore a per-protein shift but not a per-batch one, so they
    # tell the corrected matrix from the uncorrected one.
    z <- v$prot_clustering_obj$excel_order$zscore_mat
    expect_false(is.null(z))   # NULL means clustering skipped (< 2 DE features)
    ids <- rownames(z)
    qc_cols <- colnames(pre$expr_imp_single)
    expect_true(same_values(z, zscore_rows(pre$expr_imp_single[ids, qc_cols]), tolerance = 1e-10))
    expect_false(same_values(z, zscore_rows(raw$expr_imp_single[ids, qc_cols]), tolerance = 1e-8))

    # --- Exports and Shiny payload ------------------------------------------------
    files <- list.files(res$proj, recursive = TRUE, full.names = TRUE)
    expect_true(any(basename(files) == "final_results.tsv"))
    expect_true(any(grepl("\\.xlsx$", files)))
    cells_f <- find_one(res$proj, "^batch_design_cells\\.tsv$")
    contr_f <- find_one(res$proj, "^batch_design_contrasts\\.tsv$")
    cells <- utils::read.delim(cells_f, check.names = FALSE)
    expect_setequal(colnames(cells), c("Condition", "B1", "B2"))
    expect_true("Treatment_vs_Control" %in% utils::read.delim(contr_f)$Contrast)

    payload <- readRDS(find_one(res$proj, "^shiny_payload_proteomics_.*\\.rds$"))
    rn <- rownames(pre$expr_filt); cn <- colnames(pre$expr_filt)
    expect_identical(is.na(payload$expr_raw[rn, cn]), is.na(pre$expr_filt))
    expect_true(same_values(payload$expr_raw[rn, cn][obs], pre$expr_filt[obs]))
    expect_true(same_values(payload$expr_norm[rn, cn], pre$expr_imp_single[rn, cn]))

    # --- Report ------------------------------------------------------------------
    if (render) {
        html <- flat_html(find_one(res$proj, "^report_proteomics\\.html$"))
        expect_match(html, "Batch in the DE model", fixed = TRUE)
        expect_match(html, "Samples per condition and batch", fixed = TRUE)
        expect_match(html, "Batch effects were corrected with ComBat for the batch column Batch",
                     fixed = TRUE)
        expect_match(html, "with batch as a covariate", fixed = TRUE)
        expect_match(html, "included Batch as an additive batch covariate and was fitted on uncorrected values",
                     fixed = TRUE)
    }

    # --- Same sandbox, batch adjustment off: no stale diagnostics ---------------
    cfg <- res$config
    cfg$modes$proteomics$batch_correction <- list(batch_col = "Batch", method = "none")
    yaml::write_yaml(cfg, res$config_path)    # tracked as a file target -> rebuild
    store <- file.path(res$proj, "_targets_store")
    off_tgts <- c("prot_de_res", "prot_exports", if (render) "prot_report")
    de_off <- withr::with_dir(res$proj, {
        withr::local_envvar(c(MULTIOMICS_CONFIG = res$config_path))
        suppressMessages(targets::tar_make(names = tidyselect::all_of(off_tgts), store = store,
                                           callr_function = NULL, reporter = "silent"))
        targets::tar_read_raw("prot_de_res", store = store)
    })
    expect_null(de_off$batch_design)
    expect_false("batch_B2" %in% colnames(de_off$runs[[1]]$design))
    left <- list.files(res$proj, pattern = "^batch_design_.*\\.tsv$", recursive = TRUE)
    expect_length(left, 0)
    if (render) {
        html_off <- flat_html(find_one(res$proj, "^report_proteomics\\.html$"))
        expect_false(grepl("Batch in the DE model", html_off, fixed = TRUE))
        expect_false(grepl("Batch effects were corrected", html_off, fixed = TRUE))
    }

    # The shipped example data was only ever read.
    shipped <- utils::read.csv(file.path(e2e_repo_root(), "data", "example_proteomics",
                                         "proteomics", "metadata.csv"))
    expect_false("Batch" %in% colnames(shipped))
})
