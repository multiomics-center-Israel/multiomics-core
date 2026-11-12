# Tests for the before/after QC of the RNA batch correction
# (R/modules/rnaseq/01b_mod_batch_qc.R, rna_batch_correction_sentence(),
# draw_rna_sample_structure_plots(), draw_rna_top_de_heatmap())
# test-rna-batch-qc.R
#
# All data here is synthetic; the fixtures are in helper-rna-batch.R.

# The established QC file names for this fixture (200 genes, so no top-N PCA;
# the 3D PCA is switched off in bc_config()).
bc_structure_files <- function(suffix = "") {
    paste0(c("PCA_PC1.vs.PC2", "PCA_PC1.vs.PC3", "sample_distance_heatmap",
             "sample_correlation_heatmap", "samples_rna_heatmap",
             "samples_rna_heatmap_wo_col", "samples_rna_heatmap_unclustered"),
           suffix, ".png")
}

bc_de <- function(pre) {
    suppressMessages(mod_rnaseq_de(pre, list(contrasts = pre$contrasts),
                                   bc_config(enabled = FALSE)))
}

# --- Methods sentence ---------------------------------------------------------

test_that("rna_batch_correction_sentence describes the correction that was recorded", {
    expect_null(rna_batch_correction_sentence(NULL))

    pre <- make_bc_pre()
    limma_rec <- compute_rna_batch_correction(pre, bc_config(cols = c("Batch", "Run")))$record
    txt <- rna_batch_correction_sentence(limma_rec)
    expect_match(txt, "limma::removeBatchEffect", fixed = TRUE)
    expect_match(txt, "batch columns Batch and Run", fixed = TRUE)
    expect_match(txt, "group column Condition kept in the design", fixed = TRUE)
    expect_match(txt, "variance-stabilized (VST)", fixed = TRUE)
    expect_match(txt, "only for the outputs labelled batch-corrected", fixed = TRUE)
    expect_match(txt, "do not show that the remaining variation is biological", fixed = TRUE)
    expect_no_match(txt, "left uncorrected", fixed = TRUE)
    expect_no_match(txt, "ComBat", fixed = TRUE)

    # ComBat records are built by hand, so this needs no sva.
    combat_rec <- limma_rec
    combat_rec$method <- "combat"
    combat_rec$batch_cols <- "Batch"
    combat_rec$mean_only <- TRUE
    combat_rec$input_matrix <- "TMMlogCPM"
    combat_rec$n_features_unadjusted_qc <- 3L
    txt <- rna_batch_correction_sentence(combat_rec)
    expect_match(txt, "sva::ComBat (parametric empirical Bayes, batch means only)", fixed = TRUE)
    expect_match(txt, "batch column Batch)", fixed = TRUE)
    expect_match(txt, "TMM log-CPM", fixed = TRUE)
    expect_match(txt, "Up to 3 gene(s)", fixed = TRUE)
    combat_rec$mean_only <- FALSE
    expect_no_match(rna_batch_correction_sentence(combat_rec), "batch means only", fixed = TRUE)
})

test_that("the Methods sentence follows the record file, not the config", {
    dir <- withr::local_tempdir()
    record <- file.path(dir, "batch_correction.tsv")
    # What report_template.Rmd does.
    methods_text <- function() {
        rna_batch_correction_sentence(
            if (file.exists(record)) read.delim(record, stringsAsFactors = FALSE))
    }
    expect_null(methods_text())
    write_rna_batch_correction(compute_rna_batch_correction(make_bc_pre(), bc_config()), dir)
    expect_match(methods_text(), "batch column Batch)", fixed = TRUE)
    write_rna_batch_correction(list(ran = FALSE), dir)
    expect_null(methods_text())
})

# --- Plots --------------------------------------------------------------------

test_that("the established QC plots keep their names; the corrected set is suffixed", {
    skip_if_not_installed("pheatmap")
    skip_if_not_installed("patchwork")
    out_dir <- withr::local_tempdir()
    diag <- file.path(out_dir, "Diagnostic_plots")
    pre <- make_bc_pre()
    cfg <- bc_config()
    de_res <- bc_de(pre)

    qc_pre <- suppressMessages(mod_rnaseq_qc_pre(pre, cfg, out_dir))
    expect_setequal(basename(qc_pre$files),
                    c(bc_structure_files(), "rna_histograms_summary.png",
                      "raw_boxplot.png", "norm_boxplot.png"))
    expect_s3_class(qc_pre$pca_scores, "data.frame")
    top_de <- file.path(diag, "heatmap_top_DE.png")
    sumdf <- build_rnaseq_summary_df(de_res$tables, cfg$modes$rna$de)
    expect_identical(
        suppressMessages(draw_rna_top_de_heatmap(pre$expr_work, pre$meta, sumdf,
                                                 cfg$modes$rna, top_de)),
        top_de)

    originals <- list.files(diag, full.names = TRUE, recursive = TRUE)
    originals <- originals[!dir.exists(originals)]
    before <- tools::md5sum(originals)

    corrected <- compute_rna_batch_correction(pre, cfg)
    files <- suppressMessages(mod_rnaseq_batch_qc(pre, corrected, de_res, cfg, out_dir))
    expect_setequal(basename(files),
                    c(bc_structure_files("_batch_corr_limma"),
                      "heatmap_top_DE_batch_corr_limma.png",
                      "PCA_batch_before_after.png"))
    expect_true(all(file.exists(files)))
    # The uncorrected plots are not redrawn or replaced.
    expect_identical(tools::md5sum(originals), before)

    # Correction off in the same directory: only the corrected set goes.
    expect_identical(
        suppressMessages(mod_rnaseq_batch_qc(pre, list(ran = FALSE), de_res, cfg, out_dir)),
        character(0))
    expect_false(any(file.exists(files)))
    expect_identical(tools::md5sum(originals), before)
})

test_that("a change of method does not leave the other method's plots behind", {
    skip_if_not_installed("pheatmap")
    skip_if_not_installed("patchwork")
    out_dir <- withr::local_tempdir()
    diag <- file.path(out_dir, "Diagnostic_plots")
    pre <- make_bc_pre()
    cfg <- bc_config()
    corrected <- compute_rna_batch_correction(pre, cfg)

    dir.create(diag, recursive = TRUE)
    leftover <- file.path(diag, "PCA_PC1.vs.PC2_batch_corr_combat.png")
    file.create(leftover)
    suppressMessages(mod_rnaseq_batch_qc(pre, corrected, list(tables = list()), cfg, out_dir))
    expect_false(file.exists(leftover))
    expect_true(file.exists(file.path(diag, "PCA_PC1.vs.PC2_batch_corr_limma.png")))
})

test_that("draw_rna_top_de_heatmap writes nothing without DE genes", {
    pre <- make_bc_pre()
    out_file <- file.path(withr::local_tempdir(), "heatmap_top_DE.png")
    no_de <- data.frame(FeatureID = bc_ids(), pass_any_contrast = 0)
    expect_identical(
        draw_rna_top_de_heatmap(pre$expr_work, pre$meta, no_de,
                                bc_config()$modes$rna, out_file),
        character(0))
    expect_false(file.exists(out_file))
})

# --- Target tracking ----------------------------------------------------------

test_that("the pipeline tracks the batch-corrected QC plots as files", {
    skip_if_not_installed("targets")
    withr::local_package("targets")
    formats <- vapply(pipe_rnaseq(), function(t) t$settings$format, character(1))
    names(formats) <- vapply(pipe_rnaseq(), function(t) t$settings$name, character(1))
    expect_identical(formats[["rna_batch_qc_files"]], "file")
    expect_identical(formats[["rna_batch_corr_files"]], "file")
})

test_that("a deleted batch-corrected plot is redrawn, and a disabled run tracks no files", {
    skip_if_not_installed("targets")
    skip_if_not_installed("pheatmap")
    skip_if_not_installed("patchwork")
    sandbox <- withr::local_tempdir()
    withr::local_dir(sandbox)
    out_dir <- file.path(sandbox, "rna")
    input <- file.path(sandbox, "bc_qc_input.rds")
    pre <- make_bc_pre()
    cfg <- bc_config()
    args <- list(pre = pre, corrected = compute_rna_batch_correction(pre, cfg),
                 de_res = bc_de(pre), config = cfg)
    saveRDS(args, input)

    # The functions come from the session (helper.R sourced R/), as tar_make()
    # is run in this process.
    writeLines(c(
        "library(targets)",
        "list(",
        sprintf("  tar_target(input_file, %s, format = \"file\"),", deparse(input)),
        "  tar_target(x, readRDS(input_file)),",
        "  tar_target(",
        "    qc_files,",
        sprintf("    mod_rnaseq_batch_qc(x$pre, x$corrected, x$de_res, x$config, %s),",
                deparse(out_dir)),
        "    format = \"file\"",
        "  )",
        ")"
    ), "_targets.R")
    make <- function() {
        suppressMessages(targets::tar_make(callr_function = NULL, reporter = "silent"))
    }
    outdated <- function() {
        suppressMessages(targets::tar_outdated(callr_function = NULL, reporter = "silent"))
    }

    make()
    files <- targets::tar_read(qc_files)
    expect_true("PCA_batch_before_after.png" %in% basename(files))
    expect_true(all(file.exists(files)))
    expect_length(outdated(), 0)

    file.remove(files[basename(files) == "PCA_batch_before_after.png"])
    expect_true("qc_files" %in% outdated())
    make()
    expect_true(all(file.exists(files)))

    # Correction off: the target still succeeds, with nothing to track.
    args$corrected <- list(ran = FALSE)
    saveRDS(args, input)
    make()
    expect_identical(targets::tar_read(qc_files), character(0))
    expect_false(any(file.exists(files)))
})

# --- Before/after PCA panel ---------------------------------------------------

test_that("the PCA panel shows the group the correction kept, labelled as such", {
    # effects$color is the group: one row for it, then the batch.
    same <- rna_batch_pca_panel_columns(bc_config()$modes$rna, "Condition", "Batch")
    expect_identical(same$column, c("Condition", "Batch"))
    expect_identical(same$label, c("Condition (group kept by the correction)", "Batch (batch)"))

    # effects$color is another column: the kept group still comes first, and
    # the colour column is shown under its own name, not as the kept group.
    cfg <- bc_config(cols = c("Batch", "Run"))$modes$rna
    cfg$effects$color <- c("Pair", "Condition")
    other <- rna_batch_pca_panel_columns(cfg, "Condition", c("Batch", "Run"))
    expect_identical(other$column, c("Condition", "Pair", "Batch", "Run"))
    expect_identical(other$label, c("Condition (group kept by the correction)", "Pair",
                                    "Batch (batch)", "Run (batch)"))
    expect_identical(sum(grepl("kept by the correction", other$label)), 1L)
})

test_that("the panel is drawn when the kept group differs from effects$color", {
    skip_if_not_installed("pheatmap")
    skip_if_not_installed("patchwork")
    out_dir <- withr::local_tempdir()
    pre <- make_bc_pre()
    cfg <- bc_config()
    cfg$modes$rna$effects$color <- "Run"
    corrected <- compute_rna_batch_correction(pre, cfg)
    # The correction follows the contrasts, not the colour column.
    expect_identical(corrected$record$group_col, "Condition")
    files <- suppressMessages(mod_rnaseq_batch_qc(pre, corrected, bc_de(pre), cfg, out_dir))
    expect_true("PCA_batch_before_after.png" %in% basename(files))
})

test_that("rna_shared_pca_features keeps genes usable on both sides", {
    samples <- paste0("S", 1:4)
    before <- rbind(
        ok           = c(1, 2, 3, 4),
        flat_before  = c(5, 5, 5, 5),
        flat_after   = c(1, 2, 3, 4),
        one_finite   = c(1, NA, NA, Inf),
        only_before  = c(4, 3, 2, 1)
    )
    colnames(before) <- samples
    after <- before[c("ok", "flat_before", "flat_after", "one_finite"), ]
    after["flat_before", ] <- c(1, 2, 3, 4)
    after["flat_after", ] <- 7
    after <- rbind(after, only_after = c(1, 2, 3, 4))

    expect_identical(rna_shared_pca_features(before, after), "ok")
    # A gene that lost a value on one side still counts while two remain.
    after["ok", 1] <- NA
    expect_identical(rna_shared_pca_features(before, after), "ok")
})

test_that("both sides of the panel are fitted on the same genes", {
    skip_if_not_installed("patchwork")
    pre <- make_bc_pre()
    cfg <- bc_config()
    corrected <- compute_rna_batch_correction(pre, cfg)
    before <- bc_plain(pre$expr_work)
    after <- corrected$qc_expr
    # A gene that varies before and is constant after would be dropped from
    # one PCA only.
    after["G200", ] <- 0
    fitted <- list()
    # Stand in for the PCA plot to record the genes each fit receives. The
    # panel function looks it up where helper.R sourced R/ (its own
    # environment), so that is where it is replaced, and restored afterwards.
    fn_env <- environment(draw_rna_batch_pca_panel)
    original <- get("qc_pca_scatter", envir = fn_env)
    withr::defer(assign("qc_pca_scatter", original, envir = fn_env))
    assign("qc_pca_scatter", function(expr_mat, ...) {
        fitted[[length(fitted) + 1]] <<- rownames(expr_mat)
        ggplot2::ggplot()
    }, envir = fn_env)
    out_file <- file.path(withr::local_tempdir(), "panel.png")
    draw_rna_batch_pca_panel(before, after, pre$meta, cfg$modes$rna, "Condition", "Batch",
                             "limma", out_file)
    expect_length(fitted, 4)
    expect_true(all(vapply(fitted, identical, logical(1), fitted[[1]])))
    expect_false("G200" %in% fitted[[1]])
    expect_identical(length(fitted[[1]]), 199L)
})
