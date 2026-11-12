# Tests for batch correction of the RNA display matrices
# (R/domain/rnaseq/04c_batch_correction.R)
# test-rna-batch-correction.R
#
# All data here is synthetic.

# The fixtures (bc_meta(), make_bc_pre(), bc_config(), bc_shift(), ...) are in
# helper-rna-batch.R, shared with test-rna-batch-qc.R.

# --- Config -------------------------------------------------------------------

test_that("get_rna_batch_correction_config is off unless enabled, limma by default", {
    off <- list(enabled = FALSE, method = NA_character_, batch_cols = character(0),
                mean_only = NA)
    expect_identical(get_rna_batch_correction_config(list()), off)
    expect_identical(get_rna_batch_correction_config(
        list(batch_correction = list(enabled = FALSE, method = "nonsense"))), off)

    on <- function(...) get_rna_batch_correction_config(
        list(batch_correction = list(enabled = TRUE, ...)))
    expect_identical(on(batch_cols = "Batch")$method, "limma")
    expect_identical(on(batch_cols = c("Batch", "Run"))$batch_cols, c("Batch", "Run"))
    expect_identical(on(batch_cols = "Batch", combat_mean_only = TRUE)$mean_only, NA)
    expect_false(on(batch_cols = "Batch", method = "combat")$mean_only)
    expect_true(on(batch_cols = "Batch", method = "combat", combat_mean_only = TRUE)$mean_only)

    for (bad in list(character(0), c("A", "B", "C"), c("Batch", "Batch"), "")) {
        expect_error(on(batch_cols = bad), "batch_cols")
    }
    expect_error(on(batch_cols = "Batch", method = "sva"), "limma.*combat")
    expect_error(on(batch_cols = c("Batch", "Run"), method = "combat"), "exactly one")
})

test_that("validate_rna_config checks the batch_correction block", {
    cfg <- list(
        id_columns = list(gene_id = "gene"),
        files = list(counts = "counts.tsv", metadata = "meta.tsv")
    )
    expect_true(validate_rna_config(cfg))
    cfg$batch_correction <- list(enabled = TRUE, batch_cols = "Batch")
    expect_true(validate_rna_config(cfg))

    cfg$batch_correction$combat_mean_only <- TRUE
    expect_warning(validate_rna_config(cfg), "combat_mean_only")
    cfg$batch_correction$method <- "combat"
    expect_true(validate_rna_config(cfg))

    cfg$batch_correction$batch_cols <- c("Batch", "Run")
    expect_error(validate_rna_config(cfg), "exactly one")
    cfg$batch_correction <- "limma"
    expect_error(validate_rna_config(cfg), "mapping")
})

test_that("resolve_rna_de_factor_col follows the contrasts, then the DE fallback", {
    con <- data.frame(Factor = c("Genotype", "Genotype"))
    expect_identical(resolve_rna_de_factor_col(list(effects = list(color = "Group")), con),
                     "Genotype")
    expect_error(resolve_rna_de_factor_col(list(), data.frame(Factor = c("A", "B"))),
                 "same factor")
    expect_identical(resolve_rna_de_factor_col(
        list(effects = list(color = c("Group", "Batch")),
             filtering = list(group_col = "Other")), NULL), "Group")
    expect_identical(resolve_rna_de_factor_col(
        list(filtering = list(group_col = "Other")), NULL), "Other")
    expect_identical(resolve_rna_de_factor_col(list(), NULL), "Condition")
})

# --- Validation of the columns ------------------------------------------------

test_that("prepare_rna_batch_correction aligns to the matrix and makes batch categorical", {
    meta <- bc_meta()
    ids <- rev(meta$SampleID)
    out <- prepare_rna_batch_correction(meta, ids, "SampleID", "Condition", c("Batch", "Pair"))
    expect_identical(as.character(out$group), rev(meta$Condition))
    expect_identical(as.character(out$batch$Batch), rev(meta$Batch))
    expect_s3_class(out$batch$Pair, "factor")
    expect_identical(nlevels(out$batch$Pair), 3L)

    # Row names stand in for a missing sample column.
    by_rowname <- prepare_rna_batch_correction(meta[, -1], ids, "SampleID", "Condition", "Batch")
    expect_identical(by_rowname$batch$Batch, out$batch$Batch)
})

test_that("prepare_rna_batch_correction stops on unusable columns", {
    meta <- bc_meta()
    prep <- function(m, cols, ids = meta$SampleID, group = "Condition") {
        prepare_rna_batch_correction(m, ids, "SampleID", group, cols)
    }
    expect_error(prep(meta, "Nope"), "not found in metadata")
    expect_error(prep(meta, "Batch", group = "Nope"), "not found in metadata")
    expect_error(prep(meta, "Condition"), "is the group column")
    expect_error(prep(meta, "Batch", ids = c(meta$SampleID, "S99")), "not in the metadata")

    spaced <- meta
    spaced[["Batch no"]] <- spaced$Batch
    expect_error(prep(spaced, "Batch no"), "not syntactic")

    holed <- meta
    holed$Batch[3] <- NA
    expect_error(prep(holed, "Batch"), "has no value for 1 sample")
    holed$Batch[3] <- ""
    expect_error(prep(holed, "Batch"), "has no value for 1 sample")
    no_group <- meta
    no_group$Condition[5] <- NA
    expect_error(prep(no_group, "Batch"), "Group column")

    single <- meta
    single$Batch <- "B1"
    expect_error(prep(single, "Batch"), "single level")

    confounded <- meta
    confounded$Batch <- ifelse(confounded$Condition == "Ctrl", "B1", "B2")
    expect_error(prep(confounded, "Batch"), "confounded")
})

# --- The correction -----------------------------------------------------------

test_that("disabled or absent: nothing is computed and nothing is written", {
    pre <- make_bc_pre()
    expect_identical(compute_rna_batch_correction(pre, bc_config(enabled = FALSE)),
                     list(ran = FALSE))
    no_block <- bc_config()
    no_block$modes$rna$batch_correction <- NULL
    expect_identical(compute_rna_batch_correction(pre, no_block), list(ran = FALSE))

    dir <- withr::local_tempdir()
    expect_identical(write_rna_batch_correction(list(ran = FALSE), dir), character(0))
    expect_length(list.files(dir), 0)
})

test_that("limma removes the batch shift, keeps the group effect and leaves pre alone", {
    pre <- make_bc_pre()
    before <- pre
    res <- compute_rna_batch_correction(pre, bc_config())
    expect_identical(pre, before)
    expect_true(res$ran)

    orig <- bc_plain(pre$expr_work)
    expect_identical(dimnames(res$expr_work), dimnames(orig))
    expect_false(isTRUE(all.equal(res$expr_work, orig)))
    # VST has no QC matrix of its own.
    expect_identical(res$qc_expr, res$expr_work)

    sgn <- bc_batch_sign()[bc_batch_sign() != 0]
    expect_gt(bc_shift(orig, "Batch", bc_batch_genes(), sgn), 1.3)
    expect_lt(abs(bc_shift(res$expr_work, "Batch", bc_batch_genes(), sgn)), 0.1)

    # The group effect is kept where it exists, and the batch no longer leaks
    # into the group comparison of genes that have none.
    expect_equal(bc_shift(res$expr_work, "Condition", bc_group_genes()), 2, tolerance = 0.1)
    expect_gt(abs(bc_shift(orig, "Condition", bc_batch_genes(), sgn)), 0.3)
    expect_lt(abs(bc_shift(res$expr_work, "Condition", bc_batch_genes(), sgn)), 0.1)

    rec <- res$record
    expect_identical(nrow(rec), 1L)
    expect_identical(rec$method, "limma")
    expect_identical(rec$batch_cols, "Batch")
    expect_identical(rec$group_col, "Condition")
    expect_true(is.na(rec$mean_only))
    expect_true(is.na(rec$seed))
    expect_identical(rec$input_matrix, "VST")
    expect_identical(c(rec$n_samples, rec$n_features), c(12L, 200L))
    expect_identical(rec$n_features_unadjusted, 0L)
    expect_identical(rec$output_file, "rna_batch_corr_limma_VST.tsv")
})

test_that("two batch columns are both removed, and the exported matrix is the corrected one", {
    pre <- make_bc_pre()
    orig <- bc_plain(pre$expr_work)
    one <- compute_rna_batch_correction(pre, bc_config(cols = "Batch"))
    two <- compute_rna_batch_correction(pre, bc_config(cols = c("Batch", "Run")))

    sgn <- bc_run_sign()[bc_run_sign() != 0]
    expect_gt(bc_shift(one$expr_work, "Run", bc_run_genes(), sgn), 1.3)
    expect_lt(abs(bc_shift(two$expr_work, "Run", bc_run_genes(), sgn)), 0.1)
    expect_false(isTRUE(all.equal(two$expr_work, one$expr_work)))
    expect_false(isTRUE(all.equal(two$expr_work, orig)))
    expect_identical(two$record$batch_cols, "Batch;Run")

    # The class of bug to guard: a corrected result that is computed but not
    # the one written, so the export holds the uncorrected matrix.
    dir <- withr::local_tempdir()
    files <- write_rna_batch_correction(two, dir)
    expect_identical(basename(files), c("rna_batch_corr_limma_VST.tsv", "batch_correction.tsv"))
    back <- read.delim(files[1], check.names = FALSE, stringsAsFactors = FALSE)
    expect_identical(back$gene, rownames(two$expr_work))
    written <- as.matrix(back[, -1])
    rownames(written) <- back$gene
    expect_equal(written, two$expr_work, tolerance = 1e-8)
    expect_false(isTRUE(all.equal(written, orig, tolerance = 1e-8)))
})

test_that("TMMlogCPM corrects the QC matrix and the analysis matrix separately", {
    pre <- make_bc_pre(method = "TMMlogCPM")
    cfg <- bc_config()
    res <- compute_rna_batch_correction(pre, cfg)
    qc_orig <- rna_qc_expr(pre, cfg$modes$rna$qc)

    expect_identical(dimnames(res$qc_expr), dimnames(qc_orig))
    expect_false(isTRUE(all.equal(res$qc_expr, res$expr_work)))
    expect_false(isTRUE(all.equal(res$qc_expr, qc_orig)))
    sgn <- bc_batch_sign()[bc_batch_sign() != 0]
    expect_lt(abs(bc_shift(res$qc_expr, "Batch", bc_batch_genes(), sgn)), 0.1)
    expect_match(res$record$qc_matrix, "log2(TMM CPM", fixed = TRUE)
    expect_identical(res$record$output_file, "rna_batch_corr_limma_TMMlogCPM.tsv")
})

test_that("rna_constant_within_batch flags only features constant inside a batch", {
    batch <- factor(c("B1", "B1", "B2", "B2", "B3"))
    m <- rbind(
        varies   = c(1, 2, 3, 4, 5),
        flat_b1  = c(0, 0, 3, 4, 5),
        flat_b2  = c(1, 2, 7, 7, 5),
        flat_all = c(1, 1, 1, 1, 1)
    )
    expect_identical(rna_constant_within_batch(m, batch), c(FALSE, TRUE, TRUE, TRUE))
})

test_that("ComBat corrects, records its settings and leaves within-batch constants alone", {
    skip_if_not_installed("sva")
    pre <- make_bc_pre()
    # A gene with no signal in B1, as a gene without reads there would be.
    flat <- "G150"
    pre$expr_work[flat, pre$meta$Batch == "B1"] <- 0
    orig <- bc_plain(pre$expr_work)

    cfg <- bc_config(method = "combat")
    suppressMessages(expect_warning(
        res <- compute_rna_batch_correction(pre, cfg), "constant within a batch"))
    expect_identical(res$record$n_features_unadjusted, 1L)
    expect_identical(res$expr_work[flat, ], orig[flat, ])
    expect_true(all(is.finite(res$expr_work)))
    expect_identical(dimnames(res$expr_work), dimnames(orig))

    sgn <- bc_batch_sign()[bc_batch_sign() != 0]
    expect_lt(abs(bc_shift(res$expr_work, "Batch", bc_batch_genes(), sgn)), 0.2)
    expect_equal(bc_shift(res$expr_work, "Condition", bc_group_genes()), 2, tolerance = 0.2)

    expect_identical(res$record$method, "combat")
    expect_false(res$record$mean_only)
    expect_identical(res$record$seed, 1L)

    again <- suppressMessages(suppressWarnings(compute_rna_batch_correction(pre, cfg)))
    expect_identical(again$expr_work, res$expr_work)

    mean_only <- suppressMessages(suppressWarnings(compute_rna_batch_correction(
        pre, bc_config(method = "combat", combat_mean_only = TRUE))))
    expect_true(mean_only$record$mean_only)
    expect_false(isTRUE(all.equal(mean_only$expr_work, res$expr_work)))
})

test_that("ComBat refuses a single-sample batch unless mean-only is on", {
    skip_if_not_installed("sva")
    pre <- make_bc_pre()
    pre$meta$Batch[1] <- "B3"
    expect_error(
        suppressMessages(compute_rna_batch_correction(pre, bc_config(method = "combat"))),
        "single sample")
})

# --- DE is not affected -------------------------------------------------------

test_that("DE results are the same with batch correction on and off", {
    pre <- make_bc_pre()
    inputs <- list(contrasts = pre$contrasts)
    de_off <- suppressMessages(mod_rnaseq_de(pre, inputs, bc_config(enabled = FALSE)))
    compute_rna_batch_correction(pre, bc_config())
    de_on <- suppressMessages(mod_rnaseq_de(pre, inputs, bc_config()))
    expect_identical(de_on$tables, de_off$tables)
    expect_identical(de_on$design, de_off$design)
})

# --- Files --------------------------------------------------------------------

test_that("the record and the matrix always reflect the current run", {
    dir <- withr::local_tempdir()
    pre <- make_bc_pre()
    record <- file.path(dir, "batch_correction.tsv")

    limma_run <- compute_rna_batch_correction(pre, bc_config(cols = c("Batch", "Run")))
    files <- write_rna_batch_correction(limma_run, dir)
    expect_true(all(file.exists(files)))
    # The report and the fact sheet find the normalized table by this prefix.
    expect_length(list.files(dir, pattern = "^rna_norm_"), 0)
    back <- read.delim(record, stringsAsFactors = FALSE)
    expect_identical(back$method, "limma")
    expect_identical(back$batch_cols, "Batch;Run")
    expect_identical(back$output_file, basename(files[1]))

    # Another matrix behind the correction: the earlier export must not stay.
    tmm_run <- compute_rna_batch_correction(make_bc_pre(method = "TMMlogCPM"), bc_config())
    write_rna_batch_correction(tmm_run, dir)
    expect_identical(sort(list.files(dir)),
                     c("batch_correction.tsv", "rna_batch_corr_limma_TMMlogCPM.tsv"))

    # Correction turned off in the same directory: no evidence may remain.
    expect_identical(write_rna_batch_correction(list(ran = FALSE), dir), character(0))
    expect_length(list.files(dir), 0)
})

# --- Target tracking ----------------------------------------------------------

test_that("the pipeline tracks the batch-correction exports as files", {
    skip_if_not_installed("targets")
    withr::local_package("targets")
    formats <- vapply(pipe_rnaseq(), function(t) t$settings$format, character(1))
    names(formats) <- vapply(pipe_rnaseq(), function(t) t$settings$name, character(1))
    expect_identical(formats[["rna_batch_corr_files"]], "file")
})

test_that("a deleted export is rebuilt, and a disabled run tracks no files", {
    skip_if_not_installed("targets")
    sandbox <- withr::local_tempdir()
    withr::local_dir(sandbox)
    datasets <- file.path(sandbox, "Datasets")
    input <- file.path(sandbox, "bc_res.rds")
    saveRDS(compute_rna_batch_correction(make_bc_pre(), bc_config()), input)

    # The functions come from the session (helper.R sourced R/), as tar_make()
    # is run in this process.
    writeLines(c(
        "library(targets)",
        "list(",
        sprintf("  tar_target(bc_res_file, %s, format = \"file\"),", deparse(input)),
        "  tar_target(bc_res, readRDS(bc_res_file)),",
        sprintf("  tar_target(bc_files, write_rna_batch_correction(bc_res, %s), format = \"file\")",
                deparse(datasets)),
        ")"
    ), "_targets.R")
    make <- function() {
        suppressMessages(targets::tar_make(callr_function = NULL, reporter = "silent"))
    }
    outdated <- function() {
        suppressMessages(targets::tar_outdated(callr_function = NULL, reporter = "silent"))
    }

    make()
    files <- targets::tar_read(bc_files)
    expect_identical(basename(files), c("rna_batch_corr_limma_VST.tsv", "batch_correction.tsv"))
    expect_true(all(file.exists(files)))
    expect_length(outdated(), 0)

    # The record is the report's evidence of the correction: losing it must
    # make the target outdated, and the next run must write it again.
    file.remove(files[2])
    expect_true("bc_files" %in% outdated())
    make()
    expect_true(all(file.exists(files)))

    # Correction off: the target still succeeds, with nothing to track.
    saveRDS(list(ran = FALSE), input)
    make()
    expect_identical(targets::tar_read(bc_files), character(0))
    expect_length(list.files(datasets), 0)
})
