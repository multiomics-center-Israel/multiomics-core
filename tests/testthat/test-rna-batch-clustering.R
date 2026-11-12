# Tests for the clustering of the batch-corrected RNA matrix
# (R/modules/rnaseq/04b_mod_clustering_batch.R)
# test-rna-batch-clustering.R
#
# All data here is synthetic; the fixtures are in helper-rna-batch.R.

bc_clustering_config <- function(enabled = TRUE, group_col = "Condition") {
    cfg <- bc_config()
    cfg$modes$rna$clustering <- list(
        enabled   = enabled,
        group_col = group_col,
        steps = list(
            hierarchical    = list(enabled = TRUE, k = 3),
            partition       = list(enabled = FALSE),
            binary_patterns = list(enabled = FALSE)
        )
    )
    cfg
}

test_that("compare_rna_cluster_assignments lines the two runs up by gene", {
    cmp <- compare_rna_cluster_assignments(c(g2 = 2, g1 = 1, g3 = 2), c(g1 = 1, g2 = 1, g4 = 3))
    expect_identical(cmp$gene, c("g1", "g2", "g3", "g4"))
    expect_identical(cmp$cluster_uncorrected, c(1, 2, 2, NA))
    expect_identical(cmp$cluster_batch_corrected, c(1, 1, NA, 3))

    expect_null(compare_rna_cluster_assignments(NULL, c(g1 = 1)))
    expect_null(compare_rna_cluster_assignments(c(g1 = 1), NULL))
    expect_null(compare_rna_cluster_assignments(c(1, 2), c(g1 = 1)))
})

test_that("no correction ran: nothing is clustered and earlier outputs are removed", {
    out_dir <- withr::local_tempdir()
    old_dir <- file.path(out_dir, "Batch_corrected_limma", "Clustering", "Hierarchical")
    dir.create(old_dir, recursive = TRUE)
    dir.create(file.path(out_dir, "Datasets"))
    stale <- c(file.path(old_dir, "Hierarchical_DE_heatmap.png"),
               file.path(out_dir, "Datasets", "clustering_batch_comparison.tsv"))
    file.create(stale)
    # The original clustering directory is not this module's to clean.
    dir.create(file.path(out_dir, "Clustering"))
    original <- file.path(out_dir, "Clustering", "Hierarchical_DE_heatmap.png")
    file.create(original)

    pre <- make_bc_pre()
    res <- suppressMessages(mod_rnaseq_clustering_batch_corr(
        pre, list(ran = FALSE), list(tables = list()), list(objects = list()),
        bc_clustering_config(), out_dir))
    expect_identical(res$status, "not_run")
    expect_null(res$record)
    expect_identical(res$files, character(0))
    expect_false(any(file.exists(stale)))
    # The emptied directories go too, so an old one cannot stand in for a new.
    expect_false(dir.exists(file.path(out_dir, "Batch_corrected_limma")))
    expect_true(file.exists(original))
    expect_false(file.exists(file.path(out_dir, "Datasets", "clustering_batch_correction.tsv")))
})

test_that("clustering disabled: the corrected run is skipped like the original", {
    out_dir <- withr::local_tempdir()
    pre <- make_bc_pre()
    cfg <- bc_clustering_config(enabled = FALSE)
    corrected <- compute_rna_batch_correction(pre, cfg)
    res <- suppressMessages(mod_rnaseq_clustering_batch_corr(
        pre, corrected, list(tables = list()), list(objects = list()), cfg, out_dir))
    # Disabled clustering is not a skip of the corrected run: nothing was
    # attempted for either matrix, so nothing is recorded.
    expect_identical(res$status, "not_run")
    expect_null(res$comparison)
    expect_length(res$files, 0)
    expect_length(list.files(file.path(out_dir, "Datasets")), 0)
})

test_that("the corrected matrix is what is clustered, next to an untouched original", {
    skip_if_not_installed("pheatmap")
    out_dir <- withr::local_tempdir()
    pre <- make_bc_pre()
    before <- pre
    cfg <- bc_clustering_config()
    de_res <- suppressMessages(mod_rnaseq_de(pre, list(contrasts = pre$contrasts), cfg))
    corrected <- compute_rna_batch_correction(pre, cfg)

    original <- suppressMessages(mod_rnaseq_clustering(pre, de_res, cfg, out_dir))
    original_files <- original$files[file.exists(original$files)]
    expect_gt(length(original_files), 0)
    md5_before <- tools::md5sum(original_files)

    res <- suppressMessages(mod_rnaseq_clustering_batch_corr(
        pre, corrected, de_res, original, cfg, out_dir))
    expect_identical(pre, before)
    expect_identical(res$status, "completed")
    expect_identical(res$record$relation, "same")
    expect_true(is.na(res$record$reason))

    # Its own directory, same file names; the original files are not rewritten.
    corrected_root <- file.path(out_dir, "Batch_corrected_limma", "Clustering")
    clustering_files <- res$files[!startsWith(res$files, file.path(out_dir, "Datasets"))]
    expect_setequal(basename(setdiff(res$files, clustering_files)),
                    c("clustering_batch_comparison.tsv", "clustering_batch_correction.tsv"))
    expect_true(all(startsWith(clustering_files, corrected_root)))
    expect_setequal(basename(clustering_files), basename(original$files))
    expect_identical(tools::md5sum(original_files), md5_before)

    # The z-scores behind the corrected heatmap come from the corrected values.
    z_corrected <- res$clustering$excel_order$zscore_mat
    ids <- intersect(original$excel_order$ordered_ids, rownames(z_corrected))
    expect_gt(length(ids), 1)
    expect_equal(unname(z_corrected[ids, ]),
                 unname(zscore_rows(corrected$expr_work[ids, ])))
    expect_false(isTRUE(all.equal(unname(z_corrected[ids, ]),
                                  unname(original$excel_order$zscore_mat[ids, ]))))

    # Both assignments are kept, for the same genes.
    cmp <- res$comparison
    expect_named(cmp, c("gene", "cluster_uncorrected", "cluster_batch_corrected"))
    expect_setequal(cmp$gene, names(original$objects$clusters))
    expect_false(anyNA(cmp))
    written <- read.delim(file.path(out_dir, "Datasets", "clustering_batch_comparison.tsv"),
                          stringsAsFactors = FALSE)
    expect_identical(written$gene, cmp$gene)
})

# --- The clustering group and what the correction kept ------------------------

test_that("classify_rna_clustering_group tells the five relations apart", {
    meta <- bc_meta()
    classify <- function(kept, clustering, kept_col, clustering_col, batch_cols = "Batch") {
        classify_rna_clustering_group(kept, clustering, lapply(batch_cols, function(b) meta[[b]]),
                                      kept_col, clustering_col, batch_cols)
    }

    same <- classify(meta$Condition, meta$Condition, "Condition", "Condition")
    expect_identical(same, list(relation = "same", proceed = TRUE))

    # A coarser column: every kept group lies in one clustering group.
    fine <- paste(meta$Condition, meta$Run)
    determined <- classify(fine, meta$Condition, "Fine", "Condition")
    expect_identical(determined$relation, "determined_by_kept_group")
    expect_true(determined$proceed)

    batch_col <- classify(meta$Condition, meta$Batch, "Condition", "Batch")
    expect_identical(batch_col, list(relation = "batch_column", proceed = FALSE))
    second_batch <- classify(meta$Condition, meta$Run, "Condition", "Run",
                             batch_cols = c("Batch", "Run"))
    expect_identical(second_batch$relation, "batch_column")

    # Within each condition, this column is the batch under another name.
    site <- paste(meta$Condition, meta$Batch)
    confounded <- classify(meta$Condition, site, "Condition", "Site")
    expect_identical(confounded, list(relation = "confounded_with_batch", proceed = FALSE))

    # Crossed with batch: identifiable, but not in the correction's model.
    crossed <- classify(meta$Condition, meta$Run, "Condition", "Run")
    expect_identical(crossed, list(relation = "not_protected", proceed = TRUE))

    # Imbalance is not confounding: this column follows the batch except for
    # one sample, so it is heavily unbalanced yet still identifiable.
    nearly_batch <- meta$Batch
    nearly_batch[1] <- "B2"
    unbalanced <- classify(meta$Condition, nearly_batch, "Condition", "NearlyBatch")
    expect_identical(unbalanced$relation, "not_protected")
})

test_that("rna_batch_clustering_note says why it was skipped or what was not protected", {
    rec <- function(status, reason, relation, method = "limma", clustering_col = "Line") {
        rna_batch_clustering_record(status, reason, method, "Condition", clustering_col,
                                    relation, "Batch")
    }
    expect_null(rna_batch_clustering_note(NULL))
    expect_null(rna_batch_clustering_note(rec("completed", NA_character_, "same")))
    expect_null(rna_batch_clustering_note(
        rec("completed", NA_character_, "determined_by_kept_group")))

    skipped <- rna_batch_clustering_note(rec("skipped", "batch_column", "batch_column",
                                             clustering_col = "Batch"))
    expect_match(skipped, "was skipped", fixed = TRUE)
    expect_match(skipped, "'Batch' is a batch column", fixed = TRUE)
    confounded <- rna_batch_clustering_note(
        rec("skipped", "confounded_with_batch", "confounded_with_batch"))
    expect_match(confounded, "cannot be told apart from batch", fixed = TRUE)
    expect_match(confounded, "within the groups of 'Condition'", fixed = TRUE)
    expect_match(rna_batch_clustering_note(rec("skipped", "no_output", "same")),
                 "produced no output", fixed = TRUE)

    limma_note <- rna_batch_clustering_note(rec("completed", NA_character_, "not_protected"))
    expect_match(limma_note, "not explicitly protected", fixed = TRUE)
    # Identifiable or balanced does not mean preserved.
    expect_match(limma_note, "also where it is balanced across batches", fixed = TRUE)
    expect_no_match(limma_note, "ComBat", fixed = TRUE)
    combat_note <- rna_batch_clustering_note(
        rec("completed", NA_character_, "not_protected", method = "combat"))
    expect_match(combat_note, "ComBat also adjusts each batch's spread", fixed = TRUE)
})

test_that("an unusable clustering group skips the corrected clustering without failing", {
    pre <- make_bc_pre()
    pre$meta$Site <- paste(pre$meta$Condition, pre$meta$Batch)
    cases <- list(Batch = "batch_column", Site = "confounded_with_batch")
    for (group_col in names(cases)) {
        out_dir <- withr::local_tempdir()
        cfg <- bc_clustering_config(group_col = group_col)
        corrected <- compute_rna_batch_correction(pre, cfg)
        expect_message(
            res <- mod_rnaseq_clustering_batch_corr(
                pre, corrected, list(tables = list()), list(objects = list()), cfg, out_dir),
            "was skipped")
        expect_identical(res$status, "skipped")
        expect_identical(res$record$reason, cases[[group_col]])
        expect_null(res$clustering)

        # The skip is on record, with its reason, and nothing was clustered.
        record_file <- file.path(out_dir, "Datasets", "clustering_batch_correction.tsv")
        expect_identical(res$files, record_file)
        written <- read.delim(record_file, stringsAsFactors = FALSE)
        expect_identical(written$status, "skipped")
        expect_identical(written$reason, cases[[group_col]])
        expect_identical(written$kept_group_col, "Condition")
        expect_identical(written$clustering_group_col, group_col)
        expect_false(dir.exists(file.path(out_dir, "Batch_corrected_limma")))
    }
})

test_that("an unprotected clustering group runs with a warning and is recorded", {
    skip_if_not_installed("pheatmap")
    out_dir <- withr::local_tempdir()
    pre <- make_bc_pre()
    cfg <- bc_clustering_config(group_col = "Run")
    de_res <- suppressMessages(mod_rnaseq_de(pre, list(contrasts = pre$contrasts), cfg))
    corrected <- compute_rna_batch_correction(pre, cfg)
    original <- suppressMessages(mod_rnaseq_clustering(pre, de_res, cfg, out_dir))

    expect_warning(
        res <- suppressMessages(mod_rnaseq_clustering_batch_corr(
            pre, corrected, de_res, original, cfg, out_dir)),
        "not explicitly protected")
    expect_identical(res$status, "completed")
    expect_identical(res$record$relation, "not_protected")
    written <- read.delim(file.path(out_dir, "Datasets", "clustering_batch_correction.tsv"),
                          stringsAsFactors = FALSE)
    expect_identical(written$relation, "not_protected")
    expect_match(rna_batch_clustering_note(written), "not explicitly protected", fixed = TRUE)
})

# --- Stale directories and report discovery -----------------------------------

test_that("the directory of an earlier K does not survive or stand in for the new one", {
    out_dir <- withr::local_tempdir()
    root <- file.path(out_dir, "Batch_corrected_limma", "Clustering")
    old_k <- file.path(root, "Partition_clustering", "Partition_clustering_7_clusters")
    dir.create(old_k, recursive = TRUE)
    file.create(file.path(old_k, "Partition_clustering_heatmap.png"))
    file.create(file.path(old_k, "partition_clusters.tsv"))

    # Cleanup: the old K's files and its directory go.
    remove_rna_batch_clustering_outputs(out_dir)
    expect_false(dir.exists(old_k))
    expect_false(dir.exists(file.path(out_dir, "Batch_corrected_limma")))

    # Discovery: an empty directory left by anything else is not picked.
    new_k <- file.path(root, "Partition_clustering", "Partition_clustering_4_clusters")
    dir.create(old_k, recursive = TRUE)
    expect_identical(rna_partition_heatmap_dir(root), NA_character_)
    dir.create(new_k)
    file.create(file.path(new_k, "Partition_clustering_heatmap.png"))
    expect_identical(rna_partition_heatmap_dir(root), new_k)

    # A directory that still holds a file of someone else's is kept.
    foreign <- file.path(out_dir, "Clustering", "Partition_clustering")
    dir.create(foreign, recursive = TRUE)
    file.create(file.path(foreign, "kept.txt"))
    remove_rna_batch_clustering_outputs(out_dir)
    expect_true(file.exists(file.path(foreign, "kept.txt")))
    expect_false(dir.exists(root))
})

test_that("the report view of the corrected clustering does not depend on the uncorrected one", {
    rna_dir <- withr::local_tempdir()
    rec <- function(status, reason, relation) {
        rna_batch_clustering_record(status, reason, "limma", "Condition", "Run", relation, "Batch")
    }
    expect_false(rna_batch_clustering_view(NULL, rna_dir)$show)

    # A skip is shown, with its reason, although there is no figure at all.
    skipped <- rna_batch_clustering_view(rec("skipped", "batch_column", "batch_column"), rna_dir)
    expect_true(skipped$show)
    expect_match(skipped$note, "was skipped", fixed = TRUE)
    expect_identical(skipped$binary_heatmaps, character(0))

    # Completed without a figure the report shows (hierarchical only): nothing.
    completed <- rec("completed", NA_character_, "not_protected")
    expect_false(rna_batch_clustering_view(completed, rna_dir)$show)

    # Completed with binary patterns only, and no uncorrected output anywhere.
    binary_dir <- file.path(rna_dir, "Batch_corrected_limma", "Clustering", "Binary_patterns")
    dir.create(binary_dir, recursive = TRUE)
    file.create(file.path(binary_dir, "Heatmap_01.png"))
    view <- rna_batch_clustering_view(completed, rna_dir)
    expect_true(view$show)
    expect_identical(basename(view$binary_heatmaps), "Heatmap_01.png")
    expect_true(is.na(view$partition_heatmap))
    expect_match(view$note, "not explicitly protected", fixed = TRUE)

    # The template shows the section and the tab on that view, not only when
    # the uncorrected run has figures.
    rmd <- file.path(root_dir, "R", "domain", "rnaseq", "report_template.Rmd")
    skip_if_not(file.exists(rmd), "report_template.Rmd not found")
    src <- readLines(rmd, warn = FALSE)
    discover <- find_rmd_chunk_lines(rmd, "clustering-discover-rna")
    discover <- paste(src[seq(discover[1], discover[2])], collapse = "\n")
    expect_match(discover, "rna_batch_clustering_view(bc_clust_record, rna_dir)", fixed = TRUE)
    expect_match(discover, "has_part_rna || has_binary_rna || bc_clust_show_rna", fixed = TRUE)
    expect_true(any(grepl("{r clustering-batch-corrected-rna, eval=bc_clust_show_rna,", src,
                          fixed = TRUE)))
})

# --- Target tracking ----------------------------------------------------------

test_that("the pipeline tracks the corrected clustering files and keeps the clustering object", {
    skip_if_not_installed("targets")
    withr::local_package("targets")
    formats <- vapply(pipe_rnaseq(), function(t) t$settings$format, character(1))
    names(formats) <- vapply(pipe_rnaseq(), function(t) t$settings$name, character(1))
    expect_identical(formats[["rna_clustering_batch_corr_files"]], "file")
    # The canonical clustering target is still an ordinary object target.
    expect_false(identical(formats[["rna_clustering_obj"]], "file"))
})

test_that("a deleted corrected-clustering file is rebuilt, and no correction tracks no files", {
    skip_if_not_installed("targets")
    skip_if_not_installed("pheatmap")
    sandbox <- withr::local_tempdir()
    withr::local_dir(sandbox)
    out_dir <- file.path(sandbox, "rna")
    input <- file.path(sandbox, "bc_clustering_input.rds")
    pre <- make_bc_pre()
    cfg <- bc_clustering_config()
    de_res <- suppressMessages(mod_rnaseq_de(pre, list(contrasts = pre$contrasts), cfg))
    args <- list(pre = pre, corrected = compute_rna_batch_correction(pre, cfg), de_res = de_res,
                 original = suppressMessages(mod_rnaseq_clustering(pre, de_res, cfg, out_dir)),
                 config = cfg)
    saveRDS(args, input)

    # The functions come from the session (helper.R sourced R/), as tar_make()
    # is run in this process.
    writeLines(c(
        "library(targets)",
        "list(",
        sprintf("  tar_target(input_file, %s, format = \"file\"),", deparse(input)),
        "  tar_target(x, readRDS(input_file)),",
        "  tar_target(",
        "    clustering_files,",
        "    mod_rnaseq_clustering_batch_corr(x$pre, x$corrected, x$de_res, x$original,",
        sprintf("                                     x$config, %s)$files,", deparse(out_dir)),
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
    files <- targets::tar_read(clustering_files)
    record <- file.path(out_dir, "Datasets", "clustering_batch_correction.tsv")
    expect_true(record %in% files)
    expect_true(any(startsWith(files, file.path(out_dir, "Batch_corrected_limma"))))
    expect_true(all(file.exists(files)))
    expect_length(outdated(), 0)

    # The record is what the report reads: losing it must rebuild the target.
    file.remove(record)
    expect_true("clustering_files" %in% outdated())
    make()
    expect_true(all(file.exists(files)))

    args$corrected <- list(ran = FALSE)
    saveRDS(args, input)
    make()
    expect_identical(targets::tar_read(clustering_files), character(0))
    expect_false(any(file.exists(files)))
})
