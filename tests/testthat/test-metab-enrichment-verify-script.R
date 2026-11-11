# tests/testthat/test-metab-enrichment-verify-script.R
#
# Regression tests for the before/after verification script of the enrichment
# layout refactor (fixtures/metab_enrichment_layout/verify_layout_refactor.R).
#
# The script is a check, so its failure modes are the dangerous ones: a check
# that passes without having looked at anything. Three such cases were found on
# a real run and are pinned here:
#   - the enrichment tables were looked for under "Datasets/" spelled exactly,
#     and a folder spelled "datasets" on disk yielded none;
#   - with no table collected, the table checks compared two empty sets and
#     passed;
#   - plot data was compared with identical(), which compares references by
#     address and so fails on equal content from two R sessions.
#
# Everything is synthetic. No pipeline run, {targets} store or real result is
# involved: the snapshots are built by hand.

verify_env <- new.env()
sys.source(test_path("fixtures", "metab_enrichment_layout",
                     "verify_layout_refactor.R"), envir = verify_env)

# A matching before/after pair: one contrast, ORA and QEA, in the old flat
# layout and the unit layout. `datasets` sets how the before run spells its
# Datasets folder.
snapshot_pair <- function(datasets = "Datasets") {
    md5 <- c(ora = "aaaa", ora_first = "aaaa", qea = "bbbb")
    before_tsv <- stats::setNames(unname(md5), c(
        file.path(datasets, "enrichment_ora_A_vs_B_results.tsv"),
        file.path(datasets, "enrichment_ora_results.tsv"),
        file.path(datasets, "enrichment_qea_results.tsv")))
    after_tsv <- stats::setNames(unname(md5), c(
        "Enrichment/ORA/combined/contrasts/A_vs_B/de/results.tsv",
        "Datasets/enrichment_ora_results.tsv",
        "Enrichment/QEA/combined/all_samples/results.tsv"))
    before_png <- c("Enrichment/enrichment_ora_dotplot_A_vs_B.png",
                    "Enrichment/enrichment_ora_dotplot.png",
                    "Enrichment/enrichment_qea_barplot.png")
    after_png <- c("Enrichment/ORA/combined/contrasts/A_vs_B/de/dotplot.png",
                   "Enrichment/enrichment_ora_dotplot.png",
                   "Enrichment/QEA/combined/all_samples/barplot.png")
    shared <- list(
        results = list(ora = list(table = data.frame(pathway = "P1", FDR = 0.01))),
        plots = list(plots = list(ora_dotplot = list(
            data = data.frame(pathway = "P1", FDR = 0.01),
            labels = list(title = "ORA")))),
        result_names = c("ora", "plots", "files"),
        by_contrast_names = list(ora_by_contrast = list(
            "A vs B" = c("res", "plots", "files"))),
        mummichog_md5 = stats::setNames(character(0), character(0)),
        payload_names = c("expr_norm", "enrichment_ora"),
        payload_tables = list(enrichment_ora = data.frame(pathway = "P1")),
        slide_titles = "Enrichment: Enrichment Ora Dotplot a vs b")
    list(
        before = c(shared, list(files = c(names(before_tsv), before_png),
                                tsv_md5 = before_tsv)),
        after = c(shared, list(
            files = c(names(after_tsv), after_png,
                      "Enrichment/enrichment_manifest.csv"),
            tsv_md5 = after_tsv)))
}

# Write a pair to a temporary folder, compare, and return the per-check results.
compare_pair <- function(pair) {
    d <- withr::local_tempdir(.local_envir = parent.frame())
    saveRDS(pair$before, file.path(d, "before.rds"))
    saveRDS(pair$after, file.path(d, "after.rds"))
    out <- NULL
    printed <- utils::capture.output(out <- verify_env$compare_enrichment_layout(d))
    list(all = isTRUE(as.logical(out)), checks = attr(out, "checks"),
         printed = printed)
}

# ---- the fixture itself -------------------------------------------------------

test_that("a matching pair of snapshots passes every check", {
    res <- compare_pair(snapshot_pair())
    expect_true(res$all)
    expect_true(all(res$checks))
    expect_identical(names(res$checks),
                     c("results", "plots", "fields", "tsv_content", "tsv_complete",
                       "png", "manifest", "mummichog", "payload_keys",
                       "payload_tables", "slides"))
})

# ---- case of the output folders -----------------------------------------------

test_that("the output folders are recognised whatever their case on disk", {
    canon <- verify_env$.el_canon_paths
    expect_identical(
        canon(c("datasets/enrichment_ora_results.tsv", "DATASETS/x.tsv",
                "enrichment/ORA/combined/contrasts/A/de/results.tsv",
                "Mummichog_Pinned/A/input.tsv", "Datasets", "datasets\\a.tsv")),
        c("Datasets/enrichment_ora_results.tsv", "Datasets/x.tsv",
          "Enrichment/ORA/combined/contrasts/A/de/results.tsv",
          "mummichog_pinned/A/input.tsv", "Datasets", "Datasets/a.tsv"))
    # Only the first component, and only a whole component, is rewritten.
    expect_identical(
        canon(c("Network/datasets/a.tsv", "datasets_old/a.tsv",
                "Enrichment/ORA/combined/contrasts/datasets/de/results.tsv")),
        c("Network/datasets/a.tsv", "datasets_old/a.tsv",
          "Enrichment/ORA/combined/contrasts/datasets/de/results.tsv"))
})

test_that("enrichment tables are found under a lower-case datasets folder", {
    is_tsv <- verify_env$.el_is_enrichment_tsv
    expect_identical(
        is_tsv(c("datasets/enrichment_ora_A_vs_B_results.tsv",
                 "Datasets/enrichment_qea_results.tsv",
                 "enrichment/QEA/combined/all_samples/results.tsv",
                 "datasets/final_results.tsv",
                 "datasets/sub/enrichment_ora_results.tsv",
                 "Enrichment/enrichment_manifest.csv",
                 "mummichog_pinned/A/enrichment_x.tsv")),
        c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
})

test_that("the inventory hashes the tables of a lower-case datasets folder", {
    mode_dir <- file.path(withr::local_tempdir(), "metabolomics")
    dir.create(file.path(mode_dir, "datasets"), recursive = TRUE)
    dir.create(file.path(mode_dir, "enrichment"))
    dir.create(file.path(mode_dir, "mummichog_pinned", "A"), recursive = TRUE)
    writeLines("ora", file.path(mode_dir, "datasets", "enrichment_ora_A_results.tsv"))
    writeLines("other", file.path(mode_dir, "datasets", "final_results.tsv"))
    writeLines("png", file.path(mode_dir, "enrichment", "enrichment_ora_dotplot_A.png"))
    writeLines("mmc", file.path(mode_dir, "mummichog_pinned", "A", "input.tsv"))
    # A contrast whose name starts with "." has a hidden unit directory.
    hidden <- file.path(mode_dir, "enrichment", "ORA", "combined", "contrasts",
                        ".treated", "de")
    dir.create(hidden, recursive = TRUE)
    writeLines("hidden", file.path(hidden, "results.tsv"))

    inv <- verify_env$.el_inventory(mode_dir)

    expect_setequal(names(inv$tsv_md5),
                    c("Datasets/enrichment_ora_A_results.tsv",
                      "Enrichment/ORA/combined/contrasts/.treated/de/results.tsv"))
    expect_identical(
        inv$tsv_md5[["Datasets/enrichment_ora_A_results.tsv"]],
        unname(tools::md5sum(file.path(mode_dir, "datasets",
                                       "enrichment_ora_A_results.tsv"))))
    expect_identical(names(inv$mummichog_md5), "mummichog_pinned/A/input.tsv")
    expect_setequal(inv$files, c("Datasets/enrichment_ora_A_results.tsv",
                                 "Datasets/final_results.tsv",
                                 "Enrichment/enrichment_ora_dotplot_A.png",
                                 "Enrichment/ORA/combined/contrasts/.treated/de/results.tsv",
                                 "mummichog_pinned/A/input.tsv"))
})

test_that("a before run with a lower-case datasets folder still compares", {
    res <- compare_pair(snapshot_pair(datasets = "datasets"))
    expect_true(res$all)
})

# ---- no vacuous pass on the table checks ----------------------------------------

test_that("the table checks fail when neither snapshot collected a table", {
    pair <- snapshot_pair()
    for (side in c("before", "after")) {
        pair[[side]]$tsv_md5 <- stats::setNames(character(0), character(0))
        pair[[side]]$files <- pair[[side]]$files[!grepl("\\.tsv$", pair[[side]]$files)]
    }
    res <- compare_pair(pair)

    expect_false(res$checks[["tsv_content"]])
    expect_false(res$checks[["tsv_complete"]])
    expect_false(res$all)
    expect_true(any(grepl("Take both snapshots again", res$printed)))
    # The checks that do not depend on the tables are unaffected. (The plot-file
    # check reads the contrast folders off the table paths, so it is not one.)
    expect_true(all(res$checks[c("results", "plots", "fields", "manifest")]))
})

test_that("the table checks fail when only one snapshot collected tables", {
    for (side in c("before", "after")) {
        pair <- snapshot_pair()
        pair[[side]]$tsv_md5 <- stats::setNames(character(0), character(0))
        res <- compare_pair(pair)
        expect_false(res$checks[["tsv_content"]], info = side)
        expect_false(res$checks[["tsv_complete"]], info = side)
    }
})

test_that("a snapshot that hashed only some of the tables it lists is refused", {
    # The shape of a snapshot taken by the case-sensitive version: the tables
    # are in the file list, under the lower-case folder, but were not hashed.
    pair <- snapshot_pair(datasets = "datasets")
    pair$before$tsv_md5 <- pair$before$tsv_md5[0]
    expect_false(verify_env$.el_tsv_collected(pair$before))
    expect_true(verify_env$.el_tsv_collected(pair$after))

    partly <- snapshot_pair()
    partly$before$tsv_md5 <- partly$before$tsv_md5[1]
    expect_false(verify_env$.el_tsv_collected(partly$before))
    res <- compare_pair(partly)
    expect_false(res$checks[["tsv_content"]])
    expect_false(res$checks[["tsv_complete"]])
})

test_that("the table checks still catch a real difference", {
    changed <- snapshot_pair()
    changed$after$tsv_md5[["Enrichment/QEA/combined/all_samples/results.tsv"]] <- "cccc"
    res <- compare_pair(changed)
    expect_false(res$checks[["tsv_content"]])
    expect_true(res$checks[["tsv_complete"]])

    extra <- snapshot_pair()
    extra$after$tsv_md5 <- c(extra$after$tsv_md5,
        "Enrichment/GSEA/combined/ranking_by_fc/A_vs_B/results.tsv" = "dddd")
    extra$after$files <- c(extra$after$files,
        "Enrichment/GSEA/combined/ranking_by_fc/A_vs_B/results.tsv")
    res <- compare_pair(extra)
    expect_true(res$checks[["tsv_content"]])
    expect_false(res$checks[["tsv_complete"]])
})

test_that("the plot-file check fails when the before run has no plot", {
    pair <- snapshot_pair()
    for (side in c("before", "after")) {
        pair[[side]]$files <- pair[[side]]$files[!grepl("\\.png$", pair[[side]]$files)]
    }
    expect_false(compare_pair(pair)$checks[["png"]])
})

test_that("taking a snapshot of a folder without enrichment tables is an error", {
    skip_if_not_installed("targets")
    mode_dir <- file.path(withr::local_tempdir(), "metabolomics")
    dir.create(file.path(mode_dir, "Datasets"), recursive = TRUE)
    writeLines("other", file.path(mode_dir, "Datasets", "final_results.tsv"))
    store <- withr::local_tempdir()
    # Stand in for the store read: a non-NULL result, so the call reaches the
    # inventory.
    snapshot <- verify_env$snapshot_enrichment_layout
    testthat::local_mocked_bindings(
        tar_read_raw = function(...) list(ora = NULL, files = character(0)),
        .package = "targets")

    expect_error(
        snapshot("before", mode_dir = mode_dir,
                 snap_dir = withr::local_tempdir(), store = store),
        "No enrichment result table was found")
})

# ---- plot data compared by value -------------------------------------------------

test_that("plot data equal in value compares equal where identical() does not", {
    same <- verify_env$.el_same_plots
    plot_data <- function(env) list(plots = list(dot = list(
        data = data.frame(pathway = c("P1", "P2"), FDR = c(0.01, 0.2)),
        labels = list(title = "ORA", scope = env))))
    e1 <- new.env(); assign("k", 1, envir = e1)
    e2 <- new.env(); assign("k", 1, envir = e2)
    x <- plot_data(e1); y <- plot_data(e2)

    # Two sessions never share an environment, which is all identical() sees.
    expect_false(identical(x, y))
    expect_true(same(x, y))
    expect_true(same(NULL, NULL))
})

test_that("the plot comparison is exact: no tolerance, attributes included", {
    same <- verify_env$.el_same_plots
    base <- list(dot = list(data = data.frame(FDR = c(0.01, 0.2)),
                            labels = list(title = "ORA")))

    off <- base
    off$dot$data$FDR[1] <- 0.01 + 1e-12
    # Within all.equal()'s default tolerance, so only tolerance = 0 catches it.
    expect_true(isTRUE(all.equal(base, off)))
    expect_false(same(base, off))

    attr_changed <- base
    attr(attr_changed$dot$data, "note") <- "edited"
    expect_true(isTRUE(all.equal(base, attr_changed, check.attributes = FALSE)))
    expect_false(same(base, attr_changed))

    relabelled <- base
    relabelled$dot$labels$title <- "GSEA"
    expect_false(same(base, relabelled))
    missing_plot <- base[0]
    expect_false(same(base, missing_plot))
})

test_that("the comparison reports a plot difference on check 2 only", {
    pair <- snapshot_pair()
    pair$after$plots$plots$ora_dotplot$data$FDR <- 0.01 + 1e-12
    res <- compare_pair(pair)
    expect_false(res$checks[["plots"]])
    expect_true(res$checks[["results"]])
    # The hint names the branch that differs, not its values.
    expect_true(any(grepl("differs in: plots", res$printed)))
})
