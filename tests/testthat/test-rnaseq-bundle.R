# RNA-seq omics bundle (R/domain/rnaseq/13_bundle.R): the run's finished
# objects, packaged without recomputation, minus model objects, plot paths and
# config bookkeeping.

rna_bundle_inputs <- function() {
    counts <- matrix(c(10, 20, 30, 40), nrow = 2,
                     dimnames = list(c("G1", "G2"), c("S1", "S2")))
    list(
        pre = list(expr_raw = counts, expr_filt = counts, expr_work = log2(counts + 1),
                   de_input = counts,
                   meta = data.frame(SampleID = c("S1", "S2"), group = c("A", "B")),
                   row_data = data.frame(gene_id = c("G1", "G2")),
                   info = list(filter_mode = "adaptive"),
                   contrasts = data.frame(name = "B_vs_A")),
        de_res = list(dds = structure(list(), class = "fake_dds"),
                      tables = list(B_vs_A = data.frame(FeatureID = c("G1", "G2"),
                                                         pvalue = c(0.01, 0.5))),
                      design = list(formula = "~ group")),
        pathway_res = list(annotation = data.frame(gene_id = "G1", symbol = "ABC1"),
                           pathway_results = list(B_vs_A = data.frame(pathway = "P1")),
                           plot_files = c("/home/someone/run/plot.png")),
        annot = data.frame(gene_id = "G1", symbol = "ABC1"),
        # project$dir is required: omics_bundle_info() names the run through
        # get_run_name(), which resolves the run folder under it.
        config = list(project = list(dir = tempdir(), name = "Demo", analysis_round = "R1"),
                      modes = list(rna = list(de = list(p_cutoff = 0.05))),
                      .config_path = "/home/someone/config.yaml",
                      .config_mtime = 0)
    )
}

build_parts <- function(x, ...) {
    build_rnaseq_bundle_parts(pre = x$pre, de_res = x$de_res, pathway_res = x$pathway_res,
                              annot = x$annot, config = x$config, ...)
}


test_that("the bundle keeps the run's objects and drops only model, plots and bookkeeping", {
    x <- rna_bundle_inputs()
    parts <- build_parts(x)

    # No Trinotate in this run, so no trinotate part rather than a NULL one.
    expect_named(parts, c("pre", "de", "pathway", "annotation", "config"))
    expect_identical(parts$pre, x$pre)
    expect_identical(parts$annotation, x$annot)

    expect_named(parts$de, c("tables", "design"))
    expect_identical(parts$de$tables, x$de_res$tables)

    expect_false("plot_files" %in% names(parts$pathway))
    expect_identical(parts$pathway$pathway_results, x$pathway_res$pathway_results)
    expect_identical(parts$pathway$annotation, x$pathway_res$annotation)

    expect_named(parts$config, c("project", "modes"))
    expect_identical(parts$config$modes, x$config$modes)
})

test_that("optional parts are included only when the run produced them", {
    x <- rna_bundle_inputs()
    trin <- data.frame(gene_id = "G1", go = "GO:0008150")
    with_trin <- build_parts(x, trinotate = trin)
    expect_identical(with_trin$trinotate, trin)

    bare <- build_rnaseq_bundle_parts(pre = x$pre, de_res = x$de_res, config = x$config)
    expect_named(bare, c("pre", "de", "config"))
})

test_that("objects that are not RNA results are refused", {
    x <- rna_bundle_inputs()
    expect_error(build_rnaseq_bundle_parts(pre = x$pre, de_res = list(dds = NULL),
                                           config = x$config),
                 "no DE tables")
    expect_error(build_rnaseq_bundle_parts(pre = list(meta = x$pre$meta),
                                           de_res = x$de_res, config = x$config))
})

test_that("an RNA bundle written from these parts reads back unchanged", {
    x <- rna_bundle_inputs()
    root <- withr::local_tempdir()
    prov <- list(created_at_utc = "2026-01-31T09:30:00Z", git_commit = "abc",
                 git_dirty = FALSE, r_version = "4.5.1", renv_lock_sha256 = "x",
                 packages = list(base = "4.5.1"))

    path <- write_omics_bundle(build_parts(x), omics_bundle_info(x$config, "rnaseq"),
                               root, prov)

    # Normalised: on Windows the temp dir can carry backslashes that dirname()
    # turns into forward slashes.
    norm <- function(p) normalizePath(p, winslash = "/", mustWork = TRUE)
    expect_identical(norm(dirname(path)), norm(file.path(root, "bundles", "rnaseq")))
    b <- read_omics_bundle(path)
    expect_identical(b$parts, build_parts(x))
    expect_identical(b$manifest$project, "Demo")
    expect_identical(b$manifest$analysis_round, "R1")
    expect_identical(b$manifest$info$run_name, "Demo_R1")
})

test_that("bundle identity matches the run folder, defaults and numeric rounds included", {
    configs <- list(
        list(project = list(dir = "/p", name = "Demo", analysis_round = "R1")),
        list(project = list(dir = "/p", name = "Demo", analysis_round = 2)),
        list(project = list(dir = "/p"))
    )
    for (cfg in configs) {
        info <- omics_bundle_info(cfg, "rnaseq")
        expect_identical(paste(info$project, info$analysis_round, sep = "_"),
                         get_run_name(cfg))
        expect_identical(info$run_name, get_run_name(cfg))
    }
    expect_identical(omics_bundle_info(configs[[3]], "rnaseq")$project, "Project")
    expect_identical(omics_bundle_info(configs[[2]], "rnaseq")$analysis_round, "2")
})
