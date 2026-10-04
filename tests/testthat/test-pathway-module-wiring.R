# Regression tests for arg-drift between pathway modules and run_pathway_analysis().
# Each module runs end-to-end with a tiny synthetic DE table + an inline GMT,
# and asserts no error — catches "unused argument" and similar wiring breakage
# that CI's unit tests don't otherwise exercise (CI doesn't invoke tar_make()).
# The proteomics tests also assert which identifiers reach enrichment.

test_that("mod_rnaseq_pathway: end-to-end wiring with synthetic DE + tiny GMT", {
    skip_if_not_installed("fgsea")

    gmt <- tempfile(fileext = ".gmt")
    writeLines(c(
        "SET_A\tdesc A\tg1\tg2\tg3\tg4\tg5",
        "SET_B\tdesc B\tg6\tg7\tg8\tg9\tg10\tg11\tg12"
    ), gmt)
    on.exit(unlink(gmt), add = TRUE)

    de_res <- list(tables = list(
        contrast_A = data.frame(
            FeatureID      = paste0("g", 1:15),
            log2FoldChange = c(2.0, -1.5, 0.1, 1.2, -2.1, 0.5, 0.3, -0.7, 1.8, -1.1,
                               0.4, -0.2, 1.5, -1.8, 0.6),
            pvalue         = c(0.001, 0.01, 0.8, 0.04, 0.005, 0.2, 0.3, 0.15, 0.002, 0.03,
                               0.5, 0.7, 0.008, 0.001, 0.4),
            padj           = c(0.01, 0.05, 0.9, 0.1, 0.02, 0.3, 0.4, 0.25, 0.01, 0.08,
                               0.6, 0.8, 0.04, 0.01, 0.5),
            stat           = c(5, -3, 0.2, 2.5, -4, 1, 1.2, -1.5, 4.5, -2.8,
                               0.8, -0.3, 3.2, -4.5, 1.1),
            stringsAsFactors = FALSE
        )
    ))
    pre    <- list(expr_filt = NULL, meta = NULL)
    config <- list(modes = list(rna = list(
        annotation = list(skip_annotation = TRUE, organism = "Homo sapiens"),
        pathway    = list(
            enabled   = TRUE,
            databases = character(0),
            gmt_file  = gmt,
            min_size  = 3,
            max_size  = 50
        )
    )))
    out_dir <- tempfile("test-rna-pathway-")
    dir.create(out_dir)
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    expect_no_error(suppressMessages(suppressWarnings(
        mod_rnaseq_pathway(de_res, pre, config, out_dir)
    )))
})

test_that("mod_rnaseq_pathway resolves a list of relative GMT paths against the raw dir", {
    skip_if_not_installed("fgsea")

    # Two GMTs given relative to paths$raw, as a YAML list arrives. Both broke:
    # a scalar path check errors on a list, and an unresolved relative path
    # loads zero gene sets.
    proj <- tempfile("test-rna-gmt-proj-")
    dir.create(file.path(proj, "data", "gmt"), recursive = TRUE)
    on.exit(unlink(proj, recursive = TRUE), add = TRUE)
    writeLines("SET_A\tdesc A\tg1\tg2\tg3\tg4\tg5",
               file.path(proj, "data", "gmt", "a.gmt"))
    writeLines("SET_B\tdesc B\tg6\tg7\tg8\tg9\tg10\tg11\tg12",
               file.path(proj, "data", "gmt", "b.gmt"))

    de_res <- list(tables = list(
        contrast_A = data.frame(
            FeatureID      = paste0("g", 1:15),
            log2FoldChange = c(2.0, -1.5, 0.1, 1.2, -2.1, 0.5, 0.3, -0.7, 1.8, -1.1,
                               0.4, -0.2, 1.5, -1.8, 0.6),
            pvalue         = c(0.001, 0.01, 0.8, 0.04, 0.005, 0.2, 0.3, 0.15, 0.002, 0.03,
                               0.5, 0.7, 0.008, 0.001, 0.4),
            padj           = c(0.01, 0.05, 0.9, 0.1, 0.02, 0.3, 0.4, 0.25, 0.01, 0.08,
                               0.6, 0.8, 0.04, 0.01, 0.5),
            stat           = c(5, -3, 0.2, 2.5, -4, 1, 1.2, -1.5, 4.5, -2.8,
                               0.8, -0.3, 3.2, -4.5, 1.1),
            stringsAsFactors = FALSE
        )
    ))
    config <- list(
        project = list(dir = proj),
        paths   = list(raw = "data"),
        modes   = list(rna = list(
            annotation = list(skip_annotation = TRUE, organism = "Homo sapiens"),
            pathway    = list(
                enabled   = TRUE,
                databases = character(0),
                gmt_file  = list("gmt/a.gmt", "gmt/b.gmt"),
                min_size  = 3,
                max_size  = 50
            )
        ))
    )
    out_dir <- tempfile("test-rna-pathway-rel-")
    dir.create(out_dir)
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    res <- suppressMessages(suppressWarnings(
        mod_rnaseq_pathway(de_res, list(expr_filt = NULL, meta = NULL), config, out_dir)
    ))
    # With zero gene sets loaded the module returns an empty pathway_results.
    expect_gt(length(res$pathway_results), 0)
})

# run_proteomics_pathway() maps FeatureID to the first Genes symbol before
# enrichment, except when annotation.skip_annotation is true (#148): a
# non-model run's custom GMT is keyed on the raw Protein.Group FeatureID, and
# remapping there left every gene set with an empty intersection, silently.
# The tests below assert on the identifiers that reach fgsea, read back from
# the leading edges it returns, so reverting that gate fails them.

prot_pathway_fixture <- function(skip_annotation, gmt_lines) {
    gmt <- tempfile(fileext = ".gmt")
    writeLines(gmt_lines, gmt)
    summary_df <- data.frame(
        FeatureID                  = paste0("p", 1:15),
        Genes                      = paste0("g", 1:15),
        padj.imputs.contrast_A     = c(0.01, 0.05, 0.9, 0.1, 0.02, 0.3, 0.4, 0.25, 0.01, 0.08,
                                       0.6, 0.8, 0.04, 0.01, 0.5),
        pvalue.imputs.contrast_A   = c(0.001, 0.01, 0.8, 0.04, 0.005, 0.2, 0.3, 0.15, 0.002, 0.03,
                                       0.5, 0.7, 0.008, 0.001, 0.4),
        linearFC.imputs.contrast_A = c(4.0, -2.8, 1.1, 2.3, -4.2, 1.4, 1.2, -1.6, 3.5, -2.1,
                                       1.3, -1.1, 2.8, -3.5, 1.5),
        stringsAsFactors = FALSE,
        check.names      = FALSE
    )
    config <- list(modes = list(proteomics = list(
        de_table   = list(id_col = "FeatureID"),
        annotation = list(skip_annotation = skip_annotation, organism = "Homo sapiens"),
        pathway    = list(
            enabled   = TRUE,
            databases = character(0),
            gmt_file  = gmt,
            min_size  = 3,
            max_size  = 50
        )
    )))
    list(gmt = gmt, de_res = list(summary_df = summary_df), config = config)
}

# Every identifier in the leading edges of the fgsea tables for one contrast.
prot_fgsea_ids <- function(res, contrast = "contrast_A") {
    tabs <- res[[contrast]][grepl("_fgsea$", names(res[[contrast]]))]
    unlist(lapply(tabs, function(t) strsplit(t$leadingEdge, ",", fixed = TRUE)),
           use.names = FALSE)
}

test_that("run_proteomics_pathway keeps FeatureID when skip_annotation is true", {
    skip_if_not_installed("fgsea")

    # A custom GMT keyed on the raw FeatureID, as a non-model run has it.
    fx <- prot_pathway_fixture(TRUE, c(
        "SET_A\tdesc A\tp1\tp2\tp3\tp4\tp5",
        "SET_B\tdesc B\tp6\tp7\tp8\tp9\tp10\tp11\tp12"
    ))
    on.exit(unlink(fx$gmt), add = TRUE)
    out_dir <- withr::local_tempdir("test-prot-pathway-")

    res <- suppressMessages(suppressWarnings(
        run_proteomics_pathway(fx$de_res, list(), fx$config, out_dir)
    ))

    # Remapped to g1..g15, nothing would overlap this GMT and no fgsea table
    # would be produced at all.
    ids <- prot_fgsea_ids(res)
    expect_gt(length(ids), 0)
    expect_true(all(ids %in% paste0("p", 1:15)))
})

test_that("run_proteomics_pathway maps FeatureID to gene symbols when skip_annotation is false", {
    skip_if_not_installed("fgsea")

    fx <- prot_pathway_fixture(FALSE, c(
        "SET_A\tdesc A\tg1\tg2\tg3\tg4\tg5",
        "SET_B\tdesc B\tg6\tg7\tg8\tg9\tg10\tg11\tg12"
    ))
    on.exit(unlink(fx$gmt), add = TRUE)
    out_dir <- withr::local_tempdir("test-prot-pathway-")

    # With skip_annotation false the function also annotates the ids through
    # biomaRt/OrgDb. That lookup is not what this test is about and would reach
    # the network, so it returns no annotation here; the symbol mapping comes
    # from the Genes column regardless. Sourced functions live in the global
    # environment, not a package, so stub by assignment and restore after.
    target <- environment(run_proteomics_pathway)
    real_annotate <- get("annotate_genes_v2", envir = target)
    assign("annotate_genes_v2", function(...) NULL, envir = target)
    withr::defer(assign("annotate_genes_v2", real_annotate, envir = target))

    res <- suppressMessages(suppressWarnings(
        run_proteomics_pathway(fx$de_res, list(), fx$config, out_dir)
    ))

    ids <- prot_fgsea_ids(res)
    expect_gt(length(ids), 0)
    expect_true(all(ids %in% paste0("g", 1:15)))
})
