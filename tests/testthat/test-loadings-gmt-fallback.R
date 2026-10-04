# Loadings enrichment for gene views on organisms with no KEGG code or OrgDb.
#
# enrich_feature_list() returned NULL for those, which is every non-model
# organism -- so a run produced loadings enrichment for metabolomics only, and
# the gene views were simply empty. Nothing said that the gene sets to test
# against were configured and sitting unused. The fallback runs the same
# over-representation test on those GMTs.
#
# Two things have to hold for it to mean anything:
#
#   The query feature list must survive ID resolution. MOFA2 makes feature names
#   unique across views by appending the view name, so the resolver has to
#   accept both spellings -- otherwise the fallback runs on a partial list and
#   under-reports without saying so.
#
#   The background must be the view's own features in the query's namespace. Let
#   enricher() default it and the gene sets become the universe, which inflates
#   every p-value.
#
# All fixtures synthetic; no network, no clusterProfiler call.

# Stub by assignment into the function's own environment. with_mocked_bindings()
# errors with "No packages loaded with pkgload" here, since these are sourced
# functions rather than a package. Same shape as the one in
# test-gmt-kegg-class-exclusion.R, which is file-local there too.
local_stubs <- function(stubs, env = parent.frame()) {
    target <- environment(enrich_feature_list_gmt)
    nms <- names(stubs)

    had <- vapply(nms, exists, logical(1), envir = target, inherits = FALSE)
    old <- lapply(nms[had], get, envir = target, inherits = FALSE)
    names(old) <- nms[had]

    withr::defer({
        for (nm in nms) {
            if (nm %in% names(old)) {
                assign(nm, old[[nm]], envir = target)
            } else if (exists(nm, envir = target, inherits = FALSE)) {
                rm(list = nm, envir = target)
            }
        }
    }, envir = env)

    for (nm in nms) assign(nm, stubs[[nm]], envir = target)
    invisible(NULL)
}

harm_fixture <- function(n = 4, features = NULL) {
    gpm <- data.frame(
        gene_id    = paste0("WBGene", sprintf("%05d", seq_len(n))),
        protein_id = paste0("P", sprintf("%05d", seq_len(n))),
        stringsAsFactors = FALSE
    )
    if (is.null(features)) features <- paste0("GENE_", seq_len(n))
    list(
        gene_protein_mapping = gpm,
        inputs = list(transcriptomics = list(
            expr_work = matrix(0, nrow = length(features), ncol = 2,
                               dimnames = list(features, c("s1", "s2")))
        ))
    )
}


# ---- resolving the IDs the query is built from ------------------------------

test_that("a MOFA view-suffixed GENE_N resolves like the bare one", {
    # The regression that makes this feature worth having: MOFA2 requires unique
    # feature names across views and appends the view name to shared ones, so
    # "GENE_2" and "GENE_2_transcriptomics" are the same feature.
    harm <- harm_fixture()

    bare <- suppressMessages(
        resolve_gene_n_ids("GENE_2", harm, "transcriptomics"))
    suffixed <- suppressMessages(
        resolve_gene_n_ids("GENE_2_transcriptomics", harm, "transcriptomics"))

    expect_identical(bare, "WBGene00002")
    expect_identical(suffixed, bare)
})

test_that("a mixed list resolves both spellings and keeps its length", {
    harm <- harm_fixture()

    resolved <- suppressMessages(resolve_gene_n_ids(
        c("GENE_1", "GENE_2_transcriptomics", "GENE_3_proteomics"),
        harm, "transcriptomics"))

    expect_identical(resolved,
                     c("WBGene00001", "WBGene00002", "WBGene00003"))
})

test_that("a native id passes through and an unmapped GENE_N is dropped", {
    # Dropped rather than carried as NA: an NA would reach the query list and be
    # counted as a feature that simply failed to enrich.
    harm <- harm_fixture(n = 2)

    expect_identical(
        suppressMessages(resolve_gene_n_ids(c("WBGene00001", "already_native"),
                                            harm, "transcriptomics")),
        c("WBGene00001", "already_native"))
    expect_identical(
        suppressMessages(resolve_gene_n_ids(c("GENE_1", "GENE_99"),
                                            harm, "transcriptomics")),
        "WBGene00001")
})

test_that("metabolomics ids are returned untouched", {
    # Metabolomics uses feature_N, not GENE_N, and has no row in the mapping.
    harm <- harm_fixture()

    expect_identical(
        resolve_gene_n_ids(c("feature_1", "feature_2"), harm, "metabolomics"),
        c("feature_1", "feature_2"))
})


# ---- when the fallback declines ---------------------------------------------

gmt_cfg <- function(dir, file = "sets.gmt", mode = "rna") {
    cfg <- list(project = list(dir = dir), paths = list(raw = "."),
                modes = list())
    cfg$modes[[mode]] <- list(pathway = list(gmt_file = file))
    cfg
}

write_gmt_fixture <- function(dir, file = "sets.gmt") {
    lines <- c(
        paste("SET_A", "first set", "WBGene00001", "WBGene00002", sep = "\t"),
        paste("SET_B", "second set", "WBGene00003", "WBGene00004", sep = "\t")
    )
    path <- file.path(dir, file)
    writeLines(lines, path)
    path
}

test_that("no config means NULL, which is what this path did before", {
    expect_null(enrich_feature_list_gmt("WBGene00001", "transcriptomics",
                                        harm_fixture(), NULL))
})

test_that("an omics type with no configured GMT declines", {
    # Only the two gene views carry one; metabolomics goes through the compound
    # branch and must not be routed here.
    dir <- withr::local_tempdir()
    expect_null(enrich_feature_list_gmt("cpd1", "metabolomics",
                                        harm_fixture(), gmt_cfg(dir)))
})

test_that("an unset or missing gmt_file declines rather than erroring", {
    dir <- withr::local_tempdir()
    harm <- harm_fixture()

    unset <- gmt_cfg(dir)
    unset$modes$rna$pathway$gmt_file <- ""
    expect_null(enrich_feature_list_gmt("WBGene00001", "transcriptomics",
                                        harm, unset))

    # Configured but not on disk: says which file, and returns NULL.
    expect_message(
        res <- enrich_feature_list_gmt("WBGene00001", "transcriptomics",
                                       harm, gmt_cfg(dir, "absent.gmt")),
        "gmt_file not found")
    expect_null(res)
})


# ---- the test it actually runs ----------------------------------------------

test_that("the fallback runs ORA on the loadings list and returns the KEGG shape", {
    # Stubbed at run_multi_ora_enricher(): the point here is what this function
    # hands it and what it makes of the answer, not clusterProfiler's maths.
    dir <- withr::local_tempdir()
    write_gmt_fixture(dir)
    captured <- NULL

    local_stubs(list(run_multi_ora_enricher = function(sig_genes, universe,
                                                           term2gene, term2name = NULL,
                                                           ...) {
        captured <<- list(sig_genes = sig_genes, universe = universe,
                          term2gene = term2gene)
        data.frame(pathway = "first set", ID = "SET_A", pvalue = 0.001,
                   padj = 0.01, GeneRatio = "2/2", Count = 2L,
                   stringsAsFactors = FALSE)
    }))

    res <- suppressMessages(enrich_feature_list_gmt(
        c("WBGene00001", "WBGene00002"), "transcriptomics",
        harm_fixture(), gmt_cfg(dir)))

    # The column shape .rbind_fill() stacks the compound-view results onto.
    expect_identical(names(res),
                     c("pathway", "ID", "pvalue", "padj", "GeneRatio", "setSize"))
    expect_equal(res$setSize, 2L)
    expect_equal(res$padj, 0.01)

    # The gene sets came from the GMT, not from anywhere else.
    expect_setequal(unique(captured$term2gene$term), c("SET_A", "SET_B"))
})

test_that("the background is the view's own features, not the gene sets", {
    # Left NULL, enricher() uses the gene sets as the universe and every p-value
    # is inflated. The background must also be in the query's namespace, so it
    # goes through the same resolution -- including the MOFA suffix.
    dir <- withr::local_tempdir()
    write_gmt_fixture(dir)
    captured <- NULL

    local_stubs(list(run_multi_ora_enricher = function(sig_genes, universe,
                                                           term2gene, term2name = NULL,
                                                           ...) {
        captured <<- list(sig_genes = sig_genes, universe = universe)
        NULL
    }))

    harm <- harm_fixture(
        n = 4,
        features = c("GENE_1_transcriptomics", "GENE_2_transcriptomics",
                     "GENE_3", "GENE_4"))

    suppressMessages(enrich_feature_list_gmt(
        c("WBGene00001", "WBGene00002"), "transcriptomics", harm, gmt_cfg(dir)))

    expect_setequal(captured$universe,
                    c("WBGene00001", "WBGene00002", "WBGene00003", "WBGene00004"))
    expect_setequal(captured$sig_genes, c("WBGene00001", "WBGene00002"))
})

test_that("no usable background declines instead of testing against the gene sets", {
    # The failure this guards is silent: with universe = NULL, enricher() takes
    # its background from TERM2GENE, so the universe becomes the union of the
    # gene sets rather than what the view measured, every p-value inflates, and
    # the result looks like an ordinary enrichment table. The enricher must not
    # be reached at all.
    dir <- withr::local_tempdir()
    write_gmt_fixture(dir)
    called <- FALSE

    local_stubs(list(run_multi_ora_enricher = function(...) {
        called <<- TRUE
        NULL
    }))

    harm <- harm_fixture()

    # No preprocessed data for this view at all.
    no_inputs <- harm
    no_inputs$inputs <- list()
    expect_message(
        res <- enrich_feature_list_gmt("WBGene00001", "transcriptomics",
                                       no_inputs, gmt_cfg(dir)),
        "no usable background")
    expect_null(res)

    # Present, but carrying no expr_work.
    no_expr <- harm
    no_expr$inputs$transcriptomics <- list()
    expect_null(suppressMessages(enrich_feature_list_gmt(
        "WBGene00001", "transcriptomics", no_expr, gmt_cfg(dir))))

    # Present with expr_work, but its rows are unnamed -- resolution yields
    # nothing, which reaches run_multi_ora_enricher() as NULL just the same.
    unnamed <- harm
    unnamed$inputs$transcriptomics$expr_work <- matrix(0, nrow = 2, ncol = 2)
    expect_null(suppressMessages(enrich_feature_list_gmt(
        "WBGene00001", "transcriptomics", unnamed, gmt_cfg(dir))))

    expect_false(called)
})

test_that("an enricher that finds nothing yields NULL, not an empty frame", {
    dir <- withr::local_tempdir()
    write_gmt_fixture(dir)

    local_stubs(list(run_multi_ora_enricher = function(...) NULL))
    expect_null(suppressMessages(enrich_feature_list_gmt(
        "WBGene00001", "transcriptomics", harm_fixture(), gmt_cfg(dir))))

    local_stubs(list(run_multi_ora_enricher = function(...) {
        data.frame(pathway = character(0), ID = character(0),
                   pvalue = numeric(0), padj = numeric(0),
                   GeneRatio = character(0), Count = integer(0),
                   stringsAsFactors = FALSE)
    }))
    expect_null(suppressMessages(enrich_feature_list_gmt(
        "WBGene00001", "transcriptomics", harm_fixture(), gmt_cfg(dir))))
})


# ---- reaching the fallback from enrich_feature_list() -----------------------

test_that("enrich_feature_list routes to the GMT branch with no KEGG or OrgDb", {
    # The branch this whole change exists for. Pinned at the source: reaching it
    # through the real function needs an OrgDb to be absent, which is exactly
    # the case the KEGG path cannot be exercised in.
    body_src <- paste(deparse(body(enrich_feature_list)), collapse = " ")

    expect_true(grepl("enrich_feature_list_gmt", body_src, fixed = TRUE))
    # Resolution has to precede the guard, or the fallback gets raw GENE_N ids.
    expect_lt(regexpr("resolve_gene_n_ids", body_src, fixed = TRUE),
              regexpr("enrich_feature_list_gmt", body_src, fixed = TRUE))
})

test_that("both loadings callers can supply the config the fallback needs", {
    # config reaches enrich_feature_list() only if it is threaded through the
    # DIABLO and MOFA wrappers; without it the fallback declines silently and
    # the feature looks like it does nothing.
    for (fn in list(run_diablo_loadings_enrichment, run_mofa_weights_enrichment)) {
        expect_true("config" %in% names(formals(fn)))
    }
    expect_true("config" %in% names(formals(enrich_feature_list)))

    for (fn in list(run_diablo_loadings_enrichment, run_mofa_weights_enrichment)) {
        body_src <- paste(deparse(body(fn)), collapse = " ")
        expect_true(grepl("config = config", body_src, fixed = TRUE))
    }
})
