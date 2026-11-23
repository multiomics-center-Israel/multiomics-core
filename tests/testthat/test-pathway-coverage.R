# KEGG pathway coverage per layer (07h_pathway_coverage.R): why a cell in a
# KEGG figure is empty, and that labelling it changes no result.
#
# A layer that could not test a pathway -- no KEGG members of its omics type,
# none or too few measured, no ids mapped -- is "not testable" and drawn grey
# with a cross; a layer that could test it and reported nothing is "neutral"
# and drawn white. Nothing here touches the network: compound tables are
# injected or read from a temporary cache, and the gene-pathway getter is only
# called with download = FALSE.
#
# All fixtures are synthetic.

cov_cpd_fixture <- function() {
    data.frame(
        pathway  = c(rep("map00010", 5), rep("map00020", 3), rep("map00030", 2),
                     rep("map00040", 4)),
        compound = c(sprintf("C%05d", 1:5), sprintf("C%05d", 6:8),
                     sprintf("C%05d", 9:10), sprintf("C%05d", 11:14)),
        name     = c(rep("Glycolysis", 5), rep("Citrate cycle", 3),
                     rep("Pentose phosphate", 2), rep("Fructose metabolism", 4)),
        stringsAsFactors = FALSE
    )
}

# Measured: 4 of map00010, 3 of map00020, 1 of map00030, none of map00040.
cov_measured <- c(sprintf("C%05d", 1:4), sprintf("C%05d", 6:8), "C00009")

cov_fixture <- function(kegg_org = "rno") {
    cpd <- compound_pathway_coverage(cov_cpd_fixture(), cov_measured,
                                     compound_min_measured(10), kegg_org = kegg_org)
    list(layers = list(metabolomics = cpd),
         summary = data.frame(layer = "metabolomics", status = "mapped",
                              stringsAsFactors = FALSE))
}


# ---- one floor for every compound producer ---------------------------------

test_that("the compound floor is three, whatever the gene-set minimum", {
    expect_identical(compound_min_measured(10), 3L)
    expect_identical(compound_min_measured(2), 3L)
    expect_identical(compound_min_measured(), 3L)
})

test_that("compound ORA, compound GSEA and Multi-ORA all take the shared floor", {
    src <- function(f) paste(deparse(body(f)), collapse = " ")
    expect_true(grepl("compound_min_measured(min_gs)", src(run_compound_ora), fixed = TRUE))
    expect_true(grepl("compound_min_measured(min_gs)", src(run_compound_gsea), fixed = TRUE))
    for (f in list(run_multi_ora, .run_multi_ora_contrast_group)) {
        s <- src(f)
        expect_true(grepl("run_compound_ora(de_mapped, out_dir, compound_min_measured()",
                          s, fixed = TRUE))
        # The literal floor of two is gone.
        expect_false(grepl("run_compound_ora(de_mapped, out_dir, 2,", s, fixed = TRUE))
    }
    # Multi-ORA coverage uses the same floor for metabolomics.
    expect_identical(.multi_ora_min_measured()$metabolomics, compound_min_measured())
})

test_that("a pathway with exactly two measured compounds is not tested by compound ORA", {
    cache <- withr::local_tempdir()
    saveRDS(cov_cpd_fixture(), file.path(cache, "kegg_compound_pathways.rds"))
    # map00010: 4 measured, 2 significant; map00030: 2 measured, both significant.
    universe <- c(sprintf("C%05d", 1:4), "C00009", "C00010", sprintf("C%05d", 6:8))
    de <- data.frame(KEGG_ID = universe,
                     padj = c(0.01, 0.01, 0.9, 0.9, 0.01, 0.01, 0.9, 0.9, 0.9),
                     pvalue = c(1e-4, 1e-4, 0.5, 0.5, 1e-4, 1e-4, 0.5, 0.5, 0.5),
                     stringsAsFactors = FALSE)
    res <- suppressMessages(run_compound_ora(de, cache, min_gs = 2, max_gs = 500,
                                             pval_cutoff = 1, universe = universe))
    expect_false("map00030" %in% res$ID)
    expect_true("map00010" %in% res$ID)
})

test_that("the gene floor is the one the KEGG-from-DE producer applies", {
    expect_identical(gene_min_measured(10, 200), 5L)
    expect_identical(gene_min_measured(2, 200), 3L)
    expect_identical(gene_min_measured(10, 800), 10L)
    s <- paste(deparse(body(run_kegg_enrichment_for_omics)), collapse = " ")
    expect_true(grepl("gene_min_measured(min_gs, nrow(de_mapped))", s, fixed = TRUE))
})


# ---- classes ---------------------------------------------------------------

test_that("each pathway is classed from its measured compounds", {
    cov <- cov_fixture()$layers$metabolomics
    cls <- setNames(cov$coverage, cov$pathway_id)
    expect_identical(unname(cls["map00010"]), "testable")
    expect_identical(unname(cls["map00020"]), "testable")
    expect_identical(unname(cls["map00030"]), "below_min_measured")
    expect_identical(unname(cls["map00040"]), "no_measured_members")
    expect_identical(cov$n_in_pathway[cov$pathway_id == "map00010"], 5L)
    expect_identical(cov$n_measured[cov$pathway_id == "map00010"], 4L)
    expect_true(all(cov$n_min_required == 3L))
})

test_that("a maximum set size is honoured", {
    cov <- compound_pathway_coverage(cov_cpd_fixture(), cov_measured, 3,
                                     max_measured = 3)
    expect_identical(cov$coverage[cov$pathway_id == "map00010"], "above_max_measured")
})

test_that("a KEGG pathway with no compounds is no_kegg_members; GO is left alone", {
    tab <- cov_fixture()$layers$metabolomics
    got <- layer_coverage_for_keys(c("00010", "04110", "GO:0006915"), tab, "rno")
    expect_identical(got$coverage, c("testable", "no_kegg_members", NA))
    expect_identical(got$n_in_pathway[2], 0L)
    expect_true(is.na(got$n_measured[3]))
})

test_that("organism, map and bare spellings join to the same pathway", {
    tab <- cov_fixture()$layers$metabolomics
    got <- layer_coverage_for_keys(c("rno00030", "map00030", "00030"), tab, "rno")
    expect_identical(got$coverage, rep("below_min_measured", 3))
})

test_that("nothing mapped classes every pathway no_id_mapping", {
    cov <- compound_pathway_coverage(cov_cpd_fixture(), character(0), 3)
    done <- .finish_layer_coverage("metabolomics", cov, feature_ids = c("a", "b"),
                                   n_mapped = 0L, measured = character(0), min_req = 3L)
    expect_true(all(done$pathways$coverage == "no_id_mapping"))
    expect_identical(done$summary$status, "no_id_mapping")
    # A KEGG key absent from the table is unmapped too, not "no KEGG members".
    got <- layer_coverage_for_keys("04110", done$pathways)
    expect_identical(got$coverage, "no_id_mapping")
})

test_that("annotation adds per-layer columns and leaves a table without coverage alone", {
    df <- data.frame(norm_id = c("00010", "00040", "GO:0006915"), x = 1:3,
                     stringsAsFactors = FALSE)
    expect_identical(annotate_pathway_coverage(df, NULL), df)
    out <- annotate_pathway_coverage(df, cov_fixture(), kegg_org = "rno")
    expect_true(all(c("n_in_pathway_metabolomics", "n_measured_metabolomics",
                      "n_min_required_metabolomics", "coverage_metabolomics")
                    %in% names(out)))
    expect_identical(out$coverage_metabolomics,
                     c("testable", "no_measured_members", NA))
    # Only the layers asked for.
    expect_identical(annotate_pathway_coverage(df, cov_fixture(), layers = "proteomics"), df)
})

test_that("reclassifying keeps the counts and moves only the class", {
    cov <- reclassify_pathway_coverage(cov_fixture(), list(metabolomics = 5L))
    tab <- cov$layers$metabolomics
    expect_identical(tab$coverage[tab$pathway_id == "map00010"], "below_min_measured")
    expect_identical(tab$n_measured, cov_fixture()$layers$metabolomics$n_measured)
    expect_true(all(tab$n_min_required == 5L))
})


# ---- display state ---------------------------------------------------------

test_that("cells are scored, neutral or not testable, one rule for every class", {
    st <- pathway_cell_state(
        c(2, NA, NA, NA, NA, NA, NA),
        c("no_kegg_members", "testable", "no_kegg_members", "no_measured_members",
          "below_min_measured", "no_id_mapping", NA))
    expect_identical(st, c("scored", "neutral", rep("not_testable", 4), NA))
})

test_that("the state matrix covers only the layers with coverage", {
    m <- matrix(c(1, NA, NA, 2, NA, NA), nrow = 3,
                dimnames = list(c("00010", "00030", "00040"),
                                c("proteomics", "metabolomics")))
    st <- .coverage_state_matrix(m, list(metabolomics = c("testable", "testable",
                                                         "no_kegg_members")))
    expect_true(all(is.na(st[, "proteomics"])))
    expect_identical(unname(st[, "metabolomics"]), c("scored", "neutral", "not_testable"))
    expect_null(.coverage_state_matrix(m, list()))
})

test_that("without states the heatmap scale is exactly the old one", {
    pal <- grDevices::colorRampPalette(c("white", "red"))(50)
    m <- matrix(c(1, 3), nrow = 1)
    sc <- .coverage_heatmap_scale(m, NULL, pal, c(1, 3))
    expect_identical(sc$colors, pal)
    expect_equal(sc$breaks, seq(1, 3, length.out = 51))
    expect_identical(sc$marks, FALSE)
    expect_identical(sc$values, m)
    expect_true(is.na(sc$legend_breaks))
})

test_that("neutral cells are white, real values keep their colour, crosses only where not testable", {
    pal <- grDevices::colorRampPalette(c("#fff7f3", "#7a0177"))(50)
    zl <- c(1, 4)
    m <- matrix(c(1, NA, NA, 4), nrow = 2,
                dimnames = list(c("a", "b"), c("proteomics", "metabolomics")))
    st <- matrix(c(NA, NA, "neutral", "scored"), nrow = 2, dimnames = dimnames(m))
    st["b", "proteomics"] <- "not_testable"
    sc <- .coverage_heatmap_scale(m, st, pal, zl)

    bin <- function(v) .bincode(v, sc$breaks, right = TRUE, include.lowest = TRUE)
    # The scale's lowest real value keeps the palette's first colour ...
    expect_identical(sc$colors[bin(zl[1])], pal[1])
    expect_identical(sc$colors[bin(zl[2])], pal[50])
    # ... and the neutral cell falls in the white interval below it.
    expect_identical(sc$colors[bin(sc$values["a", "metabolomics"])], "white")
    # The not-testable cell stays NA (the figure's grey) and carries the cross.
    expect_true(is.na(sc$values["b", "proteomics"]))
    expect_identical(unname(sc$marks["b", "proteomics"]), "×")
    expect_identical(sum(nzchar(sc$marks)), 1L)
    # The palette itself is still all there, in order, after the two added bins.
    expect_identical(sc$colors[-(1:2)], pal)
})

test_that("the ORA heatmap returns the coverage states it drew", {
    m <- matrix(c(1e-3, NA, NA, NA, 1e-4, NA), nrow = 3,
                dimnames = list(c("00010", "00030", "04110"),
                                c("proteomics", "metabolomics")))
    m["04110", "proteomics"] <- 1e-2
    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)
    drawn <- plot_cross_omics_ora_heatmap(m, pathway_coverage = cov_fixture(),
                                          kegg_org = "rno")
    st <- attr(drawn, "cell_state")
    expect_identical(unname(st[, "metabolomics"][order(rownames(drawn))]),
                     c("neutral", "scored", "not_testable"))
    expect_true(all(is.na(st[, "proteomics"])))
    # And without coverage nothing is attached.
    plain <- plot_cross_omics_ora_heatmap(m)
    expect_null(attr(plain, "cell_state"))
})

test_that("the meta-analysis heatmap draws with coverage columns", {
    meta <- data.frame(pathway = c("Glycolysis", "Pentose", "Cell cycle"),
                       norm_id = c("00010", "00030", "04110"),
                       pval_proteomics = c(1e-3, 1e-2, 1e-4),
                       pval_metabolomics = c(1e-2, NA, NA),
                       n_omics = c(2L, 1L, 1L), combined_pval = c(1e-4, 1e-2, 1e-4),
                       stringsAsFactors = FALSE)
    meta <- annotate_pathway_coverage(meta, cov_fixture(), kegg_org = "rno")
    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_no_error(plot_cross_omics_pathway_heatmap(meta, c("proteomics",
                                                             "metabolomics")))
})


# ---- lookup ----------------------------------------------------------------

test_that("the lookup names why the target has no result", {
    tabs <- list(
        proteomics = data.frame(pathway = c("rno00010", "rno00030", "rno04110", "rno00040"),
                                pval = c(0.01, 0.02, 0.03, 0.04), NES = c(2, 1.5, 1.2, 1),
                                size = 20, method = "fgsea", contrast = "A_vs_B",
                                stringsAsFactors = FALSE),
        metabolomics = data.frame(ID = "map00020", pvalue = 0.2, NES = 0.5,
                                  padj = 0.4, method = "fgsea", contrast = "A_vs_B",
                                  stringsAsFactors = FALSE))
    lk <- build_cross_omics_lookup(tabs, "proteomics", "metabolomics", top_n = 4,
                                   kegg_org = "rno", pathway_coverage = cov_fixture())
    all <- lk[lk$ranking == "all", ]
    st <- setNames(all$target_status, all$norm_id)
    expect_identical(unname(st["00010"]), "tested, not reported")
    expect_identical(unname(st["00030"]), "not testable: fewer than the minimum measured")
    expect_identical(unname(st["04110"]), "not testable: no KEGG members")
    expect_identical(unname(st["00040"]), "not testable: none measured")
    expect_true(all(c("target_n_in_pathway", "target_n_measured_features",
                      "target_coverage") %in% names(lk)))

    # Without coverage the statuses and columns are what they were.
    plain <- build_cross_omics_lookup(tabs, "proteomics", "metabolomics", top_n = 4,
                                      kegg_org = "rno")
    expect_true(all(plain$target_status == "not in target results"))
    expect_false("target_coverage" %in% names(plain))
})

test_that("a testable pathway is not called tested in a contrast the target lacks", {
    st <- .coverage_lookup_status(rep("not in target results", 2),
                                  c("testable", "testable"), c(TRUE, FALSE))
    expect_identical(st, c("tested, not reported", "not in target results"))
})

test_that("lookup cells follow the same states", {
    expect_identical(.lookup_cell_state(c("tested", "tested, not reported",
                                          "not testable: no KEGG members",
                                          "not in target results")),
                     c(NA, "neutral", "not_testable", NA))
})


# ---- Multi-ORA and MultiGSEA -------------------------------------------------

test_that("the Multi-ORA summary carries coverage, and the dot plot no fake p = 1", {
    pooled <- data.frame(pathway = c("Glycolysis", "Cell cycle"),
                         ID = c("rno00010", "rno04110"), pvalue = c(1e-4, 1e-3),
                         padj = c(1e-3, 1e-2), Count = c(5L, 4L),
                         stringsAsFactors = FALSE)
    per <- list(proteomics = pooled)
    combined <- build_multi_ora_summary(pooled, per, NULL,
                                        pathway_coverage = cov_fixture())
    expect_identical(combined$coverage_metabolomics, c("testable", "no_kegg_members"))
    expect_identical(combined$n_measured_metabolomics, c(4L, 0L))

    marks <- .multi_ora_coverage_markers(combined)
    expect_identical(marks$state, c("neutral", "not_testable"))
    expect_identical(unique(marks$source), "Metabolomics")

    src <- paste(deparse(body(plot_multi_ora_dotplot)), collapse = " ")
    expect_false(grepl("is.na(matched_padj), 1", src, fixed = TRUE))
})

test_that("MultiGSEA pairs record which layer could not test a term", {
    got <- .multigsea_coverage_columns(c("00010", "04110", "GO:0006915"),
                                       has1 = c(TRUE, TRUE, TRUE),
                                       has2 = c(FALSE, FALSE, FALSE),
                                       "proteomics", "metabolomics",
                                       cov_fixture(), "rno")
    expect_true(all(is.na(got$coverage_x)))
    expect_identical(got$state_y, c("neutral", "not_testable", NA))
})


# ---- run-level coverage, from the cache only --------------------------------

local_metab_run <- function(env = parent.frame(), map_hmdb = TRUE) {
    tmp <- withr::local_tempdir(.local_envir = env)
    hmdb <- sprintf("HMDB%07d", 1:9)
    dir.create(file.path(tmp, "data"))
    # Eight of the nine map into the fixture's pathways, one to a compound in
    # none of them; without mapping, the file names none of the run's HMDB ids.
    map_ids <- if (map_hmdb) hmdb else sprintf("HMDB9%06d", 1:9)
    kegg <- c(cov_measured, "C99999")
    writeLines(c("HMDB\tKEGG", paste0(map_ids, "\t", kegg)),
               file.path(tmp, "data", "HMDB2kegg_cpd.Jan2026.v2.txt"))
    withr::local_envvar(c(PIPELINE_ROOT = tmp), .local_envir = env)
    withr::local_dir(tmp, .local_envir = env)
    list(
        tmp = tmp,
        harmonization_res = list(inputs = list(metabolomics = list(
            row_data = data.frame(feature_id = hmdb, HMDB = hmdb,
                                  stringsAsFactors = FALSE)))),
        de_results = list(metabolomics = list(de_tables = list(B_vs_A = data.frame(
            feature_id = hmdb, logFC = seq(1, -1, length.out = 9),
            P.Value = seq(0.01, 0.9, length.out = 9),
            adj.P.Val = seq(0.05, 0.95, length.out = 9),
            stringsAsFactors = FALSE)))),
        config = list(global = list(organism = "not a known organism"),
                      modes = list(multiomics = list(enrichment = list())))
    )
}

test_that("no cached compound table means no coverage, and no download", {
    run <- local_metab_run()
    expect_message(
        res <- pathway_coverage_for_run(run$de_results, run$harmonization_res,
                                        run$config,
                                        compound_cache_dir = file.path(run$tmp, "none"),
                                        gene_cache_dir = NULL),
        "no cached KEGG compound-pathway table")
    expect_null(res)
    src <- paste(deparse(body(pathway_coverage_for_run)), collapse = " ")
    expect_true(grepl(".cached_compound_pathways(", src, fixed = TRUE))
    expect_false(grepl("get_kegg_compound_pathways(", src, fixed = TRUE))
})

test_that("the gene-pathway getter reads its cache and fetches nothing when told not to", {
    dir <- withr::local_tempdir()
    expect_null(get_kegg_gene_pathways("rno", dir, download = FALSE))
    saveRDS(data.frame(pathway = "rno00010", gene = "24189"),
            file.path(dir, "kegg_gene_pathways_rno.rds"))
    expect_identical(get_kegg_gene_pathways("rno", dir, download = FALSE)$gene, "24189")
})

test_that("a run's metabolite coverage is read from the cache and summarised", {
    run <- local_metab_run()
    cache <- file.path(run$tmp, "cross", "metabolomics")
    dir.create(cache, recursive = TRUE)
    saveRDS(cov_cpd_fixture(), file.path(cache, "kegg_compound_pathways.rds"))

    res <- suppressMessages(pathway_coverage_for_run(
        run$de_results, run$harmonization_res, run$config,
        compound_cache_dir = cache, gene_cache_dir = NULL))
    s <- res$summary
    expect_identical(s$layer, "metabolomics")
    expect_identical(s$status, "mapped")
    expect_identical(s$n_features, 9L)
    expect_identical(s$n_distinct_ids, 9L)
    expect_identical(s$n_pathways_testable, 2L)

    out <- file.path(run$tmp, "cross")
    writeLines("stale", file.path(out, "kegg_pathway_coverage.csv"))
    write_pathway_coverage(res, out, extra_keys = c("rno04110", "GO:0006915"),
                           kegg_org = "rno")
    tab <- utils::read.csv(file.path(out, "kegg_pathway_coverage.csv"),
                           stringsAsFactors = FALSE, colClasses = c(norm_id = "character"))
    expect_identical(tab$coverage[tab$norm_id == "04110"], "no_kegg_members")
    expect_false("GO:0006915" %in% tab$norm_id)
    expect_true(file.exists(file.path(out, "kegg_pathway_coverage_summary.tsv")))

    # Nothing to write clears both files.
    write_pathway_coverage(NULL, out)
    expect_false(file.exists(file.path(out, "kegg_pathway_coverage.csv")))
    expect_false(file.exists(file.path(out, "kegg_pathway_coverage_summary.tsv")))
})

test_that("metabolites that map to no KEGG compound give no_id_mapping", {
    run <- local_metab_run(map_hmdb = FALSE)
    cache <- file.path(run$tmp, "metab")
    dir.create(cache)
    saveRDS(cov_cpd_fixture(), file.path(cache, "kegg_compound_pathways.rds"))
    res <- suppressMessages(pathway_coverage_for_run(
        run$de_results, run$harmonization_res, run$config,
        compound_cache_dir = cache, gene_cache_dir = NULL))
    expect_identical(res$summary$status, "no_id_mapping")
    expect_true(all(res$layers$metabolomics$coverage == "no_id_mapping"))
})

test_that("coverage can be switched off", {
    run <- local_metab_run()
    run$config$modes$multiomics$enrichment$pathway_coverage <- list(enabled = FALSE)
    expect_null(suppressMessages(pathway_coverage_for_run(
        run$de_results, run$harmonization_res, run$config,
        compound_cache_dir = run$tmp, gene_cache_dir = NULL)))
})


# ---- config and report ---------------------------------------------------------

test_that("coverage is on by default, in the reader, the validator and the template", {
    expect_true(.pathway_coverage_config(list())$enabled)
    expect_false(.pathway_coverage_config(list(modes = list(multiomics = list(
        enrichment = list(pathway_coverage = list(enabled = FALSE))))))$enabled)
    cfg <- validate_multiomics_config(list(
        integration = list(methods = c("DIABLO", "SNF")),
        condition_column = "Group"))
    expect_true(cfg$enrichment$pathway_coverage$enabled)
    tpl <- yaml::read_yaml(file.path(root_dir, "config", "templates",
                                     "multiomics_config.yaml"))
    expect_true(tpl$modes$multiomics$enrichment$pathway_coverage$enabled)
})

test_that("every coverage legend is plain and within 100 words", {
    for (fig in c("heatmap", "lookup", "multi_ora", "multigsea")) {
        txt <- pathway_coverage_legend(fig)
        n <- lengths(regmatches(txt, gregexpr("[^[:space:]]+", txt)))
        expect_lte(n, 99, label = fig)
        expect_false(grepl("**", txt, fixed = TRUE))
        expect_true(grepl("cross", txt, fixed = TRUE), info = fig)
        expect_true(grepl("none were measured", txt, fixed = TRUE), info = fig)
    }
})

test_that("the report shows coverage legends only where the tables carry coverage", {
    dir <- withr::local_tempdir()
    with_cov <- file.path(dir, "a.csv")
    without <- file.path(dir, "b.csv")
    utils::write.csv(data.frame(norm_id = "00010", coverage_metabolomics = "testable"),
                     with_cov, row.names = FALSE)
    utils::write.csv(data.frame(norm_id = "00010", pval_a = 0.1), without,
                     row.names = FALSE)
    expect_true(has_pathway_coverage(with_cov))
    expect_false(has_pathway_coverage(without))
    expect_false(has_pathway_coverage(file.path(dir, "missing.csv")))

    src <- paste(readLines(file.path(root_dir, "R", "domain", "multiomics",
                                     "report_template_multiomics.Rmd"), warn = FALSE),
                 collapse = "\n")
    for (fig in c("heatmap", "lookup", "multi_ora", "multigsea")) {
        expect_true(grepl(sprintf('pathway_coverage_legend("%s")', fig), src,
                          fixed = TRUE), info = fig)
    }
    # Each heatmap legend sits behind the gate on the meta-analysis table.
    expect_gte(lengths(regmatches(src, gregexpr(
        "has_pathway_coverage(file.path(", src, fixed = TRUE))), 4)
    expect_true(grepl("Metabolites measured (n / N in pathway)", src, fixed = TRUE))
})
