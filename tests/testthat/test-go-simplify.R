# Tests for the optional GO term simplification of the pathway step
# (R/core/09b_go_simplify.R): config resolution and validation, the clustering
# rule, the per-table simplification, file outputs and module wiring.

go_simplify_has_deps <- function() {
    requireNamespace("GO.db", quietly = TRUE) &&
        requireNamespace("GOSemSim", quietly = TRUE)
}

# Synthetic fgsea-shaped table. Three GO BP groups (death, immune, cell cycle),
# two GO CC terms, one GO term that is not significant and one non-GO row.
go_simplify_fixture <- function() {
    data.frame(
        pathway      = c("GO:0006915", "GO:0012501", "GO:0008219",
                         "GO:0006955", "GO:0002376",
                         "GO:0007049", "GO:0000278",
                         "GO:0005737", "GO:0005829",
                         "GO:0006412", "hsa04110"),
        pathway_name = c("apoptotic process", "programmed cell death", "cell death",
                         "immune response", "immune system process",
                         "cell cycle", "mitotic cell cycle",
                         "cytoplasm", "cytosol",
                         "translation", "Cell cycle (KEGG)"),
        pval = c(0.0004, 0.0001, 0.0009, 0.0002, 0.0030, 0.0005, 0.0003,
                 0.0010, 0.0020, 0.4000, 0.0001),
        padj = c(0.004, 0.002, 0.006, 0.003, 0.012, 0.005, 0.004,
                 0.007, 0.010, 0.500, 0.002),
        NES  = c(2.1, 2.3, 1.9, -2.0, -1.7, 2.2, 2.4, 1.5, 1.4, 0.3, 2.0),
        stringsAsFactors = FALSE
    )
}

# ---- settings ----------------------------------------------------------------

test_that("resolve_go_simplify_settings defaults to off, 0.7, Wang", {
    s <- resolve_go_simplify_settings(NULL)
    expect_false(s$enabled)
    expect_equal(s$threshold, 0.7)
    expect_equal(s$measure, "Wang")
    expect_null(s$orgdb)

    s <- resolve_go_simplify_settings(list(simplify_go = TRUE, simplify_threshold = 0.5))
    expect_true(s$enabled)
    expect_equal(s$threshold, 0.5)
})

test_that("resolve_go_simplify_settings finds an OrgDb only for IC measures", {
    s <- resolve_go_simplify_settings(list(simplify_go = TRUE), "Homo sapiens")
    expect_null(s$orgdb)
    s <- resolve_go_simplify_settings(list(simplify_go = TRUE, simplify_measure = "Rel"),
                                      "Homo sapiens")
    expect_equal(s$orgdb, "org.Hs.eg.db")
})

# ---- validation --------------------------------------------------------------

test_that("validate_go_simplify_config checks types, range and measure", {
    expect_silent(validate_go_simplify_config(NULL))
    expect_silent(validate_go_simplify_config(list()))
    expect_silent(validate_go_simplify_config(list(
        simplify_go = TRUE, simplify_threshold = 0.4, simplify_measure = "Lin")))

    expect_error(validate_go_simplify_config(list(simplify_go = "yes")), "simplify_go")
    expect_error(validate_go_simplify_config(list(simplify_threshold = "0.7")),
                 "simplify_threshold")
    expect_error(validate_go_simplify_config(list(simplify_threshold = 1.5)),
                 "simplify_threshold")
    expect_error(validate_go_simplify_config(list(simplify_threshold = -0.1)),
                 "simplify_threshold")
    expect_error(validate_go_simplify_config(list(simplify_measure = "TCSS")),
                 "simplify_measure")
})

test_that("the RNA and proteomics validators check the simplify keys", {
    rna <- yaml::read_yaml(file.path(root_dir, "config", "templates", "rna_config.yaml"))
    rna_cfg <- rna$modes$rna
    rna_cfg$enrichment$annotation_dir <- NULL
    rna_cfg$pathway$simplify_measure <- "bogus"
    expect_error(validate_rna_config(rna_cfg), "simplify_measure")

    prot <- yaml::read_yaml(file.path(root_dir, "config", "templates", "proteins_config.yaml"))
    prot_cfg <- prot$modes$proteomics
    prot_cfg$pathway$simplify_threshold <- 2
    expect_error(suppressWarnings(validate_proteomics_config(prot_cfg)), "simplify_threshold")
})

test_that("both templates carry the simplify keys, off by default", {
    for (f in c("rna_config.yaml", "proteins_config.yaml")) {
        cfg <- yaml::read_yaml(file.path(root_dir, "config", "templates", f))
        pw <- cfg$modes[[if (f == "rna_config.yaml") "rna" else "proteomics"]]$pathway
        expect_identical(pw$simplify_go, FALSE, info = f)
        expect_equal(pw$simplify_threshold, 0.7, info = f)
        expect_identical(pw$simplify_measure, "Wang", info = f)
    }
})

# ---- clustering rule -----------------------------------------------------------

test_that("collapse_go_terms keeps the lowest p-value term of each group", {
    ids <- c("GO:0000001", "GO:0000002", "GO:0000003", "GO:0000004")
    sim <- matrix(c(1.0, 0.9, 0.8, 0.1,
                    0.9, 1.0, 0.85, 0.1,
                    0.8, 0.85, 1.0, 0.1,
                    0.1, 0.1, 0.1, 1.0), 4, dimnames = list(ids, ids))
    pval <- setNames(c(0.03, 0.01, 0.02, 0.04), ids)
    res <- collapse_go_terms(sim, pval, threshold = 0.7)

    expect_equal(res$term, ids)
    expect_equal(length(unique(res$cluster)), 2)
    expect_equal(res$representative, c(rep("GO:0000002", 3), "GO:0000004"))
})

test_that("collapse_go_terms breaks p-value ties on adjusted p, then order", {
    ids <- c("GO:0000001", "GO:0000002", "GO:0000003")
    sim <- matrix(0.9, 3, 3, dimnames = list(ids, ids)); diag(sim) <- 1
    pval <- setNames(c(0.01, 0.01, 0.01), ids)
    expect_equal(unique(collapse_go_terms(sim, pval, setNames(c(0.05, 0.02, 0.02), ids))$representative),
                 "GO:0000002")
    expect_equal(unique(collapse_go_terms(sim, pval)$representative), "GO:0000001")
})

test_that("a higher threshold merges at least as much as a lower one", {
    ids <- sprintf("GO:%07d", 1:5)
    x <- matrix(withr::with_seed(1, runif(25)), 5); sim <- (x + t(x)) / 2; diag(sim) <- 1
    dimnames(sim) <- list(ids, ids)
    pval <- setNames(c(0.01, 0.02, 0.03, 0.04, 0.05), ids)
    n_groups <- vapply(c(0, 0.3, 0.6, 0.9),
                       function(h) length(unique(collapse_go_terms(sim, pval, threshold = h)$cluster)),
                       integer(1))
    expect_equal(n_groups[1], 5L)
    expect_true(all(diff(n_groups) <= 0))
})

# ---- one table ---------------------------------------------------------------

test_that("simplify_go_results drops grouped duplicates and leaves the rest", {
    skip_if_not(go_simplify_has_deps(), "GO.db / GOSemSim not installed")
    df <- go_simplify_fixture()
    settings <- resolve_go_simplify_settings(list(simplify_go = TRUE))
    res <- suppressMessages(simplify_go_results(df, settings))

    expect_equal(res$status, "simplified")
    expect_lt(nrow(res$simplified), nrow(df))

    # Non-GO and non-significant GO rows pass through untouched and ungrouped.
    for (id in c("hsa04110", "GO:0006412")) {
        expect_true(id %in% res$simplified$pathway, info = id)
        expect_false(id %in% res$map$pathway, info = id)
        kept <- res$simplified[res$simplified$pathway == id, names(df)]
        rownames(kept) <- NULL
        orig <- df[df$pathway == id, ]
        rownames(orig) <- NULL
        expect_equal(kept, orig, info = id)
    }

    # Every group keeps exactly its lowest p-value member, and only that one.
    for (cl in unique(res$map$go_cluster)) {
        grp <- res$map[res$map$go_cluster == cl, ]
        best <- grp$pathway[order(grp$pval, grp$padj)][1]
        expect_equal(unique(grp$representative), best)
        expect_true(best %in% res$simplified$pathway)
        expect_false(any(setdiff(grp$pathway, best) %in% res$simplified$pathway))
    }

    # Groups never mix GO ontologies.
    expect_true(all(tapply(res$map$ontology, res$map$go_cluster,
                           function(x) length(unique(x))) == 1))
    expect_setequal(unique(res$map$ontology), c("BP", "CC"))

    # The collapsed_terms column lists what each representative stands for.
    collapsed <- unlist(strsplit(res$simplified$collapsed_terms, ";", fixed = TRUE))
    expect_setequal(collapsed, res$map$pathway[!res$map$is_representative])
})

test_that("simplify_go_results with threshold 0 keeps every term", {
    skip_if_not(go_simplify_has_deps(), "GO.db / GOSemSim not installed")
    df <- go_simplify_fixture()
    settings <- resolve_go_simplify_settings(list(simplify_go = TRUE, simplify_threshold = 0))
    res <- suppressMessages(simplify_go_results(df, settings))
    expect_equal(nrow(res$simplified), nrow(df))
    expect_equal(res$simplified$pathway, df$pathway)
})

test_that("simplify_go_results skips tables without enough significant GO terms", {
    df <- data.frame(pathway = c("GO:0006915", "hsa04110"), pval = c(0.01, 0.02),
                     padj = c(0.02, 0.03), stringsAsFactors = FALSE)
    res <- simplify_go_results(df, resolve_go_simplify_settings(list(simplify_go = TRUE)))
    expect_match(res$status, "^skipped")
    expect_identical(res$simplified, df)
    expect_null(res$map)
})

# ---- whole pathway result + files ------------------------------------------

test_that("simplify_pathway_results skips non-GO tables and leaves its input alone", {
    skip_if_not(go_simplify_has_deps(), "GO.db / GOSemSim not installed")
    df <- go_simplify_fixture()
    kegg <- df[df$pathway == "hsa04110", ]
    pr <- list(c1 = list(GO_BP_fgsea = df, KEGG_fgsea = kegg))
    pr_before <- pr
    settings <- resolve_go_simplify_settings(list(simplify_go = TRUE))

    out <- suppressMessages(simplify_pathway_results(pr, settings))
    expect_identical(pr, pr_before)
    expect_equal(names(out$tables), "c1|GO_BP_fgsea")
    expect_equal(out$summary$analysis, "GO_BP_fgsea")
    expect_equal(out$summary$n_terms_full, nrow(df))
    expect_lte(out$summary$n_terms_simplified, out$summary$n_terms_full)

    dir <- withr::local_tempdir()
    files <- suppressMessages(write_simplified_pathway_results(out, dir))
    expect_setequal(basename(files), c("pathway_c1_GO_BP_fgsea_simplified.csv",
                                       "pathway_c1_GO_BP_fgsea_simplified_map.csv",
                                       "go_simplification_summary.csv"))

    txt <- describe_go_simplification(file.path(dir, "go_simplification_summary.csv"))
    expect_match(txt, "Wang measure")
    expect_match(txt, "full enrichment tables are unchanged")
})

test_that("describe_go_simplification is silent when nothing was simplified", {
    expect_null(describe_go_simplification(NULL))
    expect_null(describe_go_simplification(tempfile()))
    f <- withr::local_tempfile(fileext = ".csv")
    write.csv(data.frame(status = "skipped: no GO terms", sig_cutoff = 0.1,
                         max_terms = 500, threshold = 0.7, measure = "Wang"),
              f, row.names = FALSE)
    expect_null(describe_go_simplification(f))
})

test_that("run_go_simplification does nothing when simplify_go is off", {
    dir <- file.path(withr::local_tempdir(), "simplified")
    pr <- list(c1 = list(GO_BP_fgsea = go_simplify_fixture()))
    expect_null(run_go_simplification(pr, list(), "Homo sapiens", dir))
    expect_null(run_go_simplification(pr, list(simplify_go = FALSE), "Homo sapiens", dir))
    expect_false(dir.exists(dir))
})

# ---- module wiring: off = identical outputs ---------------------------------

test_that("mod_rnaseq_pathway output is the same with simplify_go off, absent or on", {
    skip_if_not_installed("fgsea")

    gmt <- withr::local_tempfile(fileext = ".gmt")
    writeLines(c(
        "GO:0006915\tapoptotic process\tg1\tg2\tg3\tg4\tg5",
        "GO:0012501\tprogrammed cell death\tg1\tg2\tg3\tg4\tg5\tg6",
        "GO:0007049\tcell cycle\tg6\tg7\tg8\tg9\tg10\tg11\tg12",
        "SET_X\tnot GO\tg9\tg10\tg11\tg12\tg13"
    ), gmt)
    de_res <- list(tables = list(contrast_A = data.frame(
        FeatureID      = paste0("g", 1:15),
        log2FoldChange = c(2.0, -1.5, 0.1, 1.2, -2.1, 0.5, 0.3, -0.7, 1.8, -1.1,
                           0.4, -0.2, 1.5, -1.8, 0.6),
        pvalue         = c(0.001, 0.01, 0.8, 0.04, 0.005, 0.2, 0.3, 0.15, 0.002, 0.03,
                           0.5, 0.7, 0.008, 0.001, 0.4),
        padj           = c(0.01, 0.05, 0.9, 0.1, 0.02, 0.3, 0.4, 0.25, 0.01, 0.08,
                           0.6, 0.8, 0.04, 0.01, 0.5),
        stat           = c(5, -3, 0.2, 2.5, -4, 1, 1.2, -1.5, 4.5, -2.8,
                           0.8, -0.3, 3.2, -4.5, 1.1),
        stringsAsFactors = FALSE)))

    run <- function(extra) {
        config <- list(modes = list(rna = list(
            annotation = list(skip_annotation = TRUE, organism = "Homo sapiens"),
            pathway    = c(list(enabled = TRUE, databases = character(0), gmt_file = gmt,
                                min_size = 3, max_size = 50), extra))))
        out_dir <- withr::local_tempdir()
        res <- suppressMessages(suppressWarnings(
            mod_rnaseq_pathway(de_res, list(expr_filt = NULL, meta = NULL), config, out_dir)))
        enrich <- file.path(out_dir, "Enrichment")
        csvs <- sort(list.files(enrich, pattern = "\\.csv$"))
        list(res = res$pathway_results,
             csvs = setNames(lapply(file.path(enrich, csvs), readLines), csvs),
             simplified_dir = dir.exists(file.path(enrich, "simplified")))
    }

    absent <- run(list())
    off    <- run(list(simplify_go = FALSE))
    on     <- run(list(simplify_go = TRUE))

    expect_identical(off$res, absent$res)
    expect_identical(off$csvs, absent$csvs)
    expect_false(absent$simplified_dir)
    expect_false(off$simplified_dir)

    # Turning it on adds files under simplified/ and changes nothing else.
    expect_identical(on$res, off$res)
    expect_identical(on$csvs, off$csvs)
    expect_true(on$simplified_dir)
})

# ---- methods text -------------------------------------------------------------

test_that("proteomics Methods mention GO simplification only when it ran", {
    cfg <- list(modes = list(proteomics = list(
        engine = "DIANN", scale_in = "linear",
        imputation = list(method = "none", multi_imputation = FALSE),
        de = list(method = "limma"),
        pathway = list(enabled = TRUE, method = "fgsea", databases = "GO",
                       simplify_go = TRUE))))
    f <- withr::local_tempfile(fileext = ".csv")

    # The flag alone is not evidence: no summary file, no sentence.
    down <- build_proteomics_methods_text(cfg, go_simplify_summary_file = f)$downstream
    expect_false(grepl("semantic similarity", down))

    write.csv(data.frame(status = "simplified", sig_cutoff = 0.1, max_terms = 500,
                         threshold = 0.7, measure = "Wang"), f, row.names = FALSE)
    down <- build_proteomics_methods_text(cfg, go_simplify_summary_file = f)$downstream
    expect_match(down, "grouped by semantic similarity (GOSemSim, Wang measure)", fixed = TRUE)
    expect_match(down, "cut at 0.7", fixed = TRUE)
})
