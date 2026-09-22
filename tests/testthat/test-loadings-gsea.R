# GSEA on DIABLO loadings and MOFA2 weights (07e_loadings_gsea.R).
#
# All fixtures are synthetic: made-up feature ids, sample names and values. No
# network: the KEGG branch is not exercised here, and the metabolite mapper is
# stubbed.

# Stub by assignment into the functions' own environment. with_mocked_bindings()
# errors with "No packages loaded with pkgload" here, since these are sourced
# functions rather than a package. File-local, as in test-loadings-gmt-fallback.R.
local_stubs <- function(stubs, env = parent.frame()) {
    target <- environment(run_loadings_gsea)
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

# 60 features with signed values; SET_UP holds the 15 largest, SET_DOWN the 15
# smallest, SET_MIX alternates.
synthetic_values <- function() {
    stats::setNames(seq(3, by = -0.1, length.out = 60), paste0("G", 1:60))
}

write_synthetic_gmt <- function(path) {
    lines <- c(
        paste(c("SET_UP", "Set up", paste0("G", 1:15)), collapse = "\t"),
        paste(c("SET_DOWN", "Set down", paste0("G", 46:60)), collapse = "\t"),
        paste(c("SET_MIX", "Set mixed", paste0("G", seq(2, 60, by = 4))), collapse = "\t")
    )
    writeLines(lines, path)
    path
}

gmt_config <- function(gmt_path, method = NULL) {
    cfg <- list(
        project = list(dir = tempdir()),
        params = list(seed = 7),
        global = list(organism = "Synthetic organism"),
        modes = list(
            proteomics = list(pathway = list(gmt_file = gmt_path)),
            multiomics = list(condition_column = "group",
                              enrichment = list(loadings = list(min_size = 5,
                                                                max_size = 100)))))
    if (!is.null(method)) cfg$modes$multiomics$enrichment$loadings$method <- method
    cfg
}


# ---- ranking ---------------------------------------------------------------

test_that("the ranking keeps the sign and sorts decreasing", {
    r <- collapse_loading_ranks(c(-0.5, 0.9, 0.1), c("a", "b", "c"))
    expect_identical(names(r), c("b", "c", "a"))
    expect_equal(unname(r), c(0.9, 0.1, -0.5))
})

test_that("duplicate keys keep the largest absolute value", {
    r <- collapse_loading_ranks(c(0.4, -0.9, 0.2), c("g1", "g1", "g2"))
    expect_equal(r[["g1"]], -0.9)
    expect_length(r, 2)
})

test_that("a +x / -x tie goes to the positive value, whatever the row order", {
    a <- collapse_loading_ranks(c(0.5, -0.5), c("g1", "g1"))
    b <- collapse_loading_ranks(c(-0.5, 0.5), c("g1", "g1"))
    expect_equal(a[["g1"]], 0.5)
    expect_identical(a, b)
})

test_that("missing keys and non-finite values are dropped", {
    r <- collapse_loading_ranks(c(1, 2, NA, Inf), c("a", NA, "c", "d"))
    expect_identical(names(r), "a")
    expect_identical(collapse_loading_ranks(numeric(0), character(0)), numeric(0))
})

test_that("GENE_N ids resolve in place, keeping the vector's length", {
    harm <- list(gene_protein_mapping = data.frame(
        gene_id = c("GID1", "GID2"), protein_id = c("PID1", "PID2"),
        stringsAsFactors = FALSE))
    ids <- c("GENE_2", "GENE_9", "GENE_1_proteomics", "P_OTHER")
    out <- .resolve_feature_ids_aligned(ids, harm, "proteomics")
    expect_identical(out, c("PID2", NA, "PID1", "P_OTHER"))
    expect_identical(.resolve_feature_ids_aligned(ids, harm, "metabolomics"), ids)
})


# ---- entries ---------------------------------------------------------------

diablo_fixture <- function() {
    list(
        loadings = list(
            proteomics = matrix(c(0.5, -0.2, 0.1, 0.3, -0.4, 0.0), ncol = 2,
                                dimnames = list(c("P1", "P2", "P3"),
                                                c("comp1", "comp2"))),
            Y = matrix(c(1, -1), ncol = 1, dimnames = list(c("A", "B"), "comp1"))),
        sample_scores = list(
            proteomics = matrix(c(2, 1, -1, -2, 0.5, -0.5, 0.1, -0.1), ncol = 2,
                                dimnames = list(paste0("s", 1:4),
                                                c("comp1", "comp2"))),
            Y = matrix(0, 4, 2, dimnames = list(paste0("s", 1:4),
                                                c("comp1", "comp2")))))
}

test_that("DIABLO entries skip the outcome block and keep signed loadings", {
    e <- diablo_loading_entries(diablo_fixture())
    expect_setequal(names(e), c("DIABLO_proteomics_comp1", "DIABLO_proteomics_comp2"))
    one <- e[["DIABLO_proteomics_comp1"]]
    expect_identical(one$integration, "DIABLO")
    expect_identical(one$axis, "comp1")
    expect_equal(one$values, c(P1 = 0.5, P2 = -0.2, P3 = 0.1))
    expect_equal(one$scores, c(s1 = 2, s2 = 1, s3 = -1, s4 = -2))
})

test_that("MOFA2 entries take the leading factors and their sample scores", {
    w <- matrix(1:8 / 10, ncol = 4,
                dimnames = list(c("m1", "m2"), paste0("Factor", 1:4)))
    factors <- data.frame(sample_id = c("s1", "s2"), Factor1 = c(1, -1),
                          Factor2 = c(0, 2), Factor3 = c(3, 3), Factor4 = c(1, 1))
    e <- mofa_weight_entries(list(weights = list(metabolomics = w),
                                  factors = factors), n_factors = 3)
    expect_identical(names(e), paste0("MOFA_metabolomics_Factor", 1:3))
    expect_equal(e[["MOFA_metabolomics_Factor2"]]$scores, c(s1 = 0, s2 = 2))
    expect_equal(e[["MOFA_metabolomics_Factor1"]]$values, c(m1 = 0.1, m2 = 0.2))
})


# ---- orientation -----------------------------------------------------------

test_that("the orientation names the condition that scores higher", {
    o <- axis_orientation(c(s1 = 2, s2 = 1, s3 = -1, s4 = -2),
                          c(s1 = "B", s2 = "B", s3 = "A", s4 = "A"), "comp1")
    expect_identical(o$higher, "B")
    expect_identical(o$lower, "A")
    expect_match(o$note, "positive NES points towards B")
})

test_that("a tie in the group means breaks on the condition name", {
    o <- axis_orientation(c(s1 = 1, s2 = 1), c(s1 = "Z", s2 = "A"), "Factor1")
    expect_identical(o$higher, "A")
})

test_that("with one condition or no scores the orientation is not stated", {
    one <- axis_orientation(c(s1 = 1, s2 = 2), c(s1 = "A", s2 = "A"), "comp1")
    expect_true(is.na(one$higher))
    expect_true(is.na(axis_orientation(NULL, c(s1 = "A"), "comp1")$lower))
})


# ---- config ----------------------------------------------------------------

test_that("GSEA is the default and ORA has to be asked for", {
    expect_identical(loadings_enrichment_method(list()), "gsea")
    cfg <- list(modes = list(multiomics = list(enrichment = list(
        loadings = list(method = "ORA")))))
    expect_identical(loadings_enrichment_method(cfg), "ora")
})

test_that("an unknown method falls back to GSEA with a warning", {
    cfg <- list(modes = list(multiomics = list(enrichment = list(
        loadings = list(method = "fisher")))))
    expect_warning(m <- loadings_enrichment_method(cfg), "gsea")
    expect_identical(m, "gsea")
})

test_that("the validator fills GSEA and the size bounds, and rejects a typo", {
    validated <- suppressMessages(suppressWarnings(validate_multiomics_config(
        list(integration = list(methods = "SNF")))))
    expect_identical(validated$enrichment$loadings$method, "gsea")
    expect_equal(validated$enrichment$loadings$min_size, 10)
    expect_equal(validated$enrichment$loadings$max_size, 500)

    expect_error(suppressMessages(suppressWarnings(validate_multiomics_config(
        list(integration = list(methods = "SNF"),
             enrichment = list(loadings = list(method = "fisher")))))),
        "loadings.method")
})


# ---- dispatch --------------------------------------------------------------

test_that("the default runs GSEA and records it", {
    out_dir <- withr::local_tempdir()
    calls <- character(0)
    local_stubs(list(
        run_loadings_gsea = function(...) { calls <<- c(calls, "gsea"); list() },
        run_diablo_loadings_enrichment = function(...) { calls <<- c(calls, "ora"); NULL },
        run_mofa_weights_enrichment = function(...) { calls <<- c(calls, "ora"); NULL }))

    suppressMessages(run_loadings_enrichment(
        list(diablo_results = list(), mofa_results = list()), list(),
        list(global = list(organism = "Synthetic organism")), out_dir))

    expect_identical(calls, "gsea")
    expect_identical(readLines(file.path(out_dir, "loadings_enrichment_method.txt")),
                     "gsea")
})

test_that("method ora runs the ORA path and records it", {
    out_dir <- withr::local_tempdir()
    calls <- character(0)
    local_stubs(list(
        run_loadings_gsea = function(...) { calls <<- c(calls, "gsea"); list() },
        run_diablo_loadings_enrichment = function(...) { calls <<- c(calls, "diablo_ora"); NULL },
        run_mofa_weights_enrichment = function(...) { calls <<- c(calls, "mofa_ora"); NULL }))

    cfg <- list(global = list(organism = "Synthetic organism"),
                modes = list(multiomics = list(enrichment = list(
                    loadings = list(method = "ora")))))
    suppressMessages(run_loadings_enrichment(
        list(diablo_results = list(x = 1), mofa_results = list(x = 1)), list(),
        cfg, out_dir))

    expect_identical(calls, c("diablo_ora", "mofa_ora"))
    expect_identical(readLines(file.path(out_dir, "loadings_enrichment_method.txt")),
                     "ora")
})


# ---- scoring ---------------------------------------------------------------

test_that("GMT-based loadings GSEA is reproducible and signed", {
    skip_if_not_installed("fgsea")
    gmt <- write_synthetic_gmt(file.path(withr::local_tempdir(), "sets.gmt"))
    cfg <- gmt_config(gmt)

    run <- function() suppressMessages(suppressWarnings(gene_loadings_gsea(
        synthetic_values(), "proteomics", harmonization_res = list(), config = cfg,
        kegg_org = NULL, org_db = NULL, min_size = 5, max_size = 100, seed = 7)))
    a <- run()
    b <- run()

    expect_identical(a, b)
    expect_gt(a$NES[a$ID == "SET_UP"], 0)
    expect_lt(a$NES[a$ID == "SET_DOWN"], 0)
    expect_identical(a$pathway[a$ID == "SET_UP"], "Set up")
})

test_that("metabolite loadings are scored on their signed values", {
    skip_if_not_installed("fgsea")
    ids <- paste0("feature_", 1:40)
    vals <- stats::setNames(seq(2, by = -0.1, length.out = 40), ids)
    cpd <- sprintf("C%05d", 1:40)
    local_stubs(list(map_metabolite_ids_to_kegg = function(de_tables, harmonization_res) {
        data.frame(feature_id = ids, KEGG_CPD = cpd, stringsAsFactors = FALSE)
    }))
    cpd_pathways <- data.frame(
        pathway = rep(c("map90001", "map90002"), each = 8),
        compound = c(cpd[1:8], cpd[33:40]),
        name = rep(c("Pathway high", "Pathway low"), each = 8),
        stringsAsFactors = FALSE)

    res <- suppressMessages(metabolite_loadings_gsea(
        vals, harmonization_res = list(), cache_dir = withr::local_tempdir(),
        min_size = 3, max_size = 100, seed = 7, cpd_pathways = cpd_pathways))

    expect_gt(res$NES[res$ID == "map90001"], 0)
    expect_lt(res$NES[res$ID == "map90002"], 0)
})


# ---- plot and files --------------------------------------------------------

gsea_table <- function() {
    data.frame(pathway = c("Set up", "Set down", "Set flat"),
               ID = c("SET_UP", "SET_DOWN", "SET_FLAT"),
               NES = c(2.4, -2.1, 0.3), pvalue = c(0.001, 0.004, 0.8),
               padj = c(0.003, 0.006, 0.8), stringsAsFactors = FALSE)
}

test_that("the NES plot carries dashed guides at -2 and +2", {
    skip_if_not_installed("ggplot2")
    p <- plot_loadings_gsea_nes(gsea_table(), "DIABLO_proteomics_comp1",
                                subtitle = "B scores higher than A on comp1")
    expect_true(inherits(p, "ggplot"))

    vlines <- Filter(function(l) inherits(l$geom, "GeomVline"), p$layers)
    guides <- Filter(function(l) identical(l$aes_params$linetype, "dashed"), vlines)
    expect_length(guides, 1)
    expect_setequal(guides[[1]]$data$xintercept, c(-2, 2))
    expect_identical(guides[[1]]$aes_params$colour, "red")
})

test_that("an empty table draws nothing", {
    expect_null(plot_loadings_gsea_nes(gsea_table()[0, ], "x"))
})

test_that("entries are written with their orientation, replacing a previous run", {
    skip_if_not_installed("ggplot2")
    out_dir <- withr::local_tempdir()
    writeLines("old", file.path(out_dir, "DIABLO_proteomics_comp9_gsea.csv"))
    writeLines("keep", file.path(out_dir, "DIABLO_proteomics_comp1_enrichment.csv"))
    local_stubs(list(gene_loadings_gsea = function(...) gsea_table()))

    entries <- diablo_loading_entries(diablo_fixture())
    ctx <- list(harmonization_res = list(), config = list(), kegg_org = NULL,
                org_db = NULL, min_size = 5, max_size = 100, seed = 7,
                exclude_classes = NULL, cpd_pathways = NULL,
                conditions = c(s1 = "B", s2 = "B", s3 = "A", s4 = "A"),
                cache_dir = out_dir)
    combined <- suppressMessages(run_loadings_gsea_entries(
        entries, ctx, out_dir, "diablo_loadings_gsea_all.csv"))

    expect_false(file.exists(file.path(out_dir, "DIABLO_proteomics_comp9_gsea.csv")))
    expect_true(file.exists(file.path(out_dir, "DIABLO_proteomics_comp1_enrichment.csv")))
    for (lab in c("DIABLO_proteomics_comp1", "DIABLO_proteomics_comp2")) {
        expect_true(file.exists(file.path(out_dir, paste0(lab, "_gsea.csv"))))
        expect_true(file.exists(file.path(out_dir, paste0(lab, "_gsea_nes.png"))))
    }
    expect_true(file.exists(file.path(out_dir, "diablo_loadings_gsea_all.csv")))
    comp1 <- combined[combined$axis == "comp1", ]
    expect_true(all(comp1$higher_scoring_group == "B"))
    expect_true(all(comp1$integration == "DIABLO"))
})
