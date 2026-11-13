# tests/testthat/test-metab-enrichment-provenance.R
#
# Provenance of the metabolomics ORA and GSEA results
# (R/domain/metabolomics/06r_enrichment_provenance.R).
#
# ORA and GSEA test COMPOUNDS, and several measured features can map to one
# compound. These tests cover:
#   1. the statistics did not move: table, ranks, contrast and method equal a
#      frozen copy of the functions as they were before provenance was added
#      (fixtures/metab_enrichment_provenance/baseline_1f145fb.R);
#   2. which step chose each feature's compound id;
#   3. ORA: universe, foreground, significance and overlap kept apart, with
#      compound-based counts;
#   4. GSEA: the representative feature that was actually kept, ties included,
#      and the members fgsea used as opposed to the members it was given;
#   5. pathways: repeated ids, equal labels across libraries, repeated members;
#   6. a provenance that fails or disagrees leaves the result alone and is
#      refused where provenance is required.
#
# Everything is synthetic and lives in a temporary folder.

prov_env <- new.env()
sys.source(test_path("fixtures", "metab_enrichment_provenance",
                     "baseline_1f145fb.R"), envir = prov_env)

# ---- the synthetic case ---------------------------------------------------------
#
# 40 plain features g01..g40, one each for compounds C00001..C00040, plus:
#   m1a m1b m1c  three features of C00101. m1b has the largest |logFC| but is
#                not significant; m1c is significant with a smaller one.
#   t1a t1b      two features of C00102 with the SAME logFC (tie, same sign).
#   o1a o1b      two features of C00103 with +x and -x (tie, opposite signs).
#   x1           KEGG cell "C00104;C00105": the first id is used.
#   h1           no KEGG cell; HMDB0000201 maps to C00106 through the file.
#   n1           name only: dropped, because the data carries KEGG ids.
#   k1           KEGG cell C00107 AND an HMDB id mapping elsewhere: KEGG wins.
#   z1           C00109, but its logFC is missing.
#   d1           C00108, measured but absent from the DE table.
#   e1           in the DE table only, not measured.
#
# Libraries (every line has a description):
#   libA  PW1 PW2 PW3 PW4 PW5, PW6 twice, PW7      libB  PW7 again, PW8
#   PW2 lists HMDB0000201 and C00106: the same compound twice after translation.
#   PW3 lists C00108 and C00109 (measured, not ranked) and C00500 (not measured).
#   PW4 has one matching member. PW5 is HMDB0000201 + C00106 only.
make_prov_case <- function(reverse_ties = FALSE, env = parent.frame()) {
    dir <- withr::local_tempdir(.local_envir = env)

    base_feat <- sprintf("g%02d", 1:40)
    base_kegg <- sprintf("C%05d", 1:40)
    special <- c("m1a", "m1b", "m1c", "t1a", "t1b", "o1a", "o1b",
                 "x1", "h1", "n1", "k1", "z1", "d1")
    row_data <- data.frame(
        feature_id = c(base_feat, special),
        Name = c(paste0("cmpd", 1:40), "multiA", "multiB", "multiC", "tieA", "tieB",
                 "oppA", "oppB", "Ambig", "HmdbOnly", "Citrate", "Both", "NoFC", "NoDE"),
        HMDB = c(rep(NA_character_, 40), NA, NA, NA, NA, NA, NA, NA,
                 NA, "HMDB0000201", NA, "HMDB0000207", NA, NA),
        KEGG = c(base_kegg, "C00101", "C00101", "C00101", "C00102", "C00102",
                 "C00103", "C00103", "C00104;C00105", NA, NA, "cpd:C00107",
                 "C00109", "C00108"),
        stringsAsFactors = FALSE)
    expr <- withr::with_seed(1, matrix(
        stats::rnorm(nrow(row_data) * 6, mean = 10), nrow = nrow(row_data),
        dimnames = list(row_data$feature_id, paste0("S", 1:6))))

    i <- 1:40
    de <- data.frame(
        feature_id = c(base_feat, "m1a", "m1b", "m1c", "t1a", "t1b", "o1a", "o1b",
                       "x1", "h1", "n1", "k1", "z1", "e1"),
        logFC = c((-1)^i * (0.6 + 0.05 * i), 0.5, -2.013, 1.013, 1.523, 1.523,
                  1.733, -1.733, 0.913, -1.113, 0.9, 1.313, NA, 0.777),
        stringsAsFactors = FALSE)
    significant <- c(base_feat[1:6], "m1c", "t1a", "x1", "n1", "z1", "e1")
    de$P.Value   <- ifelse(de$feature_id %in% significant, 0.001, 0.5)
    de$adj.P.Val <- de$P.Value
    if (reverse_ties) {
        swap <- function(a, b) {
            ia <- match(a, de$feature_id); ib <- match(b, de$feature_id)
            de[c(ia, ib), ] <<- de[c(ib, ia), ]
        }
        swap("t1a", "t1b")
        swap("o1a", "o1b")
        rownames(de) <- NULL
    }

    gmt_line <- function(id, desc, members) paste(c(id, desc, members), collapse = "\t")
    lib_a <- file.path(dir, "libA.gmt")
    lib_b <- file.path(dir, "libB.gmt")
    writeLines(c(
        gmt_line("PW1", "First",  c(base_kegg[1:10], "C00101", "C00104")),
        gmt_line("PW2", "Second", c(base_kegg[11:20], "C00102", "C00103",
                                    "HMDB0000201", "C00106")),
        gmt_line("PW3", "Third",  c(base_kegg[21:30], "C00108", "C00109", "C00500")),
        gmt_line("PW4", "Small",  c("C00031", "C00777")),
        gmt_line("PW5", "Dup",    c("HMDB0000201", "C00106")),
        gmt_line("PW6", "Repeat", c("C00001", "C00002")),
        gmt_line("PW6", "Repeat", base_kegg[31:35]),
        gmt_line("PW7", "Shared", base_kegg[1:5])
    ), lib_a)
    writeLines(c(
        gmt_line("PW7", "Shared", c(base_kegg[36:40], "C00101")),
        gmt_line("PW8", "Other",  base_kegg[12:18])
    ), lib_b)
    mapping <- file.path(dir, "hmdb_kegg.tsv")
    writeLines(c("HMDB\tKEGG", "HMDB0000201\tC00106", "HMDB0000207\tC00999"), mapping)
    missing_gmt <- file.path(dir, "absent.gmt")

    list(
        pre = list(row_data = row_data, expr_raw = expr),
        de_res = list(de_tables = list("A vs B" = de)),
        config = list(
            project = list(dir = dir), paths = list(raw = "data"),
            modes = list(metabolomics = list(
                de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5),
                enrichment = list(run_enrichment = TRUE,
                                  gmt_file = list(lib_a, lib_b, missing_gmt),
                                  mapping_file = mapping)))),
        files = list(lib_a = lib_a, lib_b = lib_b, mapping = mapping,
                     missing = missing_gmt)
    )
}

# Silence the functions' progress messages and fgsea's own warnings, but let a
# provenance warning through: an unexpected one must show up as a warning.
keep_provenance_warnings <- function(expr) {
    withCallingHandlers(suppressMessages(expr), warning = function(w) {
        if (!grepl("provenance", conditionMessage(w), fixed = TRUE)) {
            invokeRestart("muffleWarning")
        }
    })
}

run_ora <- function(case, fun = run_metabolomics_ora) {
    keep_provenance_warnings(
        fun(case$pre, case$de_res, case$config, contrast = "A vs B"))
}

# fgsea draws random numbers; the same seed on both sides makes the two runs
# comparable to the last digit.
run_gsea <- function(case, fun = run_metabolomics_gsea, seed = 11) {
    keep_provenance_warnings(withr::with_seed(
        seed, fun(case$pre, case$de_res, case$config, contrast = "A vs B")))
}

row_of <- function(df, ...) {
    cond <- list(...)
    keep <- rep(TRUE, nrow(df))
    for (nm in names(cond)) keep <- keep & !is.na(df[[nm]]) & df[[nm]] == cond[[nm]]
    df[keep, , drop = FALSE]
}

# Replace a function where the pipeline functions live, for one test.
local_prov_stub <- function(name, fn, env = parent.frame()) {
    target <- environment(map_compounds_for_enrichment)
    old <- get(name, envir = target, inherits = FALSE)
    withr::defer(assign(name, old, envir = target), envir = env)
    assign(name, fn, envir = target)
    invisible(NULL)
}

# ---- 1. the statistics did not move ------------------------------------------------

test_that("read_gmt_list() returns what it returned before, plus records on request", {
    case <- make_prov_case()
    for (f in c(case$files$lib_a, case$files$lib_b, case$files$missing)) {
        expect_identical(read_gmt_list(f), prov_env$baseline_read_gmt_list(f))
        expect_identical(
            read_gmt_list(f, include_descriptions = TRUE),
            prov_env$baseline_read_gmt_list(f, include_descriptions = TRUE))
        with_records <- read_gmt_list(f, include_descriptions = TRUE,
                                      include_records = TRUE)
        expect_identical(with_records[c("sets", "descriptions")],
                         prov_env$baseline_read_gmt_list(f, include_descriptions = TRUE))
    }

    rec <- read_gmt_list(case$files$lib_a, include_descriptions = TRUE,
                         include_records = TRUE)
    ids <- vapply(rec$records, function(r) r$pathway_id, character(1))
    expect_identical(ids, c("PW1", "PW2", "PW3", "PW4", "PW5", "PW6", "PW6", "PW7"))
    expect_identical(vapply(rec$records, function(r) r$line, integer(1)), 1:8)
    # A repeated id keeps its last line in `sets`; both lines are on record.
    expect_identical(rec$sets$PW6, rec$records[[7]]$members)
    expect_identical(rec$records[[6]]$members, c("C00001", "C00002"))
    expect_identical(read_gmt_list(case$files$missing, include_descriptions = TRUE,
                                   include_records = TRUE)$records, list())
})

test_that("map_compounds_for_enrichment() returns the same mapping as before", {
    case <- make_prov_case()
    args <- list(case$pre$row_data, case$pre$expr_raw, case$files$mapping)
    now  <- suppressMessages(do.call(map_compounds_for_enrichment, args))
    base <- suppressMessages(do.call(prov_env$baseline_map_compounds_for_enrichment, args))

    expect_identical(names(now), c(names(base), "feature_compounds"))
    expect_identical(now[names(base)], base)

    # Also without a mapping file, and without the KEGG restriction.
    for (extra in list(list(NULL), list(case$files$mapping, annotated_only = FALSE))) {
        a <- c(list(case$pre$row_data, case$pre$expr_raw), extra)
        expect_identical(
            suppressMessages(do.call(map_compounds_for_enrichment, a))[names(base)],
            suppressMessages(do.call(prov_env$baseline_map_compounds_for_enrichment, a)))
    }
})

test_that("ORA returns the same table, contrast and method as before", {
    for (reverse in c(FALSE, TRUE)) {
        case <- make_prov_case(reverse_ties = reverse)
        now  <- run_ora(case)
        base <- run_ora(case, prov_env$baseline_run_metabolomics_ora)

        expect_identical(names(base), c("table", "contrast", "method"))
        expect_identical(names(now), c(names(base), "provenance"))
        expect_identical(now[names(base)], base)
        expect_identical(now$provenance$status, "ok")
    }
})

test_that("GSEA returns the same table, ranks, contrast and method as before", {
    skip_if_not_installed("fgsea")
    for (reverse in c(FALSE, TRUE)) {
        case <- make_prov_case(reverse_ties = reverse)
        now  <- run_gsea(case)
        base <- run_gsea(case, prov_env$baseline_run_metabolomics_gsea)

        expect_identical(names(base), c("table", "ranks", "contrast", "method"))
        expect_identical(names(now), c(names(base), "provenance"))
        expect_identical(now$ranks, base$ranks)
        expect_identical(now$table, base$table)
        expect_identical(now[c("contrast", "method")], base[c("contrast", "method")])
        expect_identical(now$provenance$status, "ok")
    }
})

test_that("the early exits still return NULL", {
    case <- make_prov_case()
    off <- case$config
    off$modes$metabolomics$enrichment$run_enrichment <- FALSE
    expect_null(suppressMessages(run_metabolomics_ora(case$pre, case$de_res, off)))
    expect_null(suppressMessages(run_metabolomics_gsea(case$pre, case$de_res, off)))

    none <- case$de_res
    none$de_tables[[1]]$P.Value <- 0.5
    expect_null(suppressMessages(run_metabolomics_ora(case$pre, none, case$config)))
})

# ---- 2. which step chose each compound id -------------------------------------------

test_that("id_source is the step that made the final choice", {
    case <- make_prov_case()
    fc <- suppressMessages(map_compounds_for_enrichment(
        case$pre$row_data, case$pre$expr_raw, case$files$mapping))$feature_compounds

    expect_s3_class(fc, "data.frame")
    expect_identical(names(fc), c("feature_id", "compound_id", "id_source",
                                  "upstream_origin", "excluded_reason", "kegg_cell",
                                  "kegg_candidates", "n_kegg_candidates", "ambiguous"))
    expect_identical(fc$feature_id, case$pre$row_data$feature_id)

    # A KEGG cell wins, and where it came from is not known to this stage.
    g <- row_of(fc, feature_id = "g01")
    expect_identical(c(g$compound_id, g$id_source, g$upstream_origin),
                     c("C00001", "kegg_column", "not_tracked"))
    expect_false(g$ambiguous)

    # ... also over an HMDB id that maps elsewhere.
    k <- row_of(fc, feature_id = "k1")
    expect_identical(c(k$compound_id, k$id_source), c("C00107", "kegg_column"))
    expect_identical(k$kegg_cell, "cpd:C00107")

    # No KEGG cell: the mapping file this stage read supplies the id.
    h <- row_of(fc, feature_id = "h1")
    expect_identical(c(h$compound_id, h$id_source, h$upstream_origin),
                     c("C00106", "hmdb_mapping", "enrichment_mapping_file"))
    expect_identical(h$n_kegg_candidates, 0L)

    # Several ids in one cell: the first is used, all are on record.
    x <- row_of(fc, feature_id = "x1")
    expect_identical(x$compound_id, "C00104")
    expect_identical(x$kegg_candidates, "C00104;C00105")
    expect_identical(x$n_kegg_candidates, 2L)
    expect_true(x$ambiguous)

    # A name-only feature is excluded when the data carries KEGG ids.
    n <- row_of(fc, feature_id = "n1")
    expect_true(is.na(n$compound_id))
    expect_true(is.na(n$id_source))
    expect_identical(n$excluded_reason, "not_kegg_annotated")

    expect_identical(sum(is.na(fc$excluded_reason)), 52L)
})

test_that("names and matrix row names are recorded as the source when they are used", {
    expr <- matrix(1:12, nrow = 3, dimnames = list(c("f1", "f2", "f3"), paste0("S", 1:4)))

    named <- data.frame(feature_id = c("f1", "f2", "f3"),
                        Name = c("Glucose", NA, "Citrate"), stringsAsFactors = FALSE)
    fc <- suppressMessages(map_compounds_for_enrichment(named, expr))$feature_compounds
    expect_identical(fc$compound_id, c("Glucose", NA, "Citrate"))
    expect_identical(fc$id_source, c("name_column", NA, "name_column"))
    expect_identical(fc$upstream_origin, c("annotation", NA, "annotation"))
    expect_identical(fc$excluded_reason, c(NA, "empty_id", NA))
    expect_false(any(fc$ambiguous))

    bare <- data.frame(feature_id = c("f1", "f2", "f3"), stringsAsFactors = FALSE)
    fc <- suppressMessages(map_compounds_for_enrichment(bare, expr))$feature_compounds
    expect_identical(fc$compound_id, c("f1", "f2", "f3"))
    expect_identical(unique(fc$id_source), "matrix_rowname")
})

test_that("describe_metab_feature_compounds() refuses vectors that do not line up", {
    expect_error(
        describe_metab_feature_compounds(
            feature_id = c("f1", "f2"), compound_id = "C00001", valid = TRUE,
            id_source = "kegg_column", restricted = FALSE, kegg_cell = "C00001"),
        "one entry per feature")
    expect_error(
        describe_metab_feature_compounds(
            feature_id = NULL, compound_id = "C00001", valid = TRUE,
            id_source = "kegg_column", restricted = FALSE, kegg_cell = "C00001"),
        "one entry per feature")
})

# ---- 3. ORA ---------------------------------------------------------------------------

test_that("ORA provenance has the agreed elements and passes its own checks", {
    case <- make_prov_case()
    res <- run_ora(case)
    p <- res$provenance

    expect_identical(names(p), c("schema_version", "method", "contrast", "status",
                                 "diagnostics", "error", "reference_files",
                                 "feature_compounds", "pathway_sources",
                                 "pathway_members", "thresholds", "compounds",
                                 "features"))
    expect_identical(p$schema_version, 1L)
    expect_identical(c(p$method, p$contrast, p$status), c("ora", "A vs B", "ok"))
    expect_null(p$error)
    expect_true(all(p$diagnostics$passed))
    expect_setequal(p$diagnostics$check,
                    c("universe_size", "universe_is_measured_compounds",
                      "foreground_size", "foreground_from_significant_features",
                      "library_sources", "pathway_inputs", "pathway_size", "overlap",
                      "overlap_compounds"))
    expect_identical(p$thresholds, list(pvalue_column = "P.Value", p_cutoff = 0.05,
                                        logfc_cutoff = log2(1.5)))
    expect_identical(require_metab_enrichment_provenance(res), p)
})

test_that("ORA counts are compound-based and agree with the table", {
    case <- make_prov_case()
    res <- run_ora(case)
    p <- res$provenance
    cmp <- p$compounds

    # 40 plain compounds + C00101-C00104, C00106-C00109. Not C00105 (second id
    # of an ambiguous cell) and not the name-only feature.
    expect_identical(sum(cmp$in_universe), 48L)
    expect_identical(unique(res$table$bg_size), 48L)
    expect_identical(sort(cmp$compound_id[cmp$in_foreground]),
                     sort(c(sprintf("C%05d", 1:6), "C00101", "C00102", "C00104")))
    expect_identical(unique(res$table$fg_size), 9L)

    # Three features, one compound, counted once; one of them is significant.
    multi <- row_of(cmp, compound_id = "C00101")
    expect_identical(c(multi$n_features, multi$n_significant_features), c(3L, 1L))
    expect_true(multi$in_foreground)
    expect_identical(sum(cmp$n_features), 52L)
})

test_that("ORA keeps significance, universe and foreground apart for each feature", {
    case <- make_prov_case()
    f <- run_ora(case)$provenance$features
    flags <- function(id) {
        r <- row_of(f, feature_id = id)
        c(de = r$in_de_table, sig = r$significant, universe = r$compound_in_universe,
          foreground = r$compound_in_foreground)
    }

    # The significant feature of a foreground compound ...
    expect_identical(flags("m1c"), c(de = TRUE, sig = TRUE, universe = TRUE, foreground = TRUE))
    # ... and its siblings, which are not significant but share the compound.
    expect_identical(flags("m1b"), c(de = TRUE, sig = FALSE, universe = TRUE, foreground = TRUE))
    expect_identical(flags("m1a"), c(de = TRUE, sig = FALSE, universe = TRUE, foreground = TRUE))
    # In the universe without being significant.
    expect_identical(flags("g20"), c(de = TRUE, sig = FALSE, universe = TRUE, foreground = FALSE))
    # Measured, absent from the DE table: universe only.
    expect_identical(flags("d1"), c(de = FALSE, sig = FALSE, universe = TRUE, foreground = FALSE))
    # Significant but without a usable compound: in neither set.
    expect_identical(flags("n1"), c(de = TRUE, sig = TRUE, universe = FALSE, foreground = FALSE))
    expect_identical(flags("e1"), c(de = TRUE, sig = TRUE, universe = FALSE, foreground = FALSE))
    # A missing logFC is not significant, whatever its p-value.
    expect_identical(flags("z1"), c(de = TRUE, sig = FALSE, universe = TRUE, foreground = FALSE))

    expect_true(is.na(row_of(f, feature_id = "e1")$compound_id))
    expect_identical(nrow(f), 54L)
})

test_that("ORA pathway sizes count every listed occurrence, as the test does", {
    case <- make_prov_case()
    res <- run_ora(case)
    m <- res$provenance$pathway_members
    tbl <- res$table

    expect_identical(names(m), c("pathway", "position", "member_as_listed",
                                 "member_after_translation", "compound_id",
                                 "occurrence", "in_input_set",
                                 "in_effective_statistic", "in_foreground"))
    # For ORA the two are the same thing.
    expect_identical(m$in_effective_statistic, m$in_input_set)

    for (r in seq_len(nrow(tbl))) {
        rows <- m[m$pathway == tbl$pathway[r] & m$in_effective_statistic, ]
        expect_identical(nrow(rows), as.integer(tbl$pathway_size[r]))
        expect_identical(length(unique(rows$compound_id[rows$in_foreground])),
                         as.integer(tbl$overlap[r]))
    }

    # PW2 lists one compound twice (an HMDB id and its KEGG id): two rows, both
    # counted in the size, so 14 and not 13.
    dup <- row_of(m, pathway = "PW2 - Second", compound_id = "C00106")
    expect_identical(dup$member_as_listed, c("HMDB0000201", "C00106"))
    expect_identical(dup$member_after_translation, c("C00106", "C00106"))
    expect_identical(dup$occurrence, 1:2)
    expect_identical(dup$in_effective_statistic, c(TRUE, TRUE))
    expect_identical(row_of(tbl, pathway = "PW2 - Second")$pathway_size, 14L)
    expect_identical(row_of(tbl, pathway = "PW2 - Second")$overlap, 1L)
    # PW5 is that compound twice and nothing else: ORA tests it with size 2.
    expect_identical(row_of(tbl, pathway = "PW5 - Dup")$pathway_size, 2L)

    # A listed member that is not measured is on record and not counted.
    out <- row_of(m, pathway = "PW3 - Third", member_as_listed = "C00500")
    expect_true(is.na(out$compound_id))
    expect_false(out$in_input_set)
    expect_true(is.na(out$occurrence))

    expect_identical(row_of(tbl, pathway = "PW1 - First")$overlap, 8L)
})

# ---- 4. GSEA ---------------------------------------------------------------------------

test_that("GSEA provenance has the agreed elements and passes its own checks", {
    skip_if_not_installed("fgsea")
    case <- make_prov_case()
    res <- run_gsea(case)
    p <- res$provenance

    expect_identical(names(p), c("schema_version", "method", "contrast", "status",
                                 "diagnostics", "error", "reference_files",
                                 "feature_compounds", "pathway_sources",
                                 "pathway_members", "ranked", "features"))
    expect_identical(c(p$method, p$contrast, p$status), c("gsea", "A vs B", "ok"))
    expect_true(all(p$diagnostics$passed))
    expect_setequal(p$diagnostics$check,
                    c("ranked_follows_ranks", "ranked_universe", "representative_logfc",
                      "representative_is_largest_abs_logfc", "library_sources",
                      "pathway_inputs", "pathway_size", "leading_edge"))
    expect_identical(require_metab_enrichment_provenance(res), p)
})

test_that("`ranked` follows the final ranks vector entry by entry", {
    skip_if_not_installed("fgsea")
    case <- make_prov_case()
    res <- run_gsea(case)
    ranked <- res$provenance$ranked

    expect_identical(names(ranked), c("rank", "compound_id", "logFC",
                                      "representative_feature", "n_features",
                                      "tied_abs_logfc", "tied_opposite_sign"))
    # The ranked universe: compounds with a DE row and a logFC. C00108 (not in
    # the DE table) and C00109 (missing logFC) are measured but not ranked.
    expect_identical(nrow(ranked), 46L)
    expect_identical(ranked$compound_id, names(res$ranks))
    expect_identical(ranked$logFC, unname(res$ranks))
    expect_identical(ranked$rank, seq_along(res$ranks))
    expect_false(any(c("C00108", "C00109") %in% ranked$compound_id))
    expect_false(anyDuplicated(ranked$representative_feature) > 0)

    # Each representative is a real DE row of that compound with that logFC.
    de <- case$de_res$de_tables[[1]]
    expect_identical(de$logFC[match(ranked$representative_feature, de$feature_id)],
                     ranked$logFC)
})

test_that("the representative is the feature the deduplication kept", {
    skip_if_not_installed("fgsea")
    case <- make_prov_case()
    ranked <- run_gsea(case)$provenance$ranked

    # Largest |logFC| of three; it is not the significant one.
    multi <- row_of(ranked, compound_id = "C00101")
    expect_identical(multi$representative_feature, "m1b")
    expect_identical(multi$logFC, -2.013)
    expect_identical(multi$n_features, 3L)
    expect_false(multi$tied_abs_logfc)
    expect_false(multi$tied_opposite_sign)

    plain <- row_of(ranked, compound_id = "C00001")
    expect_identical(c(plain$representative_feature, as.character(plain$n_features)),
                     c("g01", "1"))
    expect_false(plain$tied_abs_logfc)
})

test_that("a tie goes to the earlier DE row, and is recorded as a tie", {
    skip_if_not_installed("fgsea")
    case_forward  <- make_prov_case()
    case_backward <- make_prov_case(reverse_ties = TRUE)
    forward  <- run_gsea(case_forward)
    backward <- run_gsea(case_backward)
    rep_of <- function(res, compound) {
        row_of(res$provenance$ranked, compound_id = compound)
    }

    # Same sign: the logFC cannot tell the two apart, the row order decides.
    expect_identical(rep_of(forward, "C00102")$representative_feature, "t1a")
    expect_identical(rep_of(backward, "C00102")$representative_feature, "t1b")
    expect_identical(forward$ranks[["C00102"]], backward$ranks[["C00102"]])
    for (res in list(forward, backward)) {
        expect_true(rep_of(res, "C00102")$tied_abs_logfc)
        expect_false(rep_of(res, "C00102")$tied_opposite_sign)
    }

    # Opposite signs: the row order decides the SIGN the compound is ranked
    # with. Provenance names the feature that was used in each case.
    expect_identical(rep_of(forward, "C00103")$representative_feature, "o1a")
    expect_identical(forward$ranks[["C00103"]], 1.733)
    expect_identical(rep_of(backward, "C00103")$representative_feature, "o1b")
    expect_identical(backward$ranks[["C00103"]], -1.733)
    for (res in list(forward, backward)) {
        expect_true(rep_of(res, "C00103")$tied_abs_logfc)
        expect_true(rep_of(res, "C00103")$tied_opposite_sign)
        expect_identical(res$provenance$status, "ok")
    }
})

test_that("every feature has a role, and only representatives are ranked", {
    skip_if_not_installed("fgsea")
    case <- make_prov_case()
    res <- run_gsea(case)
    f <- res$provenance$features
    role <- function(id) row_of(f, feature_id = id)$role

    expect_identical(names(f), c("feature_id", "compound_id", "in_de_table", "logFC", "role"))
    expect_identical(role("m1b"), "representative")
    expect_identical(c(role("m1a"), role("m1c")), rep("not_representative", 2))
    expect_identical(c(role("t1a"), role("t1b")), c("representative", "not_representative"))
    expect_identical(c(role("n1"), role("e1")), rep("no_compound", 2))
    expect_identical(role("z1"), "missing_logfc")
    expect_identical(role("d1"), "not_in_de_table")
    expect_false(row_of(f, feature_id = "d1")$in_de_table)

    expect_identical(sum(f$role == "representative"), length(res$ranks))
    expect_setequal(f$feature_id[f$role == "representative"],
                    res$provenance$ranked$representative_feature)
    expect_identical(nrow(f), 54L)
})

test_that("GSEA keeps the occurrences given to fgsea apart from the members it used", {
    skip_if_not_installed("fgsea")
    case <- make_prov_case()
    res <- run_gsea(case)
    m <- res$provenance$pathway_members
    tbl <- res$table

    expect_identical(names(m), c("pathway", "position", "member_as_listed",
                                 "member_after_translation", "compound_id",
                                 "occurrence", "in_input_set",
                                 "in_effective_statistic", "measured_compound_id",
                                 "measured", "in_leading_edge"))
    expect_true(all(m$in_input_set[m$in_effective_statistic]))

    # fgsea's size is the number of DISTINCT ranked members.
    for (r in seq_len(nrow(tbl))) {
        rows <- m[m$pathway == tbl$pathway[r], ]
        expect_identical(sum(rows$in_effective_statistic), as.integer(tbl$pathway_size[r]))
        expect_identical(sum(rows$in_effective_statistic),
                         length(unique(rows$compound_id[rows$in_input_set])))
    }

    # PW2: both occurrences were given to fgsea, one counted. 14 given, 13 used.
    pw2 <- m[m$pathway == "PW2 - Second", ]
    dup <- row_of(pw2, compound_id = "C00106")
    expect_identical(dup$member_as_listed, c("HMDB0000201", "C00106"))
    expect_identical(dup$in_input_set, c(TRUE, TRUE))
    expect_identical(dup$in_effective_statistic, c(TRUE, FALSE))
    expect_identical(sum(pw2$in_input_set), 14L)
    expect_identical(row_of(tbl, pathway = "PW2 - Second")$pathway_size, 13L)

    # PW5 is one compound listed twice. It passed the pipeline's own filter (two
    # matches) and fgsea returned nothing for it: given, never used.
    pw5 <- m[m$pathway == "PW5 - Dup", ]
    expect_identical(pw5$in_input_set, c(TRUE, TRUE))
    expect_identical(pw5$in_effective_statistic, c(FALSE, FALSE))
    expect_false(any(pw5$in_leading_edge))
    expect_false("PW5 - Dup" %in% tbl$pathway)
})

test_that("measured members that were not ranked never count towards the score", {
    skip_if_not_installed("fgsea")
    case <- make_prov_case()
    m <- run_gsea(case)$provenance$pathway_members
    pw3 <- m[m$pathway == "PW3 - Third", ]

    for (id in c("C00108", "C00109")) {
        r <- row_of(pw3, member_as_listed = id)
        expect_true(r$measured)
        expect_identical(r$measured_compound_id, id)
        expect_true(is.na(r$compound_id))
        expect_identical(c(r$in_input_set, r$in_effective_statistic, r$in_leading_edge),
                         c(FALSE, FALSE, FALSE))
    }
    absent <- row_of(pw3, member_as_listed = "C00500")
    expect_false(absent$measured)
    expect_false(absent$in_input_set)
    expect_identical(sum(pw3$in_effective_statistic), 10L)
})

test_that("the leading edge maps to compounds fgsea used, and to their representatives", {
    skip_if_not_installed("fgsea")
    case <- make_prov_case()
    res <- run_gsea(case)
    m <- res$provenance$pathway_members
    ranked <- res$provenance$ranked
    membership <- metab_pathway_feature_membership(res$provenance)

    expect_true(all(m$in_effective_statistic[m$in_leading_edge]))
    for (r in seq_len(nrow(res$table))) {
        pw <- res$table$pathway[r]
        edge <- strsplit(res$table$leading_edge[r], ";", fixed = TRUE)[[1]]
        expect_setequal(m$compound_id[m$pathway == pw & m$in_leading_edge], edge)

        # The features of the leading edge are the representatives, one each.
        feats <- membership$feature_id[membership$pathway == pw & membership$in_leading_edge]
        expect_setequal(feats, ranked$representative_feature[match(edge, ranked$compound_id)])
    }
    expect_true(sum(m$in_leading_edge) > 0)
})

# ---- 5. pathways ------------------------------------------------------------------------

test_that("every GMT line is on record with what became of it", {
    skip_if_not_installed("fgsea")
    case <- make_prov_case()
    ora  <- run_ora(case)$provenance$pathway_sources
    gsea <- run_gsea(case)$provenance$pathway_sources

    expect_identical(names(ora), c("gmt_file", "library", "line", "pathway_id",
                                   "description", "pathway", "n_listed", "fate"))
    for (s in list(ora, gsea)) {
        expect_identical(s$library, c(rep("libA", 8), rep("libB", 2)))
        expect_identical(s$line, c(1:8, 1:2))
        expect_identical(s$pathway_id, c("PW1", "PW2", "PW3", "PW4", "PW5", "PW6",
                                         "PW6", "PW7", "PW7", "PW8"))
        expect_identical(s$pathway[1], "PW1 - First")
        expect_identical(s$n_listed, c(12L, 14L, 13L, 2L, 2L, 2L, 5L, 5L, 6L, 7L))
        expect_identical(unique(s$gmt_file[s$library == "libB"]), case$files$lib_b)
    }

    expect_identical(ora$fate, c(
        "tested", "tested", "tested", "too_few_matches", "tested",
        "replaced_same_id_in_file", "tested", "replaced_same_label_later_library",
        "tested", "tested"))
    # The same lines for GSEA, except the set fgsea returned nothing for.
    expect_identical(gsea$fate, replace(ora$fate, 5, "not_returned"))
})

test_that("a label shared by two libraries is tested with the later library's members", {
    case <- make_prov_case()
    res <- run_ora(case)
    m <- res$provenance$pathway_members
    shared <- m[m$pathway == "PW7 - Shared", ]

    expect_identical(shared$member_as_listed, c(sprintf("C%05d", 36:40), "C00101"))
    expect_identical(row_of(res$table, pathway = "PW7 - Shared")$library, "libB")
    expect_identical(row_of(res$table, pathway = "PW7 - Shared")$pathway_size, 6L)
    # The repeated id keeps its second line.
    expect_identical(m$member_as_listed[m$pathway == "PW6 - Repeat"], sprintf("C%05d", 31:35))
})

test_that("the reference files that were read are recorded with their checksums", {
    case <- make_prov_case()
    ref <- run_ora(case)$provenance$reference_files

    expect_identical(names(ref), c("role", "path", "md5", "read"))
    expect_identical(ref$role, c("gmt", "gmt", "gmt", "mapping"))
    expect_identical(ref$path, c(case$files$lib_a, case$files$lib_b,
                                 case$files$missing, case$files$mapping))
    expect_identical(ref$read, c(TRUE, TRUE, FALSE, TRUE))
    expect_identical(ref$md5[c(1, 2, 4)], unname(tools::md5sum(
        c(case$files$lib_a, case$files$lib_b, case$files$mapping))))
    # A configured library that does not exist was skipped, and says so.
    expect_true(is.na(ref$md5[3]))
})

# ---- feature membership for display --------------------------------------------------------

test_that("ORA membership lists every feature of a pathway compound, with separate roles", {
    case <- make_prov_case()
    res <- run_ora(case)
    mb <- metab_pathway_feature_membership(res$provenance)

    expect_identical(names(mb), c("pathway", "pathway_id", "library", "compound_id",
                                  "feature_id", "member_occurrences",
                                  "in_effective_statistic", "compound_in_universe",
                                  "feature_significant", "compound_in_foreground",
                                  "compound_in_overlap"))

    multi <- row_of(mb, pathway = "PW1 - First", compound_id = "C00101")
    expect_identical(multi$feature_id, c("m1a", "m1b", "m1c"))
    expect_identical(multi$feature_significant, c(FALSE, FALSE, TRUE))
    expect_identical(multi$compound_in_foreground, rep(TRUE, 3))
    expect_identical(multi$compound_in_overlap, rep(TRUE, 3))
    expect_identical(unique(multi$pathway_id), "PW1")
    expect_identical(unique(multi$library), "libA")
    expect_identical(unique(row_of(mb, pathway = "PW7 - Shared")$library), "libB")

    # One row for the compound listed twice, with the count kept.
    dup <- row_of(mb, pathway = "PW2 - Second", compound_id = "C00106")
    expect_identical(dup$feature_id, "h1")
    expect_identical(dup$member_occurrences, 2L)
    expect_false(dup$feature_significant)
    expect_false(dup$compound_in_overlap)

    # The compounds in the overlap are the ones the table counted.
    for (r in seq_len(nrow(res$table))) {
        rows <- mb[mb$pathway == res$table$pathway[r], ]
        expect_identical(length(unique(rows$compound_id[rows$compound_in_overlap])),
                         as.integer(res$table$overlap[r]))
        expect_true(all(rows$compound_in_universe))
    }
})

test_that("GSEA membership never credits a feature that did not enter the score", {
    skip_if_not_installed("fgsea")
    case <- make_prov_case()
    res <- run_gsea(case)
    mb <- metab_pathway_feature_membership(res$provenance)

    expect_identical(names(mb), c("pathway", "pathway_id", "library", "compound_id",
                                  "feature_id", "member_occurrences",
                                  "in_effective_statistic", "compound_ranked",
                                  "feature_is_representative", "contributes_to_score",
                                  "in_leading_edge"))

    # Hard rules, for every row.
    expect_true(all(mb$compound_ranked[mb$contributes_to_score]))
    expect_true(all(mb$feature_is_representative[mb$contributes_to_score]))
    expect_true(all(mb$in_effective_statistic[mb$contributes_to_score]))
    expect_true(all(mb$contributes_to_score[mb$in_leading_edge]))

    # Three features of one compound: all shown, one contributes.
    multi <- row_of(mb, pathway = "PW1 - First", compound_id = "C00101")
    expect_identical(multi$feature_id, c("m1a", "m1b", "m1c"))
    expect_identical(multi$feature_is_representative, c(FALSE, TRUE, FALSE))
    expect_identical(multi$contributes_to_score, c(FALSE, TRUE, FALSE))
    expect_identical(multi$in_leading_edge[c(1, 3)], c(FALSE, FALSE))

    # Measured, not ranked: shown, never contributing.
    for (id in c("C00108", "C00109")) {
        r <- row_of(mb, pathway = "PW3 - Third", compound_id = id)
        expect_identical(nrow(r), 1L)
        expect_identical(c(r$compound_ranked, r$feature_is_representative,
                           r$contributes_to_score, r$in_leading_edge), rep(FALSE, 4))
    }
    expect_identical(row_of(mb, pathway = "PW3 - Third", compound_id = "C00108")$feature_id, "d1")

    # A set fgsea returned nothing for: its feature is shown and credited nothing.
    pw5 <- row_of(mb, pathway = "PW5 - Dup")
    expect_identical(pw5$feature_id, "h1")
    expect_identical(pw5$member_occurrences, 2L)
    expect_true(pw5$compound_ranked)
    expect_false(pw5$contributes_to_score)

    # What contributes is exactly what fgsea sized.
    for (r in seq_len(nrow(res$table))) {
        rows <- mb[mb$pathway == res$table$pathway[r], ]
        expect_identical(sum(rows$contributes_to_score), as.integer(res$table$pathway_size[r]))
    }
})

# ---- 6. provenance that fails or disagrees ----------------------------------------------------

test_that("a provenance that cannot be built leaves the ORA result as it was", {
    case <- make_prov_case()
    base <- run_ora(case, prov_env$baseline_run_metabolomics_ora)
    local_prov_stub("describe_metab_feature_compounds", function(...) stop("boom"))

    expect_warning(res <- run_ora(case), "provenance could not be built")
    expect_identical(res[names(base)], base)

    p <- res$provenance
    expect_identical(names(p), c("schema_version", "method", "contrast", "status",
                                 "diagnostics", "error"))
    expect_identical(c(p$method, p$contrast, p$status), c("ora", "A vs B", "failed"))
    expect_identical(p$error$step, "build")
    expect_match(p$error$message, "boom")
    expect_identical(nrow(p$diagnostics), 0L)

    expect_error(require_metab_enrichment_provenance(res),
                 "Provenance is required for metabolomics ORA 'A vs B'.*could not be built.*boom")
    expect_error(metab_pathway_feature_membership(p), "cannot be used")
})

test_that("a provenance that cannot be built leaves the GSEA result as it was", {
    skip_if_not_installed("fgsea")
    case <- make_prov_case()
    base <- run_gsea(case, prov_env$baseline_run_metabolomics_gsea)
    local_prov_stub("describe_metab_feature_compounds", function(...) stop("boom"))

    expect_warning(res <- run_gsea(case), "provenance could not be built")
    expect_identical(res$table, base$table)
    expect_identical(res$ranks, base$ranks)
    expect_identical(res$provenance$status, "failed")
    expect_error(require_metab_enrichment_provenance(res), "could not be built")
})

test_that("a provenance that disagrees with the result is marked, kept and refused", {
    case <- make_prov_case()
    base <- run_ora(case, prov_env$baseline_run_metabolomics_ora)
    local_prov_stub(".metab_prov_sources_match", function(...) FALSE)

    expect_warning(res <- run_ora(case), "does not agree with the result")
    expect_identical(res[names(base)], base)

    p <- res$provenance
    expect_identical(p$status, "inconsistent")
    expect_null(p$error)
    expect_identical(p$diagnostics$check[!p$diagnostics$passed], "library_sources")
    # The content is kept, for whoever has to find out what disagreed.
    expect_true(is.data.frame(p$pathway_members) && nrow(p$pathway_members) > 0)

    expect_error(require_metab_enrichment_provenance(res),
                 "does not agree with the result it describes \\(failed: library_sources\\)")
    expect_error(metab_pathway_feature_membership(p), "cannot be used")
})

test_that("provenance that is absent or of another version is refused, never assumed", {
    case <- make_prov_case()
    legacy <- run_ora(case, prov_env$baseline_run_metabolomics_ora)
    expect_error(require_metab_enrichment_provenance(legacy), "there is no provenance record")
    expect_error(require_metab_enrichment_provenance(NULL), "there is no provenance record")
    expect_error(require_metab_enrichment_provenance(legacy, label = "the payload"),
                 "Provenance is required for the payload")

    res <- run_ora(case)
    expect_identical(validate_metab_enrichment_provenance(res$provenance), character(0))

    newer <- res
    newer$provenance$schema_version <- 2L
    expect_error(require_metab_enrichment_provenance(newer), "schema version is 2")

    partial <- res
    partial$provenance$compounds <- NULL
    expect_error(require_metab_enrichment_provenance(partial), "it lacks: compounds")
})
