# Tests for the per-run results fact sheet.
#
# The sheet pairs each headline number with the artifact it can be checked
# against, so the things worth pinning are: that a value never appears without
# a source, that a missing artifact drops rows instead of failing the run, and
# that the numbers match what the files actually contain.
#
# Base R only:
#   testthat::test_file("tests/testthat/test-results-fact-sheet.R")
#
# Test map:
#   F1 fact() / bind_facts() / fmt_range() building blocks
#   F2 de_summary_counts is read under either of the two schemas in the repo
#   F3 build_rnaseq_fact_sheet on a synthetic results directory
#   F4 missing artifacts degrade to fewer rows, never to an error
#   F5 every row carries a non-empty source

`%||%` <- function(a, b) if (is.null(a)) b else a

find_repo_file <- function(rel) {
    dir <- normalizePath(".", mustWork = FALSE)
    for (i in 1:8) {
        cand <- file.path(dir, rel)
        if (file.exists(cand)) return(cand)
        dir <- dirname(dir)
    }
    stop("Could not locate ", rel, " from working dir ", getwd())
}

source(find_repo_file("R/core/00_paths.R"))
source(find_repo_file("R/core/01_io.R"))
source(find_repo_file("R/core/results_fact_sheet.R"))
source(find_repo_file("R/domain/rnaseq/12_fact_sheet.R"))


# =============================================================================
# F1 — building blocks
# =============================================================================
test_that("F1 fact() drops rows with no value rather than writing a blank", {
    expect_null(fact("x", NULL, "f.tsv"))
    expect_null(fact("x", NA, "f.tsv"))
    r <- fact("x", 61, "f.tsv")
    expect_equal(r$value, "61")
    expect_equal(names(r), c("claim", "value", "source_file"))
})

test_that("F1 bind_facts drops NULLs and always returns the three columns", {
    s <- bind_facts(list(fact("a", 1, "f"), NULL, fact("b", 2, "g")), NULL)
    expect_equal(nrow(s), 2)
    expect_equal(names(s), c("claim", "value", "source_file"))
    expect_equal(nrow(bind_facts(NULL)), 0)
    expect_equal(names(bind_facts()), c("claim", "value", "source_file"))
})

test_that("F1 fmt_range reports min, max and median and ignores non-finite values", {
    expect_equal(fmt_range(c(1.68, 2.61, 351)), "1.68 to 351, median 2.61")
    expect_equal(fmt_range(c(1, NA, Inf, 3, 2)), "1 to 3, median 2")
    expect_null(fmt_range(c(NA, Inf)))
})


# =============================================================================
# F2 — the two de_summary_counts schemas in this repo
# =============================================================================
write_run <- function(root, de_schema = c("contrast", "Name", "none")) {
    de_schema <- match.arg(de_schema)
    ds <- file.path(root, "Datasets")
    dir.create(ds, recursive = TRUE, showWarnings = FALSE)

    counts <- data.frame(gene = paste0("g", 1:4),
                         S_1 = c(100, 10, 0, 40), S_2 = c(120, 12, 0, 44),
                         N_1 = c(20, 11, 30, 42), N_2 = c(24, 9, 34, 38),
                         check.names = FALSE)
    write.table(counts, file.path(ds, "rna_counts_filtered.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)
    write.table(counts, file.path(ds, "rna_norm_TMMlogCPM.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)

    final <- data.frame(
        Gene = paste0("g", 1:4),
        Mean.S = c(110, 11, 0, 42), Mean.N = c(22, 10, 32, 40),
        log2FC.S_vs_N = c(2.32, 0.14, -Inf, 0.07),
        linearFC.S_vs_N = c(5.0, 1.1, -1000, 1.05),
        pvalue.S_vs_N = c(1e-8, 0.4, 1e-5, 0.9),
        padj.S_vs_N = c(1e-6, 0.6, 1e-4, 0.95),
        pass_any_contrast = c(1, NA, 1, NA),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    write.table(final, file.path(ds, "final_results.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)

    if (de_schema == "contrast") {
        d <- data.frame(contrast = c("S_vs_N", "any"), up = c(1, 0),
                        down = c(1, 0), total = c(2, 2))
        write.table(d, file.path(root, "de_summary_counts.tsv"), sep = "\t",
                    quote = FALSE, row.names = FALSE)
    } else if (de_schema == "Name") {
        d <- data.frame(Name = c("S_vs_N", "pass_any"), up = c(1, 0),
                        down = c(1, 0), any = c(2, 2))
        write.table(d, file.path(ds, "de_summary_counts.tsv"), sep = "\t",
                    quote = FALSE, row.names = FALSE)
    }
    root
}

fixture_config <- function() {
    list(modes = list(rna = list(
        de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5, deseq_mode = "default"),
        effects = list(samples = "SampleID", color = "grp")
    )))
}
fixture_meta <- function() {
    data.frame(SampleID = c("S_1", "S_2", "N_1", "N_2"),
               grp = c("S", "S", "N", "N"), stringsAsFactors = FALSE)
}
fixture_contrasts <- function() {
    data.frame(Contrast_name = "S_vs_N", Factor = "grp",
               Numerator = "S", Denominator = "N", stringsAsFactors = FALSE)
}

test_that("F2 both de_summary_counts schemas produce the same DE row", {
    for (schema in c("contrast", "Name")) {
        root <- write_run(file.path(withr::local_tempdir(), "rna"), schema)
        d <- .read_de_summary_counts(root, create_legacy_output_dirs(root, create = FALSE))
        expect_equal(names(d)[1], "contrast", info = schema)
        expect_equal(nrow(d), 1, info = schema)   # the all-contrasts row is dropped
        expect_equal(d$contrast, "S_vs_N", info = schema)
        expect_equal(d$total, 2, info = schema)
    }
})

test_that("F2 no de_summary_counts file at all yields NULL, not an error", {
    root <- write_run(file.path(withr::local_tempdir(), "rna"), "none")
    expect_null(.read_de_summary_counts(root, create_legacy_output_dirs(root, create = FALSE)))
})


# =============================================================================
# F3 — the assembled sheet
# =============================================================================
test_that("F3 the sheet carries the numbers the files actually contain", {
    root <- write_run(file.path(withr::local_tempdir(), "rna"), "contrast")
    s <- build_rnaseq_fact_sheet(root, fixture_config(),
                                 pre = list(meta = fixture_meta()),
                                 inputs = list(contrasts = fixture_contrasts()))
    val <- function(claim) s$value[s$claim == claim]

    expect_equal(val("contrasts tested"), "S_vs_N")
    expect_equal(val("libraries analysed"), "4")
    expect_equal(val("genes after expression filtering"), "4")
    expect_equal(val("genes receiving an adjusted p-value"), "4")
    expect_equal(val("differentially expressed genes, S_vs_N"), "2 total (1 up, 1 down)")
    expect_equal(val("genes at raw p < 0.05, S_vs_N"), "2")
    expect_equal(val("significance rule"),
                 "adjusted p <= 0.05 and |linear fold change| >= 1.5")
    # g3 has Mean.S == 0 and passes, so it is the one on-off gene of the two
    expect_match(val("differentially expressed genes with a zero group mean"), "^1 of 2")
    # and the fold-change range must then cover only g1
    expect_match(val("|linear fold change| among differentially expressed genes detected in both groups, S_vs_N"),
                 "^5 to 5")
})

test_that("F3 PCA and correlation rows are recomputed, and say so", {
    root <- write_run(file.path(withr::local_tempdir(), "rna"), "contrast")
    s <- build_rnaseq_fact_sheet(root, fixture_config(),
                                 pre = list(meta = fixture_meta()),
                                 inputs = list(contrasts = fixture_contrasts()))
    pc <- s[s$claim == "variance explained by PC1 and PC2", ]
    expect_equal(nrow(pc), 1)
    expect_match(pc$source_file, "rna_norm_TMMlogCPM\\.tsv \\(recomputed\\)")
    expect_match(pc$value, "^\\d+\\.\\d{2}% and \\d+\\.\\d{2}%$")
    expect_true(any(s$claim == "mean correlation within group and between groups"))
})

test_that("F3 the collection label drops the contrast prefix", {
    root <- file.path(withr::local_tempdir(), "rna")
    write_run(root, "contrast")
    enr <- file.path(root, "Enrichment")
    dir.create(enr, showWarnings = FALSE)
    fg <- data.frame(pathway = c("p1", "p2"), padj = c(0.01, 0.9), NES = c(2, -1))
    write.csv(fg, file.path(enr, "pathway_S_vs_N_KEGG_spalangia_fgsea.csv"), row.names = FALSE)
    ora <- data.frame(pathway = c("p1", "p2"), padj = c(1, 1), overlap = c(1, 2))
    write.csv(ora, file.path(enr, "pathway_S_vs_N_KEGG_spalangia_ora_up.csv"), row.names = FALSE)

    s <- build_rnaseq_fact_sheet(root, fixture_config(),
                                 pre = list(meta = fixture_meta()),
                                 inputs = list(contrasts = fixture_contrasts()))
    claim <- grep("^gene sets tested", s$claim, value = TRUE)
    expect_equal(claim, "gene sets tested, ranked and over-representation hits (KEGG_spalangia, S_vs_N)")
    expect_match(s$value[s$claim == claim],
                 "2 tested; 1 ranked at adjusted p <= 0.05; 0 over-represented")
})


# =============================================================================
# F4 / F5 — degrade gracefully, and never state a value without a source
# =============================================================================
test_that("F4 an empty results directory yields a sheet, not an error", {
    root <- file.path(withr::local_tempdir(), "rna")
    dir.create(root, recursive = TRUE)
    s <- build_rnaseq_fact_sheet(root, fixture_config())
    expect_s3_class(s, "data.frame")
    # config-derived rows survive because they need no artifact
    expect_true("significance rule" %in% s$claim)
    expect_false(any(grepl("^genes ", s$claim)))
})

test_that("F4 a run without enrichment simply omits those rows", {
    root <- write_run(file.path(withr::local_tempdir(), "rna"), "contrast")
    s <- build_rnaseq_fact_sheet(root, fixture_config(),
                                 pre = list(meta = fixture_meta()),
                                 inputs = list(contrasts = fixture_contrasts()))
    expect_false(any(grepl("^gene sets tested", s$claim)))
    expect_true("libraries analysed" %in% s$claim)
})

test_that("F5 every row states where its value came from", {
    root <- write_run(file.path(withr::local_tempdir(), "rna"), "contrast")
    s <- build_rnaseq_fact_sheet(root, fixture_config(),
                                 pre = list(meta = fixture_meta()),
                                 inputs = list(contrasts = fixture_contrasts()))
    expect_true(nrow(s) > 10)
    expect_true(all(nzchar(s$claim)))
    expect_true(all(nzchar(s$value)))
    expect_true(all(nzchar(s$source_file)))
    expect_false(any(duplicated(s$claim)))
})


# =============================================================================
# F6 — claims that must stay true: per contrast, per source, per run
# =============================================================================

# Two contrasts, and a padj column with NAs where the raw p-value is present --
# the shape DESeq2's independent filtering produces.
write_two_contrast_run <- function(root) {
    ds <- file.path(root, "Datasets")
    dir.create(ds, recursive = TRUE, showWarnings = FALSE)

    counts <- data.frame(gene = paste0("g", 1:4),
                         S_1 = c(100, 10, 5, 40), S_2 = c(120, 12, 6, 44),
                         N_1 = c(20, 11, 30, 42), N_2 = c(24, 9, 34, 38),
                         check.names = FALSE)
    write.table(counts, file.path(ds, "rna_counts_filtered.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)
    write.table(counts, file.path(ds, "rna_norm_TMMlogCPM.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)

    final <- data.frame(
        Gene = paste0("g", 1:4),
        Mean.S = c(110, 11, 5.5, 42), Mean.N = c(22, 10, 32, 40),
        linearFC.S_vs_N = c(5.0, 1.1, -6.0, 1.05),
        pvalue.S_vs_N   = c(1e-8, 0.4, 1e-5, 0.9),
        padj.S_vs_N     = c(1e-6, NA,  1e-4, NA),      # filtered out for 2 genes
        S_vs_N_pass     = c(1, NA, 1, NA),
        linearFC.X_vs_N = c(1.02, 9.0, 1.01, 1.03),
        pvalue.X_vs_N   = c(0.7, 1e-9, 0.8, 0.6),
        padj.X_vs_N     = c(0.9, 1e-7, 0.95, 0.85),
        X_vs_N_pass     = c(NA, 1, NA, NA),
        pass_any_contrast = c(1, 1, 1, NA),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    write.table(final, file.path(ds, "final_results.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)
    root
}

two_contrast_inputs <- function() {
    list(contrasts = data.frame(
        Contrast_name = c("S_vs_N", "X_vs_N"), Factor = "grp",
        Numerator = c("S", "X"), Denominator = "N", stringsAsFactors = FALSE))
}

test_that("F6 every contrast gets its own rows, not just the first", {
    root <- write_two_contrast_run(file.path(withr::local_tempdir(), "rna"))
    s <- build_rnaseq_fact_sheet(root, fixture_config(),
                                 pre = list(meta = fixture_meta()),
                                 inputs = two_contrast_inputs())
    val <- function(claim) s$value[s$claim == claim]

    expect_equal(val("genes at raw p < 0.05, S_vs_N"), "2")
    expect_equal(val("genes at raw p < 0.05, X_vs_N"), "1")
    expect_equal(val("smallest adjusted p-value, S_vs_N"), "1e-06")
    expect_equal(val("smallest adjusted p-value, X_vs_N"), "1e-07")
    expect_false(any(duplicated(s$claim)))
})

test_that("F6 the by-chance comparator counts genes that have a raw p-value", {
    root <- write_two_contrast_run(file.path(withr::local_tempdir(), "rna"))
    s <- build_rnaseq_fact_sheet(root, fixture_config(),
                                 pre = list(meta = fixture_meta()),
                                 inputs = two_contrast_inputs())

    # All 4 genes carry a raw p-value for S_vs_N; only 2 survived independent
    # filtering into padj. Counting the 2 would understate chance by half.
    expect_match(s$source_file[s$claim == "genes expected at raw p < 0.05 by chance, S_vs_N"],
                 "0.05 x 4 genes with a raw p-value", fixed = TRUE)
})

test_that("F6 the fold-change range uses that contrast's own pass flag", {
    root <- write_two_contrast_run(file.path(withr::local_tempdir(), "rna"))
    s <- build_rnaseq_fact_sheet(root, fixture_config(),
                                 pre = list(meta = fixture_meta()),
                                 inputs = two_contrast_inputs())
    val <- function(claim) s$value[s$claim == claim]

    # S_vs_N passes g1 (|FC| 5) and g3 (|FC| 6); X_vs_N passes only g2 (|FC| 9).
    # Under pass_any_contrast both rows would have covered g1, g2 and g3.
    expect_match(val("|linear fold change| among differentially expressed genes detected in both groups, S_vs_N"),
                 "^5 to 6")
    expect_match(val("|linear fold change| among differentially expressed genes detected in both groups, X_vs_N"),
                 "^9 to 9")
})

test_that("F6 a legacy padj_cutoff config is reported, not silently replaced by 0.05", {
    root <- write_run(file.path(withr::local_tempdir(), "rna"), "contrast")
    cfg <- list(modes = list(rna = list(
        de = list(padj_cutoff = 0.1, linear_fc_cutoff = 1.5),   # no p_cutoff
        effects = list(samples = "SampleID", color = "grp")
    )))
    s <- build_rnaseq_fact_sheet(root, cfg, pre = list(meta = fixture_meta()),
                                 inputs = list(contrasts = fixture_contrasts()))

    expect_equal(s$value[s$claim == "significance rule"],
                 "adjusted p <= 0.1 and |linear fold change| >= 1.5")
})

test_that("F6 the fresh Datasets DE counts win over a stale copy at the root", {
    root <- write_run(file.path(withr::local_tempdir(), "rna"), "Name")   # Datasets copy
    stale <- data.frame(contrast = c("S_vs_N", "any"), up = c(99, 0),
                        down = c(99, 0), total = c(198, 198))
    write.table(stale, file.path(root, "de_summary_counts.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)

    d <- .read_de_summary_counts(root, create_legacy_output_dirs(root, create = FALSE))
    expect_equal(d$total, 2)   # not the 198 left behind by a previous run
})

test_that("F6 provenance rows point outside the mode directory", {
    run <- withr::local_tempdir()
    root <- file.path(run, "rna")
    write_run(root, "contrast")
    info <- file.path(run, "execution_info")
    dir.create(info, recursive = TRUE, showWarnings = FALSE)
    writeLines("2026-01-01 00:00:00", file.path(info, "timestamp.txt"))

    s <- build_rnaseq_fact_sheet(root, fixture_config(),
                                 pre = list(meta = fixture_meta()),
                                 inputs = list(contrasts = fixture_contrasts()),
                                 run_dir = run)

    # source_file is relative to the mode directory, and execution_info is that
    # directory's sibling -- "execution_info/..." would name nothing.
    expect_equal(s$source_file[s$claim == "run produced at"],
                 "../execution_info/timestamp.txt")
    expect_true(all(grepl("^\\.\\./execution_info/",
                          s$source_file[grepl("execution_info", s$source_file)])))
})

test_that("F6 the nested local-enrichment layout is found too", {
    root <- file.path(withr::local_tempdir(), "rna")
    write_run(root, "contrast")
    unit <- file.path(root, "Enrichment", "GSEA", "GO_BP", "ranking_by_fc", "S_vs_N")
    dir.create(unit, recursive = TRUE, showWarnings = FALSE)
    write.csv(data.frame(ID = c("p1", "p2", "p3"), p.adjust = c(0.01, 0.02, 0.9)),
              file.path(unit, "results.csv"), row.names = FALSE)

    s <- build_rnaseq_fact_sheet(root, fixture_config(),
                                 pre = list(meta = fixture_meta()),
                                 inputs = list(contrasts = fixture_contrasts()))
    claim <- grep("^gene sets tested and hit", s$claim, value = TRUE)

    expect_length(claim, 1)
    expect_match(claim, "GO_BP, ranking_by_fc, S_vs_N", fixed = TRUE)
    expect_equal(s$value[s$claim == claim], "3 tested; 2 significant")
})

test_that("F6 missing values in the normalised matrix do not silently drop the structure rows", {
    root <- file.path(withr::local_tempdir(), "rna")
    write_run(root, "contrast")
    ds <- file.path(root, "Datasets")
    norm <- read.delim(file.path(ds, "rna_norm_TMMlogCPM.tsv"), check.names = FALSE)
    norm[2, "N_1"] <- NA                       # what preprocessed input can carry through
    write.table(norm, file.path(ds, "rna_norm_TMMlogCPM.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)

    s <- build_rnaseq_fact_sheet(root, fixture_config(),
                                 pre = list(meta = fixture_meta()),
                                 inputs = list(contrasts = fixture_contrasts()))
    pc <- s[s$claim == "variance explained by PC1 and PC2", ]

    expect_equal(nrow(pc), 1)                  # prcomp() would have aborted
    expect_match(pc$source_file, "3 of 4 genes")
    expect_false(any(grepl("NA", s$value[s$claim == "pearson correlation between libraries"])))
})
