# KO-space pathview: ID transforms (utils/build_spalangia_ko_map.R) and the
# feature -> KO join used by generate_per_omic_union_pathview().
# All fixtures are synthetic; nothing here touches KEGG or the network.

repo_root <- normalizePath(if (dir.exists("R")) "." else "../..")
source(file.path(repo_root, "utils", "build_spalangia_ko_map.R"))

# A three-row eggNOG-mapper annotation, with the emapper banner line, the
# "#query" header, one query without any KO, and one query carrying two.
write_fake_eggnog <- function(path) {
    writeLines(c(
        "## emapper-2.1.13",
        paste("#query", "seed_ortholog", "KEGG_ko", "PFAMs", sep = "\t"),
        paste("evm.model.ctg1_np1.1", "9999.XP_1", "ko:K00001,ko:K00002", "PF00001", sep = "\t"),
        paste("BRK_g1.t1", "9999.XP_2", "ko:K00001", "PF00002", sep = "\t"),
        paste("BRK_g2.t1", "9999.XP_3", "-", "-", sep = "\t")
    ), path)
    path
}

test_that("read_eggnog_ko_map strips the ko: prefix and drops '-'", {
    f <- write_fake_eggnog(withr::local_tempfile(fileext = ".txt"))
    ko <- read_eggnog_ko_map(f)

    expect_named(ko, c("evm.model.ctg1_np1.1", "BRK_g1.t1", "BRK_g2.t1"))
    expect_equal(ko[["evm.model.ctg1_np1.1"]], c("K00001", "K00002"))
    expect_equal(ko[["BRK_g1.t1"]], "K00001")
    expect_length(ko[["BRK_g2.t1"]], 0)
})

test_that("read_eggnog_ko_map rejects a file without the emapper header", {
    f <- withr::local_tempfile(fileext = ".txt")
    writeLines(c("query\tKEGG_ko", "g1\tko:K00001"), f)
    expect_error(read_eggnog_ko_map(f), "#query")
})

test_that("RNA feature ids map onto the eggNOG query key", {
    ids <- c("evm.TU.ctg1_np1.1", "BRK_g1", "BRK_g2.t1", "tRNA_7")
    expect_equal(
        rna_feature_to_eggnog_key(ids),
        # evm gene -> model naming; a bare BRK gene gains its first transcript;
        # anything already transcript-shaped is left alone.
        c("evm.model.ctg1_np1.1", "BRK_g1.t1", "BRK_g2.t1", "tRNA_7")
    )
})

test_that("protein groups split into per-member eggNOG keys", {
    expect_equal(
        protein_group_to_eggnog_keys("BRK_g1.t1|Fake_species"),
        "BRK_g1.t1"
    )
    expect_equal(
        protein_group_to_eggnog_keys("BRK_g2.t1|Fake_species;BRK_g1.t1|Fake_species"),
        c("BRK_g2.t1", "BRK_g1.t1")
    )
})

test_that("resolve_features_to_ko emits one long row per (feature, KO)", {
    ko <- read_eggnog_ko_map(write_fake_eggnog(withr::local_tempfile(fileext = ".txt")))

    rna <- resolve_features_to_ko(
        c("evm.TU.ctg1_np1.1", "BRK_g1", "BRK_g2", "BRK_g404"),
        rna_feature_to_eggnog_key, ko, "transcriptomics"
    )
    expect_equal(names(rna), c("omics", "feature_id", "KO"))
    # Two KOs for the evm gene, one for BRK_g1; BRK_g2 has no KO and BRK_g404 is
    # absent from the annotation, so neither contributes a row.
    expect_equal(nrow(rna), 3L)
    expect_setequal(unique(rna$feature_id), c("evm.TU.ctg1_np1.1", "BRK_g1"))
    expect_true(all(rna$omics == "transcriptomics"))

    # The first group member that carries a KO represents the group.
    prot <- resolve_features_to_ko(
        "BRK_g2.t1|Fake_species;BRK_g1.t1|Fake_species",
        protein_group_to_eggnog_keys, ko, "proteomics"
    )
    expect_equal(prot$KO, "K00001")
})

test_that("resolve_features_to_ko returns an empty frame when nothing maps", {
    ko <- read_eggnog_ko_map(write_fake_eggnog(withr::local_tempfile(fileext = ".txt")))
    out <- resolve_features_to_ko("no_such_feature", rna_feature_to_eggnog_key,
                                  ko, "transcriptomics")
    expect_equal(nrow(out), 0L)
    expect_equal(names(out), c("omics", "feature_id", "KO"))
})

test_that("aggregate_log2fc_by_ko averages features sharing a KO", {
    ko_map <- data.frame(
        omics = c(rep("transcriptomics", 3), "proteomics"),
        feature_id = c("g1", "g2", "g3", "p1"),
        KO = c("K00001", "K00001", "K00002", "K00003"),
        stringsAsFactors = FALSE
    )
    de <- data.frame(
        feature_id = c("g1", "g2", "g3", "g4"),
        log2fc = c(1, 3, -2, 10),
        stringsAsFactors = FALSE
    )

    fc <- aggregate_log2fc_by_ko(de, ko_map, "transcriptomics")
    expect_equal(fc[["K00001"]], 2)   # mean of 1 and 3
    expect_equal(fc[["K00002"]], -2)
    expect_false("K00003" %in% names(fc))  # other layer
    expect_length(fc, 2L)              # g4 has no KO
})

test_that("aggregate_log2fc_by_ko ignores non-finite log2FC and blank KOs", {
    ko_map <- data.frame(
        omics = rep("proteomics", 3),
        feature_id = c("p1", "p2", "p3"),
        KO = c("K00001", "", "K00002"),
        stringsAsFactors = FALSE
    )
    de <- data.frame(
        feature_id = c("p1", "p2", "p3"),
        log2fc = c(NA_real_, 1, Inf),
        stringsAsFactors = FALSE
    )
    expect_null(aggregate_log2fc_by_ko(de, ko_map, "proteomics"))
    expect_null(aggregate_log2fc_by_ko(de, ko_map, "metabolomics"))
    expect_null(aggregate_log2fc_by_ko(NULL, ko_map, "proteomics"))
})

test_that("load_feature_ko_map reads the configured TSV and drops blank KOs", {
    raw_dir <- withr::local_tempdir()
    dir.create(file.path(raw_dir, "data"), showWarnings = FALSE)
    ko_rel <- "ko.tsv"
    write.table(
        data.frame(
            omics = c("transcriptomics", "proteomics", "proteomics"),
            feature_id = c("g1", "p1", "p2"),
            KO = c("K00001", "K00002", ""),
            stringsAsFactors = FALSE
        ),
        file.path(raw_dir, "data", ko_rel),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    config <- list(
        project = list(dir = raw_dir),
        paths = list(raw = "data"),
        modes = list(multiomics = list(enrichment = list(
            pathview = list(species = "ko", ko_map = ko_rel)
        )))
    )

    ko_map <- load_feature_ko_map(config)
    expect_equal(names(ko_map), c("omics", "feature_id", "KO"))
    expect_equal(nrow(ko_map), 2L)
    expect_setequal(ko_map$KO, c("K00001", "K00002"))
})

test_that("load_feature_ko_map is a no-op without the config key or file", {
    config <- list(project = list(dir = tempdir()), paths = list(raw = "data"))
    expect_null(load_feature_ko_map(config))

    config$modes <- list(multiomics = list(enrichment = list(
        pathview = list(ko_map = "definitely_not_here.tsv")
    )))
    expect_message(expect_null(load_feature_ko_map(config)), "ko_map not found")
})

test_that("generate_per_omic_union_pathview is a no-op without a KO map", {
    skip_if_not_installed("pathview")
    config <- list(
        global = list(organism = "Spalangia cameroni"),
        project = list(dir = tempdir()),
        paths = list(raw = "data")
    )
    # No KEGG code for the organism and no ko_map configured: nothing to draw,
    # and in particular no KEGG request is made.
    expect_null(generate_per_omic_union_pathview(list(), list(), config,
                                                 withr::local_tempdir()))
})
