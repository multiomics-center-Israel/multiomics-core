# KO-space pathview: the ID transforms in utils/build_feature_ko_map.R, the
# feature -> KO join, the pathway selector and the colour-key placement used by
# generate_per_omic_union_pathview().
#
# All fixtures are synthetic; nothing here touches KEGG or the network.

repo_root <- normalizePath(if (dir.exists("R")) "." else "../..")
source(file.path(repo_root, "utils", "build_feature_ko_map.R"))

# A small eggNOG-mapper annotation, with the emapper banner line, the "#query"
# header, one query without any KO, and one query carrying two.
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


# ---- reading the annotation -------------------------------------------------

test_that("read_eggnog_ko_map strips the ko: prefix and drops '-'", {
    path <- write_fake_eggnog(withr::local_tempfile(fileext = ".tsv"))
    ko <- read_eggnog_ko_map(path)

    expect_named(ko, c("evm.model.ctg1_np1.1", "BRK_g1.t1", "BRK_g2.t1"))
    expect_equal(ko[["evm.model.ctg1_np1.1"]], c("K00001", "K00002"))
    expect_equal(ko[["BRK_g1.t1"]], "K00001")
    expect_length(ko[["BRK_g2.t1"]], 0)
})

test_that("read_eggnog_ko_map rejects a file without the emapper header", {
    path <- withr::local_tempfile(fileext = ".tsv")
    writeLines(c("gene\tko", "g1\tK00001"), path)
    expect_error(read_eggnog_ko_map(path), "#query")
})


# ---- feature id -> annotation key -------------------------------------------

test_that("RNA feature ids map onto the eggNOG query key", {
    ids <- c("evm.TU.ctg1_np1.1", "BRK_g1", "BRK_g2.t1", "tRNA_7")
    expect_identical(
        rna_feature_to_eggnog_key(ids),
        # evm gene -> model naming; a bare BRK gene gains its first transcript;
        # anything already suffixed, or of another shape, is left alone.
        c("evm.model.ctg1_np1.1", "BRK_g1.t1", "BRK_g2.t1", "tRNA_7")
    )
})

test_that("protein groups split into per-member eggNOG keys", {
    expect_identical(
        protein_group_to_eggnog_keys("BRK_g1.t1|Fake_tag"),
        "BRK_g1.t1"
    )
    expect_identical(
        protein_group_to_eggnog_keys("BRK_g2.t1|Fake_tag;BRK_g1.t1|Fake_tag"),
        c("BRK_g2.t1", "BRK_g1.t1")
    )
})

test_that("resolve_features_to_ko emits one long row per (feature, KO)", {
    path <- write_fake_eggnog(withr::local_tempfile(fileext = ".tsv"))
    ko_by_query <- read_eggnog_ko_map(path)

    rna <- resolve_features_to_ko(
        c("evm.TU.ctg1_np1.1", "BRK_g1", "BRK_g2", "BRK_g404"),
        rna_feature_to_eggnog_key, ko_by_query, "transcriptomics"
    )
    # Two KOs for the evm gene, one for BRK_g1; BRK_g2 has no KO and BRK_g404 is
    # not in the annotation at all, so neither contributes a row.
    expect_equal(nrow(rna), 3L)
    expect_setequal(unique(rna$feature_id), c("evm.TU.ctg1_np1.1", "BRK_g1"))
    expect_true(all(rna$omics == "transcriptomics"))

    prot <- resolve_features_to_ko(
        "BRK_g2.t1|Fake_tag;BRK_g1.t1|Fake_tag",
        protein_group_to_eggnog_keys, ko_by_query, "proteomics"
    )
    # The first group member has no KO, so the second one represents the group.
    expect_equal(prot$KO, "K00001")
})

test_that("resolve_features_to_ko returns an empty frame when nothing maps", {
    path <- write_fake_eggnog(withr::local_tempfile(fileext = ".tsv"))
    empty <- resolve_features_to_ko("no_such_gene", rna_feature_to_eggnog_key,
                                    read_eggnog_ko_map(path), "transcriptomics")
    expect_equal(nrow(empty), 0L)
    # The columns must survive an empty result, or rbind() downstream drops them.
    expect_named(empty, c("omics", "feature_id", "KO"))
    expect_type(empty$KO, "character")
})

test_that("read_feature_ids names the column it could not find", {
    path <- withr::local_tempfile(fileext = ".tsv")
    utils::write.table(data.frame(Geneid = c("g1", "g2"), S1 = c(1, 2)),
                       path, sep = "\t", quote = FALSE, row.names = FALSE)

    expect_identical(read_feature_ids(path, "Geneid"), c("g1", "g2"))
    expect_error(read_feature_ids(path, "Protein.Group"), "Protein.Group")
})


# ---- feature log2FC -> KO node ----------------------------------------------

test_that("aggregate_log2fc_by_ko averages features sharing a KO", {
    de <- data.frame(feature_id = c("g1", "g2", "g3"),
                     log2fc = c(1, 3, -2), stringsAsFactors = FALSE)
    ko_map <- data.frame(
        omics = "transcriptomics",
        feature_id = c("g1", "g2", "g3"),
        KO = c("K00001", "K00001", "K00002"),
        stringsAsFactors = FALSE
    )

    fc <- aggregate_log2fc_by_ko(de, ko_map, "transcriptomics")
    # A KO box stands for an ortholog group, so the two paralogues share it and
    # the node carries their mean rather than whichever row came first.
    expect_equal(fc[["K00001"]], 2)
    expect_equal(fc[["K00002"]], -2)
})

test_that("aggregate_log2fc_by_ko ignores non-finite log2FC and other layers", {
    de <- data.frame(feature_id = c("g1", "g2"),
                     log2fc = c(NA_real_, Inf), stringsAsFactors = FALSE)
    ko_map <- data.frame(omics = "transcriptomics",
                         feature_id = c("g1", "g2"), KO = c("K00001", "K00002"),
                         stringsAsFactors = FALSE)

    expect_null(aggregate_log2fc_by_ko(de, ko_map, "transcriptomics"))
    # Nothing in the map belongs to these layers.
    expect_null(aggregate_log2fc_by_ko(de, ko_map, "proteomics"))
    expect_null(aggregate_log2fc_by_ko(NULL, ko_map, "transcriptomics"))
})


# ---- the configured KO map --------------------------------------------------

test_that("load_feature_ko_map reads the configured TSV and drops blank KOs", {
    raw_dir <- withr::local_tempdir()
    dir.create(file.path(raw_dir, "data"), showWarnings = FALSE)
    utils::write.table(
        data.frame(
            omics = c("transcriptomics", "proteomics", "proteomics"),
            feature_id = c("g1", "p1", "p2"),
            KO = c("K00001", "K00002", ""),
            stringsAsFactors = FALSE
        ),
        file.path(raw_dir, "data", "ko.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    config <- list(
        project = list(dir = raw_dir),
        paths = list(raw = "data"),
        modes = list(multiomics = list(enrichment = list(
            pathview = list(species = "ko", ko_map = "ko.tsv")
        )))
    )

    ko_map <- load_feature_ko_map(config)
    expect_named(ko_map, c("omics", "feature_id", "KO"))
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


# ---- which pathways get a map -----------------------------------------------
#
# The selector normalizes first and tests the normalized value, so every
# spelling of one pathway number selects it once, while an identifier that
# normalization does not recognise -- including another organism's accession --
# is still rejected. These pin that contract; the renderer composes the same two
# #201 helpers inline.

test_that("every KEGG spelling of one pathway selects it, foreign prefixes do not", {
    pathways <- c("00010", "map00010", "ko00010", "hsa00010",
                  "hsa00010 Glycolysis / Gluconeogenesis",
                  "mmu00010", "GO:0006915", "My custom gene set")

    keys <- normalize_pathway_join_key(pathways, "hsa")
    keep <- is_kegg_pathway_accession(keys, "hsa")

    expect_identical(keep, c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
    # All five accepted spellings are one map, so one map is rendered.
    expect_identical(unique(keys[keep]), "00010")
})

test_that("in KO space only the species-neutral prefixes are accepted", {
    pathways <- c("map00010", "ko00010", "00010", "hsa00010")

    keys <- normalize_pathway_join_key(pathways, NULL)
    keep <- is_kegg_pathway_accession(keys, NULL)

    # A run with no KEGG code of its own has no organism whose accessions it
    # could legitimately claim, so hsa00010 in the table is foreign here too.
    expect_identical(keep, c(TRUE, TRUE, TRUE, FALSE))
})


# ---- colour key placement ---------------------------------------------------

test_that("pick_key_position avoids the inked corners", {
    skip_if_not_installed("png")

    # A blank template with ink in both right-hand corners. bottomleft is the
    # only quiet candidate, and topleft is never a candidate.
    img <- matrix(1, nrow = 100, ncol = 100)
    img[1:28, 73:100] <- 0
    img[73:100, 73:100] <- 0
    path <- withr::local_tempfile(fileext = ".png")
    png::writePNG(img, path)

    expect_identical(pick_key_position(path), "bottomleft")
})

test_that("pick_key_position falls back to the default corner", {
    # A first run has not downloaded the template yet, and an unreadable file
    # must not take the map down with it.
    expect_identical(pick_key_position(file.path(tempdir(), "not_here.png")),
                     "topright")

    bad <- withr::local_tempfile(fileext = ".png")
    writeLines("not a png", bad)
    expect_identical(pick_key_position(bad), "topright")
})


# ---- the renderer's guard rails ---------------------------------------------

test_that("generate_per_omic_union_pathview is a no-op without a KO map", {
    skip_if_not_installed("pathview")
    config <- list(
        # An organism get_kegg_organism() cannot resolve to a KEGG code.
        global = list(organism = "Unlisted nonmodel species"),
        project = list(dir = tempdir()),
        paths = list(raw = "data")
    )
    # No KEGG code for the organism and no ko_map configured: nothing to draw,
    # and in particular no KEGG request is made.
    expect_null(generate_per_omic_union_pathview(list(), list(), config,
                                                 withr::local_tempdir()))
})

test_that("enrichment.pathview.species accepts only 'ko'", {
    skip_if_not_installed("pathview")
    config <- list(
        global = list(organism = "Unlisted nonmodel species"),
        project = list(dir = tempdir()),
        paths = list(raw = "data"),
        modes = list(multiomics = list(enrichment = list(
            # global.organism is the only place an organism is named, so a code
            # here is ignored rather than quietly overriding it.
            pathview = list(species = "hsa")
        )))
    )
    expect_message(
        expect_null(generate_per_omic_union_pathview(list(), list(), config,
                                                     withr::local_tempdir())),
        "accepts only 'ko'"
    )
})
