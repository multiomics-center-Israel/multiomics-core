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

test_that("a run with no KEGG code accepts only the species-neutral prefixes", {
    pathways <- c("map00010", "ko00010", "00010", "hsa00010")

    keys <- normalize_pathway_join_key(pathways, NULL)
    keep <- is_kegg_pathway_accession(keys, NULL)

    # Such a run has no organism whose accessions it could legitimately claim,
    # so hsa00010 in its table is foreign.
    expect_identical(keep, c(TRUE, TRUE, TRUE, FALSE))
})

test_that("the render space reads the same organism registry as pathway identity", {
    # get_kegg_organism() knows six exact species names; the pathway identity
    # helpers resolve through get_organism_info() as well. Resolving the render
    # space with the narrow one pushed an organism KEGG does cover into KO mode,
    # where it then needed a ko_map to draw anything.
    expect_null(get_kegg_organism("Saccharomyces cerevisiae"))
    expect_identical(resolve_kegg_org_code("Saccharomyces cerevisiae"), "sce")

    expect_null(get_kegg_organism("Arabidopsis thaliana"))
    expect_identical(resolve_kegg_org_code("Arabidopsis thaliana"), "ath")

    # And an organism neither knows still resolves to nothing, which is what
    # puts a genuinely code-less run into KO mode.
    expect_null(resolve_kegg_org_code("Unlisted nonmodel species"))
})

test_that("forcing KO rendering does not make the run's own accessions foreign", {
    # Rendering space and identity space are separate questions: species = "ko"
    # says draw the reference artwork, not "this run has no organism". Dropping
    # the code here would reject the "<accession> <name>" form the run's own
    # gene sets carry, and the forced mode would select nothing.
    labelled <- "hsa00010 Glycolysis / Gluconeogenesis"

    expect_true(is_kegg_pathway_accession(
        normalize_pathway_join_key(labelled, "hsa"), "hsa"))
    expect_false(is_kegg_pathway_accession(
        normalize_pathway_join_key(labelled, NULL), NULL))
})


# ---- one contrast's pathways, one contrast's fold changes -------------------
#
# The ORA inputs span every contrast, so selecting pathways globally and then
# reading log2FC from whichever DE table happened to come first rendered a
# pathway found in one contrast using another contrast's values. Selection is
# keyed on the `contrast` column the ORA rows already carry, and the DE table is
# then looked up by that name.

ora_rows <- function(contrast, pathway, pvalue) {
    data.frame(pathway = pathway, pvalue = pvalue, contrast = contrast,
               stringsAsFactors = FALSE)
}

test_that(".kegg_hits_by_contrast keeps each contrast's pathways to itself", {
    tables <- list(
        ora_rows("A", c("hsa00010", "hsa00020"), c(0.001, 0.30)),
        ora_rows("B", "hsa00030", 0.002)
    )

    hits <- .kegg_hits_by_contrast(tables, "hsa")

    expect_setequal(names(hits), c("A", "B"))
    # 00020 is above the cutoff, so A keeps only 00010.
    expect_identical(hits[["A"]], "00010")
    expect_identical(hits[["B"]], "00030")
    # The pathway only B found must not appear under A.
    expect_false("00030" %in% hits[["A"]])
})

test_that(".kegg_hits_by_contrast drops rows it cannot attribute to a contrast", {
    no_contrast <- data.frame(pathway = "hsa00010", pvalue = 0.001,
                              stringsAsFactors = FALSE)
    blank <- ora_rows("", "hsa00020", 0.001)

    expect_length(.kegg_hits_by_contrast(list(no_contrast), "hsa"), 0L)
    expect_length(.kegg_hits_by_contrast(list(blank), "hsa"), 0L)
})

test_that(".de_table_for_contrast never falls back to the first table", {
    tables <- list(A = data.frame(feature_id = "g1", log2fc = 1),
                   B = data.frame(feature_id = "g1", log2fc = -5))

    # The value that identifies the table is the one under that contrast's name.
    expect_equal(.de_table_for_contrast(tables, "B")$log2fc, -5)
    # A contrast with no DE table means the layer is absent for it, which is a
    # skipped layer -- not a licence to reach for tables[[1]].
    expect_null(.de_table_for_contrast(tables, "C"))
    expect_null(.de_table_for_contrast(NULL, "A"))
    expect_null(.de_table_for_contrast(list(data.frame(x = 1)), "A"))
})

test_that("a pathway from one contrast cannot pick up another's log2FC", {
    # End to end over the two pure pieces the renderer composes: B's pathway
    # resolves to B's table, and the value it would carry is B's.
    tables <- list(
        ora_rows("A", "hsa00010", 0.001),
        ora_rows("B", "hsa00030", 0.001)
    )
    de <- list("A" = data.frame(feature_id = "g1", log2fc = 2,
                                stringsAsFactors = FALSE),
               "B" = data.frame(feature_id = "g1", log2fc = -7,
                                stringsAsFactors = FALSE))

    hits <- .kegg_hits_by_contrast(tables, "hsa")
    owner <- names(hits)[vapply(hits, function(k) "00030" %in% k, logical(1))]

    expect_identical(owner, "B")
    expect_equal(.de_table_for_contrast(de, owner)$log2fc, -7)
})

test_that("contrast names that differ only syntactically still match", {
    # extract_de_tables() and the ORA rows are named from the same contrast
    # strings, but one of them having been through make.names() must not read as
    # a different contrast.
    tables <- list("cond A vs B" = data.frame(feature_id = "g1", log2fc = 3))
    expect_equal(.de_table_for_contrast(tables, "cond.A.vs.B")$log2fc, 3)
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


# ---- the union artifact is not the ">= 2 omics" artifact --------------------

test_that("the union renderer writes its own PDF, not the supported one", {
    # "supported" is a claim the report repeats: enriched in two or more omics
    # layers. This renderer unions single-layer hits on purpose, so sharing that
    # filename presented one layer's evidence under the other's promise. The
    # separation is the fix, and this is what stops it being undone.
    union_body <- paste(deparse(body(generate_per_omic_union_pathview)), collapse = " ")
    supported_body <- paste(deparse(body(generate_multi_ora_pathview)), collapse = " ")

    expect_true(grepl("multi_ora_pathview_union.pdf", union_body, fixed = TRUE))
    expect_false(grepl("multi_ora_pathview_supported.pdf", union_body, fixed = TRUE))
    expect_true(grepl("multi_ora_pathview_supported.pdf", supported_body, fixed = TRUE))
    expect_false(grepl("multi_ora_pathview_union.pdf", supported_body, fixed = TRUE))
})


# ---- the renderer's guard rails ---------------------------------------------

test_that("generate_per_omic_union_pathview is a no-op without a KO map", {
    skip_if_not_installed("pathview")
    config <- list(
        # An organism neither KEGG registry resolves to a code.
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
