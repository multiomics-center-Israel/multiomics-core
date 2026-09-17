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

test_that("a leading Gene: prefix is dropped before the lookup", {
    # Some GFF-derived counts matrices carry "Gene:" where the proteome, and so
    # the eggNOG query column, does not.
    expect_identical(
        rna_feature_to_eggnog_key(c("Gene:EHI_012345", "Gene:evm.TU.ctg1_np1.1",
                                    "Gene:BRK_g1")),
        c("EHI_012345", "evm.model.ctg1_np1.1", "BRK_g1.t1")
    )
})

test_that("the KO map keeps the feature id the DE tables join on", {
    path <- withr::local_tempfile(fileext = ".tsv")
    writeLines(c(
        "## emapper-2.1.13",
        paste("#query", "seed_ortholog", "KEGG_ko", "PFAMs", sep = "\t"),
        paste("EHI_012345", "9999.XP_1", "ko:K00001", "PF00001", sep = "\t")
    ), path)

    out <- resolve_features_to_ko("Gene:EHI_012345", rna_feature_to_eggnog_key,
                                  read_eggnog_ko_map(path), "transcriptomics")

    # The prefix is stripped for the lookup only: the map has to key on the id
    # the DE table actually carries, or the join finds nothing.
    expect_identical(out$feature_id, "Gene:EHI_012345")
    expect_identical(out$KO, "K00001")
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
            pathview = list(ko_map = "ko.tsv")
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
# is still rejected. These pin that contract. Identity is organism-aware even
# though the maps are always drawn in KO space -- the two are separate questions.

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

test_that("pathway identity reads the full organism registry", {
    # get_kegg_organism() knows six exact species names; the pathway identity
    # helpers resolve through get_organism_info() as well. Resolving with the
    # narrow one left an organism KEGG does cover unable to recognise its own
    # pathway accessions.
    expect_null(get_kegg_organism("Saccharomyces cerevisiae"))
    expect_identical(resolve_kegg_org_code("Saccharomyces cerevisiae"), "sce")

    expect_null(get_kegg_organism("Arabidopsis thaliana"))
    expect_identical(resolve_kegg_org_code("Arabidopsis thaliana"), "ath")

    # And an organism neither knows still resolves to nothing, which is what
    # leaves a genuinely code-less run accepting only map/ko prefixes.
    expect_null(resolve_kegg_org_code("Unlisted nonmodel species"))
})

test_that("drawing in KO space does not make the run's own accessions foreign", {
    # The maps are always KEGG reference maps, but an organism that has a code
    # still spells its own pathways with it. Dropping the code from identity
    # would reject the "<accession> <name>" form its gene sets carry, and the
    # run would select nothing to draw.
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
# keyed on the `contrast` column the ORA rows carry, grouped by the canonical
# key the cross-omics module already uses, and the DE table is looked up by that
# same key -- the spellings differ between exports.

ora_rows <- function(contrast, pathway, pvalue) {
    data.frame(pathway = pathway, pvalue = pvalue, contrast = contrast,
               stringsAsFactors = FALSE)
}

test_that("one biological contrast gets one key, however it is spelled", {
    # The documented case: a proteomics ORA export drops the spaces its DE
    # tables keep, and make.names() alone leaves those as two contrasts.
    expect_identical(normalize_contrast_key("1.56ppmvs.0ppm"),
                     normalize_contrast_key("1.56ppm vs. 0ppm"))
    expect_false(identical(make.names("1.56ppmvs.0ppm"),
                           make.names("1.56ppm vs. 0ppm")))

    # The other spellings this pipeline produces agree too: RNA's underscores,
    # metabolomics' dash, and the bare "vs".
    expect_identical(normalize_contrast_key(c("A vs. B", "A_vs_B", "a - b", "A vs B")),
                     rep("avsb", 4L))
    expect_length(unique(normalize_contrast_key(
        c("1.56ppm vs. 0ppm", "1.56ppm_vs_0ppm", "1.56ppm - 0ppm", "1.56ppmvs.0ppm")
    )), 1L)

    # Two genuinely different contrasts stay apart.
    expect_false(identical(normalize_contrast_key("A vs. B"),
                           normalize_contrast_key("A vs. C")))
    # Style is what gets dropped, not content: stripping every non-alphanumeric
    # character made these two doses one key, which is how a contrast would
    # have been rendered under another's fold changes.
    expect_false(identical(normalize_contrast_key("1.56ppm vs. 0ppm"),
                           normalize_contrast_key("15.6ppm vs. 0ppm")))
})

test_that(".kegg_hits_by_contrast keeps each contrast's pathways to itself", {
    tables <- list(
        ora_rows("A vs. B", c("hsa00010", "hsa00020"), c(0.001, 0.30)),
        ora_rows("A vs. C", "hsa00030", 0.002)
    )

    hits <- .kegg_hits_by_contrast(tables, "hsa")

    expect_setequal(names(hits), c("avsb", "avsc"))
    # 00020 is above the cutoff, so the first contrast keeps only 00010.
    expect_identical(hits[["avsb"]]$pathways, "00010")
    expect_identical(hits[["avsc"]]$pathways, "00030")
    # The pathway only the second contrast found must not appear under the first.
    expect_false("00030" %in% hits[["avsb"]]$pathways)
    # The readable spelling survives for headings and filenames.
    expect_identical(hits[["avsb"]]$label, "A vs. B")
})

test_that(".kegg_hits_by_contrast merges the spellings of one contrast", {
    # The same contrast arriving from two exports is one entry, not two.
    tables <- list(
        ora_rows("1.56ppm vs. 0ppm", "hsa00010", 0.001),
        ora_rows("1.56ppmvs.0ppm", "hsa00020", 0.001)
    )

    hits <- .kegg_hits_by_contrast(tables, "hsa")

    expect_length(hits, 1L)
    expect_setequal(hits[[1]]$pathways, c("00010", "00020"))
})

test_that(".kegg_hits_by_contrast drops rows it cannot attribute to a contrast", {
    no_contrast <- data.frame(pathway = "hsa00010", pvalue = 0.001,
                              stringsAsFactors = FALSE)
    blank <- ora_rows("", "hsa00020", 0.001)

    expect_length(.kegg_hits_by_contrast(list(no_contrast), "hsa"), 0L)
    expect_length(.kegg_hits_by_contrast(list(blank), "hsa"), 0L)
})

test_that(".de_table_for_contrast never falls back to the first table", {
    tables <- list("A vs. B" = data.frame(feature_id = "g1", log2fc = 1),
                   "A vs. C" = data.frame(feature_id = "g1", log2fc = -5))

    expect_equal(.de_table_for_contrast(tables, "avsc")$log2fc, -5)
    # A contrast with no DE table means the layer is absent for it, which is a
    # skipped layer -- not a licence to reach for tables[[1]].
    expect_null(.de_table_for_contrast(tables, "avsd"))
    expect_null(.de_table_for_contrast(NULL, "avsb"))
    expect_null(.de_table_for_contrast(list(data.frame(x = 1)), "avsb"))
})

test_that(".de_table_for_contrast matches across the export spellings", {
    # The DE tables keep the spaces the ORA export drops.
    tables <- list("1.56ppm vs. 0ppm" = data.frame(feature_id = "g1", log2fc = 3))
    expect_equal(
        .de_table_for_contrast(tables, normalize_contrast_key("1.56ppmvs.0ppm"))$log2fc,
        3
    )
})

test_that("two distinct contrasts cannot write one output filename", {
    # "A-B" is a comparison and "A B" is a group whose name has a space: two
    # canonical contrasts, and make.names() turns both into "A.B", so one
    # contrast's maps would have overwritten the other's.
    k1 <- normalize_contrast_key("A-B")
    k2 <- normalize_contrast_key("A B")

    expect_false(identical(k1, k2))
    expect_identical(make.names("A-B"), make.names("A B"))   # the collision
    expect_false(identical(.contrast_out_key(k1), .contrast_out_key(k2)))
})

test_that(".contrast_out_key is one-to-one over canonical keys", {
    # The key alphabet is [a-z0-9.], so replacing the dots cannot merge two.
    keys <- normalize_contrast_key(c("1.56ppm vs. 0ppm", "15.6ppm vs. 0ppm",
                                     "A vs. B", "A vs. C", "A-B", "A B"))
    out <- .contrast_out_key(keys)

    expect_length(unique(out), length(unique(keys)))
    expect_false(any(grepl("[^a-z0-9_]", out)))
})

test_that(".de_table_for_contrast treats an ambiguous key as a miss", {
    # Two table names collapsing to one key is not a match to guess at.
    tables <- list("A vs. B" = data.frame(feature_id = "g1", log2fc = 1),
                   "A_vs_B"  = data.frame(feature_id = "g1", log2fc = 9))
    expect_null(.de_table_for_contrast(tables, "avsb"))
})

test_that("a pathway from one contrast cannot pick up another's log2FC", {
    # End to end over the two pure pieces the renderer composes, with the two
    # sides spelled as the real exports spell them.
    tables <- list(
        ora_rows("1.56ppmvs.0ppm", "hsa00010", 0.001),
        ora_rows("15.6ppmvs.0ppm", "hsa00030", 0.001)
    )
    de <- list("1.56ppm vs. 0ppm" = data.frame(feature_id = "g1", log2fc = 2,
                                               stringsAsFactors = FALSE),
               "15.6ppm vs. 0ppm" = data.frame(feature_id = "g1", log2fc = -7,
                                               stringsAsFactors = FALSE))

    hits <- .kegg_hits_by_contrast(tables, "hsa")
    owner <- names(hits)[vapply(hits, function(h) "00030" %in% h$pathways, logical(1))]

    expect_length(owner, 1L)
    expect_equal(.de_table_for_contrast(de, owner)$log2fc, -7)
    # And the other contrast still resolves to its own value, not this one.
    other <- names(hits)[vapply(hits, function(h) "00010" %in% h$pathways, logical(1))]
    expect_equal(.de_table_for_contrast(de, other)$log2fc, 2)
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


# ---- a PDF is written whole or not at all -----------------------------------

test_that(".compile_pathview_pdf leaves no partial file behind on failure", {
    skip_if_not_installed("png")

    out <- withr::local_tempdir()
    pdf_path <- file.path(out, "multi_ora_pathview_union.pdf")

    good <- file.path(out, "ok.png")
    png::writePNG(matrix(1, 8, 8), good)
    bad <- file.path(out, "broken.png")
    writeLines("not a png", bad)

    before <- grDevices::dev.cur()
    # The report links this PDF as a download and cannot tell a truncated one
    # from a whole one, so a failure part-way through must leave nothing --
    # here the first page renders and the second does not.
    expect_message(
        expect_null(.compile_pathview_pdf(c(good, bad), pdf_path)),
        "PDF compilation failed"
    )
    expect_false(file.exists(pdf_path))
    # And the device it opened must not be left behind for the next plot.
    expect_identical(grDevices::dev.cur(), before)
})

test_that(".compile_pathview_pdf captions pages and ignores a mismatched set", {
    skip_if_not_installed("png")

    out <- withr::local_tempdir()
    good <- file.path(out, "ok.png")
    png::writePNG(matrix(1, 8, 8), good)

    # One pathway enriched in two contrasts renders two similar maps, and a
    # downloaded PDF has no filename to tell them apart, so each page is
    # captioned with its contrast.
    labelled <- file.path(out, "labelled.pdf")
    expect_identical(
        .compile_pathview_pdf(c(good, good), labelled, c("A vs. B", "A vs. C")),
        labelled
    )
    expect_gt(file.size(labelled), 0)

    # A label set that does not line up with the pages is dropped rather than
    # recycled onto the wrong maps.
    unlabelled <- file.path(out, "unlabelled.pdf")
    expect_identical(
        .compile_pathview_pdf(c(good, good), unlabelled, "only one label"),
        unlabelled
    )
})

test_that(".compile_pathview_pdf writes the file when every page reads", {
    skip_if_not_installed("png")

    out <- withr::local_tempdir()
    pdf_path <- file.path(out, "multi_ora_pathview_union.pdf")
    good <- file.path(out, "ok.png")
    png::writePNG(matrix(1, 8, 8), good)

    expect_identical(.compile_pathview_pdf(good, pdf_path), pdf_path)
    expect_true(file.exists(pdf_path))
    expect_gt(file.size(pdf_path), 0)

    # Nothing to compile is not a failure to report, just nothing.
    expect_null(.compile_pathview_pdf(character(0), pdf_path))
})


# ---- one run's maps are not another's ---------------------------------------

test_that("clear_multi_ora_pathview_outputs removes only what Multi-ORA owns", {
    out <- withr::local_tempdir()
    pv <- file.path(out, "pathview")
    dir.create(pv)

    mine <- c(file.path(pv, "ko00010.multi_ora_A.vs.B.multi.png"),
              file.path(pv, "ko00020.multi_ora_A.vs.B.png"),
              file.path(out, "multi_ora_pathview_supported.pdf"),
              file.path(out, "multi_ora_pathview_union.pdf"),
              file.path(out, "multi_ora_pathview_union.yaml"))
    # Another renderer's overlays, and the download cache: the blank template
    # and the KGML the report reads pathway titles from.
    theirs <- c(file.path(pv, "ko00010.metab_top.png"),
                file.path(pv, "ko00030.prot_top.png"),
                file.path(pv, "ko00010.png"),
                file.path(pv, "ko00010.xml"))
    for (f in c(mine, theirs)) writeLines("x", f)

    expect_message(clear_multi_ora_pathview_outputs(out), "cleared")

    expect_false(any(file.exists(mine)))
    expect_true(all(file.exists(theirs)))
})

test_that("clear_multi_ora_pathview_outputs is quiet on a fresh directory", {
    out <- withr::local_tempdir()
    expect_length(clear_multi_ora_pathview_outputs(out), 0L)
})

test_that("Multi-ORA clears stale maps before deciding it has too few layers", {
    # The rerun that drops an omics layer produces no Multi-ORA at all, and is
    # exactly the rerun whose previous maps would otherwise stay on the page.
    out <- withr::local_tempdir()
    pv <- file.path(out, "pathview")
    dir.create(pv)
    stale <- file.path(pv, "ko00010.multi_ora_avsb.multi.png")
    writeLines("x", stale)
    writeLines("x", file.path(out, "multi_ora_pathview_union.pdf"))

    expect_null(suppressMessages(
        run_multi_ora(list(transcriptomics = list()), NULL, list(), out)
    ))
    expect_false(file.exists(stale))
    expect_false(file.exists(file.path(out, "multi_ora_pathview_union.pdf")))
})

test_that("the renderer's default top_n matches the config validator's", {
    # The validator fills enrichment.pathview.top_n with 5 when it is absent, so
    # a renderer default of anything else means two different answers to one
    # question depending on which path the config took.
    validated <- suppressMessages(suppressWarnings(validate_multiomics_config(
        list(integration = list(methods = "SNF"))
    )))
    expect_equal(
        validated$enrichment$pathview$top_n,
        eval(formals(generate_per_omic_union_pathview)$top_n)
    )
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

test_that("the renderer is a no-op without a KO map, whatever the organism", {
    skip_if_not_installed("pathview")
    no_code <- list(
        global = list(organism = "Unlisted nonmodel species"),
        project = list(dir = tempdir()), paths = list(raw = "data")
    )
    # A KO map is the only route onto a reference map, so without one there is
    # nothing to draw -- and in particular no KEGG request is made.
    expect_message(
        expect_null(generate_per_omic_union_pathview(list(), list(), no_code,
                                                     withr::local_tempdir())),
        "no feature-to-KO map configured"
    )

    # Having a KEGG code does not provide that route: this renderer never draws
    # in an organism's native gene space, because project feature ids are not
    # KEGG gene ids.
    with_code <- no_code
    with_code$global$organism <- "human"
    expect_message(
        expect_null(generate_per_omic_union_pathview(list(), list(), with_code,
                                                     withr::local_tempdir())),
        "no feature-to-KO map configured"
    )
})

test_that("run_pathview: false stops the renderer before it reads anything", {
    skip_if_not_installed("pathview")
    config <- list(
        global = list(organism = "Unlisted nonmodel species"),
        project = list(dir = tempdir()), paths = list(raw = "data"),
        modes = list(multiomics = list(enrichment = list(
            pathview = list(run_pathview = FALSE, ko_map = "ko.tsv")
        )))
    )
    # Even with a KO map configured, the switch wins -- and it is read before
    # the map, so a disabled run does not go looking for files.
    expect_message(
        expect_null(generate_per_omic_union_pathview(list(), list(), config,
                                                     withr::local_tempdir())),
        "disabled by enrichment.pathview.run_pathview"
    )
})

test_that("the renderer draws only in KO space", {
    # The native-gene branch passed project feature ids to pathview as
    # gene.idtype = "KEGG", which they are not. It is gone, and this is what
    # stops it coming back by accident.
    body_src <- paste(deparse(body(generate_per_omic_union_pathview)), collapse = " ")

    expect_true(grepl('species = "ko"', body_src, fixed = TRUE))
    # No organism code reaches pathview, and the gene-protein bridge that only
    # the native path needed is no longer consulted.
    expect_false(grepl("species = pv_species", body_src, fixed = TRUE))
    expect_false(grepl("gene_protein_mapping", body_src, fixed = TRUE))
})
