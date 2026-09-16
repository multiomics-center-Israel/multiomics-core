# tests/testthat/test-cross-omics-pathway-identity.R
#
# Cross-omics enrichment joins layers that do not agree on where pathway identity
# lives: clusterProfiler-derived layers put the readable description in `pathway`
# and the accession in `ID`, compound ORA does the same with map##### accessions,
# and a table bound from custom GMT collections carries the gene-set name in
# `pathway` with no `ID` at all. pathway_join_key() resolves identity per row and
# normalizes only genuine KEGG accessions, so hsa00010 and map00010 are one
# pathway while GO, PFAM and custom names are left exactly as they were.
#
# All data below is synthetic.

kegg_gene_layer <- function(pvals = c(0.001, 0.02)) {
    data.frame(
        pathway = c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)"),
        ID      = c("hsa00010", "hsa00020"),
        pvalue  = pvals,
        stringsAsFactors = FALSE
    )
}

compound_layer <- function(pvals = c(0.005, 0.04)) {
    data.frame(
        pathway = c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)"),
        ID      = c("map00010", "map00020"),
        pvalue  = pvals,
        stringsAsFactors = FALSE
    )
}


# ---- join key -------------------------------------------------------------

test_that("the key is resolved per row, not per column", {
    # An ID column that only some rows fill: bind_rows() across collections
    # produces exactly this, and a whole-column rule would key the NA rows on NA.
    df <- data.frame(
        pathway = c("Glycolysis / Gluconeogenesis", "Zinc finger domain", "Ribosome"),
        ID      = c("hsa00010", NA, ""),
        stringsAsFactors = FALSE
    )

    expect_equal(
        pathway_join_key(df, kegg_org = "hsa"),
        c("00010", "Zinc finger domain", "Ribosome")
    )
})

test_that("Description is the last resort", {
    df <- data.frame(Description = c("Ribosome", "Spliceosome"),
                     stringsAsFactors = FALSE)
    expect_equal(pathway_join_key(df), c("Ribosome", "Spliceosome"))
})

test_that("a row with nothing to key on yields NA rather than a made-up key", {
    df <- data.frame(pathway = c("Ribosome", NA), ID = c(NA, NA),
                     stringsAsFactors = FALSE)
    expect_equal(pathway_join_key(df), c("Ribosome", NA_character_))
})


# ---- KEGG detection is organism-aware, not shape-only ----------------------

test_that("the KEGG forms this pipeline produces all reduce to the map number", {
    ids <- c("00010", "map00010", "ko00010", "hsa00010")
    expect_equal(normalize_pathway_join_key(ids, kegg_org = "hsa"),
                 rep("00010", 4))
})

test_that("an organism code that is not this run's is left alone", {
    # Shape alone would strip mmu00010 too. Only the run's own organism counts,
    # so a custom term that merely looks like an accession survives intact.
    expect_equal(normalize_pathway_join_key("mmu00010", kegg_org = "hsa"),
                 "mmu00010")
    expect_equal(normalize_pathway_join_key("mmu00010", kegg_org = "mmu"),
                 "00010")
})

test_that("GO, PFAM, InterPro and custom names are byte-identical", {
    ids <- c("GO:0006915", "PF00089", "IPR001", "Glycolysis_custom",
             "HALLMARK_APOPTOSIS", "ko:K00844")
    expect_identical(normalize_pathway_join_key(ids, kegg_org = "hsa"), ids)
    expect_false(any(is_kegg_pathway_accession(ids, kegg_org = "hsa")))
})

test_that("a non-KEGG identifier keeps even its surrounding whitespace", {
    # "Byte-identical" has to mean exactly that. Trimming decides whether a cell
    # counts as empty; it must not become a quiet edit to identifiers this
    # function does not understand.
    ids <- c("  MY_CUSTOM_SET  ", " GO:0006915", "PF00089 ", "\tInterPro_domain")
    expect_identical(normalize_pathway_join_key(ids, kegg_org = "hsa"), ids)
})

test_that("a padded KEGG accession is still recognised and normalized", {
    # Detection runs on a trimmed copy, so padding does not hide a real accession.
    expect_equal(
        normalize_pathway_join_key(c(" map00010 ", "hsa00010\t", "\n00010"),
                                   kegg_org = "hsa"),
        rep("00010", 3)
    )
})

test_that("the join key carries a padded custom value through unaltered", {
    df <- data.frame(pathway = c("  MY_CUSTOM_SET  ", " map00010 "),
                     stringsAsFactors = FALSE)

    expect_identical(pathway_join_key(df, kegg_org = "hsa"),
                     c("  MY_CUSTOM_SET  ", "00010"))
})

test_that("a KEGGREST-style key carrying its name still joins on the accession", {
    # fetch_kegg_via_rest() (R/core/09_enrichment.R) names each gene set
    # "<accession> <readable name>". That is the non-model KEGG fallback -- the
    # path this join most needs to work -- and an end-anchored rule would leave
    # those keys unnormalized and unable to meet map00010.
    ids <- c("hsa00010 Glycolysis / Gluconeogenesis",
             "ko00020 Citrate cycle (TCA cycle)",
             "map00030 Pentose phosphate pathway")

    expect_equal(normalize_pathway_join_key(ids, kegg_org = "hsa"),
                 c("00010", "00020", "00030"))
})

test_that("a name-carrying key joins the bare accession form", {
    tables <- list(
        rna          = data.frame(pathway = "hsa00010 Glycolysis / Gluconeogenesis",
                                  pvalue = 0.01, stringsAsFactors = FALSE),
        metabolomics = data.frame(pathway = "Glycolysis / Gluconeogenesis",
                                  ID = "map00010", pvalue = 0.02,
                                  stringsAsFactors = FALSE)
    )

    merged <- merge_pathway_pvalues(tables, "00010", names(tables),
                                    kegg_org = "hsa")

    expect_equal(nrow(merged), 1)
    expect_false(anyNA(merged$pval_rna))
    expect_false(anyNA(merged$pval_metabolomics))
})

test_that("a leading token that is not an accession is left alone", {
    # Only a genuine accession at the head counts; prose stays prose.
    ids <- c("Glycolysis hsa00010", "HALLMARK_00010 set", "IPR001 domain")
    expect_identical(normalize_pathway_join_key(ids, kegg_org = "hsa"), ids)
})

test_that("an organism only the core registry knows still yields its KEGG code", {
    # get_kegg_organism() knows six species; get_organism_info() knows those plus
    # yeast, Arabidopsis, chicken, pig, cow and Giardia. Giardia is in the core
    # registry precisely because it has a KEGG code and no OrgDb -- the non-model
    # shape this work is for -- so the join key has to reach the wider registry.
    skip_if_not(exists("get_organism_info", mode = "function"),
                "source R/core/11_annotation.R for this test")

    expect_null(get_kegg_organism("Giardia lamblia"))
    expect_equal(resolve_kegg_org_code("Giardia lamblia"), "gla")
    expect_equal(resolve_kegg_org_code("yeast"), "sce")

    # And the six the local map already had keep resolving.
    expect_equal(resolve_kegg_org_code("human"), "hsa")

    # An organism neither registry knows has no code at all, not NA.
    expect_null(resolve_kegg_org_code("Nonexistent organism"))
})

test_that("a Giardia-style key normalizes once the wider registry supplies the code", {
    skip_if_not(exists("get_organism_info", mode = "function"),
                "source R/core/11_annotation.R for this test")

    kegg_org <- resolve_kegg_org_code("Giardia lamblia")

    expect_equal(
        normalize_pathway_join_key(
            c("gla00010", "gla00010 Glycolysis / Gluconeogenesis", "map00010"),
            kegg_org = kegg_org
        ),
        rep("00010", 3)
    )
})

test_that("with no KEGG organism, only the species-neutral forms normalize", {
    expect_equal(
        normalize_pathway_join_key(c("map00010", "ko00010", "00010", "hsa00010"),
                                   kegg_org = NULL),
        c("00010", "00010", "00010", "hsa00010")
    )
})


# ---- the join -------------------------------------------------------------

test_that("KEGG prefix variants of one pathway join across layers", {
    tables <- list(rna = kegg_gene_layer(), metabolomics = compound_layer())

    merged <- merge_pathway_pvalues(tables, c("00010", "00020"),
                                    names(tables), kegg_org = "hsa")

    expect_equal(nrow(merged), 2)
    expect_true(all(c("pval_rna", "pval_metabolomics") %in% names(merged)))
    # Both layers landed on both rows; before the normalization one side keyed on
    # hsa##### and the other on map#####, so every cell of one column was NA.
    expect_false(anyNA(merged$pval_rna))
    expect_false(anyNA(merged$pval_metabolomics))
})

test_that("two custom collections with different names stay separate", {
    # Neither name is KEGG-like, so nothing is rewritten and nothing merges.
    tables <- list(
        rna = data.frame(pathway = c("Glycolysis_setA", "Glycolysis_setB"),
                         pvalue = c(0.01, 0.02), stringsAsFactors = FALSE),
        proteomics = data.frame(pathway = c("Glycolysis_setA", "Glycolysis_setB"),
                                pvalue = c(0.03, 0.04), stringsAsFactors = FALSE)
    )

    merged <- merge_pathway_pvalues(tables, c("Glycolysis_setA", "Glycolysis_setB"),
                                    names(tables), kegg_org = "hsa")

    expect_equal(nrow(merged), 2)
    expect_setequal(merged$norm_id, c("Glycolysis_setA", "Glycolysis_setB"))
})

test_that("contrasts still collapse to one row per key by minimum p-value", {
    # Three contrasts of one pathway in one layer: the pre-existing aggregate()
    # behaviour, which the change to a normalized key must not alter.
    tables <- list(
        rna = data.frame(
            pathway = rep("Glycolysis / Gluconeogenesis", 3),
            ID      = rep("hsa00010", 3),
            pvalue  = c(0.05, 0.001, 0.02),
            stringsAsFactors = FALSE
        ),
        proteomics = data.frame(pathway = "Glycolysis / Gluconeogenesis",
                                ID = "hsa00010", pvalue = 0.01,
                                stringsAsFactors = FALSE)
    )

    merged <- merge_pathway_pvalues(tables, "00010", names(tables),
                                    kegg_org = "hsa")

    expect_equal(nrow(merged), 1)
    expect_equal(merged$pval_rna, 0.001)
})

test_that("the merge cannot multiply rows", {
    # Both layers carry several rows per key -- across contrasts and across KEGG
    # prefix variants. aggregate() reduces each layer to unique keys first, so the
    # result has exactly one row per requested key rather than their product.
    tables <- list(
        rna = data.frame(
            pathway = rep("Glycolysis / Gluconeogenesis", 4),
            ID      = c("hsa00010", "hsa00010", "map00010", "00010"),
            pvalue  = c(0.01, 0.02, 0.03, 0.04),
            stringsAsFactors = FALSE
        ),
        metabolomics = data.frame(
            pathway = rep("Glycolysis / Gluconeogenesis", 3),
            ID      = c("map00010", "map00010", "ko00010"),
            pvalue  = c(0.05, 0.06, 0.07),
            stringsAsFactors = FALSE
        )
    )

    merged <- merge_pathway_pvalues(tables, "00010", names(tables),
                                    kegg_org = "hsa")

    expect_equal(nrow(merged), 1)
    expect_equal(merged$pval_rna, 0.01)
    expect_equal(merged$pval_metabolomics, 0.05)
})


# ---- display --------------------------------------------------------------

test_that("the readable label is attached and identity stays available", {
    tables <- list(rna = kegg_gene_layer(), metabolomics = compound_layer())
    meta <- data.frame(norm_id = c("00010", "00020"), stringsAsFactors = FALSE)

    out <- attach_pathway_display_names(meta, tables, kegg_org = "hsa")

    expect_equal(out$pathway,
                 c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)"))
    expect_equal(out$norm_id, c("00010", "00020"))
    # The readable column reads first in the exported table.
    expect_equal(names(out)[1], "pathway")
})

test_that("a readable name beats an accession-only label for the same key", {
    # The GMT-keyed layer can only offer "map00010" as its label. It is listed
    # first, and must still not become the display name for that pathway.
    tables <- list(
        custom = data.frame(pathway = "map00010", pvalue = 0.01,
                            stringsAsFactors = FALSE),
        rna    = data.frame(pathway = "Glycolysis / Gluconeogenesis",
                            ID = "hsa00010", pvalue = 0.02,
                            stringsAsFactors = FALSE)
    )
    meta <- data.frame(norm_id = "00010", stringsAsFactors = FALSE)

    out <- attach_pathway_display_names(meta, tables, kegg_org = "hsa")

    expect_equal(out$pathway, "Glycolysis / Gluconeogenesis")
})

test_that("pathway_name wins over pathway for display, and neither is mutated", {
    tables <- list(rna = data.frame(
        pathway      = "GO:0006915",
        pathway_name = "apoptotic process",
        pvalue       = 0.01,
        stringsAsFactors = FALSE
    ))
    meta <- data.frame(norm_id = "GO:0006915", stringsAsFactors = FALSE)

    out <- attach_pathway_display_names(meta, tables, kegg_org = "hsa")

    expect_equal(out$pathway, "apoptotic process")
    # The source table is untouched -- no display column written back into it.
    expect_equal(tables$rna$pathway, "GO:0006915")
    expect_equal(colnames(tables$rna),
                 c("pathway", "pathway_name", "pvalue"))
})

test_that("a row with no pathway_name falls back to its own pathway", {
    # This is the shape the compound-ORA report table hits: add_pathway_names()
    # fills pathway_name for the rows it could resolve and leaves the rest. A
    # whole-table column choice renders those rows blank; the label has to be
    # resolved per row, which is why the report calls this helper.
    df <- data.frame(
        pathway      = c("map00010", "map00020", "map00030"),
        pathway_name = c("Glycolysis / Gluconeogenesis", NA, "   "),
        stringsAsFactors = FALSE
    )

    expect_equal(
        pathway_display_label(df),
        c("Glycolysis / Gluconeogenesis", "map00020", "map00030")
    )
})

test_that("a key no layer names falls back to the key itself", {
    tables <- list(rna = kegg_gene_layer())
    meta <- data.frame(norm_id = "99999", stringsAsFactors = FALSE)

    out <- attach_pathway_display_names(meta, tables, kegg_org = "hsa")

    expect_equal(out$pathway, "99999")
})

test_that("colliding display labels are separated by their key", {
    # Both figures position rows by label. Two keys sharing a name would stack on
    # one axis slot in the dot plot and make pheatmap error on duplicate rownames.
    out <- disambiguate_pathway_labels(
        c("Glycolysis", "Glycolysis", "Citrate cycle"),
        c("00010", "00051", "00020")
    )

    expect_equal(out, c("Glycolysis (00010)", "Glycolysis (00051)", "Citrate cycle"))
    expect_false(anyDuplicated(out) > 0)
})

test_that("labels that do not collide are untouched", {
    labs <- c("Glycolysis", "Citrate cycle")
    expect_identical(disambiguate_pathway_labels(labs, c("00010", "00020")), labs)
})

test_that("truncation happens before disambiguation, so the key survives", {
    # The key suffix sits at the end. Disambiguating first and truncating second
    # cuts it back off and rebuilds the collision -- which then reaches
    # factor(levels = ) in the dot plot, and duplicated levels are an error.
    long <- paste0(strrep("Glycolysis and gluconeogenesis ", 3))
    labs <- truncate_pathway_label(c(long, long), 45)
    out <- disambiguate_pathway_labels(labs, c("00010", "00051"))

    expect_false(anyDuplicated(out) > 0)
    expect_true(all(grepl("(00010)", out[1], fixed = TRUE),
                    grepl("(00051)", out[2], fixed = TRUE)))
})

test_that("truncate_pathway_label lands on the stated width and leaves short labels alone", {
    expect_equal(nchar(truncate_pathway_label(strrep("a", 80), 45)), 45)
    expect_identical(truncate_pathway_label("Glycolysis", 45), "Glycolysis")
})

test_that("labels still collide-proof with no keys to fall back on", {
    out <- disambiguate_pathway_labels(c("Glycolysis", "Glycolysis"), NULL)
    expect_false(anyDuplicated(out) > 0)
})

test_that("competing readable names resolve deterministically without duplicating rows", {
    # KEGG names differ between sources -- clusterProfiler appends the organism.
    # First readable name wins; the row count must not change either way.
    tables <- list(
        rna          = data.frame(pathway = "Glycolysis - Homo sapiens (human)",
                                  ID = "hsa00010", pvalue = 0.01,
                                  stringsAsFactors = FALSE),
        metabolomics = data.frame(pathway = "Glycolysis / Gluconeogenesis",
                                  ID = "map00010", pvalue = 0.02,
                                  stringsAsFactors = FALSE)
    )
    meta <- data.frame(norm_id = "00010", stringsAsFactors = FALSE)

    out <- attach_pathway_display_names(meta, tables, kegg_org = "hsa")

    expect_equal(nrow(out), 1)
    expect_equal(out$pathway, "Glycolysis - Homo sapiens (human)")
})
