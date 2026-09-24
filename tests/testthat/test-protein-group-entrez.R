# tests/testthat/test-protein-group-entrez.R
#
# Protein groups ("P1;P2") used to be passed whole to the UniProt -> Entrez
# lookup, where no group string is ever a key, so every multi-protein group fell
# out of the KEGG enrichment, the loadings enrichment and the pathway maps.
#
# A group is one measurement and now gets one representative gene: its first
# annotated gene symbol when that maps, otherwise the first of its accessions
# that resolves. A single-accession protein is mapped exactly as before, with an
# isoform tried as its canonical accession only when it has no entry of its own.
#
# Synthetic accessions, symbols and a fake lookup throughout: the org.Db itself
# is not what is under test.

fake_lookup <- function(table) {
    calls <- new.env()
    calls$keys <- list()
    fn <- function(keys) {
        calls$keys[[length(calls$keys) + 1]] <- keys
        stats::setNames(unname(table[keys]), keys)
    }
    list(fn = fn, calls = calls)
}

entrez_of <- function(res) stats::setNames(res$ENTREZID, res$feature_id)

test_that("a group is split into its members, leading protein first", {
    m <- protein_group_members(c("P1;P2; P3", "P4", NA, "", "P1;P2; P3"))
    expect_identical(m$feature_id, c(rep("P1;P2; P3", 3), "P4"))
    expect_identical(m$accession, c("P1", "P2", "P3", "P4"))
    expect_identical(m$position, c(1L, 2L, 3L, 1L))
})

test_that("empty members from stray separators are dropped", {
    m <- protein_group_members("P1;;P2;")
    expect_identical(m$accession, c("P1", "P2"))
    expect_identical(m$position, c(1L, 3L))
})

test_that("the canonical accession drops only a trailing isoform number", {
    expect_identical(canonical_uniprot_accession(c("P12345-2", "P12345", "A0A0B1C2D3", "cRAP-P1")),
                     c("P12345", "P12345", "A0A0B1C2D3", "cRAP-P1"))
})

test_that("a group whose members come from different genes gets one gene, not all", {
    lk <- fake_lookup(c(P1 = "101", P2 = "102"))
    res <- map_protein_groups_to_entrez("P1;P2", lk$fn)
    expect_identical(res$ENTREZID, "101")
    expect_identical(res$source, "accession")
    expect_identical(res$matched_key, "P1")
})

test_that("an unmapped leading accession falls back to the next one that resolves", {
    lk <- fake_lookup(c(P1 = NA_character_, P2 = "102", P3 = "103"))
    res <- map_protein_groups_to_entrez("P1;P2;P3", lk$fn)
    expect_identical(res$ENTREZID, "102")
    expect_identical(res$matched_key, "P2")
})

test_that("an isoform without its own entry maps through its canonical accession", {
    lk <- fake_lookup(c(P9 = "109", P8 = "108", `P8-2` = "208"))
    res <- map_protein_groups_to_entrez(c("P9-3", "P8-2"), lk$fn)
    # P8-2 has an entry of its own and keeps it; P9-3 falls back to P9.
    expect_identical(res$feature_id, c("P8-2", "P9-3"))
    expect_identical(res$ENTREZID, c("208", "109"))
    expect_identical(res$matched_key, c("P8-2", "P9"))
    expect_identical(res$source, c("accession", "canonical_accession"))
})

test_that("a group takes its first gene symbol when it maps, over its accessions", {
    acc <- fake_lookup(c(P1 = "101", P2 = "102"))
    sym <- fake_lookup(c(GENEB = "202"))
    res <- map_protein_groups_to_entrez("P1;P2", acc$fn,
                                        symbols = c(`P1;P2` = "GENEB"),
                                        symbol_lookup = sym$fn)
    expect_identical(res$ENTREZID, "202")
    expect_identical(res$source, "gene_symbol")
    expect_identical(res$matched_key, "GENEB")
})

test_that("a group falls back to its accessions when its symbol is missing or unmapped", {
    acc <- fake_lookup(c(P1 = "101", P3 = "103"))
    sym <- fake_lookup(c(KNOWN = "900"))
    res <- map_protein_groups_to_entrez(
        c("P1;P2", "P3;P4", "P5;P1"), acc$fn,
        symbols = c(`P1;P2` = "UNKNOWN", `P3;P4` = NA, `P5;P1` = ""),
        symbol_lookup = sym$fn)
    expect_identical(entrez_of(res), c(`P1;P2` = "101", `P3;P4` = "103", `P5;P1` = "101"))
    expect_true(all(res$source == "accession"))
})

test_that("a symbol lookup that finds nothing leaves the accession fallback working", {
    acc <- fake_lookup(c(P1 = "101"))
    nothing <- function(keys) stats::setNames(rep(NA_character_, length(keys)), keys)
    res <- map_protein_groups_to_entrez("P1;P2", acc$fn,
                                        symbols = c(`P1;P2` = "GENEA"),
                                        symbol_lookup = nothing)
    expect_identical(res$ENTREZID, "101")
})

test_that("a single-accession protein maps exactly as before, whatever its Genes says", {
    # Regression: the gene annotation is read for protein groups only. P1's
    # Genes entry names a different gene; its UniProt mapping must still win.
    acc <- fake_lookup(c(P1 = "101", `P2-2` = NA_character_, P2 = "102"))
    sym <- fake_lookup(c(OTHER = "999"))
    res <- map_protein_groups_to_entrez(c("P1", "P2-2"), acc$fn,
                                        symbols = c(P1 = "OTHER", `P2-2` = "OTHER"),
                                        symbol_lookup = sym$fn)
    expect_identical(entrez_of(res), c(P1 = "101", `P2-2` = "102"))
    expect_identical(res$source, c("accession", "canonical_accession"))
    expect_length(sym$calls$keys, 0)
})

test_that("every feature gets at most one row, however its groups overlap", {
    acc <- fake_lookup(c(P1 = "101", P2 = "102", P3 = "103", `P1-2` = NA_character_))
    sym <- fake_lookup(c(GENEA = "101", GENEB = "102"))
    ids <- c("P1;P2", "P2;P1", "P2", "P1-2;P3", "P1", "P3;P1;P2")
    res <- map_protein_groups_to_entrez(ids, acc$fn,
                                        symbols = c(`P1;P2` = "GENEA", `P2;P1` = "GENEB"),
                                        symbol_lookup = sym$fn)
    expect_false(anyDuplicated(res$feature_id) > 0)
    expect_setequal(res$feature_id, ids)
    # Several features may share a gene -- collapsing by gene is the caller's job.
    expect_identical(entrez_of(res)[c("P1;P2", "P1", "P1-2;P3")],
                     c(`P1;P2` = "101", P1 = "101", `P1-2;P3` = "101"))
})

test_that("a feature maps the same whether alone or among other features", {
    acc <- fake_lookup(c(P1 = "101", P2 = "102", P3 = "103"))
    sym <- fake_lookup(c(GENEB = "102"))
    symbols <- c(`P1;P2` = "GENEB", `P3;P1` = NA)
    ids <- c("P1", "P2", "P1;P2", "P3;P1", "P2;P3")
    together <- map_protein_groups_to_entrez(ids, acc$fn, symbols, sym$fn)
    for (id in ids) {
        alone <- map_protein_groups_to_entrez(id, acc$fn, symbols, sym$fn)
        expected <- together[together$feature_id == id, , drop = FALSE]
        rownames(expected) <- NULL
        expect_identical(alone, expected, info = id)
    }
})

test_that("the result does not depend on input order", {
    acc <- fake_lookup(c(P1 = "101", P2 = "102", P3 = "103"))
    sym <- fake_lookup(c(GENEC = "103"))
    symbols <- c(`P3;P1` = "GENEC")
    a <- map_protein_groups_to_entrez(c("P3;P1", "P1;P2", "P2"), acc$fn, symbols, sym$fn)
    b <- map_protein_groups_to_entrez(c("P2", "P1;P2", "P3;P1"), acc$fn, symbols, sym$fn)
    expect_identical(a, b)
})

test_that("a group none of whose members map is left out", {
    lk <- fake_lookup(c(P1 = NA_character_, P2 = NA_character_, P5 = "105"))
    res <- map_protein_groups_to_entrez(c("P1;P2", "P5"), lk$fn)
    expect_identical(res$feature_id, "P5")
})

test_that("nothing mapped gives an empty frame with the expected columns", {
    lk <- fake_lookup(c(P1 = NA_character_))
    res <- map_protein_groups_to_entrez("P1", lk$fn)
    expect_identical(nrow(res), 0L)
    expect_identical(names(res), c("feature_id", "ENTREZID", "source", "matched_key"))
    expect_identical(nrow(map_protein_groups_to_entrez(character(0), lk$fn)), 0L)
})

test_that("each key, canonical forms included, is looked up once in one call", {
    acc <- fake_lookup(c(P1 = "101", P2 = "102"))
    sym <- fake_lookup(c(GENEA = NA_character_))
    map_protein_groups_to_entrez(c("P1;P2", "P2;P1", "P1-2"), acc$fn,
                                 symbols = c(`P1;P2` = "GENEA", `P2;P1` = "GENEA"),
                                 symbol_lookup = sym$fn)
    expect_length(acc$calls$keys, 1)
    expect_setequal(acc$calls$keys[[1]], c("P1", "P2", "P1-2"))
    expect_false(anyDuplicated(acc$calls$keys[[1]]) > 0)
    expect_length(sym$calls$keys, 1)
    expect_identical(sym$calls$keys[[1]], "GENEA")
})

test_that("gene symbols are joined to features by value, first gene only", {
    pre <- list(
        expr_work = matrix(0, 3, 1, dimnames = list(c("P1;P2", "P3", "P4;P5"), "S1")),
        row_data = data.frame(Protein.Group = c("P4;P5", "P1;P2", "P3"),
                              Genes = c("", "GENEA; GENEB", "GENE3"),
                              stringsAsFactors = FALSE))
    sy <- protein_group_gene_symbols(pre)
    expect_identical(sy[["P1;P2"]], "GENEA")
    expect_identical(sy[["P3"]], "GENE3")
    expect_true(is.na(sy[["P4;P5"]]))
})

test_that("no gene symbols are returned when the annotation cannot be aligned", {
    expr <- matrix(0, 2, 1, dimnames = list(c("P1;P2", "P3"), "S1"))
    no_id_col <- list(expr_work = expr,
                      row_data = data.frame(Genes = c("GENEA", "GENE3"),
                                            stringsAsFactors = FALSE))
    expect_null(suppressMessages(protein_group_gene_symbols(no_id_col)))

    no_gene_col <- list(expr_work = expr,
                        row_data = data.frame(Protein.Group = c("P1;P2", "P3"),
                                              stringsAsFactors = FALSE))
    expect_null(protein_group_gene_symbols(no_gene_col))
    expect_null(protein_group_gene_symbols(list(expr_work = expr)))
})
