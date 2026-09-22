# tests/testthat/test-protein-group-entrez.R
#
# Protein groups ("P1;P2") used to be passed whole to the UniProt -> Entrez
# lookup, where no group string is ever a key, so every multi-protein group fell
# out of the KEGG enrichment, the loadings enrichment and the pathway maps. The
# lookup now goes member by member: a group maps to every gene its members come
# from, an isoform with no entry of its own is tried as its canonical accession,
# and a gene measured as its own protein is not claimed again through another
# group's secondary member.
#
# Synthetic accessions and a fake lookup throughout: the org.Db itself is not
# what is under test.

fake_lookup <- function(table) {
    calls <- new.env()
    calls$keys <- list()
    fn <- function(keys) {
        calls$keys[[length(calls$keys) + 1]] <- keys
        stats::setNames(unname(table[keys]), keys)
    }
    list(fn = fn, calls = calls)
}

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

test_that("a group whose members come from different genes maps to all of them", {
    lk <- fake_lookup(c(P1 = "101", P2 = "102"))
    res <- map_protein_groups_to_entrez("P1;P2", lk$fn)
    expect_identical(res$ENTREZID, c("101", "102"))
    expect_identical(res$member_rank, c(1L, 2L))
    expect_identical(res$matched_accession, c("P1", "P2"))
})

test_that("two members of one gene count once", {
    lk <- fake_lookup(c(P1 = "101", P2 = "101"))
    res <- map_protein_groups_to_entrez("P1;P2", lk$fn)
    expect_identical(res$ENTREZID, "101")
    expect_identical(res$matched_accession, "P1")
})

test_that("an unmapped leading member leaves the next one first", {
    lk <- fake_lookup(c(P1 = NA_character_, P2 = "102", P3 = "103"))
    res <- map_protein_groups_to_entrez("P1;P2;P3", lk$fn)
    expect_identical(res$ENTREZID, c("102", "103"))
    expect_identical(res$member_rank, c(1L, 2L))
})

test_that("an isoform without its own entry maps through its canonical accession", {
    lk <- fake_lookup(c(P9 = "109", P8 = "108", `P8-2` = "208"))
    res <- map_protein_groups_to_entrez(c("P9-3", "P8-2"), lk$fn)
    res <- res[order(res$feature_id), ]
    # P8-2 has an entry of its own and keeps it; P9-3 falls back to P9.
    expect_identical(res$feature_id, c("P8-2", "P9-3"))
    expect_identical(res$ENTREZID, c("208", "109"))
    expect_identical(res$matched_accession, c("P8-2", "P9"))
    expect_identical(res$via_canonical, c(FALSE, TRUE))
})

test_that("a gene measured on its own is not claimed again through a secondary member", {
    lk <- fake_lookup(c(P1 = "101", P2 = "102", P3 = NA_character_))
    res <- map_protein_groups_to_entrez(c("P1", "P2;P1", "P3;P1"), lk$fn)
    by_group <- split(res$ENTREZID, res$feature_id)
    expect_identical(by_group[["P1"]], "101")
    # P1 is a secondary member of P2;P1, and gene 101 is P1's own: dropped.
    expect_identical(by_group[["P2;P1"]], "102")
    # In P3;P1 the leading member is unmapped, so P1 is that group's first gene
    # and stays -- two features sharing a first gene is left to the callers,
    # which dedupe by gene as they always have.
    expect_identical(by_group[["P3;P1"]], "101")
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
    expect_identical(names(res), c("feature_id", "ENTREZID", "matched_accession",
                                   "member_rank", "via_canonical"))
    expect_identical(nrow(map_protein_groups_to_entrez(character(0), lk$fn)), 0L)
})

test_that("each key, canonical forms included, is looked up once in one call", {
    lk <- fake_lookup(c(P1 = "101", P2 = "102"))
    map_protein_groups_to_entrez(c("P1;P2", "P2;P1", "P1-2"), lk$fn)
    expect_length(lk$calls$keys, 1)
    expect_setequal(lk$calls$keys[[1]], c("P1", "P2", "P1-2"))
    expect_false(anyDuplicated(lk$calls$keys[[1]]) > 0)
})

test_that("the result does not depend on input order", {
    lk <- fake_lookup(c(P1 = "101", P2 = "102", P3 = "103"))
    a <- map_protein_groups_to_entrez(c("P3;P1", "P1;P2", "P2"), lk$fn)
    b <- map_protein_groups_to_entrez(c("P2", "P1;P2", "P3;P1"), lk$fn)
    expect_identical(a, b)
})
