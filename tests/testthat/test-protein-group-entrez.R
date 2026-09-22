# tests/testthat/test-protein-group-entrez.R
#
# Protein groups ("P1;P2") used to be passed whole to the UniProt -> Entrez
# lookup, where no group string is ever a key, so every multi-protein group fell
# out of the KEGG enrichment, the loadings enrichment and the pathway maps. The
# lookup now goes member by member and a group takes its first member that maps.
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

test_that("a group takes the leading member when it maps", {
    lk <- fake_lookup(c(P1 = "101", P2 = "102"))
    res <- map_protein_groups_to_entrez("P1;P2", lk$fn)
    expect_identical(res$ENTREZID, "101")
    expect_identical(res$matched_accession, "P1")
})

test_that("a group falls back to the next member when the leading one does not map", {
    lk <- fake_lookup(c(P1 = NA_character_, P2 = "102", P3 = "103"))
    res <- map_protein_groups_to_entrez(c("P1;P2;P3", "P3"), lk$fn)
    res <- res[order(res$feature_id), ]
    expect_identical(res$feature_id, c("P1;P2;P3", "P3"))
    expect_identical(res$ENTREZID, c("102", "103"))
    expect_identical(res$matched_accession, c("P2", "P3"))
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
    expect_identical(names(res), c("feature_id", "ENTREZID", "matched_accession"))
    expect_identical(nrow(map_protein_groups_to_entrez(character(0), lk$fn)), 0L)
})

test_that("each accession is looked up once, in one call", {
    lk <- fake_lookup(c(P1 = "101", P2 = "102"))
    map_protein_groups_to_entrez(c("P1;P2", "P2;P1", "P1"), lk$fn)
    expect_length(lk$calls$keys, 1)
    expect_setequal(lk$calls$keys[[1]], c("P1", "P2"))
    expect_false(anyDuplicated(lk$calls$keys[[1]]) > 0)
})

test_that("member order decides, not input order", {
    lk <- fake_lookup(c(P1 = "101", P2 = "102"))
    a <- map_protein_groups_to_entrez(c("P2;P1", "P1;P2"), lk$fn)
    b <- map_protein_groups_to_entrez(c("P1;P2", "P2;P1"), lk$fn)
    expect_identical(a, b)
    expect_identical(a$ENTREZID[a$feature_id == "P2;P1"], "102")
})
