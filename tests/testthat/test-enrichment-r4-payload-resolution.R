# tests/testthat/test-enrichment-r4-payload-resolution.R
#
# §17 R4 — Compact Payload Integration. The R2 manifest hierarchy
# (subgroup/group_type/group_label) + the R3 result reuse resolve, at Shiny read
# time, to the exact ORA table via resolve_ora_selection() using ONLY the compact
# whole-table nodes in enrichment$ora — no duplicated per-subgroup storage. The
# validator enforces this end-to-end. payload_version is unaffected (stays "2.0").

# Build a small but realistic compact `enrichment` block whose manifest is
# produced by the REAL R2 code (.build_ora_manifest_rows), so manifest counts and
# the whole-table ora nodes are consistent BY CONSTRUCTION.
.mk_enrichment_r4 <- function() {
    cw  <- data.frame(Cluster = c("up", "up", "down"), ID = c("a", "b", "c"),
                      geneID = c("g1/g2", "g3", "g4"),
                      p.adjust = c(1e-3, 2e-3, 3e-3), stringsAsFactors = FALSE)
    cwo <- data.frame(Cluster = c("all", "all"), ID = c("a", "d"),
                      geneID = c("g1/g2", "g5"),
                      p.adjust = c(1e-3, 2e-3), stringsAsFactors = FALSE)
    pt  <- data.frame(Cluster = c("1", "1", "2"), ID = c("a", "b", "c"),
                      geneID = c("g1", "g2", "g3"),
                      p.adjust = c(1e-3, 2e-3, 3e-3), stringsAsFactors = FALSE)
    bp  <- data.frame(Cluster = c("0011", "0101"), ID = c("a", "b"),
                      geneID = c("g1", "g2"),
                      p.adjust = c(1e-3, 2e-3), stringsAsFactors = FALSE)
    adf <- data.frame(Cluster = "all_genes", ID = "a", geneID = "g1",
                      p.adjust = 1e-3, stringsAsFactors = FALSE)

    ora <- list(KEGG = list(
        contrasts_with_direction    = list(c1 = cw),
        contrasts_without_direction = list(c1 = cwo),
        partition                   = pt,
        binary_patterns             = bp,
        all_DE                      = adf))

    manifest <- rbind(
        .build_ora_manifest_rows("contrasts", "c1", "KEGG", cw, FALSE, NULL,
                                 "KEGG_contrasts_c1_ora"),
        .build_ora_manifest_rows("contrasts_wo_direction", "c1", "KEGG", cwo, FALSE, NULL,
                                 "KEGG_contrasts_wo_direction_c1_ora"),
        .build_ora_manifest_rows("partition", "k", "KEGG", pt, FALSE,
                                 setNames(c("1", "1", "2"), c("g1", "g2", "g3")),
                                 "KEGG_partition_k_ora"),
        .build_ora_manifest_rows("binary_patterns", "best", "KEGG", bp, FALSE,
                                 setNames(c("0011", "0101"), c("g1", "g2")),
                                 "KEGG_binary_patterns_best_ora"),
        .build_ora_manifest_rows("all_DE", "any_contrast", "KEGG", adf, FALSE, NULL,
                                 "KEGG_all_DE_any_contrast_ora"))

    list(available = TRUE, gene_index = paste0("g", 1:8), config = list(),
         manifest = manifest, ora = ora, gsea = list(),
         gsea_leading_edge = list(), gsea_rankings = list(),
         pathway_membership = list())
}

.sel_n <- function(enr, ...) {
    s <- resolve_ora_selection(enr, ...)
    if (is.null(s)) NA_integer_ else nrow(s)
}

# ---------------------------------------------------------------------------
# resolve_ora_selection — UP / DOWN / UP+DOWN / All
# ---------------------------------------------------------------------------

test_that("UP / DOWN filter the ONE with_direction table; UP+DOWN is that whole table", {
    enr <- .mk_enrichment_r4()
    up <- resolve_ora_selection(enr, "KEGG", "contrasts", "c1", "up",      "direction")
    dn <- resolve_ora_selection(enr, "KEGG", "contrasts", "c1", "down",    "direction")
    ud <- resolve_ora_selection(enr, "KEGG", "contrasts", "c1", "up_down", "direction")
    expect_true(all(up$Cluster == "up"));   expect_equal(nrow(up), 2L)
    expect_true(all(dn$Cluster == "down")); expect_equal(nrow(dn), 1L)
    expect_setequal(unique(ud$Cluster), c("up", "down"))
    expect_equal(nrow(ud), 3L)                       # whole with_direction table
    expect_equal(nrow(up) + nrow(dn), nrow(ud))      # reuse, not a recomputed union
})

test_that("All resolves to the DISTINCT without_direction node (UP+DOWN != All)", {
    enr <- .mk_enrichment_r4()
    all_sel <- resolve_ora_selection(enr, "KEGG", "contrasts", "c1", "all", "direction")
    ud      <- resolve_ora_selection(enr, "KEGG", "contrasts", "c1", "up_down", "direction")
    expect_true(all(all_sel$Cluster == "all"))
    expect_equal(nrow(all_sel), 2L)                  # whole without_direction table
    # different node, different content => never collapsed / never a union
    expect_false(isTRUE(all.equal(all_sel, ud)))
    expect_false(identical(nrow(all_sel), nrow(ud)))
})

# ---------------------------------------------------------------------------
# Compare Clusters / Compare Patterns
# ---------------------------------------------------------------------------

test_that("a single cluster/pattern filters; Compare returns the whole multi-group table", {
    enr <- .mk_enrichment_r4()
    c1 <- resolve_ora_selection(enr, "KEGG", "partition", "1", NA_character_, "partition")
    cc <- resolve_ora_selection(enr, "KEGG", "partition", "compare", NA_character_, "partition")
    expect_true(all(c1$Cluster == "1")); expect_equal(nrow(c1), 2L)
    expect_setequal(unique(cc$Cluster), c("1", "2")); expect_equal(nrow(cc), 3L)

    p0011 <- resolve_ora_selection(enr, "KEGG", "binary_patterns", "0011", NA_character_, "pattern")
    cp    <- resolve_ora_selection(enr, "KEGG", "binary_patterns", "compare", NA_character_, "pattern")
    expect_true(all(p0011$Cluster == "0011")); expect_equal(nrow(p0011), 1L)
    expect_setequal(unique(cp$Cluster), c("0011", "0101")); expect_equal(nrow(cp), 2L)
})

test_that("all_DE resolves to its single whole table", {
    enr <- .mk_enrichment_r4()
    a <- resolve_ora_selection(enr, "KEGG", "all_DE", "any_contrast", NA_character_, NA_character_)
    expect_equal(nrow(a), 1L)
})

# ---------------------------------------------------------------------------
# compact-storage invariant: no duplicated per-subgroup tables
# ---------------------------------------------------------------------------

test_that("enrichment$ora holds ONLY whole-table nodes (no up/down/compare copies)", {
    enr <- .mk_enrichment_r4()
    expect_setequal(names(enr$ora$KEGG),
                    c("contrasts_with_direction", "contrasts_without_direction",
                      "partition", "binary_patterns", "all_DE"))
    # contrasts node keyed only by contrast (one whole table), never by subgroup
    expect_equal(names(enr$ora$KEGG$contrasts_with_direction), "c1")
    expect_false(any(c("up", "down", "up_down", "compare") %in%
                     names(enr$ora$KEGG$contrasts_with_direction)))
    # partition/binary are a single table each, not per-cluster copies
    expect_s3_class(enr$ora$KEGG$partition, "data.frame")
    expect_s3_class(enr$ora$KEGG$binary_patterns, "data.frame")
})

# ---------------------------------------------------------------------------
# validator — accepts a consistent block, rejects inconsistencies
# ---------------------------------------------------------------------------

test_that("validate_enrichment_payload accepts an R2/R4-consistent block", {
    enr <- .mk_enrichment_r4()
    expect_invisible(validate_enrichment_payload(enr, strict = TRUE))
    expect_true(validate_enrichment_payload(enr, strict = TRUE))
})

test_that("validator rejects a manifest count that does not match the resolved table", {
    enr <- .mk_enrichment_r4()
    # inflate the UP subgroup count (2 -> 99): resolution now disagrees
    up_row <- enr$manifest$analysis == "ORA" & enr$manifest$group == "contrasts" &
              enr$manifest$subgroup == "up"
    enr$manifest$n_significant[up_row] <- 99L
    expect_error(validate_enrichment_payload(enr, strict = TRUE), "n_significant")
})

test_that("validator rejects an out-of-vocabulary subgroup / group_type", {
    enr <- .mk_enrichment_r4()
    bad_sg <- enr; bad_sg$manifest$subgroup[1] <- "sideways"
    expect_error(validate_enrichment_payload(bad_sg, strict = TRUE), "subgroup")

    bad_gt <- .mk_enrichment_r4(); bad_gt$manifest$group_type[1] <- "phylum"
    expect_error(validate_enrichment_payload(bad_gt, strict = TRUE), "group_type")
})

test_that("validator tolerates an evaluated-but-empty subgroup (0 filtered rows)", {
    # add an empty 'up' scenario: a contrast whose with_direction table has only down
    cw_downonly <- data.frame(Cluster = c("down", "down"), ID = c("c", "d"),
                              geneID = c("g4", "g5"), p.adjust = c(1e-3, 2e-3),
                              stringsAsFactors = FALSE)
    man <- .build_ora_manifest_rows("contrasts", "c9", "KEGG", cw_downonly, FALSE, NULL,
                                    "KEGG_contrasts_c9_ora")
    enr <- list(available = TRUE, gene_index = paste0("g", 1:6), config = list(),
                manifest = man,
                ora = list(KEGG = list(contrasts_with_direction = list(c9 = cw_downonly),
                                       contrasts_without_direction = list())),
                gsea = list(), gsea_leading_edge = list(),
                gsea_rankings = list(), pathway_membership = list())
    # up row is status 'empty' (0 rows) -> validator must accept
    expect_true(validate_enrichment_payload(enr, strict = TRUE))
    expect_equal(.sel_n(enr, "KEGG", "contrasts", "c9", "up", "direction"), 0L)
    expect_equal(.sel_n(enr, "KEGG", "contrasts", "c9", "down", "direction"), 2L)
})
