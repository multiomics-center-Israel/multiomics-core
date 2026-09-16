# tests/testthat/test-multigsea-pathway-identity.R
#
# MultiGSEA reads the same per-omics frames the cross-omics join does, so it
# faces the same identity problem #201 fixed there: a gene layer keyed on
# hsa00010 and a compound layer keyed on map00010 were two unrelated pathways.
# .multigsea_term_ids() now delegates to pathway_join_key(), so both sides key on
# the normalized accession while GO, PFAM, InterPro and custom names are left
# exactly as they are.
#
# Identity and display stay apart. resolve_term() shows the readable name where
# the table has one, and recovers it from an "<accession> <name>" identifier
# where it does not -- normalizing the key must not turn that label into a bare
# map number.
#
# All data below is synthetic.

mg_frame <- function(pathway, name = NULL, id = NULL, padj = NULL) {
    df <- data.frame(pathway = pathway, stringsAsFactors = FALSE)
    if (!is.null(name)) df$pathway_name <- name
    if (!is.null(id)) df$ID <- id
    df$padj <- if (is.null(padj)) rep(0.01, length(pathway)) else padj
    df
}


# ---- identity -------------------------------------------------------------

test_that("the KEGG forms collapse to one term id under the run's organism", {
    ids <- c("hsa00010", "map00010", "ko00010", "00010")
    expect_equal(.multigsea_term_ids(mg_frame(ids), kegg_org = "hsa"),
                 rep("00010", 4))
})

test_that("an accession carrying its own name yields the same id", {
    # The non-model fallback names sets "<accession> <readable name>".
    expect_equal(
        .multigsea_term_ids(
            mg_frame(c("gla00010 Glycolysis / Gluconeogenesis", "map00010")),
            kegg_org = "gla"
        ),
        c("00010", "00010")
    )
})

test_that("GO, PFAM, InterPro and custom names are byte-identical", {
    ids <- c("GO:0006915", "PF00089", "IPR001", "HALLMARK_APOPTOSIS")
    expect_identical(.multigsea_term_ids(mg_frame(ids), kegg_org = "hsa"), ids)
})

test_that("another organism's prefix is left alone", {
    expect_identical(.multigsea_term_ids(mg_frame("mmu00010"), kegg_org = "hsa"),
                     "mmu00010")
})

test_that("term keeps its precedence over the other columns", {
    df <- mg_frame("map00020")
    df$term <- "GO:0006915"
    expect_identical(.multigsea_term_ids(df, kegg_org = "hsa"), "GO:0006915")
})

test_that("a frame with no identifier column still returns NULL", {
    expect_null(.multigsea_term_ids(data.frame(padj = c(0.1, 0.2))))
})

test_that("two omics keyed on different KEGG forms now share a term", {
    rna  <- mg_frame("hsa00010", name = "Glycolysis / Gluconeogenesis")
    metab <- mg_frame("map00010", name = "Glycolysis / Gluconeogenesis")

    t1 <- .multigsea_term_ids(rna, kegg_org = "hsa")
    t2 <- .multigsea_term_ids(metab, kegg_org = "hsa")

    expect_equal(union(t1, t2), "00010")
    expect_false(is.na(match(t1, t2)))
})


# ---- the name map lines up with the ids it is looked up by ----------------

test_that("term names are keyed by the same ids term_ids returns", {
    # The coupling that matters: plot_multigsea_combined() does
    # .multigsea_term_names(list(res))[term_ids], so a key built from a different
    # column would silently look up NA for every row.
    df <- mg_frame("hsa00010", name = "Glycolysis / Gluconeogenesis")

    keys <- .multigsea_term_ids(df, kegg_org = "hsa")
    nms  <- .multigsea_term_names(list(df), kegg_org = "hsa")

    expect_true(all(keys %in% names(nms)))
    expect_identical(unname(nms[keys]), "Glycolysis / Gluconeogenesis")
})

test_that("a frame with pathway, ID and pathway_name resolves by the id in use", {
    # term_ids prefers ID here, so the name map has to be keyed on ID as well --
    # keying it on pathway is the mismatch this change removes.
    df <- mg_frame("Glycolysis / Gluconeogenesis", name = "Glycolysis",
                   id = "hsa00010")

    keys <- .multigsea_term_ids(df, kegg_org = "hsa")
    nms  <- .multigsea_term_names(list(df), kegg_org = "hsa")

    expect_identical(keys, "00010")
    expect_identical(unname(nms[keys]), "Glycolysis")
})


# ---- display is preserved, shape by shape ---------------------------------

test_that("an <accession> <name> identifier still displays without the accession", {
    # Before this change resolve_term() found no mapping and stripped the leading
    # accession off the term itself. Normalizing the key removes the prefix it
    # stripped, so the name map now has to carry the readable remainder -- and it
    # must be the remainder, not the whole string with the accession back in it.
    df <- mg_frame("gla00010 Glycolysis / Gluconeogenesis")

    keys <- .multigsea_term_ids(df, kegg_org = "gla")
    nms  <- .multigsea_term_names(list(df), kegg_org = "gla")

    expect_identical(keys, "00010")
    expect_identical(unname(nms[keys]), "Glycolysis / Gluconeogenesis")
})

test_that("a bare accession with no name displays as it does today", {
    # Stripping leaves nothing, so the identifier itself is the label -- the same
    # rule resolve_term() applies with its `if (nchar(t) == 0) t <- term`.
    df <- mg_frame("hsa00010")

    nms <- .multigsea_term_names(list(df), kegg_org = "hsa")

    expect_identical(unname(nms[["00010"]]), "hsa00010")
})

test_that("a non-KEGG identifier with no name is unchanged as a label", {
    df <- mg_frame(c("GO:0006915", "HALLMARK_APOPTOSIS"))

    nms <- .multigsea_term_names(list(df), kegg_org = "hsa")

    expect_identical(unname(nms[["GO:0006915"]]), "GO:0006915")
    expect_identical(unname(nms[["HALLMARK_APOPTOSIS"]]), "HALLMARK_APOPTOSIS")
})

test_that("a real pathway_name still wins over anything in the identifier", {
    df <- mg_frame("gla00010 Glycolysis / Gluconeogenesis",
                   name = "Glycolysis")

    nms <- .multigsea_term_names(list(df), kegg_org = "gla")

    expect_identical(unname(nms[["00010"]]), "Glycolysis")
})

test_that("the first name seen wins across frames, as before", {
    a <- mg_frame("hsa00010", name = "Glycolysis / Gluconeogenesis")
    b <- mg_frame("map00010", name = "Glycolysis")

    nms <- .multigsea_term_names(list(a, b), kegg_org = "hsa")

    expect_length(nms, 1L)
    expect_identical(unname(nms[["00010"]]), "Glycolysis / Gluconeogenesis")
})


# ---- a blank label must not reserve a normalized key ----------------------
#
# This collision only exists because of normalization: hsa00010 and map00010 did
# not share a key before, so a layer with no usable name could not block a layer
# that had one.

test_that("an NA label does not lock out a readable one from another layer", {
    rna   <- data.frame(pathway = NA_character_, ID = "hsa00010",
                        stringsAsFactors = FALSE)
    metab <- mg_frame("map00010", name = "Glycolysis / Gluconeogenesis")

    nms <- .multigsea_term_names(list(rna, metab), kegg_org = "hsa")

    expect_identical(unname(nms[["00010"]]), "Glycolysis / Gluconeogenesis")
})

test_that("an empty label does not lock out a readable one from another layer", {
    rna   <- data.frame(pathway = "", ID = "hsa00010", stringsAsFactors = FALSE)
    metab <- mg_frame("map00010", name = "Glycolysis / Gluconeogenesis")

    nms <- .multigsea_term_names(list(rna, metab), kegg_org = "hsa")

    expect_identical(unname(nms[["00010"]]), "Glycolysis / Gluconeogenesis")
})

test_that("a whitespace-only label does not lock out a readable one", {
    rna   <- data.frame(pathway = "   \t ", ID = "hsa00010",
                        stringsAsFactors = FALSE)
    metab <- mg_frame("map00010", name = "Glycolysis / Gluconeogenesis")

    nms <- .multigsea_term_names(list(rna, metab), kegg_org = "hsa")

    expect_identical(unname(nms[["00010"]]), "Glycolysis / Gluconeogenesis")
})

test_that("a key no layer can name is simply absent from the map", {
    # resolve_term() and the combined panel both fall back to the id itself, so
    # an absent key is the right outcome -- not a blank string stored under it.
    rna <- data.frame(pathway = NA_character_, ID = "hsa00010",
                      stringsAsFactors = FALSE)

    nms <- .multigsea_term_names(list(rna), kegg_org = "hsa")

    expect_false("00010" %in% names(nms))
})

test_that("pathway_name missing on only some rows does not blank those rows", {
    # The column exists, so a whole-column branch would hand back NA for row 1
    # even though its identifier carries the readable text.
    df <- mg_frame(c("gla00010 Glycolysis / Gluconeogenesis",
                     "gla00020 Citrate cycle"),
                   name = c(NA, "Citrate cycle (TCA cycle)"))

    nms <- .multigsea_term_names(list(df), kegg_org = "gla")

    expect_identical(unname(nms[["00010"]]), "Glycolysis / Gluconeogenesis")
    expect_identical(unname(nms[["00020"]]), "Citrate cycle (TCA cycle)")
})

test_that("a frame with ID on some rows only falls back per row", {
    # pathway_join_key() decides identity row by row, so row 2 keys on `pathway`
    # and normalizes to 00020. Its readable fallback has to behave like an
    # identifier-only row too -- a whole-column ID branch would leave the
    # accession sitting in the label as "hsa00020 Citrate cycle".
    df <- data.frame(
        ID      = c("hsa00010", NA),
        pathway = c("Glycolysis / Gluconeogenesis", "hsa00020 Citrate cycle"),
        stringsAsFactors = FALSE
    )

    nms <- .multigsea_term_names(list(df), kegg_org = "hsa")

    expect_identical(unname(nms[["00010"]]), "Glycolysis / Gluconeogenesis")
    expect_identical(unname(nms[["00020"]]), "Citrate cycle")
})

test_that("a blank ID counts as absent for the readable fallback too", {
    # Same usable-value rule as pathway_join_key() uses to skip a blank ID.
    df <- data.frame(
        ID      = c("hsa00010", "   "),
        pathway = c("Glycolysis / Gluconeogenesis", "hsa00020 Citrate cycle"),
        stringsAsFactors = FALSE
    )

    nms <- .multigsea_term_names(list(df), kegg_org = "hsa")

    expect_identical(unname(nms[["00020"]]), "Citrate cycle")
})

test_that("an empty pathway_name on one row falls back to that row's identifier", {
    df <- mg_frame(c("hsa00010", "hsa00020"),
                   name = c("   ", "Citrate cycle (TCA cycle)"))

    nms <- .multigsea_term_names(list(df), kegg_org = "hsa")

    expect_identical(unname(nms[["00010"]]), "hsa00010")
    expect_identical(unname(nms[["00020"]]), "Citrate cycle (TCA cycle)")
})


# ---- no organism configured ----------------------------------------------

test_that("with no KEGG organism the behaviour is what it is today", {
    # kegg_org NULL leaves the species-prefixed forms alone, which is how every
    # existing caller that does not pass one behaves.
    ids <- c("hsa00010", "GO:0006915", "HALLMARK_APOPTOSIS")
    expect_identical(.multigsea_term_ids(mg_frame(ids)), ids)
    expect_identical(.multigsea_term_ids(mg_frame(ids), kegg_org = NULL), ids)
})

test_that("map and ko still normalize without an organism", {
    expect_equal(.multigsea_term_ids(mg_frame(c("map00010", "ko00010"))),
                 rep("00010", 2))
})
