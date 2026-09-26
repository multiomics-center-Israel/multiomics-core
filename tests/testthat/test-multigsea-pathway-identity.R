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


# ---- the identity object ---------------------------------------------------
#
# One object carries the key and the decision that produced it. Everything that
# needs to know where a row's identity came from asks this rather than guessing
# from the columns, which is what repeatedly put wrong labels on rows.

test_that("the identity object reports source, raw and key for each rung", {
    df <- data.frame(
        term        = c("GO:0006915", NA, NA, NA),
        ID          = c("hsa99999",   "hsa00010", NA, NA),
        pathway     = c("ignored",    "ignored",  "hsa00020 Citrate cycle", NA),
        Description = c("ignored",    "ignored",  "ignored", "hsa00030 PPP"),
        stringsAsFactors = FALSE
    )

    id <- .multigsea_identity(df, kegg_org = "hsa")

    expect_identical(id$source, c("term", "ID", "pathway", "Description"))
    expect_identical(id$raw, c("GO:0006915", "hsa00010",
                               "hsa00020 Citrate cycle", "hsa00030 PPP"))
    expect_identical(id$key, c("GO:0006915", "00010", "00020", "00030"))
})

test_that("blank and NA candidates fall through row by row", {
    df <- data.frame(
        ID      = c(NA, "   ", ""),
        pathway = c("hsa00010", "hsa00020", NA),
        stringsAsFactors = FALSE
    )

    id <- .multigsea_identity(df, kegg_org = "hsa")

    expect_identical(id$source, c("pathway", "pathway", NA))
    expect_identical(id$key, c("00010", "00020", NA))
})

test_that("a table with no candidate column at all yields NULL", {
    expect_null(.multigsea_identity(data.frame(padj = c(0.1, 0.2))))
})

test_that("term_ids is exactly the identity object's key", {
    df <- data.frame(ID = c("hsa00010", NA),
                     pathway = c("x", "map00020"), stringsAsFactors = FALSE)

    expect_identical(.multigsea_term_ids(df, kegg_org = "hsa"),
                     .multigsea_identity(df, kegg_org = "hsa")$key)
})

test_that("without term, the key agrees with pathway_join_key()", {
    # The shared rungs must not drift from #201's ladder. term is MultiGSEA's own
    # addition in front of it, so a frame without term has to match exactly.
    df <- data.frame(
        ID          = c("hsa00010", NA, NA),
        pathway     = c("readable", "map00020", NA),
        Description = c("ignored", "ignored", "ko00030 PPP"),
        stringsAsFactors = FALSE
    )

    expect_identical(.multigsea_identity(df, kegg_org = "hsa")$key,
                     pathway_join_key(df, kegg_org = "hsa"))
})

test_that("term takes precedence over a conflicting ID", {
    df <- data.frame(term = "GO:0006915", ID = "hsa00010",
                     stringsAsFactors = FALSE)

    id <- .multigsea_identity(df, kegg_org = "hsa")

    expect_identical(id$source, "term")
    expect_identical(id$key, "GO:0006915")
})


# ---- labels follow the source the identity actually used -------------------

test_that("a term carrying its own name keeps that name", {
    df <- data.frame(term = "gla00010 Glycolysis / Gluconeogenesis",
                     stringsAsFactors = FALSE)

    nms <- .multigsea_term_names(list(df), kegg_org = "gla")

    expect_identical(unname(nms[["00010"]]), "Glycolysis / Gluconeogenesis")
})

test_that("an ID-only row still displays its accession", {
    # Nothing beside the ID to name it, and normalization has taken the prefix
    # off the key -- so the raw value is the only thing left that says anything.
    df <- data.frame(ID = "hsa00010", padj = 0.01, stringsAsFactors = FALSE)

    nms <- .multigsea_term_names(list(df), kegg_org = "hsa")

    expect_identical(unname(nms[["00010"]]), "hsa00010")
})

test_that("a conflicting term and ID label from term, not from the ID", {
    # The row keyed on term, so the Description beside the ID describes a
    # different pathway and must not be attached to it.
    df <- data.frame(term = "GO:0006915", ID = "hsa00010",
                     Description = "Glycolysis / Gluconeogenesis",
                     stringsAsFactors = FALSE)

    nms <- .multigsea_term_names(list(df), kegg_org = "hsa")

    expect_identical(unname(nms[["GO:0006915"]]), "GO:0006915")
    expect_false("Glycolysis / Gluconeogenesis" %in% unname(nms))
})


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

test_that("a row with no name column still contributes its own raw identifier", {
    # This asserted the opposite until the identity object landed, and the old
    # expectation was the bug: with the key normalized to 00010 there is nothing
    # left to recover downstream, so the row has to carry its accession here.
    rna <- data.frame(pathway = NA_character_, ID = "hsa00010",
                      stringsAsFactors = FALSE)

    nms <- .multigsea_term_names(list(rna), kegg_org = "hsa")

    expect_identical(unname(nms[["00010"]]), "hsa00010")
})

test_that("a row with no usable identifier at all contributes nothing", {
    # The only way a key is absent from the map now: there was no key to begin
    # with. A blank string is never stored under one.
    unidentifiable <- data.frame(ID = NA_character_, pathway = NA_character_,
                                 stringsAsFactors = FALSE)

    expect_length(.multigsea_term_names(list(unidentifiable), kegg_org = "hsa"), 0L)

    # And such a row does not disturb an identifiable one beside it.
    mixed <- data.frame(ID = c(NA_character_, NA_character_),
                        pathway = c(NA_character_, "hsa00010"),
                        stringsAsFactors = FALSE)

    nms <- .multigsea_term_names(list(mixed), kegg_org = "hsa")

    expect_identical(names(nms), "00010")
    expect_true(all(.multigsea_usable_label(nms)))
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


# ---- duplicates one omic can now hold --------------------------------------
#
# Normalizing the identity makes hsa00010 and map00010 the same pathway, so a
# bound per-omic table can hold two rows for it where it held two pathways
# before. The reduction rule is #201's: keep the minimum adjusted p-value.

test_that("two prefix variants in one omic collapse to a single row", {
    df <- data.frame(pathway = c("hsa00010", "map00010", "hsa00020"),
                     padj = c(0.20, 0.001, 0.03), stringsAsFactors = FALSE)
    terms <- .multigsea_term_ids(df, kegg_org = "hsa")

    keep <- .multigsea_collapse_duplicate_terms(df, terms)

    expect_equal(terms[keep], c("00010", "00020"))
    expect_false(anyDuplicated(terms[keep]) > 0)
})

test_that("the retained row is the most significant of the duplicates", {
    df <- data.frame(pathway = c("hsa00010", "map00010"),
                     padj = c(0.20, 0.001), stringsAsFactors = FALSE)
    terms <- .multigsea_term_ids(df, kegg_org = "hsa")

    keep <- .multigsea_collapse_duplicate_terms(df, terms)

    expect_equal(df$padj[keep], 0.001)
})

test_that("a duplicate with no p-value does not displace one that has a value", {
    df <- data.frame(pathway = c("hsa00010", "map00010"),
                     padj = c(NA_real_, 0.04), stringsAsFactors = FALSE)
    terms <- .multigsea_term_ids(df, kegg_org = "hsa")

    keep <- .multigsea_collapse_duplicate_terms(df, terms)

    expect_equal(df$padj[keep], 0.04)
})

test_that("the collapse reduces on the same column the panel scores on", {
    # No padj here, so p.adjust is what MultiGSEA would score on and what the
    # reduction has to read.
    df <- data.frame(pathway = c("hsa00010", "map00010"),
                     p.adjust = c(0.30, 0.002), stringsAsFactors = FALSE)
    terms <- .multigsea_term_ids(df, kegg_org = "hsa")

    keep <- .multigsea_collapse_duplicate_terms(df, terms)

    expect_equal(df$p.adjust[keep], 0.002)
})

test_that("rows with no key survive the collapse rather than folding together", {
    df <- data.frame(pathway = c(NA_character_, NA_character_, "hsa00010"),
                     padj = c(0.01, 0.02, 0.03), stringsAsFactors = FALSE)
    terms <- .multigsea_term_ids(df, kegg_org = "hsa")

    expect_length(.multigsea_collapse_duplicate_terms(df, terms), 3L)
})


test_that("pairwise matching does not depend on row order among duplicates", {
    skip_if_not_installed("ggplot2")

    # Same four rows, two of which are prefix variants of one pathway, in two
    # different orders. match() used to take whichever came first.
    rows <- data.frame(
        pathway = c("hsa00010", "map00010", "hsa00020", "hsa00030"),
        padj    = c(0.20, 0.001, 0.01, 0.02),
        stringsAsFactors = FALSE
    )

    run_pair <- function(idx) {
        out_dir <- withr::local_tempdir()
        suppressMessages(suppressWarnings(
            .save_multigsea_pair_plot(rows[idx, , drop = FALSE],
                                      rows[idx, , drop = FALSE],
                                      "transcriptomics", "metabolomics",
                                      out_dir = out_dir, kegg_org = "hsa")
        ))
        # term is read as character on purpose: a normalized KEGG key is all
        # digits, and read.csv's default type conversion would turn "00010" into
        # the number 10 before the assertions below ever see it.
        written <- utils::read.csv(
            file.path(out_dir, "multigsea_transcriptomics_vs_metabolomics.csv"),
            colClasses = c(term = "character"),
            stringsAsFactors = FALSE
        )
        out <- written[order(written$term), c("term", "x", "y")]
        rownames(out) <- NULL          # so the two orders compare on content alone
        out
    }

    forward <- run_pair(seq_len(4))
    reversed <- run_pair(rev(seq_len(4)))

    expect_equal(nrow(forward), 3L)
    expect_setequal(forward$term, c("00010", "00020", "00030"))
    # Row order changes nothing, and the collapsed pathway kept its best p-value.
    expect_equal(forward, reversed)
    expect_equal(forward$x[forward$term == "00010"], -log10(0.001))
})


test_that("a name on a discarded duplicate still reaches the plot label", {
    skip_if_not_installed("ggplot2")

    # The row that survives the collapse is the most significant one, and here it
    # is the one with no name. Resolving names after the collapse would lose
    # "Glycolysis" permanently and label the pathway with its bare key.
    rows <- data.frame(
        pathway      = c("hsa00010",   "map00010", "hsa00020", "hsa00030"),
        pathway_name = c("Glycolysis", NA,         "TCA",      "PPP"),
        padj         = c(0.20,         0.001,      0.01,       0.02),
        stringsAsFactors = FALSE
    )

    out_dir <- withr::local_tempdir()
    suppressMessages(suppressWarnings(
        .save_multigsea_pair_plot(rows, rows, "transcriptomics", "metabolomics",
                                  out_dir = out_dir, kegg_org = "hsa")
    ))
    written <- utils::read.csv(
        file.path(out_dir, "multigsea_transcriptomics_vs_metabolomics.csv"),
        colClasses = c(term = "character"), stringsAsFactors = FALSE
    )

    # Scoring kept the most significant row...
    expect_equal(written$x[written$term == "00010"], -log10(0.001))
    # ...and display still found the name the discarded duplicate carried.
    expect_identical(written$label[written$term == "00010"], "Glycolysis")
})


test_that("a row with no identity does not take the pairwise plot down with it", {
    skip_if_not_installed("ggplot2")

    # An unkeyable row used to reach resolve_term() as NA, where
    # nchar(NA_character_) is NA and `if (nchar(t) == 0)` errors. The main
    # pairwise loop has no tryCatch around it, so that aborted every MultiGSEA
    # panel, not just this pair.
    rows <- data.frame(
        pathway = c("hsa00010", "hsa00020", "hsa00030", NA),
        padj    = c(0.01, 0.02, 0.03, 0.04),
        stringsAsFactors = FALSE
    )

    out_dir <- withr::local_tempdir()
    expect_no_error(suppressMessages(suppressWarnings(
        .save_multigsea_pair_plot(rows, rows, "transcriptomics", "metabolomics",
                                  out_dir = out_dir, kegg_org = "hsa")
    )))

    written <- utils::read.csv(
        file.path(out_dir, "multigsea_transcriptomics_vs_metabolomics.csv"),
        colClasses = c(term = "character"), stringsAsFactors = FALSE
    )
    expect_setequal(written$term, c("00010", "00020", "00030"))
})


test_that("an omic left with no keyed rows is not plotted against imputed zeros", {
    skip_if_not_installed("ggplot2")

    # Dropping the unkeyable rows can empty one side entirely. union() cannot
    # see that -- the other omic alone carries the pair past the >= 3 check --
    # and match() then scores the empty side as all-NA, i.e. a column of zeros
    # that reads as a real "no enrichment anywhere" result.
    keyed <- data.frame(
        pathway = c("hsa00010", "hsa00020", "hsa00030"),
        padj    = c(0.01, 0.02, 0.03),
        stringsAsFactors = FALSE
    )
    unkeyable <- data.frame(
        pathway = c(NA_character_, NA_character_, NA_character_),
        padj    = c(0.01, 0.02, 0.03),
        stringsAsFactors = FALSE
    )

    out_dir <- withr::local_tempdir()
    expect_no_error(suppressMessages(suppressWarnings(
        .save_multigsea_pair_plot(keyed, unkeyable, "transcriptomics", "metabolomics",
                                  out_dir = out_dir, kegg_org = "hsa")
    )))

    expect_length(list.files(out_dir), 0L)
})


# ---- the synthesized label must not undercut resolve_term() ----------------
#
# Installing a label in the map shadows resolve_term()'s own fallback, so
# anything the synthesized label fails to strip is text that used to be removed.

test_that("the fallback and the synthesized label share one implementation", {
    # resolve_term()'s fallback and the name map used to strip separately, which
    # is how they drifted: the map installed a label the fallback would have
    # shortened, and shadowed it. Both now call this helper, so agreement is by
    # construction -- this pins it so a future edit to either cannot re-split it.
    shapes <- c("GO:0006915~Apoptosis", "gla00010 Glycolysis / Gluconeogenesis")

    fallback <- .multigsea_readable_from_identifier(shapes)
    expect_identical(fallback, c("Apoptosis", "Glycolysis / Gluconeogenesis"))

    synthesized <- vapply(shapes, function(s) {
        nms <- .multigsea_term_names(
            list(data.frame(pathway = s, stringsAsFactors = FALSE)),
            kegg_org = "gla"
        )
        unname(nms[[1]])
    }, character(1), USE.NAMES = FALSE)

    expect_identical(synthesized, fallback)
})

test_that("a GO identifier with a tilde name keeps only the name", {
    df <- data.frame(pathway = "GO:0006915~Apoptosis", stringsAsFactors = FALSE)

    nms <- .multigsea_term_names(list(df), kegg_org = "hsa")

    expect_identical(unname(nms[["GO:0006915~Apoptosis"]]), "Apoptosis")
})

test_that("another species' accession is still stripped for display", {
    # mmu00010 is correctly NOT the same pathway as hsa00010 for joining, so the
    # key keeps its prefix -- but resolve_term() strips any two or three letter
    # prefix for display, and the label has to match that.
    df <- data.frame(pathway = "mmu00010 Glycolysis / Gluconeogenesis",
                     stringsAsFactors = FALSE)

    keys <- .multigsea_term_ids(df, kegg_org = "hsa")
    nms  <- .multigsea_term_names(list(df), kegg_org = "hsa")

    expect_identical(keys, "mmu00010 Glycolysis / Gluconeogenesis")
    expect_identical(unname(nms[[keys]]), "Glycolysis / Gluconeogenesis")
})


# ---- a stated name outranks one synthesized from an identifier -------------

test_that("an explicit pathway_name beats a bare identifier from an earlier omic", {
    # "GO:0006915" is not a KEGG accession, so ranking on that alone called it
    # readable and let it win on input order. It is still not a name.
    rna   <- data.frame(pathway = "GO:0006915", padj = 0.01,
                        stringsAsFactors = FALSE)
    prot  <- data.frame(pathway = "GO:0006915", pathway_name = "Apoptotic process",
                        padj = 0.02, stringsAsFactors = FALSE)

    nms <- .multigsea_term_names(list(rna, prot), kegg_org = "hsa")

    expect_identical(unname(nms[["GO:0006915"]]), "Apoptotic process")
})

test_that("the same holds for a bare PFAM accession", {
    rna  <- data.frame(pathway = "PF00089", padj = 0.01, stringsAsFactors = FALSE)
    prot <- data.frame(pathway = "PF00089", pathway_name = "Trypsin",
                       padj = 0.02, stringsAsFactors = FALSE)

    nms <- .multigsea_term_names(list(rna, prot), kegg_org = "hsa")

    expect_identical(unname(nms[["PF00089"]]), "Trypsin")
})


# ---- co-significance counts omics, not rows --------------------------------

test_that("duplicate rows within one omic do not make a pathway co-significant", {
    # Both rows are the same omic and the same pathway. Counting rows would call
    # this co-significant across two layers; it is one layer twice.
    summary_df <- data.frame(
        term  = c("00010", "00010"),
        omic  = c("transcriptomics", "transcriptomics"),
        padj  = c(0.001, 0.002),
        stringsAsFactors = FALSE
    )

    expect_length(.multigsea_cosignificant_terms(summary_df), 0L)
})

test_that("the same pathway significant in two distinct omics is co-significant", {
    summary_df <- data.frame(
        term  = c("00010", "00010"),
        omic  = c("transcriptomics", "metabolomics"),
        padj  = c(0.001, 0.002),
        stringsAsFactors = FALSE
    )

    expect_equal(.multigsea_cosignificant_terms(summary_df), "00010")
})

test_that("a pathway below threshold in only one omic is not co-significant", {
    summary_df <- data.frame(
        term  = c("00010", "00010"),
        omic  = c("transcriptomics", "metabolomics"),
        padj  = c(0.001, 0.900),
        stringsAsFactors = FALSE
    )

    expect_length(.multigsea_cosignificant_terms(summary_df), 0L)
})


# ---- ID + Description, with no pathway column ------------------------------

test_that("an ID and Description table contributes its readable Description", {
    # A supported clusterProfiler shape. Normalizing the id means the accession
    # cannot be recovered downstream either, so without this the label would be
    # the bare key.
    df <- data.frame(ID = "hsa00010",
                     Description = "Glycolysis / Gluconeogenesis",
                     stringsAsFactors = FALSE)

    keys <- .multigsea_term_ids(df, kegg_org = "hsa")
    nms  <- .multigsea_term_names(list(df), kegg_org = "hsa")

    expect_identical(keys, "00010")
    expect_identical(unname(nms[["00010"]]), "Glycolysis / Gluconeogenesis")
})

test_that("a row keyed on Description reads its readable text from Description", {
    # pathway_join_key() falls all the way through to Description here, so the
    # display side has to follow the same row-wise ladder. Choosing the fallback
    # by "does this row key on ID" alone left this row nameless, and after
    # normalization the accession is gone from the key too -- the panel would
    # have shown a bare 00010.
    df <- data.frame(
        ID          = NA_character_,
        pathway     = NA_character_,
        Description = "gla00010 Glycolysis / Gluconeogenesis",
        stringsAsFactors = FALSE
    )

    keys <- .multigsea_term_ids(df, kegg_org = "gla")
    nms  <- .multigsea_term_names(list(df), kegg_org = "gla")

    expect_identical(keys, "00010")
    expect_identical(unname(nms[["00010"]]), "Glycolysis / Gluconeogenesis")
})

test_that("the identity and display ladders agree row by row in one frame", {
    # One frame, one row per rung of ID -> pathway -> Description.
    df <- data.frame(
        ID          = c("hsa00010", NA, NA),
        pathway     = c("Glycolysis / Gluconeogenesis", "hsa00020 Citrate cycle", NA),
        Description = c(NA, NA, "hsa00030 Pentose phosphate pathway"),
        stringsAsFactors = FALSE
    )

    keys <- .multigsea_term_ids(df, kegg_org = "hsa")
    nms  <- .multigsea_term_names(list(df), kegg_org = "hsa")

    expect_identical(keys, c("00010", "00020", "00030"))
    expect_identical(unname(nms[["00010"]]), "Glycolysis / Gluconeogenesis")
    expect_identical(unname(nms[["00020"]]), "Citrate cycle")
    expect_identical(unname(nms[["00030"]]), "Pentose phosphate pathway")
})

test_that("Description fills in per row where pathway is blank", {
    df <- data.frame(
        ID          = c("hsa00010", "hsa00020"),
        pathway     = c("Glycolysis / Gluconeogenesis", NA),
        Description = c("ignored", "Citrate cycle (TCA cycle)"),
        stringsAsFactors = FALSE
    )

    nms <- .multigsea_term_names(list(df), kegg_org = "hsa")

    # `pathway` keeps its precedence where it says something.
    expect_identical(unname(nms[["00010"]]), "Glycolysis / Gluconeogenesis")
    expect_identical(unname(nms[["00020"]]), "Citrate cycle (TCA cycle)")
})


# ---- one label per normalized term across omics ----------------------------

test_that("bare accessions from two omics resolve to a single label", {
    # Two layers, no names anywhere. Resolved per omic they would label the same
    # pathway "hsa00010" and "map00010" and take a y-axis row each.
    rna   <- data.frame(pathway = "hsa00010", padj = 0.01, stringsAsFactors = FALSE)
    metab <- data.frame(pathway = "map00010", padj = 0.02, stringsAsFactors = FALSE)

    nms <- .multigsea_term_names(list(rna, metab), kegg_org = "hsa")

    expect_length(nms, 1L)
    expect_identical(names(nms), "00010")
})

test_that("a readable name from either omic beats an accession-only label", {
    # The accession-only layer is first, and must not fix the display name.
    rna   <- data.frame(pathway = "hsa00010", padj = 0.01, stringsAsFactors = FALSE)
    metab <- mg_frame("map00010", name = "Glycolysis / Gluconeogenesis")

    nms <- .multigsea_term_names(list(rna, metab), kegg_org = "hsa")

    expect_identical(unname(nms[["00010"]]), "Glycolysis / Gluconeogenesis")
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

# ---- outputs of an earlier run -------------------------------------------

# A MultiGSEA directory as an earlier run left it, beside Multi-ORA's output and
# a file nothing here wrote.
stale_multigsea_dir <- function(envir = parent.frame()) {
    out <- withr::local_tempdir(.local_envir = envir)
    old <- c(file.path(out, c("multigsea_transcriptomics_vs_proteomics.png",
                              "multigsea_transcriptomics_vs_proteomics.pdf",
                              "multigsea_transcriptomics_vs_proteomics.csv",
                              "multigsea_combined_enrichment.png")),
             file.path(out, "per_contrast", "A_vs_B",
                       "multigsea_transcriptomics_vs_proteomics.png"))
    kept <- c(file.path(out, "multi_ora", "multi_ora_pooled_barplot.png"),
              file.path(out, "multi_ora", "per_contrast", "A_vs_B", "multi_ora_summary.csv"),
              file.path(out, "notes.txt"))
    for (f in c(old, kept)) {
        dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
        writeLines("x", f)
    }
    list(out = out, old = old, kept = kept)
}

test_that("clear_multigsea_outputs removes only what MultiGSEA wrote", {
    d <- stale_multigsea_dir()
    expect_message(clear_multigsea_outputs(d$out), "cleared")
    expect_false(any(file.exists(d$old)))
    expect_false(dir.exists(file.path(d$out, "per_contrast")))
    # Multi-ORA shares the directory and keeps its own per_contrast/.
    expect_true(all(file.exists(d$kept)))
    expect_true(dir.exists(file.path(d$out, "multi_ora", "per_contrast")))
})

test_that("clear_multigsea_outputs is quiet on a fresh or absent directory", {
    expect_length(clear_multigsea_outputs(withr::local_tempdir()), 0L)
    expect_length(clear_multigsea_outputs(NULL), 0L)
    expect_length(clear_multigsea_outputs(file.path(withr::local_tempdir(), "none")), 0L)
})

test_that("a MultiGSEA run that draws nothing still clears the previous run's figures", {
    # The report globs these files, so an early return must not leave the last
    # run's figures to be shown as this run's.
    no_results <- stale_multigsea_dir()
    expect_null(suppressMessages(run_multigsea_plots(NULL, list(), no_results$out)))
    expect_false(any(file.exists(no_results$old)))
    expect_true(all(file.exists(no_results$kept)))

    switched_off <- stale_multigsea_dir()
    cfg <- list(modes = list(multiomics = list(enrichment = list(
        multigsea = list(run_multigsea = FALSE)))))
    expect_null(suppressMessages(run_multigsea_plots(
        list(per_omics = list(transcriptomics = data.frame())), cfg, switched_off$out)))
    expect_false(any(file.exists(switched_off$old)))
    expect_true(all(file.exists(switched_off$kept)))

    one_layer <- stale_multigsea_dir()
    expect_null(suppressMessages(run_multigsea_plots(
        list(per_omics = list(transcriptomics = data.frame())),
        list(global = list(organism = "Homo sapiens")), one_layer$out)))
    expect_false(any(file.exists(one_layer$old)))
    expect_true(all(file.exists(one_layer$kept)))
})

test_that("the pipeline clears MultiGSEA outputs when it skips or fails the step", {
    # Without cross-omics enrichment the target returns before calling
    # run_multigsea_plots(), and an error part-way through is caught here after
    # some pairs are written, so the cleanup has to run on both paths too.
    src <- paste(readLines(testthat::test_path(
        "..", "..", "R", "pipeline", "multiomics", "00_pipe_multiomics.R")),
        collapse = "\n")
    block <- regmatches(src, regexpr(
        "(?s)multiomics_multigsea,.*?\n        \\),", src, perl = TRUE))
    expect_length(block, 1)
    skip <- regmatches(block, regexpr(
        "(?s)Skipping MultiGSEA plots.*?return\\(NULL\\)", block, perl = TRUE))
    expect_match(skip, "clear_multigsea_outputs(mg_dir)", fixed = TRUE)
    # And when the step fails part-way, so no partial set is left behind.
    err <- regmatches(block, regexpr(
        "(?s)MultiGSEA plots failed.*?NULL\n", block, perl = TRUE))
    expect_match(err, "clear_multigsea_outputs(mg_dir)", fixed = TRUE)
    # And it clears the directory the step writes to.
    expect_match(block, "out_dir = mg_dir", fixed = TRUE)
})
