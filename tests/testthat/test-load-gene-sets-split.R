# tests/testthat/test-load-gene-sets-split.R
#
# Unit tests for load_gene_sets() keeping several custom GMTs apart
# (R/core/09_enrichment.R). One GMT stays the historical "custom" collection;
# several GMTs become one collection each, named after the file, so every
# source is scored and FDR-corrected on its own instead of being pooled.
#
# Also covers add_pathway_names() preferring the names a collection already
# carries, which is what keeps a custom GO GMT naming its own terms once the
# collection is called "GO_<something>" rather than "custom".

write_gmt_file <- function(lines, name) {
    path <- file.path(tempfile("gmtdir"), name)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(lines, path)
    path
}

test_that("a single GMT stays one collection named 'custom'", {
    gmt <- write_gmt_file(
        paste(c("PF00089", "Trypsin", paste0("g", 1:4)), collapse = "\t"),
        "PFAM_spalangia.gmt")
    on.exit(unlink(dirname(gmt), recursive = TRUE), add = TRUE)

    gs <- load_gene_sets(organism = "Spalangia cameroni", gmt_file = gmt)

    expect_equal(names(gs), "custom")
    expect_equal(gs$custom[["PF00089"]], paste0("g", 1:4))
})

test_that("several GMTs become one collection each, named after the file", {
    kegg <- write_gmt_file(
        paste(c("map00500", "Starch and sucrose metabolism", paste0("g", 1:5)),
              collapse = "\t"),
        "KEGG_spalangia.gmt")
    pfam <- write_gmt_file(
        paste(c("PF00089", "Trypsin", paste0("g", 4:9)), collapse = "\t"),
        "PFAM_spalangia.gmt")
    on.exit(unlink(c(dirname(kegg), dirname(pfam)), recursive = TRUE), add = TRUE)

    gs <- load_gene_sets(organism = "Spalangia cameroni",
                         gmt_file = list(kegg, pfam))

    expect_setequal(names(gs), c("KEGG_spalangia", "PFAM_spalangia"))
    expect_equal(gs$KEGG_spalangia[["map00500"]], paste0("g", 1:5))
    expect_equal(gs$PFAM_spalangia[["PF00089"]], paste0("g", 4:9))
})

test_that("each split collection keeps its own descriptions", {
    kegg <- write_gmt_file(
        paste(c("map00500", "Starch and sucrose metabolism", paste0("g", 1:5)),
              collapse = "\t"),
        "KEGG_spalangia.gmt")
    ipr <- write_gmt_file(
        paste(c("IPR001254", "Serine proteases, trypsin domain", paste0("g", 2:6)),
              collapse = "\t"),
        "InterPro_spalangia.gmt")
    on.exit(unlink(c(dirname(kegg), dirname(ipr)), recursive = TRUE), add = TRUE)

    gs <- load_gene_sets(organism = "Spalangia cameroni", gmt_file = c(kegg, ipr))

    expect_equal(attr(gs$KEGG_spalangia, "descriptions")[["map00500"]],
                 "Starch and sucrose metabolism")
    expect_equal(attr(gs$InterPro_spalangia, "descriptions")[["IPR001254"]],
                 "Serine proteases, trypsin domain")
})

test_that("same basename in different directories still yields distinct names", {
    a <- write_gmt_file(paste(c("A", "set a", "g1", "g2"), collapse = "\t"),
                        "pathways.gmt")
    b <- write_gmt_file(paste(c("B", "set b", "g3", "g4"), collapse = "\t"),
                        "pathways.gmt")
    on.exit(unlink(c(dirname(a), dirname(b)), recursive = TRUE), add = TRUE)

    gs <- load_gene_sets(organism = "Spalangia cameroni", gmt_file = c(a, b))

    expect_length(gs, 2)
    expect_false(anyDuplicated(names(gs)) > 0)
    expect_true("pathways" %in% names(gs))
})

test_that("a missing path is skipped without dropping the readable ones", {
    ok <- write_gmt_file(paste(c("A", "set a", "g1", "g2"), collapse = "\t"),
                         "KEGG_spalangia.gmt")
    on.exit(unlink(dirname(ok), recursive = TRUE), add = TRUE)

    expect_warning(
        gs <- load_gene_sets(organism = "Spalangia cameroni",
                             gmt_file = c(ok, "/nonexistent/PFAM.gmt")),
        "not found")

    # One surviving path is a single GMT again, so it keeps the "custom" name.
    expect_equal(names(gs), "custom")
    expect_equal(gs$custom$A, c("g1", "g2"))
})

test_that("add_pathway_names prefers the collection's own GO descriptions", {
    gs <- list("GO:0000002" = c("g1", "g2"))
    attr(gs, "descriptions") <- c("GO:0000002" = "mitochondrial genome maintenance [BP]")
    df <- data.frame(pathway = "GO:0000002", padj = 0.01)

    out <- add_pathway_names(df, "GO_spalangia", gs)

    expect_equal(out$pathway_name, "mitochondrial genome maintenance [BP]")
})
