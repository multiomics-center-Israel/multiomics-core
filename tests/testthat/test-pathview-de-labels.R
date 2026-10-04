# Naming the changed genes on a KEGG map, and not drawing a map twice.
#
# A KEGG box can stand for several genes but prints one name -- an EC number,
# or the first gene of the node -- so a box coloured because one member changed
# can carry another member's name. de_gene_node_labels() decides which boxes
# get the changed genes' symbols written above them. All fixtures synthetic;
# no pathview call, no network.

nodes_fixture <- function() {
    data.frame(
        kegg.names = c("11", "21", "31", "41"),
        labels = c("GeneA1", "GeneB", "1.2.3.4", "GeneD"),
        all.mapped = c("11,12,13", "21", "31,32", ""),
        x = c(100, 200, 300, 400), y = c(50, 60, 70, 80),
        width = 46, height = 17,
        stringsAsFactors = FALSE
    )
}
symbols_fixture <- c("11" = "GeneA1", "12" = "GeneA2", "13" = "GeneA3",
                     "21" = "GeneB", "31" = "GeneC1", "32" = "GeneC2")


test_that("a box is labelled with the member that changed, not the one printed", {
    lab <- de_gene_node_labels(nodes_fixture(), de_entrez = c("12", "31"),
                               symbols = symbols_fixture)

    # GeneA2 changed but the box prints GeneA1; GeneC1 changed on an EC box.
    expect_equal(lab$label, c("GeneA2", "GeneC1"))
    expect_equal(lab$x, c(100, 300))
})

test_that("several changed members on one box are listed together", {
    lab <- de_gene_node_labels(nodes_fixture(), c("12", "13"), symbols_fixture)
    expect_equal(lab$label, "GeneA2/GeneA3")
})

test_that("a box already printing the changed gene is left alone on symbol maps", {
    lab <- de_gene_node_labels(nodes_fixture(), "21", symbols_fixture)
    expect_equal(nrow(lab), 0)

    # On an EC-number map the printed label is never a symbol, so label it.
    lab_ec <- de_gene_node_labels(nodes_fixture(), "21", symbols_fixture,
                                  skip_if_shown = FALSE)
    expect_equal(lab_ec$label, "GeneB")
})

test_that("a gene with no symbol falls back to its Entrez id", {
    lab <- de_gene_node_labels(nodes_fixture(), "32", c("32" = NA_character_))
    expect_equal(lab$label, "32")
})

test_that("nothing changed, or no node table, gives no labels", {
    expect_equal(nrow(de_gene_node_labels(nodes_fixture(), character(0),
                                          symbols_fixture)), 0)
    expect_equal(nrow(de_gene_node_labels(NULL, "11", symbols_fixture)), 0)
    expect_equal(nrow(de_gene_node_labels(nodes_fixture()[, c("x", "y")], "11",
                                          symbols_fixture)), 0)
})


# ---- maps already drawn by another renderer ----------------------------------

test_that("maps other renderers drew are found per set and contrast", {
    d <- withr::local_tempdir()
    file.create(file.path(d, c(
        "xyz00100.metab_top.png", "xyz00200.prot_top.png",
        "xyz00300.multi_ora_A_vs_B.png", "xyz00400.multi_ora_C_vs_D.png",
        "xyz00500.png", "xyz00600.gsea_pair_a_vs_b.png")))

    first <- .pathview_already_drawn(d, "xyz", "A_vs_B", first_contrast = TRUE)
    expect_setequal(names(first), c("00100", "00200", "00300"))
    expect_equal(first[["00300"]], "multi-omics maps")

    # The per-omics maps show the first contrast only; the blank KEGG template
    # and this renderer's own maps never count.
    other <- .pathview_already_drawn(d, "xyz", "C_vs_D", first_contrast = FALSE)
    expect_equal(names(other), "00400")
})
