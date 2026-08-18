# tests/testthat/test-multi-ora-collections.R
#
# Unit tests for the collection-aware Multi-ORA plotting helpers in
# R/domain/multiomics/07b_multigsea_plots.R:
#   - classify_pathway_collection()
#   - select_top_ora_per_collection()
#   - summarize_ora_collections()
#   - plot_multi_ora_barplot()
#
# All data is synthetic: term ids follow the real namespaces (GO:#######,
# map#####, PF#####) but the p-values and gene ids are made up.

make_pooled_ora <- function(n_go = 100, n_kegg = 5, n_pfam = 2) {
    go <- data.frame(
        pathway = sprintf("go term %d", seq_len(n_go)),
        ID = sprintf("GO:%07d", seq_len(n_go)),
        # GO wins every head-to-head: this is the situation that made the plot
        # look GO-only.
        pvalue = seq_len(n_go) * 1e-8,
        stringsAsFactors = FALSE
    )
    kegg <- data.frame(
        pathway = sprintf("kegg map %d", seq_len(n_kegg)),
        ID = sprintf("map%05d", seq_len(n_kegg)),
        pvalue = seq_len(n_kegg) * 1e-3,
        stringsAsFactors = FALSE
    )
    pfam <- data.frame(
        pathway = sprintf("pfam domain %d", seq_len(n_pfam)),
        ID = sprintf("PF%05d", seq_len(n_pfam)),
        pvalue = seq_len(n_pfam) * 1e-2,
        stringsAsFactors = FALSE
    )
    rbind(go, kegg, pfam)
}

test_that("classify_pathway_collection recognises each gene-set namespace", {
    ids <- c("GO:0006955", "map00010", "ko00010", "hsa04110", "cel04110",
             "PF00001", "IPR000001", "weird_id", "", NA)
    expect_equal(
        classify_pathway_collection(ids),
        c("GO", "KEGG", "KEGG", "KEGG", "KEGG",
          "Pfam", "InterPro", "Other", "Other", "Other")
    )
})

test_that("select_top_ora_per_collection keeps every collection in the top slots", {
    ora <- make_pooled_ora()
    ora$collection <- classify_pathway_collection(ora$ID)

    top <- select_top_ora_per_collection(ora, top_n = 20)

    expect_equal(nrow(top), 20)
    expect_setequal(unique(top$collection), c("GO", "KEGG", "Pfam"))
    # Pfam only has 2 rows, KEGG 5; the rest of the slots go to GO.
    expect_equal(sum(top$collection == "Pfam"), 2)
    expect_equal(sum(top$collection == "KEGG"), 5)
    expect_equal(sum(top$collection == "GO"), 13)
    # Within a collection the strongest terms are the ones taken.
    expect_equal(top$ID[top$collection == "KEGG"], sprintf("map%05d", 1:5))
    # Rows come back in global p-value order.
    expect_false(is.unsorted(top$pvalue))
})

test_that("select_top_ora_per_collection is a plain top-n for a single collection", {
    ora <- make_pooled_ora(n_go = 40, n_kegg = 0, n_pfam = 0)
    ora$collection <- classify_pathway_collection(ora$ID)

    top <- select_top_ora_per_collection(ora, top_n = 20)

    expect_equal(top$ID, sprintf("GO:%07d", 1:20))
})

test_that("select_top_ora_per_collection returns every row when there are fewer than top_n", {
    ora <- make_pooled_ora(n_go = 3, n_kegg = 2, n_pfam = 1)
    ora <- ora[order(-ora$pvalue), ]
    ora$collection <- classify_pathway_collection(ora$ID)

    top <- select_top_ora_per_collection(ora, top_n = 20)

    expect_equal(nrow(top), 6)
    expect_false(is.unsorted(top$pvalue))
})

test_that("summarize_ora_collections separates tested-but-flat from never-tested", {
    t2g <- data.frame(
        term = c(rep("GO:0000001", 3), rep("map00010", 3),
                 rep("PF00001", 3), rep("PF00002", 3)),
        gene = paste0("g", 1:12),
        stringsAsFactors = FALSE
    )
    # PF00002 never reached the test (no hit gene / size filter); PF00001 was
    # tested and came out flat.
    tested <- data.frame(
        ID = c("GO:0000001", "map00010", "PF00001"),
        pvalue = c(1e-6, 1e-3, 0.8),
        stringsAsFactors = FALSE
    )
    pooled <- data.frame(
        pathway = c("go term", "kegg map"),
        ID = c("GO:0000001", "map00010"),
        pvalue = c(1e-6, 1e-3),
        stringsAsFactors = FALSE
    )
    attr(pooled, "all_tested") <- tested

    res <- summarize_ora_collections(t2g, pooled)

    expect_setequal(res$collection, c("GO", "KEGG", "Pfam"))
    pfam <- res[res$collection == "Pfam", ]
    expect_equal(pfam$n_sets, 2)
    expect_equal(pfam$n_tested, 1)
    expect_equal(pfam$n_p05, 0)
    expect_equal(pfam$n_reported, 0)
    expect_equal(res$n_reported[res$collection == "GO"], 1)
})

test_that("summarize_ora_collections tolerates a NULL ORA result", {
    t2g <- data.frame(term = c("GO:0000001", "map00010"),
                      gene = c("g1", "g2"), stringsAsFactors = FALSE)

    res <- summarize_ora_collections(t2g, NULL)

    expect_equal(sum(res$n_sets), 2)
    expect_equal(sum(res$n_tested), 0)
    expect_equal(sum(res$n_reported), 0)
})

test_that("plot_multi_ora_barplot draws every collection present in the table", {
    ora <- make_pooled_ora()
    png_path <- tempfile(fileext = ".png")
    grDevices::png(png_path, width = 1000, height = 700, res = 120)
    drawn <- withCallingHandlers(
        plot_multi_ora_barplot(ora, "synthetic pooled ORA"),
        warning = function(w) invokeRestart("muffleWarning")
    )
    grDevices::dev.off()

    expect_setequal(drawn, c("GO", "KEGG", "Pfam"))
    expect_true(file.exists(png_path))
    expect_gt(file.size(png_path), 0)
})

test_that("plot_multi_ora_barplot leaves a single-collection table untagged", {
    ora <- make_pooled_ora(n_go = 0, n_kegg = 8, n_pfam = 0)
    png_path <- tempfile(fileext = ".png")
    grDevices::png(png_path, width = 1000, height = 700, res = 120)
    drawn <- plot_multi_ora_barplot(ora, "synthetic compound ORA")
    grDevices::dev.off()

    expect_equal(drawn, "KEGG")
})

test_that("plot_multi_ora_barplot draws a searched-but-flat collection anyway", {
    # A collection whose sets were all tested and all came out flat has nothing
    # to draw; the run must still be able to say it was searched.
    ora <- make_pooled_ora(n_go = 5, n_kegg = 3, n_pfam = 0)
    attr(ora, "all_tested") <- rbind(
        ora[, c("ID", "pvalue")],
        data.frame(ID = sprintf("PF%05d", 1:4), pvalue = rep(0.9, 4),
                   stringsAsFactors = FALSE)
    )

    png_path <- tempfile(fileext = ".png")
    grDevices::png(png_path, width = 1000, height = 700, res = 120)
    drawn <- plot_multi_ora_barplot(ora, "synthetic pooled ORA")
    grDevices::dev.off()

    # Pfam is reported by the caption, not by a bar, so it stays out of the
    # drawn-collections return value.
    expect_setequal(drawn, c("GO", "KEGG"))
    expect_true(file.exists(png_path))
    expect_gt(file.size(png_path), 0)
})

test_that("plot_multi_ora_barplot tolerates a table with no all_tested attribute", {
    ora <- make_pooled_ora(n_go = 4, n_kegg = 2, n_pfam = 1)

    png_path <- tempfile(fileext = ".png")
    grDevices::png(png_path, width = 1000, height = 700, res = 120)
    drawn <- plot_multi_ora_barplot(ora, "synthetic pooled ORA")
    grDevices::dev.off()

    expect_setequal(drawn, c("GO", "KEGG", "Pfam"))
})
