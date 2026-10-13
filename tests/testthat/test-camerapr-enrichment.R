# tests/testthat/test-camerapr-enrichment.R
#
# cameraPR next to fgsea in run_pathway_analysis() (R/core/09_enrichment.R).
# Both must test the same sets with the same members, so a difference between
# their tables is a difference between the methods and not between inputs.

camera_fixture <- function() {
    withr::with_seed(3, {
        genes <- paste0("g", 1:200)
        stat <- rnorm(200)
        stat[1:15] <- stat[1:15] + 2
    })
    sets <- list(
        UP_SET    = c(genes[1:15], "not_ranked_1"),
        MIXED     = genes[40:70],
        SMALL     = genes[100:103],          # below min_size once intersected
        BIG       = genes[1:190],            # above max_size
        UNRANKED  = paste0("x", 1:20),       # no overlap at all
        EDGE_MIN  = c(genes[150:159], "not_ranked_2")  # exactly min_size after intersection
    )
    list(de = data.frame(FeatureID = genes, stat = stat, log2FoldChange = stat / 2,
                         stringsAsFactors = FALSE),
         gene_sets = list(CUSTOM = sets))
}

test_that("fgsea and cameraPR test identical sets with identical sizes", {
    skip_if_not_installed("fgsea")
    fx <- camera_fixture()
    res <- suppressMessages(run_pathway_analysis(
        list(c1 = fx$de), fx$gene_sets, method = c("fgsea", "camerapr"),
        min_size = 10, max_size = 100))
    fg <- res$c1$CUSTOM_fgsea
    cam <- res$c1$CUSTOM_camerapr
    expect_setequal(fg$pathway, cam$pathway)
    m <- match(fg$pathway, cam$pathway)
    expect_equal(as.integer(fg$size), cam$NGenes[m])
    expect_setequal(cam$pathway, c("UP_SET", "MIXED", "EDGE_MIN"))
    expect_true(all(c("pathway", "NGenes", "Direction", "PValue", "FDR", "inter_gene_cor")
                    %in% names(cam)))
    expect_equal(cam$Direction[cam$pathway == "UP_SET"], "Up")
})

test_that("non-finite ranks are dropped from both methods", {
    skip_if_not_installed("fgsea")
    fx <- camera_fixture()
    fx$de$stat[c(1, 2, 3)] <- c(NA, Inf, -Inf)
    res <- suppressMessages(run_pathway_analysis(
        list(c1 = fx$de), fx$gene_sets, method = c("fgsea", "camerapr"),
        min_size = 10, max_size = 100))
    fg <- res$c1$CUSTOM_fgsea
    cam <- res$c1$CUSTOM_camerapr
    # UP_SET had 15 ranked members; three of them are no longer finite.
    expect_equal(cam$NGenes[cam$pathway == "UP_SET"], 12L)
    expect_equal(as.integer(fg$size[fg$pathway == "UP_SET"]), 12L)
})

test_that("the set filter reproduces fgsea's size rule", {
    ranks <- setNames(c(3, 2, 1, 0, -1), c("a", "b", "c", "d", "e"))
    gs <- list(s1 = c("a", "a", "b", "z"), s2 = c("c", "d", "e"), s3 = "z")
    kept <- .filter_sets_to_ranks(gs, ranks, min_size = 2, max_size = 3)
    expect_named(kept, c("s1", "s2"))
    expect_equal(kept$s1, c("a", "b"))
})

test_that("cameraPR reports an FDR for a single set", {
    ranks <- setNames(c(5, 4, 3, -1, -2, -3, 0.5, -0.5), letters[1:8])
    out <- run_camerapr(ranks, list(only = c("a", "b", "c")))
    expect_equal(nrow(out), 1L)
    expect_true("FDR" %in% names(out))
    expect_equal(out$FDR, out$PValue)
    expect_equal(out$inter_gene_cor, 0.01)
})

test_that("an invalid inter-gene correlation is rejected", {
    ranks <- setNames(c(5, 4, 3, -1, -2), letters[1:5])
    sets <- list(s = c("a", "b"))
    expect_error(run_camerapr(ranks, sets, inter_gene_cor = 1), "strictly between -1 and 1")
    expect_error(run_camerapr(ranks, sets, inter_gene_cor = -1.5), "strictly between -1 and 1")
    expect_error(run_camerapr(ranks, sets, inter_gene_cor = NA_real_), "strictly between -1 and 1")
    expect_error(run_camerapr(ranks, sets, inter_gene_cor = c(0.1, 0.2)), "strictly between -1 and 1")
    expect_error(run_camerapr(c(ranks, f = NA), sets), "must be finite")
})

test_that("method = 'both' keeps its meaning and adds no cameraPR table", {
    skip_if_not_installed("fgsea")
    fx <- camera_fixture()
    res <- suppressMessages(run_pathway_analysis(
        list(c1 = fx$de), fx$gene_sets, method = "both", min_size = 10, max_size = 100))
    expect_false(any(grepl("_camerapr$", names(res$c1))))
    expect_true("CUSTOM_fgsea" %in% names(res$c1))
})

test_that("cameraPR results get a dotplot without an NES axis", {
    skip_if_not_installed("ggplot2")
    cam <- data.frame(pathway = paste0("S", 1:4), pathway_name = paste0("Set ", 1:4),
                      NGenes = c(10L, 12L, 15L, 20L), Direction = c("Up", "Down", "Up", "Down"),
                      PValue = c(0.001, 0.01, 0.2, 0.5), FDR = c(0.004, 0.02, 0.27, 0.5),
                      inter_gene_cor = 0.01, stringsAsFactors = FALSE)
    out_dir <- withr::local_tempdir("test-camerapr-plot-")
    files <- suppressMessages(generate_pathway_plots(list(c1 = list(DB_camerapr = cam)), out_dir))
    expect_length(files, 1L)
    expect_true(file.exists(files[[1]]))
})
