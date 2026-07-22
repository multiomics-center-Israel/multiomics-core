# Tests for write_pca_companions().
#
# Regression tests for the bug where the proteomics report's "With Sample Names"
# and interactive PCA tabs were always empty: the QC module wrote only the plain
# PCA PNGs, never the `*_labeled.png` and `*_scores.csv` companions the report
# template looks for. write_pca_companions() must produce both, named exactly as
# the report expects, with a scores CSV that carries sample + PC + metadata.

# Synthetic, fully reproducible PCA input: 20 features x 6 samples, two groups.
make_pca_inputs <- function() {
    withr::with_seed(1, {
        expr <- matrix(rnorm(20 * 6), nrow = 20,
                       dimnames = list(paste0("F", 1:20), paste0("S", 1:6)))
    })
    meta <- data.frame(
        sample_id = paste0("S", 1:6),
        group     = rep(c("ctrl", "trt"), each = 3),
        stringsAsFactors = FALSE
    )
    cfg <- list(effects = list(samples = "sample_id", color = "group", shape = NULL))
    list(expr = expr, meta = meta, cfg = cfg)
}

test_that("write_pca_companions writes labeled PNG and scores CSV with report names", {
    inp <- make_pca_inputs()
    out_dir <- withr::local_tempdir()

    p <- qc_pca_scatter(inp$expr, inp$meta, inp$cfg, pcs = c(1, 2), out_file = NULL)
    written <- write_pca_companions(p, out_dir, pcs = c(1, 2))

    f_scores  <- file.path(out_dir, "PCA_PC1.vs.PC2_scores.csv")
    f_labeled <- file.path(out_dir, "PCA_PC1.vs.PC2_labeled.png")
    expect_setequal(written, c(f_scores, f_labeled))
    expect_true(file.exists(f_scores))
    expect_true(file.exists(f_labeled))

    # The interactive chunk needs sample, PC columns, and the color column.
    sc <- read.csv(f_scores, stringsAsFactors = FALSE)
    expect_true(all(c("sample", "PC1", "PC2", "group") %in% names(sc)))
    expect_equal(nrow(sc), 6)
})

test_that("write_pca_companions honors the requested PC pair in file names", {
    inp <- make_pca_inputs()
    out_dir <- withr::local_tempdir()

    p <- qc_pca_scatter(inp$expr, inp$meta, inp$cfg, pcs = c(1, 3), out_file = NULL)
    write_pca_companions(p, out_dir, pcs = c(1, 3))

    expect_true(file.exists(file.path(out_dir, "PCA_PC1.vs.PC3_scores.csv")))
    expect_true(file.exists(file.path(out_dir, "PCA_PC1.vs.PC3_labeled.png")))
    sc <- read.csv(file.path(out_dir, "PCA_PC1.vs.PC3_scores.csv"))
    expect_true(all(c("PC1", "PC3") %in% names(sc)))
})

test_that("write_pca_companions warns and no-ops on a plot without scores", {
    out_dir <- withr::local_tempdir()
    bare <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
        ggplot2::geom_point()
    expect_warning(res <- write_pca_companions(bare, out_dir, pcs = c(1, 2)))
    expect_length(res, 0)
    expect_length(list.files(out_dir), 0)
})
