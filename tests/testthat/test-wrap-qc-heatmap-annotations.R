# Which sample annotation bars wrap_qc_heatmap() draws. With annot_cols the
# clustered QC heatmap uses the configured heatmap_annotations, like its
# sibling heatmaps; without them it keeps the color + shape bars. Synthetic data.

heatmap_fixture <- function() {
    m <- matrix(c(1, 5, 2, 6, 3, 7, 4, 8, 9, 2, 8, 1), nrow = 3,
                dimnames = list(paste0("g", 1:3), paste0("S", 1:4)))
    meta <- data.frame(SampleID = paste0("S", 1:4), group = c("A", "A", "B", "B"),
                       treatment = c("x", "y", "x", "y"), replicate = c("1", "2", "1", "2"))
    cfg <- list(effects = list(samples = "SampleID", color = "group", shape = "replicate"))
    list(m = m, meta = meta, cfg = cfg)
}

# Names of the annotation bars pheatmap actually drew.
annotation_drawn <- function(...) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    ph <- suppressMessages(wrap_qc_heatmap(...))
    lay <- ph$gtable$layout
    ph$gtable$grobs[[which(lay$name == "col_annotation_names")]]$label
}

test_that("annot_cols sets the annotation bars", {
    f <- heatmap_fixture()
    skip_if_not_installed("pheatmap")
    expect_identical(annotation_drawn(f$m, f$meta, f$cfg, stage = "rna",
                                      annot_cols = c("group", "treatment")),
                     c("group", "treatment"))
})

test_that("without annot_cols the color and shape columns are used", {
    f <- heatmap_fixture()
    skip_if_not_installed("pheatmap")
    expect_identical(annotation_drawn(f$m, f$meta, f$cfg, stage = "rna"),
                     c("group", "replicate"))
})
