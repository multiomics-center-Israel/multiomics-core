# tests/testthat/test-proteomics-heatmap-annotation.R
#
# The proteomics heatmap sample annotation (effects$heatmap_annotations):
#   - order_annotation_levels(), build_proteomics_sample_annotation(),
#     build_sample_annotation_colors(), annotation_for_pheatmap(),
#     proteomics_heatmap_annotation()                              (07g)
#   - wrap_qc_heatmap() drawing the resolver's annotation        (08_qc.R)
#
# Synthetic data only.

test_that("level order: leading number, then no/yes, then first appearance", {
    expect_equal(levels(order_annotation_levels(c("24h", "0h", "4h", "24h", "0h"))),
                 c("0h", "4h", "24h"))
    expect_equal(levels(order_annotation_levels(c(24, 0, 4))), c("0", "4", "24"))
    expect_equal(levels(order_annotation_levels(c("yes", "no", "yes"))), c("no", "yes"))
    expect_equal(levels(order_annotation_levels(c("Present", "Absent"))), c("Absent", "Present"))
    expect_equal(levels(order_annotation_levels(c("ctrl", "trt", "ctrl"))), c("ctrl", "trt"))
    expect_equal(levels(order_annotation_levels(c("b", "a", NA))), c("b", "a"))
})

test_that("sample annotation follows effects$heatmap_annotations, with the EDA colours", {
    meta <- data.frame(SampleID = paste0("S", 1:6),
                       group = c("g24", "g24", "g0", "g0", "g4", "g4"),
                       time = c("24h", "24h", "0h", "0h", "4h", "4h"),
                       heme = c("yes", "no", "no", "no", "yes", "no"),
                       stringsAsFactors = FALSE)
    p_cfg <- list(effects = list(samples = "SampleID", color = "group",
                                 heatmap_annotations = list("time", "heme")))
    annot <- build_proteomics_sample_annotation(meta, p_cfg)
    expect_equal(colnames(annot), c("time", "heme"))
    expect_equal(rownames(annot), meta$SampleID)
    expect_equal(levels(annot$time), c("0h", "4h", "24h"))
    expect_equal(levels(annot$heme), c("no", "yes"))

    cols <- build_sample_annotation_colors(annot)
    expect_equal(cols$time, c("0h" = "#DADAEB", "4h" = "#9E9AC8", "24h" = "#54278F"))
    expect_equal(cols$heme, c(no = "#BDBDBD", yes = "#B8864B"))

    # A column without a matching palette is left to pheatmap's defaults.
    p_cfg$effects$heatmap_annotations <- list("group")
    expect_null(build_sample_annotation_colors(build_proteomics_sample_annotation(meta, p_cfg)))

    # Missing columns: warn and keep the rest; none left -> the colour column.
    p_cfg$effects$heatmap_annotations <- list("time", "missing_col")
    expect_warning(a2 <- build_proteomics_sample_annotation(meta, p_cfg), "missing_col")
    expect_equal(colnames(a2), "time")
    p_cfg$effects$heatmap_annotations <- list("missing_col")
    expect_warning(a3 <- build_proteomics_sample_annotation(meta, p_cfg), "colour column")
    expect_equal(a3, build_heatmap_annotation_col(meta, p_cfg))
})

test_that("a QC heatmap draws exactly the resolver's annotation when the key is set", {
    meta <- data.frame(SampleID = paste0("S", 1:6),
                       group = rep(c("g0", "g4", "g24"), each = 2),
                       time = rep(c("0h", "4h", "24h"), each = 2),
                       heme = c("no", "no", "no", "yes", "no", "yes"),
                       stringsAsFactors = FALSE)
    p_cfg <- list(effects = list(samples = "SampleID", color = "group",
                                 heatmap_annotations = list("time", "heme")))
    x <- withr::with_seed(3, matrix(rnorm(60), nrow = 10,
                                    dimnames = list(paste0("P", 1:10), rev(meta$SampleID))))

    hm_ann <- proteomics_heatmap_annotation(meta, p_cfg, colnames(x))
    # pheatmap draws the first column at the bottom and lists it last, so the
    # columns arrive reversed: time ends up on top and first in the legend.
    expected <- build_proteomics_sample_annotation(meta, p_cfg)[colnames(x), c("heme", "time"),
                                                                drop = FALSE]
    expect_equal(hm_ann$annot, expected)

    # Capture what wrap_qc_heatmap() hands to pheatmap.
    # The R/ files are sourced, not a package, so the stub replaces the binding
    # where wrap_qc_heatmap() looks it up and is restored afterwards.
    seen <- NULL
    env <- environment(wrap_qc_heatmap)
    original <- get("plot_heatmap_core", envir = env)
    withr::defer(assign("plot_heatmap_core", original, envir = env))
    assign("plot_heatmap_core", function(expr_mat, annotation_col = NULL, ...) {
        seen <<- annotation_col
        invisible(NULL)
    }, envir = env)
    wrap_qc_heatmap(x, meta, p_cfg, stage = "proteomics",
                    annot = hm_ann$annot, annot_colors = hm_ann$colors)
    expect_equal(seen, expected)

    # Without the key the helper returns NULL and callers keep their annotation,
    # and the default annotation is not reordered.
    default_annot <- data.frame(a = 1:2, b = 3:4)
    expect_identical(annotation_for_pheatmap(default_annot,
                                             list(effects = list(samples = "SampleID"))),
                     default_annot)
    p_cfg$effects$heatmap_annotations <- NULL
    expect_null(proteomics_heatmap_annotation(meta, p_cfg, colnames(x)))
})

test_that("without the key the annotation is the existing condition annotation", {
    meta <- data.frame(SampleID = paste0("S", 1:4), group = c("a", "a", "b", "b"),
                       stringsAsFactors = FALSE)
    p_cfg <- list(effects = list(samples = "SampleID", color = "group"))
    expect_equal(build_proteomics_sample_annotation(meta, p_cfg),
                 build_heatmap_annotation_col(meta, p_cfg))
})

