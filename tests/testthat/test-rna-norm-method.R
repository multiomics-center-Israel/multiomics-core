# tests/testthat/test-rna-norm-method.R
#
# The RNA QC plots, report, pipeline summary, PowerPoint and commentary used to
# label normalized values "TMMlogCPM" whatever ran. A tximport run cannot be
# TMMlogCPM at all (normalize_counts() refuses it), so a VST run carried TMM
# titles. Every reader-facing consumer now takes the method from
# resolve_rna_norm_method() / find_rna_norm_file() and the wording from
# describe_rna_norm_method().

synthetic_pre <- function(method = NULL, info_method = NULL) {
    mat <- matrix(c(1, 2, 3, 4), nrow = 2,
                  dimnames = list(c("G1", "G2"), c("S1", "S2")))
    if (!is.null(method)) attr(mat, "method") <- method
    list(expr_work = mat, info = list(expr_work_method = info_method))
}

test_that("resolve_rna_norm_method reads the method that ran, not the config", {
    expect_equal(resolve_rna_norm_method(synthetic_pre("VST")), "VST")
    expect_equal(resolve_rna_norm_method(synthetic_pre("TMMlogCPM")), "TMMlogCPM")
    expect_equal(resolve_rna_norm_method(synthetic_pre("preprocessed")), "preprocessed")
    # attribute wins over pre$info; pre$info is the fallback
    expect_equal(resolve_rna_norm_method(synthetic_pre("VST", "TMMlogCPM")), "VST")
    expect_equal(resolve_rna_norm_method(synthetic_pre(NULL, "VST")), "VST")
    expect_equal(resolve_rna_norm_method(synthetic_pre()), "unknown")
    expect_equal(resolve_rna_norm_method(NULL), "unknown")
})

test_that("describe_rna_norm_method gives method-specific wording", {
    vst <- describe_rna_norm_method("VST")
    expect_equal(vst$label, "VST")
    expect_match(vst$sentence, "variance-stabilizing", fixed = TRUE)
    expect_false(grepl("TMM", vst$sentence, fixed = TRUE))
    expect_false(grepl("CPM", vst$value_label, fixed = TRUE))
    expect_false(grepl("prior", paste(vst$tags, collapse = " "), fixed = TRUE))

    tmm <- describe_rna_norm_method("TMMlogCPM", prior_count = 2)
    expect_equal(tmm$label, "TMMlogCPM")
    expect_match(tmm$sentence, "prior count = 2", fixed = TRUE)

    pp <- describe_rna_norm_method("preprocessed")
    expect_match(pp$sentence, "no further normalization", fixed = TRUE)

    # Unrecognised or missing methods must not borrow another method's text
    for (m in list("unknown", NULL, NA_character_, "")) {
        d <- describe_rna_norm_method(m)
        expect_equal(d$label, "unknown")
        expect_false(grepl("TMM|VST|CPM", d$sentence))
        expect_equal(d$value_label, "Normalized expression")
    }
    other <- describe_rna_norm_method("quantile")
    expect_equal(other$label, "quantile")
    expect_false(grepl("TMM|VST|CPM", other$sentence))
})

test_that("find_rna_norm_file returns the single export and its method", {
    ds <- withr::local_tempdir()
    expect_null(find_rna_norm_file(ds))

    file.create(file.path(ds, "rna_norm_VST.tsv"))
    hit <- find_rna_norm_file(ds)
    expect_equal(hit$method, "VST")
    expect_equal(basename(hit$path), "rna_norm_VST.tsv")

    # Two candidates are ambiguous; do not guess which run wrote which
    file.create(file.path(ds, "rna_norm_TMMlogCPM.tsv"))
    expect_null(find_rna_norm_file(ds))
})

test_that("write_rnaseq_datasets_legacy removes another method's stale export", {
    out <- withr::local_tempdir()
    dirs <- create_legacy_output_dirs(out)
    file.create(file.path(dirs$datasets, "rna_norm_TMMlogCPM.tsv"))

    pre <- synthetic_pre("VST")
    pre$expr_filt <- pre$expr_work
    write_rnaseq_datasets_legacy(pre, config = list(), dirs = dirs)

    expect_equal(sort(list.files(dirs$datasets, pattern = "^rna_norm_")),
                 "rna_norm_VST.tsv")
    expect_equal(find_rna_norm_file(dirs$datasets)$method, "VST")
})

test_that("the RNA report no longer hard-codes TMMlogCPM outside the TMM branch", {
    candidates <- c(
        testthat::test_path("..", "..", "R", "domain", "rnaseq", "report_template.Rmd"),
        "R/domain/rnaseq/report_template.Rmd"
    )
    f <- candidates[file.exists(candidates)][1]
    if (is.na(f)) stop("R/domain/rnaseq/report_template.Rmd not found from the test directory")
    src <- readLines(f, warn = FALSE)

    expect_false(any(grepl("rna_norm_TMMlogCPM.tsv", src, fixed = TRUE)))
    expect_false(any(grepl("Expression (TMMlogCPM)", src, fixed = TRUE)))
    expect_false(any(grepl("normalization$method %||%", src, fixed = TRUE)))
})

test_that("QC plots carry the axis label they are given", {
    mat <- matrix(c(5, 6, 7, 8, 5.5, 6.5, 7.5, 8.5), nrow = 4,
                  dimnames = list(paste0("G", 1:4), c("S1", "S2")))
    meta <- data.frame(SampleID = c("S1", "S2"), Group = c("A", "B"))
    cfg <- list(effects = list(samples = "SampleID", color = "Group"))
    lbl <- describe_rna_norm_method("VST")$value_label

    p_dens <- qc_omic_density(mat, meta, cfg, x_label = lbl)
    expect_equal(p_dens$labels$x, lbl)
    # Default unchanged for the other omics
    expect_equal(qc_omic_density(mat, meta, cfg)$labels$x, "Intensity")

    p_box <- norm_boxplot(mat, meta, cfg, y_label = lbl)
    expect_equal(p_box$labels$y, lbl)
})
