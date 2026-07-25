# tests/testthat/test-metab-de-input-log.R
#
# chosen_norm = "none" lets an already-processed table be fed straight into the
# pipeline. When that table is ALREADY log-scaled (preprocessing.input_is_log =
# true), the non-limma DE methods must NOT compute fold change as
# log2(mean_B / mean_A) on a linear matrix — they must take the difference of
# group means on the log scale (the correct log2 fold change).
#
# These tests pin:
#   1. de_two_group()'s two FC modes directly (unit), and
#   2. run_metabolomics_de() routing FC to the log-difference path when
#      input_is_log = true (integration).
#
# All fixtures are synthetic (no filesystem, no real data).

# =============================================================================
# Unit: de_two_group FC modes
# =============================================================================
test_that("de_two_group: mat_for_fc = NULL -> log2FC is the difference of log means", {
    mat <- matrix(c(2, 2.2, 4, 4.2), nrow = 1,
                  dimnames = list("F1", c("S1", "S2", "S3", "S4")))
    condition <- factor(c("A", "A", "B", "B"))

    res <- de_two_group(mat, condition, "B - A", mat_for_fc = NULL,
                        test_fn = function(b, a) stats::t.test(b, a, var.equal = FALSE))

    # mean_B - mean_A = 4.1 - 2.1 = 2.0  (already the log2 fold change)
    expect_equal(res$logFC[1], 2.0)
})

test_that("de_two_group: mat_for_fc provided -> log2FC = log2(mean_B / mean_A) on that (linear) matrix", {
    mat <- matrix(c(2, 2.2, 4, 4.2), nrow = 1,          # log-scale test matrix
                  dimnames = list("F1", c("S1", "S2", "S3", "S4")))
    fc  <- matrix(c(100, 100, 110, 110), nrow = 1,      # linear-scale FC matrix
                  dimnames = list("F1", c("S1", "S2", "S3", "S4")))
    condition <- factor(c("A", "A", "B", "B"))

    res <- de_two_group(mat, condition, "B - A", mat_for_fc = fc,
                        test_fn = function(b, a) stats::t.test(b, a, var.equal = FALSE))

    expect_equal(res$logFC[1], log2(110 / 100))
})


# =============================================================================
# Integration: run_metabolomics_de honors preprocessing.input_is_log
# =============================================================================
make_de_fixture <- function() {
    samples <- c("S1", "S2", "S3", "S4")
    # expr_work / expr_log: log-scale test matrix (A ~ 2.1, B ~ 4.1 -> log2FC 2.0)
    W <- matrix(c(2, 2.2, 4, 4.2), nrow = 1, dimnames = list("F1", samples))
    # expr_filt: a DIFFERENT (linear) matrix -> log2(110/100) ~ 0.1375 if used
    F <- matrix(c(100, 100, 110, 110), nrow = 1, dimnames = list("F1", samples))
    meta <- data.frame(sample_id = samples,
                       group = c("A", "A", "B", "B"),
                       stringsAsFactors = FALSE)
    pre <- list(expr_work = W, expr_log = W, expr_filt = F, meta = meta,
                row_data = data.frame(feature_id = "F1", stringsAsFactors = FALSE),
                info = list(normalization = list(scaling = "none")))
    contrasts <- data.frame(Contrast_name = "B - A", Numerator = "B",
                            Denominator = "A", stringsAsFactors = FALSE)
    list(pre = pre, contrasts = contrasts, samples = samples)
}

make_de_config <- function(input_is_log, transform = "none") {
    list(modes = list(metabolomics = list(
        de            = list(method = "t_test", condition_column = "group"),
        effects       = list(samples = "sample_id"),
        preprocessing = list(chosen_norm = "none", input_is_log = input_is_log,
                             transform = transform),
        qc            = list()
    )))
}

test_that("run_metabolomics_de: input_is_log = true -> FC is the log-scale mean difference", {
    fx  <- make_de_fixture()
    cfg <- make_de_config(input_is_log = TRUE, transform = "none")

    res <- suppressMessages(suppressWarnings(
        run_metabolomics_de(fx$pre, cfg, fx$contrasts)
    ))

    lfc <- res$de_tables[["B - A"]]$logFC[1]
    # Difference of log means (2.0), NOT log2(110/100) from expr_filt.
    expect_equal(lfc, 2.0)
})

test_that("run_metabolomics_de: input_is_log = false (default) -> FC is the linear log2 ratio + caution warning", {
    fx  <- make_de_fixture()
    cfg <- make_de_config(input_is_log = FALSE, transform = "none")

    expect_warning(
        res <- suppressMessages(run_metabolomics_de(fx$pre, cfg, fx$contrasts)),
        "already be log-scaled|already-log|log-scaled"
    )

    lfc <- res$de_tables[["B - A"]]$logFC[1]
    # log2(mean_B / mean_A) on expr_filt = log2(110 / 100).
    expect_equal(lfc, log2(110 / 100))
})

test_that("run_metabolomics_de: input_is_log = true with transform != none warns about double-log", {
    fx  <- make_de_fixture()
    cfg <- make_de_config(input_is_log = TRUE, transform = "log2")

    expect_warning(
        suppressMessages(run_metabolomics_de(fx$pre, cfg, fx$contrasts)),
        "double-log"
    )
})
