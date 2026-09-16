# Regression tests: the fold change and the p-value must come from ONE matrix.
#
# Metabolomics and lipidomics DE used to run the test on the normalised,
# transformed, pre-scaling matrix (pre$expr_log) while taking the fold change
# from the raw filtered matrix (pre$expr_filt) via the mat_for_fc argument.
# Those two matrices differ by the sample normalisation, so the reported
# direction could contradict the test that produced the p-value -- on the Elah
# Pick lanes it did so for 72 of 575 lipids (9 at padj <= 0.05) and 9 of 662
# metabolites.
#
# This is the assertion test-proda-de.R already makes for proteomics ("reporting
# a crude group-mean ratio instead of the fitted coefficient"), carried over to
# the two domains that actually had the defect.

# The fixture is built so the two matrices disagree on EVERY feature: the
# normalised data halve in group B, while a 4x sample-loading inflation makes
# the raw intensities double in group B. Old behaviour: logFC = +1 everywhere.
# Correct behaviour: logFC = -1 everywhere, agreeing with the test statistic.
make_fc_fixture <- function(transform = "glog10", n_feat = 40, seed = 11) {
    set.seed(seed)
    samples <- c(sprintf("A%d", 1:5), sprintf("B%d", 1:5))
    meta <- data.frame(sample_id = samples,
                       Group = rep(c("A", "B"), each = 5),
                       stringsAsFactors = FALSE)
    idx_B <- which(meta$Group == "B")

    # Sample-normalised linear intensities; group B is EXACTLY 2x LOWER, so the
    # true fold change is exactly log2(0.5) = -1 and the unit conversion can be
    # asserted without sampling slack.
    idx_A <- which(meta$Group == "A")
    lin <- matrix(stats::rlnorm(n_feat * length(samples), meanlog = 12, sdlog = 0.1),
                  nrow = n_feat,
                  dimnames = list(sprintf("F%03d", seq_len(n_feat)), samples))
    lin[, idx_B] <- lin[, idx_A] * 0.5

    # Raw, un-normalised intensities: a 4x sample loading in group B flips the
    # crude ratio to 2x HIGHER.
    raw <- lin
    raw[, idx_B] <- raw[, idx_B] * 4

    expr_log <- transform_metab(lin, method = transform, pseudocount = 1)
    list(meta = meta, raw = raw, expr_log = expr_log,
         expr_work = scale_metab(expr_log, method = "auto"),
         transform = transform, idx_A = idx_A, idx_B = idx_B)
}

metab_pre <- function(f) {
    list(expr_raw = f$raw, expr_filt = f$raw, expr_log = f$expr_log,
         expr_work = f$expr_work, meta = f$meta, row_data = NULL,
         info = list(normalization = list(sample_norm = "none",
                                          transform = f$transform,
                                          scaling = "auto")))
}

metab_cfg <- function(f, method) {
    list(modes = list(metabolomics = list(
        effects       = list(samples = "sample_id", color = "Group"),
        preprocessing = list(transform = f$transform, scaling = "auto"),
        de            = list(method = method, p_cutoff = 0.05, linear_fc_cutoff = 1)
    )))
}

lipid_cfg <- function(f, method) {
    list(modes = list(lipidomics = list(
        effects       = list(samples = "sample_id", color = "Group"),
        normalization = list(transform = f$transform, scaling = "auto"),
        de            = list(method = method, p_cutoff = 0.05, linear_fc_cutoff = 1,
                             contrasts = "B - A")
    )))
}

ctr_table <- data.frame(Contrast_name = "B_vs_A", Numerator = "B",
                        Denominator = "A", stringsAsFactors = FALSE)

# The invariant, stated once: the sign of the reported fold change is the sign
# of the group difference on the matrix that was tested.
expect_fc_from_tested_matrix <- function(de, f) {
    tested_diff <- rowMeans(f$expr_log[, f$idx_B, drop = FALSE]) -
                   rowMeans(f$expr_log[, f$idx_A, drop = FALSE])
    tested_diff <- tested_diff[de$feature_id]
    expect_equal(sign(de$logFC), unname(sign(tested_diff)))
    # and it is NOT the raw-matrix ratio, which points the other way here
    expect_true(all(sign(de$logFC) != sign(log2(
        rowMeans(f$raw[de$feature_id, f$idx_B, drop = FALSE]) /
        rowMeans(f$raw[de$feature_id, f$idx_A, drop = FALSE])))))
}

for (m in c("limma", "t_test", "t_test_equal", "wilcoxon")) {
    test_that(paste0("metabolomics DE [", m,
                     "]: logFC comes from the tested matrix, not expr_filt"), {
        f <- make_fc_fixture()
        res <- suppressMessages(
            run_metabolomics_de(metab_pre(f), metab_cfg(f, m), ctr_table))
        de <- res$de_tables[["B_vs_A"]]
        expect_fc_from_tested_matrix(de, f)
        # glog10 differences are reported in log2 units: a true 2x decrease
        # reads as -1, not as the raw glog10 difference of -0.301. The residual
        # is glog10's own curvature near its offset, not a unit error.
        expect_equal(de$logFC, rep(-1, nrow(de)), tolerance = 5e-3)
    })

    test_that(paste0("lipidomics DE [", m,
                     "]: logFC comes from the tested matrix, not expr_filt"), {
        f <- make_fc_fixture(transform = "log2")
        pre <- metab_pre(f)
        pre$info$normalization$transform <- "log2"
        res <- suppressMessages(
            run_lipidomics_de(pre, lipid_cfg(f, m)))
        de <- res$de_tables[["B_vs_A"]]
        expect_fc_from_tested_matrix(de, f)
        # log2 is already log2 units; the residual is the +1 pseudocount.
        expect_equal(de$logFC, rep(-1, nrow(de)), tolerance = 5e-3)
    })
}

test_that("metabolomics/lipidomics DE: sign(logFC) agrees with sign(statistic)", {
    # Wilcoxon's W is unsigned, so this half of the invariant covers the three
    # methods whose statistic carries a direction.
    for (m in c("limma", "t_test", "t_test_equal")) {
        f <- make_fc_fixture()
        de <- suppressMessages(
            run_metabolomics_de(metab_pre(f), metab_cfg(f, m), ctr_table)
        )$de_tables[["B_vs_A"]]
        expect_equal(sign(de$logFC), sign(de$statistic),
                     info = paste("metabolomics", m))
    }
})

test_that("fc_to_log2_units converts each transform into log2 units", {
    f <- make_fc_fixture()
    cond <- factor(f$meta$Group)
    # glog10 / log10: a difference of d decades is d * log2(10) in log2 units
    tbl <- data.frame(feature_id = rownames(f$expr_log),
                      logFC = -0.3010, AveExpr = 12)
    for (tf in c("glog10", "log10")) {
        out <- fc_to_log2_units(tbl, tf, f$expr_log, cond, "B - A")
        expect_equal(out$logFC, rep(-1, nrow(tbl)), tolerance = 1e-3)
        expect_equal(out$AveExpr, rep(12 * log2(10), nrow(tbl)), tolerance = 1e-9)
    }
    # log2: already log2, left alone
    expect_equal(fc_to_log2_units(tbl, "log2", f$expr_log, cond, "B - A")$logFC,
                 rep(-0.3010, nrow(tbl)))
    # none: linear matrix, so the log2 ratio of the group means is taken on it
    lin <- f$expr_log  # stand-in linear matrix; only its group means matter
    out <- fc_to_log2_units(tbl, "none", exp(lin), cond, "B - A")
    expect_equal(unname(sign(out$logFC)),
                 unname(sign(rowMeans(exp(lin)[, f$idx_B, drop = FALSE]) -
                             rowMeans(exp(lin)[, f$idx_A, drop = FALSE]))))
})

test_that("metabolomics DE: the chosen_norm='none' scale caveat fires only when the transform says 'none'", {
    # The mat_raw block this commit removes carried a warning: with
    # chosen_norm = "none" the non-limma methods computed the fold change as a
    # linear ratio of the RAW filtered matrix, which is wrong if that table
    # already arrives log-scaled. Taking the fold change from the tested matrix
    # in its declared unit answers that for every transform except one --
    # "none", where the label means "no transform applied" but
    # fc_to_log2_units() reads it as "the values are linear". The warning is
    # kept for exactly that lane, and now covers limma too.
    warnings_of <- function(expr) {
        ws <- character()
        withCallingHandlers(
            suppressMessages(expr),
            warning = function(w) {
                ws <<- c(ws, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
        ws
    }

    f   <- make_fc_fixture(transform = "log2")
    pre <- metab_pre(f)
    pre$info$normalization$transform <- "log2"
    cfg <- metab_cfg(f, "t_test_equal")
    cfg$modes$metabolomics$preprocessing$chosen_norm <- "none"

    ws <- warnings_of(run_metabolomics_de(pre, cfg, ctr_table))
    expect_false(any(grepl("chosen_norm = 'none'", ws, fixed = TRUE)))

    pre$info$normalization$transform              <- "none"
    cfg$modes$metabolomics$preprocessing$transform <- "none"
    for (m in c("limma", "t_test", "t_test_equal", "wilcoxon")) {
        cfg$modes$metabolomics$de$method <- m
        ws <- warnings_of(run_metabolomics_de(pre, cfg, ctr_table))
        expect_true(any(grepl("chosen_norm = 'none' with transform = 'none'",
                              ws, fixed = TRUE)),
                    info = m)
    }
})
