# tests/testthat/test-de-summary-imput-padj-gate.R
#
# Regression test: the multi-imputation consensus call (pass.imputs) must stay
# consistent with the reported summary padj (padj.imputs). A feature may clear
# the min_no_passed vote yet have its q-quantile padj land just above the
# cutoff; in that case it must NOT be flagged significant, otherwise the volcano
# (which colours on padj.imputs) and the headline count disagree (the 48-vs-46
# bug). See R/domain/proteomics/05_de_summary.R.

# Build a runs_de_tables structure: a list of `n_rep` imputations, each a list
# keyed by contrast, each a data.frame with FeatureID/logFC/P.Value/adj.P.Val.
make_runs <- function(padj_by_feature, contrast = "A_vs_B") {
    ids <- names(padj_by_feature)
    n_rep <- length(padj_by_feature[[1]])
    lapply(seq_len(n_rep), function(n) {
        tbl <- data.frame(
            FeatureID = ids,
            logFC     = rep(1, length(ids)),          # |log2FC| >= 0, never binds
            P.Value   = vapply(ids, function(id) padj_by_feature[[id]][n], numeric(1)),
            adj.P.Val = vapply(ids, function(id) padj_by_feature[[id]][n], numeric(1)),
            stringsAsFactors = FALSE
        )
        setNames(list(tbl), contrast)
    })
}

make_cfg <- function(n_rep = 10, min_passed = 8) {
    list(modes = list(proteomics = list(
        de = list(use_fdr_for_pass1 = TRUE, p_cutoff = 0.05, linear_fc_cutoff = 1),
        imputation = list(multi_imputation = TRUE, no_repetitions = n_rep, min_no_passed = min_passed),
        de_table = list(id_col = "FeatureID")
    )))
}

test_that("consensus pass.imputs is gated on the reported padj.imputs", {
    # F1: significant in every imputation -> passes both vote and padj gate
    # F2: the real Scara3 padj vector -> 8/10 votes, but q=0.8 quantile = 0.0501
    # F3: never significant -> fails both
    padj_by_feature <- list(
        F1 = rep(0.001, 10),
        F2 = c(0.0025, 0.0034, 0.0101, 0.0132, 0.0169,
               0.0336, 0.0347, 0.0427, 0.0798, 0.0843),
        F3 = rep(0.5, 10)
    )

    out <- summarize_limma_mult_imputation(make_runs(padj_by_feature), make_cfg())

    pass_col <- grep("^pass\\.imputs\\.", names(out), value = TRUE)
    padj_col <- grep("^padj\\.imputs\\.", names(out), value = TRUE)
    expect_length(pass_col, 1)
    expect_length(padj_col, 1)

    pass <- setNames(out[[pass_col]], out$FeatureID)
    padj <- setNames(out[[padj_col]], out$FeatureID)

    # F2 clears 8/10 votes but its summary padj (0.0501) is over the cutoff,
    # so the gate must drop it.
    expect_gt(padj[["F2"]], 0.05)
    expect_true(is.na(pass[["F2"]]))

    # Sanity: the genuinely-significant feature survives, the null one does not.
    expect_equal(pass[["F1"]], 1)
    expect_true(is.na(pass[["F3"]]))

    # The invariant the fix guarantees: every called feature has padj <= cutoff.
    called <- !is.na(pass) & pass == 1
    expect_true(all(padj[called] <= 0.05))
})
