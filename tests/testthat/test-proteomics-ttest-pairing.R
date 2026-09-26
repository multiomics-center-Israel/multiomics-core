# Tests for the paired-design guard in run_ttest_de()
# test-proteomics-ttest-pairing.R
#
# run_ttest_de() used to resolve the paired mode per contrast as
#   is_paired <- paired && !is.null(pairing_col) && pairing_col %in% colnames(meta)
# so a configured paired design whose pairing column was missing ran the
# ordinary unpaired test with no warning, while the run's Methods text still
# said "paired". Silently swapping the statistical test is not an acceptable
# fallback, so both failures now abort. resolve_de_block() already takes this
# line for de$block_col.

tt_meta <- function() {
    data.frame(
        SampleID  = c("S1", "S2", "S3", "S4", "S5", "S6"),
        Condition = rep(c("ctl", "trt"), times = 3),
        Subject   = rep(c("P1", "P2", "P3"), each = 2),
        stringsAsFactors = FALSE
    )
}

tt_expr <- function(meta) {
    # Deterministic values: the guard under test runs before any test statistic,
    # and a seeded draw would still add noise to a failure message for nothing.
    #
    # Columns alternate ctl, trt. Each row's within-subject differences are
    # deliberately non-constant -- with a constant difference t.test(paired =
    # TRUE) has zero variance and errors, which would make the paired run return
    # all-NA p-values and prove nothing about the mode that ran.
    matrix(
        c(10, 12, 11, 13, 12, 15,
          20, 19, 21, 20, 22, 23,
           5,  9,  6, 10,  7, 13),
        nrow = 3, byrow = TRUE,
        dimnames = list(c("P0001", "P0002", "P0003"), meta$SampleID)
    )
}

tt_prot_tbl <- function() {
    data.frame(
        Protein.Group = c("P0001", "P0002", "P0003"),
        Genes         = c("GENE1", "GENE2", "GENE3"),
        stringsAsFactors = FALSE
    )
}

tt_contrasts <- function() {
    data.frame(
        Contrast_name = "trt_vs_ctl",
        Factor        = "Condition",
        Numerator     = "trt",
        Denominator   = "ctl",
        stringsAsFactors = FALSE
    )
}

tt_cfg <- function(...) {
    de <- utils::modifyList(list(p_adjust_method = "BH"), list(...))
    list(modes = list(proteomics = list(
        effects     = list(samples = "SampleID", color = "Condition"),
        id_columns  = list(protein_id = "Protein.Group", protein_annot = "Genes"),
        de          = de
    )))
}

run_tt <- function(meta = tt_meta(), ...) {
    run_ttest_de(
        expr_imp     = tt_expr(meta),
        meta         = meta,
        contrasts_df = tt_contrasts(),
        prot_tbl     = tt_prot_tbl(),
        cfg          = tt_cfg(...)
    )
}

test_that("a paired design without a pairing column aborts", {
    expect_error(run_tt(paired = TRUE), "de$pairing_col is not set", fixed = TRUE)
    expect_error(run_tt(paired = TRUE, pairing_col = ""),
                 "de$pairing_col is not set", fixed = TRUE)
})

test_that("a pairing column absent from the metadata aborts", {
    expect_error(run_tt(paired = TRUE, pairing_col = "Donor"),
                 "not in the sample metadata", fixed = TRUE)
    # The message names the column asked for and what is actually there, so the
    # fix is obvious without opening the sample sheet.
    expect_error(run_tt(paired = TRUE, pairing_col = "Donor"), "Donor", fixed = TRUE)
    expect_error(run_tt(paired = TRUE, pairing_col = "Donor"), "Subject", fixed = TRUE)
})

test_that("a valid paired design runs and is not downgraded", {
    res <- run_tt(paired = TRUE, pairing_col = "Subject")
    expect_type(res, "list")
    expect_named(res$de_tables, "trt_vs_ctl")
    expect_equal(nrow(res$de_tables[["trt_vs_ctl"]]), 3L)

    # The paired test genuinely answers differently from the unpaired one on
    # this fixture, which is what made the silent downgrade worth refusing.
    paired_p <- res$de_tables[["trt_vs_ctl"]]$P.Value
    unpaired_p <- run_tt()$de_tables[["trt_vs_ctl"]]$P.Value
    expect_false(anyNA(paired_p))
    expect_false(anyNA(unpaired_p))
    expect_false(isTRUE(all.equal(paired_p, unpaired_p)))
})

test_that("an unpaired run ignores the pairing column entirely", {
    # paired unset means the pairing column is irrelevant, present or not, and
    # must not start raising errors.
    expect_no_error(run_tt(pairing_col = "Donor"))
    expect_no_error(run_tt(paired = FALSE, pairing_col = "Donor"))
})

test_that("the fewer-than-2-common-pairs guard still applies", {
    # Kept from the original implementation: a pairing column that is present
    # but shares fewer than two pairs across the two groups is still refused,
    # now on top of the up-front validation rather than instead of it.
    meta <- tt_meta()
    meta$Subject <- c("P1", "P1", "P2", "P9", "P3", "P8")
    expect_error(run_tt(meta = meta, paired = TRUE, pairing_col = "Subject"),
                 "fewer than 2 common pairs", fixed = TRUE)
})
