# tests/testthat/test-de-percontrast-limma.R
#
# Unit tests for the per-contrast (two-group vs control) limma path:
#   - auto_generate_control_contrasts()      (R/modules/proteomics/02_mod_de.R)
#   - run_limma_percontrast_proteomics()     (R/domain/proteomics/05_de_summary.R)
#
# Synthetic 3-condition design: control "Ctrl" + tests "A", "B", 3 replicates each.
# Values are log2 and fully deterministic (no RNG in the fit); only tiny fixed
# per-replicate offsets create within-group variance.

make_percontrast_fixture <- function() {
    samples <- c("C1", "C2", "C3", "A1", "A2", "A3", "B1", "B2", "B3")
    groups  <- c("Ctrl", "Ctrl", "Ctrl", "A", "A", "A", "B", "B", "B")

    # Proteins (rows). test_only_A is at the imputation floor in Ctrl and B.
    expr <- rbind(
        up_in_A     = c(10.00, 10.10,  9.90, 13.00, 13.10, 12.90, 10.00, 10.10,  9.90),
        test_only_A = c( 6.00,  6.05,  5.95, 14.00, 14.10, 13.90,  6.00,  6.05,  5.95),
        flat        = c(10.00, 10.10,  9.90, 10.00, 10.05,  9.95, 10.00, 10.10,  9.90),
        small_up_A  = c(10.00, 10.02,  9.98, 10.30, 10.32, 10.28, 10.00, 10.02,  9.98)
    )
    colnames(expr) <- samples

    # Observed pre-imputation mask: all TRUE except test_only_A in Ctrl and B.
    observed <- matrix(TRUE, nrow = nrow(expr), ncol = ncol(expr),
                       dimnames = dimnames(expr))
    observed["test_only_A", c("C1", "C2", "C3", "B1", "B2", "B3")] <- FALSE

    meta <- data.frame(
        SampleName = samples,
        Group      = groups,
        stringsAsFactors = FALSE
    )

    prot_tbl <- data.frame(
        Protein.Group             = rownames(expr),
        Protein.Names             = paste0(rownames(expr), "_HUMAN"),
        Genes                     = rownames(expr),
        First.Protein.Description = paste0("desc ", rownames(expr)),
        stringsAsFactors = FALSE
    )

    contrasts_df <- data.frame(
        Contrast_name = c("A_vs_Ctrl", "B_vs_Ctrl"),
        Factor        = "Group",
        Numerator     = c("A", "B"),
        Denominator   = c("Ctrl", "Ctrl"),
        stringsAsFactors = FALSE
    )

    cfg <- list(modes = list(proteomics = list(
        effects    = list(samples = "SampleName", color = "Group"),
        id_columns = list(
            protein_id    = "Protein.Group",
            sample_col    = "SampleName",
            protein_annot = c("Protein.Group", "Protein.Names", "Genes",
                              "First.Protein.Description")
        ),
        de_table   = list(id_col = "FeatureID"),
        de         = list(p_adjust_method = "BH", p_cutoff = 0.05,
                          linear_fc_cutoff = 1.5, use_adj_for_pass1 = TRUE),
        imputation = list(no_repetitions = 1, min_no_passed = 1,
                          multi_imputation = FALSE)
    )))

    list(expr = expr, observed = observed, meta = meta, prot_tbl = prot_tbl,
         contrasts_df = contrasts_df, cfg = cfg)
}

test_that("config validation accepts de.method = limma_percontrast", {
    cfg <- make_percontrast_fixture()$cfg$modes$proteomics
    cfg$de$method <- "limma_percontrast"
    cfg$de$control_condition <- "Ctrl"
    expect_error(validate_proteomics_config(cfg), NA)
})

test_that("auto_generate_control_contrasts references every non-control level", {
    meta <- make_percontrast_fixture()$meta
    cdf <- auto_generate_control_contrasts(meta, "Group", "Ctrl")

    expect_setequal(cdf$Contrast_name, c("A_vs_Ctrl", "B_vs_Ctrl"))
    expect_true(all(cdf$Denominator == "Ctrl"))
    expect_setequal(cdf$Numerator, c("A", "B"))
    expect_true(all(cdf$Factor == "Group"))
})

test_that("auto_generate_control_contrasts errors on an absent control", {
    meta <- make_percontrast_fixture()$meta
    expect_error(
        auto_generate_control_contrasts(meta, "Group", "Nonexistent"),
        "not a level"
    )
})

test_that("each per-contrast table is full-length and row-aligned to the matrix", {
    fx <- make_percontrast_fixture()
    res <- run_limma_percontrast_proteomics(
        expr_imp = fx$expr, observed = fx$observed, meta = fx$meta,
        contrasts_df = fx$contrasts_df, prot_tbl = fx$prot_tbl, cfg = fx$cfg
    )

    expect_named(res$de_tables, c("A_vs_Ctrl", "B_vs_Ctrl"))
    for (de in res$de_tables) {
        expect_equal(nrow(de), nrow(fx$expr))
        expect_equal(de$FeatureID, rownames(fx$expr))
    }
})

test_that("a test-only protein shows a large positive logFC (no compression)", {
    fx <- make_percontrast_fixture()
    res <- run_limma_percontrast_proteomics(
        expr_imp = fx$expr, observed = fx$observed, meta = fx$meta,
        contrasts_df = fx$contrasts_df, prot_tbl = fx$prot_tbl, cfg = fx$cfg
    )
    de_A <- res$de_tables[["A_vs_Ctrl"]]

    lfc <- de_A$logFC[de_A$FeatureID == "test_only_A"]
    expect_false(is.na(lfc))
    # Present in A, floor in Ctrl -> should be strongly positive, not compressed.
    expect_gt(lfc, 3)
})

test_that("a protein at the floor in both groups of a contrast is NA there, tested elsewhere", {
    fx <- make_percontrast_fixture()
    res <- run_limma_percontrast_proteomics(
        expr_imp = fx$expr, observed = fx$observed, meta = fx$meta,
        contrasts_df = fx$contrasts_df, prot_tbl = fx$prot_tbl, cfg = fx$cfg
    )
    de_A <- res$de_tables[["A_vs_Ctrl"]]
    de_B <- res$de_tables[["B_vs_Ctrl"]]

    # test_only_A: floor in both Ctrl and B -> dropped (NA) in B_vs_Ctrl ...
    expect_true(is.na(de_B$logFC[de_B$FeatureID == "test_only_A"]))
    expect_true(is.na(de_B$adj.P.Val[de_B$FeatureID == "test_only_A"]))
    # ... but observed in A -> tested (non-NA) in A_vs_Ctrl.
    expect_false(is.na(de_A$logFC[de_A$FeatureID == "test_only_A"]))
})

test_that("the summary honours linear_fc_cutoff = 1.5 (blocks a significant sub-1.5 protein)", {
    fx <- make_percontrast_fixture()
    res <- run_limma_percontrast_proteomics(
        expr_imp = fx$expr, observed = fx$observed, meta = fx$meta,
        contrasts_df = fx$contrasts_df, prot_tbl = fx$prot_tbl, cfg = fx$cfg
    )

    summ <- summarize_limma_mult_imputation(list(res$de_tables), fx$cfg)

    ca       <- normalize_contrast_name("A_vs_Ctrl")
    pass_col <- paste0("pass.imputs.", ca)
    padj_col <- paste0("padj.imputs.", ca)
    fc_col   <- paste0("linearFC.imputs.", ca)

    idx <- function(id) match(id, summ$FeatureID)

    # up_in_A: large FC + significant -> passes.
    expect_equal(summ[[pass_col]][idx("up_in_A")], 1)
    # flat: no change -> does not pass.
    expect_true(is.na(summ[[pass_col]][idx("flat")]))

    # small_up_A: tight, clearly significant, but linear FC < 1.5 -> blocked by the
    # fold-change gate, NOT by the p-value. Assert both facts to isolate the gate.
    expect_true(is.na(summ[[pass_col]][idx("small_up_A")]))
    expect_lt(summ[[padj_col]][idx("small_up_A")], 0.05)
    expect_lt(abs(summ[[fc_col]][idx("small_up_A")]), 1.5)
})
