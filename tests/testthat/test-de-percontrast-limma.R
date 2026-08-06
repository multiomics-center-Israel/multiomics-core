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

test_that("min_val-imputed matrix (no NAs) is still filtered by the per-contrast floor", {
    # Mirrors the real Serge input: values are floored to a repeated minimum, and
    # the observed mask carries no NA information (all TRUE). The per-contrast floor
    # detection must still drop floor-in-both proteins.
    samples <- c("C1", "C2", "C3", "A1", "A2", "A3", "B1", "B2", "B3")
    groups  <- c("Ctrl", "Ctrl", "Ctrl", "A", "A", "A", "B", "B", "B")
    FL <- 6.0
    expr <- rbind(
        up_in_A   = c(10.0, 10.1, 9.9, 13.0, 13.1, 12.9, 10.0, 10.1,  9.9),
        onlyA     = c(  FL,   FL,  FL, 14.0, 14.1, 13.9,   FL,   FL,   FL),
        all_floor = c(  FL,   FL,  FL,   FL,   FL,   FL,   FL,   FL,   FL)
    )
    colnames(expr) <- samples
    observed <- matrix(TRUE, nrow(expr), ncol(expr), dimnames = dimnames(expr))

    meta <- data.frame(SampleName = samples, Group = groups, stringsAsFactors = FALSE)
    prot_tbl <- data.frame(
        Protein.Group = rownames(expr), Protein.Names = rownames(expr),
        Genes = rownames(expr), First.Protein.Description = rownames(expr),
        stringsAsFactors = FALSE
    )
    contrasts_df <- data.frame(
        Contrast_name = c("A_vs_Ctrl", "B_vs_Ctrl"), Factor = "Group",
        Numerator = c("A", "B"), Denominator = c("Ctrl", "Ctrl"),
        stringsAsFactors = FALSE
    )
    cfg <- make_percontrast_fixture()$cfg

    res <- suppressMessages(run_limma_percontrast_proteomics(
        expr_imp = expr, observed = observed, meta = meta,
        contrasts_df = contrasts_df, prot_tbl = prot_tbl, cfg = cfg
    ))
    de_A <- res$de_tables[["A_vs_Ctrl"]]
    de_B <- res$de_tables[["B_vs_Ctrl"]]

    # all_floor: at the floor in both groups of both contrasts -> NA everywhere.
    expect_true(is.na(de_A$logFC[de_A$FeatureID == "all_floor"]))
    expect_true(is.na(de_B$logFC[de_B$FeatureID == "all_floor"]))
    # onlyA: floor in Ctrl and B -> dropped in B_vs_Ctrl, tested (large +) in A_vs_Ctrl.
    expect_true(is.na(de_B$logFC[de_B$FeatureID == "onlyA"]))
    expect_false(is.na(de_A$logFC[de_A$FeatureID == "onlyA"]))
    expect_gt(de_A$logFC[de_A$FeatureID == "onlyA"], 3)
    # up_in_A: above floor in both groups -> tested in both contrasts.
    expect_false(is.na(de_A$logFC[de_A$FeatureID == "up_in_A"]))
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

test_that("a contrast Factor that is not the grouping column is rejected", {
    fx <- make_percontrast_fixture()
    bad <- fx$contrasts_df
    bad$Factor <- "SomeOtherColumn"   # configured group_col is "Group"

    expect_error(
        run_limma_percontrast_proteomics(
            expr_imp = fx$expr, observed = fx$observed, meta = fx$meta,
            contrasts_df = bad, prot_tbl = fx$prot_tbl, cfg = fx$cfg
        ),
        "Factor must equal the grouping"
    )

    # A missing required contrast column is also caught, not silently ignored.
    dropped <- fx$contrasts_df[, setdiff(names(fx$contrasts_df), "Denominator")]
    expect_error(
        run_limma_percontrast_proteomics(
            expr_imp = fx$expr, observed = fx$observed, meta = fx$meta,
            contrasts_df = dropped, prot_tbl = fx$prot_tbl, cfg = fx$cfg
        ),
        "missing required column"
    )
})

test_that("fdrtool_correction is applied per contrast when enabled", {
    skip_if_not_installed("fdrtool")

    # fdrtool estimates an empirical null from the moderated t-statistics, so it
    # needs a realistically sized protein set (the tiny shared fixture has too few
    # to fit). Build a 200-protein, 3-vs-3 matrix with a handful of true positives;
    # seeded so the fit - and thus the assertion - is deterministic.
    n_prot  <- 400L
    samples <- c("C1", "C2", "C3", "A1", "A2", "A3")
    groups  <- c("Ctrl", "Ctrl", "Ctrl", "A", "A", "A")
    expr <- withr::with_seed(42, {
        m <- matrix(rnorm(n_prot * length(samples), mean = 10, sd = 0.3),
                    nrow = n_prot, ncol = length(samples))
        m[1:30, 4:6] <- m[1:30, 4:6] + 2   # 30 proteins up in group A
        m
    })
    rownames(expr) <- paste0("P", seq_len(n_prot))
    colnames(expr) <- samples
    observed <- matrix(TRUE, n_prot, length(samples), dimnames = dimnames(expr))

    meta <- data.frame(SampleName = samples, Group = groups, stringsAsFactors = FALSE)
    prot_tbl <- data.frame(
        Protein.Group = rownames(expr), Protein.Names = rownames(expr),
        Genes = rownames(expr), First.Protein.Description = rownames(expr),
        stringsAsFactors = FALSE
    )
    contrasts_df <- data.frame(
        Contrast_name = "A_vs_Ctrl", Factor = "Group",
        Numerator = "A", Denominator = "Ctrl", stringsAsFactors = FALSE
    )
    cfg_plain <- make_percontrast_fixture()$cfg
    cfg_fdr   <- cfg_plain
    cfg_fdr$modes$proteomics$de$fdrtool_correction <- TRUE

    res_plain <- suppressMessages(run_limma_percontrast_proteomics(
        expr_imp = expr, observed = observed, meta = meta,
        contrasts_df = contrasts_df, prot_tbl = prot_tbl, cfg = cfg_plain
    ))
    res_fdr <- suppressMessages(suppressWarnings(run_limma_percontrast_proteomics(
        expr_imp = expr, observed = observed, meta = meta,
        contrasts_df = contrasts_df, prot_tbl = prot_tbl, cfg = cfg_fdr
    )))

    p_plain <- res_plain$de_tables[["A_vs_Ctrl"]]$P.Value
    p_fdr   <- res_fdr$de_tables[["A_vs_Ctrl"]]$P.Value

    # Same tested/untested rows, but the correction actually changes the p-values.
    expect_equal(is.na(p_plain), is.na(p_fdr))
    expect_false(any(is.na(p_fdr)))
    expect_false(isTRUE(all.equal(p_plain, p_fdr)))
})

test_that("fdrtool_correction errors early when the package is unavailable", {
    skip_if(requireNamespace("fdrtool", quietly = TRUE),
            "fdrtool is installed; cannot test the missing-package path")

    fx <- make_percontrast_fixture()
    cfg_fdr <- fx$cfg
    cfg_fdr$modes$proteomics$de$fdrtool_correction <- TRUE

    expect_error(
        run_limma_percontrast_proteomics(
            expr_imp = fx$expr, observed = fx$observed, meta = fx$meta,
            contrasts_df = fx$contrasts_df, prot_tbl = fx$prot_tbl, cfg = cfg_fdr
        ),
        "fdrtool"
    )
})

test_that("load_omics_inputs makes contrasts optional only for limma_percontrast + control_condition", {
    tmp_dir <- tempfile("percontrast_loader_")
    dir.create(file.path(tmp_dir, "data"), recursive = TRUE)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

    # Minimal on-disk inputs: a protein table and a metadata table, no contrasts file.
    utils::write.csv(
        data.frame(Protein.Group = c("P1", "P2"),
                   C1 = c(10, 11), A1 = c(12, 9),
                   check.names = FALSE),
        file.path(tmp_dir, "data", "protein.csv"), row.names = FALSE
    )
    utils::write.csv(
        data.frame(SampleName = c("C1", "A1"), Group = c("Ctrl", "A")),
        file.path(tmp_dir, "data", "meta.csv"), row.names = FALSE
    )

    base_cfg <- function() list(
        project = list(dir = tmp_dir),
        paths   = list(raw = "data"),
        modes   = list(proteomics = list(
            files   = list(protein = "protein.csv", metadata = "meta.csv"),
            effects = list(samples = "SampleName", color = "Group")
        ))
    )

    # Contrasts is still required for the default method -> hard error.
    cfg_default <- base_cfg()
    expect_error(
        load_omics_inputs(cfg_default, mode = "proteomics"),
        "contrasts"
    )

    # limma_percontrast + control_condition -> contrasts no longer required.
    cfg_pc <- base_cfg()
    cfg_pc$modes$proteomics$de <- list(method = "limma_percontrast",
                                       control_condition = "Ctrl")
    inputs <- load_omics_inputs(cfg_pc, mode = "proteomics")
    expect_null(inputs$contrasts)
    expect_true(!is.null(inputs$protein))

    # limma_percontrast WITHOUT control_condition -> contrasts required again.
    cfg_pc_no_ctrl <- base_cfg()
    cfg_pc_no_ctrl$modes$proteomics$de <- list(method = "limma_percontrast")
    expect_error(
        load_omics_inputs(cfg_pc_no_ctrl, mode = "proteomics"),
        "contrasts"
    )
})
