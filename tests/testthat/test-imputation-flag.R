# tests/testthat/test-imputation-flag.R
#
# imputed_above_measured.<contrast> (R/domain/proteomics/06_outputs_legacy.R):
# per contrast, TRUE when the imputed cells of the two groups average above the
# measured ones. Left-censored imputation assumes a missing value is low; a
# protein near the bottom of the sample distribution can be filled with values
# above what was measured, and then "not detected" reads as "higher". All values
# below are synthetic.

flag_meta <- function() {
    data.frame(SampleID  = c("S_1", "S_2", "NS_1", "NS_2"),
               condition = c("S", "S", "NS", "NS"),
               stringsAsFactors = FALSE)
}

flag_contrasts <- function() {
    data.frame(Contrast_name = "S_vs_NS", Factor = "condition",
               Numerator = "S", Denominator = "NS",
               stringsAsFactors = FALSE)
}

flag_config <- function() {
    list(modes = list(proteomics = list(
        de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5, use_adj_for_pass1 = TRUE),
        imputation = list(multi_imputation = FALSE, no_repetitions = 1, min_no_passed = 1),
        de_table = list(id_col = "FeatureID"),
        effects = list(samples = "SampleID")
    )))
}

# p_high : S never measured, filled at 18 while NS measured at 15 -> TRUE
# p_low  : one S cell filled at 9, below the measured 10-12        -> FALSE
# p_full : fully measured                                          -> NA
# p_none : nothing measured in either group                        -> NA
flag_matrices <- function() {
    ids <- c("p_high", "p_low", "p_full", "p_none")
    measured <- matrix(c(NA, NA, 15, 15,
                         10, NA, 12, 12,
                         20, 20, 21, 21,
                         NA, NA, NA, NA),
                       nrow = 4, byrow = TRUE,
                       dimnames = list(ids, flag_meta()$SampleID))
    model <- measured
    model["p_high", c("S_1", "S_2")] <- c(18, 18.4)
    model["p_low", "S_2"] <- 9
    model["p_none", ] <- c(17, 17.2, 17.1, 17.3)
    list(measured = measured, model = model)
}

test_that("imputed values above the measured ones are flagged, below are not", {
    m <- flag_matrices()
    pre <- list(expr_filt = m$measured, meta = flag_meta())
    out <- flag_imputed_above_measured_proteomics(pre, m$model, flag_contrasts(),
                                                  flag_config())

    expect_identical(colnames(out), "S_vs_NS")
    expect_identical(rownames(out), rownames(m$measured))
    expect_true(out["p_high", "S_vs_NS"])
    expect_false(out["p_low", "S_vs_NS"])
})

test_that("rows with nothing imputed or nothing measured are NA, not FALSE", {
    m <- flag_matrices()
    pre <- list(expr_filt = m$measured, meta = flag_meta())
    out <- flag_imputed_above_measured_proteomics(pre, m$model, flag_contrasts(),
                                                  flag_config())
    expect_true(is.na(out["p_full", "S_vs_NS"]))
    expect_true(is.na(out["p_none", "S_vs_NS"]))
})

test_that("only the contrast's own samples are pooled", {
    # A third group with a low measured value must not pull p_high's measured
    # mean down into a different answer for S_vs_NS.
    meta <- rbind(flag_meta(),
                  data.frame(SampleID = "X_1", condition = "X", stringsAsFactors = FALSE))
    measured <- cbind(flag_matrices()$measured, X_1 = c(30, 30, 30, 30))
    model <- cbind(flag_matrices()$model, X_1 = c(30, 30, 30, 30))
    contrasts <- rbind(flag_contrasts(),
                       data.frame(Contrast_name = "X_vs_NS", Factor = "condition",
                                  Numerator = "X", Denominator = "NS",
                                  stringsAsFactors = FALSE))
    out <- flag_imputed_above_measured_proteomics(
        list(expr_filt = measured, meta = meta), model, contrasts, flag_config())

    expect_true(out["p_high", "S_vs_NS"])
    # X_vs_NS holds no imputed cell for p_high (X measured, NS measured)
    expect_true(is.na(out["p_high", "X_vs_NS"]))
})

test_that("a model matrix missing a feature drops the flag with a warning", {
    m <- flag_matrices()
    pre <- list(expr_filt = m$measured, meta = flag_meta())
    expect_warning(
        out <- flag_imputed_above_measured_proteomics(pre, m$model[-1, ],
                                                      flag_contrasts(), flag_config()),
        "imputed_above_measured")
    expect_null(out)
})

flag_summary <- function() {
    runs <- list(list(S_vs_NS = data.frame(
        FeatureID     = c("p_high", "p_low", "p_full", "p_none"),
        Protein.Names = paste0("PROT", 1:4, "_TEST"),
        Genes         = paste0("gene", 1:4),
        First.Protein.Description = paste("test protein", 1:4),
        logFC     = c(3, -1.5, -1, 0.1),
        P.Value   = c(0.01, 0.2, 0.03, 0.9),
        adj.P.Val = c(0.04, 0.5, 0.06, 0.95),
        stringsAsFactors = FALSE
    )))
    summarize_limma_mult_imputation(runs, flag_config())
}

flag_final_results <- function() {
    m <- flag_matrices()
    # The pipeline filter never lets an all-missing row through, so the
    # end-to-end fixture leaves p_none out.
    m <- lapply(m, function(x) x[c("p_high", "p_low", "p_full"), , drop = FALSE])
    pre <- list(expr_filt = m$measured, expr_imp_single = m$model,
                meta = flag_meta(), row_data = NULL)
    sdf <- flag_summary()
    sdf <- sdf[match(rownames(m$measured), sdf$FeatureID), , drop = FALSE]
    build_final_results_proteomics(
        pre = pre, summary_df = sdf, contrasts_df = flag_contrasts(),
        feature_id_col = "FeatureID", config = flag_config(), expr_model = m$model)
}

test_that("final results carry the flag between padj and upDown", {
    fr <- flag_final_results()
    nm <- names(fr)
    expect_true("imputed_above_measured.S_vs_NS" %in% nm)
    expect_equal(nm[which(nm == "padj.imputs.S_vs_NS") + 1L],
                 "imputed_above_measured.S_vs_NS")
    expect_equal(nm[which(nm == "imputed_above_measured.S_vs_NS") + 1L],
                 "upDown.imputs.S_vs_NS")

    # The exported column is a readable label: the flag text or blank.
    flag <- setNames(fr$imputed_above_measured.S_vs_NS, fr$FeatureID)
    expect_identical(unname(flag["p_high"]), "FLAG: Imputation value higher than measured")
    expect_identical(unname(flag[c("p_low", "p_full")]), c("", ""))
})

test_that("the flag does not change the significance call", {
    fr <- flag_final_results()
    # p_high passes (adj p 0.04, |log2FC| 3) and stays "up" although flagged.
    expect_identical(fr$upDown.imputs.S_vs_NS[fr$FeatureID == "p_high"], "up")
})

test_that("the Excel writer classifies the flag as a DE statistic", {
    skip_if_not_installed("openxlsx")
    fr <- flag_final_results()
    out_dir <- withr::local_tempdir()
    files <- write_final_results_excels_legacy_generic(
        final_results = fr, config = flag_config(), out_dir = out_dir,
        mode = "proteomics", id_col = "FeatureID",
        expr_for_de = flag_matrices()$model[c("p_high", "p_low", "p_full"), ],
        with_cutoffs = FALSE
    )
    f_de <- files[length(files)]
    skip_if(is.na(f_de) || !file.exists(f_de), "DE workbook was not written")

    top <- openxlsx::readWorkbook(f_de, sheet = "Results", colNames = FALSE,
                                  rows = 1:3, skipEmptyRows = FALSE)
    hdr_row <- which(apply(top, 1, function(r) identical(as.character(r[[1]]), "FeatureID")))
    expect_length(hdr_row, 1L)
    hdr <- unlist(top[hdr_row, ], use.names = FALSE)

    # Matched DE columns get the short stat name; an unmatched one would keep
    # its full name and be moved next to the ID as an annotation.
    expect_true("imputed_above_measured" %in% hdr)
    expect_false(any(grepl("^imputed_above_measured\\.", hdr)))
    expect_equal(hdr[which(hdr == "padj") + 1L], "imputed_above_measured")
})

# =============================================================================
# de.exclude_imputed_above_measured: flagged rows left untested per contrast
# (exclude_imputed_above_measured_proteomics(), 05e_imputed_above_measured.R)
# =============================================================================

# Three groups, two contrasts. p_high is flagged in S_vs_NS (S filled at 18
# above the measured 15) but not in X_vs_NS, where nothing is imputed.
excl_fixture <- function() {
    meta <- rbind(flag_meta(),
                  data.frame(SampleID = c("X_1", "X_2"), condition = "X",
                             stringsAsFactors = FALSE))
    m <- flag_matrices()
    m <- lapply(m, function(x) x[c("p_high", "p_low", "p_full"), , drop = FALSE])
    x_cols <- matrix(c(30, 30, 25, 25, 22, 22), nrow = 3, byrow = TRUE,
                     dimnames = list(rownames(m$measured), c("X_1", "X_2")))
    measured <- cbind(m$measured, x_cols)
    model <- cbind(m$model, x_cols)
    contrasts <- rbind(flag_contrasts(),
                       data.frame(Contrast_name = "X_vs_NS", Factor = "condition",
                                  Numerator = "X", Denominator = "NS",
                                  stringsAsFactors = FALSE))
    list(pre = list(expr_filt = measured, meta = meta), model = model,
         contrasts = contrasts)
}

excl_run <- function(p) {
    tab <- function(cn, pv) data.frame(
        FeatureID = c("p_high", "p_low", "p_full"), Contrast = cn,
        logFC = c(3, -1.5, -1), P.Value = pv,
        adj.P.Val = p.adjust(pv, method = "BH"), stringsAsFactors = FALSE)
    list(S_vs_NS = tab("S_vs_NS", p$S_vs_NS), X_vs_NS = tab("X_vs_NS", p$X_vs_NS))
}

excl_p <- list(S_vs_NS = c(0.001, 0.02, 0.03), X_vs_NS = c(0.004, 0.01, 0.2))

test_that("a flagged row loses its p-values in its own contrast only", {
    fx <- excl_fixture()
    runs <- list(excl_run(excl_p))
    out <- exclude_imputed_above_measured_proteomics(
        runs, list(fx$model), fx$pre, fx$contrasts, flag_config())
    s <- out$runs_de_tables[[1]]$S_vs_NS
    x <- out$runs_de_tables[[1]]$X_vs_NS

    expect_true(is.na(s$P.Value[s$FeatureID == "p_high"]))
    expect_true(is.na(s$adj.P.Val[s$FeatureID == "p_high"]))
    # The fold change is still the fitted one.
    expect_identical(s$logFC, runs[[1]]$S_vs_NS$logFC)
    # Not flagged in X_vs_NS: that table is untouched.
    expect_identical(x, runs[[1]]$X_vs_NS)

    expect_identical(out$excluded$n_runs, 1L)
    expect_identical(unname(out$excluded$n_excluded["p_high", ]), c(1L, 0L))
    expect_identical(sum(out$excluded$n_excluded), 1L)
})

test_that("the adjustment is recomputed over the rows still tested", {
    fx <- excl_fixture()
    out <- exclude_imputed_above_measured_proteomics(
        list(excl_run(excl_p)), list(fx$model), fx$pre, fx$contrasts, flag_config())
    s <- out$runs_de_tables[[1]]$S_vs_NS
    kept <- s$FeatureID != "p_high"
    expect_equal(s$adj.P.Val[kept], p.adjust(excl_p$S_vs_NS[-1], method = "BH"))
    expect_equal(s$P.Value[kept], excl_p$S_vs_NS[-1])

    # The configured method is the one recomputed, not BH by default.
    cfg <- flag_config()
    cfg$modes$proteomics$de$p_adjust_method <- "bonferroni"
    out_b <- exclude_imputed_above_measured_proteomics(
        list(excl_run(excl_p)), list(fx$model), fx$pre, fx$contrasts, cfg)
    s_b <- out_b$runs_de_tables[[1]]$S_vs_NS
    expect_equal(s_b$adj.P.Val[kept], p.adjust(excl_p$S_vs_NS[-1], method = "bonferroni"))
})

test_that("an excluded row cannot pass, and its log2FC is still reported", {
    fx <- excl_fixture()
    out <- exclude_imputed_above_measured_proteomics(
        list(excl_run(excl_p)), list(fx$model), fx$pre, fx$contrasts, flag_config())
    before <- summarize_limma_mult_imputation(list(excl_run(excl_p)), flag_config())
    after <- summarize_limma_mult_imputation(out$runs_de_tables, flag_config())
    row <- function(df) df$FeatureID == "p_high"

    expect_identical(before$pass.imputs.S_vs_NS[row(before)], 1)
    expect_true(is.na(after$pass.imputs.S_vs_NS[row(after)]))
    expect_true(is.na(after$padj.imputs.S_vs_NS[row(after)]))
    expect_true(is.na(after$pvalue.imputs.S_vs_NS[row(after)]))
    expect_equal(after$log2FC.imputs.S_vs_NS, before$log2FC.imputs.S_vs_NS)
    # Still tested, and still passing, in the other contrast.
    expect_identical(after$pass.imputs.X_vs_NS[row(after)], 1)

    # final_results: blank upDown for the excluded row.
    # The fixture carries no annotation columns; their absence warns.
    fr <- suppressWarnings(build_final_results_proteomics(
        pre = list(expr_filt = fx$pre$expr_filt, expr_imp_single = fx$model,
                   meta = fx$pre$meta, row_data = NULL),
        summary_df = after, contrasts_df = fx$contrasts, feature_id_col = "FeatureID",
        config = flag_config(), expr_model = fx$model))
    expect_false(identical(fr$upDown.imputs.S_vs_NS[fr$FeatureID == "p_high"], "up"))
    expect_identical(fr$upDown.imputs.X_vs_NS[fr$FeatureID == "p_high"], "up")
})

test_that("under multiple imputation each run is judged on its own draw", {
    fx <- excl_fixture()
    # Run 2 filled p_high's S cells below the measured NS values: not flagged.
    model2 <- fx$model
    model2["p_high", c("S_1", "S_2")] <- c(13, 13.2)
    p2 <- list(S_vs_NS = c(0.003, 0.02, 0.03), X_vs_NS = excl_p$X_vs_NS)
    runs <- list(excl_run(excl_p), excl_run(p2))

    out <- exclude_imputed_above_measured_proteomics(
        runs, list(fx$model, model2), fx$pre, fx$contrasts, flag_config())
    expect_true(is.na(out$runs_de_tables[[1]]$S_vs_NS$P.Value[1]))
    expect_identical(out$runs_de_tables[[2]]$S_vs_NS, runs[[2]]$S_vs_NS)
    expect_identical(out$excluded$n_runs, 2L)
    expect_identical(out$excluded$n_excluded["p_high", "S_vs_NS"], 1L)

    # Applied before the summary: the quantile is over the run that tested it,
    # and the run that did not test it counts as not passing.
    cfg <- flag_config()
    cfg$modes$proteomics$imputation <- list(multi_imputation = TRUE,
                                            no_repetitions = 2, min_no_passed = 1)
    sm <- summarize_limma_mult_imputation(out$runs_de_tables, cfg)
    r <- sm$FeatureID == "p_high"
    expect_equal(sm$pvalue.imputs.S_vs_NS[r], 0.003)
    expect_identical(sm$sum.pass.S_vs_NS[r], 1)

    # Excluded in one of two runs: it keeps a summarised p-value, so enrichment
    # keeps it. Excluded in both: it goes.
    expect_length(imputed_above_measured_excluded_ids(out$excluded, "S_vs_NS"), 0L)
    both <- exclude_imputed_above_measured_proteomics(
        runs, list(fx$model, fx$model), fx$pre, fx$contrasts, flag_config())
    expect_identical(imputed_above_measured_excluded_ids(both$excluded, "S_vs_NS"), "p_high")
    expect_length(imputed_above_measured_excluded_ids(both$excluded, "X_vs_NS"), 0L)
    expect_length(imputed_above_measured_excluded_ids(NULL, "S_vs_NS"), 0L)
})

test_that("the exclusion stops rather than silently skipping when the flag cannot be computed", {
    fx <- excl_fixture()
    expect_error(suppressWarnings(exclude_imputed_above_measured_proteomics(
        list(excl_run(excl_p)), list(fx$model[-1, ]), fx$pre, fx$contrasts,
        flag_config())), "exclude_imputed_above_measured")
    expect_error(exclude_imputed_above_measured_proteomics(
        list(excl_run(excl_p)), list(fx$model, fx$model), fx$pre, fx$contrasts,
        flag_config()), "imputed matrices")
})

# --- mod_proteomics_de() end to end -------------------------------------------

# Low-abundance proteins measured in one group only: perseus-like imputation
# fills the other group from the low tail of each sample (around 16 here),
# above the measured value of 10, so they are flagged where that group is in
# the contrast.
excl_de_inputs <- function(exclude = NULL, n_rep = 3L) {
    samples <- c(paste0("A_", 1:3), paste0("B_", 1:3), paste0("C_", 1:3))
    group <- rep(c("A", "B", "C"), each = 3)
    expr <- withr::with_seed(11, matrix(rnorm(40 * 9, mean = 20, sd = 2), nrow = 40,
                                        dimnames = list(paste0("P", 1:40), samples)))
    expr["P1", ] <- c(10, 10.2, 9.9, NA, NA, NA, 10.1, 10, 10.3)
    expr["P2", ] <- c(10.4, 10.1, 10, NA, NA, NA, 9.8, 10.2, 10)
    expr["P3", ] <- c(NA, 21, 20.5, 22, 21.4, 21.9, 20, 20.2, 21)
    meta <- data.frame(SampleName = samples, Group = group, stringsAsFactors = FALSE)
    prot_tbl <- data.frame(Protein.Group = rownames(expr), Protein.Names = rownames(expr),
                           Genes = rownames(expr),
                           First.Protein.Description = rownames(expr),
                           stringsAsFactors = FALSE)
    contrasts <- data.frame(Contrast_name = c("B_vs_A", "C_vs_A"), Factor = "Group",
                            Numerator = c("B", "C"), Denominator = "A",
                            stringsAsFactors = FALSE)
    de <- list(method = "limma", p_adjust_method = "BH", p_cutoff = 0.05,
               linear_fc_cutoff = 1.5)
    de$exclude_imputed_above_measured <- exclude
    cfg <- list(params = list(seed = 42), modes = list(proteomics = list(
        effects = list(samples = "SampleName", color = "Group"),
        id_columns = list(protein_id = "Protein.Group", sample_col = "SampleName",
                          protein_annot = c("Protein.Group", "Protein.Names", "Genes",
                                            "First.Protein.Description")),
        de_table = list(id_col = "FeatureID"),
        de = de,
        imputation = list(method = "perseus_like", width = 0.3, downshift = 1.8,
                          multi_imputation = n_rep > 1L, no_repetitions = n_rep,
                          min_no_passed = min(2L, n_rep)),
        batch_correction = list(method = "none")
    )))
    pre <- list(expr_raw = expr, expr_filt = expr, expr_work = expr,
                expr_imp_single = expr, meta = meta, row_data = prot_tbl)
    list(pre = pre, inputs = list(contrasts = contrasts), config = cfg)
}

run_excl_de <- function(x) {
    suppressMessages(mod_proteomics_de(pre = x$pre, inputs = x$inputs, config = x$config))
}

test_that("with the option off the DE result is the object it was before the option", {
    skip_if_not_installed("limma")
    x <- excl_de_inputs()
    absent <- run_excl_de(x)
    off <- run_excl_de(excl_de_inputs(exclude = FALSE))

    # The pre-option path, written out: impute, fit every run, summarise.
    imputations <- make_imputations_proteomics(x$pre$expr_filt, x$config)
    runs <- run_limma_multimp(imputations, x$pre$meta, x$inputs$contrasts,
                              x$pre$row_data, x$config)
    runs_de_tables <- lapply(runs, function(r) r$de_tables)
    expected <- list(
        method = "limma", imputations = imputations, runs = runs,
        runs_de_tables = runs_de_tables,
        summary_df = summarize_limma_mult_imputation(runs_de_tables, x$config),
        de_model = runs[[1]]$fit2, batch_design = runs[[1]]$batch_design)

    expect_identical(absent, expected)
    expect_identical(off, expected)
    expect_null(off$imputed_above_measured_excluded)
})

test_that("with the option on, each run drops exactly the rows its own draw flags", {
    skip_if_not_installed("limma")
    off <- run_excl_de(excl_de_inputs(exclude = FALSE))
    x <- excl_de_inputs(exclude = TRUE)
    on <- run_excl_de(x)

    # Same imputations and the same fit: only the testing changed.
    expect_identical(on$imputations, off$imputations)
    expect_identical(on$de_model, off$de_model)

    n_hit <- 0L
    for (i in seq_along(on$runs_de_tables)) {
        flags <- flag_imputed_above_measured_proteomics(
            x$pre, on$imputations[[i]], x$inputs$contrasts, x$config)
        for (cn in names(on$runs_de_tables[[i]])) {
            t_on <- on$runs_de_tables[[i]][[cn]]
            t_off <- off$runs_de_tables[[i]][[cn]]
            hit <- flags[[cn]][match(t_on$FeatureID, rownames(flags))] %in% TRUE
            n_hit <- n_hit + sum(hit)
            expect_identical(is.na(t_on$P.Value), hit)
            expect_identical(t_on$logFC, t_off$logFC)
            expect_identical(t_on$t, t_off$t)
            expect_identical(t_on$P.Value[!hit], t_off$P.Value[!hit])
            expect_equal(t_on$adj.P.Val, p.adjust(t_on$P.Value, method = "BH"))
            expect_identical(on$runs[[i]]$de_tables[[cn]], t_on)
        }
    }
    # The fixture is built to flag P1 and P2 in B_vs_A; C_vs_A has nothing
    # imputed for them.
    expect_gt(n_hit, 0L)
    n_ex <- on$imputed_above_measured_excluded$n_excluded
    expect_identical(unname(n_ex[c("P1", "P2"), "B_vs_A"]), c(3L, 3L))
    expect_identical(unname(n_ex[c("P1", "P2"), "C_vs_A"]), c(0L, 0L))
    sm <- on$summary_df
    expect_true(all(is.na(sm$padj.imputs.B_vs_A[sm$FeatureID %in% c("P1", "P2")])))
    expect_false(anyNA(sm$padj.imputs.C_vs_A[sm$FeatureID %in% c("P1", "P2")]))
})

test_that("enrichment leaves the excluded rows out of that contrast's table only", {
    x <- excl_de_inputs(exclude = TRUE)
    on <- run_excl_de(x)

    # Capture what reaches run_pathway_analysis(): its fgsea ranking, ORA
    # significant lists and ORA background are all built from these tables.
    # Sourced functions live in the global environment, so stub by assignment.
    target <- environment(run_proteomics_pathway)
    real_rpa <- get("run_pathway_analysis", envir = target)
    real_prep <- get("prepare_proteomics_gene_sets", envir = target)
    seen <- NULL
    assign("run_pathway_analysis", function(de_tables, ...) {
        seen <<- de_tables
        list()
    }, envir = target)
    assign("prepare_proteomics_gene_sets", function(de_tables, config) {
        list(organism = "Homo sapiens", annotation = NULL,
             gene_sets = list(SET = list(S1 = paste0("P", 1:10))))
    }, envir = target)
    withr::defer({
        assign("run_pathway_analysis", real_rpa, envir = target)
        assign("prepare_proteomics_gene_sets", real_prep, envir = target)
    })

    cfg <- x$config
    cfg$modes$proteomics$annotation <- list(skip_annotation = TRUE)
    out_dir <- withr::local_tempdir("test-excl-pathway-")
    suppressMessages(suppressWarnings(run_proteomics_pathway(on, list(), cfg, out_dir)))

    expect_false(any(c("P1", "P2") %in% seen$B_vs_A$FeatureID))
    expect_true(all(c("P1", "P2") %in% seen$C_vs_A$FeatureID))
    expect_identical(nrow(seen$C_vs_A), 40L)
    expect_identical(nrow(seen$B_vs_A), 38L)

    # Option off: nothing is dropped.
    seen <- NULL
    off <- run_excl_de(excl_de_inputs(exclude = FALSE))
    suppressMessages(suppressWarnings(run_proteomics_pathway(off, list(), cfg, out_dir)))
    expect_identical(nrow(seen$B_vs_A), 40L)
})

# --- config, Methods text and Excel glossary -----------------------------------

test_that("one resolver: absent and false are off, other modes are always off", {
    expect_false(de_excludes_imputed_above_measured(NULL))
    expect_false(de_excludes_imputed_above_measured(list()))
    expect_false(de_excludes_imputed_above_measured(list(exclude_imputed_above_measured = FALSE)))
    expect_true(de_excludes_imputed_above_measured(list(exclude_imputed_above_measured = TRUE)))
    expect_false(de_excludes_imputed_above_measured(
        list(exclude_imputed_above_measured = TRUE), mode = "rna"))
})

test_that("the Excel glossary says flagged rows were excluded only when they were", {
    entry <- function(...) {
        g <- build_provenance_notes("proteomics", ...)$glossary
        g$Meaning[g$Column == "imputed_above_measured.<contrast>"]
    }
    off <- entry(multi_imputation = FALSE)
    expect_identical(off, entry(multi_imputation = FALSE,
                                exclude_imputed_above_measured = FALSE))
    expect_match(off, "not part of the significance call", fixed = TRUE)
    expect_no_match(off, "excluded")

    on1 <- entry(multi_imputation = FALSE, exclude_imputed_above_measured = TRUE)
    expect_match(on1, "excluded flagged rows from testing", fixed = TRUE)
    expect_no_match(on1, "not part of the significance call", fixed = TRUE)
    on_m <- entry(multi_imputation = TRUE, exclude_imputed_above_measured = TRUE)
    expect_match(on_m, "in each imputation run", fixed = TRUE)
})

test_that("the Methods text gains one sentence only when the option is on", {
    cfg <- excl_de_inputs()$config
    off <- build_proteomics_methods_text(cfg, de_method = "limma")$differential
    cfg$modes$proteomics$de$exclude_imputed_above_measured <- FALSE
    expect_identical(build_proteomics_methods_text(cfg, de_method = "limma")$differential, off)
    expect_no_match(off, "left untested", fixed = TRUE)

    cfg$modes$proteomics$de$exclude_imputed_above_measured <- TRUE
    on <- build_proteomics_methods_text(cfg, de_method = "limma")$differential
    expect_match(on, "left untested", fixed = TRUE)
    expect_match(on, "separately in each imputation run", fixed = TRUE)
    expect_identical(sub(" In each contrast, proteins whose imputed.*$", "", on), off)
})
