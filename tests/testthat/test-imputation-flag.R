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

    flag <- setNames(fr$imputed_above_measured.S_vs_NS, fr$FeatureID)
    expect_identical(unname(flag[c("p_high", "p_low")]), c(TRUE, FALSE))
    expect_true(is.na(flag[["p_full"]]))
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
