# tests/testthat/test-proteomics-pass1-contract.R
#
# One answer to "did proteomics pass 1 read the adjusted p-value?" (#258).
# de_uses_adjusted_p() is the single resolver: the DE summariser, the Methods
# text and the Excel export all call it. The canonical key is use_adj_for_pass1;
# use_fdr_for_pass1 is a deprecated alias read the same way; absent means FALSE.
# The validation side (deprecation warning, conflict, rejected spellings, type)
# is in test-config-validation.R.

# ---- the resolver -----------------------------------------------------------

test_that("de_uses_adjusted_p resolves the proteomics pass-1 keys", {
    r <- function(...) de_uses_adjusted_p(list(...), "proteomics")
    expect_true(r(use_adj_for_pass1 = TRUE))
    expect_false(r(use_adj_for_pass1 = FALSE))
    expect_true(r(use_fdr_for_pass1 = TRUE))
    expect_false(r(use_fdr_for_pass1 = FALSE))
    expect_true(r(use_adj_for_pass1 = TRUE, use_fdr_for_pass1 = TRUE))
    # Absent keeps the existing default.
    expect_false(r(p_cutoff = 0.05))
    expect_false(de_uses_adjusted_p(NULL, "proteomics"))
    # The other modes' spellings are no longer read for proteomics; validation
    # refuses them, so a config cannot reach here with only those set.
    expect_false(r(use_adj = TRUE))
    expect_false(r(use_adjusted_pval = TRUE))
})

test_that("de_uses_adjusted_p is silent for proteomics however often it is called", {
    cfg <- list(use_fdr_for_pass1 = TRUE)
    expect_silent(for (i in 1:5) de_uses_adjusted_p(cfg, "proteomics"))
})

test_that("de_uses_adjusted_p is unchanged for RNA and metabolomics/lipidomics", {
    expect_true(de_uses_adjusted_p(list(use_adj = TRUE), "rna"))
    expect_true(de_uses_adjusted_p(list(p_cutoff = 0.05), "rna"))
    expect_warning(res <- de_uses_adjusted_p(list(use_adj = FALSE), "rna"), "gates on padj")
    expect_true(res)
    for (m in c("metabolomics", "lipidomics")) {
        expect_true(de_uses_adjusted_p(list(use_adjusted_pval = TRUE), m))
        expect_false(de_uses_adjusted_p(list(use_adjusted_pval = FALSE), m))
        expect_true(de_uses_adjusted_p(list(p_cutoff = 0.05), m))
        # Their fallback order still includes the neighbours' spellings.
        expect_false(de_uses_adjusted_p(list(use_adj = FALSE), m))
    }
})

# ---- the DE summariser, with multiple imputation ----------------------------

# Synthetic three-run table. F_split clears pass 1 on the raw p-value in runs 2
# and 3 (run 1 fails on fold change), but on the adjusted p-value only in run
# 2, so with min_no_passed = 2 the vote depends on the pass-1 choice. Its
# summary padj (the 2/3 quantile of 0.010, 0.040, 0.065, about 0.048) stays
# under the cutoff, so the final gate does not decide it either way. F_clear
# passes under both readings.
pass1_runs <- function(contrast = "A_vs_B") {
    p    <- list(F_split = c(0.002, 0.010, 0.030), F_clear = rep(0.001, 3))
    padj <- list(F_split = c(0.010, 0.040, 0.065), F_clear = rep(0.002, 3))
    lfc  <- list(F_split = c(0.2, 1.0, 1.0),       F_clear = rep(2, 3))
    ids <- names(p)
    lapply(1:3, function(n) {
        tbl <- data.frame(
            FeatureID = ids,
            logFC     = vapply(ids, function(id) lfc[[id]][n], numeric(1)),
            P.Value   = vapply(ids, function(id) p[[id]][n], numeric(1)),
            adj.P.Val = vapply(ids, function(id) padj[[id]][n], numeric(1)),
            stringsAsFactors = FALSE
        )
        setNames(list(tbl), contrast)
    })
}

pass1_cfg <- function(...) {
    list(modes = list(proteomics = list(
        de = c(list(method = "limma", p_cutoff = 0.05, linear_fc_cutoff = 1.5), list(...)),
        imputation = list(method = "perseus_like", width = 0.3, downshift = 1.8,
                          multi_imputation = TRUE, no_repetitions = 3, min_no_passed = 2),
        de_table = list(id_col = "FeatureID")
    )))
}

summarise <- function(...) summarize_limma_mult_imputation(pass1_runs(), pass1_cfg(...))

pass_of <- function(out, id) {
    col <- grep("^pass\\.imputs\\.", names(out), value = TRUE)
    stopifnot(length(col) == 1)
    out[[col]][out$FeatureID == id]
}

test_that("the alias gives the same pass flags and statistics as the canonical key", {
    canon_t <- summarise(use_adj_for_pass1 = TRUE)
    alias_t <- summarise(use_fdr_for_pass1 = TRUE)
    both_t  <- summarise(use_adj_for_pass1 = TRUE, use_fdr_for_pass1 = TRUE)
    expect_identical(alias_t, canon_t)
    expect_identical(both_t, canon_t)
    expect_identical(summarise(use_fdr_for_pass1 = FALSE), summarise(use_adj_for_pass1 = FALSE))

    # The fixture is not vacuous: the choice decides F_split's call.
    expect_true(is.na(pass_of(canon_t, "F_split")))
    expect_equal(pass_of(summarise(use_adj_for_pass1 = FALSE), "F_split"), 1)
    expect_equal(pass_of(canon_t, "F_clear"), 1)
})

test_that("an absent key keeps the raw-p default", {
    expect_identical(summarise(), summarise(use_adj_for_pass1 = FALSE))
})

# ---- the Methods text --------------------------------------------------------

methods_for <- function(...) {
    cfg <- pass1_cfg(...)
    cfg$modes$proteomics <- utils::modifyList(cfg$modes$proteomics, list(
        engine = "DIANN", scale_in = "linear",
        filtering = list(remove_contaminants = TRUE, contaminant_prefix = "cRAP-",
                         min_count = 3, min_groups = 1),
        normalization = list(method = "none"),
        batch_correction = list(method = "none"),
        clustering = list(enabled = FALSE)
    ))
    build_proteomics_methods_text(cfg, de_method = "limma")
}

test_that("the Methods text states the pass-1 choice the summariser applied", {
    expect_identical(methods_for(use_fdr_for_pass1 = TRUE), methods_for(use_adj_for_pass1 = TRUE))
    expect_match(methods_for(use_fdr_for_pass1 = TRUE)$consensus,
                 "its adjusted p-value met the cutoff", fixed = TRUE)
    expect_match(methods_for()$consensus, "its unadjusted p-value met the cutoff", fixed = TRUE)
})

# ---- every consumer agrees ---------------------------------------------------

test_that("summariser, Excel resolver and Methods text agree for every config state", {
    states <- list(
        canonical_true  = list(use_adj_for_pass1 = TRUE),
        canonical_false = list(use_adj_for_pass1 = FALSE),
        alias_true      = list(use_fdr_for_pass1 = TRUE),
        alias_false     = list(use_fdr_for_pass1 = FALSE),
        both_true       = list(use_adj_for_pass1 = TRUE, use_fdr_for_pass1 = TRUE),
        absent          = list()
    )
    for (nm in names(states)) {
        args <- states[[nm]]
        excel_adj <- de_uses_adjusted_p(do.call(pass1_cfg, args)$modes$proteomics$de,
                                        "proteomics")
        summary_adj <- is.na(pass_of(do.call(summarise, args), "F_split"))
        text <- do.call(methods_for, args)$consensus
        text_adj <- grepl("its adjusted p-value met the cutoff", text, fixed = TRUE)
        expect_identical(summary_adj, excel_adj, info = nm)
        expect_identical(text_adj, excel_adj, info = nm)
    }
})

# ---- single source -----------------------------------------------------------

test_that("no proteomics consumer reads a pass-1 spelling directly", {
    root <- normalizePath(testthat::test_path("..", ".."))
    files <- list.files(file.path(root, c("R/domain/proteomics", "R/modules/proteomics")),
                        pattern = "\\.R$", full.names = TRUE)
    # Validation is where the spellings are checked by design.
    files <- files[basename(files) != "90_config_validate.R"]
    pat <- paste0("(\\$|\\[\\[\\s*[\"'])",
                  "use_(adj_for_pass1|fdr_for_pass1|adj|adjusted_pval)\\b")
    hits <- character(0)
    for (f in files) {
        code <- sub("#.*$", "", readLines(f, warn = FALSE))
        bad <- grep(pat, code, perl = TRUE)
        if (length(bad)) hits <- c(hits, paste0(basename(f), ":", bad))
    }
    expect_identical(hits, character(0))
})
