# Regression contract for the upDown.<contrast> column in Final_results.
#
# upDown is populated in build_final_results_generic() only when the per-contrast
# pass column (named by get_contrast_cols(mode)) is found in summary_df AND == 1.
# rnaseq and metabolomics share the same fc/p/padj/updown naming but their DE
# steps emit DIFFERENT pass affixes:
#     rnaseq DE  -> "<contrast>_pass"   (R/domain/rnaseq/04_de_summary.R)
#     metab  DE  -> "pass.<contrast>"   (R/domain/metabolomics/03_differential.R)
#     prot   DE  -> "pass.imputs.<cn>"  (R/domain/proteomics/05_de_summary.R)
# A single shared get_contrast_cols branch can only satisfy one of rnaseq/metab,
# which is exactly the see-saw that blanked metabolomics upDown. These tests lock
# the contract for ALL THREE omics so neither side can silently regress the other.
#
# Base-R only (no Bioconductor / openxlsx); run in RStudio from the repo root:
#   testthat::test_file("tests/testthat/test-updown-contract.R")
#
# Test map:
#   T-rna   baseline: rnaseq upDown populated (pass = <cn>_pass)
#   T-prot  baseline: proteomics upDown populated (pass = pass.imputs.<cn_norm>)
#   T-metab regression lock: metabolomics upDown populated (pass = pass.<cn>)
#           -> FAILS on the regressed code, PASSES after the fix
#   T-warn  missing expected pass column -> warning (not silent), upDown blank

`%||%` <- function(a, b) if (is.null(a)) b else a

find_repo_file <- function(rel) {
    dir <- normalizePath(".", mustWork = FALSE)
    for (i in 1:8) {
        cand <- file.path(dir, rel)
        if (file.exists(cand)) return(cand)
        dir <- dirname(dir)
    }
    stop("Could not locate ", rel, " from working dir ", getwd())
}

source(find_repo_file("R/core/05_export_excel.R"))
# Source the REAL normalize_contrast_name (used by get_contrast_cols' proteomics
# path AND by our proteomics fixture below) rather than stubbing it, so the
# baseline-protecting proteomics test can't drift from the real implementation.
# 01_io.R has no top-level loads, so it sources cleanly in base R.
source(find_repo_file("R/core/01_io.R"))


# ---- Fixture --------------------------------------------------------------
# Build a summary_df whose contrast column names mirror what each omics' DE step
# ACTUALLY emits (hardcoded here, NOT via get_contrast_cols), then run it through
# build_final_results_generic(mode). If get_contrast_cols(mode) looks under a
# different pass name than the DE emits, upDown comes out blank -> test fails.
#
# Three features: F1 passes with +FC (-> "up"), F2 passes with -FC (-> "down"),
# F3 does not pass (-> "").
make_result <- function(mode, pass_col, fc_col, p_col, padj_col, updown_col,
                         cn = "trt - ctrl", include_pass = TRUE) {
    ids <- c("F1", "F2", "F3")

    sdf <- data.frame(FeatureID = ids, stringsAsFactors = FALSE, check.names = FALSE)
    sdf[[fc_col]]   <- c(49, -3, 5)            # F2 negative keeps the FC-sign warning quiet
    sdf[[p_col]]    <- c(7.9e-11, 1e-4, 0.5)
    sdf[[padj_col]] <- c(9.5e-8, 1e-3, 0.6)
    if (include_pass) {
        sdf[[pass_col]] <- c(1, 1, NA)         # F1, F2 pass; F3 does not
    }
    sdf$pass_any_contrast <- c(1, 1, NA)       # present -> no pass_any fallback warning

    expr <- matrix(c(10, 20, 30, 40, 50, 60), nrow = 3, byrow = TRUE,
                   dimnames = list(ids, c("S1", "S2")))
    contrasts_df <- data.frame(
        Contrast_name = cn, Factor = "grp",
        Numerator = "trt", Denominator = "ctrl",
        stringsAsFactors = FALSE
    )

    res <- build_final_results_generic(
        summary_df     = sdf,
        expr_df        = expr,
        contrasts_df   = contrasts_df,
        feature_id_col = "FeatureID",
        mode           = mode
    )
    list(res = res, updown_col = updown_col, ids = ids)
}

expect_updown_contract <- function(out) {
    res <- out$res
    expect_true(out$updown_col %in% names(res))
    ud <- res[[out$updown_col]][match(out$ids, res$FeatureID)]
    expect_equal(ud[1], "up")     # F1: pass + positive FC
    expect_equal(ud[2], "down")   # F2: pass + negative FC
    expect_equal(ud[3], "")       # F3: did not pass
}


# ---- T-rna : baseline (must stay populated) -------------------------------
test_that("rnaseq upDown is populated (pass = <contrast>_pass)", {
    cn <- "trt - ctrl"
    out <- make_result(
        mode      = "rna",
        pass_col  = paste0(cn, "_pass"),       # rnaseq DE convention
        fc_col    = paste0("linearFC.", cn),
        p_col     = paste0("pvalue.", cn),
        padj_col  = paste0("padj.", cn),
        updown_col = paste0("upDown.", cn),
        cn = cn
    )
    expect_updown_contract(out)
})

# ---- T-prot : baseline (must stay populated) ------------------------------
test_that("proteomics upDown is populated (pass = pass.imputs.<cn_norm>)", {
    cn   <- "trt - ctrl"
    norm <- normalize_contrast_name(cn)        # real impl: strips spaces -> "trt-ctrl"
    out <- make_result(
        mode      = "proteomics",
        pass_col  = paste0("pass.imputs.", norm),   # proteomics DE convention
        fc_col    = paste0("linearFC.imputs.", norm),
        p_col     = paste0("pvalue.imputs.", norm),
        padj_col  = paste0("padj.imputs.", norm),
        updown_col = paste0("upDown.imputs.", norm),
        cn = cn
    )
    expect_updown_contract(out)
})

# ---- T-metab : the regression lock ----------------------------------------
test_that("metabolomics upDown is populated (pass = pass.<contrast>)", {
    cn <- "trt - ctrl"
    out <- make_result(
        mode      = "metabolomics",
        pass_col  = paste0("pass.", cn),       # metabolomics DE convention (03_differential.R)
        fc_col    = paste0("linearFC.", cn),
        p_col     = paste0("pvalue.", cn),
        padj_col  = paste0("padj.", cn),
        updown_col = paste0("upDown.", cn),
        cn = cn
    )
    expect_updown_contract(out)
})


# ---- T-warn : missing expected pass column warns (no silent all-NA) -------
test_that("a missing expected pass column warns and blanks upDown (not silent)", {
    cn <- "trt - ctrl"
    # mode = metabolomics expects pass.<cn>; we omit the pass column entirely.
    expect_warning(
        out <- make_result(
            mode      = "metabolomics",
            pass_col  = paste0("pass.", cn),   # name is irrelevant; include_pass=FALSE omits it
            fc_col    = paste0("linearFC.", cn),
            p_col     = paste0("pvalue.", cn),
            padj_col  = paste0("padj.", cn),
            updown_col = paste0("upDown.", cn),
            cn = cn,
            include_pass = FALSE
        ),
        "expected pass column"
    )
    ud <- out$res[[out$updown_col]]
    expect_true(all(ud == ""))
})
