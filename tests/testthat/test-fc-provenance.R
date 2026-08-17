# Tests for the fold-change provenance columns in Final_results_{ALL,DE}.
#
# Context: the per-sample columns are raw counts (RNA-seq) or unimputed log2
# intensities (proteomics), while the reported fold change is estimated on the
# normalized/imputed matrix. Readers who eyeball the sample columns therefore
# compute a different number. These tests pin the columns that close that gap:
# log2FC next to linearFC, a <sample>.norm block, and Mean.<group>.
#
# Deliberately BASE-R ONLY (no Bioconductor, no openxlsx) so they can be run
# from RStudio on Windows before pushing:
#
#   testthat::test_file("tests/testthat/test-fc-provenance.R")
#
# Test map:
#   P1 compute_group_mean_columns: math, group resolution, absent group -> NA
#   P2 build_final_results_generic: block order and naming of the new columns
#   P3 log2FC <-> linearFC is an exact round trip (the signed-reciprocal rule)
#   P4 build_rnaseq_summary_df emits log2FC.<contrast> from log2FoldChange
#   P5 metabolomics/lipidomics are untouched (no log2FC column, no norm block)
#   P6 build_provenance_notes documents every column family a mode emits

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

source(find_repo_file("R/core/01_io.R"))                       # normalize_contrast_name
source(find_repo_file("R/core/05_export_excel.R"))
source(find_repo_file("R/domain/rnaseq/04_de_summary.R"))
source(find_repo_file("R/domain/proteomics/05_de_summary.R"))
source(find_repo_file("R/domain/proteomics/06_outputs_legacy.R"))


# ---- Shared fixtures --------------------------------------------------------

prov_meta <- function() {
    data.frame(
        SampleID  = c("S_1", "S_2", "NS_1", "NS_2"),
        condition = c("S", "S", "NS", "NS"),
        stringsAsFactors = FALSE
    )
}

prov_contrasts <- function() {
    data.frame(
        Contrast_name = "S_vs_NS",
        Factor        = "condition",
        Numerator     = "S",
        Denominator   = "NS",
        stringsAsFactors = FALSE
    )
}

# Raw counts with a deliberate depth imbalance: the S libraries are twice as
# deep, so the raw ratio and the normalized ratio disagree by exactly 1 in log2.
prov_raw <- function() {
    matrix(
        c(200, 200, 100, 100,
          400, 400, 800, 800),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("g1", "g2"), c("S_1", "S_2", "NS_1", "NS_2"))
    )
}

prov_norm <- function() {
    m <- prov_raw()
    sweep(m, 2, c(2, 2, 1, 1), "/")   # size factors: S deep, NS shallow
}


# =============================================================================
# P1 — compute_group_mean_columns
# =============================================================================
test_that("P1 group means are the arithmetic mean over that group's replicates", {
    means <- compute_group_mean_columns(
        expr          = prov_norm(),
        sample_meta   = prov_meta(),
        sample_id_col = "SampleID",
        contrasts_df  = prov_contrasts()
    )

    expect_equal(sort(colnames(means)), c("Mean.NS", "Mean.S"))
    expect_equal(rownames(means), c("g1", "g2"))
    # g1 normalized: S = 100,100 ; NS = 100,100
    expect_equal(means["g1", "Mean.S"], 100)
    expect_equal(means["g1", "Mean.NS"], 100)
    # g2 normalized: S = 200,200 ; NS = 800,800
    expect_equal(means["g2", "Mean.S"], 200)
    expect_equal(means["g2", "Mean.NS"], 800)
})

test_that("P1 a group with no samples in the matrix yields NA, not an error", {
    meta <- prov_meta()
    contrasts_df <- prov_contrasts()
    contrasts_df$Denominator <- "MISSING_GROUP"

    means <- compute_group_mean_columns(prov_norm(), meta, "SampleID", contrasts_df)
    expect_true(all(is.na(means[["Mean.MISSING_GROUP"]])))
    expect_false(any(is.na(means[["Mean.S"]])))
})

test_that("P1 means and CV resolve group membership identically", {
    # Both go through compute_group_stat_columns, so the two blocks can never
    # disagree about which samples make up a group.
    args <- list(prov_norm(), prov_meta(), "SampleID", prov_contrasts())
    means <- do.call(compute_group_mean_columns, args)
    cvs   <- do.call(compute_group_cv_columns, args)

    expect_equal(sub("^Mean\\.", "", colnames(means)),
                 sub("^CV\\.", "", colnames(cvs)))
})

test_that("P1 an ambiguous grouping column warns and returns NULL", {
    meta <- prov_meta()
    meta$genotype <- c("ko", "ko", "wt", "wt")
    contrasts_df <- rbind(
        prov_contrasts(),
        data.frame(Contrast_name = "ko_vs_wt", Factor = "genotype",
                   Numerator = "ko", Denominator = "wt", stringsAsFactors = FALSE)
    )

    expect_warning(
        res <- compute_group_mean_columns(prov_norm(), meta, "SampleID", contrasts_df),
        "single grouping column"
    )
    expect_null(res)
})


# =============================================================================
# P2 — block order and naming in build_final_results_generic
# =============================================================================
prov_summary_df <- function() {
    data.frame(
        Gene                = c("g1", "g2"),
        log2FC.S_vs_NS      = c(0, -2),
        linearFC.S_vs_NS    = c(1, -4),
        pvalue.S_vs_NS      = c(0.9, 0.001),
        padj.S_vs_NS        = c(0.95, 0.01),
        S_vs_NS_pass        = c(NA, 1),
        pass_any_contrast   = c(NA, 1),
        stringsAsFactors    = FALSE,
        check.names         = FALSE
    )
}

test_that("P2 final results carry raw, .norm, Mean and CV blocks in that order", {
    means <- compute_group_mean_columns(prov_norm(), prov_meta(), "SampleID", prov_contrasts())
    cvs   <- compute_group_cv_columns(prov_norm(), prov_meta(), "SampleID", prov_contrasts())

    fr <- build_final_results_generic(
        summary_df     = prov_summary_df(),
        expr_df        = prov_raw(),
        contrasts_df   = prov_contrasts(),
        feature_id_col = "Gene",
        mode           = "rna",
        cv_cols        = cvs,
        norm_expr      = prov_norm(),
        mean_cols      = means
    )

    nm <- names(fr)
    expect_equal(nm[1], "Gene")
    expect_equal(nm[2:5], c("S_1", "S_2", "NS_1", "NS_2"))
    expect_equal(nm[6:9], c("S_1.norm", "S_2.norm", "NS_1.norm", "NS_2.norm"))
    # Mean before CV, both before the DE statistics
    expect_true(all(which(grepl("^Mean\\.", nm)) < min(which(grepl("^CV\\.", nm)))))
    expect_true(max(which(grepl("^CV\\.", nm))) < which(nm == "log2FC.S_vs_NS"))
    # log2FC immediately precedes linearFC: the estimate, then its presentation
    expect_equal(which(nm == "linearFC.S_vs_NS"), which(nm == "log2FC.S_vs_NS") + 1L)
})

test_that("P2 the .norm block holds normalized values, not the raw ones", {
    fr <- build_final_results_generic(
        summary_df     = prov_summary_df(),
        expr_df        = prov_raw(),
        contrasts_df   = prov_contrasts(),
        feature_id_col = "Gene",
        mode           = "rna",
        norm_expr      = prov_norm()
    )

    g2 <- fr[fr$Gene == "g2", ]
    expect_equal(g2$S_1, 400)          # raw
    expect_equal(g2$S_1.norm, 200)     # raw / size factor 2
    expect_equal(g2$NS_1, 800)
    expect_equal(g2$NS_1.norm, 800)    # size factor 1
})

test_that("P2 the raw columns alone do NOT reproduce log2FC; the .norm ones do", {
    # This is the reported symptom, pinned as a test: with unequal library
    # sizes the raw ratio is off by exactly the depth ratio (1 in log2 here).
    raw <- prov_raw()
    nrm <- prov_norm()

    raw_lfc  <- log2(mean(raw["g2", c("S_1", "S_2")]) / mean(raw["g2", c("NS_1", "NS_2")]))
    norm_lfc <- log2(mean(nrm["g2", c("S_1", "S_2")]) / mean(nrm["g2", c("NS_1", "NS_2")]))

    expect_equal(raw_lfc, -1)
    expect_equal(norm_lfc, -2)
    expect_equal(norm_lfc, prov_summary_df()$log2FC.S_vs_NS[2])
    expect_false(isTRUE(all.equal(raw_lfc, norm_lfc)))
})

test_that("P2 omitting the new arguments reproduces the previous column set", {
    fr_old <- build_final_results_generic(
        summary_df     = prov_summary_df()[, setdiff(names(prov_summary_df()), "log2FC.S_vs_NS")],
        expr_df        = prov_raw(),
        contrasts_df   = prov_contrasts(),
        feature_id_col = "Gene",
        mode           = "rna"
    )
    expect_false(any(grepl("^log2FC\\.|\\.norm$|^Mean\\.", names(fr_old))))
})


# =============================================================================
# P3 — log2FC <-> linearFC round trip
# =============================================================================
test_that("P3 linearFC is the signed reciprocal of 2^log2FC, exactly", {
    lfc <- c(-2, -0.6833, -0.0001, 0, 0.0001, 0.5, 3)
    linear <- ifelse(lfc >= 0, 2^lfc, -1 * (2^-lfc))

    # Recovering log2FC from linearFC: the sign carries the direction
    back <- sign(linear) * log2(abs(linear))
    expect_equal(back, lfc)

    # A reader reading -1.61 as "log2FC = -1.61" is off by more than a factor 2
    expect_equal(round(linear[2], 2), -1.61)
    expect_false(isTRUE(all.equal(linear[2], lfc[2], tolerance = 1e-3)))
})


# =============================================================================
# P4 — build_rnaseq_summary_df emits log2FC
# =============================================================================
test_that("P4 rnaseq summary carries log2FC.<contrast> alongside linearFC", {
    de_tables <- list(S_vs_NS = data.frame(
        FeatureID      = c("g1", "g2"),
        log2FoldChange = c(-0.6833, 2),
        pvalue         = c(0.028, 1e-6),
        padj           = c(0.54, 1e-4),
        stringsAsFactors = FALSE
    ))

    sdf <- build_rnaseq_summary_df(de_tables, list(p_cutoff = 0.05, linear_fc_cutoff = 1.5))

    expect_true("log2FC.S_vs_NS" %in% names(sdf))
    # Unrounded: the column must be the model estimate itself, so that the
    # documented rule reproduces linearFC exactly rather than nearly.
    expect_equal(sdf$log2FC.S_vs_NS, c(-0.6833, 2))
    lfc <- sdf$log2FC.S_vs_NS
    expect_equal(sdf$linearFC.S_vs_NS,
                 signif(ifelse(lfc >= 0, 2^lfc, -1 * (2^-lfc)), 3))
})

test_that("P4 the log2FC -> linearFC round trip holds for every feature", {
    # Regression guard: rounding log2FC before storing it made signif(2^x, 3)
    # disagree with the reported linearFC for ~1% of features.
    lfc <- seq(-5, 5, by = 0.0007)
    de_tables <- list(A_vs_B = data.frame(
        FeatureID      = paste0("g", seq_along(lfc)),
        log2FoldChange = lfc,
        pvalue         = 0.5,
        padj           = 0.5,
        stringsAsFactors = FALSE
    ))

    sdf <- build_rnaseq_summary_df(de_tables, list(p_cutoff = 0.05, linear_fc_cutoff = 1.5))
    l <- sdf$log2FC.A_vs_B
    expect_equal(sdf$linearFC.A_vs_B, signif(ifelse(l >= 0, 2^l, -1 * (2^-l)), 3))
})

test_that("P4 get_contrast_cols exposes log2fc for rna and proteomics only", {
    expect_equal(get_contrast_cols("A_vs_B", mode = "rna")$log2fc, "log2FC.A_vs_B")
    expect_equal(get_contrast_cols("A_vs_B", mode = "proteomics")$log2fc,
                 "log2FC.imputs.A_vs_B")
    expect_null(get_contrast_cols("A_vs_B", mode = "metabolomics")$log2fc)
    expect_null(get_contrast_cols("A_vs_B", mode = "lipidomics")$log2fc)
})


# =============================================================================
# P5 — metabolomics / lipidomics are untouched
# =============================================================================
test_that("P5 a metabolomics summary_df with a log2FC column still emits none", {
    sdf <- prov_summary_df()
    names(sdf)[names(sdf) == "S_vs_NS_pass"] <- "pass.S_vs_NS"

    fr <- build_final_results_generic(
        summary_df     = sdf,
        expr_df        = prov_raw(),
        contrasts_df   = prov_contrasts(),
        feature_id_col = "Gene",
        mode           = "metabolomics"
    )

    # log2FC.S_vs_NS exists in the input but metabolomics does not declare it,
    # so it must not leak into the exported table.
    expect_false("log2FC.S_vs_NS" %in% names(fr))
    expect_true("linearFC.S_vs_NS" %in% names(fr))
})


# =============================================================================
# P7 — proteomics: log2FC.imputs and the imputed .norm block
# =============================================================================
prot_runs <- function(n_runs = 3) {
    lapply(seq_len(n_runs), function(i) {
        list(S_vs_NS = data.frame(
            FeatureID     = c("p1", "p2"),
            Protein.Names = c("PROT1_TEST", "PROT2_TEST"),
            Genes         = c("gene1", "gene2"),
            First.Protein.Description = c("test protein 1", "test protein 2"),
            # Slightly different per run, so the consensus is a genuine average
            logFC      = c(-0.7 + 0.01 * i, 2 - 0.01 * i),
            P.Value    = c(0.02, 1e-5),
            adj.P.Val  = c(0.3, 1e-3),
            stringsAsFactors = FALSE
        ))
    })
}

prot_config <- function() {
    list(modes = list(proteomics = list(
        de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5, use_adj_for_pass1 = TRUE),
        imputation = list(multi_imputation = TRUE, no_repetitions = 3, min_no_passed = 2),
        de_table = list(id_col = "FeatureID"),
        effects = list(samples = "SampleID")
    )))
}

test_that("P7 proteomics summary emits log2FC.imputs consistent with linearFC.imputs", {
    sdf <- summarize_limma_mult_imputation(prot_runs(), prot_config())

    expect_true("log2FC.imputs.S_vs_NS" %in% names(sdf))
    # It is the log2 of the consensus linear ratio, not the mean of per-run logFCs
    expect_equal(sdf$log2FC.imputs.S_vs_NS, log2(sdf$linearRatio.imputs.S_vs_NS))

    l <- sdf$log2FC.imputs.S_vs_NS
    ratio <- 2^l
    expect_equal(sdf$linearFC.imputs.S_vs_NS,
                 signif(ifelse(ratio >= 1, ratio, -1 / ratio), 3))
})

test_that("P7 proteomics final results carry the imputed block and group means", {
    meta <- prov_meta()
    observed <- matrix(c(10, 10, 12, 12,
                         8, NA, 9, 9),
                       nrow = 2, byrow = TRUE,
                       dimnames = list(c("p1", "p2"), meta$SampleID))
    imputed <- observed
    imputed["p2", "S_2"] <- 7.5   # the value limma actually saw

    pre <- list(expr_filt = observed, expr_imp_single = imputed, meta = meta,
                row_data = NULL)
    sdf <- summarize_limma_mult_imputation(prot_runs(), prot_config())
    sdf <- sdf[match(c("p1", "p2"), sdf$FeatureID), , drop = FALSE]

    fr <- build_final_results_proteomics(
        pre = pre, summary_df = sdf, contrasts_df = prov_contrasts(),
        feature_id_col = "FeatureID", config = prot_config()
    )

    expect_true(all(paste0(meta$SampleID, ".norm") %in% names(fr)))
    expect_true(all(c("Mean.S", "Mean.NS") %in% names(fr)))
    # The unimputed column keeps the gap; the .norm column carries the fill-in
    p2 <- fr[fr$FeatureID == "p2", ]
    expect_true(is.na(p2$S_2))
    expect_equal(p2$S_2.norm, 7.5)
    # Means come from the imputed matrix, so the gap does not skew the group
    expect_equal(p2$Mean.S, mean(c(8, 7.5)))
    expect_equal(which(names(fr) == "linearFC.imputs.S_vs_NS"),
                 which(names(fr) == "log2FC.imputs.S_vs_NS") + 1L)
})


# =============================================================================
# P6 — the "How to read" sheet documents what each mode actually emits
# =============================================================================
test_that("P6 provenance notes cover every column family the mode emits", {
    rna <- build_provenance_notes("rna")
    expect_true(all(c("<sample>", "<sample>.norm", "Mean.<group>", "CV.<group>",
                      "log2FC.<contrast>", "linearFC.<contrast>") %in% rna$glossary$Column))

    prot <- build_provenance_notes("proteomics")
    expect_true(all(c("<sample>", "<sample>.norm", "Mean.<group>", "CV.<group>",
                      "log2FC.imputs.<contrast>", "linearFC.imputs.<contrast>") %in%
                    prot$glossary$Column))

    # Both must state the signed-reciprocal rule and flag the approximation
    expect_true(any(grepl("signed linear fold change", rna$glossary$Meaning)))
    expect_true(any(grepl("approximation", rna$notes)))
    expect_true(any(grepl("approximation", prot$notes)))

    # Unknown modes still get a usable glossary rather than an error
    expect_s3_class(build_provenance_notes("something_else")$glossary, "data.frame")
})
