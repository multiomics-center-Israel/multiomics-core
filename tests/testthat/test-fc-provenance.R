# Tests for the fold-change provenance columns in Final_results_{ALL,DE}.
#
# Context: the per-sample columns are raw counts (RNA-seq) or unimputed log2
# intensities (proteomics), while the reported fold change is estimated on the
# normalized/imputed matrix. Readers who eyeball the sample columns therefore
# compute a different number. These tests pin the columns that close that gap:
# log2FC next to linearFC, a <sample>.norm block, and Mean.<group>.
#
# Base R and openxlsx-free throughout, except P9's last case, which needs
# DESeq2 and skips when it is unavailable — so the file still runs from
# RStudio on Windows before pushing:
#
#   testthat::test_file("tests/testthat/test-fc-provenance.R")
#
# Test map:
#   P1 compute_group_mean_columns: math, group resolution, absent group -> NA
#   P2 build_final_results_generic: block order and naming of the new columns
#   P3 log2FC <-> linearFC is an exact round trip (the signed-reciprocal rule)
#   P4 build_rnaseq_summary_df emits log2FC.<contrast> from log2FoldChange
#   P5 metabolomics/lipidomics are untouched (no log2FC column, no norm block)
#   P8 linearFC vs log2FC, each derived independently of the production code
#   P9 log2FC_from_means: model-free estimate, and DESeq2 shrinkage detection
#   P7 proteomics: log2FC.imputs and the imputed .norm block
#   P10 the analyst-facing shrinkage alert (heavy shrinkage, x = 0 stripe)
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

source(find_repo_file("R/core/01_io.R"))                        # normalize_contrast_name
source(find_repo_file("R/core/05_export_excel.R"))
source(find_repo_file("R/core/de_shrinkage_check.R"))
source(find_repo_file("R/domain/rnaseq/00a_deseq_factory.R"))   # create_deseq_dataset
source(find_repo_file("R/domain/rnaseq/01_expression.R"))       # compute_cpm
source(find_repo_file("R/domain/rnaseq/04_de_summary.R"))
source(find_repo_file("R/domain/rnaseq/05_outputs_legacy.R"))   # build_final_results_rnaseq
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
    # Stat block order: the model estimate, its model-free counterpart, then
    # the linear presentation of the model estimate.
    nm <- names(fr)
    expect_equal(nm[which(nm == "log2FC.imputs.S_vs_NS") + 1L],
                 "log2FC_from_means.S_vs_NS")
    expect_equal(nm[which(nm == "log2FC.imputs.S_vs_NS") + 2L],
                 "linearFC.imputs.S_vs_NS")
    # And it is recomputable from the two Mean cells in the same row
    expect_equal(fr$log2FC_from_means.S_vs_NS, fr$Mean.S - fr$Mean.NS)
})


# =============================================================================
# P8 — linearFC vs log2FC, each derived independently of the production code
#      These deliberately do NOT reuse the pipeline's own expression. The
#      expected values below are hand-computed constants, and the round trip is
#      written a second, different way (exp/log rather than 2^ and ifelse), so a
#      change to the production formula cannot make these tests agree with it by
#      construction.
# =============================================================================
test_that("P8 linearFC matches hand-computed constants for known log2FC values", {
    # Hand-computed, not derived from any formula in R/:
    #   log2FC  0      -> ratio 1        -> linearFC  1
    #   log2FC  1      -> ratio 2        -> linearFC  2
    #   log2FC  2      -> ratio 4        -> linearFC  4
    #   log2FC -1      -> ratio 0.5      -> 1/0.5 = 2   -> linearFC -2
    #   log2FC -2      -> ratio 0.25     -> 1/0.25 = 4  -> linearFC -4
    #   log2FC  0.585  -> ratio 1.5      -> linearFC  1.5   (log2(1.5)=0.5849625)
    #   log2FC -0.585  -> ratio 0.6667   -> linearFC -1.5
    #   log2FC -0.6833 -> ratio 0.62114  -> 1/0.62114 = 1.61 -> linearFC -1.61
    cases <- data.frame(
        log2FC   = c(0, 1, 2, -1, -2, 0.5849625, -0.5849625, -0.6833),
        expected = c(1, 2, 4, -2, -4, 1.5,       -1.5,       -1.61),
        stringsAsFactors = FALSE
    )

    de_tables <- list(A_vs_B = data.frame(
        FeatureID      = paste0("g", seq_len(nrow(cases))),
        log2FoldChange = cases$log2FC,
        pvalue         = 0.5,
        padj           = 0.5,
        stringsAsFactors = FALSE
    ))
    sdf <- build_rnaseq_summary_df(de_tables, list(p_cutoff = 0.05, linear_fc_cutoff = 1.5))

    expect_equal(sdf$log2FC.A_vs_B, cases$log2FC)
    expect_equal(sdf$linearFC.A_vs_B, cases$expected, tolerance = 1e-3)
})

test_that("P8 the two columns agree when each is derived on its own", {
    lfc <- seq(-6, 6, by = 0.013)
    de_tables <- list(A_vs_B = data.frame(
        FeatureID      = paste0("g", seq_along(lfc)),
        log2FoldChange = lfc,
        pvalue         = 0.5,
        padj           = 0.5,
        stringsAsFactors = FALSE
    ))
    sdf <- build_rnaseq_summary_df(de_tables, list(p_cutoff = 0.05, linear_fc_cutoff = 1.5))

    written_log2 <- sdf$log2FC.A_vs_B
    written_lin  <- sdf$linearFC.A_vs_B

    # Derivation A: log2FC -> linearFC, via exp/log rather than 2^ and ifelse
    ratio <- exp(written_log2 * log(2))
    derived_lin <- ratio
    below <- ratio < 1
    derived_lin[below] <- -1 / ratio[below]
    expect_equal(written_lin, signif(derived_lin, 3))

    # Derivation B: linearFC -> log2FC, independently of derivation A
    magnitude <- abs(written_lin)
    magnitude[magnitude < 1] <- 1 / magnitude[magnitude < 1]   # no-op guard
    derived_log2 <- log(magnitude, base = 2) * sign(written_lin)
    # signif(, 3) on linearFC caps the recoverable precision
    expect_equal(derived_log2, written_log2, tolerance = 1e-2)

    # And the signs must never disagree between the two columns
    expect_true(all(sign(written_lin) == sign(written_log2) |
                    written_log2 == 0))
})

test_that("P8 proteomics linearFC.imputs agrees with its own log2FC.imputs", {
    sdf <- summarize_limma_mult_imputation(prot_runs(), prot_config())

    written_log2 <- sdf$log2FC.imputs.S_vs_NS
    written_lin  <- sdf$linearFC.imputs.S_vs_NS

    ratio <- exp(written_log2 * log(2))
    derived_lin <- ratio
    below <- ratio < 1
    derived_lin[below] <- -1 / ratio[below]
    expect_equal(written_lin, signif(derived_lin, 3))
})


# =============================================================================
# P9 — log2FC_from_means: the model-free estimate, and shrinkage detection
# =============================================================================
test_that("P9 the naive log2FC is exactly the log ratio of the Mean columns", {
    means <- compute_group_mean_columns(prov_norm(), prov_meta(), "SampleID", prov_contrasts())
    naive <- compute_naive_log2fc_columns(means, prov_contrasts(), scale = "linear")

    expect_equal(colnames(naive), "S_vs_NS")
    # g1 normalized: S mean 100, NS mean 100 -> 0.  g2: 200 vs 800 -> -2.
    expect_equal(naive[["S_vs_NS"]], c(0, -2))
    # Recomputed by hand from the exported cells, which is the point of it
    expect_equal(naive[["S_vs_NS"]], log2(means$Mean.S / means$Mean.NS))
})

test_that("P9 on a log2 scale the naive log2FC is a difference, not a ratio", {
    log_means <- data.frame(Mean.S = c(10, 3), Mean.NS = c(8, 5),
                            row.names = c("p1", "p2"), check.names = FALSE)
    naive <- compute_naive_log2fc_columns(log_means, prov_contrasts(), scale = "log2")
    expect_equal(naive[["S_vs_NS"]], c(2, -2))
})

test_that("P9 a non-positive group mean gives NA, not -Inf", {
    lin_means <- data.frame(Mean.S = c(100, 0), Mean.NS = c(50, 50),
                            row.names = c("g1", "g2"), check.names = FALSE)
    naive <- compute_naive_log2fc_columns(lin_means, prov_contrasts(), scale = "linear")
    expect_equal(naive[["S_vs_NS"]][1], 1)
    expect_true(is.na(naive[["S_vs_NS"]][2]))
})

test_that("P9 a contrast with no matching Mean columns warns and is skipped", {
    means <- compute_group_mean_columns(prov_norm(), prov_meta(), "SampleID", prov_contrasts())
    other <- data.frame(Contrast_name = "X_vs_Y", Factor = "condition",
                        Numerator = "X", Denominator = "Y", stringsAsFactors = FALSE)
    expect_warning(res <- compute_naive_log2fc_columns(means, other, scale = "linear"),
                   "no group means for contrast")
    expect_null(res)
})

test_that("P9 the naive column lands next to log2FC in the exported table", {
    means <- compute_group_mean_columns(prov_norm(), prov_meta(), "SampleID", prov_contrasts())
    naive <- compute_naive_log2fc_columns(means, prov_contrasts(), scale = "linear")

    fr <- build_final_results_generic(
        summary_df     = prov_summary_df(),
        expr_df        = prov_raw(),
        contrasts_df   = prov_contrasts(),
        feature_id_col = "Gene",
        mode           = "rna",
        norm_expr      = prov_norm(),
        mean_cols      = means,
        naive_log2fc   = naive
    )

    nm <- names(fr)
    expect_equal(nm[which(nm == "log2FC.S_vs_NS") + 1L], "log2FC_from_means.S_vs_NS")
    expect_equal(nm[which(nm == "log2FC.S_vs_NS") + 2L], "linearFC.S_vs_NS")
    # Recomputable from the neighbouring Mean cells
    expect_equal(fr$log2FC_from_means.S_vs_NS, log2(fr$Mean.S / fr$Mean.NS))
})

test_that("P9 DESeq2 shrinkage is visible as a gap against log2FC_from_means", {
    skip_if_not_installed("DESeq2")

    meta <- data.frame(
        SampleID  = c("A1", "A2", "A3", "B1", "B2", "B3"),
        condition = c("A", "A", "A", "B", "B", "B"),
        stringsAsFactors = FALSE
    )
    contrasts_df <- data.frame(
        Contrast_name = "A_vs_B", Factor = "condition",
        Numerator = "A", Denominator = "B", stringsAsFactors = FALSE
    )

    # Deterministic fixture: half the genes well expressed, half at low counts.
    # Within each half, a third are unchanged and the rest move symmetrically up
    # and down, so median-of-ratios normalization has a stable baseline to sit
    # on (if every gene moved the same way, the size factors would absorb the
    # effect and nothing would be left to shrink).
    n_gene   <- 600L
    base     <- c(rep(2000, n_gene / 2), rep(8, n_gene / 2))
    true_lfc <- rep(rep(c(0, 1.5, -1.5), each = 100L), times = 2L)

    counts <- withr::with_seed(42, {
        m <- vapply(seq_len(6), function(j) {
            direction <- if (j <= 3) 0.5 else -0.5
            stats::rnbinom(n_gene, mu = base * 2^(true_lfc * direction), size = 20)
        }, numeric(n_gene))
        dimnames(m) <- list(paste0("g", seq_len(n_gene)), meta$SampleID)
        m
    })

    fit <- function(deseq_mode) {
        de_cfg <- list(method = "deseq2", deseq_mode = deseq_mode, p_cutoff = 0.05,
                       linear_fc_cutoff = 1.5, sample_col = "SampleID")
        de_res <- run_deseq2_de(counts, meta, contrasts_df, de_cfg)
        sdf <- build_rnaseq_summary_df(de_res$tables, de_cfg)
        names(sdf)[names(sdf) == "FeatureID"] <- "Gene"

        norm_counts <- deseq2_normalized_counts(de_res$dds)
        means <- compute_group_mean_columns(norm_counts, meta, "SampleID", contrasts_df)
        naive <- compute_naive_log2fc_columns(means, contrasts_df, scale = "linear")

        pre <- list(expr_filt = counts, meta = meta, row_data = NULL,
                    info = list(source_type = "matrix"))
        config <- list(modes = list(rna = list(
            de = de_cfg, excel = list(group_cv = TRUE)
        )))
        build_final_results_rnaseq(pre, sdf, contrasts_df, config = config,
                                   norm_counts = norm_counts)
    }

    fr_default <- suppressWarnings(fit("default"))
    fr_legacy  <- suppressWarnings(fit("legacy"))

    expect_true("log2FC_from_means.A_vs_B" %in% names(fr_default))

    # Only the genes that were actually made to move; unchanged genes have a
    # naive log2FC near zero, so the ratio below would be pure noise for them.
    moved <- true_lfc != 0
    high  <- which(moved & seq_len(n_gene) <= n_gene / 2)
    low   <- which(moved & seq_len(n_gene) >  n_gene / 2)

    shrink_ratio <- function(fr, idx) {
        median(abs(fr$log2FC.A_vs_B[idx]) /
               abs(fr$log2FC_from_means.A_vs_B[idx]), na.rm = TRUE)
    }

    # The naive column must recover the effect that was put in, on both halves
    expect_equal(median(abs(fr_default$log2FC_from_means.A_vs_B[high])), 1.5,
                 tolerance = 0.1)
    expect_equal(median(abs(fr_default$log2FC_from_means.A_vs_B[low])), 1.5,
                 tolerance = 0.2)

    # Default mode does not shrink: the model estimate tracks the naive one
    expect_gt(shrink_ratio(fr_default, high), 0.95)
    expect_gt(shrink_ratio(fr_default, low),  0.95)

    # betaPrior pulls estimates towards zero. The naive column is model-free and
    # does not move, so the ratio drops — and it drops further where the counts
    # are low, which is where the prior has the most say.
    expect_lt(shrink_ratio(fr_legacy, low), 0.90)
    expect_lt(shrink_ratio(fr_legacy, low), shrink_ratio(fr_legacy, high))
    expect_gt(shrink_ratio(fr_legacy, high), 0.90)
})


# =============================================================================
# P10 — the analyst-facing shrinkage alert
# =============================================================================
shrink_table <- function(model_lfc, naive_lfc, padj) {
    data.frame(
        Gene                       = paste0("g", seq_along(model_lfc)),
        log2FC.S_vs_NS             = model_lfc,
        log2FC_from_means.S_vs_NS  = naive_lfc,
        linearFC.S_vs_NS           = signif(ifelse(model_lfc >= 0, 2^model_lfc,
                                                   -1 * (2^-model_lfc)), 3),
        pvalue.S_vs_NS             = padj,
        padj.S_vs_NS               = padj,
        stringsAsFactors           = FALSE,
        check.names                = FALSE
    )
}

# 100 features with a real effect, 100 null ones. The null half exists to prove
# it cannot drag the verdict either way (its ratio is 0/0-ish by construction).
shrink_fixture <- function(factor_applied = 1) {
    naive <- c(rep(c(2, -2), each = 50L), rep(c(0.002, -0.002), each = 50L))
    model <- naive * factor_applied
    padj  <- c(rep(1e-6, 100L), rep(0.9, 100L))
    shrink_table(model, naive, padj)
}

test_that("P10 a healthy fit is not flagged", {
    chk <- check_log2fc_shrinkage(shrink_fixture(1), "S_vs_NS", mode = "rna")
    expect_equal(chk$flag, "ok")
    expect_equal(chk$median_ratio, 1)
    expect_equal(chk$n_considered, 100)   # the null half is excluded by the floor
    expect_equal(chk$frac_flat, 0)
})

test_that("P10 ordinary betaPrior-scale shrinkage is not flagged", {
    # ~0.85 is what DESeq2 betaPrior produced in the P9 fixture; that is normal
    # behaviour and must not fire the alert.
    chk <- check_log2fc_shrinkage(shrink_fixture(0.85), "S_vs_NS", mode = "rna")
    expect_equal(chk$flag, "ok")
    expect_equal(chk$median_ratio, 0.85)
})

test_that("P10 heavy shrinkage is flagged as 'shrunk'", {
    chk <- check_log2fc_shrinkage(shrink_fixture(0.3), "S_vs_NS", mode = "rna")
    expect_equal(chk$flag, "shrunk")
    expect_equal(chk$median_ratio, 0.3)
    expect_equal(chk$frac_flat, 0)        # shrunk, but not flattened to zero
})

test_that("P10 a collapse to zero is flagged as 'collapsed'", {
    chk <- check_log2fc_shrinkage(shrink_fixture(1e-5), "S_vs_NS", mode = "rna")
    expect_equal(chk$flag, "collapsed")
    expect_equal(chk$frac_flat, 1)
    # Every significant feature sits at log2FC ~ 0: the x = 0 volcano stripe
    expect_equal(chk$frac_sig_flat, 1)
    expect_equal(chk$n_significant, 100)
})

test_that("P10 the x = 0 stripe is detected even when the median ratio looks fine", {
    # Most features are estimated correctly, so the median ratio is 1, but a
    # block of significant features has been flattened. The median would miss
    # this; the stripe check must not.
    naive <- rep(c(2, -2), each = 100L)
    model <- naive
    model[1:30] <- 0                       # flattened, and significant below
    padj  <- c(rep(1e-6, 100L), rep(1e-6, 100L))

    chk <- check_log2fc_shrinkage(shrink_table(model, naive, padj), "S_vs_NS", mode = "rna")
    expect_equal(chk$median_ratio, 1)      # the median is blind to it
    expect_equal(chk$frac_sig_flat, 0.15)  # 30 of 200 significant features
    expect_equal(chk$flag, "collapsed")
})

test_that("P10 null features cannot trigger the alert on their own", {
    # This is the evm.TU.ptg000675l_np1212.2 situation: both estimates ~0, so
    # their ratio is unstable and meaningless. The floor must exclude them.
    naive <- rep(c(0.005, -0.005), each = 100L)
    model <- naive / 4.6                   # a ratio of ~0.22, from noise alone
    chk <- check_log2fc_shrinkage(shrink_table(model, naive, rep(0.99, 200L)),
                                  "S_vs_NS", mode = "rna")
    expect_equal(chk$n_considered, 0)
    expect_true(is.na(chk$median_ratio))
    expect_equal(chk$flag, "ok")
})

test_that("P10 the check is skipped, not errored, when the columns are absent", {
    df <- prov_summary_df()                # has log2FC but no log2FC_from_means
    expect_null(check_log2fc_shrinkage(df, "S_vs_NS", mode = "rna"))
    expect_null(check_log2fc_shrinkage(df, "S_vs_NS", mode = "metabolomics"))
})

test_that("P10 the alert names the contrast, the symptom and the volcano stripe", {
    chk <- check_log2fc_shrinkage(shrink_fixture(1e-5), "S_vs_NS", mode = "rna")
    expect_warning(warn_log2fc_shrinkage(chk, mode = "rna"), "S_vs_NS")
    expect_warning(warn_log2fc_shrinkage(chk, mode = "rna"), "COLLAPSED")
    expect_warning(warn_log2fc_shrinkage(chk, mode = "rna"), "x = 0")
    expect_warning(warn_log2fc_shrinkage(chk, mode = "rna"), "log2fc_shrinkage_check.tsv")

    # A healthy run reports, but does not warn
    ok <- check_log2fc_shrinkage(shrink_fixture(1), "S_vs_NS", mode = "rna")
    expect_silent(suppressMessages(warn_log2fc_shrinkage(ok, mode = "rna")))
})

test_that("P10 proteomics column naming resolves too", {
    df <- data.frame(
        FeatureID                      = paste0("p", 1:4),
        `log2FC.imputs.S_vs_NS`        = c(2, -2, 2, -2),
        `log2FC_from_means.S_vs_NS`    = c(2, -2, 2, -2),
        `padj.imputs.S_vs_NS`          = 1e-6,
        stringsAsFactors = FALSE, check.names = FALSE
    )
    chk <- check_log2fc_shrinkage(df, "S_vs_NS", mode = "proteomics")
    expect_equal(chk$median_ratio, 1)
    expect_equal(chk$flag, "ok")
})


# =============================================================================
# P6 — the "How to read" sheet documents what each mode actually emits
# =============================================================================
test_that("P6 provenance notes cover every column family the mode emits", {
    rna <- build_provenance_notes("rna")
    expect_true(all(c("<sample>", "<sample>.norm", "Mean.<group>", "CV.<group>",
                      "log2FC.<contrast>", "log2FC_from_means.<contrast>",
                      "linearFC.<contrast>") %in% rna$glossary$Column))

    prot <- build_provenance_notes("proteomics")
    expect_true(all(c("<sample>", "<sample>.norm", "Mean.<group>", "CV.<group>",
                      "log2FC.imputs.<contrast>", "log2FC_from_means.<contrast>",
                      "linearFC.imputs.<contrast>") %in% prot$glossary$Column))

    # Both must name the shrinkage comparison the naive column exists for
    expect_true(any(grepl("shrinkage", rna$notes)))
    expect_true(any(grepl("shrinkage", prot$notes)))

    # Both must state the signed-reciprocal rule and say the two log2FC columns
    # are not expected to agree to the last digit
    expect_true(any(grepl("signed linear fold change", rna$glossary$Meaning)))
    expect_true(any(grepl("agree closely|can differ", rna$notes)))
    expect_true(any(grepl("differ for two reasons", prot$notes)))

    # Unknown modes still get a usable glossary rather than an error
    expect_s3_class(build_provenance_notes("something_else")$glossary, "data.frame")
})
