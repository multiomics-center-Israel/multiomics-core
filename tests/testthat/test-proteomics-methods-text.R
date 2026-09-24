# tests/testthat/test-proteomics-methods-text.R
#
# The reader-facing Methods text for the proteomics report.
#
# The section used to be prose inline in the template, and it drifted: it
# described KNN imputation the pipeline never performs, named limma-voom and
# DESeq2 (RNA-seq methods that cannot run here), stated Benjamini-Hochberg while
# de$p_adjust_method is configurable, and defaulted imputation$method to "none"
# while the dispatcher defaults it to perseus_like. These tests exist so a claim
# can only enter the Methods section together with the condition that makes it
# true.
#
# build_proteomics_methods_text() lives in R/domain/proteomics/08_report.R and
# is loaded by helper.R along with the rest of R/.

cfg_for <- function(...) {
    prot <- utils::modifyList(
        list(
            engine = "DIANN",
            scale_in = "linear",
            filtering = list(remove_contaminants = TRUE, contaminant_prefix = "cRAP-",
                             min_count = 3, min_groups = 1),
            normalization = list(method = "none"),
            batch_correction = list(method = "none"),
            imputation = list(method = "perseus_like", width = 0.3, downshift = 1.8,
                              multi_imputation = FALSE, no_repetitions = 1,
                              min_no_passed = 1),
            de = list(method = "limma", p_cutoff = 0.05, linear_fc_cutoff = 1.5),
            clustering = list(enabled = FALSE)
        ),
        list(...)
    )
    list(modes = list(proteomics = prot))
}

blocks_for <- function(..., de_method = "limma",
                       pca_cc_file = NULL, hier_file = NULL) {
    build_proteomics_methods_text(cfg_for(...), de_method = de_method,
                                  pca_cc_file = pca_cc_file, hier_file = hier_file)
}

# A path that exists, for the two Methods sentences that describe a figure and
# are therefore gated on the figure being there rather than on the config flag
# that asked for it. The caller supplies the directory -- withr::local_tempdir()
# in the test body, so it is cleaned up with that test.
existing_file <- function(dir, name) {
    p <- file.path(dir, name)
    file.create(p)
    p
}


# =============================================================================
# Data processing
# =============================================================================

test_that("no-normalization is stated explicitly, not omitted", {
    # Silence in a Methods section reads as an oversight. normalization$method
    # defaults to "none" in the proteomics template, so this is the common case.
    txt <- blocks_for()$data_processing
    expect_match(txt, "No between-sample normalization was applied", fixed = TRUE)
    expect_false(grepl("median-centred", txt, fixed = TRUE))

    txt2 <- blocks_for(normalization = list(method = "median"))$data_processing
    expect_match(txt2, "median-centred on the log2 scale", fixed = TRUE)
    expect_false(grepl("No between-sample normalization", txt2, fixed = TRUE))
})

test_that("contaminant removal is claimed only when the flag is on", {
    # filter_contaminants() honours remove_contaminants and skips the step when
    # it is FALSE, so an unconditional sentence would describe work not done.
    on <- blocks_for()$data_processing
    expect_match(on, "were removed as contaminants", fixed = TRUE)

    off <- blocks_for(filtering = list(remove_contaminants = FALSE, min_count = 3,
                                       min_groups = 1))$data_processing
    expect_false(grepl("contaminant", off, fixed = TRUE))
})

test_that("the min-count phrase survives per-group thresholds", {
    # extract_min_count() accepts a bare number, a default plus overrides, or
    # per-group values alone. Assuming a scalar would misreport the last two.
    expect_match(methods_min_count_phrase(3), "at least 3 sample(s)", fixed = TRUE)
    expect_match(methods_min_count_phrase(NULL), "at least 1 sample", fixed = TRUE)

    both <- methods_min_count_phrase(list(default = 2, Treated = 4))
    expect_match(both, "at least 2 sample(s)", fixed = TRUE)
    expect_match(both, "Treated: 4", fixed = TRUE)

    only <- methods_min_count_phrase(list(Control = 3, Treated = 4))
    expect_match(only, "per-group thresholds", fixed = TRUE)
    expect_match(only, "Control: 3", fixed = TRUE)
})

test_that("batch correction describes the real matrix contract", {
    # 04b_batch_correction.R estimates and applies the correction on the
    # single-imputed complete matrix, restores the missingness mask to
    # expr_filt (which DE then re-imputes), and leaves QC/clustering on the
    # corrected complete matrix. "Batch correction after imputation" alone
    # would collapse all of that into something inaccurate.
    txt <- blocks_for(batch_correction = list(method = "combat"))$data_processing
    expect_match(txt, "ComBat", fixed = TRUE)
    expect_match(txt, "single-imputed complete matrix", fixed = TRUE)
    expect_match(txt, "original missingness pattern was then restored", fixed = TRUE)
    expect_match(txt, "re-imputed for differential analysis", fixed = TRUE)
    expect_match(txt, "quality-control and clustering use the corrected complete matrix",
                 fixed = TRUE)

    expect_false(grepl("Batch effects", blocks_for()$data_processing, fixed = TRUE))
})

test_that("the log2 sentence is tied to scale_in and stays loader-agnostic", {
    # The loaders differ -- DIA-NN linear uses log2(x), the preprocessed-linear
    # path log2(x + 1) -- so the wording must not commit to either form.
    txt <- blocks_for()$data_processing
    expect_match(txt, "transformed to the log2 scale", fixed = TRUE)
    expect_false(grepl("x + 1", txt, fixed = TRUE))
    expect_false(grepl("acquired using", txt, fixed = TRUE))
    expect_match(txt, "quantified using DIANN", fixed = TRUE)

    expect_false(grepl("transformed to the log2 scale",
                       blocks_for(scale_in = "log")$data_processing, fixed = TRUE))
})

test_that("an absent scale_in falls back to files$is_logtransformed", {
    # get_proteomics_expression_matrix() resolves a missing scale_in from
    # files$is_logtransformed and leaves such a table untransformed. Resolving it
    # to "linear" here would claim a log2 transformation these legacy configs
    # never received.
    legacy <- blocks_for(scale_in = NULL,
                         files = list(is_logtransformed = TRUE))$data_processing
    expect_false(grepl("transformed to the log2 scale", legacy, fixed = TRUE))

    # The same fallback in the other direction: no flag, no scale_in -> linear.
    expect_match(blocks_for(scale_in = NULL)$data_processing,
                 "transformed to the log2 scale", fixed = TRUE)
    expect_match(blocks_for(scale_in = NULL,
                            files = list(is_logtransformed = FALSE))$data_processing,
                 "transformed to the log2 scale", fixed = TRUE)
})

test_that("batch correction needs the resolved enabled state, not just a method", {
    # get_proteomics_batch_config() resolves enabled as bc$enabled, defaulting to
    # method != "none". An explicit enabled: false leaves correct_batch_proteomics()
    # returning uncorrected values, whatever the method says.
    on <- blocks_for(batch_correction = list(method = "combat"))$data_processing
    expect_match(on, "Batch effects were corrected with ComBat", fixed = TRUE)

    expect_false(grepl("Batch effects", blocks_for(batch_correction = list(
        method = "combat", enabled = FALSE))$data_processing, fixed = TRUE))
    expect_match(blocks_for(batch_correction = list(
        method = "probatch", enabled = TRUE))$data_processing,
        "corrected with proBatch", fixed = TRUE)
    expect_false(grepl("Batch effects", blocks_for(batch_correction = list(
        method = "none", enabled = TRUE))$data_processing, fixed = TRUE))
})


# =============================================================================
# Missing values
# =============================================================================

test_that("an absent imputation method resolves the way the dispatcher does", {
    # impute_proteomics() defaults to perseus_like. A Methods generator
    # defaulting to "none" would print nothing about imputation for a run that
    # imputed. This does not fix the other readers that still default "none" --
    # that divergence is tracked separately -- it only refuses to add another.
    expect_equal(methods_imputation_method(list()), "perseus_like")
    expect_equal(methods_imputation_method(list(method = "perseus")), "perseus_like")
    expect_false(is.null(blocks_for(imputation = list())$missing_values))
})

test_that("each supported imputation method gets its own description", {
    expect_match(blocks_for()$missing_values, "shifted down by 1.8", fixed = TRUE)
    expect_match(blocks_for()$missing_values, "narrowed to 0.3", fixed = TRUE)

    expect_match(blocks_for(imputation = list(method = "dep2"))$missing_values,
                 "1st percentile", fixed = TRUE)
    expect_match(blocks_for(imputation = list(method = "dep2",
                                              dep2_method = "MinProb"))$missing_values,
                 "narrow normal distribution centred", fixed = TRUE)
    expect_match(blocks_for(imputation = list(method = "qrilc"))$missing_values,
                 "left-censored data (QRILC)", fixed = TRUE)
    expect_match(blocks_for(imputation = list(method = "minval"))$missing_values,
                 "minimum-based value", fixed = TRUE)

    expect_null(blocks_for(imputation = list(method = "none"))$missing_values)
})

test_that("the retired KNN and MinProb-as-method claims cannot come back", {
    # The old text said perseus_like used "KNN for MCAR", and switched on
    # "minprob"/"knn" as if they were dispatcher methods. Neither is true:
    # 03_imputation.R has no KNN at all, and minprob is imputation$dep2_method.
    for (m in list(list(), list(method = "dep2"), list(method = "dep2", dep2_method = "MinProb"),
                   list(method = "qrilc"), list(method = "minval"))) {
        txt <- blocks_for(imputation = m)$missing_values
        expect_false(grepl("KNN", txt, fixed = TRUE))
        expect_false(grepl("MCAR", txt, fixed = TRUE))
        expect_false(grepl("MNAR", txt, fixed = TRUE))
    }
})

test_that("the measured/model-input sentence is narrow", {
    txt <- blocks_for()$missing_values
    expect_match(txt, "Where measured and model-input values are both presented, they are labelled separately",
                 fixed = TRUE)
    expect_false(grepl("throughout this report", txt, fixed = TRUE))
})


# =============================================================================
# Differential protein abundance
# =============================================================================

test_that("the reader-facing term is differential protein abundance", {
    expect_match(blocks_for()$results_oneliner, "Differential protein abundance", fixed = TRUE)
    expect_false(grepl("differential expression", blocks_for()$differential, ignore.case = TRUE))
})

test_that("each supported DE method is described from its implementation", {
    expect_match(blocks_for()$differential, "empirical Bayes shrinkage", fixed = TRUE)

    pc <- blocks_for(de = list(method = "limma_percontrast", p_cutoff = 0.05,
                               linear_fc_cutoff = 1.5),
                     de_method = "limma_percontrast")$differential
    expect_match(pc, "fitted separately", fixed = TRUE)
    expect_match(pc, "observed above the imputation floor", fixed = TRUE)

    tt <- blocks_for(de = list(method = "ttest", p_cutoff = 0.05, linear_fc_cutoff = 1.5),
                     de_method = "ttest")$differential
    expect_match(tt, "assuming equal variances", fixed = TRUE)

    we <- blocks_for(de = list(method = "welch", p_cutoff = 0.05, linear_fc_cutoff = 1.5),
                     de_method = "welch")$differential
    expect_match(we, "assuming unequal variances", fixed = TRUE)
})

test_that("a paired t-test is not described as an ordinary two-sample test", {
    txt <- blocks_for(de = list(method = "ttest", paired = TRUE, pairing_col = "Subject",
                                p_cutoff = 0.05, linear_fc_cutoff = 1.5),
                      de_method = "ttest")$differential
    expect_match(txt, "Paired t-tests", fixed = TRUE)
    expect_match(txt, "matched on Subject", fixed = TRUE)
    expect_false(grepl("assuming equal variances", txt, fixed = TRUE))
})

test_that("limma blocking is described as duplicateCorrelation, not a random block", {
    txt <- blocks_for(de = list(method = "limma", block_col = "Donor",
                                p_cutoff = 0.05, linear_fc_cutoff = 1.5))$differential
    expect_match(txt, "limma::duplicateCorrelation()", fixed = TRUE)
    expect_match(txt, "repeated measurements within Donor", fixed = TRUE)
    expect_false(grepl("random block", txt, fixed = TRUE))

    expect_false(grepl("duplicateCorrelation", blocks_for()$differential, fixed = TRUE))
})

test_that("ANOVA states the Tukey path without implying an omnibus gate", {
    txt <- blocks_for(de = list(method = "anova", p_cutoff = 0.05, linear_fc_cutoff = 1.5),
                      de_method = "anova")$differential
    expect_match(txt, "one-way ANOVA", fixed = TRUE)
    expect_match(txt, "Tukey's honest significant difference test", fixed = TRUE)
    expect_match(txt, "Those Tukey p-values were then adjusted across proteins", fixed = TRUE)
    # No claim that the omnibus p-value gates the pairwise calls: it does not.
    expect_false(grepl("omnibus", txt, fixed = TRUE))
    # And not the generic sentence as well, which would state adjustment twice.
    expect_false(grepl("P-values were adjusted with the", txt, fixed = TRUE))
})

test_that("the adjustment procedure is read from config, never hard-coded", {
    expect_match(blocks_for()$differential, "adjusted with the BH procedure", fixed = TRUE)
    txt <- blocks_for(de = list(method = "limma", p_adjust_method = "holm",
                                p_cutoff = 0.05, linear_fc_cutoff = 1.5))$differential
    expect_match(txt, "adjusted with the holm procedure", fixed = TRUE)
    expect_false(grepl("Benjamini-Hochberg", txt, fixed = TRUE))
})

test_that("fdrtool is claimed only where the implementation applies it", {
    # run_limma_proteomics() and run_limma_percontrast_proteomics() apply it;
    # the t-test and ANOVA paths never read the flag.
    on <- blocks_for(de = list(method = "limma", fdrtool_correction = TRUE,
                               p_cutoff = 0.05, linear_fc_cutoff = 1.5))$differential
    expect_match(on, "empirical null distribution", fixed = TRUE)

    for (m in c("ttest", "welch", "anova")) {
        txt <- blocks_for(de = list(method = m, fdrtool_correction = TRUE,
                                    p_cutoff = 0.05, linear_fc_cutoff = 1.5),
                          de_method = m)$differential
        expect_false(grepl("empirical null", txt, fixed = TRUE))
    }
})

test_that("the retired RNA-seq method descriptions cannot come back", {
    # limma-voom and DESeq2 cannot run in the proteomics DE module.
    for (m in c("limma", "limma_percontrast", "ttest", "welch", "anova")) {
        txt <- blocks_for(de = list(method = m, p_cutoff = 0.05, linear_fc_cutoff = 1.5),
                          de_method = m)$differential
        expect_false(grepl("voom", txt, ignore.case = TRUE))
        expect_false(grepl("DESeq2", txt, fixed = TRUE))
        expect_false(grepl("negative binomial", txt, fixed = TRUE))
    }
})

test_that("the count sentence says what the number is", {
    txt <- build_proteomics_methods_text(cfg_for(), de_method = "limma",
                                         n_included = 9000L)$differential
    expect_match(txt, "9,000 proteins passed filtering and were included", fixed = TRUE)
    expect_false(grepl("were quantified", txt, fixed = TRUE))
})

test_that("the count is not called pipeline-filtered in precomputed mode", {
    # mod_proteomics_de() returns before the filtered matrix is used, so the rows
    # counted here are the ones the upstream tables carried. Saying they passed
    # this pipeline's filtering credits it with a step it did not run.
    txt <- build_proteomics_methods_text(cfg_for(), de_method = "precomputed",
                                         n_included = 9000L)$differential
    expect_match(txt, "9,000 proteins were present in the precomputed result tables",
                 fixed = TRUE)
    expect_false(grepl("passed filtering", txt, fixed = TRUE))
    expect_false(grepl("included in the differential analysis", txt, fixed = TRUE))
})


# =============================================================================
# Precomputed DE
# =============================================================================

test_that("precomputed mode describes no internal model", {
    # mod_proteomics_de() returns method "precomputed" before any internal
    # imputation or fitting, so none of it may be claimed.
    b <- blocks_for(de_method = "precomputed")
    txt <- b$differential
    expect_match(txt, "loaded from precomputed result tables", fixed = TRUE)
    expect_match(txt, "not inferred by this pipeline", fixed = TRUE)
    for (bad in c("empirical Bayes", "limma", "t-tests", "ANOVA", "Tukey", "empirical null")) {
        expect_false(grepl(bad, txt, fixed = TRUE))
    }
    expect_null(b$consensus)
})

test_that("precomputed adjusted p-values are described conditionally", {
    # load_precomputed_proteomics_de() takes the supplied adjusted p-value and
    # falls back to p.adjust(method = "BH") -- specifically BH, not
    # de$p_adjust_method -- when the table carries none.
    txt <- blocks_for(de = list(method = "limma", p_adjust_method = "holm",
                                p_cutoff = 0.05, linear_fc_cutoff = 1.5),
                      de_method = "precomputed")$differential
    expect_match(txt, "uses the supplied adjusted p-value when available", fixed = TRUE)
    expect_match(txt, "Benjamini-Hochberg procedure", fixed = TRUE)
    expect_false(grepl("holm", txt, fixed = TRUE))
})

test_that("precomputed mode does not claim imputation was QC-only", {
    # With batch correction enabled the single-imputed matrix also carries the
    # correction, so "for visualisation only" would be false.
    txt <- blocks_for(de_method = "precomputed")$missing_values
    expect_match(txt, "not used to fit the differential-abundance model", fixed = TRUE)
    expect_match(txt, "when enabled, batch correction", fixed = TRUE)
    expect_false(grepl("only;", txt, fixed = TRUE))
})

test_that("the Results pointer distinguishes the two modes", {
    internal <- blocks_for()$results_oneliner
    pre <- blocks_for(de_method = "precomputed")$results_oneliner
    expect_match(internal, "was assessed with", fixed = TRUE)
    expect_match(pre, "loaded from precomputed tables", fixed = TRUE)
    expect_false(grepl("assessed with", pre, fixed = TRUE))
})


# =============================================================================
# Multiple-imputation consensus
# =============================================================================

multi_cfg <- function(method = "perseus_like", ...) {
    blocks_for(imputation = utils::modifyList(
        list(method = method, width = 0.3, downshift = 1.8, multi_imputation = TRUE,
             no_repetitions = 10, min_no_passed = 8), list(...)))
}

test_that("the consensus section appears only for repeated runs", {
    expect_null(blocks_for()$consensus)
    expect_null(blocks_for(imputation = list(method = "perseus_like", multi_imputation = TRUE,
                                             no_repetitions = 1, min_no_passed = 1))$consensus)
    expect_false(is.null(multi_cfg()$consensus))
})

test_that("the consensus text carries the whole contract", {
    txt <- multi_cfg()$consensus
    expect_match(txt, "repeated across 10 configured imputation runs", fixed = TRUE)
    expect_match(txt, "passed in at least 8 of the 10 runs", fixed = TRUE)
    expect_match(txt, "summarised adjusted p-value also met the cutoff", fixed = TRUE)
    expect_match(txt, "0.8 quantile of the per-run values", fixed = TRUE)
    expect_match(txt, "mean of the per-run ratios on the linear scale", fixed = TRUE)
    expect_match(txt, "logarithm of that mean, not the mean of the per-run log2 fold changes",
                 fixed = TRUE)
})

test_that("repetitions are not called independent when they are identical", {
    # make_imputations_proteomics() seeds each run, but a deterministic method
    # returns the same matrix every time.
    expect_match(multi_cfg("perseus_like")$consensus, "own seed", fixed = TRUE)
    for (m in c("minval", "none")) {
        expect_match(multi_cfg(m)$consensus, "identical by construction", fixed = TRUE)
    }
    expect_match(multi_cfg("dep2", dep2_method = "MinDet")$consensus,
                 "identical by construction", fixed = TRUE)
    expect_match(multi_cfg("dep2", dep2_method = "MinProb")$consensus,
                 "own seed", fixed = TRUE)

    for (m in c("perseus_like", "minval", "none")) {
        expect_false(grepl("independent imputations", multi_cfg(m)$consensus, fixed = TRUE))
    }
})

test_that("the per-run p-value follows use_adj_for_pass1 as the summariser reads it", {
    # summarize_limma_mult_imputation() uses isTRUE(), so an absent key means
    # the raw p-value was used.
    expect_match(multi_cfg()$consensus, "its unadjusted p-value met the cutoff", fixed = TRUE)
    txt <- build_proteomics_methods_text(
        cfg_for(imputation = list(method = "perseus_like", multi_imputation = TRUE,
                                  no_repetitions = 10, min_no_passed = 8),
                de = list(method = "limma", use_adj_for_pass1 = TRUE,
                          p_cutoff = 0.05, linear_fc_cutoff = 1.5)))$consensus
    expect_match(txt, "its adjusted p-value met the cutoff", fixed = TRUE)
})


# =============================================================================
# Quality control and software
# =============================================================================

test_that("the clustering sentence is specific to hierarchical clustering", {
    # partition and binary_patterns have different contracts; distance/linkage
    # describe the hierarchical step only.
    expect_false(grepl("linkage", blocks_for()$quality_control, fixed = TRUE))

    d <- withr::local_tempdir()
    hm <- existing_file(d, "Hierarchical_DE_heatmap.png")
    txt <- blocks_for(clustering = list(
        enabled = TRUE,
        steps = list(hierarchical = list(enabled = TRUE, distance = "manhattan",
                                         linkage = "average"))),
        hier_file = hm)$quality_control
    expect_match(txt, "When hierarchical clustering was enabled", fixed = TRUE)
    expect_match(txt, "manhattan distance and average linkage", fixed = TRUE)

    off <- blocks_for(clustering = list(
        enabled = TRUE,
        steps = list(hierarchical = list(enabled = FALSE))),
        hier_file = hm)$quality_control
    expect_false(grepl("linkage", off, fixed = TRUE))
})

test_that("the clustering sentence needs the heatmap, not just the flags", {
    # mod_proteomics_clustering() returns before the hierarchical step when
    # fewer than two DE features are in the matrix, and says so with a message
    # rather than an error. On such a run the flags are on and no heatmap
    # exists, so the flags alone must not produce the claim.
    d <- withr::local_tempdir()
    on_cfg <- list(enabled = TRUE,
                   steps = list(hierarchical = list(enabled = TRUE,
                                                    distance = "euclidean",
                                                    linkage = "complete")))
    expect_false(grepl("linkage", blocks_for(clustering = on_cfg)$quality_control,
                       fixed = TRUE))
    expect_false(grepl("linkage",
                       blocks_for(clustering = on_cfg,
                                  hier_file = file.path(d, "absent.png"))$quality_control,
                       fixed = TRUE))
    expect_match(blocks_for(clustering = on_cfg,
                            hier_file = existing_file(d, "Hierarchical_DE_heatmap.png")
                            )$quality_control,
                 "euclidean distance and complete linkage", fixed = TRUE)
})

test_that("the complete-case PCA keeps its sensitivity-view framing", {
    d <- withr::local_tempdir()
    txt <- blocks_for(pca_cc_file = existing_file(d, "PCA_robust.png"))$quality_control
    expect_match(txt, "observed in every included sample is shown as a sensitivity view",
                 fixed = TRUE)
    for (bad in c("pre-imputation", "imputation-free", "before imputation")) {
        expect_false(grepl(bad, txt, fixed = TRUE))
    }
})

test_that("the sensitivity view is claimed only when the panel exists", {
    # 01_mod_qc_pre.R skips the panel with no missingness flags, with fewer than
    # three complete-case proteins, or when the plot call fails -- and the
    # template additionally requires the file. A Methods claim for a figure the
    # report does not contain sends the reader looking for nothing.
    d <- withr::local_tempdir()
    expect_false(grepl("sensitivity view", blocks_for()$quality_control, fixed = TRUE))
    expect_false(grepl("sensitivity view",
                       blocks_for(pca_cc_file = file.path(d, "absent.png"))$quality_control,
                       fixed = TRUE))
    # The PCA sentence itself is unconditional; only the panel claim is gated.
    expect_match(blocks_for()$quality_control,
                 "Principal component analysis was computed", fixed = TRUE)
})

test_that("pathway and PPI sentences survive, still gated on their flags", {
    # These were reader-facing Methods content before this change and are kept.
    expect_null(blocks_for()$downstream)

    pw <- blocks_for(pathway = list(enabled = TRUE, method = "ora",
                                    databases = c("GO", "KEGG")))$downstream
    expect_match(pw, "over-representation analysis", fixed = TRUE)
    expect_match(pw, "against GO, KEGG", fixed = TRUE)
    expect_false(grepl("gene-set enrichment", pw, fixed = TRUE))

    both <- blocks_for(pathway = list(enabled = TRUE))$downstream
    expect_match(both, "gene-set enrichment analysis and over-representation analysis",
                 fixed = TRUE)

    ppi <- blocks_for(ppi = list(enabled = TRUE))$downstream
    expect_match(ppi, "STRING database", fixed = TRUE)
    expect_false(grepl("Pathway enrichment", ppi, fixed = TRUE))
})

test_that("the version identifier comes from provenance, never a literal", {
    plain <- blocks_for()$software
    expect_match(plain, "Analyses were performed in R", fixed = TRUE)
    expect_false(grepl("v1.0", plain, fixed = TRUE))
    expect_false(grepl("commit", plain, fixed = TRUE))

    d <- tempfile("methods-prov-"); dir.create(file.path(d, "execution_info"), recursive = TRUE)
    on.exit(unlink(d, recursive = TRUE), add = TRUE)
    writeLines("0123456789abcdef0123456789abcdef01234567",
               file.path(d, "execution_info", "git_commit.txt"))
    expect_match(build_proteomics_methods_text(cfg_for(), run_dir = d)$software,
                 "(commit 01234567)", fixed = TRUE)
})


# =============================================================================
# The template consumes the generator rather than duplicating it
# =============================================================================

test_that("the report template has one source of Methods wording", {
    candidates <- c(
        testthat::test_path("..", "..", "R", "domain", "proteomics",
                            "report_template_proteomics.Rmd"),
        "R/domain/proteomics/report_template_proteomics.Rmd"
    )
    f <- candidates[file.exists(candidates)][1]
    skip_if(is.na(f), "proteomics report template not found")
    src <- paste(readLines(f, warn = FALSE), collapse = "\n")

    expect_match(src, "build_proteomics_methods_text(", fixed = TRUE)
    # Both override keys keep working; this PR does not change that contract.
    expect_match(src, "report_cfg$stats_section_text", fixed = TRUE)
    expect_match(src, "report_cfg$stats_methods_text", fixed = TRUE)
    # The duplicated Copy-ready block and its hard-coded version are gone.
    expect_false(grepl("Copy-ready version", src, fixed = TRUE))
    expect_false(grepl("multiomics-core pipeline (v1.0)", src, fixed = TRUE))
    # And the second switch() that drifted from the Methods one.
    expect_false(grepl("limma-voom", src, fixed = TRUE))
    expect_false(grepl("negative binomial generalized linear models", src, fixed = TRUE))
})

test_that("the statistics override replaces statistics only", {
    # Before the generator existed, the pathway and PPI Methods sentences were
    # emitted from the Methods chunk with no override gate. Folding every block
    # into the else-branch of stats_section_text would have deleted them for any
    # project that supplies its own statistics text.
    f <- c(testthat::test_path("..", "..", "R", "domain", "proteomics",
                               "report_template_proteomics.Rmd"),
           "R/domain/proteomics/report_template_proteomics.Rmd")
    f <- f[file.exists(f)][1]
    skip_if(is.na(f), "proteomics report template not found")
    lines <- readLines(f, warn = FALSE)

    start <- grep("^```\\{r methods-statistical-analysis[ ,}]", lines)
    expect_length(start, 1L)
    end <- start + which(grepl("^```\\s*$", lines[(start + 1):length(lines)]))[1]
    chunk <- lines[(start + 1):(end - 1)]

    # The override's else-branch closes before the blocks that are not
    # statistics, so those are emitted on both paths.
    else_close <- grep("^\\}\\s*$", chunk)
    expect_gt(length(else_close), 0L)
    tail_txt <- paste(chunk[(else_close[1] + 1):length(chunk)], collapse = "\n")
    for (blk in c("quality_control", "downstream", "software")) {
        expect_match(tail_txt, sprintf("methods_blocks$%s", blk), fixed = TRUE)
    }
    # while the statistical blocks stay inside it.
    head_txt <- paste(chunk[1:else_close[1]], collapse = "\n")
    for (blk in c("differential", "consensus")) {
        expect_match(head_txt, sprintf("methods_blocks$%s", blk), fixed = TRUE)
    }
})

test_that("the Results pointer is not emitted when Methods is hidden", {
    # report$show_methods gates the Methods chunks but not this one, so the
    # one-liner's "see Methods" would point at a section the report omits.
    f <- c(testthat::test_path("..", "..", "R", "domain", "proteomics",
                               "report_template_proteomics.Rmd"),
           "R/domain/proteomics/report_template_proteomics.Rmd")
    f <- f[file.exists(f)][1]
    skip_if(is.na(f), "proteomics report template not found")
    lines <- readLines(f, warn = FALSE)

    start <- grep("^```\\{r de-stats-description[ ,}]", lines)
    expect_length(start, 1L)
    end <- start + which(grepl("^```\\s*$", lines[(start + 1):length(lines)]))[1]
    chunk <- paste(lines[(start + 1):(end - 1)], collapse = "\n")

    expect_match(chunk, "isTRUE(show_methods)", fixed = TRUE)
    expect_match(chunk, ".ds_blocks$results_oneliner", fixed = TRUE)
    # and the full description is what replaces it, not silence.
    expect_match(chunk, ".ds_blocks$differential", fixed = TRUE)

    # The pointer itself lives only in the one-liner.
    expect_match(blocks_for()$results_oneliner, "see Methods", fixed = TRUE)
    expect_false(grepl("see Methods", blocks_for()$differential, fixed = TRUE))
})
