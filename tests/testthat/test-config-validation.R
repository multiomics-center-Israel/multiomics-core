# Tests for config validation functions
# test-config-validation.R

# --- Proteomics Config Mock ---

create_mock_proteomics_config <- function() {
    list(
        engine = "maxquant",
        scale_in = "log2",
        id_columns = list(
            protein_id = "Protein.IDs",
            sample_col = "SampleID"
        ),
        imputation = list(
            method = "perseus_like",
            no_repetitions = 5,
            min_no_passed = 3,
            width = 0.3,
            downshift = 1.8
        ),
        de = list(
            method = "limma",
            p_cutoff = 0.05,
            linear_fc_cutoff = 1.5,
            p_adjust_method = "BH"
        ),
        pathway = list(
            gsea_ranking = "stat"
        ),
        ppi = list(
            enabled = TRUE,
            significance_threshold = 0.05,
            string_score_threshold = 400,
            lfc_threshold = 0,
            active_subnetwork = FALSE,
            complex_analysis = FALSE
        ),
        advanced_stats = list(
            enabled = TRUE,
            bootstrap_n = 1000,
            ci_level = 0.95,
            run_robust_regression = FALSE,
            run_power_analysis = TRUE
        )
    )
}

test_that("validate_proteomics_config passes with valid config", {
    cfg <- create_mock_proteomics_config()
    expect_true(validate_proteomics_config(cfg))
})

test_that("validate_proteomics_config errors on missing protein_id", {
    cfg <- create_mock_proteomics_config()
    cfg$id_columns$protein_id <- NULL
    expect_error(validate_proteomics_config(cfg))
})

test_that("validate_proteomics_config errors on missing sample_col", {
    cfg <- create_mock_proteomics_config()
    cfg$id_columns$sample_col <- NULL
    expect_error(validate_proteomics_config(cfg))
})

test_that("validate_proteomics_config errors on invalid imputation method", {
    cfg <- create_mock_proteomics_config()
    cfg$imputation$method <- "invalid_method"
    expect_error(validate_proteomics_config(cfg))
})

test_that("validate_proteomics_config errors on p_cutoff > 1", {
    cfg <- create_mock_proteomics_config()
    cfg$de$p_cutoff <- 1.5
    expect_error(validate_proteomics_config(cfg))
})

test_that("validate_proteomics_config errors on p_cutoff < 0", {
    cfg <- create_mock_proteomics_config()
    cfg$de$p_cutoff <- -0.1
    expect_error(validate_proteomics_config(cfg))
})

test_that("validate_proteomics_config passes with valid PPI config", {
    cfg <- create_mock_proteomics_config()
    expect_true(validate_proteomics_config(cfg))
})

test_that("validate_proteomics_config errors on PPI significance_threshold > 1", {
    cfg <- create_mock_proteomics_config()
    cfg$ppi$significance_threshold <- 1.5
    expect_error(validate_proteomics_config(cfg))
})

test_that("validate_proteomics_config errors on PPI string_score_threshold > 1000", {
    cfg <- create_mock_proteomics_config()
    cfg$ppi$string_score_threshold <- 1500
    expect_error(validate_proteomics_config(cfg))
})

test_that("validate_proteomics_config passes with valid advanced_stats config", {
    cfg <- create_mock_proteomics_config()
    expect_true(validate_proteomics_config(cfg))
})

test_that("validate_proteomics_config errors on bootstrap_n < 100", {
    cfg <- create_mock_proteomics_config()
    cfg$advanced_stats$bootstrap_n <- 50
    expect_error(validate_proteomics_config(cfg))
})

test_that("validate_proteomics_config errors on ci_level > 0.999", {
    cfg <- create_mock_proteomics_config()
    cfg$advanced_stats$ci_level <- 1.0
    expect_error(validate_proteomics_config(cfg))
})

test_that("validate_proteomics_config allows NULL optional sections", {
    cfg <- create_mock_proteomics_config()
    cfg$ppi <- NULL
    cfg$advanced_stats <- NULL
    cfg$pathway <- NULL
    expect_true(validate_proteomics_config(cfg))
})

test_that("validate_proteomics_config errors on min_no_passed > no_repetitions", {
    cfg <- create_mock_proteomics_config()
    cfg$imputation$min_no_passed <- 10
    expect_error(validate_proteomics_config(cfg))
})

test_that("validate_proteomics_config errors on invalid DE method", {
    cfg <- create_mock_proteomics_config()
    cfg$de$method <- "invalid_de"
    expect_error(validate_proteomics_config(cfg))
})

test_that("validate_proteomics_config errors on invalid p_adjust_method", {
    cfg <- create_mock_proteomics_config()
    cfg$de$p_adjust_method <- "invalid_padj"
    expect_error(validate_proteomics_config(cfg))
})

test_that("validate_proteomics_config errors on invalid scale_in", {
    cfg <- create_mock_proteomics_config()
    cfg$scale_in <- "log10"
    expect_error(validate_proteomics_config(cfg))
})

# --- RNA-seq Config ---

create_mock_rna_config <- function() {
    list(
        id_columns = list(
            gene_id = "gene_id"
        ),
        files = list(
            counts = "counts.csv",
            metadata = "metadata.csv"
        ),
        normalization = list(
            method = "TMMlogCPM"
        ),
        filtering = list(
            group_col = "Condition",
            threshold = 10
        ),
        de = list(
            padj_cutoff = 0.05,
            linear_fc_cutoff = 1.5
        )
    )
}

test_that("validate_rna_config passes with valid config", {
    cfg <- create_mock_rna_config()
    expect_true(validate_rna_config(cfg))
})

test_that("validate_rna_config errors on missing counts file", {
    cfg <- create_mock_rna_config()
    cfg$files$counts <- NULL
    expect_error(validate_rna_config(cfg))
})

test_that("validate_rna_config errors on invalid normalization method", {
    cfg <- create_mock_rna_config()
    cfg$normalization$method <- "invalid_norm"
    expect_error(validate_rna_config(cfg))
})

# --- Enrichment block validation (the `enrichment` section) ---

test_that("validate_rna_config accepts a valid enrichment block", {
    cfg <- create_mock_rna_config()
    cfg$enrichment <- list(
        enabled        = TRUE,
        annotation_dir = "data/Func_annot_data_expanded",
        databases      = c("GO_BP", "GO_MF", "GO_CC"),
        workers        = 4,
        pvalue_cutoff  = 0.05,
        padj_method    = "fdr",
        plots          = list(dotplot = TRUE, shared_genes = FALSE),
        gsea           = list(per_pathway_artifacts = TRUE)
    )
    expect_true(validate_rna_config(cfg))
})

test_that("validate_rna_config rejects an out-of-range enrichment pvalue_cutoff", {
    cfg <- create_mock_rna_config()
    cfg$enrichment <- list(pvalue_cutoff = 2)
    expect_error(validate_rna_config(cfg))
})

test_that("validate_rna_config rejects an unknown enrichment database", {
    cfg <- create_mock_rna_config()
    cfg$enrichment <- list(databases = c("GO_BP", "NOT_A_DB"))
    expect_error(validate_rna_config(cfg))
})

test_that("validate_rna_config rejects a non-boolean enrichment plot toggle", {
    cfg <- create_mock_rna_config()
    cfg$enrichment <- list(plots = list(dotplot = "yes"))
    expect_error(validate_rna_config(cfg))
})

# --- Multi-omics config ---

create_mock_multiomics_config <- function() {
    list(
        integration = list(
            methods = c("DIABLO", "SNF"),
            diablo = list(ncomp = 2),
            snf = list(K = 3, alpha = 0.5, T = 20, n_clusters = 2)
        ),
        feature_selection = list(method = "variance", top_n = 500),
        condition_column = "Group"
    )
}

test_that("validate_multiomics_config returns the config section, not TRUE", {
    # validate_config() assigns the return value back onto
    # config$modes$multiomics; returning TRUE wiped the whole section.
    cfg <- create_mock_multiomics_config()
    res <- validate_multiomics_config(cfg)
    expect_type(res, "list")
    expect_equal(res$condition_column, "Group")
    expect_equal(res$integration$methods, c("DIABLO", "SNF"))
})

test_that("validate_config keeps config$modes$multiomics a list", {
    config <- list(
        paths  = list(raw = "data", out = "outputs"),
        params = list(seed = 1),
        modes  = list(multiomics = create_mock_multiomics_config())
    )
    validated <- validate_config(config)
    expect_type(validated$modes$multiomics, "list")
    expect_equal(
        validated$modes$multiomics$feature_selection$top_n, 500
    )
})
