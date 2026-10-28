# tests/testthat/test-proteomics-de-batch-covariate.R
#
# batch_correction$de_mode = "covariate": batch as a fixed additive term in the
# pooled limma design, fitted on uncorrected values.
#   - resolve_de_batch_covariate()  (R/domain/proteomics/05_de_summary.R)
#   - build_proteomics_de_design()  (R/domain/proteomics/05_de_summary.R)
#   - run_limma_proteomics()        (R/domain/proteomics/05_de_summary.R)
#   - validate_proteomics_config()  (R/domain/proteomics/90_config_validate.R)
#
# All data is synthetic.

# Two conditions crossed with two batches, three replicates per cell. Batch B2
# carries a large additive shift; the first 30 proteins carry a small, real
# treatment effect that the batch shift buries unless batch is modelled.
make_batch_fixture <- function(n_prot = 300L, batch_shift = 3, effect = 0.8) {
    samples <- sprintf("S%02d", 1:12)
    group   <- rep(c("Ctrl", "Trt"), times = 6)
    batch   <- rep(c("B1", "B2"), each = 6)

    expr <- withr::with_seed(7, {
        m <- matrix(rnorm(n_prot * 12, mean = 20, sd = 0.3), nrow = n_prot)
        m[, batch == "B2"] <- m[, batch == "B2"] + batch_shift
        m[1:30, group == "Trt"] <- m[1:30, group == "Trt"] + effect
        m
    })
    rownames(expr) <- paste0("P", seq_len(n_prot))
    colnames(expr) <- samples

    meta <- data.frame(SampleName = samples, Group = group, Batch = batch,
                       stringsAsFactors = FALSE)
    prot_tbl <- data.frame(
        Protein.Group = rownames(expr), Protein.Names = rownames(expr),
        Genes = rownames(expr), First.Protein.Description = rownames(expr),
        stringsAsFactors = FALSE
    )
    contrasts_df <- data.frame(Contrast_name = "Trt_vs_Ctrl", Factor = "Group",
                               Numerator = "Trt", Denominator = "Ctrl",
                               stringsAsFactors = FALSE)
    list(expr = expr, meta = meta, prot_tbl = prot_tbl, contrasts_df = contrasts_df)
}

batch_cfg <- function(de_mode = "covariate", batch_col = "Batch", block_col = NULL) {
    list(modes = list(proteomics = list(
        effects    = list(samples = "SampleName", color = "Group"),
        id_columns = list(
            protein_id    = "Protein.Group",
            sample_col    = "SampleName",
            protein_annot = c("Protein.Group", "Protein.Names", "Genes",
                              "First.Protein.Description")
        ),
        de_table   = list(id_col = "FeatureID"),
        de         = list(method = "limma", p_adjust_method = "BH", block_col = block_col),
        batch_correction = list(method = "none", batch_col = batch_col, de_mode = de_mode)
    )))
}

run_fixture <- function(fx, cfg) {
    run_limma_proteomics(expr_imp = fx$expr, meta = fx$meta,
                         contrasts_df = fx$contrasts_df, prot_tbl = fx$prot_tbl,
                         cfg = cfg)
}

# --- resolver -------------------------------------------------------------

test_that("resolve_de_batch_covariate is NULL unless de_mode is covariate", {
    meta <- make_batch_fixture()$meta
    expect_null(resolve_de_batch_covariate(meta, batch_cfg("corrected")$modes$proteomics))
    no_bc <- batch_cfg()$modes$proteomics
    no_bc$batch_correction <- NULL
    expect_null(resolve_de_batch_covariate(meta, no_bc))
})

test_that("resolve_de_batch_covariate returns batch aligned to the metadata rows", {
    meta <- make_batch_fixture()$meta
    b <- resolve_de_batch_covariate(meta, batch_cfg()$modes$proteomics)
    expect_s3_class(b, "factor")
    expect_equal(as.character(b), meta$Batch)
})

test_that("resolve_de_batch_covariate refuses unusable batch columns", {
    meta <- make_batch_fixture()$meta
    expect_error(resolve_de_batch_covariate(meta, batch_cfg(batch_col = NULL)$modes$proteomics),
                 "no batch column")
    expect_error(resolve_de_batch_covariate(meta, batch_cfg(batch_col = "Plate")$modes$proteomics),
                 "not in the sample metadata")

    na_meta <- meta
    na_meta$Batch[3] <- NA
    expect_error(resolve_de_batch_covariate(na_meta, batch_cfg()$modes$proteomics),
                 "missing values")

    one <- meta
    one$Batch <- "B1"
    expect_error(resolve_de_batch_covariate(one, batch_cfg()$modes$proteomics),
                 "single level")

    expect_error(resolve_de_batch_covariate(meta, batch_cfg(block_col = "Batch")$modes$proteomics),
                 "both de\\$block_col")
})

test_that("an unknown de_mode is refused by the shared resolver", {
    expect_error(get_proteomics_batch_config(batch_cfg("both")$modes$proteomics),
                 "de_mode")
})

# --- design ---------------------------------------------------------------

test_that("without batch the design is one column per condition, as before", {
    g <- factor(c("Ctrl", "Ctrl", "Trt", "Trt"))
    d <- build_proteomics_de_design(g)
    expect_equal(colnames(d), c("Ctrl", "Trt"))
    expect_equal(unname(d[, "Trt"]), c(0, 0, 1, 1))
})

test_that("with batch the condition columns keep their names and batch is treatment-coded", {
    fx <- make_batch_fixture()
    d <- build_proteomics_de_design(factor(fx$meta$Group), factor(fx$meta$Batch))
    expect_equal(colnames(d), c("Ctrl", "Trt", "batch_B2"))
    expect_equal(unname(d[, "batch_B2"]), as.numeric(fx$meta$Batch == "B2"))
})

test_that("a batch column that clashes with a condition name never renames the condition", {
    g <- factor(c("batch_B2", "batch_B2", "X", "X", "batch_B2", "X"))
    b <- factor(c("B1", "B1", "B1", "B2", "B2", "B2"))
    d <- build_proteomics_de_design(g, b)
    expect_equal(colnames(d)[1:2], c("batch_B2", "X"))
    expect_false(anyDuplicated(colnames(d)) > 0)
})

test_that("batch confounded with condition is refused, not fitted", {
    g <- factor(c("Ctrl", "Ctrl", "Ctrl", "Trt", "Trt", "Trt"))
    b <- factor(c("B1", "B1", "B1", "B2", "B2", "B2"))
    expect_error(build_proteomics_de_design(g, b), "confounded")
})

test_that("a design with no residual degrees of freedom is refused", {
    g <- factor(c("Ctrl", "Trt", "Ctrl", "Trt"))
    b <- factor(c("B1", "B1", "B2", "B3"))
    expect_error(build_proteomics_de_design(g, b), "no residual degrees of freedom")
})

# --- the fit --------------------------------------------------------------

test_that("modelling batch recovers an effect that the batch shift buries", {
    skip_if_not_installed("limma")
    fx <- make_batch_fixture()

    plain   <- run_fixture(fx, batch_cfg("corrected"))
    covaried <- run_fixture(fx, batch_cfg("covariate"))

    expect_equal(colnames(covaried$design), c("Ctrl", "Trt", "batch_B2"))
    expect_equal(colnames(plain$design), c("Ctrl", "Trt"))
    # One more coefficient, one fewer residual df: limma is charged for batch.
    expect_equal(unique(covaried$fit2$df.residual), unique(plain$fit2$df.residual) - 1)

    n_true_hits <- function(res) {
        de <- res$de_tables[["Trt_vs_Ctrl"]]
        sum(de$adj.P.Val[de$FeatureID %in% paste0("P", 1:30)] <= 0.05)
    }
    expect_gt(n_true_hits(covaried), n_true_hits(plain))

    # The contrast still estimates the treatment effect, not a batch artefact.
    de <- covaried$de_tables[["Trt_vs_Ctrl"]]
    expect_equal(mean(de$logFC[1:30]), 0.8, tolerance = 0.1)
})

test_that("tables keep their shape and contrast names with batch in the model", {
    skip_if_not_installed("limma")
    fx <- make_batch_fixture()
    res <- run_fixture(fx, batch_cfg("covariate"))
    expect_named(res$de_tables, "Trt_vs_Ctrl")
    expect_equal(res$de_tables[["Trt_vs_Ctrl"]]$FeatureID, rownames(fx$expr))
    expect_equal(colnames(res$contrast_matrix), "Trt_vs_Ctrl")
    expect_equal(unname(res$contrast_matrix["batch_B2", 1]), 0)
})

# --- config validation ----------------------------------------------------

test_that("validation accepts covariate mode with pooled limma", {
    cfg <- batch_cfg()$modes$proteomics
    expect_true(validate_proteomics_config(cfg))
})

test_that("validation refuses covariate mode where batch cannot enter the model", {
    for (m in c("limma_percontrast", "ttest", "welch", "anova")) {
        cfg <- batch_cfg()$modes$proteomics
        cfg$de$method <- m
        expect_error(validate_proteomics_config(cfg), "only with de\\$method 'limma'")
    }

    pre <- batch_cfg()$modes$proteomics
    pre$files <- list(de_table = list(Trt_vs_Ctrl = "de.tsv"))
    expect_error(validate_proteomics_config(pre), "precomputed")

    expect_error(validate_proteomics_config(batch_cfg(batch_col = NULL)$modes$proteomics),
                 "no batch column")
    expect_error(validate_proteomics_config(batch_cfg(block_col = "Batch")$modes$proteomics),
                 "Model it once")
})
