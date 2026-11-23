# tests/testthat/test-proteomics-covariate-removal.R
#
# The DE matrix with de$covariates removed (R/domain/proteomics/05e_covariate_removal.R):
#   - remove_covariate_effects(): refitting gives a zero covariate slope and the
#     adjusted group means equal the model's group coefficients
#   - resolve_proteomics_display_matrix(): the matrix the sample heatmaps and the
#     "after" PCA panel use is limma::removeBatchEffect() output; unchanged
#     without covariates
#   - build_group_mean_adjusted_proteomics(): Mean.adj differences reproduce the
#     fitted log2FC of run_limma_proteomics() (single imputation draw)
#
# Synthetic, seeded data only.

make_removal_fixture <- function(n_prot = 60L) {
    samples <- sprintf("S%02d", 1:12)
    group <- rep(c("Ctrl", "TrtA", "TrtB"), each = 4)
    # Covariate deliberately unbalanced across groups, so removing it moves the
    # raw group means and the test can tell adjusted from unadjusted.
    depth <- c(0.9, 1.1, 0.7, 1.3, -0.2, 0.1, -0.5, 0.3, -0.8, -0.4, -1.0, -0.6)
    expr <- withr::with_seed(7, {
        m <- matrix(rnorm(n_prot * 12, mean = 20, sd = 0.4), nrow = n_prot)
        slopes <- rnorm(n_prot, mean = 0.8, sd = 0.3)
        m <- m + outer(slopes, depth)
        m[1:10, group == "TrtA"] <- m[1:10, group == "TrtA"] + 1
        m[11:20, group == "TrtB"] <- m[11:20, group == "TrtB"] - 1
        m
    })
    rownames(expr) <- paste0("P", seq_len(n_prot))
    colnames(expr) <- samples
    meta <- data.frame(SampleName = samples, Group = group, depth = depth,
                       stringsAsFactors = FALSE)
    prot_tbl <- data.frame(Protein.Group = rownames(expr), Protein.Names = rownames(expr),
                           Genes = rownames(expr), First.Protein.Description = rownames(expr),
                           stringsAsFactors = FALSE)
    contrasts_df <- data.frame(
        Contrast_name = c("TrtA_vs_Ctrl", "TrtB_vs_Ctrl", "TrtB_vs_TrtA"),
        Factor = "Group", Numerator = c("TrtA", "TrtB", "TrtB"),
        Denominator = c("Ctrl", "Ctrl", "TrtA"), stringsAsFactors = FALSE)
    list(expr = expr, meta = meta, prot_tbl = prot_tbl, contrasts_df = contrasts_df)
}

removal_cfg <- function(covariates = list("depth")) {
    list(
        params = list(seed = 1),
        modes = list(proteomics = list(
            effects = list(samples = "SampleName", color = "Group"),
            id_columns = list(protein_id = "Protein.Group", sample_col = "SampleName",
                              protein_annot = c("Protein.Group", "Protein.Names", "Genes",
                                                "First.Protein.Description")),
            de_table = list(id_col = "FeatureID"),
            de = list(method = "limma", p_adjust_method = "BH", p_cutoff = 0.05,
                      linear_fc_cutoff = 1.5, covariates = covariates),
            imputation = list(method = "perseus_like", multi_imputation = FALSE),
            batch_correction = list(method = "none")
        ))
    )
}

test_that("refitting the adjusted matrix gives a zero covariate slope and the model's group means", {
    fx <- make_removal_fixture()
    g <- factor(fx$meta$Group)
    cov_c <- cbind(depth = fx$meta$depth - mean(fx$meta$depth))
    adj <- remove_covariate_effects(fx$expr, cov_c, g)

    design <- cbind(stats::model.matrix(~ 0 + g), depth = cov_c[, 1])
    refit <- limma::lmFit(adj, design)
    expect_lt(max(abs(refit$coefficients[, "depth"])), 1e-8)

    original <- limma::lmFit(fx$expr, design)
    for (lv in levels(g)) {
        adj_mean <- rowMeans(adj[, g == lv, drop = FALSE])
        expect_equal(unname(adj_mean), unname(original$coefficients[, paste0("g", lv)]),
                     tolerance = 1e-10)
    }
    # The raw group means are not preserved when the covariate differs by group.
    raw_ctrl <- rowMeans(fx$expr[, g == "Ctrl"])
    expect_gt(max(abs(raw_ctrl - rowMeans(adj[, g == "Ctrl"]))), 0.1)
})

test_that("the display matrix is removeBatchEffect output with the resolver's covariates", {
    fx <- make_removal_fixture()
    config <- removal_cfg()
    pre <- list(meta = fx$meta, expr_filt = fx$expr, expr_imp_single = fx$expr)
    de_res <- list(imputations = list(fx$expr), covariate_design = list(notes = character()))

    disp <- resolve_proteomics_display_matrix(pre, de_res, config)
    expected <- limma::removeBatchEffect(
        fx$expr, covariates = cbind(fx$meta$depth - mean(fx$meta$depth)),
        design = stats::model.matrix(~ 0 + factor(fx$meta$Group)))
    expect_equal(disp$mat, expected, tolerance = 1e-10)
    expect_equal(disp$before, fx$expr)
    expect_equal(colnames(disp$covariates), "depth")
    expect_match(disp$label, "covariate depth removed", fixed = TRUE)
    expect_true(all(disp$observed))
})

test_that("without a fitted covariate the display matrix is the DE matrix unchanged", {
    fx <- make_removal_fixture()
    pre <- list(meta = fx$meta, expr_filt = fx$expr, expr_imp_single = fx$expr + 5)
    de_res <- list(imputations = list(fx$expr), covariate_design = NULL)
    disp <- resolve_proteomics_display_matrix(pre, de_res, removal_cfg(covariates = NULL))
    expect_equal(disp$mat, fx$expr)
    expect_null(disp$covariates)
    expect_identical(disp$label, "DE matrix")
})

test_that("Mean.adj differences reproduce the fitted log2FC for every contrast", {
    fx <- make_removal_fixture()
    config <- removal_cfg()
    fit <- run_limma_proteomics(fx$expr, fx$meta, fx$contrasts_df, fx$prot_tbl, config)
    expect_false(is.null(fit$covariate_design))

    pre <- list(meta = fx$meta, expr_filt = fx$expr, expr_imp_single = fx$expr)
    de_res <- list(imputations = list(fx$expr), covariate_design = fit$covariate_design)
    adj_means <- build_group_mean_adjusted_proteomics(pre, de_res, fx$contrasts_df, config)
    expect_setequal(colnames(adj_means), c("Mean.adj.Ctrl", "Mean.adj.TrtA", "Mean.adj.TrtB"))

    for (i in seq_len(nrow(fx$contrasts_df))) {
        cn <- fx$contrasts_df$Contrast_name[i]
        tab <- fit$de_tables[[cn]]
        diff <- adj_means[tab$FeatureID, paste0("Mean.adj.", fx$contrasts_df$Numerator[i])] -
            adj_means[tab$FeatureID, paste0("Mean.adj.", fx$contrasts_df$Denominator[i])]
        expect_equal(unname(diff), unname(tab$logFC), tolerance = 1e-8)
    }
})

test_that("no Mean.adj columns without a fitted covariate", {
    fx <- make_removal_fixture()
    pre <- list(meta = fx$meta, expr_filt = fx$expr, expr_imp_single = fx$expr)
    de_res <- list(imputations = list(fx$expr), covariate_design = NULL)
    expect_null(build_group_mean_adjusted_proteomics(pre, de_res, fx$contrasts_df,
                                                     removal_cfg(covariates = NULL)))
})

test_that("the Excel notes name the design and the adjusted columns only with covariates", {
    plain <- build_provenance_notes("proteomics", multi_imputation = FALSE)
    expect_false("Mean.adj.<group>" %in% plain$glossary$Column)

    with_cov <- build_provenance_notes("proteomics", multi_imputation = FALSE,
                                       de_covariates = "depth",
                                       de_design = "~ 0 + Group + depth")
    cols <- with_cov$glossary$Column
    expect_true("Mean.adj.<group>" %in% cols)
    expect_equal(which(cols == "Mean.adj.<group>"), which(cols == "Mean.<group>") + 1L)
    expect_true(any(grepl("~ 0 + Group + depth", with_cov$notes, fixed = TRUE)))
    expect_match(with_cov$glossary$Meaning[cols == "Mean.<group>"], "Not adjusted", fixed = TRUE)
})
