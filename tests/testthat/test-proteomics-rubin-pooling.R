# tests/testthat/test-proteomics-rubin-pooling.R
#
# pool_limma_imputations() combines per-imputation limma fits with Rubin's
# rules (05d_rubin_pooled_stat.R). The fits here are real eBayes MArrayLM
# objects on a small synthetic matrix, so every quantity is checked against
# what limma itself produced.

# Three-vs-three toy data with missing values; `identical_draws` fills every
# draw with the same constant, otherwise each draw is a fresh random fill.
make_toy_runs <- function(m = 3, n_prot = 40, seed = 11, identical_draws = FALSE) {
    samples <- paste0("S", 1:6)
    withr::with_seed(seed, {
        base <- matrix(rnorm(n_prot * 6, mean = 20), n_prot,
                       dimnames = list(paste0("P", seq_len(n_prot)), samples))
        base[1:10, 4:6] <- base[1:10, 4:6] + 1.5
        na_mask <- matrix(FALSE, n_prot, 6, dimnames = dimnames(base))
        na_mask[sample(length(na_mask), 30)] <- TRUE
        draws <- lapply(seq_len(m), function(i) {
            x <- base
            x[na_mask] <- if (identical_draws) 17 else rnorm(sum(na_mask), 17, 0.5)
            x
        })
    })
    grp <- factor(rep(c("A", "B"), each = 3))
    design <- stats::model.matrix(~ 0 + grp)
    colnames(design) <- c("A", "B")
    cm <- limma::makeContrasts(B - A, levels = design)
    colnames(cm) <- "B_vs_A"
    runs <- lapply(draws, function(x) {
        fit2 <- limma::eBayes(limma::contrasts.fit(limma::lmFit(x, design), cm))
        tt <- limma::topTable(fit2, coef = "B_vs_A", sort.by = "none", number = Inf)
        list(
            expr_imp = x, design = design, contrast_matrix = cm, fit2 = fit2,
            de_tables = list(B_vs_A = data.frame(
                FeatureID = rownames(x), Genes = paste0("G", seq_len(nrow(x))),
                Contrast = "B_vs_A", tt[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")],
                stringsAsFactors = FALSE, row.names = NULL))
        )
    })
    list(runs = runs, observed = !na_mask)
}

test_that("the pooled statistics match a hand computation of Rubin's rules", {
    toy <- make_toy_runs(m = 4)
    res <- pool_limma_imputations(toy$runs, "B_vs_A", observed = toy$observed)
    k <- which(rowSums(!toy$observed) > 0)[1]

    m  <- 4
    Q  <- vapply(toy$runs, function(r) r$fit2$coefficients[k, 1], numeric(1))
    U  <- vapply(toy$runs, function(r) r$fit2$s2.post[k] * r$fit2$stdev.unscaled[k, 1]^2,
                 numeric(1))
    dfc <- mean(vapply(toy$runs, function(r) r$fit2$df.total[k], numeric(1)))
    W  <- mean(U); B <- var(Q); Tt <- W + (1 + 1 / m) * B
    lam <- (1 + 1 / m) * B / Tt
    r   <- (1 + 1 / m) * B / W
    df_old <- (m - 1) / lam^2
    df_obs <- (dfc + 1) / (dfc + 3) * dfc * (1 - lam)
    df_br  <- df_old * df_obs / (df_old + df_obs)
    t_hand <- mean(Q) / sqrt(Tt)

    expect_true(res$pooled[k])
    expect_equal(res$pooled_logFC[k], mean(Q))
    expect_equal(res$W[k], W)
    expect_equal(res$B[k], B)
    expect_equal(res$T[k], Tt)
    expect_equal(res$t_pooled[k], t_hand)
    expect_equal(res$df_com[k], dfc)
    expect_equal(res$df_BR[k], df_br)
    expect_equal(res$riv[k], r)
    expect_equal(res$fmi[k], (r + 2 / (df_br + 3)) / (r + 1))
    expect_equal(res$z_pooled[k], limma::zscoreT(t_hand, df = df_br))
})

test_that("identical imputations reproduce the single-fit moderated t", {
    toy <- make_toy_runs(m = 3, identical_draws = TRUE)
    res <- pool_limma_imputations(toy$runs, "B_vs_A")
    fit <- toy$runs[[1]]$fit2
    expect_false(any(res$pooled))
    expect_equal(res$t_pooled, unname(fit$t[, 1]))
    expect_equal(res$df_BR, as.numeric(fit$df.total))
    expect_equal(res$z_pooled, limma::zscoreT(unname(fit$t[, 1]), df = fit$df.total))
    expect_true(all(res$fmi == 0))
})

test_that("the pooled t is not the mean of the per-draw t-statistics", {
    toy <- make_toy_runs(m = 5)
    res <- pool_limma_imputations(toy$runs, "B_vs_A")
    mean_t <- rowMeans(sapply(toy$runs, function(r) r$fit2$t[, 1]))
    imputed <- res$pooled
    expect_true(any(imputed))
    expect_false(isTRUE(all.equal(res$t_pooled[imputed], unname(mean_t[imputed]))))
    # The between-draw term only adds variance.
    expect_true(all(res$T >= res$W))
})

test_that("one imputation gives the single-fit statistic, unpooled", {
    toy <- make_toy_runs(m = 1)
    res <- pool_limma_imputations(toy$runs, "B_vs_A")
    expect_false(any(res$pooled))
    expect_equal(res$t_pooled, unname(toy$runs[[1]]$fit2$t[, 1]))
    expect_equal(unique(res$m), 1L)
})

test_that("fits with different protein rows are rejected", {
    toy <- make_toy_runs(m = 2)
    n <- nrow(toy$runs[[2]]$fit2)
    toy$runs[[2]]$fit2 <- toy$runs[[2]]$fit2[rev(seq_len(n)), ]
    expect_error(pool_limma_imputations(toy$runs, "B_vs_A"), "different protein rows")
})

test_that("DE methods without a limma fit are rejected with a clear message", {
    toy <- make_toy_runs(m = 2)
    toy$runs[[1]]$fit2 <- NULL
    expect_error(pool_limma_imputations(toy$runs, "B_vs_A"), "de.method 'limma'")
})

test_that("fdrtool's p-value replacement does not change the pooled z", {
    skip_if_not_installed("fdrtool")
    toy <- make_toy_runs(m = 3)
    before <- pool_limma_imputations(toy$runs, "B_vs_A")
    # As run_limma_proteomics() applies it: only fit2$p.value is rewritten.
    corrected <- lapply(toy$runs, function(r) {
        r$fit2$p.value[, 1] <- fdrtool::fdrtool(r$fit2$t[, 1], statistic = "normal",
                                                plot = FALSE, verbose = FALSE)$pval
        r
    })
    after <- pool_limma_imputations(corrected, "B_vs_A")
    expect_identical(after$z_pooled, before$z_pooled)
})

test_that("observed counts follow the contrast sides of the design", {
    toy <- make_toy_runs(m = 2)
    res <- pool_limma_imputations(toy$runs, "B_vs_A", observed = toy$observed)
    expect_equal(res$n_observed_denominator, unname(rowSums(toy$observed[, 1:3])))
    expect_equal(res$n_observed_numerator, unname(rowSums(toy$observed[, 4:6])))
    expect_equal(res$frac_imputed, unname(1 - rowSums(toy$observed) / 6))
    # Without a mask the counts are NA rather than invented.
    expect_true(all(is.na(pool_limma_imputations(toy$runs, "B_vs_A")$n_observed)))
})

test_that("pooling leaves the DE summary unchanged", {
    toy <- make_toy_runs(m = 3)
    cfg <- list(modes = list(proteomics = list(
        de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5, use_adj_for_pass1 = TRUE),
        imputation = list(multi_imputation = TRUE, no_repetitions = 3, min_no_passed = 1),
        de_table = list(id_col = "FeatureID"),
        effects = list(samples = "SampleID")
    )))
    runs_de_tables <- lapply(toy$runs, function(r) r$de_tables)
    before <- summarize_limma_mult_imputation(runs_de_tables, cfg)
    invisible(pool_limma_imputations(toy$runs, "B_vs_A", observed = toy$observed))
    after <- summarize_limma_mult_imputation(lapply(toy$runs, function(r) r$de_tables), cfg)
    expect_identical(after, before)
})

test_that("the collapse keeps the legacy padj rule by default and most-observed on request", {
    cfg <- list(modes = list(proteomics = list(de_table = list(id_col = "FeatureID"))))
    ann <- data.frame(FeatureID = c("p1", "p2", "p3"), Genes = c("GA", "GA;GX", "GB"),
                      stringsAsFactors = FALSE)
    tbl <- data.frame(FeatureID = c("p1", "p2", "p3"), padj = c(0.5, 0.01, 0.2),
                      n_observed = c(6L, 3L, 5L), AveExpr = c(18, 22, 20),
                      stringsAsFactors = FALSE)

    legacy <- map_proteins_to_gene_symbols(tbl, ann, cfg)
    expect_equal(legacy$ProteinID[legacy$FeatureID == "GA"], "p2")

    obs <- map_proteins_to_gene_symbols(tbl, ann, cfg, priority = "observed")
    expect_equal(obs$ProteinID[obs$FeatureID == "GA"], "p1")

    # Ties on observed counts fall to the higher average expression.
    tbl$n_observed <- c(6L, 6L, 5L)
    tie <- map_proteins_to_gene_symbols(tbl, ann, cfg, priority = "observed")
    expect_equal(tie$ProteinID[tie$FeatureID == "GA"], "p2")
})

test_that("the observed collapse prefers a protein with a finite statistic", {
    cfg <- list(modes = list(proteomics = list(de_table = list(id_col = "FeatureID"))))
    ann <- data.frame(FeatureID = c("p1", "p2", "p3"), Genes = c("GA", "GA", "GB"),
                      stringsAsFactors = FALSE)
    tbl <- data.frame(FeatureID = c("p1", "p2", "p3"), stat = c(NA, 1.5, NaN),
                      n_observed = c(6L, 3L, 5L), AveExpr = c(18, 22, 20),
                      stringsAsFactors = FALSE)
    obs <- map_proteins_to_gene_symbols(tbl, ann, cfg, priority = "observed")
    # p1 has more observed values but no z; the gene keeps p2.
    expect_equal(obs$ProteinID[obs$FeatureID == "GA"], "p2")
    # A gene whose only protein has no z is still kept, as before.
    expect_equal(obs$ProteinID[obs$FeatureID == "GB"], "p3")
})

test_that("build_pooled_ranking_tables ranks on z_pooled and collapses per gene", {
    toy <- make_toy_runs(m = 3)
    cfg <- list(modes = list(proteomics = list(de_table = list(id_col = "FeatureID"))))
    pre <- list(expr_filt = ifelse(toy$observed, 1, NA_real_))
    de_res <- list(method = "limma", runs = toy$runs, summary_df = NULL)
    tabs <- build_pooled_ranking_tables(de_res, pre, cfg)
    expect_named(tabs, "B_vs_A")
    tb <- tabs$B_vs_A
    expect_identical(tb$stat, tb$z_pooled)
    expect_false(any(duplicated(tb$FeatureID)))
    expect_true(all(tb$FeatureID %in% paste0("G", 1:40)))
    expect_true(all(is.na(tb$legacy_stat)))
})

test_that("the pooled module skips unless limma ran and the option is on", {
    off <- list(modes = list(proteomics = list(pathway = list(enabled = TRUE))))
    expect_null(suppressMessages(mod_proteomics_pathway_pooled(list(method = "limma"), NULL, off, tempdir())))
    on <- list(modes = list(proteomics = list(pathway = list(
        enabled = TRUE, pooled_ranking = list(enabled = TRUE)))))
    expect_message(res <- mod_proteomics_pathway_pooled(list(method = "ttest"), NULL, on, tempdir()),
                   "needs de.method 'limma'")
    expect_null(res)
})
