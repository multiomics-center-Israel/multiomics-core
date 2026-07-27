# Wiring tests for the metabolomics TabPFN-vs-RF classifier benchmark POC
# (R/domain/metabolomics/04c_classifier_benchmark.R).
#
# The RF path runs with no Python. The TabPFN path is skipped unless a tabpfn
# venv is available, so CI stays green without the Python dependency.

# Build a tiny, seeded, 2-class synthetic `pre` satisfying assert_pre_contract.
make_synthetic_pre <- function(n_per_class = 12, n_feat = 30, n_signal = 8, seed = 1) {
    withr::with_seed(seed, {
        n <- 2 * n_per_class
        samples <- paste0("S", seq_len(n))
        feats   <- paste0("feat_", seq_len(n_feat))
        grp     <- rep(c("A", "B"), each = n_per_class)

        mat <- matrix(stats::rnorm(n_feat * n, sd = 1),
                      nrow = n_feat, dimnames = list(feats, samples))
        # Shift the first n_signal features in group B so a classifier has signal.
        mat[seq_len(n_signal), grp == "B"] <- mat[seq_len(n_signal), grp == "B"] + 2.5

        meta <- data.frame(sample_id = samples, sample_type = grp,
                           stringsAsFactors = FALSE)
        row_data <- data.frame(feature_id = feats, Name = feats,
                               stringsAsFactors = FALSE)
        list(expr_raw = mat, expr_filt = mat, expr_work = mat,
             meta = meta, row_data = row_data)
    })
}

make_config <- function(top_n = 10, k = 5, seed = 1, tabpfn_enabled = TRUE) {
    list(modes = list(metabolomics = list(
        effects   = list(samples = "sample_id", color = "sample_type"),
        benchmark = list(top_n = top_n, k = k, seed = seed,
                         tabpfn = list(enabled = tabpfn_enabled))
    )))
}

tabpfn_python <- function() {
    p <- Sys.getenv("TABPFN_PYTHON", "")
    if (nzchar(p)) return(p)
    if (.Platform$OS.type == "windows") "envs/tabpfn/Scripts/python.exe"
    else "envs/tabpfn/bin/python"
}


test_that("RF baseline path: benchmark returns a valid AUC without Python", {
    skip_if(
        !requireNamespace("ranger", quietly = TRUE) &&
            !requireNamespace("randomForest", quietly = TRUE),
        "no RF backend installed"
    )

    pre    <- make_synthetic_pre()
    # feature_sel_res = NULL exercises the univariate t-test feature fallback.
    config <- make_config(tabpfn_enabled = FALSE)

    bench <- suppressMessages(
        compute_metab_classifier_benchmark(pre, NULL, config)
    )

    expect_type(bench, "list")
    expect_false(is.null(bench$rf))
    expect_true(is.finite(bench$rf$auc))
    expect_gte(bench$rf$auc, 0)
    expect_lte(bench$rf$auc, 1)
    # Signal is present, so a working classifier should beat chance.
    expect_gt(bench$rf$auc, 0.5)
    expect_null(bench$tabpfn)  # disabled in config
    expect_equal(bench$n_features, 10L)
    expect_true(all(c("fpr", "tpr") %in% names(bench$rf$roc_data)))
})

test_that("benchmark is reproducible for a fixed seed", {
    skip_if(
        !requireNamespace("ranger", quietly = TRUE) &&
            !requireNamespace("randomForest", quietly = TRUE),
        "no RF backend installed"
    )
    pre    <- make_synthetic_pre()
    config <- make_config(tabpfn_enabled = FALSE)

    a <- suppressMessages(compute_metab_classifier_benchmark(pre, NULL, config))
    b <- suppressMessages(compute_metab_classifier_benchmark(pre, NULL, config))
    expect_equal(a$rf$auc, b$rf$auc)
})

test_that("format_classifier_benchmark: 2-row table, NA for skipped models", {
    skip_if(
        !requireNamespace("ranger", quietly = TRUE) &&
            !requireNamespace("randomForest", quietly = TRUE),
        "no RF backend installed"
    )
    pre    <- make_synthetic_pre()
    config <- make_config(tabpfn_enabled = FALSE)
    bench  <- suppressMessages(compute_metab_classifier_benchmark(pre, NULL, config))

    tbl <- format_classifier_benchmark(bench)
    expect_equal(nrow(tbl), 2L)
    expect_setequal(tbl$model, c("RF", "TabPFN"))
    expect_true(is.na(tbl$AUC[tbl$model == "TabPFN"]))  # TabPFN disabled -> NA
    expect_false(is.na(tbl$AUC[tbl$model == "RF"]))
})

test_that("non-2-class problems are skipped, not errored", {
    pre <- make_synthetic_pre()
    pre$meta$sample_type <- rep(c("A", "B", "C"), length.out = nrow(pre$meta))
    config <- make_config()

    expect_null(suppressMessages(
        compute_metab_classifier_benchmark(pre, NULL, config)
    ))
})

test_that("uses the RF importance ranking when feature selection is supplied", {
    skip_if(
        !requireNamespace("ranger", quietly = TRUE) &&
            !requireNamespace("randomForest", quietly = TRUE),
        "no RF backend installed"
    )
    pre <- make_synthetic_pre()
    # Rank the signal features (feat_1..feat_8) at the top.
    imp <- data.frame(
        feature_id = paste0("feat_", c(1:8, 9:30)),
        importance = c(seq(30, 23), seq(22, 1)),
        stringsAsFactors = FALSE
    )
    fs <- list(rf = list(importance_df = imp), plsda = NULL)
    config <- make_config(top_n = 8, tabpfn_enabled = FALSE)

    bench <- suppressMessages(compute_metab_classifier_benchmark(pre, fs, config))
    expect_equal(bench$feature_source, "rf_importance")
    expect_equal(bench$n_features, 8L)
    expect_gt(bench$rf$auc, 0.5)
})

test_that("TabPFN path runs end-to-end when a venv is available", {
    py <- tabpfn_python()
    skip_if(!file.exists(py), "no tabpfn venv (set TABPFN_PYTHON / run setup-tabpfn-venv.sh)")
    skip_if(
        !requireNamespace("ranger", quietly = TRUE) &&
            !requireNamespace("randomForest", quietly = TRUE),
        "no RF backend installed"
    )

    pre    <- make_synthetic_pre()
    config <- make_config(tabpfn_enabled = TRUE)
    bench  <- suppressMessages(compute_metab_classifier_benchmark(pre, NULL, config))

    expect_false(is.null(bench$tabpfn))
    expect_true(is.finite(bench$tabpfn$auc))
    expect_gte(bench$tabpfn$auc, 0)
    expect_lte(bench$tabpfn$auc, 1)
})
