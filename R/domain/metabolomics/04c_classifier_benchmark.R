# R/domain/metabolomics/04c_classifier_benchmark.R
#
# POC harness: does a tabular foundation model (TabPFN) beat a Random Forest
# baseline at classifying samples from metabolomics abundances, after feature
# selection? Both models are scored under ONE shared, seeded k-fold CV so the two
# AUCs are directly comparable.
#
# Scope: a standalone benchmark — NOT wired into the {targets} DAG, reports, or
# the Shiny payload. It reuses the existing feature-selection output
# (metab_feature_sel_res) and the pre-processing contract (expr_work, meta).
#
# The TabPFN side runs as an ISOLATED Python subprocess via processx (never
# imported into R), reusing the pinned-venv pattern established by mummichog
# (see 06c_mummichog_pinned.R). The main R pipeline keeps zero Python deps, and
# a missing/disabled TabPFN venv degrades gracefully to the RF-only result.
#
# Data safety: the TabPFN model runs LOCALLY only. We never call the hosted
# TabPFN API/client — that would send abundances off-machine.
#
# Reuses: assert_pre_contract, assert_numeric_matrix, %||%, and the generic
# venv-path helpers .mmc_exists_nofollow / .mmc_abs_keep_symlink (defined in
# 06c; resolved at call time, so source order within the folder is irrelevant).


# ==== venv resolution =======================================================

#' Platform-appropriate default path to the pinned TabPFN venv's Python
#'
#' Mirrors `.mmc_default_python()`: `python -m venv` puts the interpreter under
#' `bin/` on Unix and `Scripts/` on Windows.
#' @return Relative path to the venv Python for the current OS.
#' @noRd
.clf_default_python <- function() {
    if (.Platform$OS.type == "windows") {
        "envs/tabpfn/Scripts/python.exe"
    } else {
        "envs/tabpfn/bin/python"
    }
}

#' Is a usable TabPFN interpreter available?
#'
#' Soft preflight (returns a logical rather than aborting) so the benchmark can
#' fall back to RF-only when the venv is absent or `import tabpfn` fails. Uses
#' the same non-canonicalising path handling as mummichog so a venv symlink is
#' invoked as itself.
#'
#' @param python Path to the candidate interpreter.
#' @return TRUE if the interpreter exists and imports tabpfn; FALSE otherwise.
#' @noRd
.clf_tabpfn_available <- function(python) {
    if (!nzchar(python) || !.mmc_exists_nofollow(python)) return(FALSE)
    python <- .mmc_abs_keep_symlink(python)
    probe <- tryCatch(
        processx::run(
            python, c("-c", "import tabpfn"),
            error_on_status = FALSE, timeout = 120
        ),
        error = function(e) NULL
    )
    !is.null(probe) && isTRUE(probe$status == 0)
}


# ==== data preparation ======================================================

#' Resolve the 2-class response and feature matrix for the benchmark
#'
#' Mirrors `run_metabolomics_rf()`'s label resolution and matrix prep exactly so
#' the benchmark trains on the same design the feature-selection step did.
#'
#' @param pre    List from preprocess_metabolomics() (pre contract).
#' @param config Full pipeline config.
#' @return list(X = samples x features matrix, resp = 0/1 integer,
#'   condition = factor, levels) or NULL when there are not exactly 2 classes.
#' @noRd
.clf_prep_matrix <- function(pre, config) {
    cfg <- config$modes$metabolomics
    bench_cfg <- cfg$benchmark %||% list()

    condition_col <- bench_cfg$condition_column %||%
                     cfg$rf$condition_column %||%
                     cfg$de$condition_column %||%
                     cfg$effects$color %||% "sample_type"
    sample_col <- cfg$effects$samples %||% "sample_id"

    mat  <- pre$expr_work
    meta <- pre$meta
    assert_numeric_matrix(mat, "metab_expr_work")

    meta <- meta[match(colnames(mat), meta[[sample_col]]), , drop = FALSE]
    condition <- factor(meta[[condition_col]])

    lvls <- levels(condition)
    if (length(lvls) != 2) {
        message("classifier benchmark: need exactly 2 levels in '",
                condition_col, "' — found ", length(lvls), "; skipping")
        return(NULL)
    }

    resp <- as.integer(condition == lvls[2])

    # samples x features, median-impute NAs, drop zero-variance columns
    X <- t(mat)
    for (j in seq_len(ncol(X))) {
        nas <- is.na(X[, j])
        if (any(nas)) X[nas, j] <- stats::median(X[!nas, j], na.rm = TRUE)
    }
    col_vars <- apply(X, 2, stats::var, na.rm = TRUE)
    X <- X[, col_vars > 0, drop = FALSE]

    list(X = X, resp = resp, condition = condition, levels = lvls)
}

#' Pick the top-N feature columns to feed the classifiers
#'
#' Foundation models such as TabPFN 2.5 are happiest at <=~50 features, so
#' selection-first is essential. Preference order: reuse the RF importance
#' ranking from feature selection, else the PLS-DA VIP ranking, else a quick
#' univariate Welch t-test fallback (so the benchmark still runs when neither
#' `rf:` nor `plsda:` is configured — the metabolomics template ships neither).
#'
#' @param X               samples x features matrix (from .clf_prep_matrix).
#' @param resp            0/1 integer response.
#' @param feature_sel_res Result of mod_metabolomics_feature_selection() or NULL.
#' @param top_n           Number of features to keep.
#' @return list(X = reduced matrix, source = which ranking was used).
#' @noRd
.clf_select_features <- function(X, resp, feature_sel_res, top_n) {
    ranked <- NULL
    src <- NULL

    imp_df <- feature_sel_res$rf$importance_df
    vip_df <- feature_sel_res$plsda$vip_df
    if (!is.null(imp_df) && "feature_id" %in% names(imp_df)) {
        ranked <- as.character(imp_df$feature_id)
        src <- "rf_importance"
    } else if (!is.null(vip_df) && "feature_id" %in% names(vip_df)) {
        ranked <- as.character(vip_df$feature_id)
        src <- "plsda_vip"
    }

    if (is.null(ranked)) {
        # Univariate Welch t-test ranking (ascending p-value).
        pvals <- apply(X, 2, function(col) {
            tryCatch(
                stats::t.test(col[resp == 1], col[resp == 0])$p.value,
                error = function(e) NA_real_
            )
        })
        ranked <- names(sort(pvals, na.last = NA))
        src <- "univariate_ttest"
    }

    ranked <- ranked[ranked %in% colnames(X)]
    keep <- utils::head(ranked, top_n)
    if (length(keep) < 2) {
        # Ranking did not line up with the matrix — fall back to all features.
        keep <- colnames(X)
        src <- paste0(src, "+fallback_all")
    }

    list(X = X[, keep, drop = FALSE], source = src)
}

#' Assign stratified k-fold CV labels
#'
#' Each class is spread evenly across folds so no fold trains on a single class.
#' k is capped at the smallest class size; for tiny cohorts this approaches
#' leave-one-out. Seeded via withr for reproducibility.
#'
#' @param resp 0/1 integer response.
#' @param k    Requested number of folds.
#' @param seed RNG seed.
#' @return Integer fold vector (1..k_eff) aligned to `resp`, or NULL if a class
#'   has fewer than 2 members (CV not possible).
#' @noRd
.clf_make_folds <- function(resp, k, seed) {
    min_class <- min(table(resp))
    if (min_class < 2) {
        message("classifier benchmark: smallest class has ", min_class,
                " sample(s) — cannot cross-validate; skipping")
        return(NULL)
    }
    k_eff <- max(2L, min(as.integer(k), min_class))

    withr::with_seed(seed, {
        folds <- integer(length(resp))
        for (cls in unique(resp)) {
            idx <- which(resp == cls)
            folds[idx] <- sample(rep(seq_len(k_eff), length.out = length(idx)))
        }
        folds
    })
}


# ==== ROC helper ============================================================

#' Build a ROC data.frame + AUC from held-out predictions
#'
#' Deliberately does NOT auto-flip to AUC >= 0.5: a model doing worse than
#' chance is a real result the benchmark should surface. Both classifiers use
#' the same response encoding and the same `direction`, so the comparison is
#' fair. Falls back to a manual sweep when pROC is unavailable.
#'
#' @param resp     0/1 integer response (held-out samples only).
#' @param prob_pos Predicted P(class = 1) for the same samples.
#' @return list(roc_data = data.frame(fpr, tpr), auc) or NULL on failure.
#' @noRd
.clf_roc_from_predictions <- function(resp, prob_pos) {
    if (length(unique(resp)) < 2) return(NULL)

    if (requireNamespace("pROC", quietly = TRUE)) {
        roc_obj <- tryCatch(
            pROC::roc(resp, prob_pos, quiet = TRUE, levels = c(0, 1),
                      direction = "<"),
            error = function(e) NULL
        )
        if (is.null(roc_obj)) return(NULL)
        return(list(
            roc_data = data.frame(
                fpr = 1 - roc_obj$specificities,
                tpr = roc_obj$sensitivities
            ),
            auc = as.numeric(pROC::auc(roc_obj))
        ))
    }

    # Manual ROC + trapezoidal AUC.
    thr <- sort(unique(c(-Inf, prob_pos, Inf)), decreasing = TRUE)
    n_pos <- sum(resp == 1); n_neg <- sum(resp == 0)
    tpr <- vapply(thr, function(t) sum(prob_pos >= t & resp == 1) / max(n_pos, 1), numeric(1))
    fpr <- vapply(thr, function(t) sum(prob_pos >= t & resp == 0) / max(n_neg, 1), numeric(1))
    ord <- order(fpr, tpr)
    fpr <- fpr[ord]; tpr <- tpr[ord]
    auc <- sum(diff(fpr) * (utils::head(tpr, -1) + utils::tail(tpr, -1)) / 2)
    list(roc_data = data.frame(fpr = fpr, tpr = tpr), auc = auc)
}


# ==== Random Forest CV baseline =============================================

#' Cross-validated Random Forest AUC
#'
#' Trains a probability forest on each CV training split and pools the held-out
#' positive-class probabilities into a single ROC. Uses explicit CV (not OOB) so
#' the score is computed exactly like the TabPFN side. Prefers ranger, falls
#' back to randomForest.
#'
#' @param X       samples x features matrix.
#' @param resp    0/1 integer response.
#' @param folds   Integer fold vector aligned to rows of X.
#' @param seed    RNG seed.
#' @param n_trees Number of trees per forest.
#' @return list(roc_data, auc, method) or NULL when no RF backend / ROC fails.
compute_rf_cv_auc <- function(X, resp, folds, seed = 1234, n_trees = 500) {
    yf <- factor(resp, levels = c(0, 1))
    oof <- rep(NA_real_, length(resp))

    use_ranger <- requireNamespace("ranger", quietly = TRUE)
    use_rf     <- requireNamespace("randomForest", quietly = TRUE)
    if (!use_ranger && !use_rf) {
        message("classifier benchmark (RF): no RF package available — skipping")
        return(NULL)
    }

    withr::with_seed(seed, {
        for (f in sort(unique(folds))) {
            tr <- folds != f
            te <- folds == f
            if (length(unique(resp[tr])) < 2) next  # single-class train split

            if (use_ranger) {
                fit <- tryCatch(
                    ranger::ranger(x = X[tr, , drop = FALSE], y = yf[tr],
                                   num.trees = n_trees, probability = TRUE,
                                   seed = seed),
                    error = function(e) NULL
                )
                if (is.null(fit)) next
                pred <- stats::predict(fit, data = X[te, , drop = FALSE])$predictions
                oof[te] <- pred[, "1"]
            } else {
                Xtr <- X[tr, , drop = FALSE]; Xte <- X[te, , drop = FALSE]
                safe <- paste0("F", seq_len(ncol(Xtr)))
                colnames(Xtr) <- safe; colnames(Xte) <- safe
                fit <- tryCatch(
                    randomForest::randomForest(x = Xtr, y = yf[tr],
                                               ntree = n_trees),
                    error = function(e) NULL
                )
                if (is.null(fit)) next
                oof[te] <- stats::predict(fit, newdata = Xte, type = "prob")[, "1"]
            }
        }
    })

    mask <- !is.na(oof)
    if (sum(mask) < 2) {
        message("classifier benchmark (RF): too few held-out predictions — skipping")
        return(NULL)
    }
    roc <- .clf_roc_from_predictions(resp[mask], oof[mask])
    if (is.null(roc)) return(NULL)

    list(roc_data = roc$roc_data, auc = roc$auc,
         method = if (use_ranger) "ranger" else "randomForest")
}


# ==== TabPFN CV via isolated subprocess =====================================

#' Cross-validated TabPFN AUC (isolated Python subprocess)
#'
#' Writes the feature matrix and the (label, fold) assignments to TSVs, runs
#' `scripts/tabpfn_cv.py` in the pinned venv, and reads back the pooled-held-out
#' AUC + ROC points. Uses the SAME `folds` as the RF baseline for a fair
#' comparison. Degrades gracefully to NULL (RF-only benchmark) when disabled in
#' config or when the venv / tabpfn import is unavailable — unlike mummichog,
#' which hard-stops, because this is an optional benchmark, not a pipeline stage.
#'
#' @param X      samples x features matrix.
#' @param resp   0/1 integer response.
#' @param folds  Integer fold vector aligned to rows of X.
#' @param config Full pipeline config.
#' @param seed   RNG seed passed through to numpy/torch.
#' @return list(roc_data, auc, method) or NULL when skipped/unavailable.
run_tabpfn_cv_auc <- function(X, resp, folds, config, seed = 1234) {
    bench_cfg  <- config$modes$metabolomics$benchmark %||% list()
    tabpfn_cfg <- bench_cfg$tabpfn %||% list()
    if (identical(tabpfn_cfg$enabled, FALSE)) {
        message("classifier benchmark (TabPFN): disabled in config — skipping")
        return(NULL)
    }

    python <- Sys.getenv("TABPFN_PYTHON", .clf_default_python())
    if (!.clf_tabpfn_available(python)) {
        message("classifier benchmark (TabPFN): interpreter '", python,
                "' not found or cannot import tabpfn — running RF only. ",
                "Build the venv with `bash setup-tabpfn-venv.sh` and set TABPFN_PYTHON.")
        return(NULL)
    }
    python <- .mmc_abs_keep_symlink(python)

    script <- file.path(getwd(), "scripts", "tabpfn_cv.py")
    if (!file.exists(script)) {
        message("classifier benchmark (TabPFN): scripts/tabpfn_cv.py not found — skipping")
        return(NULL)
    }

    work <- file.path(tempdir(), paste0("tabpfn_cv_", as.integer(seed)))
    dir.create(work, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)

    x_file   <- file.path(work, "X.tsv")
    lab_file <- file.path(work, "labels.tsv")
    out_file <- file.path(work, "result.json")

    sample_ids <- rownames(X) %||% as.character(seq_len(nrow(X)))
    utils::write.table(
        data.frame(sample_id = sample_ids, X, check.names = FALSE),
        x_file, sep = "\t", quote = FALSE, row.names = FALSE
    )
    utils::write.table(
        data.frame(sample_id = sample_ids, label = resp, fold = folds),
        lab_file, sep = "\t", quote = FALSE, row.names = FALSE
    )

    res <- tryCatch(
        processx::run(
            command = python,
            args = c(script, "--x", x_file, "--labels", lab_file,
                     "--out", out_file, "--seed", as.character(seed)),
            env = c("current", MPLBACKEND = "Agg"),
            error_on_status = FALSE, timeout = 3600
        ),
        error = function(e) NULL
    )
    if (is.null(res) || !file.exists(out_file)) {
        message("classifier benchmark (TabPFN): subprocess produced no output — skipping")
        return(NULL)
    }

    parsed <- tryCatch(jsonlite::fromJSON(out_file), error = function(e) NULL)
    if (is.null(parsed) || !is.null(parsed$error) || is.null(parsed$auc)) {
        msg <- if (!is.null(parsed$error)) parsed$error else "no AUC returned"
        message("classifier benchmark (TabPFN): ", msg, " — skipping")
        return(NULL)
    }

    list(
        roc_data = data.frame(fpr = as.numeric(parsed$fpr),
                              tpr = as.numeric(parsed$tpr)),
        auc = as.numeric(parsed$auc),
        method = "tabpfn"
    )
}


# ==== orchestrator ==========================================================

#' Benchmark TabPFN against a Random Forest baseline on metabolomics data
#'
#' Reuses the existing feature-selection ranking to reduce the feature space,
#' then scores both a Random Forest and TabPFN under one shared, seeded k-fold
#' cross-validation and reports their pooled-held-out AUCs side by side. TabPFN
#' is optional: when its venv is absent or disabled the RF baseline still runs.
#'
#' @param pre             List from preprocess_metabolomics() (pre contract).
#' @param feature_sel_res Result of mod_metabolomics_feature_selection(), or NULL
#'   to fall back to a univariate feature ranking.
#' @param config          Full pipeline config. Reads
#'   `modes$metabolomics$benchmark`: `top_n` (default `rf$top_n` or 20),
#'   `k` (default 5), `seed` (default `rf$seed` or 1234), and `tabpfn$enabled`.
#' @return list(rf, tabpfn, n_samples, n_features, feature_source, k, seed,
#'   levels) — rf/tabpfn are each list(roc_data, auc, method) or NULL — or NULL
#'   when the data cannot be benchmarked (not 2-class, or CV impossible).
#' @examples
#' \dontrun{
#' pre <- targets::tar_read(metab_pre)
#' fs  <- targets::tar_read(metab_feature_sel_res)
#' bench <- compute_metab_classifier_benchmark(pre, fs, config)
#' format_classifier_benchmark(bench)
#' }
compute_metab_classifier_benchmark <- function(pre, feature_sel_res, config) {
    assert_pre_contract(pre, stage = "metabolomics")

    bench_cfg <- config$modes$metabolomics$benchmark %||% list()
    top_n <- bench_cfg$top_n %||% config$modes$metabolomics$rf$top_n %||% 20
    k     <- bench_cfg$k    %||% 5
    seed  <- bench_cfg$seed %||% config$modes$metabolomics$rf$seed %||% 1234
    n_trees <- config$modes$metabolomics$rf$n_trees %||% 500

    prep <- .clf_prep_matrix(pre, config)
    if (is.null(prep)) return(NULL)

    sel <- .clf_select_features(prep$X, prep$resp, feature_sel_res, top_n)
    X <- sel$X

    folds <- .clf_make_folds(prep$resp, k, seed)
    if (is.null(folds)) return(NULL)
    k_eff <- length(unique(folds))

    message("classifier benchmark: ", nrow(X), " samples x ", ncol(X),
            " features (", sel$source, "), ", k_eff, "-fold CV, seed ", seed)

    rf_res     <- compute_rf_cv_auc(X, prep$resp, folds, seed = seed, n_trees = n_trees)
    tabpfn_res <- run_tabpfn_cv_auc(X, prep$resp, folds, config, seed = seed)

    list(
        rf             = rf_res,
        tabpfn         = tabpfn_res,
        n_samples      = nrow(X),
        n_features     = ncol(X),
        feature_source = sel$source,
        k              = k_eff,
        seed           = seed,
        levels         = prep$levels
    )
}


#' Format a classifier benchmark as a compact AUC table
#'
#' @param bench Result of compute_metab_classifier_benchmark().
#' @return A data.frame with one row per model (RF, TabPFN) — models that were
#'   skipped report `NA` AUC — or NULL when the benchmark itself was skipped.
format_classifier_benchmark <- function(bench) {
    if (is.null(bench)) return(NULL)
    row <- function(model, res) {
        data.frame(
            model      = model,
            n_samples  = bench$n_samples,
            n_features = bench$n_features,
            k          = bench$k,
            AUC        = if (is.null(res)) NA_real_ else round(res$auc, 4),
            stringsAsFactors = FALSE
        )
    }
    rbind(row("RF", bench$rf), row("TabPFN", bench$tabpfn))
}
