# Tests for the self-contained limma rotation set test
# (run_metabolomics_selfcontained) — Fix 1.
#
# The point of a self-contained test is that it detects a set that changed
# *at all*, so it stays powered when a large fraction of features move — the
# regime where the competitive methods (ORA/GSEA/ssGSEA) return nothing.

# ---- shared synthetic fixture -----------------------------------------------

# Build a small metabolomics `pre` + `config` with one coordinately-moving
# pathway (MOVER) and two flat pathways (FLAT1/FLAT2). Data generation is
# seeded; fry() itself is deterministic.
.make_selfcontained_fixture <- function(enabled = TRUE, method = "fry") {
    n_per_grp <- 5L
    groups  <- rep(c("B", "A"), each = n_per_grp)   # B = numerator
    samples <- paste0("s", seq_along(groups))

    n_feat   <- 30L
    kegg_ids <- sprintf("C%05d", 10001:(10000 + n_feat))
    feat_ids <- paste0("F", seq_len(n_feat))

    mat <- withr::with_seed(1, matrix(
        stats::rnorm(n_feat * length(groups), mean = 10, sd = 0.3),
        nrow = n_feat, dimnames = list(feat_ids, samples)
    ))
    # Movers: first 6 features strongly UP in group B (coordinated change).
    movers <- 1:6
    mat[movers, groups == "B"] <- mat[movers, groups == "B"] + 3

    meta <- data.frame(sample_id = samples, sample_type = groups,
                       stringsAsFactors = FALSE)
    row_data <- data.frame(feature_id = feat_ids, KEGG = kegg_ids,
                           stringsAsFactors = FALSE)
    rownames(row_data) <- feat_ids

    pre <- list(
        expr_work = mat, expr_log = mat, expr_filt = mat,
        meta = meta, row_data = row_data,
        info = list(normalization = list(scaling = "none"))
    )

    gmt <- tempfile(fileext = ".gmt")
    writeLines(c(
        paste(c("MOVER", "Moving pathway", kegg_ids[movers]),  collapse = "\t"),
        paste(c("FLAT1", "Flat pathway 1", kegg_ids[7:12]),    collapse = "\t"),
        paste(c("FLAT2", "Flat pathway 2", kegg_ids[13:18]),   collapse = "\t")
    ), gmt)

    config <- list(modes = list(metabolomics = list(
        de      = list(condition_column = "sample_type"),
        effects = list(samples = "sample_id", color = "sample_type"),
        enrichment = list(
            run_enrichment = TRUE,
            gmt_file       = gmt,
            mapping_file   = NULL,
            selfcontained  = list(enabled = enabled, method = method)
        )
    )))

    list(pre = pre, config = config, gmt = gmt, movers = movers)
}


test_that("run_metabolomics_selfcontained recovers a coordinately-moving pathway", {
    skip_if_not_installed("limma")
    skip_if_not_installed("withr")

    fx <- .make_selfcontained_fixture()
    on.exit(unlink(fx$gmt), add = TRUE)

    res <- run_metabolomics_selfcontained(fx$pre, fx$config, "B - A",
                                          contrast_label = "B_vs_A")

    expect_false(is.null(res))
    expect_true(is.data.frame(res$table))
    expect_identical(res$method, "limma_fry")
    expect_identical(res$contrast, "B_vs_A")

    # Schema slots into the shared enrichment result shape.
    expect_true(all(c("pathway", "n_hits", "direction", "PValue", "FDR",
                      "PValue_mixed", "FDR_mixed", "library") %in%
                    colnames(res$table)))

    # The moving pathway is recovered, up-regulated, with all 6 members matched,
    # and is more significant than the flat pathways.
    mover <- res$table[grepl("^MOVER", res$table$pathway), , drop = FALSE]
    expect_equal(nrow(mover), 1L)
    expect_lt(mover$FDR, 0.05)
    expect_identical(mover$direction, "Up")
    expect_equal(mover$n_hits, length(fx$movers))

    flat <- res$table[grepl("^FLAT", res$table$pathway), , drop = FALSE]
    expect_true(all(flat$FDR > mover$FDR))
})


test_that("run_metabolomics_selfcontained is deterministic (fry, no seed)", {
    skip_if_not_installed("limma")
    skip_if_not_installed("withr")

    fx <- .make_selfcontained_fixture()
    on.exit(unlink(fx$gmt), add = TRUE)

    r1 <- run_metabolomics_selfcontained(fx$pre, fx$config, "B - A")
    r2 <- run_metabolomics_selfcontained(fx$pre, fx$config, "B - A")
    expect_equal(r1$table, r2$table)
})


test_that("run_metabolomics_selfcontained is opt-in and skips when disabled", {
    skip_if_not_installed("limma")
    skip_if_not_installed("withr")

    # selfcontained.enabled = FALSE -> NULL (no-op)
    fx_off <- .make_selfcontained_fixture(enabled = FALSE)
    on.exit(unlink(fx_off$gmt), add = TRUE)
    expect_null(run_metabolomics_selfcontained(fx_off$pre, fx_off$config, "B - A"))

    # run_enrichment = FALSE -> NULL even when selfcontained.enabled = TRUE
    fx_on <- .make_selfcontained_fixture(enabled = TRUE)
    on.exit(unlink(fx_on$gmt), add = TRUE)
    fx_on$config$modes$metabolomics$enrichment$run_enrichment <- FALSE
    expect_null(run_metabolomics_selfcontained(fx_on$pre, fx_on$config, "B - A"))
})


test_that("run_metabolomics_selfcontained falls back to fry on an unknown method", {
    skip_if_not_installed("limma")
    skip_if_not_installed("withr")

    fx <- .make_selfcontained_fixture(method = "bogus")
    on.exit(unlink(fx$gmt), add = TRUE)

    expect_warning(
        res <- run_metabolomics_selfcontained(fx$pre, fx$config, "B - A"),
        "unknown method"
    )
    expect_identical(res$method, "limma_fry")
})


test_that("run_metabolomics_selfcontained returns NULL for an unbuildable contrast", {
    skip_if_not_installed("limma")
    skip_if_not_installed("withr")

    fx <- .make_selfcontained_fixture()
    on.exit(unlink(fx$gmt), add = TRUE)

    # "C" is not a level of the condition column -> makeContrasts fails -> NULL.
    expect_warning(
        res <- run_metabolomics_selfcontained(fx$pre, fx$config, "C - A"),
        "cannot build contrast"
    )
    expect_null(res)
})
