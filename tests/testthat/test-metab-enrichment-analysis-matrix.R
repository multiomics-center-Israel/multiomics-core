# Tests for the enrichment analysis-matrix fix: QEA (globaltest) and ssGSEA now
# start from the SAME upstream analysis matrix, scale, and biological-sample set
# as DE (via .metab_de_matrix_condition()), instead of pre$expr_raw. Also covers
# the QEA `standardize` option and its wiring into globaltest::gt().
#
# Covers both QEA and ssGSEA, hence the file name. Each test carries its own
# skip_if_not_installed() so a missing globaltest does not skip the ssGSEA tests
# (or vice versa), and the pure-R helper tests always run.

# ---- shared fixture ----------------------------------------------------------

# One coordinately-moving pathway (MOVER) + two flat pathways (FLAT1/FLAT2).
# expr_raw is deliberately DIFFERENT from expr_work/expr_log when expr_raw_noise
# is set, so a test can prove the fixed code never reads expr_raw. `hetero`
# inflates per-feature scale heterogeneity so globaltest's `standardize` bites.
.make_enrichment_fixture <- function(with_qc = FALSE, standardize = FALSE,
                                     shuffle_meta = FALSE, expr_raw_noise = FALSE,
                                     hetero = FALSE) {
    groups  <- rep(c("B", "A"), each = 5L)   # B = numerator
    samples <- paste0("s", seq_along(groups))
    n_feat   <- 30L
    kegg_ids <- sprintf("C%05d", 10001:(10000 + n_feat))
    feat_ids <- paste0("F", seq_len(n_feat))

    mat <- withr::with_seed(7, matrix(
        stats::rnorm(n_feat * length(groups), mean = 10, sd = 0.3),
        nrow = n_feat, dimnames = list(feat_ids, samples)
    ))
    movers <- 1:6
    mat[movers, groups == "B"] <- mat[movers, groups == "B"] + 3
    if (hetero) {
        # Feature scales span ~1 to ~1e4 so standardize=TRUE (equal weighting)
        # produces a materially different globaltest result than FALSE.
        mat <- mat * 10^seq(0, 4, length.out = n_feat)
    }

    meta <- data.frame(sample_id = samples, sample_type = groups,
                       stringsAsFactors = FALSE)

    if (with_qc) {
        # QC sample with extreme values: if biological filtering did NOT drop it,
        # it would grossly change results (and add a 3rd group level).
        qc <- matrix(1e6, nrow = n_feat, ncol = 1,
                     dimnames = list(feat_ids, "QC_1"))
        mat  <- cbind(mat, qc)
        meta <- rbind(meta, data.frame(sample_id = "QC_1", sample_type = "QC"))
    }
    if (shuffle_meta) {
        # Prove alignment is by sample id, not row position.
        meta <- meta[withr::with_seed(3, sample(nrow(meta))), , drop = FALSE]
    }

    row_data <- data.frame(feature_id = feat_ids, KEGG = kegg_ids,
                           stringsAsFactors = FALSE)
    rownames(row_data) <- feat_ids

    expr_raw <- if (expr_raw_noise) {
        withr::with_seed(99, matrix(stats::rnorm(length(mat), 1e3, 500),
                                    nrow = nrow(mat), dimnames = dimnames(mat)))
    } else {
        mat
    }

    pre <- list(expr_raw = expr_raw, expr_work = mat, expr_log = mat,
                expr_filt = mat, meta = meta, row_data = row_data,
                info = list(normalization = list(scaling = "none")))

    gmt <- tempfile(fileext = ".gmt")
    writeLines(c(
        paste(c("MOVER", "Moving pathway", kegg_ids[movers]),  collapse = "\t"),
        paste(c("FLAT1", "Flat pathway 1", kegg_ids[7:12]),    collapse = "\t"),
        paste(c("FLAT2", "Flat pathway 2", kegg_ids[13:18]),   collapse = "\t")
    ), gmt)

    config <- list(modes = list(metabolomics = list(
        de      = list(condition_column = "sample_type"),
        effects = list(samples = "sample_id", color = "sample_type"),
        enrichment = list(run_enrichment = TRUE, gmt_file = gmt,
                          mapping_file = NULL,
                          qea = list(standardize = standardize))
    )))
    list(pre = pre, config = config, gmt = gmt)
}


# ---- .metab_de_matrix_condition() (pure R, always runs) ---------------------

test_that(".metab_de_matrix_condition selects expr_work vs expr_log by scaling", {
    work <- matrix(1:12, nrow = 3, dimnames = list(c("f1","f2","f3"),
                                                    c("s1","s2","s3","s4")))
    logm <- work * 10  # distinct so we can tell which was chosen
    meta <- data.frame(sample_id = c("s1","s2","s3","s4"),
                       sample_type = c("A","A","B","B"),
                       stringsAsFactors = FALSE)
    mk <- function(scaling) list(
        expr_work = work, expr_log = logm, meta = meta,
        info = list(normalization = list(scaling = scaling))
    )
    config <- list(modes = list(metabolomics = list(
        de = list(condition_column = "sample_type"),
        effects = list(samples = "sample_id", color = "sample_type")
    )))

    # NOTE: the expected branch set mirrors run_metabolomics_de() in
    # 03_differential.R (the parity contract). Keep both in sync.
    for (sc in c("auto", "pareto", "range")) {
        expect_equal(.metab_de_matrix_condition(mk(sc), config)$mat, logm, info = sc)
    }
    for (sc in c("none", "center")) {
        expect_equal(.metab_de_matrix_condition(mk(sc), config)$mat, work, info = sc)
    }
    out <- .metab_de_matrix_condition(mk("none"), config)
    expect_s3_class(out$condition, "factor")
    expect_equal(as.character(out$condition), meta$sample_type)
})

test_that(".metab_de_matrix_condition drops QC/blanks and matches by sample id", {
    work <- matrix(1:15, nrow = 3,
                   dimnames = list(c("f1","f2","f3"),
                                   c("s1","s2","s3","s4","QC1")))
    # Shuffled row order: alignment must be by id, not position.
    meta <- data.frame(sample_id = c("QC1","s3","s1","s4","s2"),
                       sample_type = c("QC","B","A","B","A"),
                       stringsAsFactors = FALSE)
    pre <- list(expr_work = work, expr_log = work, meta = meta,
                info = list(normalization = list(scaling = "none")))
    config <- list(modes = list(metabolomics = list(
        de = list(condition_column = "sample_type"),
        effects = list(samples = "sample_id", color = "sample_type")
    )))

    out <- .metab_de_matrix_condition(pre, config)
    expect_false("QC1" %in% colnames(out$mat))
    expect_equal(ncol(out$mat), 4L)
    expect_setequal(levels(out$condition), c("A", "B"))
    # condition must line up with the (id-matched) retained columns
    expect_equal(as.character(out$condition),
                 out$meta$sample_type[match(colnames(out$mat), out$meta$sample_id)])
})


# ---- QEA (needs globaltest) --------------------------------------------------

test_that("run_metabolomics_qea runs on the DE matrix and returns a valid table", {
    skip_if_not_installed("globaltest")
    skip_if_not_installed("withr")

    fx <- .make_enrichment_fixture()
    on.exit(unlink(fx$gmt), add = TRUE)

    res <- run_metabolomics_qea(fx$pre, fx$config)
    expect_false(is.null(res))
    expect_true(all(c("pathway", "raw_p", "hits", "FDR") %in% colnames(res$table)))
    expect_identical(res$method, "globaltest_qea")

    mover <- res$table[grepl("^MOVER", res$table$pathway), , drop = FALSE]
    flat  <- res$table[grepl("^FLAT",  res$table$pathway), , drop = FALSE]
    expect_equal(nrow(mover), 1L)
    expect_true(all(mover$raw_p <= flat$raw_p))
})

test_that("run_metabolomics_qea does not read pre$expr_raw anymore", {
    skip_if_not_installed("globaltest")
    skip_if_not_installed("withr")

    fx_clean <- .make_enrichment_fixture(expr_raw_noise = FALSE)
    fx_noise <- .make_enrichment_fixture(expr_raw_noise = TRUE)
    on.exit(unlink(c(fx_clean$gmt, fx_noise$gmt)), add = TRUE)

    m <- function(r) r$table[order(r$table$pathway), c("pathway", "raw_p")]
    # expr_work/expr_log are identical between the two fixtures; only expr_raw
    # differs. Identical results prove the statistic no longer uses expr_raw.
    expect_equal(m(run_metabolomics_qea(fx_noise$pre, fx_noise$config)),
                 m(run_metabolomics_qea(fx_clean$pre, fx_clean$config)))
})

test_that("run_metabolomics_qea excludes QC/blank samples", {
    skip_if_not_installed("globaltest")
    skip_if_not_installed("withr")

    fx_no <- .make_enrichment_fixture(with_qc = FALSE)
    fx_qc <- .make_enrichment_fixture(with_qc = TRUE, shuffle_meta = TRUE)
    on.exit(unlink(c(fx_no$gmt, fx_qc$gmt)), add = TRUE)

    m <- function(r) r$table[order(r$table$pathway), c("pathway", "raw_p")]
    # The extreme QC sample would change results (and add a 3rd level) if kept;
    # biological filtering drops it, so results match the no-QC fixture.
    expect_equal(m(run_metabolomics_qea(fx_qc$pre, fx_qc$config)),
                 m(run_metabolomics_qea(fx_no$pre, fx_no$config)))
})

test_that("run_metabolomics_qea threads standardize through to globaltest::gt", {
    skip_if_not_installed("globaltest")
    skip_if_not_installed("withr")

    # Scale-heterogeneous features so equal-weighting (standardize=TRUE) yields a
    # materially different result than FALSE — proving the flag reaches gt().
    fx_f <- .make_enrichment_fixture(standardize = FALSE, hetero = TRUE)
    fx_t <- .make_enrichment_fixture(standardize = TRUE,  hetero = TRUE)
    on.exit(unlink(c(fx_f$gmt, fx_t$gmt)), add = TRUE)

    p_false <- run_metabolomics_qea(fx_f$pre, fx_f$config)$table
    p_true  <- run_metabolomics_qea(fx_t$pre, fx_t$config)$table
    key <- function(t) t$raw_p[order(t$pathway)]
    expect_false(isTRUE(all.equal(key(p_false), key(p_true))))
})

test_that("run_metabolomics_qea rejects a non-logical standardize", {
    # standardize is validated only after the globaltest availability guard, so
    # this needs globaltest present to reach the check.
    skip_if_not_installed("globaltest")
    skip_if_not_installed("withr")
    fx <- .make_enrichment_fixture(standardize = "yes")
    on.exit(unlink(fx$gmt), add = TRUE)
    expect_error(run_metabolomics_qea(fx$pre, fx$config),
                 "standardize must be TRUE or FALSE")
})

test_that(".run_globaltest forwards standardize to globaltest::gt", {
    skip_if_not_installed("globaltest")
    # Direct seam check on scale-heterogeneous data.
    set.seed(1)
    X <- matrix(stats::rnorm(8 * 6), nrow = 8)
    X <- sweep(X, 2, 10^(0:5), `*`)              # heterogeneous column scales
    colnames(X) <- paste0("c", 1:6)
    response <- factor(rep(c("A", "B"), each = 4))
    subsets  <- list(set1 = c("c1", "c2", "c3"))
    r_f <- .run_globaltest(response, X, subsets, standardize = FALSE)
    r_t <- .run_globaltest(response, X, subsets, standardize = TRUE)
    p_f <- globaltest::result(r_f)[, "p-value"]
    p_t <- globaltest::result(r_t)[, "p-value"]
    expect_false(isTRUE(all.equal(p_f, p_t)))
})


# ---- ssGSEA (needs GSVA) -----------------------------------------------------

test_that("run_metabolomics_ssgsea scores only biological samples (QC excluded)", {
    skip_if_not_installed("GSVA")
    skip_if_not_installed("withr")

    fx <- .make_enrichment_fixture(with_qc = TRUE, shuffle_meta = TRUE)
    on.exit(unlink(fx$gmt), add = TRUE)

    res <- run_metabolomics_ssgsea(fx$pre, fx$config)
    skip_if(is.null(res) || is.null(res$scores), "ssGSEA produced no scores")

    bio_ids <- paste0("s", 1:10)
    expect_false("QC_1" %in% colnames(res$scores))
    expect_equal(ncol(res$scores), length(bio_ids))
    expect_setequal(colnames(res$scores), bio_ids)
})

test_that("run_metabolomics_ssgsea does not read pre$expr_raw anymore", {
    skip_if_not_installed("GSVA")
    skip_if_not_installed("withr")

    fx_clean <- .make_enrichment_fixture(expr_raw_noise = FALSE)
    fx_noise <- .make_enrichment_fixture(expr_raw_noise = TRUE)
    on.exit(unlink(c(fx_clean$gmt, fx_noise$gmt)), add = TRUE)

    r_clean <- run_metabolomics_ssgsea(fx_clean$pre, fx_clean$config)
    r_noise <- run_metabolomics_ssgsea(fx_noise$pre, fx_noise$config)
    skip_if(is.null(r_clean$scores) || is.null(r_noise$scores),
            "ssGSEA produced no scores")

    # Only expr_raw differs between the fixtures; identical scores prove ssGSEA
    # now derives from expr_work, not expr_raw.
    expect_equal(r_noise$scores, r_clean$scores)
})
