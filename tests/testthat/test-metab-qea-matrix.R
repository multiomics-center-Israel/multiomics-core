# Tests for the QEA/ssGSEA matrix fix: both now run on the same analysis matrix
# as DE (log-normalized, biological samples only) via .metab_de_matrix_condition(),
# instead of pre$expr_raw. Also covers the QEA `standardize` option.

# ---- .metab_de_matrix_condition() -------------------------------------------

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

    # Variance scalings use the pre-scaling matrix (expr_log)
    for (sc in c("auto", "pareto", "range")) {
        out <- .metab_de_matrix_condition(mk(sc), config)
        expect_equal(out$mat, logm, info = sc)
    }
    # "none"/"center" use expr_work
    for (sc in c("none", "center")) {
        out <- .metab_de_matrix_condition(mk(sc), config)
        expect_equal(out$mat, work, info = sc)
    }
    out <- .metab_de_matrix_condition(mk("none"), config)
    expect_s3_class(out$condition, "factor")
    expect_equal(as.character(out$condition), meta$sample_type)
})

test_that(".metab_de_matrix_condition excludes QC/blank samples", {
    work <- matrix(1:15, nrow = 3,
                   dimnames = list(c("f1","f2","f3"),
                                   c("s1","s2","s3","s4","QC1")))
    meta <- data.frame(sample_id = c("s1","s2","s3","s4","QC1"),
                       sample_type = c("A","A","B","B","QC"),
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
})


# ---- QEA integration (needs globaltest) -------------------------------------

.make_qea_fixture <- function(with_qc = FALSE, standardize = FALSE) {
    groups  <- rep(c("B", "A"), each = 5L)
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

    meta <- data.frame(sample_id = samples, sample_type = groups,
                       stringsAsFactors = FALSE)

    if (with_qc) {
        # Add a QC sample that biological filtering must drop; its values are
        # deliberately extreme so it would perturb results if NOT excluded.
        qc_col <- matrix(rep(50, n_feat), ncol = 1,
                         dimnames = list(feat_ids, "QC1"))
        mat  <- cbind(mat, qc_col)
        meta <- rbind(meta, data.frame(sample_id = "QC1", sample_type = "QC"))
    }

    row_data <- data.frame(feature_id = feat_ids, KEGG = kegg_ids,
                           stringsAsFactors = FALSE)
    rownames(row_data) <- feat_ids

    pre <- list(expr_work = mat, expr_log = mat, expr_filt = mat,
                meta = meta, row_data = row_data,
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

test_that("run_metabolomics_qea runs on the DE matrix and returns a valid table", {
    skip_if_not_installed("globaltest")
    skip_if_not_installed("withr")

    fx <- .make_qea_fixture()
    on.exit(unlink(fx$gmt), add = TRUE)

    res <- run_metabolomics_qea(fx$pre, fx$config)
    expect_false(is.null(res))
    expect_true(all(c("pathway", "raw_p", "hits", "FDR") %in% colnames(res$table)))
    expect_identical(res$method, "globaltest_qea")

    # The moving pathway is more strongly associated than the flat pathways.
    mover <- res$table[grepl("^MOVER", res$table$pathway), , drop = FALSE]
    flat  <- res$table[grepl("^FLAT",  res$table$pathway), , drop = FALSE]
    expect_equal(nrow(mover), 1L)
    expect_true(all(mover$raw_p <= flat$raw_p))
})

test_that("run_metabolomics_qea excludes QC/blank samples (matrix now matches DE)", {
    skip_if_not_installed("globaltest")
    skip_if_not_installed("withr")

    fx_no <- .make_qea_fixture(with_qc = FALSE)
    fx_qc <- .make_qea_fixture(with_qc = TRUE)
    on.exit(unlink(c(fx_no$gmt, fx_qc$gmt)), add = TRUE)

    res_no <- run_metabolomics_qea(fx_no$pre, fx_no$config)
    res_qc <- run_metabolomics_qea(fx_qc$pre, fx_qc$config)

    # An extreme QC sample would change results if it were included; because
    # biological filtering drops it, the pathway p-values match.
    m <- function(r) r$table[order(r$table$pathway), c("pathway", "raw_p")]
    expect_equal(m(res_qc), m(res_no))
})

test_that("run_metabolomics_qea honors the standardize option", {
    skip_if_not_installed("globaltest")
    skip_if_not_installed("withr")

    fx <- .make_qea_fixture(standardize = TRUE)
    on.exit(unlink(fx$gmt), add = TRUE)

    res <- run_metabolomics_qea(fx$pre, fx$config)
    expect_false(is.null(res))
    expect_true("MOVER - Moving pathway" %in% res$table$pathway)
})
