# tests/testthat/test-metab-expr-norm-missing-mask.R
#
# Schema 2.1 adds `expr_norm_missing` to the Shiny payload contract.
#
# `expr_norm` is declared NA-free, and the metabolomics builder enforces that by
# replacing missing cells with row medians. Those filled values are
# indistinguishable from measurements, so any pairwise-complete analysis
# downstream computes on a partly synthetic sample at the wrong effective n.
# The mask records what was filled, and the three states must stay
# distinguishable:
#
#   logical matrix, any TRUE  -> missingness existed and was recorded
#   logical matrix, all FALSE -> schema >= 2.1 verified the data complete
#   NULL / absent             -> pre-2.1 payload; provenance unknown
#
# A NULL in the second slot would be indistinguishable from the third, which is
# the specific regression these tests exist to prevent.


make_payload_inputs <- function(expr_work) {
    meta <- data.frame(
        SampleID  = colnames(expr_work),
        treatment = rep(c("ctrl", "trt"), length.out = ncol(expr_work)),
        stringsAsFactors = FALSE
    )
    row_data <- data.frame(
        feature_id = rownames(expr_work),
        Name       = paste0("compound_", seq_len(nrow(expr_work))),
        stringsAsFactors = FALSE
    )
    pre <- list(
        expr_raw  = expr_work,
        expr_filt = expr_work,
        expr_log  = expr_work,
        expr_work = expr_work,
        meta      = meta,
        row_data  = row_data
    )
    config <- list(modes = list(metabolomics = list(
        effects      = list(samples = "SampleID", color = "treatment"),
        de           = list(condition_column = "treatment", p_cutoff = 0.05),
        preprocessing = list(chosen_norm = "none")
    )))
    list(pre = pre, config = config, inputs = list(contrasts = NULL))
}

toy_expr <- function() {
    m <- matrix(
        c(1, 2, 3, 4,
          5, 6, 7, 8,
          9, 10, 11, 12),
        nrow = 3, byrow = TRUE,
        dimnames = list(c("F1", "F2", "F3"), paste0("S", 1:4))
    )
    m
}

build <- function(expr_work) {
    fx <- make_payload_inputs(expr_work)
    suppressWarnings(suppressMessages(
        build_shiny_payload_metabolomics(
            pre = fx$pre, de_res = NULL, inputs = fx$inputs, config = fx$config,
            include_legacy = FALSE
        )
    ))
}


# ---- the three states --------------------------------------------------------

test_that("missing values are recorded in the mask before the fill overwrites them", {
    m <- toy_expr()
    m["F2", "S3"] <- NA
    m["F3", "S1"] <- NA

    payload <- build(m)

    expect_true(is.matrix(payload$expr_norm_missing))
    expect_true(is.logical(payload$expr_norm_missing))
    expect_identical(payload$expr_norm_missing, !is.finite(m))
    expect_equal(sum(payload$expr_norm_missing), 2L)
    expect_identical(dimnames(payload$expr_norm_missing), dimnames(m))

    # The fill still happens -- expr_norm keeps its NA-free guarantee.
    expect_false(anyNA(payload$expr_norm))
})

test_that("a complete matrix yields an all-FALSE mask, NOT NULL", {
    # This is the state rev5 exists to make expressible: "verified complete" has
    # to be distinguishable from "we don't know", and NULL cannot carry both.
    payload <- build(toy_expr())

    expect_true("expr_norm_missing" %in% names(payload))
    expect_false(is.null(payload$expr_norm_missing))
    expect_true(is.matrix(payload$expr_norm_missing))
    expect_true(is.logical(payload$expr_norm_missing))
    expect_false(any(payload$expr_norm_missing))
    expect_identical(dim(payload$expr_norm_missing), dim(toy_expr()))
})

test_that("the key is present in names() either way", {
    # `payload$key <- NULL` deletes the entry rather than blanking it, so a
    # naive is.null() check would pass while the key had silently vanished.
    for (m in list(toy_expr(), {x <- toy_expr(); x[1, 1] <- NA; x})) {
        payload <- build(m)
        expect_true("expr_norm_missing" %in% names(payload))
    }
})


# ---- non-finite values -------------------------------------------------------

test_that("-Inf is captured by the mask even though the fill leaves it behind", {
    # transform_metab() emits -Inf when log-transforming a non-positive value.
    # is.na() does not catch it and the row-median fill does not replace it, so
    # only an is.finite()-based mask keeps it out of downstream analysis.
    m <- toy_expr()
    m["F1", "S2"] <- -Inf

    payload <- build(m)

    expect_true(payload$expr_norm_missing["F1", "S2"])
    expect_equal(sum(payload$expr_norm_missing), 1L)

    # Restoring turns it into NA, so every consumer -- correlation and plot
    # alike -- sees one shared definition of an observed value.
    restored <- restore_missing_values(payload$expr_norm,
                                       payload$expr_norm_missing)
    expect_true(is.na(restored["F1", "S2"]))
    expect_false(any(is.infinite(restored)))
})

test_that("NaN is captured too", {
    m <- toy_expr()
    m["F3", "S4"] <- NaN

    payload <- build(m)
    expect_true(payload$expr_norm_missing["F3", "S4"])
})


# ---- schema and contract validation -----------------------------------------

test_that("the payload declares schema 2.1", {
    payload <- build(toy_expr())
    expect_equal(payload$payload_version, "2.1")
})

test_that("a freshly built payload passes strict contract validation", {
    payload <- build(toy_expr())
    expect_silent(
        assert_shiny_payload_contract(payload, strict = TRUE,
                                      context = "metabolomics")
    )
})

test_that("a malformed mask is rejected by the contract", {
    payload <- build(toy_expr())

    numeric_mask <- payload
    numeric_mask$expr_norm_missing <- matrix(
        0, nrow = nrow(payload$expr_norm), ncol = ncol(payload$expr_norm),
        dimnames = dimnames(payload$expr_norm)
    )
    expect_error(
        assert_shiny_payload_contract(numeric_mask, strict = TRUE),
        "logical matrix"
    )

    na_mask <- payload
    na_mask$expr_norm_missing[1, 1] <- NA
    expect_error(
        assert_shiny_payload_contract(na_mask, strict = TRUE),
        "must not contain NA"
    )

    wrong_dim <- payload
    wrong_dim$expr_norm_missing <- t(payload$expr_norm_missing)
    expect_error(
        assert_shiny_payload_contract(wrong_dim, strict = TRUE),
        "must match expr_norm"
    )

    renamed <- payload
    colnames(renamed$expr_norm_missing) <- paste0("X", seq_len(ncol(payload$expr_norm)))
    expect_error(
        assert_shiny_payload_contract(renamed, strict = TRUE),
        "dimnames must match"
    )
})

test_that("a metabolomics payload at >= 2.1 must carry the mask", {
    # Without this the version bump would promise a provenance guarantee that
    # nothing actually enforces.
    payload <- build(toy_expr())
    payload["expr_norm_missing"] <- list(NULL)

    expect_error(
        assert_shiny_payload_contract(payload, strict = TRUE),
        "must carry expr_norm_missing"
    )
})

test_that("the >= 2.1 gate compares versions numerically, not lexically", {
    # "2.9" < "2.10" is FALSE as strings, so a string gate would stop firing at
    # 2.10. numeric_version() orders them correctly.
    payload <- build(toy_expr())
    payload["expr_norm_missing"] <- list(NULL)

    payload$payload_version <- "2.10"
    expect_error(assert_shiny_payload_contract(payload, strict = TRUE),
                 "must carry expr_norm_missing")

    payload$payload_version <- "3.0"
    expect_error(assert_shiny_payload_contract(payload, strict = TRUE),
                 "must carry expr_norm_missing")
})

test_that("the requirement does not fire for other omics", {
    # rnaseq and proteomics perform no fill, so a NULL mask there means
    # "not applicable", not "provenance lost".
    payload <- build(toy_expr())
    payload["expr_norm_missing"] <- list(NULL)
    payload$payload_source <- "rnaseq"

    expect_silent(assert_shiny_payload_contract(payload, strict = TRUE))
})


# ---- backward compatibility with v2.0 payloads on disk -----------------------

test_that("a pre-2.1 payload is readable and degrades honestly", {
    payload <- build(toy_expr())

    # Simulate one written before the key existed.
    legacy <- payload
    legacy["expr_norm_missing"] <- list(NULL)
    legacy$payload_version <- "2.0"

    # The load path never asserts, so opening an old payload still works.
    f <- tempfile(fileext = ".rds")
    on.exit(unlink(f), add = TRUE)
    saveRDS(legacy, f)
    reloaded <- suppressMessages(load_shiny_payload_metabolomics(f, verbose = FALSE))
    expect_equal(reloaded$payload_version, "2.0")

    # The >= 2.1 requirement does not fire on a 2.0 payload...
    msgs <- capture_warnings(assert_shiny_payload_contract(legacy, strict = FALSE))
    expect_false(any(grepl("must carry expr_norm_missing", msgs)))

    # ...and restoring is a no-op, so the GUI degrades rather than crashing.
    expect_identical(
        restore_missing_values(legacy$expr_norm, legacy$expr_norm_missing),
        legacy$expr_norm
    )
})

test_that("a payload with no version at all is treated as legacy", {
    payload <- build(toy_expr())
    payload["expr_norm_missing"] <- list(NULL)
    payload["payload_version"] <- list(NULL)

    msgs <- capture_warnings(assert_shiny_payload_contract(payload, strict = FALSE))
    expect_false(any(grepl("must carry expr_norm_missing", msgs)))
})
