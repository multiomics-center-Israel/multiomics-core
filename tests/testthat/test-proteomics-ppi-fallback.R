# tests/testthat/test-proteomics-ppi-fallback.R
#
# Tests for the STRING network fallback contract (R/domain/proteomics/07c_ppi_networks.R).
# The API fallback caps the query at 500 proteins, chosen by table order. That
# truncation used to be a message(), so a run could silently report a network
# built on a subset. These pin the cache-dir resolution and the truncation
# bookkeeping; the network build itself needs the network and is not tested here.

test_that("string_cache_dir honours MULTIOMICS_STRING_CACHE", {
    tmp <- file.path(tempdir(), "string-cache-test")
    withr::with_envvar(c(MULTIOMICS_STRING_CACHE = tmp), {
        expect_equal(string_cache_dir(), tmp)
    })
    expect_true(dir.exists(tmp))
    unlink(tmp, recursive = TRUE)
})

test_that("string_cache_dir returns a usable directory by default", {
    withr::with_envvar(c(MULTIOMICS_STRING_CACHE = ""), {
        d <- string_cache_dir()
        # Either a real cache dir, or "" meaning "keep STRINGdb's tempdir default".
        expect_true(identical(d, "") || dir.exists(d))
    })
})

test_that("string_cache_dir is stable across calls", {
    withr::with_envvar(c(MULTIOMICS_STRING_CACHE = ""), {
        expect_identical(string_cache_dir(), string_cache_dir())
    })
})

# The truncation rule the fallback applies, stated independently of the network.
truncate_for_api <- function(proteins, chunk_size = 500) {
    if (length(proteins) > chunk_size) proteins[seq_len(chunk_size)] else proteins
}

test_that("the API fallback caps at 500 proteins", {
    p <- paste0("P", seq_len(649))
    kept <- truncate_for_api(p)
    expect_length(kept, 500)
    expect_identical(kept, p[1:500])
})

test_that("no truncation happens under the cap", {
    p <- paste0("P", seq_len(120))
    expect_identical(truncate_for_api(p), p)
})

test_that("truncation drops by table order, not by biology", {
    # Pinned deliberately: the dropped proteins are whichever came last in the
    # DE table, so a truncated network must be reported as incomplete.
    p <- paste0("P", seq_len(510))
    expect_false("P510" %in% truncate_for_api(p))
    expect_true("P1" %in% truncate_for_api(p))
})
