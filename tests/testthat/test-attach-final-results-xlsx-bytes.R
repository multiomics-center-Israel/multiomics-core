# Regression guard for attach_final_results_xlsx_bytes().
# Confirms:
#   - both keys are populated as raw vectors when the two xlsx paths are present
#   - keys stay NULL when xlsx_files is NULL / empty / no pattern match
#   - missing file on disk falls back to NULL (no crash)
#   - existing payload keys (incl. de_final_table) are NOT clobbered

make_xlsx_fixture <- function(dir, p_tag = "0.05") {
    all_path <- file.path(dir, sprintf("Final_results_ALL_P_%s.xlsx", p_tag))
    de_path  <- file.path(dir, sprintf("Final_results_DE_P_%s.xlsx",  p_tag))
    # The helper does not parse the xlsx; raw bytes are enough for the test.
    writeBin(as.raw(c(0x50, 0x4b, 0x03, 0x04, 1:10)), all_path)
    writeBin(as.raw(c(0x50, 0x4b, 0x03, 0x04, 11:20)), de_path)
    list(all = all_path, de = de_path)
}

test_that("attach_final_results_xlsx_bytes populates both keys from real files", {
    tmp <- tempfile("xlsx_fixture_"); dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    paths <- make_xlsx_fixture(tmp)

    payload <- list(de_final_table = data.frame(x = 1))
    payload <- attach_final_results_xlsx_bytes(payload, c(paths$all, paths$de))

    expect_true(is.raw(payload$all_final_xlsx))
    expect_true(is.raw(payload$de_final_xlsx))
    expect_equal(length(payload$all_final_xlsx), file.info(paths$all)$size)
    expect_equal(length(payload$de_final_xlsx),  file.info(paths$de)$size)
    expect_s3_class(payload$de_final_table, "data.frame")  # untouched
})

test_that("attach_final_results_xlsx_bytes handles unrelated paths gracefully", {
    payload <- list()
    payload <- attach_final_results_xlsx_bytes(
        payload,
        c("/tmp/some_other.tsv", "/tmp/Final_results_summary.csv")
    )
    expect_null(payload$all_final_xlsx)
    expect_null(payload$de_final_xlsx)
})

test_that("attach_final_results_xlsx_bytes is NULL-safe", {
    p1 <- attach_final_results_xlsx_bytes(list(), NULL)
    p2 <- attach_final_results_xlsx_bytes(list(), character(0))
    expect_null(p1$all_final_xlsx); expect_null(p1$de_final_xlsx)
    expect_null(p2$all_final_xlsx); expect_null(p2$de_final_xlsx)
})

test_that("attach_final_results_xlsx_bytes tolerates missing files on disk", {
    payload <- list()
    payload <- attach_final_results_xlsx_bytes(
        payload,
        c("/nonexistent/Final_results_ALL_P_0.05.xlsx",
          "/nonexistent/Final_results_DE_P_0.05.xlsx")
    )
    expect_null(payload$all_final_xlsx)
    expect_null(payload$de_final_xlsx)
})

test_that("attach_final_results_xlsx_bytes picks first match when multiple present", {
    tmp <- tempfile("xlsx_fixture_"); dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    p1 <- make_xlsx_fixture(tmp, p_tag = "0.05")
    # Second pair with a different p_tag in the same vector
    p2 <- make_xlsx_fixture(tmp, p_tag = "0.01")

    payload <- attach_final_results_xlsx_bytes(
        list(),
        c(p1$all, p2$all, p1$de, p2$de)
    )
    expect_true(is.raw(payload$all_final_xlsx))
    expect_true(is.raw(payload$de_final_xlsx))
})
