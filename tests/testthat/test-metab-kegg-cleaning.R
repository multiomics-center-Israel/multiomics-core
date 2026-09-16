# Tests for clean_kegg_chemspider(): the KEGG annotation column can arrive with
# ChemSpider IDs (CSID…) mixed in with real KEGG compound IDs (C#####). The
# cleaner keeps only valid KEGG ids in `KEGG` and routes CSID ids to a separate
# `ChemSpider` column.

test_that("clean_kegg_chemspider splits KEGG and ChemSpider correctly", {
    row_data <- data.frame(
        feature_id = paste0("F", 1:6),
        KEGG = c("C00031", "CSID144529", "cpd:C00267", "", NA, "foo"),
        Name = paste0("m", 1:6),
        stringsAsFactors = FALSE
    )

    out <- clean_kegg_chemspider(row_data)

    # Valid KEGG kept as-is; embedded form extracted to the bare C##### id;
    # CSID / empty / junk -> NA in KEGG.
    expect_equal(out$KEGG, c("C00031", NA, "C00267", NA, NA, NA))
    # CSID routed to ChemSpider; everything else NA there.
    expect_equal(out$ChemSpider, c(NA, "CSID144529", NA, NA, NA, NA))
    # Other columns and row count preserved.
    expect_equal(out$feature_id, row_data$feature_id)
    expect_equal(out$Name, row_data$Name)
    expect_equal(nrow(out), 6L)
})

test_that("clean_kegg_chemspider places ChemSpider immediately after KEGG", {
    row_data <- data.frame(
        feature_id = "F1", mz = 100, KEGG = "C00031", RT = 1.2,
        stringsAsFactors = FALSE
    )
    out <- clean_kegg_chemspider(row_data)
    cols <- colnames(out)
    expect_equal(cols[match("KEGG", cols) + 1L], "ChemSpider")
    # No columns lost.
    expect_setequal(cols, c("feature_id", "mz", "KEGG", "RT", "ChemSpider"))
})

test_that("clean_kegg_chemspider is a no-op without a KEGG column", {
    row_data <- data.frame(feature_id = c("F1", "F2"), Name = c("a", "b"),
                           stringsAsFactors = FALSE)
    expect_identical(clean_kegg_chemspider(row_data), row_data)
    expect_null(clean_kegg_chemspider(NULL))
})

test_that("clean_kegg_chemspider preserves a pre-existing ChemSpider column", {
    row_data <- data.frame(
        feature_id = c("F1", "F2", "F3"),
        KEGG       = c("C00031", "CSID999", "C00267"),
        ChemSpider = c("CSID111", NA, ""),   # F1 already has a curated CSID
        stringsAsFactors = FALSE
    )
    out <- clean_kegg_chemspider(row_data)
    # F1: keep curated CSID111 (not clobbered); F2: routed from KEGG; F3: none.
    expect_equal(out$ChemSpider, c("CSID111", "CSID999", NA))
    expect_equal(out$KEGG, c("C00031", NA, "C00267"))
})

test_that("clean_kegg_chemspider reproduces the 438 -> 315 KEGG / 123 CSID split", {
    kegg_ids <- sprintf("C%05d", 1:315)
    csid_ids <- sprintf("CSID%d", 1:123)
    blanks   <- rep(NA_character_, 20)               # genuinely unannotated
    vals     <- c(kegg_ids, csid_ids, blanks)

    row_data <- data.frame(
        feature_id = paste0("F", seq_along(vals)),
        KEGG = vals,
        stringsAsFactors = FALSE
    )
    expect_equal(sum(!is.na(row_data$KEGG) & nzchar(row_data$KEGG)), 438L)  # before

    out <- clean_kegg_chemspider(row_data)
    expect_equal(sum(!is.na(out$KEGG)), 315L)         # after: only real KEGG
    expect_equal(sum(!is.na(out$ChemSpider)), 123L)   # CSID routed out
    expect_true(all(grepl("^C[0-9]{5}$", out$KEGG[!is.na(out$KEGG)])))
})
