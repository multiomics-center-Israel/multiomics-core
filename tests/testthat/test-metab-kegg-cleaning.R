# Tests for clean_kegg_chemspider(): the KEGG annotation column can arrive with
# ChemSpider IDs (CSID…) mixed in with real KEGG compound IDs (C#####). The
# cleaner keeps only valid KEGG ids in `KEGG` and routes CSID ids to a separate
# `ChemSpider` column.

test_that("clean_kegg_chemspider splits KEGG and ChemSpider correctly", {
    row_data <- data.frame(
        feature_id = paste0("F", 1:6),
        KEGG = c("C00031", "CSID900000", "cpd:C00267", "", NA, "foo"),
        Name = paste0("m", 1:6),
        stringsAsFactors = FALSE
    )

    out <- clean_kegg_chemspider(row_data)

    # Valid KEGG kept as-is; embedded form extracted to the bare C##### id;
    # CSID / empty / junk -> NA in KEGG.
    expect_equal(out$KEGG, c("C00031", NA, "C00267", NA, NA, NA))
    # CSID routed to ChemSpider; everything else NA there.
    expect_equal(out$ChemSpider, c(NA, "CSID900000", NA, NA, NA, NA))
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

test_that("clean_kegg_chemspider makes a non-empty KEGG count mean real KEGG coverage", {
    # Synthetic: 5 real KEGG ids, 3 ChemSpider ids and 2 unannotated features.
    kegg_ids <- sprintf("C%05d", 1:5)
    csid_ids <- sprintf("CSID%d", 900001:900003)
    blanks   <- rep(NA_character_, 2)
    vals     <- c(kegg_ids, csid_ids, blanks)

    row_data <- data.frame(
        feature_id = paste0("F", seq_along(vals)),
        KEGG = vals,
        stringsAsFactors = FALSE
    )
    # Before: a naive non-empty count also counts the ChemSpider ids.
    expect_equal(sum(!is.na(row_data$KEGG) & nzchar(row_data$KEGG)), 8L)

    out <- clean_kegg_chemspider(row_data)
    expect_equal(sum(!is.na(out$KEGG)), 5L)           # after: only real KEGG
    expect_equal(sum(!is.na(out$ChemSpider)), 3L)     # CSID routed out
    expect_true(all(grepl("^C[0-9]{5}$", out$KEGG[!is.na(out$KEGG)])))
})

test_that("a KEGG id available from HMDB survives a CSID being routed out of KEGG", {
    skip_if_not_installed("readr")
    # F1 carries a ChemSpider id in the KEGG column but has an HMDB id the
    # mapping knows. F2 already has a real KEGG id, which the fill must not
    # replace. F3 is blank and filled from a mapping value that is not a bare
    # C##### id, which the second clean has to sanitize.
    map_file <- tempfile(fileext = ".tsv")
    on.exit(unlink(map_file), add = TRUE)
    writeLines(c("HMDB\tKEGG",
                 "HMDB0000001\tC00022",
                 "HMDB0000002\tC99999",
                 "HMDB0000003\tcpd:C00158"), map_file)

    row_data <- data.frame(
        feature_id = c("F1", "F2", "F3"),
        HMDB       = c("HMDB0000001", "HMDB0000002", "HMDB0000003"),
        KEGG       = c("CSID900001", "C00031", NA),
        stringsAsFactors = FALSE
    )

    # The order mod_met_raw() uses: clean, fill, clean.
    out <- suppressMessages(clean_kegg_chemspider(
        add_kegg_from_hmdb(clean_kegg_chemspider(row_data), map_file)))
    expect_equal(out$KEGG, c("C00022", "C00031", "C00158"))
    expect_equal(out$ChemSpider, c("CSID900001", NA, NA))

    # Filling first loses F1's KEGG id: the CSID blocks the lookup and is then
    # routed out, leaving the cell empty. Pinned so the test above cannot pass
    # for a reason other than the order.
    fill_first <- suppressMessages(clean_kegg_chemspider(
        add_kegg_from_hmdb(row_data, map_file)))
    expect_true(is.na(fill_first$KEGG[1]))
})

test_that("mod_met_raw cleans the KEGG column before and after the HMDB fill", {
    src  <- readLines(file.path(root_dir, "R", "modules", "metabolomics",
                                "00_mod_preprocessing.R"), warn = FALSE)
    code <- sub("#.*$", "", src)
    clean_at <- grep("clean_kegg_chemspider\\(", code)
    fill_at  <- grep("add_kegg_from_hmdb\\(", code)

    expect_length(fill_at, 1L)
    expect_true(any(clean_at < fill_at))
    expect_true(any(clean_at > fill_at))
})

test_that("clean_kegg_chemspider rejects a KEGG-like id with too many digits", {
    row_data <- data.frame(
        feature_id = c("F1", "F2", "F3"),
        KEGG       = c("C000311", "C000311;C00267", "cpd:C00031"),
        stringsAsFactors = FALSE
    )
    out <- suppressMessages(clean_kegg_chemspider(row_data))
    # "C000311" is not C00031; the next valid id in the cell is used instead.
    expect_equal(out$KEGG, c(NA, "C00267", "C00031"))
    expect_message(clean_kegg_chemspider(row_data), "1 other dropped")
})

test_that("clean_kegg_chemspider counts only CSIDs it moves out of KEGG", {
    row_data <- data.frame(
        feature_id = c("F1", "F2"),
        KEGG       = c("CSID900001", "C00031"),
        ChemSpider = c(NA, "CSID900002"),     # F2 already carries a curated id
        stringsAsFactors = FALSE
    )
    expect_message(first <- clean_kegg_chemspider(row_data), "1 ChemSpider routed")
    expect_equal(first$ChemSpider, c("CSID900001", "CSID900002"))

    # A second pass, as mod_met_raw() runs after the HMDB fill, routes nothing.
    expect_message(clean_kegg_chemspider(first), "0 ChemSpider routed")
})

test_that("clean_kegg_chemspider extracts a KEGG id from a line-separated cell", {
    row_data <- data.frame(
        feature_id = c("F1", "F2"),
        KEGG       = c("label\nC00031", "C00031\nC00267"),
        stringsAsFactors = FALSE
    )
    out <- suppressMessages(clean_kegg_chemspider(row_data))
    expect_equal(out$KEGG, c("C00031", "C00031"))
})

test_that("clean_kegg_chemspider cleans an alternative KEGG column in place", {
    row_data <- data.frame(
        feature_id = c("F1", "F2", "F3"),
        KEGG_ID    = c("CSID900011", "cpd:C00022", NA),
        Name       = c("m1", "m2", "m3"),
        stringsAsFactors = FALSE
    )
    out <- suppressMessages(clean_kegg_chemspider(row_data))

    # Cleaned under its own name; no KEGG column is invented.
    expect_equal(out$KEGG_ID, c(NA, "C00022", NA))
    expect_false("KEGG" %in% colnames(out))
    expect_equal(out$ChemSpider, c("CSID900011", NA, NA))
    cols <- colnames(out)
    expect_equal(cols[match("KEGG_ID", cols) + 1L], "ChemSpider")
    expect_setequal(cols, c("feature_id", "KEGG_ID", "Name", "ChemSpider"))
})

test_that("clean_kegg_chemspider cleans KEGG and an alternative column together", {
    row_data <- data.frame(
        feature_id = c("F1", "F2", "F3"),
        KEGG       = c("C00031", "CSID900021", NA),
        `KEGG ID`  = c("CSID900022", "C00267", "CSID900023"),
        ChemSpider = c(NA, NA, "CSID900099"),   # F3 already carries a curated id
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_message(out <- clean_kegg_chemspider(row_data), "2 ChemSpider routed")

    expect_equal(out$KEGG, c("C00031", NA, NA))
    expect_equal(out[["KEGG ID"]], c(NA, "C00267", NA))
    # F1 routed from the alternative column; F2 from KEGG, the first in order;
    # F3 keeps its curated value.
    expect_equal(out$ChemSpider, c("CSID900022", "CSID900021", "CSID900099"))
    expect_setequal(colnames(out), colnames(row_data))
})
