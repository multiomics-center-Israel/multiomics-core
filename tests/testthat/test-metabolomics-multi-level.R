# tests/testthat/test-metabolomics-multi-level.R
#
# Unit and integration tests for the "multi_level" metabolomics input format.
#
# Strategy:
#   - All unit tests operate on pre-built R objects (no Excel, no filesystem I/O).
#     They directly exercise merge_level_parsed() and parse_multi_level() using
#     synthetic parse results constructed inline.
#   - One integration test writes a pair of real XLSX files to a temp directory
#     and exercises read_multi_level_dir() end-to-end.
#
# Functions under test (all in R/domain/metabolomics/00_inputs.R):
#   merge_level_parsed()           — pure merge logic
#   parse_multi_level()            — per-level dispatch + sample validation
#   read_multi_level_dir()         — directory enumeration (integration only)
#   validate_metabolomics_config() — config validation guards


# ==============================================================================
# Test fixtures — synthetic parse results (no Excel dependency)
# ==============================================================================

# Helper: build a minimal parse result like parse_cd_raw() would return.
.make_parsed <- function(features, samples, extra_annot = list()) {
    n_feat <- length(features)
    n_samp <- length(samples)

    mat <- matrix(
        seq_len(n_feat * n_samp) * 10,
        nrow     = n_feat,
        ncol     = n_samp,
        dimnames = list(features, samples)
    )
    storage.mode(mat) <- "double"

    rd <- data.frame(
        feature_id = features,
        Name       = paste0("Compound_", features),
        stringsAsFactors = FALSE
    )
    for (col in names(extra_annot)) rd[[col]] <- extra_annot[[col]]

    sm <- data.frame(
        cd_column = paste0("Area: ", samples, ".raw (F1)"),
        sample_id = samples,
        stringsAsFactors = FALSE
    )

    list(expr_raw = mat, row_data = rd, sample_ids = samples, sample_map = sm)
}

.samples <- c("S1", "S2", "S3")

.lv1 <- .make_parsed(
    features = c("Glucose_mz180.0634_rt1.23", "Fructose_mz180.0634_rt2.10"),
    samples  = .samples
)
.lv2 <- .make_parsed(
    features = c("Alanine_mz89.0477_rt0.45", "Valine_mz117.0790_rt0.88"),
    samples  = .samples
)


# ==============================================================================
# 1. Contract conformance
# ==============================================================================

test_that("merge_level_parsed returns all required contract keys", {
    result <- merge_level_parsed(
        parsed_levels = list(Level_1 = .lv1, Level_2 = .lv2),
        level_names   = c("Level_1", "Level_2")
    )
    expect_named(result, c("expr_raw", "row_data", "sample_ids", "sample_map"))
    expect_true(is.matrix(result$expr_raw))
    expect_true(is.numeric(result$expr_raw))
    expect_true(is.data.frame(result$row_data))
    expect_true(is.character(result$sample_ids))
})

test_that("merge_level_parsed produces correct dimensions", {
    result <- merge_level_parsed(
        parsed_levels = list(Level_1 = .lv1, Level_2 = .lv2),
        level_names   = c("Level_1", "Level_2")
    )
    expect_equal(nrow(result$expr_raw), 4L)   # 2 features per level
    expect_equal(ncol(result$expr_raw), 3L)   # 3 samples
    expect_equal(nrow(result$row_data), 4L)
    expect_equal(result$sample_ids, .samples)
})


# ==============================================================================
# 2. Source_File injection
# ==============================================================================

test_that("Source_File is populated correctly for each level", {
    result <- merge_level_parsed(
        parsed_levels = list(Level_1 = .lv1, Level_2 = .lv2),
        level_names   = c("Level_1", "Level_2")
    )
    rd <- result$row_data
    expect_true("Source_File" %in% colnames(rd))

    lv1_rows <- startsWith(rd$feature_id, "Level_1__")
    lv2_rows <- startsWith(rd$feature_id, "Level_2__")
    expect_true(all(rd$Source_File[lv1_rows] == "Level_1"))
    expect_true(all(rd$Source_File[lv2_rows] == "Level_2"))
})


# ==============================================================================
# 3. Cross-level unique IDs via Level__ prefix
# ==============================================================================

test_that("feature IDs are prefixed and unique even when raw IDs collide across levels", {
    shared_feats <- c("Feature_A_mz100.0000_rt1.00", "Feature_B_mz200.0000_rt2.00")
    lv_a <- .make_parsed(shared_feats, .samples)
    lv_b <- .make_parsed(shared_feats, .samples)

    result <- merge_level_parsed(
        parsed_levels = list(Level_1 = lv_a, Level_2 = lv_b),
        level_names   = c("Level_1", "Level_2")
    )

    ids <- rownames(result$expr_raw)
    expect_equal(length(ids), length(unique(ids)))   # all unique
    expect_true(all(startsWith(ids[1:2], "Level_1__")))
    expect_true(all(startsWith(ids[3:4], "Level_2__")))
})

test_that("feature_id in row_data matches rownames of expr_raw", {
    result <- merge_level_parsed(
        parsed_levels = list(Level_1 = .lv1, Level_2 = .lv2),
        level_names   = c("Level_1", "Level_2")
    )
    expect_identical(rownames(result$expr_raw), result$row_data$feature_id)
})


# ==============================================================================
# 4. feature_id_orig traceability (mandatory)
# ==============================================================================

test_that("feature_id_orig is present and equals the pre-prefix feature_id", {
    result <- merge_level_parsed(
        parsed_levels = list(Level_1 = .lv1, Level_2 = .lv2),
        level_names   = c("Level_1", "Level_2")
    )
    rd <- result$row_data
    expect_true("feature_id_orig" %in% colnames(rd))

    # Strip level prefix and compare to feature_id_orig
    reconstructed <- sub("^[^_]+__", "", rd$feature_id)
    expect_identical(rd$feature_id_orig, reconstructed)
})


# ==============================================================================
# 5. Deterministic column ordering
# ==============================================================================

test_that("row_data first three columns are always feature_id, Source_File, feature_id_orig", {
    result <- merge_level_parsed(
        parsed_levels = list(Level_1 = .lv1, Level_2 = .lv2),
        level_names   = c("Level_1", "Level_2")
    )
    expect_equal(colnames(result$row_data)[1:3],
                 c("feature_id", "Source_File", "feature_id_orig"))
})

test_that("column order is identical regardless of level list order", {
    result_ab <- merge_level_parsed(
        parsed_levels = list(Level_1 = .lv1, Level_2 = .lv2),
        level_names   = c("Level_1", "Level_2")
    )
    result_ba <- merge_level_parsed(
        parsed_levels = list(Level_2 = .lv2, Level_1 = .lv1),
        level_names   = c("Level_2", "Level_1")
    )
    expect_identical(colnames(result_ab$row_data), colnames(result_ba$row_data))
})


# ==============================================================================
# 6. Column mismatch — missing annotation columns filled with NA
# ==============================================================================

test_that("annotation column present in level 1 but absent in level 2 is NA for level 2", {
    lv1_extra <- .make_parsed(
        features    = c("G_mz180.0634_rt1.23"),
        samples     = .samples,
        extra_annot = list(Formula = "C6H12O6")
    )
    lv2_plain <- .make_parsed(c("A_mz89.0477_rt0.45"), .samples)

    result <- merge_level_parsed(
        parsed_levels = list(Level_1 = lv1_extra, Level_2 = lv2_plain),
        level_names   = c("Level_1", "Level_2")
    )

    rd <- result$row_data
    expect_true("Formula" %in% colnames(rd))
    expect_equal(rd$Formula[rd$Source_File == "Level_1"], "C6H12O6")
    expect_true(is.na(rd$Formula[rd$Source_File == "Level_2"]))
})

test_that("annotation column appearing only in level 2 is NA for level 1", {
    lv1_plain <- .make_parsed(c("G_mz180.0634_rt1.00"), .samples)
    lv2_extra <- .make_parsed(
        features    = c("A_mz89.0477_rt0.45"),
        samples     = .samples,
        extra_annot = list(Adducts = "[M+H]+")
    )

    result <- merge_level_parsed(
        parsed_levels = list(Level_1 = lv1_plain, Level_2 = lv2_extra),
        level_names   = c("Level_1", "Level_2")
    )

    expect_true("Adducts" %in% colnames(result$row_data))
    expect_true(is.na(result$row_data$Adducts[result$row_data$Source_File == "Level_1"]))
})


# ==============================================================================
# 7. Sample mismatch error
# ==============================================================================

test_that("sample mismatch between levels raises an informative error", {
    # Replicate the validation logic from parse_multi_level() so we can test
    # it without needing a full config / parse_cd_raw call
    check_sample_consistency <- function(p_list, lv_names) {
        ref <- sort(p_list[[1]]$sample_ids)
        for (i in seq_along(p_list)[-1]) {
            lv <- sort(p_list[[i]]$sample_ids)
            if (!identical(ref, lv)) {
                stop(sprintf("sample mismatch between '%s' and '%s'",
                             lv_names[1], lv_names[i]))
            }
        }
    }

    lv_a <- .make_parsed(c("F1_mz100.0000_rt1.00"), c("S1", "S2"))
    lv_b <- .make_parsed(c("F2_mz200.0000_rt2.00"), c("S1", "S3"))  # S3 ≠ S2

    expect_error(
        check_sample_consistency(
            list(Level_A = lv_a, Level_B = lv_b),
            c("Level_A", "Level_B")
        ),
        regexp = "sample mismatch"
    )
})


# ==============================================================================
# 8. Single-level input
# ==============================================================================

test_that("single-level merge produces the same structure as multi-level", {
    result_single <- merge_level_parsed(
        parsed_levels = list(Level_1 = .lv1),
        level_names   = "Level_1"
    )
    result_multi <- merge_level_parsed(
        parsed_levels = list(Level_1 = .lv1, Level_2 = .lv2),
        level_names   = c("Level_1", "Level_2")
    )

    expect_named(result_single, names(result_multi))
    expect_equal(nrow(result_single$expr_raw), 2L)  # only .lv1 features
    expect_equal(colnames(result_single$row_data)[1:3],
                 colnames(result_multi$row_data)[1:3])
})


# ==============================================================================
# 9. sample_map handling
# ==============================================================================

test_that("sample_map is deduplicated when levels share identical mappings", {
    result <- merge_level_parsed(
        parsed_levels = list(Level_1 = .lv1, Level_2 = .lv2),
        level_names   = c("Level_1", "Level_2")
    )
    sm <- result$sample_map
    expect_false(is.null(sm))
    expect_equal(nrow(sm), length(.samples))   # one row per sample, no dupes
    expect_equal(nrow(sm), nrow(unique(sm)))
})

test_that("sample_map is NULL when all levels return NULL (processed_wide path)", {
    lv1 <- .lv1; lv1$sample_map <- NULL
    lv2 <- .lv2; lv2$sample_map <- NULL

    result <- merge_level_parsed(
        parsed_levels = list(Level_1 = lv1, Level_2 = lv2),
        level_names   = c("Level_1", "Level_2")
    )
    expect_null(result$sample_map)
})


# ==============================================================================
# 10. validate_metabolomics_config — multi_level guards (strict indexing)
# ==============================================================================

test_that("validate accepts a correct multi_level config", {
    cfg <- list(
        input = list(format = "multi_level", level_format = "cd_raw"),
        files = list(data_dir = "metabolomics/levels/", data = NULL),
        normalization = NULL
    )
    expect_true(validate_metabolomics_config(cfg))
})

test_that("validate rejects multi_level config missing data_dir", {
    cfg <- list(
        input = list(format = "multi_level", level_format = "cd_raw"),
        files = list(data_dir = NULL, data = NULL),
        normalization = NULL
    )
    expect_error(validate_metabolomics_config(cfg), regexp = "data_dir is required")
})

test_that("validate rejects multi_level config missing level_format", {
    cfg <- list(
        input = list(format = "multi_level", level_format = NULL),
        files = list(data_dir = "metabolomics/levels/", data = NULL),
        normalization = NULL
    )
    expect_error(validate_metabolomics_config(cfg), regexp = "input.level_format")
})

test_that("validate rejects level_format value that is not a per-file format", {
    cfg <- list(
        input = list(format = "multi_level", level_format = "long"),
        files = list(data_dir = "metabolomics/levels/", data = NULL),
        normalization = NULL
    )
    expect_error(validate_metabolomics_config(cfg))
})

test_that("validate rejects config with both data and data_dir set (mutual exclusivity)", {
    cfg <- list(
        input = list(format = "multi_level", level_format = "cd_raw"),
        files = list(data = "file.xlsx", data_dir = "metabolomics/levels/"),
        normalization = NULL
    )
    expect_error(validate_metabolomics_config(cfg), regexp = "mutually exclusive")
})

test_that("partial matching cannot cause a false positive: data_dir present, data absent", {
    # This is the core regression test for the strict-indexing fix.
    # Before the fix, cfg$files$data would partially match data_dir,
    # making has_data = TRUE and triggering the mutual exclusivity error.
    cfg <- list(
        input = list(format = "multi_level", level_format = "cd_raw"),
        files = list(data_dir = "metabolomics/levels/"),   # data key entirely absent
        normalization = NULL
    )
    expect_true(validate_metabolomics_config(cfg))
})

test_that("partial matching cannot cause a false positive: data present, data_dir absent", {
    # Symmetric case: cfg$files$data_dir must not partially match data.
    cfg <- list(
        input = list(format = "cd_raw"),
        files = list(data = "metabolomics/export.xlsx"),   # data_dir key entirely absent
        normalization = NULL
    )
    expect_true(validate_metabolomics_config(cfg))
})

test_that("validate still rejects unknown format values", {
    cfg <- list(
        input = list(format = "not_a_real_format"),
        files = list(data = "x.xlsx"),
        normalization = NULL
    )
    expect_error(validate_metabolomics_config(cfg))
})

test_that("validate accepts existing single-file formats unchanged", {
    cfg_cd <- list(
        input = list(format = "cd_raw"),
        files = list(data = "export.xlsx"),
        normalization = NULL
    )
    expect_true(validate_metabolomics_config(cfg_cd))

    cfg_pw <- list(
        input = list(format = "processed_wide"),
        files = list(data = "data.csv"),
        normalization = NULL
    )
    expect_true(validate_metabolomics_config(cfg_pw))
})


# ==============================================================================
# 11. Integration test — real XLSX files in a temp directory
#     (Only test that touches the filesystem; requires readxl + writexl)
# ==============================================================================

test_that("read_multi_level_dir reads real XLSX files from a temp directory", {
    skip_if_not_installed("writexl")
    skip_if_not_installed("readxl")

    tmp_dir <- tempfile("multi_level_test_")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

    df1 <- data.frame(
        Name        = c("Glucose", "Fructose"),
        `m/z`       = c(180.0634, 180.0634),
        `RT [min]`  = c(1.23, 2.10),
        check.names = FALSE
    )
    df2 <- data.frame(
        Name        = c("Alanine"),
        `m/z`       = c(89.0477),
        `RT [min]`  = c(0.45),
        check.names = FALSE
    )

    writexl::write_xlsx(df1, file.path(tmp_dir, "Level_1.xlsx"))
    writexl::write_xlsx(df2, file.path(tmp_dir, "Level_2.xlsx"))

    result <- read_multi_level_dir(tmp_dir, pattern = "\\.xlsx$")

    expect_length(result, 2L)
    expect_named(result, c("Level_1", "Level_2"))
    expect_true(is.data.frame(result[["Level_1"]]$data_df))
    expect_true(is.data.frame(result[["Level_2"]]$data_df))
    expect_equal(nrow(result[["Level_1"]]$data_df), 2L)
    expect_equal(nrow(result[["Level_2"]]$data_df), 1L)
    expect_equal(result[["Level_1"]]$level_name, "Level_1")
    expect_equal(result[["Level_2"]]$level_name, "Level_2")
})

test_that("read_multi_level_dir errors on an empty directory", {
    tmp_dir <- tempfile("empty_multi_level_")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

    expect_error(
        read_multi_level_dir(tmp_dir, pattern = "\\.xlsx$"),
        regexp = "No files matching"
    )
})
