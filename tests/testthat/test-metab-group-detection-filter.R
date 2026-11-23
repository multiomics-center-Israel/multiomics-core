# tests/testthat/test-metab-group-detection-filter.R
#
# preprocessing.min_detected_per_group keeps a metabolite only when, in at
# least one group, at least min_samples samples are detected (non-NA and non-zero).
# All data are synthetic.

make_gd_fixture <- function() {
  ids  <- sprintf("S%d", 1:8)
  meta <- data.frame(sample = ids,
                     group  = rep(c("A", "B"), each = 4),
                     stringsAsFactors = FALSE)
  mat <- rbind(
    all_on    = c(5, 5, 5, 5,  5, 5, 5, 5),
    a_three   = c(5, 5, 5, 0,  0, 0, 0, 0),   # 3/4 in A -> keep
    b_three   = c(0, 0, 0, 0,  5, NA, 5, 5),  # 3/4 in B -> keep
    split_two = c(5, 5, 0, 0,  5, 5, 0, NA),  # 2/4 in each -> drop
    zeros     = c(0, 0, 0, 0,  0, 0, 0, 0),   # zeros are not detections
    nas       = rep(NA_real_, 8)
  )
  colnames(mat) <- ids
  list(mat = mat, meta = meta)
}

test_that("keeps a feature detected in >= n samples of at least one group", {
  fx  <- make_gd_fixture()
  out <- filter_by_group_detection(fx$mat, fx$meta, "sample", "group", min_n = 3)
  expect_setequal(rownames(out$mat), c("all_on", "a_three", "b_three"))
  expect_setequal(out$dropped_features, c("split_two", "zeros", "nas"))
  expect_identical(colnames(out$mat), colnames(fx$mat))
})

test_that("zeros and NA both count as not detected", {
  fx  <- make_gd_fixture()
  out <- filter_by_group_detection(fx$mat, fx$meta, "sample", "group", min_n = 1)
  expect_false("zeros" %in% rownames(out$mat))
  expect_false("nas" %in% rownames(out$mat))
  expect_true("split_two" %in% rownames(out$mat))
})

test_that("ignored groups do not rescue a feature", {
  fx <- make_gd_fixture()
  fx$meta$group[5:8] <- "Pool"
  out <- filter_by_group_detection(fx$mat, fx$meta, "sample", "group",
                                   min_n = 3, ignore_groups = "Pool")
  expect_false("b_three" %in% rownames(out$mat))
  expect_true("a_three" %in% rownames(out$mat))
})

test_that("works with a single feature row", {
  fx  <- make_gd_fixture()
  out <- filter_by_group_detection(fx$mat["a_three", , drop = FALSE], fx$meta,
                                   "sample", "group", min_n = 3)
  expect_identical(rownames(out$mat), "a_three")
})

test_that("a missing group column stops with the column list", {
  fx <- make_gd_fixture()
  expect_error(
    filter_by_group_detection(fx$mat, fx$meta, "sample", "condition", min_n = 3),
    "no such column")
})

test_that("resolver: absent block means off, a valid block resolves", {
  expect_null(resolve_group_detection_cfg(list(feat_missing_threshold = 0.5)))
  expect_null(resolve_group_detection_cfg(NULL))
  gd <- resolve_group_detection_cfg(list(
    min_detected_per_group = list(min_samples = 3, group_col = "group",
                                  ignore_groups = c("Pool"))))
  expect_identical(gd$min_samples, 3L)
  expect_identical(gd$group_col, "group")
  expect_identical(gd$ignore_groups, "Pool")
  gd2 <- resolve_group_detection_cfg(list(
    min_detected_per_group = list(min_samples = 2, group_col = "group")))
  expect_identical(gd2$ignore_groups, character(0))
})

test_that("resolver rejects bad settings", {
  bad <- function(gd) resolve_group_detection_cfg(list(min_detected_per_group = gd))
  expect_error(bad(3), "block of settings")
  expect_error(bad(list(min_samples = 0, group_col = "group")))
  expect_error(bad(list(min_samples = 2.5, group_col = "group")), "whole number")
  expect_error(bad(list(min_samples = 3)))
  expect_error(bad(list(min_samples = 3, group_col = "group", min_n = 3)), "unknown key")
})

test_that("validate_metabolomics_config checks the block", {
  base <- list(input = list(format = "processed_wide"),
               files = list(data = "x.xlsx"),
               preprocessing = list(min_detected_per_group =
                                      list(min_samples = 3, group_col = "group")))
  expect_true(validate_metabolomics_config(base))
  base$preprocessing$min_detected_per_group$min_samples <- "three"
  expect_error(validate_metabolomics_config(base))
})

test_that("mod_met_filtered applies the filter only when configured", {
  fx  <- make_gd_fixture()
  raw <- list(expr_raw = fx$mat, meta = fx$meta, sample_col = "sample",
              row_data = data.frame(id = rownames(fx$mat),
                                    row.names = rownames(fx$mat)))
  cfg_off <- list(modes = list(metabolomics = list(preprocessing = list(
    feat_missing_threshold = 1.0))))
  cfg_on  <- cfg_off
  cfg_on$modes$metabolomics$preprocessing$min_detected_per_group <-
    list(min_samples = 3, group_col = "group")

  off <- suppressMessages(mod_met_filtered(raw, cfg_off))
  on  <- suppressMessages(mod_met_filtered(raw, cfg_on))
  expect_identical(nrow(off$mat), nrow(fx$mat))
  expect_setequal(rownames(on$mat), c("all_on", "a_three", "b_three"))
  expect_identical(rownames(on$row_data), rownames(on$mat))
  expect_true(all(c("split_two", "zeros") %in% on$dropped_features))
})

test_that("the block survives a YAML round trip", {
  # yaml reads a bare `n:` key as FALSE, which is why the key is min_samples.
  txt <- paste("min_detected_per_group:",
               "  min_samples: 3",
               "  group_col: \"group\"",
               "  ignore_groups: [\"Pool\"]", sep = "\n")
  gd <- resolve_group_detection_cfg(yaml::yaml.load(txt))
  expect_identical(gd$min_samples, 3L)
  expect_identical(gd$ignore_groups, "Pool")
})
