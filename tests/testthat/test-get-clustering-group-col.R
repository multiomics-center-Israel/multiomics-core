# tests/testthat/test-get-clustering-group-col.R
#
# Unit tests for get_clustering_group_col() — strict config validation

test_that("get_clustering_group_col returns valid column", {
  cfg <- list(clustering = list(group_col = "treatment"))
  meta <- data.frame(sample_id = c("s1", "s2"), treatment = c("A", "B"))

  expect_equal(get_clustering_group_col(cfg, meta), "treatment")
})

test_that("get_clustering_group_col errors on missing group_col", {
  cfg <- list(clustering = list(enabled = TRUE))
  meta <- data.frame(sample_id = c("s1", "s2"), treatment = c("A", "B"))

  expect_error(
    get_clustering_group_col(cfg, meta),
    "clustering\\$group_col is required but missing or empty"
  )
})

test_that("get_clustering_group_col errors on empty string", {
  cfg <- list(clustering = list(group_col = ""))
  meta <- data.frame(sample_id = c("s1", "s2"), treatment = c("A", "B"))

  expect_error(
    get_clustering_group_col(cfg, meta),
    "clustering\\$group_col is required but missing or empty"
  )
})

test_that("get_clustering_group_col errors when column not in meta", {
  cfg <- list(clustering = list(group_col = "nonexistent"))
  meta <- data.frame(sample_id = c("s1", "s2"), treatment = c("A", "B"))

  expect_error(
    get_clustering_group_col(cfg, meta),
    "not found in metadata columns"
  )
})

test_that("get_n_groups_from_effects returns 0 when group_col missing", {
  cfg <- list(clustering = list(enabled = TRUE))
  pre <- list(meta = data.frame(sample_id = "s1", treatment = "A"))

  expect_equal(get_n_groups_from_effects(pre, cfg), 0L)
})

test_that("get_n_groups_from_effects counts distinct groups", {
  cfg <- list(clustering = list(group_col = "treatment"))
  pre <- list(meta = data.frame(
    sample_id = paste0("s", 1:6),
    treatment = c("A", "A", "B", "B", "C", "C")
  ))

  expect_equal(get_n_groups_from_effects(pre, cfg), 3L)
})
