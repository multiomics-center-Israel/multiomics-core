test_that("resolve_group_col walks the canonical chain", {
    expect_equal(
        resolve_group_col(list(filtering = list(group_col = "treatment"),
                               de_table  = list(group_col = "ignored"),
                               effects   = list(color = "ignored2"))),
        "treatment"
    )
    expect_equal(
        resolve_group_col(list(de_table = list(group_col = "Condition"),
                               effects  = list(color = "ignored"))),
        "Condition"
    )
    expect_equal(
        resolve_group_col(list(effects = list(color = c("Group", "batch")))),
        "Group"
    )
    expect_equal(resolve_group_col(list()), "Condition")
})

test_that("resolve_group_col is strict — empty/NA fall through", {
    cfg <- list(filtering = list(group_col = ""),
                de_table  = list(group_col = NA_character_),
                effects   = list(color = "treatment"))
    expect_equal(resolve_group_col(cfg), "treatment")
})

test_that("resolve_group_col uses only the first element of effects$color", {
    cfg <- list(effects = list(color = c("primary", "secondary", "tertiary")))
    expect_equal(resolve_group_col(cfg), "primary")
})

test_that("assert_design_viable stops on group with n=1", {
    meta <- data.frame(
        sample    = paste0("S", 1:5),
        Condition = c("ctrl", "ctrl", "ctrl", "trt", "trt")
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    meta1 <- data.frame(sample = "S1", Condition = "ctrl")
    expect_error(assert_design_viable(meta1, cfg), "n >= 2 per group")
})

test_that("assert_design_viable stops naming all underpowered groups when multiple", {
    meta <- data.frame(
        sample    = paste0("S", 1:2),
        Condition = c("ctrl", "trt")
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_error(assert_design_viable(meta, cfg), "Groups with n < 2")
})

test_that("assert_design_viable warns severe at n=2 (single composite warning)", {
    meta <- data.frame(
        sample    = paste0("S", 1:4),
        Condition = c("ctrl", "ctrl", "trt", "trt")
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_warning(assert_design_viable(meta, cfg), "n=2: severe")
})

test_that("assert_design_viable warns low-power at 3 <= n < 5", {
    meta <- data.frame(
        sample    = paste0("S", 1:6),
        Condition = rep(c("ctrl", "trt"), each = 3)
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_warning(assert_design_viable(meta, cfg), "low statistical power")
})

test_that("assert_design_viable warns on global imbalance separately", {
    meta <- data.frame(
        sample    = paste0("S", 1:23),
        Condition = c(rep("ctrl", 20), rep("trt", 3))
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_warning(assert_design_viable(meta, cfg), "imbalanced")
})

test_that("assert_design_viable is silent at n=5 per group (boundary)", {
    meta <- data.frame(
        sample    = paste0("S", 1:10),
        Condition = rep(c("ctrl", "trt"), each = 5)
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_silent(assert_design_viable(meta, cfg))
})

test_that("assert_design_viable is silent at n=5/5/5 three-group design", {
    meta <- data.frame(
        sample    = paste0("S", 1:15),
        Condition = rep(c("a", "b", "c"), each = 5)
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_silent(assert_design_viable(meta, cfg))
})

test_that("assert_design_viable stops when group column is missing from meta", {
    meta <- data.frame(sample = paste0("S", 1:6),
                       treatment = rep(c("ctrl", "trt"), each = 3))
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_error(
        assert_design_viable(meta, cfg),
        "not found in metadata"
    )
})

test_that("assert_design_viable stops on NA values in the group column", {
    meta <- data.frame(
        sample    = paste0("S", 1:6),
        Condition = c("ctrl", "ctrl", NA, "trt", "trt", "trt"),
        stringsAsFactors = FALSE
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_error(assert_design_viable(meta, cfg), "Assign a group label")
})

test_that("assert_design_viable stops on empty-string values in the group column", {
    meta <- data.frame(
        sample    = paste0("S", 1:6),
        Condition = c("ctrl", "ctrl", "", "trt", "trt", "trt"),
        stringsAsFactors = FALSE
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_error(assert_design_viable(meta, cfg), "Assign a group label")
})
