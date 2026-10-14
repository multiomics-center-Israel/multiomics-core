# Tests for the retired-config-key warning in validate_config()
# (R/core/04_config.R: .retired_config_keys, find_retired_config_keys(),
# warn_retired_config_keys()).

retired_key_warnings <- function(expr) {
    msgs <- character(0)
    withCallingHandlers(expr, warning = function(w) {
        msgs <<- c(msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    grep("no longer read", msgs, value = TRUE)
}

minimal_config <- function() {
    list(paths = list(raw = "data", out = "outputs"), params = list(seed = 1),
         modes = list())
}

test_that("a config with retired keys warns once, listing each key and its reason", {
    cfg <- minimal_config()
    cfg$output <- list(output_dir = "outputs", report_format = "html")
    cfg$design <- list(condition_column = "condition", batch_column = NULL)
    cfg$modes$lipidomics <- list(commentary = list(retry_delay = 2))

    w <- retired_key_warnings(validate_config(cfg))
    expect_length(w, 1)
    for (k in c("output.output_dir", "output.report_format", "design.batch_column",
                "modes.lipidomics.commentary.retry_delay")) {
        expect_match(w, paste0(k, " is no longer read ("), fixed = TRUE)
    }
    expect_match(w, "outputs are written under paths.out", fixed = TRUE)
    expect_match(w, "; remove it", fixed = TRUE)
    expect_false(grepl("condition_column", w))
})

test_that("a key set to null still counts as present", {
    cfg <- minimal_config()
    cfg$design <- list(batch_column = NULL)
    expect_equal(find_retired_config_keys(cfg)$key, "design.batch_column")
})

test_that("a clean config does not warn", {
    expect_length(retired_key_warnings(validate_config(minimal_config())), 0)
    expect_equal(nrow(find_retired_config_keys(minimal_config())), 0)
})

test_that("a retired key is never an error", {
    cfg <- minimal_config()
    cfg$output <- list(report_format = "html")
    expect_no_error(suppressWarnings(validate_config(cfg)))
})

test_that("the retired-key list is well formed", {
    rk <- .retired_config_keys
    expect_named(rk, c("key", "reason", "since"))
    expect_false(anyDuplicated(rk$key) > 0)
    expect_true(all(nzchar(rk$reason)) && all(nzchar(rk$since)))
})

test_that("the shipped templates validate without retired-key warnings", {
    tpl_dir <- file.path(root_dir, "config", "templates")
    for (f in list.files(tpl_dir, pattern = "\\.ya?ml$", full.names = TRUE)) {
        cfg <- yaml::read_yaml(f)
        expect_equal(nrow(find_retired_config_keys(cfg)), 0, info = basename(f))
        w <- tryCatch(retired_key_warnings(validate_config(cfg)),
                      error = function(e) paste("error:", conditionMessage(e)))
        expect_identical(w, character(0), info = basename(f))
    }
})
