# Tests for input path resolution
# test-paths.R

test_that("resolve_input_path resolves relative paths under raw/", {
    config <- list(project = list(dir = "/proj"), paths = list(raw = "data"))
    expect_equal(resolve_input_path(config, "metabolics/kegg.gmt"),
                 file.path("/proj", "data", "metabolics/kegg.gmt"))
})

test_that("resolve_input_path defaults the raw directory to 'data'", {
    config <- list(project = list(dir = "/proj"))
    expect_equal(resolve_input_path(config, "x.gmt"),
                 file.path("/proj", "data", "x.gmt"))
})

test_that("resolve_input_path leaves absolute paths unchanged", {
    config <- list(project = list(dir = "/proj"), paths = list(raw = "data"))
    expect_equal(resolve_input_path(config, "/abs/x.gmt"), "/abs/x.gmt")
    expect_equal(resolve_input_path(config, "C:/Users/me/x.gmt"), "C:/Users/me/x.gmt")
    expect_equal(resolve_input_path(config, "C:\\Users\\me\\x.gmt"), "C:\\Users\\me\\x.gmt")
    expect_equal(resolve_input_path(config, "~/x.gmt"), "~/x.gmt")
    expect_equal(resolve_input_path(config, "\\\\host\\share\\x.gmt"),
                 "\\\\host\\share\\x.gmt")
})

test_that("resolve_input_path handles NULL and mixed vectors", {
    config <- list(project = list(dir = "/proj"), paths = list(raw = "data"))
    expect_null(resolve_input_path(config, NULL))
    res <- resolve_input_path(config, c("/abs/a.gmt", "rel/b.gmt"))
    expect_equal(res, c("/abs/a.gmt", file.path("/proj", "data", "rel/b.gmt")))
})
