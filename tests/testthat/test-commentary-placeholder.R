# tests/testthat/test-commentary-placeholder.R
#
# Regression for #142: R/domain/multiomics/13_commentary.R used to define its
# own create_placeholder_commentary(figure_id, error_message). The pipeline
# sources domain/ after services/, so that copy replaced the services one and
# every services fallback, which passes `reason =`, failed with "unused
# argument". helper.R sources R/ alphabetically (services last), so the rest
# of the suite never saw the shadowing copy. These tests load the files in the
# order _targets.R uses instead.

# Files under R/ in the order _targets.R sources them.
pipeline_source_order <- function() {
    layer_files <- function(layer) {
        sort(list.files(file.path(root_dir, "R", layer), pattern = "\\.R$",
                        full.names = TRUE, recursive = TRUE))
    }
    unlist(lapply(c("core", "services", "domain", "modules", "pipeline"), layer_files))
}

# Files that assign `name` at top level, in pipeline source order.
files_defining <- function(name) {
    defines <- function(f) {
        any(vapply(parse(f, keep.source = FALSE), function(e) {
            is.call(e) && is.name(e[[1]]) && as.character(e[[1]]) %in% c("<-", "=") &&
                identical(e[[2]], as.name(name))
        }, logical(1)))
    }
    Filter(defines, pipeline_source_order())
}

test_that("the layer order mirrored here is the one _targets.R uses", {
    targets_src <- readLines(file.path(root_dir, "_targets.R"), warn = FALSE)
    expect_true(any(grepl("c(core_files, service_files, domain_files, module_files)",
                          targets_src, fixed = TRUE)))
})

test_that("a missing figure gets a placeholder when files load in pipeline order", {
    env <- new.env(parent = globalenv())
    for (f in files_defining("create_placeholder_commentary")) {
        source(f, local = env)
    }

    out_dir <- withr::local_tempdir()
    figures_tbl <- data.frame(figure_id = "fig_1",
                              filepath = file.path(out_dir, "missing.png"))

    res <- suppressMessages(
        env$generate_all_commentary(figures_tbl, config = list(), output_dir = out_dir)
    )

    expect_named(res, "fig_1")
    expect_equal(res$fig_1$backend, "placeholder")
    expect_match(res$fig_1$what_is_this, "Figure file not found", fixed = TRUE)
})

test_that("create_placeholder_commentary has exactly one definition under R/", {
    expect_equal(basename(files_defining("create_placeholder_commentary")),
                 "12_commentary.R")
})
