# tests/testthat/test-execution-info-cue.R
#
# execution_info_files writes git_commit.txt and timestamp.txt into the run
# folder. Its inputs are the config and paths only, so without a cue a
# code-only rerun skipped it and the folder kept naming the previous run's
# commit. These tests read _targets.R as code (not text) and pin the cue.

# Every tar_target(<name>, ...) call in _targets.R, found by walking the parsed
# expressions, so formatting and comments cannot fool the check.
find_tar_target <- function(name) {
    exprs <- parse(file.path(root_dir, "_targets.R"), keep.source = FALSE)
    hits <- list()
    walk <- function(e) {
        if (is.call(e)) {
            fn <- e[[1]]
            if (is.symbol(fn) && identical(as.character(fn), "tar_target") &&
                length(e) >= 2 && is.symbol(e[[2]]) &&
                identical(as.character(e[[2]]), name)) {
                hits[[length(hits) + 1L]] <<- e
            }
            # Only calls can hold a tar_target(); checking is.call() also steps
            # over empty arguments such as x[, 1] without evaluating them.
            args <- as.list(e)[-1]
            for (i in seq_along(args)) if (is.call(args[[i]])) walk(args[[i]])
        }
    }
    for (e in exprs) walk(e)
    hits
}

test_that("execution_info_files is defined exactly once", {
    expect_length(find_tar_target("execution_info_files"), 1L)
})

test_that("execution_info_files reruns on every tar_make()", {
    call <- find_tar_target("execution_info_files")[[1]]
    cue <- as.list(call)$cue
    expect_false(is.null(cue))
    expect_identical(as.character(cue[[1]]), "tar_cue")
    expect_identical(as.list(cue)$mode, "always")
})
