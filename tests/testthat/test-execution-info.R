# tests/testthat/test-execution-info.R
#
# execution_info_files writes git_commit.txt and timestamp.txt into the run
# folder. Its inputs are the config and paths only, so without a cue a
# code-only rerun skipped it and the folder kept naming the previous run's
# commit. The first tests read _targets.R as code (not text) and pin the cue;
# the rest exercise write_execution_info() and the run-status record on a
# synthetic run folder.

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

# A throwaway working dir holding a fake R/pipeline tree, because
# write_execution_info() snapshots R/pipeline relative to the working dir.
local_exec_fixture <- function(pipeline_files = c("00_pipe_a.R", "00_pipe_b.R"),
                               env = parent.frame()) {
    wd <- withr::local_tempdir(.local_envir = env)
    dir.create(file.path(wd, "R", "pipeline", "mode"), recursive = TRUE)
    for (f in pipeline_files) writeLines("# synthetic", file.path(wd, "R", "pipeline", "mode", f))
    writeLines("# synthetic targets", file.path(wd, "_targets.R"))
    writeLines("project: {name: demo}", file.path(wd, "config.yaml"))
    withr::local_dir(wd, .local_envir = env)
    wd
}

# The fixture is not a git repo, so `git rev-parse` warns; that path is what
# the stale-commit test exercises, and the warning itself is not under test.
write_fixture_info <- function(wd) {
    suppressWarnings(write_execution_info(config = list(project = list(name = "demo")),
                         run_dir = file.path(wd, "run"),
                         config_path = file.path(wd, "config.yaml"),
                         targets_file = "_targets.R"))
}

test_that("every written file is returned, timestamp and sessionInfo included", {
    wd <- local_exec_fixture()
    files <- write_fixture_info(wd)
    exec <- file.path(wd, "run", "execution_info")
    expect_true(all(file.exists(files)))
    expect_true(all(c("timestamp.txt", "sessionInfo.txt", "config_used.yaml",
                      "_targets_used.R") %in% basename(files)))
    # Every tracked file the folder holds is in the return value, so a change
    # to any of them reaches dependents such as rna_fact_sheet.
    on_disk <- list.files(exec, recursive = TRUE, full.names = TRUE)
    on_disk <- setdiff(normalizePath(on_disk),
                       normalizePath(file.path(exec, "run_status.txt")))
    expect_setequal(on_disk, normalizePath(files))
})

test_that("a rerun drops snapshots of pipeline files that no longer exist", {
    wd <- local_exec_fixture()
    write_fixture_info(wd)
    snap <- file.path(wd, "run", "execution_info", "pipeline_files_used")
    expect_true(file.exists(file.path(snap, "00_pipe_b.R")))

    file.remove(file.path(wd, "R", "pipeline", "mode", "00_pipe_b.R"))
    files <- write_fixture_info(wd)
    expect_false(file.exists(file.path(snap, "00_pipe_b.R")))
    expect_true(file.exists(file.path(snap, "00_pipe_a.R")))
    expect_false(any(basename(files) == "00_pipe_b.R"))
})

test_that("a stale git_commit.txt is removed when no commit can be read", {
    wd <- local_exec_fixture()   # a temp dir outside any git repo
    exec <- file.path(wd, "run", "execution_info")
    dir.create(exec, recursive = TRUE)
    writeLines("0000000000000000000000000000000000000000", file.path(exec, "git_commit.txt"))
    files <- write_fixture_info(wd)
    expect_false(file.exists(file.path(exec, "git_commit.txt")))
    expect_false("git_commit.txt" %in% basename(files))
})

test_that("files other writers put in execution_info are left alone", {
    wd <- local_exec_fixture()
    exec <- file.path(wd, "run", "execution_info")
    dir.create(exec, recursive = TRUE)
    writeLines("kept", file.path(exec, "other_writer.txt"))
    write_fixture_info(wd)
    expect_identical(readLines(file.path(exec, "other_writer.txt")), "kept")
})

test_that("run status starts as started and moves only from started", {
    wd <- local_exec_fixture()
    run <- file.path(wd, "run")
    write_fixture_info(wd)
    expect_identical(read_run_status(run), "started")

    expect_true(write_run_status(run, "completed", only_from = "started"))
    expect_identical(read_run_status(run), "completed")

    # A later failure must not relabel a record its own run did not write.
    expect_false(write_run_status(run, "failed", only_from = "started"))
    expect_identical(read_run_status(run), "completed")

    # A new run resets the record to started.
    write_fixture_info(wd)
    expect_identical(read_run_status(run), "started")
})

test_that("read_run_status returns NA when nothing usable is recorded", {
    d <- withr::local_tempdir()
    expect_identical(read_run_status(d), NA_character_)
    dir.create(file.path(d, "execution_info"))
    writeLines("not a status", file.path(d, "execution_info", "run_status.txt"))
    expect_identical(read_run_status(d), NA_character_)
    expect_error(write_run_status(d, "finished"))
})
