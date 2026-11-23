# tests/testthat/test-enrichment-gsea-progress.R
#
# Live GSEA progress reporting and the per-job timing table. Workers announce
# each finished job with a `gsea_job_done` immediateCondition; a handler in the
# main process prints progress. These tests lock in the report cadence and
# content, the timing table, and that reporting leaves the results untouched.
# The GSEA call itself is stubbed (see local_stub_gsea_progress_worker), so no
# enrichment package is needed except where stated.

# A `gsea_job_done` condition as .with_gsea_job_timing() signals it.
make_gsea_done_cond <- function(i, status = "ok", duration_sec = 42) {
    structure(
        class = c("gsea_job_done", "immediateCondition", "condition"),
        list(message = "GSEA job done", call = NULL,
             contrast = paste0("c", i), db_name = "KEGG", ranking_method = "fc",
             duration_sec = duration_sec, status = status))
}

# Replace the GSEA worker factory where run_gsea_all() resolves it (globalenv),
# restored on test exit. `fail_contrast` makes that contrast's jobs return a
# caught error; every other job returns a successful empty result.
local_stub_gsea_progress_worker <- function(fail_contrast = character(0),
                                            env = parent.frame()) {
    target <- globalenv()
    orig <- get(".make_gsea_worker", envir = target)
    stub <- function(...) {
        worker <- function(job) list(
            ranking_method = job$ranking_method, contrast = job$contrast,
            db_name = job$db_name,
            gsea_result = if (job$contrast %in% fail_contrast) {
                structure(list(message = "stub failure"), class = "gsea_error")
            } else data.frame())
        # Keep the test frame out of the closure, so a multisession worker is
        # sent only what the stub uses.
        environment(worker) <- list2env(list(fail_contrast = fail_contrast),
                                        parent = globalenv())
        worker
    }
    assign(".make_gsea_worker", stub, envir = target)
    withr::defer(assign(".make_gsea_worker", orig, envir = target), envir = env)
}

# n_contrasts contrasts x 1 ranking x 1 database = n_contrasts jobs.
make_progress_fixture <- function(n_contrasts) {
    ranks <- setNames(c(2, 1), c("g1", "g2"))
    contrasts <- setNames(rep(list(ranks), n_contrasts), paste0("c", seq_len(n_contrasts)))
    list(
        ranked = list(fc = contrasts),
        lt = list(KEGG = list(
            TERM2GENE = data.frame(term = "p1", gene = c("g1", "g2"), stringsAsFactors = FALSE),
            TERM2NAME = data.frame(term = "p1", name = "P1", stringsAsFactors = FALSE))))
}

collect_messages <- function(expr) {
    msgs <- character(0)
    value <- withCallingHandlers(expr, message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
    })
    list(value = value, messages = msgs)
}

test_that(".format_gsea_duration picks a readable unit", {
    expect_identical(.format_gsea_duration(42.4), "42s")
    expect_identical(.format_gsea_duration(725), "12m 05s")
    expect_identical(.format_gsea_duration(7380), "2h 03m")
    expect_identical(.format_gsea_duration(NA_real_), "NA")
})

test_that("the progress handler reports every `every` jobs and on the last one", {
    ticks <- 0
    handler <- .make_gsea_progress_handler(
        total = 12L, every = 5L, clock = function() { ticks <<- ticks + 60; ticks })
    got <- collect_messages(for (i in 1:12) handler(make_gsea_done_cond(i)))$messages

    expect_length(got, 3L)
    expect_match(got[1], "[GSEA] 5/12 (41.7%) | elapsed: ", fixed = TRUE)
    expect_match(got[2], "[GSEA] 10/12 (83.3%)", fixed = TRUE)
    expect_match(got[3], "[GSEA] 12/12 (100.0%)", fixed = TRUE)
    # Each report names the job that just finished, with its own duration.
    expect_match(got[1], "Completed: c5 | KEGG | fc | 42s", fixed = TRUE)
    expect_match(got[3], "Completed: c12 | KEGG | fc | 42s", fixed = TRUE)
})

test_that("the progress handler counts a re-delivered job once and flags failures", {
    handler <- .make_gsea_progress_handler(total = 2L, every = 1L, clock = function() 0)
    got <- collect_messages({
        handler(make_gsea_done_cond(1))
        handler(make_gsea_done_cond(1))
        handler(make_gsea_done_cond(2, status = "failed"))
    })$messages

    expect_length(got, 2L)
    expect_match(got[1], "[GSEA] 1/2 ", fixed = TRUE)
    expect_false(grepl("FAILED", got[1], fixed = TRUE))
    expect_match(got[2], "[GSEA] 2/2 ", fixed = TRUE)
    expect_match(got[2], "| FAILED", fixed = TRUE)
})

test_that("the timing wrapper adds duration and status and leaves the result alone", {
    worker <- function(job) list(ranking_method = "fc", contrast = job$contrast,
                                 db_name = "KEGG", gsea_result = job$result)
    timed <- .with_gsea_job_timing(worker)

    signalled <- list()
    run <- function(job) withCallingHandlers(
        timed(job), gsea_job_done = function(cond) signalled[[length(signalled) + 1]] <<- cond)

    ok <- run(list(contrast = "c1", result = data.frame(x = 1)))
    expect_identical(ok[c("ranking_method", "contrast", "db_name", "gsea_result")],
                     worker(list(contrast = "c1", result = data.frame(x = 1))))
    expect_identical(ok$status, "ok")
    expect_true(is.numeric(ok$duration_sec) && ok$duration_sec >= 0)

    failed <- run(list(contrast = "c2",
                       result = structure(list(message = "boom"), class = "gsea_error")))
    expect_identical(failed$status, "failed")
    expect_identical(run(list(contrast = "c3", result = NULL))$status, "failed")

    expect_length(signalled, 3L)
    expect_s3_class(signalled[[1]], "immediateCondition")
    expect_identical(signalled[[2]]$contrast, "c2")
    expect_identical(signalled[[2]]$status, "failed")

    # With no handler registered the signal is a no-op.
    expect_no_error(timed(list(contrast = "c4", result = data.frame())))
})

test_that("the timing wrapper captures only the worker", {
    timed <- .with_gsea_job_timing(function(job) job)
    expect_identical(ls(environment(timed)), "worker")
})

test_that("run_gsea_all reports progress, then announces the writing stage", {
    local_stub_gsea_progress_worker()
    fx <- make_progress_fixture(7L)

    got <- collect_messages(
        run_gsea_all(fx$ranked, fx$lt, output_dir = NULL, workers = 1))$messages

    progress <- grep("[GSEA] ", got, fixed = TRUE)
    expect_length(progress, 2L)
    expect_match(got[progress[1]], "[GSEA] 5/7 (71.4%)", fixed = TRUE)
    expect_match(got[progress[2]], "[GSEA] 7/7 (100.0%)", fixed = TRUE)
    expect_match(got[progress[2]], "Completed: c7 | KEGG | fc | ", fixed = TRUE)

    complete <- grep("GSEA computations complete: 7 jobs in ", got, fixed = TRUE)
    expect_length(complete, 1L)
    expect_match(got[complete], "writing results and generating plots", fixed = TRUE)
    # Progress comes before the hand-over line, per-unit lines after it.
    per_unit <- grep("no gene set within the size limits", got, fixed = TRUE)
    expect_length(per_unit, 7L)
    expect_true(max(progress) < complete)
    expect_true(complete < min(per_unit))
})

test_that("run_gsea_all reports progress through multisession workers", {
    skip_if_not_installed("future")
    skip_if_not_installed("future.apply")
    local_stub_gsea_progress_worker()
    fx <- make_progress_fixture(6L)

    got <- collect_messages(
        run_gsea_all(fx$ranked, fx$lt, output_dir = NULL, workers = 2))$messages

    progress <- grep("[GSEA] ", got, fixed = TRUE, value = TRUE)
    expect_length(progress, 2L)
    expect_match(progress[1], "[GSEA] 5/6 ", fixed = TRUE)
    expect_match(progress[2], "[GSEA] 6/6 (100.0%)", fixed = TRUE)
})

test_that("run_gsea_all writes one timing row per job, with failures marked", {
    local_stub_gsea_progress_worker(fail_contrast = "c2")
    fx <- make_progress_fixture(3L)
    out_dir <- withr::local_tempdir()

    got <- collect_messages(
        run_gsea_all(fx$ranked, fx$lt, output_dir = out_dir, workers = 1))

    timing_file <- file.path(out_dir, "gsea_job_timings.tsv")
    expect_true(file.exists(timing_file))
    timings <- utils::read.delim(timing_file, stringsAsFactors = FALSE)
    expect_identical(colnames(timings),
                     c("contrast", "database", "ranking_method", "duration_sec", "status"))
    expect_identical(timings$contrast, c("c1", "c2", "c3"))
    expect_identical(timings$database, rep("KEGG", 3L))
    expect_identical(timings$ranking_method, rep("fc", 3L))
    expect_identical(timings$status, c("ok", "failed", "ok"))
    expect_true(is.numeric(timings$duration_sec) && all(timings$duration_sec >= 0))

    expect_true(any(grepl("3 jobs in .* \\(1 failed\\)", got$messages)))
    # The manifest still records the failure exactly as before.
    expect_identical(got$value$manifest$status, c("empty", "failed", "empty"))
})

test_that("timings stay out of run_gsea_all's return value", {
    local_stub_gsea_progress_worker()
    fx <- make_progress_fixture(2L)
    out <- suppressMessages(run_gsea_all(fx$ranked, fx$lt, output_dir = NULL, workers = 1))
    expect_named(out, c("results", "plot_files", "manifest", "index"), ignore.order = TRUE)
})

test_that("progress reporting does not change GSEA results", {
    skip_if_not_installed("fgsea")
    ranked <- setNames(c(seq(10, 3, length.out = 8), seq(2, -8, length.out = 22)),
                       paste0("g", 1:30))
    lt <- list(KEGG = list(
        TERM2GENE = data.frame(term = "p1", gene = paste0("g", 1:8), stringsAsFactors = FALSE),
        TERM2NAME = data.frame(term = "p1", name = "P1", stringsAsFactors = FALSE)))
    jobs <- list(list(ranking_method = "fc", contrast = "c1", db_name = "KEGG",
                      ranked = ranked,
                      gene_sets = term2gene_to_gene_sets(lt$KEGG$TERM2GENE, lt$KEGG$TERM2NAME),
                      max_size = 8L))
    worker <- .make_gsea_worker(0.05, "fdr", engine = "fgsea", min_size = 2, seed = 111L)

    plain <- suppressWarnings(run_enrichment_jobs(jobs, worker, workers = 1L, seed = 111L))
    timed <- suppressWarnings(suppressMessages(withCallingHandlers(
        run_enrichment_jobs(jobs, .with_gsea_job_timing(worker), workers = 1L, seed = 111L),
        gsea_job_done = .make_gsea_progress_handler(1L))))

    expect_false(inherits(plain[[1]]$gsea_result, "gsea_error"))
    expect_identical(timed[[1]]$gsea_result, plain[[1]]$gsea_result)
    expect_identical(timed[[1]]$status, "ok")
})
