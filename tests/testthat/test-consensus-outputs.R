# tests/testthat/test-consensus-outputs.R
#
# The report reads the consensus files by name and shows them as the current
# run's, so every path that does not produce a complete set clears the last
# run's first: a disabled step, too few methods, an error part-way through, and
# the pipeline skipping the module altogether. Stability writes beside it and
# must survive all of these.

# Stub by assignment into the function's own environment: these are sourced
# functions, not a package, so with_mocked_bindings() has nothing to rebind.
# As in test-loadings-gsea.R.
local_stubs <- function(stubs, env = parent.frame()) {
    target <- environment(mod_multiomics_consensus)
    nms <- names(stubs)
    had <- vapply(nms, exists, logical(1), envir = target, inherits = FALSE)
    old <- lapply(nms[had], get, envir = target, inherits = FALSE)
    names(old) <- nms[had]
    withr::defer({
        for (nm in nms) {
            if (nm %in% names(old)) {
                assign(nm, old[[nm]], envir = target)
            } else if (exists(nm, envir = target, inherits = FALSE)) {
                rm(list = nm, envir = target)
            }
        }
    }, envir = env)
    for (nm in nms) assign(nm, stubs[[nm]], envir = target)
    invisible(NULL)
}

# An outer consensus/ directory as an earlier run left it: the inner consensus
# outputs the report reads, and the stability outputs beside them.
stale_consensus <- function(envir = parent.frame()) {
    outer <- withr::local_tempdir(.local_envir = envir)
    inner <- file.path(outer, "consensus")
    old <- file.path(inner, c("sample_clusters_comparison.csv", "method_ari_matrix.csv",
                              "method_nmi_matrix.csv", file.path("plots", "meta_pca.png")))
    kept <- file.path(outer, "stability", c("feature_stability.csv", "stability.png"))
    for (f in c(old, kept)) {
        dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
        writeLines("x", f)
    }
    list(outer = outer, inner = inner, old = old, kept = kept)
}

no_stability <- function(consensus = list()) {
    list(modes = list(multiomics = list(consensus = consensus,
                                        stability = list(run_stability = FALSE))))
}

test_that("clear_consensus_outputs removes the inner directory only", {
    d <- stale_consensus()
    expect_message(clear_consensus_outputs(d$inner), "cleared")
    expect_false(dir.exists(d$inner))
    expect_true(all(file.exists(d$kept)))
    # Handed the outer directory, it refuses rather than take stability with it.
    expect_error(clear_consensus_outputs(d$outer), "stability")
    expect_true(all(file.exists(d$kept)))
})

test_that("clear_consensus_outputs is quiet on a fresh or absent directory", {
    expect_length(clear_consensus_outputs(NULL), 0L)
    expect_length(clear_consensus_outputs(file.path(withr::local_tempdir(), "none")), 0L)
})

test_that("a consensus run that compares nothing still clears the last run's", {
    # Fewer than two methods: the early return comes after the cleanup.
    one <- stale_consensus()
    expect_null(suppressMessages(run_integration_consensus(
        list(snf = list(clusters = 1), diablo = NULL, mofa = NULL),
        mae = NULL, config = list(), out_dir = one$inner)))
    expect_false(any(file.exists(one$old)))
    expect_true(all(file.exists(one$kept)))

    # Method comparison switched off inside the step.
    off <- stale_consensus()
    cfg <- list(modes = list(multiomics = list(consensus = list(compare_methods = FALSE))))
    expect_null(suppressMessages(run_integration_consensus(
        list(snf = list(), diablo = list()), mae = NULL, config = cfg, out_dir = off$inner)))
    expect_false(any(file.exists(off$old)))
    expect_true(all(file.exists(off$kept)))
})

test_that("the consensus module clears when disabled and when the step fails", {
    disabled <- stale_consensus()
    res <- suppressMessages(mod_multiomics_consensus(
        harmonization_res = list(), integration_res = list(),
        config = no_stability(list(run_consensus = FALSE)), out_dir = disabled$outer))
    expect_null(res$consensus_results)
    expect_false(any(file.exists(disabled$old)))
    expect_true(all(file.exists(disabled$kept)))

    # A failure after some files were written leaves none of them behind.
    failed <- stale_consensus()
    local_stubs(list(run_integration_consensus = function(integration_results, mae,
                                                          config, out_dir) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        writeLines("partial", file.path(out_dir, "method_ari_matrix.csv"))
        stop("boom")
    }))
    res <- suppressWarnings(suppressMessages(mod_multiomics_consensus(
        harmonization_res = list(), integration_res = list(),
        config = no_stability(), out_dir = failed$outer)))
    expect_null(res$consensus_results)
    expect_false(dir.exists(failed$inner))
    expect_true(all(file.exists(failed$kept)))
})

test_that("the pipeline clears consensus outputs when it never calls the module", {
    src <- paste(readLines(testthat::test_path(
        "..", "..", "R", "pipeline", "multiomics", "00_pipe_multiomics.R")),
        collapse = "\n")
    block <- regmatches(src, regexpr(
        "(?s)multiomics_consensus,.*?\n        \\),", src, perl = TRUE))
    expect_length(block, 1)
    skip <- regmatches(block, regexpr(
        "(?s)Skipping consensus analysis.*?return\\(NULL\\)", block, perl = TRUE))
    expect_match(skip, 'clear_consensus_outputs(file.path(multiomics_out_dir,', fixed = TRUE)
    expect_match(skip, '"consensus", "consensus"))', fixed = TRUE)
    # And that is the directory the module writes its consensus into.
    expect_match(block, 'out_dir = file.path(multiomics_out_dir, "consensus")', fixed = TRUE)
})
