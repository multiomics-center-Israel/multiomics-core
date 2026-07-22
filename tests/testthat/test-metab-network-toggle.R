# tests/testthat/test-metab-network-toggle.R
#
# mod_metabolomics_network() is opt-in via enrichment$run_network, mirroring
# enrichment$run_enrichment: it is OFF by default. When unset or false the step
# short-circuits with a NULL return before any DE/KEGG work; only an explicit
# true lets it proceed.

test_that("network step is skipped by default (run_network unset)", {
    cfg <- list(modes = list(metabolomics = list(enrichment = list())))
    expect_message(
        res <- mod_metabolomics_network(de_res = NULL, pre = NULL,
                                        config = cfg, out_dir = tempdir()),
        "run_network"
    )
    expect_null(res)
})

test_that("network step is skipped when run_network = false", {
    cfg <- list(modes = list(metabolomics = list(
        enrichment = list(run_network = FALSE)
    )))
    expect_message(
        res <- mod_metabolomics_network(de_res = NULL, pre = NULL,
                                        config = cfg, out_dir = tempdir()),
        "run_network"
    )
    expect_null(res)
})

test_that("run_network = true proceeds past the toggle to the DE guard", {
    # Enabled -> the toggle must NOT fire; the function proceeds to its own
    # downstream guard (no DE tables here) and returns NULL without emitting
    # the toggle-skip message.
    cfg <- list(modes = list(metabolomics = list(
        enrichment = list(run_network = TRUE)
    )))
    msgs <- capture_messages(
        res <- mod_metabolomics_network(de_res = NULL, pre = NULL,
                                        config = cfg, out_dir = tempdir())
    )
    expect_false(any(grepl("run_network", msgs)))
    expect_null(res)
})
