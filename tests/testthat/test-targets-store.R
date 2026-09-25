# Regression tests for per-run targets stores. Before this, run.R used the
# store named in the shared _targets.yaml, so a run of one project built in,
# and `--fresh` destroyed, whichever store the file last pointed at.

cfg_for <- function(name, round = "A01", store = NULL) {
    list(project = list(name = name, analysis_round = round,
                        targets_store = store))
}

test_that("targets_store_for: store is named from project name and round", {
    expect_equal(targets_store_for(cfg_for("My Project-X", "A03")),
                 "_targets_my_project_x_a03")
    expect_equal(targets_store_for(list(project = list(name = "Solo"))),
                 "_targets_solo")
})

test_that("targets_store_for: two projects never share a store", {
    expect_false(identical(targets_store_for(cfg_for("ProjA")),
                           targets_store_for(cfg_for("ProjB"))))
    expect_false(identical(targets_store_for(cfg_for("ProjA", "A01")),
                           targets_store_for(cfg_for("ProjA", "A02"))))
})

test_that("targets_store_for: project$targets_store overrides the derived name", {
    expect_equal(targets_store_for(cfg_for("ProjA", store = "_targets_legacy")),
                 "_targets_legacy")
    expect_equal(targets_store_for(cfg_for("ProjA", store = "")),
                 "_targets_proja_a01")
})

test_that("targets_store_for: a config without a project name is an error", {
    expect_error(targets_store_for(list(project = list())), "project\\$name")
})

test_that("use_targets_store: writes a per-run yaml and points TAR_CONFIG at it", {
    root <- tempfile("test-targets-store-")
    dir.create(root)
    old <- Sys.getenv("TAR_CONFIG", unset = NA)
    on.exit({
        if (is.na(old)) Sys.unsetenv("TAR_CONFIG") else Sys.setenv(TAR_CONFIG = old)
        unlink(root, recursive = TRUE)
    }, add = TRUE)

    store <- use_targets_store(cfg_for("ProjA"), root = root)
    yaml_path <- file.path(root, "_targets_proja_a01.yaml")

    expect_equal(store, "_targets_proja_a01")
    expect_equal(Sys.getenv("TAR_CONFIG"), yaml_path)
    expect_equal(yaml::read_yaml(yaml_path)$main$store, "_targets_proja_a01")
    expect_equal(yaml::read_yaml(yaml_path)$main$script, "_targets.R")
})

test_that("assert_targets_store_owner: a store stamped by another project is refused", {
    store <- tempfile("_targets_shared_")
    on.exit(unlink(store, recursive = TRUE), add = TRUE)

    stamp_targets_store(store, cfg_for("ProjA"))

    expect_true(assert_targets_store_owner(store, cfg_for("ProjA")))
    expect_error(assert_targets_store_owner(store, cfg_for("ProjB")),
                 "belongs to 'ProjA A01'")
    expect_error(assert_targets_store_owner(store, cfg_for("ProjA", "A02")),
                 "belongs to")
})

test_that("assert_targets_store_owner: missing or unstamped stores are accepted", {
    store <- tempfile("_targets_new_")
    expect_true(assert_targets_store_owner(store, cfg_for("ProjA")))

    dir.create(store)
    on.exit(unlink(store, recursive = TRUE), add = TRUE)
    expect_true(assert_targets_store_owner(store, cfg_for("ProjA")))
})
