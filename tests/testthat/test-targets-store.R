# Regression tests for per-run targets stores. Before this, run.R used the
# store named in the shared _targets.yaml, so a run of one project built in,
# and `--fresh` destroyed, whichever store the file last pointed at.

cfg_for <- function(name, round = "A01", store = NULL) {
    list(project = list(name = name, analysis_round = round,
                        targets_store = store))
}

# Keeps TAR_CONFIG as the test found it.
local_tar_config <- function(env = parent.frame()) {
    old <- Sys.getenv("TAR_CONFIG", unset = NA)
    withr::defer(restore_tar_config(old), envir = env)
    invisible(old)
}

test_that("targets_store_for: store is named from project name and round", {
    store <- targets_store_for(cfg_for("My Project-X", "A03"))
    expect_match(store, "^_targets_my_project_x_a03_[0-9a-f]{12}$")
    expect_match(targets_store_for(list(project = list(name = "Solo"))),
                 "^_targets_solo_[0-9a-f]{12}$")
    # Same config, same store: a rerun finds its cache.
    expect_identical(targets_store_for(cfg_for("My Project-X", "A03")), store)
})

test_that("targets_store_for: two projects never share a store", {
    expect_false(identical(targets_store_for(cfg_for("ProjA")),
                           targets_store_for(cfg_for("ProjB"))))
    expect_false(identical(targets_store_for(cfg_for("ProjA", "A01")),
                           targets_store_for(cfg_for("ProjA", "A02"))))
})

test_that("targets_store_for: names that slug alike still get their own store", {
    punct <- vapply(list(cfg_for("My Project-X"), cfg_for("My Project_X"),
                         cfg_for("My Project X"), cfg_for("My.Project.X")),
                    targets_store_for, character(1))
    expect_length(unique(punct), 4L)

    # Distinct on a case-insensitive filesystem too.
    cased <- vapply(list(cfg_for("ProjA"), cfg_for("proja"), cfg_for("PROJA"),
                         cfg_for("ProjA", "a01")),
                    targets_store_for, character(1))
    expect_length(unique(tolower(cased)), 4L)

    # Where the name ends and the round begins is part of the identity.
    split <- vapply(list(cfg_for("A B", "C"), cfg_for("A", "B C"),
                         cfg_for("A_B", NULL), cfg_for("A", "B")),
                    targets_store_for, character(1))
    expect_length(unique(tolower(split)), 4L)
})

test_that("targets_store_for: a config without a usable project name is an error", {
    expect_error(targets_store_for(list(project = list())), "project\\$name")
    expect_error(targets_store_for(cfg_for("   ")), "project\\$name")
    expect_error(targets_store_for(cfg_for(c("A", "B"))), "single value")
})

test_that("targets_store_for: project$targets_store overrides the derived name", {
    expect_equal(targets_store_for(cfg_for("ProjA", store = "_targets_legacy")),
                 "_targets_legacy")
    expect_equal(targets_store_for(cfg_for("ProjA", store = "_targets")), "_targets")
})

test_that("targets_store_for: an unsafe or malformed override is refused", {
    for (bad in list("", "   ", NA_character_, c("_targets_a", "_targets_b"), 1,
                     "/tmp/_targets_x", "C:/_targets_x", "~/_targets_x",
                     ".", "..", "../_targets_x", "_targets/../R", "_targets_x/sub",
                     "_targets..x", "sub\\_targets_x", "R", "data", "config",
                     "targets_x", "_targets x")) {
        expect_error(targets_store_for(cfg_for("ProjA", store = bad)),
                     "targets_store", info = paste(format(bad), collapse = " "))
    }
})

test_that("use_targets_store: writes a per-run yaml and points TAR_CONFIG at it", {
    root <- withr::local_tempdir()
    local_tar_config()

    store <- use_targets_store(cfg_for("ProjA"), root = root)
    yaml_path <- file.path(root, paste0(store, ".yaml"))

    expect_identical(store, targets_store_for(cfg_for("ProjA")))
    expect_equal(Sys.getenv("TAR_CONFIG"), yaml_path)
    expect_equal(yaml::read_yaml(yaml_path)$main$store, store)
    expect_equal(yaml::read_yaml(yaml_path)$main$script, "_targets.R")
})

test_that("use_targets_store: a refused store never becomes the active one", {
    root <- withr::local_tempdir()
    local_tar_config()
    Sys.setenv(TAR_CONFIG = "previous.yaml")

    store <- targets_store_for(cfg_for("ProjA"))
    stamp_targets_store(file.path(root, store), cfg_for("Someone else"))
    # A project whose derived store someone else stamped: use the override path.
    cfg <- cfg_for("ProjA", store = store)

    expect_error(use_targets_store(cfg, root = root), "Refusing to use it")
    expect_equal(Sys.getenv("TAR_CONFIG"), "previous.yaml")
    expect_false(file.exists(file.path(root, paste0(store, ".yaml"))))
})

test_that("restore_tar_config: a run leaves TAR_CONFIG as it found it, even on error", {
    root <- withr::local_tempdir()
    local_tar_config()
    # The shape run_pipeline() uses: save, register the restore, then switch.
    run <- function(fail) {
        old <- Sys.getenv("TAR_CONFIG", unset = NA)
        on.exit(restore_tar_config(old), add = TRUE)
        use_targets_store(cfg_for("ProjA"), root = root)
        if (fail) stop("pipeline failed")
        Sys.getenv("TAR_CONFIG")
    }

    Sys.unsetenv("TAR_CONFIG")
    expect_match(run(fail = FALSE), "_targets_proja_a01_")
    expect_true(is.na(Sys.getenv("TAR_CONFIG", unset = NA)))
    expect_error(run(fail = TRUE), "pipeline failed")
    expect_true(is.na(Sys.getenv("TAR_CONFIG", unset = NA)))

    Sys.setenv(TAR_CONFIG = "mine.yaml")
    run(fail = FALSE)
    expect_equal(Sys.getenv("TAR_CONFIG"), "mine.yaml")
    expect_error(run(fail = TRUE), "pipeline failed")
    expect_equal(Sys.getenv("TAR_CONFIG"), "mine.yaml")
})

test_that("run_pipeline() restores TAR_CONFIG and checks the store before using it", {
    src <- readLines(file.path(root_dir, "run.R"), warn = FALSE)
    start <- grep("^run_pipeline <- function", src)
    body <- src[start:(start + 60L)]
    save <- grep("old_tar_config <- Sys.getenv(\"TAR_CONFIG\", unset = NA)", body, fixed = TRUE)
    restore <- grep("on.exit(restore_tar_config(old_tar_config), add = TRUE)", body, fixed = TRUE)
    use <- grep("store <- use_targets_store(cfg_tmp)", body, fixed = TRUE)
    destroy <- grep("tar_destroy(ask = FALSE)", body, fixed = TRUE)
    expect_length(c(save, restore, use, destroy), 4L)
    expect_true(save < restore && restore < use && use < destroy)
    # No other Sys.setenv(TAR_CONFIG ...) in the runner.
    expect_false(any(grepl("Sys.setenv(TAR_CONFIG", src, fixed = TRUE)))
})

test_that("assert_targets_store_owner: a store stamped by another run is refused", {
    store <- withr::local_tempdir()
    stamp_targets_store(store, cfg_for("ProjA"))

    expect_true(assert_targets_store_owner(store, cfg_for("ProjA")))
    expect_error(assert_targets_store_owner(store, cfg_for("ProjB")),
                 "belongs to project 'ProjA' \\(round 'A01'\\)")
    expect_error(assert_targets_store_owner(store, cfg_for("ProjA", "A02")), "belongs to")
    expect_error(assert_targets_store_owner(store, cfg_for("ProjA", NULL)), "belongs to")
})

test_that("assert_targets_store_owner: identities that paste() alike are told apart", {
    store <- withr::local_tempdir()
    stamp_targets_store(store, cfg_for("A B", "C"))
    expect_true(assert_targets_store_owner(store, cfg_for("A B", "C")))
    expect_error(assert_targets_store_owner(store, cfg_for("A", "B C")), "belongs to")
    expect_error(assert_targets_store_owner(store, cfg_for("A B C", NULL)), "belongs to")
    expect_error(assert_targets_store_owner(store, cfg_for("a b", "C")), "belongs to")
})

test_that("assert_targets_store_owner: a malformed stamp fails closed", {
    store <- withr::local_tempdir()
    stamp <- targets_store_stamp_path(store)
    dir.create(dirname(stamp), recursive = TRUE)
    for (content in list("ProjA A01", "stamp_version: 2\nidentity: 'v1|5:ProjA|3:A01'",
                         "stamp_version: 1", "stamp_version: 1\nidentity: [a, b]",
                         "stamp_version: 1\nidentity: 'ProjA A01'", "{ not: [valid")) {
        writeLines(content, stamp)
        expect_error(assert_targets_store_owner(store, cfg_for("ProjA")),
                     "stamp .* is unreadable", info = content)
    }
})

test_that("assert_targets_store_owner: missing, empty or real unstamped stores are accepted", {
    missing <- file.path(withr::local_tempdir(), "_targets_new")
    expect_true(assert_targets_store_owner(missing, cfg_for("ProjA")))

    empty <- withr::local_tempdir()
    expect_true(assert_targets_store_owner(empty, cfg_for("ProjA")))

    legacy <- withr::local_tempdir()
    dir.create(file.path(legacy, "meta"))
    writeLines("name|type", file.path(legacy, "meta", "meta"))
    expect_true(assert_targets_store_owner(legacy, cfg_for("ProjA")))
})

test_that("assert_targets_store_owner: a folder that is not a store is refused", {
    folder <- withr::local_tempdir()
    writeLines("x", file.path(folder, "results.tsv"))
    expect_error(assert_targets_store_owner(folder, cfg_for("ProjA")),
                 "is not a targets store")
})
