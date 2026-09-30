# tests/testthat/test-run-name.R
#
# get_run_name(): the run directory name without its "Results_" prefix, used to
# name the Shiny payload files. Pure path logic, no files are touched.

make_run_config <- function(name = "Noam_Ziv", round = "A07") {
    list(project = list(dir = "proj", name = name, analysis_round = round))
}

test_that("get_run_name is the run directory name without the Results_ prefix", {
    cfg <- make_run_config()
    expect_identical(basename(get_run_out_dir(cfg)), "Results_Noam_Ziv_A07")
    expect_identical(get_run_name(cfg), "Noam_Ziv_A07")
})

test_that("get_run_name uses the same defaults as get_run_out_dir", {
    cfg <- list(project = list(dir = "proj"))
    expect_identical(get_run_name(cfg), "Project_Analysis")
})

test_that("get_run_name removes only the prefix that get_run_out_dir adds", {
    cfg <- make_run_config(name = "Results_X")
    expect_identical(basename(get_run_out_dir(cfg)), "Results_Results_X_A07")
    expect_identical(get_run_name(cfg), "Results_X_A07")
})
