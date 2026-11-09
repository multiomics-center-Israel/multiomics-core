# tests/testthat/test-contrasts-used.R
#
# write_contrasts_used() (R/core/04_config.R) copies each mode's files$contrasts
# into execution_info, so a results folder carries which groups each contrast
# compares. Contrast names are free text; only Numerator / Denominator name the
# arms. All tables below are synthetic.

contrasts_fixture <- function() {
    root <- withr::local_tempdir(.local_envir = parent.frame())
    dir.create(file.path(root, "data", "study"), recursive = TRUE)
    # The contrast name deliberately differs from its group levels.
    writeLines(c("Contrast_name,Factor,Numerator,Denominator",
                 "T1_vs_base,group,drugA_4h,ctrl_0h"),
               file.path(root, "data", "study", "contrasts.csv"))
    writeLines(c("Contrast_name\tFactor\tNumerator\tDenominator",
                 "KO_vs_WT\tgenotype\tKO\tWT"),
               file.path(root, "data", "study", "contrasts_rna.tsv"))
    exec_dir <- file.path(root, "out", "execution_info")
    dir.create(exec_dir, recursive = TRUE)
    list(root = root, exec_dir = exec_dir)
}

test_that("each mode's contrasts file is copied verbatim, named by mode", {
    fx <- contrasts_fixture()
    config <- list(
        project = list(dir = fx$root), paths = list(raw = "data"),
        modes = list(
            proteomics = list(files = list(contrasts = "study/contrasts.csv")),
            rna        = list(files = list(contrasts = "study/contrasts_rna.tsv"))
        ))

    out <- write_contrasts_used(config, fx$exec_dir)

    expect_setequal(basename(out),
                    c("contrasts_used_proteomics.csv", "contrasts_used_rna.tsv"))
    expect_identical(
        readLines(file.path(fx$exec_dir, "contrasts_used_proteomics.csv")),
        readLines(file.path(fx$root, "data", "study", "contrasts.csv")))
    expect_identical(
        readLines(file.path(fx$exec_dir, "contrasts_used_rna.tsv")),
        readLines(file.path(fx$root, "data", "study", "contrasts_rna.tsv")))
})

test_that("an absolute contrasts path is honoured", {
    fx <- contrasts_fixture()
    abs_path <- file.path(fx$root, "data", "study", "contrasts.csv")
    config <- list(project = list(dir = "/nonexistent"), paths = list(raw = "data"),
                   modes = list(metabolomics = list(files = list(contrasts = abs_path))))
    out <- write_contrasts_used(config, fx$exec_dir)
    expect_identical(basename(out), "contrasts_used_metabolomics.csv")
})

test_that("modes without a contrasts file write nothing and do not fail", {
    fx <- contrasts_fixture()
    config <- list(project = list(dir = fx$root), paths = list(raw = "data"),
                   modes = list(rna = list(files = list(contrasts = "")),
                                proteomics = list(files = list())))
    expect_identical(write_contrasts_used(config, fx$exec_dir), character(0))
    expect_length(list.files(fx$exec_dir), 0L)
})

test_that("a missing contrasts file warns and is skipped", {
    fx <- contrasts_fixture()
    config <- list(project = list(dir = fx$root), paths = list(raw = "data"),
                   modes = list(proteomics = list(files = list(contrasts = "study/none.csv"))))
    expect_warning(out <- write_contrasts_used(config, fx$exec_dir), "not found")
    expect_identical(out, character(0))
})

test_that("write_execution_info includes the snapshot", {
    fx <- contrasts_fixture()
    run_dir <- file.path(fx$root, "out")
    config <- list(project = list(dir = fx$root), paths = list(raw = "data"),
                   modes = list(proteomics = list(files = list(contrasts = "study/contrasts.csv"))))
    files <- write_execution_info(config, run_dir)
    expect_true(file.path(run_dir, "execution_info", "contrasts_used_proteomics.csv") %in% files)
})
