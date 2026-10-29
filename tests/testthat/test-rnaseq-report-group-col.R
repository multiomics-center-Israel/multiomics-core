# tests/testthat/test-rnaseq-report-group-col.R
#
# Regression test for #74: the RNA report's heatmap chunk chose its grouping
# column with its own chain (de_table$group_col, then effects$color[1] twice),
# which skipped filtering$group_col. When that differed from effects$color,
# the heatmap annotated samples by a different variable than the rest of the
# report. Every top-level `group_col <-` in the template must now agree with
# resolve_group_col().
#
# The technical summary, both settings tables and the expression box-plot read
# filtering$group_col directly, so a config whose filtering block was commented
# out failed to render with "missing value where TRUE/FALSE needed". They now
# resolve it the same way, and no chunk reads the raw key any more.

rna_report_template_lines <- function() {
    candidates <- c(
        testthat::test_path("..", "..", "R", "domain", "rnaseq", "report_template.Rmd"),
        "R/domain/rnaseq/report_template.Rmd"
    )
    f <- candidates[file.exists(candidates)][1]
    if (is.na(f)) stop("R/domain/rnaseq/report_template.Rmd not found from the test directory")
    readLines(f, warn = FALSE)
}

rna_report_group_col_lines <- function() {
    grep("^(group_col|grp_col|group_col_exp) <- ", rna_report_template_lines(), value = TRUE)
}

test_that("every group_col assignment in the RNA report follows resolve_group_col()", {
    lines <- rna_report_group_col_lines()
    # Technical summary, design table, heatmap and expression box-plot.
    expect_gte(length(lines), 4L)

    cfgs <- list(
        filtering_vs_color = list(filtering = list(group_col = "treatment"),
                                  effects   = list(color = c("batch", "treatment"))),
        filtering_vs_de_table = list(filtering = list(group_col = "treatment"),
                                     de_table  = list(group_col = "Condition"),
                                     effects   = list(color = "batch")),
        # filtering block absent (commented out in the YAML)
        no_filtering = list(effects = list(color = c("treatment", "batch")))
    )
    for (nm in names(cfgs)) {
        for (ln in lines) {
            env <- new.env(parent = globalenv())
            env$rna_cfg <- cfgs[[nm]]
            eval(parse(text = ln), envir = env)
            var <- sub(" <- .*", "", ln)
            expect_equal(env[[var]], "treatment", info = paste(nm, "::", ln))
        }
    }
})

test_that("the RNA report never reads filtering$group_col directly", {
    src <- rna_report_template_lines()
    expect_length(grep("filtering$group_col", src, fixed = TRUE), 0L)
})
