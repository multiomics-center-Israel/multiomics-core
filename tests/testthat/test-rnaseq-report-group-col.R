# tests/testthat/test-rnaseq-report-group-col.R
#
# Regression test for #74: the RNA report's heatmap chunk chose its grouping
# column with its own chain (de_table$group_col, then effects$color[1] twice),
# which skipped filtering$group_col. When that differed from effects$color,
# the heatmap annotated samples by a different variable than the rest of the
# report. Every top-level `group_col <-` in the template must now agree with
# resolve_group_col().

rna_report_group_col_lines <- function() {
    candidates <- c(
        testthat::test_path("..", "..", "R", "domain", "rnaseq", "report_template.Rmd"),
        "R/domain/rnaseq/report_template.Rmd"
    )
    f <- candidates[file.exists(candidates)][1]
    if (is.na(f)) stop("R/domain/rnaseq/report_template.Rmd not found from the test directory")
    src <- readLines(f, warn = FALSE)
    grep("^group_col <- ", src, value = TRUE)
}

test_that("every group_col assignment in the RNA report follows resolve_group_col()", {
    lines <- rna_report_group_col_lines()
    # The design-table chunk and the heatmap chunk both assign it.
    expect_gte(length(lines), 2L)

    cfgs <- list(
        filtering_vs_color = list(filtering = list(group_col = "treatment"),
                                  effects   = list(color = c("batch", "treatment"))),
        filtering_vs_de_table = list(filtering = list(group_col = "treatment"),
                                     de_table  = list(group_col = "Condition"),
                                     effects   = list(color = "batch"))
    )
    for (nm in names(cfgs)) {
        for (ln in lines) {
            env <- new.env(parent = globalenv())
            env$rna_cfg <- cfgs[[nm]]
            eval(parse(text = ln), envir = env)
            expect_equal(env$group_col, "treatment", info = paste(nm, "::", ln))
        }
    }
})
