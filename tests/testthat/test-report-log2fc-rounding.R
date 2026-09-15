# tests/testthat/test-report-log2fc-rounding.R
#
# The proteomics report reconstructed log2 fold changes from the linearFC
# column, which is written with signif(x, 3). Anything in [1.4950, 1.5049] is
# stored as 1.50, and log2(1.50) = 0.5849625 is EXACTLY the log2(1.5) threshold,
# so features genuinely below 1.5-fold were counted as passing.
#
# On a real run this inflated the interactive volcano from 222 to 227:
# PLS3, SERINC3, BLTP3A, CKAP2 and ISYNA1, whose true ratios were 1.4957 to
# 1.4996. The tables, which use log2FC directly, reported 222 throughout.
#
# resolve_log2fc() lives in R/core/02_validation.R and is shared by the report,
# the PowerPoint deck, the executive summary, the static volcano tables and the
# pathway ranking, so they all count the same features.

# The five real features, with their true log2FC and their stored linearFC.
borderline <- data.frame(
    Gene                      = c("PLS3", "SERINC3", "BLTP3A", "CKAP2", "ISYNA1"),
    `log2FC.imputs.C_vs_V`    = c(-0.58329, -0.58461, 0.58085, 0.58300, 0.58202),
    `linearFC.imputs.C_vs_V`  = c(-1.5, -1.5, 1.5, 1.5, 1.5),
    check.names = FALSE, stringsAsFactors = FALSE
)

# The same features as the DE step writes them, all significant on padj.
borderline_summary <- function() {
    data.frame(
        FeatureID                = paste0("P", 1:5),
        `log2FC.imputs.C_vs_V`   = borderline$`log2FC.imputs.C_vs_V`,
        `linearFC.imputs.C_vs_V` = borderline$`linearFC.imputs.C_vs_V`,
        `pvalue.imputs.C_vs_V`   = rep(1e-4, 5),
        `padj.imputs.C_vs_V`     = rep(1e-3, 5),
        check.names = FALSE, stringsAsFactors = FALSE
    )
}

# Locate a repo file from either the repo root or tests/testthat.
repo_file <- function(...) {
    candidates <- c(testthat::test_path("..", "..", ...), file.path(...))
    candidates[file.exists(candidates)][1]
}

test_that("the stored log2FC is preferred over the rounded linearFC", {
    got <- resolve_log2fc(borderline, "C_vs_V")
    expect_equal(got, borderline$`log2FC.imputs.C_vs_V`)
})

test_that("all five borderline features fall BELOW the 1.5-fold threshold", {
    cut <- log2(1.5)
    expect_true(all(abs(resolve_log2fc(borderline, "C_vs_V")) < cut))
    # and their true ratios really are under 1.5
    expect_true(all(2^abs(borderline$`log2FC.imputs.C_vs_V`) < 1.5))
})

test_that("reconstructing from linearFC would wrongly pass all five", {
    # This is the bug, pinned so it cannot come back.
    cut <- log2(1.5)
    from_linear <- signed_fc_to_log2(borderline$`linearFC.imputs.C_vs_V`)
    expect_true(all(abs(from_linear) >= cut))
    expect_equal(sum(abs(from_linear) >= cut), 5L)
})

test_that("the fix changes the count by exactly the borderline features", {
    cut <- log2(1.5)
    correct <- sum(abs(resolve_log2fc(borderline, "C_vs_V")) >= cut)
    buggy <- sum(abs(signed_fc_to_log2(borderline$`linearFC.imputs.C_vs_V`)) >= cut)
    expect_equal(correct, 0L)
    expect_equal(buggy, 5L)
})

test_that("the linearFC fallback still works when log2FC is absent", {
    df <- borderline
    df$`log2FC.imputs.C_vs_V` <- NULL
    got <- resolve_log2fc(df, "C_vs_V")
    expect_equal(got, signed_fc_to_log2(df$`linearFC.imputs.C_vs_V`))
})

test_that("un-suffixed column names are also resolved", {
    df <- data.frame(`log2FC.C_vs_V` = c(1, -1), check.names = FALSE)
    expect_equal(resolve_log2fc(df, "C_vs_V"), c(1, -1))
})

test_that("a contrast with no fold-change column yields NA, not an error", {
    df <- data.frame(Gene = c("A", "B"), stringsAsFactors = FALSE)
    got <- resolve_log2fc(df, "C_vs_V")
    expect_length(got, 2)
    expect_true(all(is.na(got)))
})

test_that("features comfortably past the threshold are unaffected", {
    df <- data.frame(
        `log2FC.imputs.C_vs_V`   = c(2.0, -2.0),
        `linearFC.imputs.C_vs_V` = c(4.0, -4.0),
        check.names = FALSE
    )
    expect_true(all(abs(resolve_log2fc(df, "C_vs_V")) >= log2(1.5)))
})

test_that("the executive summary fallback count uses the stored log2FC", {
    cfg <- list(de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5))
    stats <- get_de_summary_stats_proteomics(borderline_summary(), cfg)
    expect_equal(stats$C_vs_V$n_sig, 0)
})

test_that("the static volcano tables carry the stored log2FC", {
    tabs <- qc_post_tables_from_summary(borderline_summary(), cfg = list(), use_adj = TRUE)
    expect_equal(tabs$C_vs_V$logFC, borderline$`log2FC.imputs.C_vs_V`)
})

test_that("the pathway DE table carries the stored log2FC", {
    config <- list(modes = list(proteomics = list(de_table = list(id_col = "FeatureID"))))
    tbl <- extract_de_table_for_pathway(borderline_summary(), "C_vs_V", config)
    expect_equal(tbl$log2FoldChange, borderline$`log2FC.imputs.C_vs_V`)
})

test_that("the pipeline summary fallback counts from the stored log2FC", {
    config <- list(
        project = list(name = "test"),
        modes   = list(proteomics = list(de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5)))
    )

    # No pass.imputs columns, so the summary recomputes from padj and log2FC.
    stats <- collect_proteomics_pipeline_stats(config, pre = list(),
                                               de_res = list(summary_df = borderline_summary()),
                                               pathway_res = NULL)
    expect_named(stats$de$contrasts, "C_vs_V")   # not "imputs.C_vs_V"
    expect_equal(stats$de$contrasts$C_vs_V$total, 0)

    # A table from before the log2FC column still falls back to linearFC.
    old_sdf <- borderline_summary()
    old_sdf$`log2FC.imputs.C_vs_V` <- NULL
    stats_old <- collect_proteomics_pipeline_stats(config, pre = list(),
                                                   de_res = list(summary_df = old_sdf),
                                                   pathway_res = NULL)
    expect_equal(stats_old$de$contrasts$C_vs_V$total, 5)
})

test_that("the pipeline summary prefers log2FoldChange in per-contrast tables", {
    config <- list(
        project = list(name = "test"),
        modes   = list(proteomics = list(de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5)))
    )
    tbl <- data.frame(
        padj           = rep(1e-3, 5),
        linearFC       = borderline$`linearFC.imputs.C_vs_V`,
        log2FoldChange = borderline$`log2FC.imputs.C_vs_V`
    )
    stats <- collect_proteomics_pipeline_stats(config, pre = list(),
                                               de_res = list(tables = list(C_vs_V = tbl)),
                                               pathway_res = NULL)
    expect_equal(stats$de$contrasts$C_vs_V$total, 0)
})

test_that("the report template uses the core helper instead of its own copy", {
    f <- repo_file("R", "domain", "proteomics", "report_template_proteomics.Rmd")
    skip_if(is.na(f), "report template not found from the test working directory")

    src <- readLines(f, warn = FALSE)
    expect_identical(grep("resolve_log2fc <- function", src, fixed = TRUE), integer(0))
    expect_true(any(grepl("resolve_log2fc(", src, fixed = TRUE)))
    expect_identical(grep("signed_fc_to_log2(lfc_vals)", src, fixed = TRUE), integer(0))
    expect_identical(grep("signed_fc_to_log2(as.numeric(df[[lfc_col]]))", src, fixed = TRUE),
                     integer(0))
})

test_that("pptx, summary, volcano and pathway code read log2FC through the helper", {
    rel_paths <- c(
        file.path("R", "domain", "proteomics", "11_powerpoint.R"),
        file.path("R", "domain", "proteomics", "09_executive_summary.R"),
        file.path("R", "domain", "proteomics", "07b_pathway.R"),
        file.path("R", "modules", "proteomics", "03_mod_qc_post.R")
    )
    for (rel in rel_paths) {
        f <- repo_file(rel)
        skip_if(is.na(f), paste(rel, "not found from the test working directory"))
        src <- readLines(f, warn = FALSE)
        expect_identical(grep("signed_fc_to_log2(", src, fixed = TRUE), integer(0),
                         label = paste("signed_fc_to_log2() lines in", rel))
        expect_true(any(grepl("resolve_log2fc(", src, fixed = TRUE)),
                    label = paste("resolve_log2fc() used in", rel))
    }
})

test_that("the pipeline summary and AI commentary no longer rebuild log2FC inline", {
    checks <- list(
        list(rel = file.path("R", "domain", "proteomics", "10_pipeline_summary.R"),
             inline = "log2(abs(fc_vals)) * sign(fc_vals)"),
        list(rel = file.path("R", "domain", "proteomics", "09b_commentary.R"),
             inline = "log2(pmax(abs(lfc_vals), 1e-10))")
    )
    for (chk in checks) {
        f <- repo_file(chk$rel)
        skip_if(is.na(f), paste(chk$rel, "not found from the test working directory"))
        src <- readLines(f, warn = FALSE)
        expect_identical(grep(chk$inline, src, fixed = TRUE), integer(0),
                         label = paste("inline log2 rebuild in", chk$rel))
        expect_true(any(grepl("resolve_log2fc(", src, fixed = TRUE)),
                    label = paste("resolve_log2fc() used in", chk$rel))
    }
})
