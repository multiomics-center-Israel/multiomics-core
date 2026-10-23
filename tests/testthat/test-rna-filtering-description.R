# RNA-seq reports must describe the gene filter that ran, not the config.
#
# describe_rna_filtering() is the one reading of preprocess_rna()'s info that
# the HTML report and the pipeline summary share. These tests pin its wording
# for every path (adaptive, fixed, deseq2_only, none, preprocessed, missing
# info), check that preprocess_rna() records the mode that actually ran, and
# check that the report template no longer carries the old, wrong wording.
# All fixtures are synthetic.

info_for <- function(mode, filter_norm = "CPM", threshold = NA, group = "Condition", ...) {
    c(list(filter_mode = mode, filter_norm = filter_norm,
           filter_threshold = threshold, filter_group_col = group), list(...))
}

test_that("adaptive: KDE valley, applied threshold, range, fallback and the half-of-group rule", {
    d <- describe_rna_filtering(
        info_for("adaptive", threshold = 0.734),
        list(auto_filter = list(min = 0.5, max = 2, fallback = 1)))

    expect_identical(d$mode, "adaptive")
    expect_true(d$applied)
    expect_identical(d$matrix, "CPM")
    expect_identical(d$group_col, "Condition")
    expect_equal(d$threshold, 0.734)
    expect_identical(d$threshold_label, "0.734 CPM")
    expect_match(d$methods, "kernel density estimate", fixed = TRUE)
    expect_match(d$methods, "first local minimum", fixed = TRUE)
    expect_match(d$methods, "0.5–2 CPM", fixed = TRUE)
    expect_match(d$methods, "the threshold is 2^v", fixed = TRUE)
    # The fallback is returned unclamped, so it must not read as constrained.
    expect_match(d$methods, "the fallback of 1 CPM is used as configured", fixed = TRUE)
    expect_no_match(d$methods, "constrained", fixed = TRUE)
    expect_match(d$methods, "applied in this run was 0.734 CPM", fixed = TRUE)
    expect_match(d$methods, "at least half (rounded up) of the samples with a non-missing value", fixed = TRUE)
    expect_match(d$methods, "'Condition'", fixed = TRUE)
    # No scan of candidate thresholds runs, and the rule is not "all samples".
    expect_no_match(d$methods, "candidate|optimal|ALL samples|complete", ignore.case = TRUE)
    expect_match(d$de_note, "independent filtering", fixed = TRUE)
})

test_that("adaptive on tximport input is described as TPM, not CPM", {
    d <- describe_rna_filtering(info_for("adaptive", filter_norm = "TPM (tximport)", threshold = 1.2))

    expect_identical(d$matrix, "TPM (tximport abundance)")
    expect_identical(d$threshold_label, "1.2 TPM")
    expect_match(d$methods, "log2(TPM + 0.1)", fixed = TRUE)
    expect_no_match(d$methods, "CPM", fixed = TRUE)
    expect_match(d$rule, "TPM ≥ the threshold", fixed = TRUE)
})

test_that("fixed: the fixed threshold and the same group rule, no density wording", {
    d <- describe_rna_filtering(info_for("fixed", threshold = 2))

    expect_true(d$applied)
    expect_identical(d$threshold_label, "2 CPM")
    expect_match(d$methods, "fixed threshold of 2 CPM", fixed = TRUE)
    expect_match(d$methods, "at least half (rounded up)", fixed = TRUE)
    expect_no_match(d$methods, "density", fixed = TRUE)
})

test_that("deseq2_only: zero-total genes removed, kept apart from DESeq2 independent filtering", {
    d <- describe_rna_filtering(info_for("deseq2_only", group = NA_character_))

    expect_true(d$applied)
    expect_true(is.na(d$threshold))
    expect_true(is.na(d$matrix))
    expect_true(is.na(d$group_col))
    expect_null(d$rule)
    expect_identical(d$threshold_label, "N/A (no threshold-based filtering)")
    expect_match(d$methods, "total count of zero across all samples", fixed = TRUE)
    expect_match(d$methods, "no expression threshold", fixed = TRUE)
    expect_no_match(d$methods, "DESeq2", fixed = TRUE)
    expect_match(d$de_note, "does not change which genes were kept", fixed = TRUE)
})

test_that("none: no filter; the VST zero-count drop is mentioned only when it removed genes", {
    plain <- describe_rna_filtering(info_for("none", filter_norm = "none", group = NA_character_,
                                             expr_work_method = "TMMlogCPM",
                                             n_genes_raw = 10, n_genes_filt = 10))
    expect_false(plain$applied)
    expect_identical(plain$threshold_label, "N/A (no threshold-based filtering)")
    expect_identical(plain$methods, "No gene filtering was applied; all genes were kept.")

    vst <- describe_rna_filtering(info_for("none", filter_norm = "none", group = NA_character_,
                                           expr_work_method = "VST",
                                           n_genes_raw = 10, n_genes_filt = 8))
    expect_match(vst$methods, "dropped before the variance-stabilizing transformation", fixed = TRUE)
})

test_that("preprocessed input: filtering skipped, no DESeq2 note", {
    d <- describe_rna_filtering(info_for("preprocessed", filter_norm = "preprocessed",
                                         group = NA_character_))

    expect_false(d$applied)
    expect_identical(d$short, "skipped (preprocessed input)")
    expect_match(d$methods, "preprocessed expression matrix", fixed = TRUE)
    expect_null(d$de_note)
})

test_that("the config mode is not read: info decides", {
    d <- describe_rna_filtering(info_for("adaptive", threshold = 1),
                                list(mode = "fixed", fixed_threshold = 5))
    expect_identical(d$mode, "adaptive")
    expect_identical(d$threshold_label, "1 CPM")
})

test_that("missing info says 'not recorded' instead of guessing", {
    d <- describe_rna_filtering(NULL, list(mode = "adaptive"))

    expect_identical(d$mode, "unknown")
    expect_false(d$applied)
    expect_identical(d$threshold_label, "N/A (not recorded for this report)")
    expect_match(d$methods, "not recorded", fixed = TRUE)

    # An info list from before filter_mode was recorded still resolves the
    # paths filter_norm identifies on its own.
    old <- describe_rna_filtering(list(filter_norm = "preprocessed", filter_threshold = NA))
    expect_identical(old$mode, "preprocessed")
})

# ---- preprocess_rna() records the mode that ran -----------------------------

# Five probe genes on top of a deterministic background that keeps TMM
# well-behaved. At a 10 CPM threshold every non-zero cell passes and every zero
# fails, so the expected kept set follows from the group rule alone.
probe_genes <- c("g_hi", "g_A", "g_half", "g_one", "g_zero")
n_background <- 12L

filter_fixture <- function() {
    samples <- paste0("S", 1:6)
    probes <- rbind(
        g_hi   = c(1000, 1000, 1000, 1000, 1000, 1000),
        g_A    = c( 500,  500,  500,    0,    0,    0),
        g_half = c( 500,  500,    0,    0,    0,    0),
        g_one  = c( 500,    0,    0,    0,    0,    0),
        g_zero = c(   0,    0,    0,    0,    0,    0)
    )
    background <- outer(seq_len(n_background), seq_len(6), function(i, j) 100 + 20 * i + 3 * j)
    rownames(background) <- paste0("bg", seq_len(n_background))
    mat <- rbind(probes, background)
    colnames(mat) <- samples
    counts <- data.frame(gene_id = rownames(mat), mat, check.names = FALSE,
                         stringsAsFactors = FALSE)
    metadata <- data.frame(SampleID = samples, Condition = rep(c("A", "B"), each = 3),
                           stringsAsFactors = FALSE)
    list(counts = counts, metadata = metadata)
}

filter_config <- function(filtering, out_dir = tempdir(), input = NULL) {
    list(
        paths = list(out = out_dir),
        modes = list(rna = list(
            id_columns = list(sample_col = "SampleID", gene_id = "gene_id"),
            effects    = list(samples = "SampleID"),
            input      = input,
            filtering  = c(list(group_col = "Condition"), filtering)
        ))
    )
}

run_pre <- function(filtering, ...) {
    suppressWarnings(suppressMessages(
        preprocess_rna(filter_fixture(), filter_config(filtering, ...), verbose = FALSE)
    ))
}

test_that("fixed mode records mode, matrix, threshold and group, and keeps the half-of-group gene", {
    skip_if_not_installed("edgeR")
    pre <- run_pre(list(mode = "fixed", fixed_threshold = 10))

    expect_identical(pre$info$filter_mode, "fixed")
    expect_identical(pre$info$filter_norm, "CPM")
    expect_equal(pre$info$filter_threshold, 10)
    expect_identical(pre$info$filter_group_col, "Condition")
    # g_half: 2 of 3 samples of group A pass, so it is kept, as the report now
    # says. g_one: 1 of 3 does not reach the rounded-up half.
    kept <- rownames(pre$expr_filt)
    expect_true(all(c("g_hi", "g_A", "g_half") %in% kept))
    expect_false(any(c("g_one", "g_zero") %in% kept))
})

test_that("deseq2_only and none record their mode with no threshold or group", {
    skip_if_not_installed("edgeR")
    d2 <- run_pre(list(mode = "deseq2_only"))
    expect_identical(d2$info$filter_mode, "deseq2_only")
    expect_true(is.na(d2$info$filter_threshold))
    expect_true(is.na(d2$info$filter_group_col))
    expect_identical(setdiff(c(probe_genes, paste0("bg", seq_len(n_background))),
                             rownames(d2$expr_filt)), "g_zero")

    nn <- run_pre(list(mode = "none"))
    expect_identical(nn$info$filter_mode, "none")
    expect_identical(nrow(nn$expr_filt), length(probe_genes) + n_background)
})

test_that("an unrecognised mode is recorded as the adaptive filter that ran", {
    skip_if_not_installed("edgeR")
    out <- withr::local_tempdir()
    pre <- run_pre(list(mode = "fixedd"), out_dir = out)

    expect_identical(pre$info$filter_mode, "adaptive")
    expect_false(is.na(pre$info$filter_threshold))
    expect_identical(describe_rna_filtering(pre$info)$mode, "adaptive")
})

test_that("preprocessed input records filtering as skipped whatever the mode says", {
    pre <- run_pre(list(mode = "fixed", fixed_threshold = 10),
                   input = list(format = "preprocessed"))

    expect_identical(pre$info$filter_mode, "preprocessed")
    expect_identical(nrow(pre$expr_filt), length(probe_genes) + n_background)
})

# ---- consumers ----------------------------------------------------------------

test_that("the pipeline summary's filtering step uses the recorded description", {
    info <- info_for("deseq2_only", group = NA_character_)
    pre <- list(info = info,
                expr_raw = matrix(1, 5, 2), expr_filt = matrix(1, 4, 2),
                meta = data.frame(Condition = c("A", "B")))
    config <- list(project = list(name = "P"),
                   modes = list(rna = list(filtering = list(group_col = "Condition",
                                                            mode = "adaptive"))))

    stats <- collect_pipeline_stats(config, pre, de_res = list())
    expect_identical(stats$methods$filtering, describe_rna_filtering(info, config$modes$rna$filtering))

    body <- generate_summary_body_r(stats)
    expect_match(body, "total count of zero across all samples", fixed = TRUE)
    expect_no_match(body, "automatic threshold", fixed = TRUE)
    expect_no_match(body, "Adaptive density-based", fixed = TRUE)
})

rna_report_template_lines <- function() {
    candidates <- c(
        testthat::test_path("..", "..", "R", "domain", "rnaseq", "report_template.Rmd"),
        "R/domain/rnaseq/report_template.Rmd"
    )
    f <- candidates[file.exists(candidates)][1]
    if (is.na(f)) stop("R/domain/rnaseq/report_template.Rmd not found from the test directory")
    readLines(f, warn = FALSE, encoding = "UTF-8")
}

test_that("the report template reads the recorded filter and drops the old wording", {
    src <- paste(rna_report_template_lines(), collapse = "\n")

    expect_match(src, "\nparams:\n  rna_info: NULL\n", fixed = TRUE)
    expect_match(src, "describe_rna_filtering(params$rna_info, rna_cfg$filtering)", fixed = TRUE)
    for (stale in c("ALL samples", "at least one complete", "candidate CPM threshold",
                    "For each candidate threshold", "adaptive CPM threshold approach",
                    "adaptive CPM-based gene filtering", "auto_filter$min")) {
        expect_false(grepl(stale, src, fixed = TRUE), info = stale)
    }
})

# ---- pre$info reaches the rendered report ------------------------------------

RNA_ROOT <- normalizePath(if (dir.exists("R")) "." else "../..", mustWork = FALSE)
RNA_TEMPLATE <- file.path(RNA_ROOT, "R", "domain", "rnaseq", "report_template.Rmd")

# Pull one named chunk (```{r <label>...} ... ```) out of an Rmd, verbatim.
extract_rna_chunk <- function(lines, label) {
    starts <- grep(sprintf("^```\\{r %s[,}]", label), lines)
    if (length(starts) != 1L) stop("expected exactly one '", label, "' chunk, found ", length(starts))
    ends <- grep("^```\\s*$", lines)
    lines[starts:ends[ends > starts][1]]
}

# Knit a harness made of the template's own params block, its rna_filt line and
# its filtering chunks, so the test cannot drift from what the report renders.
# A plot file is always present, as a stale one from an earlier run would be.
render_filtering_section <- function(dir, params = list()) {
    tmpl <- readLines(RNA_TEMPLATE, warn = FALSE, encoding = "UTF-8")
    front_end <- which(tmpl == "---")[2]
    p_start <- grep("^params:", tmpl[seq_len(front_end)])
    p_lines <- tmpl[p_start:(front_end - 1L)]

    plot_path <- file.path(dir, "filtering_threshold_qc.png")
    # Any valid PNG will do. Large enough for the default margins, which a tiny
    # device rejects with "figure margins too large".
    grDevices::png(plot_path, width = 400, height = 300)
    graphics::plot.new()
    grDevices::dev.off()

    rmd <- c(
        "---", 'title: "filtering render smoke test"',
        # Not self-contained, so the plot shows up as its own file name rather
        # than as a data URI that is hard to tell from the theme's assets.
        "output:", "  html_document:", "    self_contained: false",
        p_lines, "---", "",
        "```{r setup-harness, include=FALSE}",
        "rna_cfg <- list(filtering = list(auto_filter = list(min = 0.5, max = 2, fallback = 1)))",
        sprintf("filter_plot <- %s", deparse(plot_path)),
        "figure_legend <- function(text) cat(sprintf(\"<p class='figure-legend'>%s</p>\\n\", text))",
        grep("^rna_filt <- ", tmpl, value = TRUE),
        grep("^rna_de <- ", tmpl, value = TRUE),
        "```", "",
        extract_rna_chunk(tmpl, "filtering-strategy"), "",
        extract_rna_chunk(tmpl, "filter-plot"), "",
        extract_rna_chunk(tmpl, "legend-filter-plot")
    )
    path <- file.path(dir, "smoke.Rmd")
    writeLines(rmd, path, useBytes = TRUE)
    out <- rmarkdown::render(path, output_file = "smoke.html", output_dir = dir,
                             params = params, envir = new.env(parent = globalenv()),
                             quiet = TRUE)
    gsub("[[:space:]]+", " ", paste(readLines(out, warn = FALSE, encoding = "UTF-8"), collapse = " "))
}

skip_unless_rna_renderable <- function() {
    skip_if_not_installed("rmarkdown")
    skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
    skip_if_not(file.exists(RNA_TEMPLATE), "report template not found")
}

test_that("the report renders the recorded adaptive filter and its threshold plot", {
    skip_unless_rna_renderable()
    html <- render_filtering_section(
        withr::local_tempdir(),
        list(rna_info = info_for("adaptive", threshold = 0.734)))

    expect_match(html, "applied in this run was 0.734 CPM", fixed = TRUE)
    expect_match(html, "at least half (rounded up) of the samples with a non-missing value", fixed = TRUE)
    expect_match(html, "Threshold applied: 0.734 CPM.", fixed = TRUE)
    expect_match(html, "filtering_threshold_qc.png", fixed = TRUE)
})

test_that("the report hides a stale threshold plot for a non-adaptive run", {
    skip_unless_rna_renderable()
    html <- render_filtering_section(
        withr::local_tempdir(),
        list(rna_info = info_for("fixed", threshold = 2)))

    expect_match(html, "fixed threshold of 2 CPM", fixed = TRUE)
    expect_no_match(html, "figure-legend", fixed = TRUE)
    expect_no_match(html, "filtering_threshold_qc.png", fixed = TRUE)
})

test_that("without params$rna_info the report says 'not recorded' and shows no plot", {
    skip_unless_rna_renderable()
    html <- render_filtering_section(withr::local_tempdir())

    expect_match(html, "not recorded for this report", fixed = TRUE)
    expect_no_match(html, "figure-legend", fixed = TRUE)
    expect_no_match(html, "filtering_threshold_qc.png", fixed = TRUE)
})

test_that("render_rnaseq_report() forwards pre$info to every render as params$rna_info", {
    skip_if_not_installed("rmarkdown")
    # The renderer falls back to the template path under getwd().
    withr::local_dir(RNA_ROOT)
    run_dir <- file.path(withr::local_tempdir(), "rna")
    seen <- list()
    local_mocked_bindings(
        render = function(input, output_file, params = NULL, ...) {
            seen[[length(seen) + 1L]] <<- params
            file.create(output_file)
            output_file
        },
        .package = "rmarkdown"
    )
    info <- info_for("fixed", threshold = 2)

    suppressMessages(render_rnaseq_report(run_dir, config = list(project = list(name = "P")),
                                          pre_info = info))

    expect_gte(length(seen), 1L)
    for (p in seen) expect_identical(p$rna_info, info)
})

test_that("the rna_report target passes rna_pre$info to the renderer", {
    src <- paste(deparse(body(pipe_rnaseq)), collapse = " ")
    expect_match(src, "pre_info = rna_pre$info", fixed = TRUE)
})

test_that("the DESeq2 independent-filtering note shows only when DESeq2 ran", {
    skip_unless_rna_renderable()
    note <- "independent filtering does not change which genes were kept"
    deseq <- render_filtering_section(withr::local_tempdir(),
        list(rna_info = info_for("adaptive", threshold = 0.734),
             de_info = list(method = "deseq2", design = "~Group")))
    voom <- render_filtering_section(withr::local_tempdir(),
        list(rna_info = info_for("adaptive", threshold = 0.734),
             de_info = list(method = "limma_voom", design = "~Group")))
    expect_match(deseq, note, fixed = TRUE)
    expect_no_match(voom, note, fixed = TRUE)
})
