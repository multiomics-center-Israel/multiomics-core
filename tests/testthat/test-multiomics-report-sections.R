# tests/testthat/test-multiomics-report-sections.R
#
# The multiomics report drops sections that are switched off, empty, or a copy of
# another tab. A markdown heading written straight into the template renders
# whether or not the chunks under it run, so every heading these rules govern has
# to be emitted conditionally. These tests read the template as text: rendering
# it needs a full run directory, and what they pin is exactly that no governed
# heading is left static.

template_lines <- function() {
    f <- file.path(root_dir, "R", "domain", "multiomics",
                   "report_template_multiomics.Rmd")
    skip_if_not(file.exists(f), "multiomics report template not found")
    readLines(f, warn = FALSE)
}

# Headings outside code chunks. Chunk bodies are skipped because they emit their
# own headings with cat(), which is already conditional on the chunk running.
static_headings <- function(src) {
    in_chunk <- FALSE
    out <- character(0)
    for (line in src) {
        if (grepl("^```\\{", line)) { in_chunk <- TRUE; next }
        if (in_chunk && grepl("^```\\s*$", line)) { in_chunk <- FALSE; next }
        if (!in_chunk && grepl("^#{1,6} ", line)) out <- c(out, line)
    }
    out
}

test_that("switched-off sections leave no static heading behind", {
    headings <- static_headings(template_lines())
    governed <- c("Data Harmonization", "Sample Concordance Across Omics",
                  "Mechanistic & Causal Inference", "Protein-Metabolite Correlations",
                  "COSMOS Causal Network", "TF Activity", "Pathway Activity (PROGENy)",
                  "Regulatory / Co-expression Network")
    for (h in governed) {
        expect_false(any(grepl(h, headings, fixed = TRUE)),
                     info = paste("static heading still present:", h))
    }
})

test_that("each whole-section switch is read and gates its heading", {
    src <- paste(template_lines(), collapse = "\n")
    for (sw in c("harmonization", "sample_concordance", "mechanistic")) {
        expect_true(grepl(sprintf('show_section("%s")', sw), src, fixed = TRUE),
                    info = sw)
    }
    expect_true(grepl('`r if (show_harmonization) "# Data Harmonization"`', src, fixed = TRUE))
    # Sample Concordance shows only with content: its heading follows the group
    # distance result, which folds in the switch and is computed above it.
    expect_true(grepl('`r if (has_group_dist) "# Sample Concordance Across Omics"`',
                      src, fixed = TRUE))
    expect_true(grepl("has_group_dist <- show_sample_concordance &&", src, fixed = TRUE))
    setup_at <- as.integer(regexpr("```{r group-dist-setup", src, fixed = TRUE))
    heading_at <- as.integer(regexpr('"# Sample Concordance Across Omics"', src, fixed = TRUE))
    expect_gt(setup_at, 0)
    expect_lt(setup_at, heading_at)
    expect_true(grepl('`r if (show_mechanistic) "# Mechanistic & Causal Inference {.tabset}"`',
                      src, fixed = TRUE))
})

test_that("RNA-protein tabs appear only with results, otherwise one top-level note", {
    src <- template_lines()
    headings <- static_headings(src)
    for (h in c("Concordance Distribution", "DE Scatter (All Contrasts)",
                "Per-Contrast DE Details", "Top Proteins by RNA Agreement")) {
        expect_false(any(grepl(h, headings, fixed = TRUE)), info = h)
        # Per-contrast details need a contrast with files, the rest the results.
        gate <- if (h == "Per-Contrast DE Details") "has_rna_prot_contrasts" else "has_rna_prot"
        expect_true(any(grepl(sprintf('`r if (%s) "## %s', gate, h), src,
                              fixed = TRUE)), info = h)
    }
    expect_true(any(grepl('`r if (has_rna_prot && show_translation_eff) "## Translation Efficiency',
                          src, fixed = TRUE)))
    # Without results: the same top-level heading, no tabset, and one note that
    # says whether the RNA layer is missing or the analysis produced nothing.
    expect_true(any(grepl(paste0('`r if (has_rna_prot) "# RNA-Protein Correlation {.tabset}" ',
                                 'else "# RNA-Protein Correlation"`'), src, fixed = TRUE)))
    expect_true(any(grepl("rna-prot-note, eval=!has_rna_prot", src, fixed = TRUE)))
    expect_false(any(grepl("## RNA-Protein Analyses", src, fixed = TRUE)))
})

test_that("a single-contrast run gets one tab, not a combined and a per-contrast copy", {
    src <- template_lines()
    headings <- static_headings(src)
    expect_false(any(grepl("Per-Contrast", headings, fixed = TRUE)))
    expect_false(any(grepl("All Contrasts (", headings, fixed = TRUE)))
    for (chunk in c("multigsea-per-contrast", "multi-ora-per-contrast",
                    "enrichment-per-contrast")) {
        line <- grep(sprintf("^```\\{r %s,", chunk), src, value = TRUE)
        expect_length(line, 1)
        expect_true(grepl("!single_contrast", line, fixed = TRUE), info = chunk)
    }
})

test_that("the shipped config template leaves the new switches on", {
    cfg <- yaml::read_yaml(file.path(root_dir, "config", "templates",
                                     "multiomics_config.yaml"))
    sections <- (cfg$modes$multiomics$report %||% list())$sections %||% list()
    show_section <- function(name) !identical(sections[[name]], FALSE)
    for (sw in c("harmonization", "sample_concordance", "mechanistic")) {
        expect_true(show_section(sw), info = sw)
    }
})


# ---- the contrast count reads outputs as well as the config ------------------

contrast_dirs <- function(root, names) {
    for (n in names) dir.create(file.path(root, "per_contrast", n), recursive = TRUE)
    root
}

test_that("one configured contrast and one output is a single-contrast run", {
    enr <- contrast_dirs(withr::local_tempdir(), "A_vs_B")
    lay <- report_contrast_layout(list("A_vs_B"), enr, withr::local_tempdir())
    expect_true(lay$single_contrast)
    expect_identical(lay$n_contrasts, 1L)
    expect_identical(lay$label, "A vs B")
})

test_that("more contrasts in the outputs than in the config is not single-contrast", {
    # The design lists one contrast, but per-omics DE tables brought two: the
    # combined view pools both, so it must not be labelled as the one.
    enr <- contrast_dirs(withr::local_tempdir(), c("A_vs_B", "C_vs_D"))
    lay <- report_contrast_layout(list("A_vs_B"), enr, withr::local_tempdir())
    expect_false(lay$single_contrast)
    expect_identical(lay$n_contrasts, 2L)

    # The same from MultiGSEA's outputs alone.
    mg <- contrast_dirs(withr::local_tempdir(), c("A_vs_B", "C_vs_D", "E_vs_F"))
    lay <- report_contrast_layout(list("A_vs_B"), withr::local_tempdir(), mg)
    expect_false(lay$single_contrast)
    expect_identical(lay$n_contrasts, 3L)
})

test_that("Multi-ORA's own per-contrast outputs count too", {
    # MultiGSEA produced nothing per contrast, but Multi-ORA, which runs on its
    # own, wrote two: the Multi-ORA per-contrast section must not be hidden.
    mg <- withr::local_tempdir()
    contrast_dirs(file.path(mg, "multi_ora"), c("A_vs_B", "C_vs_D"))
    lay <- report_contrast_layout(list("A_vs_B"), withr::local_tempdir(), mg)
    expect_false(lay$single_contrast)
    expect_identical(lay$n_contrasts, 2L)

    # And a lone Multi-ORA directory names a single-contrast run.
    mg1 <- withr::local_tempdir()
    contrast_dirs(file.path(mg1, "multi_ora"), "P_vs_Q")
    expect_identical(report_contrast_layout(NULL, withr::local_tempdir(), mg1)$label,
                     "P vs Q")
})

test_that("more contrasts in the config than in the outputs is not single-contrast", {
    enr <- contrast_dirs(withr::local_tempdir(), "A_vs_B")
    lay <- report_contrast_layout(list("A_vs_B", "C_vs_D"), enr, withr::local_tempdir())
    expect_false(lay$single_contrast)
    expect_identical(lay$n_contrasts, 2L)
})

test_that("the lone contrast is named from the outputs first, then the config", {
    # Output directory wins over a differently spelled design entry.
    enr <- contrast_dirs(withr::local_tempdir(), "Treated_vs_Control")
    expect_identical(report_contrast_layout(list("treated - control"), enr,
                                            withr::local_tempdir())$label,
                     "Treated vs Control")
    # No enrichment output: MultiGSEA's directory names it.
    mg <- contrast_dirs(withr::local_tempdir(), "X_vs_Y")
    expect_identical(report_contrast_layout(NULL, withr::local_tempdir(), mg)$label,
                     "X vs Y")
    # No outputs at all: the design's own spelling.
    expect_identical(report_contrast_layout(list("A_vs_B"), withr::local_tempdir(),
                                            withr::local_tempdir())$label, "A_vs_B")
    # Nothing known anywhere.
    none <- report_contrast_layout(NULL, withr::local_tempdir(), withr::local_tempdir())
    expect_true(none$single_contrast)
    expect_identical(none$label, "Results")
})

test_that("the report takes its contrast layout from the resolver", {
    src <- paste(template_lines(), collapse = "\n")
    expect_true(grepl("report_contrast_layout(cfg$design$contrasts, enrichment_dir,",
                      src, fixed = TRUE))
    expect_true(grepl("single_contrast <- contrast_layout$single_contrast", src, fixed = TRUE))
    expect_true(grepl("single_contrast_label <- contrast_layout$label", src, fixed = TRUE))
})


# ---- RNA-protein per-contrast tabs follow real artifacts ---------------------

# The report's own setup chunk, evaluated as the report evaluates it.
rna_prot_setup <- function(rna_prot_dir, has_rna_prot = TRUE) {
    src <- paste(template_lines(), collapse = "\n")
    chunk <- regmatches(src, regexpr(
        "(?s)```\\{r rna-prot-per-contrast-setup[^}]*\\}\n.*?\n```", src, perl = TRUE))
    expect_length(chunk, 1)
    code <- sub("^```\\{r[^}]*\\}\n", "", sub("\n```$", "", chunk))
    env <- new.env()
    env$rna_prot_dir <- rna_prot_dir
    env$has_rna_prot <- has_rna_prot
    eval(parse(text = code), envir = env)
    env
}

test_that("the per-contrast RNA-protein tabset needs a contrast with files", {
    d <- withr::local_tempdir()
    # No per_contrast/ at all: the top-level results exist, but no DE join did.
    expect_false(rna_prot_setup(d)$has_rna_prot_contrasts)

    # One contrast with a table, one that joined nothing.
    dir.create(file.path(d, "per_contrast", "A_vs_B", "tables"), recursive = TRUE)
    file.create(file.path(d, "per_contrast", "A_vs_B", "tables",
                          "rna_protein_de_concordance.csv"))
    dir.create(file.path(d, "per_contrast", "C_vs_D"), recursive = TRUE)
    env <- rna_prot_setup(d)
    expect_true(env$has_rna_prot_contrasts)
    expect_identical(basename(env$rna_prot_contrast_dirs), "A_vs_B")

    # Without RNA-protein results the section is closed regardless.
    expect_false(rna_prot_setup(d, has_rna_prot = FALSE)$has_rna_prot_contrasts)

    src <- paste(template_lines(), collapse = "\n")
    expect_true(grepl('`r if (has_rna_prot_contrasts) "## Per-Contrast DE Details {.tabset}"`',
                      src, fixed = TRUE))
    expect_true(grepl("rna-prot-per-contrast, eval=has_rna_prot_contrasts", src, fixed = TRUE))
})


# ---- cross-omics heatmap legends --------------------------------------------

# The quoted text of one chunk, joined: what its figure_legend() calls can print.
chunk_text <- function(src, label) {
    start <- grep(sprintf("^```\\{r %s,", label), src)
    expect_length(start, 1)
    end <- start + which(grepl("^```\\s*$", src[(start + 1):length(src)]))[1]
    body <- src[(start + 1):(end - 1)]
    body <- body[!grepl("^\\s*#", body)]
    quoted <- regmatches(body, gregexpr('"[^"]*"', body))
    paste(gsub('"', "", unlist(quoted)), collapse = " ")
}

test_that("the cross-omics heatmap legends stay within three lines", {
    src <- template_lines()
    for (label in c("enrichment-heatmap-legend", "enrichment-ora-heatmap-legend")) {
        txt <- chunk_text(src, label)
        n_words <- length(strsplit(trimws(txt), "\\s+")[[1]])
        # 120 words is the ceiling asked for; 400 characters is about three lines
        # at the legend's width and font size.
        expect_lte(n_words, 120, label = label)
        expect_lte(nchar(txt), 400, label = label)
    }
})

test_that("GO figures say why metabolomics is absent from them", {
    src <- paste(template_lines(), collapse = "\n")
    expect_true(grepl("collection_note <- function(coll)", src, fixed = TRUE))
    expect_true(grepl('grepl("^GO", coll)', src, fixed = TRUE))
    # Both run-level loops and the per-contrast helper place the note.
    expect_gte(lengths(regmatches(src, gregexpr("collection_note(coll)", src,
                                                fixed = TRUE))), 3)
})
