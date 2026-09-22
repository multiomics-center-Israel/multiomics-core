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
    expect_true(grepl('`r if (show_sample_concordance) "# Sample Concordance Across Omics"`',
                      src, fixed = TRUE))
    expect_true(grepl('`r if (show_mechanistic) "# Mechanistic & Causal Inference {.tabset}"`',
                      src, fixed = TRUE))
})

test_that("RNA-protein tabs are replaced by one note when there is no RNA layer", {
    src <- template_lines()
    headings <- static_headings(src)
    for (h in c("Concordance Distribution", "DE Scatter (All Contrasts)",
                "Per-Contrast DE Details", "Top Proteins by RNA Agreement")) {
        expect_false(any(grepl(h, headings, fixed = TRUE)), info = h)
        expect_true(any(grepl(sprintf('`r if (has_transcriptomics) "## %s', h), src,
                              fixed = TRUE)), info = h)
    }
    expect_true(any(grepl('`r if (!has_transcriptomics) "## RNA-Protein Analyses"`',
                          src, fixed = TRUE)))
    expect_true(any(grepl("rna-prot-absent-note, eval=!has_transcriptomics", src,
                          fixed = TRUE)))
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

# ---- pathway maps --------------------------------------------------------

test_that("pathway maps have their own top-level section after cross-omics enrichment", {
    src <- template_lines()
    headings <- static_headings(src)
    top <- grep("^# ", headings, value = TRUE)
    pv <- which(top == "# Pathway Maps (Pathview) {.tabset}")
    expect_length(pv, 1)
    expect_identical(top[pv - 1], "# Cross-Omics Enrichment {.tabset}")

    # Every pathview chunk sits inside that section, none back in MultiGSEA.
    at <- function(pattern) grep(pattern, src)[1]
    section_start <- at("^# Pathway Maps \\(Pathview\\)")
    section_end <- at("^`r if \\(show_harmonization\\)")
    for (chunk in c("pathview-setup", "per-omics-pathview-setup", "pathview-all-pdf",
                    "pathview-overview", "pathview-index-table", "pathview-plots",
                    "per-omics-pathview-metab", "per-omics-pathview-prot")) {
        line <- at(sprintf("^```\\{r %s,", chunk))
        expect_true(line > section_start && line < section_end, info = chunk)
    }
    expect_false(any(grepl("Multi-Omics Pathway Maps", headings, fixed = TRUE)))
})
