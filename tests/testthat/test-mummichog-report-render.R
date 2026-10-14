# tests/testthat/test-mummichog-report-render.R
#
# Render smoke test for the metabolomics report's Mummichog section.
#
# The chunks are EXTRACTED FROM THE REAL TEMPLATE at test time (never copied
# here), so this cannot drift from what the pipeline actually renders. It checks
# that the section survives knitting in the three shapes it can take:
#   * full      -- ORA plot + table, GSEA scatter + table, evidence drill-down
#   * ORA only  -- GSEA and evidence genuinely unavailable
#   * nothing   -- no mummichog results at all
#
# Adapted from the integration branch of #266; the section contract differs
# there (that branch swapped the primary plot for the GSEA scatter), so the
# assertions below pin THIS contract: the bubble plot stays primary and GSEA is
# strictly additive.
#
# Skipped when rmarkdown/pandoc are unavailable, like the report module itself.

METAB_ROOT <- normalizePath(if (dir.exists("R")) "." else "../..",
                            mustWork = FALSE)
METAB_TEMPLATE <- file.path(METAB_ROOT, "R", "pipeline", "metabolomics",
                            "templates", "report_metabolomics.Rmd")

# Pull one named chunk (```{r <label>...} ... ```) out of an Rmd, verbatim.
extract_chunk <- function(lines, label) {
  starts <- grep(sprintf("^```\\{r %s[,}]", label), lines)
  if (length(starts) != 1L) {
    stop("expected exactly one '", label, "' chunk, found ", length(starts))
  }
  ends <- grep("^```\\s*$", lines)
  end  <- ends[ends > starts[1]][1]
  lines[starts[1]:end]
}

# A synthetic section list shaped exactly like build_mummichog_report_sections()
# returns, so the harness exercises the template's rendering contract without
# needing fgsea or a mummichog run.
render_fixture_sections <- function(with_gsea = TRUE, with_evidence = TRUE) {
  pw <- data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    pathway      = c("Fatty acid degradation", "Glycolysis"),
    overlap_size = c(2, 1),
    pathway_size = c(2, 4),
    "p-value"    = c(0.01, 0.30),
    "overlap_EmpiricalCompounds (id)" = c("E196,E3", "E2")
  )

  gsea <- NULL
  gsea_plot <- NULL
  if (with_gsea) {
    gtab <- data.frame(
      check.names = FALSE, stringsAsFactors = FALSE,
      "Pathway"                         = c("Fatty acid degradation", "Glycolysis"),
      "Pathway size (model compounds)"  = c(2, 4),
      "Detected ECs"                    = c(2, 2),
      "Tested size"                     = c(2, 2),
      "ES"                              = c(0.8, -0.4),
      "NES"                             = c(1.9, -1.1),
      "P.Value"                         = c(0.002, 0.400),
      "padj"                            = c(0.004, 0.400),
      "Leading-edge EmpiricalCompounds" = c("E196; E3", "E2"),
      "ES defaulted by reference implementation" = c(FALSE, FALSE),
      "ES fallback reason"                       = c(NA_character_, NA_character_)
    )
    gsea <- list(
      table = gtab, n_ec = 24L, n_pathways = 2L,
      n_es_defaulted = 0L, es_fallback_reasons = character(0),
      metric = "t_statistic",
      metric_label = "signed t statistic (DE method 'limma', column 'statistic')",
      de_method = "limma",
      nes_note = mmc_gsea_nes_note(),
      deviations = character(0),
      params = list(n_perm = 100L, min_size = 1, max_size = Inf,
                    gsea_param = 1, seed = 123L),
      reference = list(repo = "xia-lab/MetaboAnalystR",
                       commit = .MMC_GSEA_REF_COMMIT))
    gsea_plot <- plot_mummichog_gsea_scatter(gtab, title = "GSEA",
                                             p_cutoff = 0.05)
  }

  evidence <- NULL
  if (with_evidence) {
    evidence <- list(
      pathway_summary = data.frame(
        check.names = FALSE, stringsAsFactors = FALSE,
        "Pathway"               = c("Fatty acid degradation", "Glycolysis"),
        "Overlap"               = c(2, 1),
        "Detected pathway size" = c(2, 4),
        "Enrichment ratio"      = c(1, 0.25),
        "p.value"               = c(0.01, 0.30),
        "Supporting ECs"        = c(2, 1),
        "Supporting features"   = c(3, 1),
        "ECs Match"             = c(1, 1),
        "ECs Conflict"          = c(1, 0),
        "ECs Mixed"             = c(0, 0),
        "ECs Not assessed"      = c(0, 0),
        "feature-EC links Match"        = c(1, 1),
        "feature-EC links Conflict"     = c(1, 0),
        "feature-EC links Not assessed" = c(1, 0)),
      ec_table = data.frame(
        check.names = FALSE, stringsAsFactors = FALSE,
        "Pathway"                       = c("Fatty acid degradation",
                                            "Fatty acid degradation", "Glycolysis"),
        "EmpiricalCompound"             = c("E196", "E3", "E2"),
        "# Features"                    = c(2, 1, 1),
        "Pathway-matching candidate(s)" = c("C00020", "C00186", "C00031"),
        "Candidate name(s)"             = c("AMP", "(S)-Lactate", "D-Glucose"),
        "Candidate KEGG ID(s)"          = c("C00020", "C00186", "C00031"),
        "Original annotation"           = c("AMP", "D-Glucose", "D-glucose"),
        "Annotation confidence"         = c("Level 1", "Level 1", "Level 2"),
        "Agreement"                     = c("Match", "Conflict", "Match"),
        "n_match"                       = c(1, 0, 1),
        "n_conflict"                    = c(0, 1, 0),
        "n_not_assessed"                = c(1, 0, 0)),
      feature_table = data.frame(
        check.names = FALSE, stringsAsFactors = FALSE,
        "Pathway"               = c("Fatty acid degradation",
                                    "Fatty acid degradation",
                                    "Fatty acid degradation", "Glycolysis"),
        "EmpiricalCompound"     = c("E196", "E196", "E3", "E2"),
        "Feature"               = c("feat_1", "feat_2", "feat_4", "feat_3"),
        "m/z"                   = c(348.07, 370.05, 89.02, 181.07),
        "RT"                    = c(3.1, 3.1, 0.8, 1.2),
        "Adduct/ion"            = c("M+H[1+]", "M+Na[1+]", "M-H[-]", "M+H[1+]"),
        "Feature p-value"       = c(0.002, 0.010, 0.001, 0.030),
        "Feature statistic"     = c(4.1, 2.6, -5.5, -2.0),
        "Significant"           = c(TRUE, TRUE, TRUE, TRUE),
        "Original annotation"   = c("AMP", NA, "D-Glucose", "D-glucose"),
        "Annotation ID"         = c("C00020", NA, "C00031", NA),
        "Annotation confidence" = c("Level 1", NA, "Level 1", "Level 2"),
        "Agreement"             = c("Match", "Not assessed", "Conflict", "Match"))
    )
  }

  list("HL_vs_LL" = list(
    title     = "Mummichog pathway analysis - test (test model)",
    subtitle  = "HL vs LL, all features",
    plot      = plot_mummichog_bubble(pw, title = "ORA", p_cutoff = 0.05),
    table     = build_mummichog_pathway_table(pw),
    gsea      = gsea,
    gsea_plot = gsea_plot,
    evidence  = evidence,
    slug      = "HL_vs_LL"))
}

# Build and knit a harness Rmd made of the template's own setup + mummichog
# chunks. Returns the rendered HTML as a single string.
render_mummichog_section <- function(sections, dir) {
  tmpl <- readLines(METAB_TEMPLATE, warn = FALSE)
  chunks <- c(
    extract_chunk(tmpl, "setup"),
    "",
    # The template's own unpack chunk, so `mummi_sections` and `has_gsea` are
    # derived exactly as the real report derives them.
    extract_chunk(tmpl, "unpack-results"),
    "",
    "# Mummichog Pathway Analysis (m/z-based) {.tabset}",
    "",
    extract_chunk(tmpl, "mummichog-methodology"),
    "",
    extract_chunk(tmpl, "mummichog-per-contrast")
  )
  rmd <- c(
    "---", 'title: "mummichog render smoke test"',
    "output:", "  html_document:", "    self_contained: true",
    "params:", "  config: NULL", "  mummichog_sections: NULL",
    "  pre: NULL", "  qc_res: NULL", "  de_res: NULL", "  rf_res: NULL",
    "  plsda_res: NULL", "  enrichment_res: NULL", "  network_res: NULL",
    "  clustering_res: NULL", "  commentary_file: NULL", "---", "",
    chunks
  )
  path <- file.path(dir, "smoke.Rmd")
  writeLines(rmd, path)

  cfg <- list(
    project = list(name = "Test", analyst = "Tester", analysis_round = "A01"),
    modes = list(metabolomics = list(
      organism = "Coelastrella",
      enrichment = list(mummichog = list(
        enabled = TRUE, p_cutoff = 0.05, tolerance_ppm = 10,
        n_permutations = 100))))
  )

  # Mirror the report module: render in a child of an environment holding the
  # pipeline's R functions, which is where `%||%` and the plot builders live.
  fns <- new.env(parent = globalenv())
  for (f in list.files(file.path(METAB_ROOT, "R"), pattern = "\\.[Rr]$",
                       recursive = TRUE, full.names = TRUE)) {
    try(sys.source(f, envir = fns), silent = TRUE)
  }

  out <- rmarkdown::render(
    input = path, output_file = "smoke.html", output_dir = dir,
    params = list(config = cfg, mummichog_sections = sections),
    envir = new.env(parent = fns), quiet = TRUE
  )
  # Collapse whitespace: pandoc hard-wraps prose, so a phrase we assert on can
  # otherwise straddle a line break.
  gsub("[[:space:]]+", " ", paste(readLines(out, warn = FALSE), collapse = " "))
}

skip_unless_renderable <- function() {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("DT")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  skip_if_not(file.exists(METAB_TEMPLATE), "report template not found")
}

# Character offset of a pandoc section anchor, or NA when the tab is absent.
# Tab ORDER is part of this section's contract -- the GSEA scatter has to be the
# first plot a reader meets -- so the tests assert on positions, not just on
# presence.
anchor_at <- function(html, id) {
  m <- regexpr(sprintf('id="%s"', id), html, fixed = TRUE)
  if (m[[1]] == -1L) NA_integer_ else as.integer(m[[1]])
}

# The HTML between two anchors, i.e. one tab's own content. Lets a test assert
# what a tab does NOT contain, which is how "no ORA plot under a GSEA heading"
# is checked.
tab_body <- function(html, from, to) {
  a <- anchor_at(html, from)
  b <- anchor_at(html, to)
  expect_false(is.na(a), info = paste("missing tab:", from))
  expect_false(is.na(b), info = paste("missing tab:", to))
  substr(html, a, b)
}


test_that("the mummichog section renders with GSEA and evidence", {
  skip_unless_renderable()

  html <- render_mummichog_section(render_fixture_sections(),
                                   withr::local_tempdir())

  # The five per-contrast tabs, matched on their section ids so prose that
  # merely mentions a tab name cannot stand in for the tab itself.
  expect_match(html, 'id="gsea-plot"', fixed = TRUE)
  expect_match(html, 'id="gsea-results-table"', fixed = TRUE)
  expect_match(html, 'id="ora-bubble-plot"', fixed = TRUE)
  expect_match(html, 'id="ora-results-table"', fixed = TRUE)
  expect_match(html, 'id="supporting-evidence"', fixed = TRUE)
  expect_match(html, 'id="pathway-drill-down"', fixed = TRUE)

  # The contract: GSEA is the PRIMARY visualization. Its tab precedes every ORA
  # tab, so the first plot a reader meets in a contrast is the NES scatter.
  expect_lt(anchor_at(html, "gsea-plot"),
            anchor_at(html, "gsea-results-table"))
  expect_lt(anchor_at(html, "gsea-results-table"),
            anchor_at(html, "ora-bubble-plot"))
  expect_lt(anchor_at(html, "ora-bubble-plot"),
            anchor_at(html, "ora-results-table"))
  expect_lt(anchor_at(html, "ora-results-table"),
            anchor_at(html, "supporting-evidence"))

  # ... and the plot in that first tab is the scatter, described by its own
  # NES / -log10(GSEA p) legend, not the enrichment-ratio bubble.
  gsea_tab <- tab_body(html, "gsea-plot", "gsea-results-table")
  expect_match(gsea_tab, "Mummichog GSEA over the full ranked EmpiricalCompound list")
  # The two axes named separately, so a reflow of the legend prose between them
  # cannot break the assertion that this is the NES-vs-p scatter.
  expect_match(gsea_tab, "x = NES", fixed = TRUE)
  expect_match(gsea_tab, "-log10(raw GSEA p-value)", fixed = TRUE)
  expect_false(grepl("Each bubble is a pathway", gsea_tab, fixed = TRUE))

  # ORA is preserved, demoted in position only: its plot and legend still render
  # under their own heading.
  expect_match(tab_body(html, "ora-bubble-plot", "ora-results-table"),
               "Each bubble is a pathway")

  # Both analyses are named and distinguished in the methodology block.
  expect_match(html, "complete ranked EmpiricalCompound list with no significance cutoff")

  # NES is described conservatively: parity only, no direction of any kind.
  expect_match(html, "retained to reproduce the pinned MetaboAnalystR result")
  expect_match(html, "should not be interpreted as biological up/down direction")

  # The pinned reference is stated in the report.
  expect_match(html, "MetaboAnalystR")
  expect_match(html, substr(.MMC_GSEA_REF_COMMIT, 1, 12))

  # The evidence terminology the biologist has to see, and the one it must not.
  expect_match(html, "Pathway-matching candidate")
  expect_match(html, "not.{0,40}Best guess", perl = TRUE)
  expect_match(html, "<details>")
  expect_match(html, "feat_1")
})

test_that("the mummichog section renders without GSEA or evidence", {
  skip_unless_renderable()

  html <- render_mummichog_section(
    render_fixture_sections(with_gsea = FALSE, with_evidence = FALSE),
    withr::local_tempdir())

  # The ORA tabs are unaffected by GSEA's absence: ORA results stay available
  # under their own headings.
  expect_match(html, 'id="ora-bubble-plot"', fixed = TRUE)
  expect_match(html, 'id="ora-results-table"', fixed = TRUE)
  expect_match(html, "Each bubble is a pathway")

  # The GSEA Plot tab is rendered even with nothing to plot, so the primary tab
  # always answers what GSEA found. Its results table is the tab that goes away.
  expect_match(html, 'id="gsea-plot"', fixed = TRUE)
  expect_false(grepl('id="gsea-results-table"', html, fixed = TRUE))
  expect_lt(anchor_at(html, "gsea-plot"), anchor_at(html, "ora-bubble-plot"))

  # And what it carries is the explanation, NOT the ORA bubble wearing a GSEA
  # heading. This is the substitution the section must never make: a reader on
  # the GSEA tab must not be shown an enrichment-ratio plot.
  gsea_tab <- tab_body(html, "gsea-plot", "ora-bubble-plot")
  expect_match(gsea_tab, "No GSEA results for this contrast")
  expect_false(grepl("Each bubble is a pathway", gsea_tab, fixed = TRUE))
  expect_false(grepl("x = NES", gsea_tab, fixed = TRUE))
  # It points the reader at where the ORA results actually are.
  expect_match(gsea_tab, "ORA Bubble Plot")

  # Evidence is simply absent, and the methodology says why rather than leaving
  # a reader to infer it from a missing tab.
  expect_false(grepl('id="supporting-evidence"', html, fixed = TRUE))
  expect_match(html, "has no results in this report")
  # And it must NOT assert a cause: a NULL GSEA result means either that the
  # analysis could not run or that it ran and found nothing, and the section
  # cannot tell them apart. Naming only the setup causes would report a
  # legitimate negative result as a misconfiguration.
  expect_match(html, "no pathway passed the size filter")
})

test_that("a contrast without GSEA says so when another contrast has it", {
  skip_unless_renderable()

  # Mixed outcome: the methodology block takes the has_gsea branch and promises
  # two analyses per contrast, so a contrast that returned NULL must say why
  # rather than just missing a tab.
  secs <- render_fixture_sections()
  absent <- render_fixture_sections(with_gsea = FALSE, with_evidence = FALSE)[[1]]
  absent$slug <- "XX_vs_YY"
  secs[["XX_vs_YY"]] <- absent

  html <- render_mummichog_section(secs, withr::local_tempdir())

  # The contrast that has it still renders its real GSEA section ...
  expect_match(html, "Mummichog GSEA over the full ranked EmpiricalCompound list")
  # ... and the one that does not is explained, without asserting a cause it
  # cannot know.
  expect_match(html, "No GSEA results for this contrast")
  expect_match(html, "no pathway passed the size filter")
})

test_that("the mummichog section renders with no results at all", {
  skip_unless_renderable()

  html <- render_mummichog_section(list(), withr::local_tempdir())
  expect_true(nzchar(html))
  expect_false(grepl('id="gsea-plot"', html, fixed = TRUE))
  expect_false(grepl('id="ora-bubble-plot"', html, fixed = TRUE))
  expect_false(grepl('id="supporting-evidence"', html, fixed = TRUE))
})
