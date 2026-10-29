# tests/testthat/test-mummichog-report-render.R
#
# Render smoke test for the metabolomics report's Mummichog section.
#
# The chunks are EXTRACTED FROM THE REAL TEMPLATE at test time (never copied
# here), so this cannot drift from what the pipeline actually renders. It checks
# that the section survives knitting in the shapes it can take:
#   * full      -- Summary, ORA plot, GSEA plot, Pathway Evidence, Advanced
#   * ORA only  -- GSEA and evidence genuinely unavailable
#   * nothing   -- no mummichog results at all
#
# The per-contrast contract pinned below: Summary -> ORA plot -> GSEA plot ->
# Pathway Evidence -> Advanced / Downloads, with the full tables under Advanced
# rather than as tabs of their own.
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
render_fixture_sections <- function(with_gsea = TRUE, with_evidence = TRUE,
                                    nominal = TRUE, gsea_only = FALSE) {
  pw <- data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    pathway      = c("Fatty acid degradation", "Glycolysis"),
    overlap_size = c(2, 1),
    pathway_size = c(2, 4),
    "p-value"    = if (nominal) c(0.01, 0.30) else c(0.20, 0.30),
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
      # gsea_only: Glycolysis becomes nominal in GSEA while its ORA p stays 0.30,
      # although an ORA evidence row exists for it.
      "P.Value"                         = if (!nominal) c(0.2, 0.4) else
                                            c(0.002, if (gsea_only) 0.010 else 0.400),
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

  tbl <- build_mummichog_pathway_table(pw)
  summary <- build_mummichog_summary(
    tbl, gsea_table = gsea$table, evidence = evidence,
    pathway_ids = c("Fatty acid degradation" = "xx00071",
                    "Glycolysis" = "xx00010"),
    library_ids = c("xx00071", "xx00010"), p_cutoff = 0.05)

  list("HL_vs_LL" = list(
    title     = "Mummichog pathway analysis - test (test model)",
    subtitle  = "HL vs LL, all features",
    plot      = plot_mummichog_bubble(pw, title = "ORA", p_cutoff = 0.05),
    table     = tbl,
    gsea      = gsea,
    gsea_plot = gsea_plot,
    evidence  = evidence,
    summary   = summary,
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
# Tab ORDER is part of this section's contract -- the Summary is what a reader
# meets first -- so the tests assert on positions, not just on presence.
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

# The data of the tables in one tab: the htmlwidget JSON payloads, without the
# tab's prose. The Summary legend defines every row label whether or not a row
# carries it, so a label assertion made on the whole tab cannot fail (or, negated,
# cannot pass); made on the payload it checks the rows.
tab_table_data <- function(tab) {
  m <- regmatches(tab, gregexpr(
    '<script[^>]*data-for="htmlwidget-[^"]*"[^>]*>.*?</script>', tab, perl = TRUE))[[1]]
  expect_gt(length(m), 0L, label = "number of table widgets in the tab")
  paste(m, collapse = " ")
}


test_that("the mummichog section renders with GSEA and evidence", {
  skip_unless_renderable()

  html <- render_mummichog_section(render_fixture_sections(),
                                   withr::local_tempdir())

  # The five per-contrast tabs, matched on their section ids so prose that
  # merely mentions a tab name cannot stand in for the tab itself, in order.
  tabs <- c("summary", "ora-plot", "gsea-plot", "pathway-evidence",
            "advanced-downloads")
  at <- vapply(tabs, function(id) anchor_at(html, id), integer(1))
  expect_false(anyNA(at), info = paste(tabs[is.na(at)], collapse = ", "))
  expect_true(all(diff(at) > 0))

  # The full tables are no longer tabs of their own.
  expect_false(grepl('id="gsea-results-table"', html, fixed = TRUE))
  expect_false(grepl('id="ora-results-table"', html, fixed = TRUE))
  expect_false(grepl('id="supporting-evidence"', html, fixed = TRUE))

  # The Summary: the top-ranked pathways with their compact cells and labels.
  summary_tab <- tab_body(html, "summary", "ora-plot")
  expect_match(summary_tab, "Top 20 Mummichog pathways for this contrast", fixed = TRUE)
  # pageLength 20 for this table only, so the whole top 20 fits one page.
  expect_match(summary_tab, '"pageLength":20', fixed = TRUE)
  summary_rows <- tab_table_data(summary_tab)
  expect_match(summary_tab, "Fatty acid degradation", fixed = TRUE)
  expect_match(summary_rows, "Mummichog-only candidate", fixed = TRUE)
  expect_match(summary_rows, "Not tested", fixed = TRUE)
  expect_match(summary_tab, "2 ECs / 3 features", fixed = TRUE)
  expect_match(summary_tab, "no adjusted pathway p-value is reported", fixed = TRUE)
  # Ranked, but not nominal: listed, and never called a candidate.
  expect_match(summary_tab, "Glycolysis", fixed = TRUE)
  expect_match(summary_rows, "No nominal Mummichog signal", fixed = TRUE)

  # Each plot tab carries its own plot, legend and one-line count.
  ora_tab <- tab_body(html, "ora-plot", "gsea-plot")
  expect_match(ora_tab, "Each bubble is a pathway")
  expect_match(ora_tab, "One pathway reached nominal p &lt; 0.05 in Mummichog ORA",
               fixed = TRUE)
  gsea_tab <- tab_body(html, "gsea-plot", "pathway-evidence")
  expect_match(gsea_tab, "Mummichog GSEA over the full ranked EmpiricalCompound list")
  expect_match(gsea_tab, "x = NES", fixed = TRUE)
  expect_match(gsea_tab, "-log10(raw GSEA p-value)", fixed = TRUE)
  expect_match(gsea_tab, "One pathway reached nominal p &lt; 0.05 in Mummichog GSEA; 1 remained significant after BH correction",
               fixed = TRUE)
  expect_false(grepl("Each bubble is a pathway", gsea_tab, fixed = TRUE))

  # Pathway Evidence: Summary pathways only, Pathway -> EC -> features.
  ev_tab <- tab_body(html, "pathway-evidence", "advanced-downloads")
  expect_match(ev_tab, "pathway-matching candidate", fixed = TRUE)
  expect_match(ev_tab, "<details>")
  expect_match(ev_tab, "E196 — Match", fixed = TRUE)
  expect_match(ev_tab, "Mummichog pathway candidate", fixed = TRUE)
  expect_match(ev_tab, "feat_1")
  expect_false(grepl("E2 — Match", ev_tab, fixed = TRUE))

  # Advanced keeps every full table, and names the exports that were written.
  expect_match(html, "`mummichog_gsea_table`|<code>mummichog_gsea_table</code>")
  expect_match(html, "<code>mummichog_evidence_features</code>", fixed = TRUE)
  expect_match(html, "Full Mummichog ORA table", fixed = TRUE)
  expect_match(html, "Full Mummichog GSEA table", fixed = TRUE)
  expect_match(html, "Pathway evidence summary (all pathways)", fixed = TRUE)
  expect_match(html, "Feature-level evidence", fixed = TRUE)

  # Methodology: both analyses named, NES described conservatively.
  expect_match(html, "complete ranked EC list, without a significance cutoff")
  expect_match(html, "unadjusted DE p-value below the configured Mummichog cutoff")
  expect_match(html, "retained to reproduce the pinned MetaboAnalystR result")
  expect_match(html, "should not be interpreted as biological up/down direction")
  expect_false(grepl("toward one end", html, fixed = TRUE))

  # Provenance and the evidence rules survive, in the collapsed block.
  expect_match(html, "Technical details / Reproducibility", fixed = TRUE)
  expect_match(html, "MetaboAnalystR")
  expect_match(html, substr(.MMC_GSEA_REF_COMMIT, 1, 12))
  expect_match(html, "not.{0,40}Best guess", perl = TRUE)
})

test_that("the mummichog section renders without GSEA or evidence", {
  skip_unless_renderable()

  html <- render_mummichog_section(
    render_fixture_sections(with_gsea = FALSE, with_evidence = FALSE),
    withr::local_tempdir())

  # ORA is unaffected by GSEA's absence.
  expect_match(html, 'id="ora-plot"', fixed = TRUE)
  expect_match(html, "Each bubble is a pathway")

  # The GSEA plot tab is rendered even with nothing to plot, and what it
  # carries is the explanation, NOT the ORA bubble wearing a GSEA heading.
  gsea_tab <- tab_body(html, "gsea-plot", "advanced-downloads")
  expect_match(gsea_tab, "No GSEA results for this contrast")
  expect_false(grepl("Each bubble is a pathway", gsea_tab, fixed = TRUE))
  expect_false(grepl("x = NES", gsea_tab, fixed = TRUE))
  expect_match(gsea_tab, "ORA plot")

  # The Summary falls back to ORA and says GSEA did not run.
  expect_match(tab_body(html, "summary", "ora-plot"), "Not run", fixed = TRUE)

  # Evidence is simply absent, and the methodology says why GSEA is.
  expect_false(grepl('id="pathway-evidence"', html, fixed = TRUE))
  # Advanced names only exports that exist for this contrast.
  expect_match(html, "<code>mummichog_summary</code>", fixed = TRUE)
  expect_false(grepl("mummichog_gsea_table", html, fixed = TRUE))
  expect_false(grepl("mummichog_evidence_", html, fixed = TRUE))
  expect_match(html, "has no results in this report")
  # It must NOT assert a cause: a NULL GSEA result means either that the
  # analysis could not run or that it ran and found nothing.
  expect_match(html, "no pathway passed the size filter")
})

test_that("a GSEA-only pathway shows no ORA evidence as its support", {
  skip_unless_renderable()

  html <- render_mummichog_section(render_fixture_sections(gsea_only = TRUE),
                                   withr::local_tempdir())
  # In the Summary, labelled as a GSEA-only signal ...
  summary_tab <- tab_body(html, "summary", "ora-plot")
  expect_match(summary_tab, "Glycolysis", fixed = TRUE)
  expect_match(summary_tab, "GSEA-only signal", fixed = TRUE)
  # ... and absent from the default Pathway Evidence view, though its evidence
  # row still exists under Advanced.
  ev_tab <- tab_body(html, "pathway-evidence", "advanced-downloads")
  expect_match(ev_tab, "Fatty acid degradation", fixed = TRUE)
  expect_false(grepl("Glycolysis", ev_tab, fixed = TRUE))
})

test_that("a contrast with no nominal pathway says so above its ranked rows", {
  skip_unless_renderable()

  html <- render_mummichog_section(render_fixture_sections(nominal = FALSE),
                                   withr::local_tempdir())
  # Pandoc escapes "<" in prose, so the message is matched in its HTML form.
  summary_tab <- tab_body(html, "summary", "ora-plot")
  expect_match(summary_tab,
               "No Mummichog pathway reached nominal p &lt; 0.05 in ORA or GSEA for this contrast.",
               fixed = TRUE)
  expect_match(summary_tab, "not significant results", fixed = TRUE)
  # The ranked rows are still listed, none as a candidate.
  summary_rows <- tab_table_data(summary_tab)
  expect_match(summary_tab, "Fatty acid degradation", fixed = TRUE)
  expect_match(summary_rows, "No nominal Mummichog signal", fixed = TRUE)
  expect_false(grepl("Mummichog-only candidate", summary_rows, fixed = TRUE))
  expect_match(tab_body(html, "pathway-evidence", "advanced-downloads"),
               "No Summary pathway has a nominal ORA result with ORA-derived evidence",
               fixed = TRUE)
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
  expect_false(grepl('id="summary"', html, fixed = TRUE))
  expect_false(grepl('id="gsea-plot"', html, fixed = TRUE))
  expect_false(grepl('id="ora-plot"', html, fixed = TRUE))
  expect_false(grepl('id="pathway-evidence"', html, fixed = TRUE))
})
