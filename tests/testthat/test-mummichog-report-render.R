# tests/testthat/test-mummichog-report-render.R
#
# Render smoke test for the metabolomics report's Mummichog section.
#
# The chunks are EXTRACTED FROM THE REAL TEMPLATE at test time (never copied
# here), so this cannot drift from what the pipeline actually renders. It checks
# that the section survives knitting in all three shapes it can take:
#   * full — ORA table + GSEA scatter/table + supporting-evidence drill-down
#   * ORA only — evidence/GSEA files genuinely unavailable
#   * nothing at all — no mummichog results
#
# Skipped when rmarkdown/pandoc are unavailable, like the report module itself.

# testthat runs with the working directory at tests/testthat; locate the repo
# root the same way helper.R does.
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
# returns, so the harness exercises the template's rendering contract.
render_fixture_sections <- function(with_gsea = TRUE, with_evidence = TRUE) {
  pw <- data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    pathway      = c("Fatty acid degradation", "Glycolysis"),
    overlap_size = c(2, 1),
    pathway_size = c(2, 4),
    "p-value"    = c(0.01, 0.30),
    "overlap_EmpiricalCompounds (id)" = c("E196,E3", "E2")
  )
  ora_plot <- plot_mummichog_bubble(pw, title = "ORA", p_cutoff = 0.05)

  gsea <- NULL
  gsea_plot <- NULL
  if (with_gsea) {
    gtab <- data.frame(
      check.names = FALSE, stringsAsFactors = FALSE,
      "Pathway"                         = c("Fatty acid degradation", "Glycolysis"),
      "Pathway size (model compounds)"  = c(2, 4),
      "Detected ECs"                    = c(2, 2),
      "Hits used in GSEA"               = c(2, 2),
      "ES"                              = c(0.8, -0.4),
      "NES"                             = c(1.9, -1.1),
      "P.Value"                         = c(0.002, 0.400),
      "padj"                            = c(0.004, 0.400),
      "Direction"                       = c("toward HL", "toward LL"),
      "Leading-edge EmpiricalCompounds" = c("E196;E3", "E2")
    )
    gsea <- list(table = gtab, n_ec = 24L, n_perm = 1000L, seed = 42L,
                 metric = "moderated_t",
                 metric_label = "moderated t statistic (DE column 'statistic')",
                 direction_note = mmc_gsea_direction_note(
                   "HL_vs_LL", "moderated t statistic"))
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
        "Match"                 = c(1, 1),
        "Conflict"              = c(1, 0),
        "Not assessed"          = c(0, 0)),
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
        "Agreement"                     = c("Match", "Conflict", "Match")),
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
    title     = "Mummichog pathway analysis — test (test model)",
    subtitle  = "HL vs LL, all features",
    plot      = gsea_plot %||% ora_plot,
    plot_kind = if (is.null(gsea_plot)) "ora_bubble" else "gsea_scatter",
    ora_plot  = ora_plot,
    table     = build_mummichog_pathway_table(pw),
    gsea      = gsea,
    gsea_plot = gsea_plot,
    evidence  = evidence,
    slug      = "HL_vs_LL"))
}

# Build and knit a harness Rmd made of the template's own setup + mummichog
# chunks. Returns the rendered HTML as a single string.
render_mummichog_section <- function(sections, dir) {
  tmpl  <- readLines(METAB_TEMPLATE, warn = FALSE)
  chunks <- c(
    extract_chunk(tmpl, "setup"),
    "",
    # The template's own unpack chunk, so the section flags the mummichog
    # section reads (mummi_sections, has_gsea, has_mummi_evidence) are derived
    # exactly as the real report derives them.
    extract_chunk(tmpl, "unpack-results"),
    "",
    "# Mummichog Pathway Analysis (m/z-based) {.tabset}",
    "",
    extract_chunk(tmpl, "mummichog-methodology"),
    "",
    extract_chunk(tmpl, "mummichog-evidence-helpers"),
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

  # Mirror the report module: render in a child of an environment that has the
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


test_that("the mummichog report section renders with GSEA and evidence", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("DT")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  skip_if_not(file.exists(METAB_TEMPLATE), "report template not found")

  html <- render_mummichog_section(render_fixture_sections(),
                                   withr::local_tempdir())

  # the four per-contrast tabs (matched on their section ids, so prose that
  # merely mentions a tab name cannot stand in for the tab itself)
  expect_match(html, 'id="pathway-plot"', fixed = TRUE)
  expect_match(html, 'id="ora-results-table"', fixed = TRUE)
  expect_match(html, 'id="gsea-results-table"', fixed = TRUE)
  expect_match(html, 'id="supporting-evidence"', fixed = TRUE)
  expect_match(html, 'id="pathway-drill-down"', fixed = TRUE)

  # the GSEA plot is the section's pathway plot, described honestly
  expect_match(html, "Mummichog GSEA over the full ranked EmpiricalCompound list")
  expect_match(html, "re-encodes significance")

  # the evidence terminology the biologist has to see, and the one it must not
  expect_match(html, "Pathway-matching candidate")
  expect_match(html, "not.{0,40}Best guess", perl = TRUE)
  expect_match(html, "Annotation confidence")
  expect_match(html, "Not assessed")

  # the drill-down really is collapsible, with the features nested inside
  expect_match(html, "<details>")
  expect_match(html, "Underlying measured features")
  expect_match(html, "feat_1")
  # NES direction is expressed in the contrast's own terms
  expect_match(html, "higher moderated t statistic in HL relative to LL")
})

test_that("the mummichog report section renders without GSEA or evidence", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("DT")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  skip_if_not(file.exists(METAB_TEMPLATE), "report template not found")

  html <- render_mummichog_section(
    render_fixture_sections(with_gsea = FALSE, with_evidence = FALSE),
    withr::local_tempdir())

  expect_match(html, 'id="pathway-plot"', fixed = TRUE)
  expect_match(html, 'id="ora-results-table"', fixed = TRUE)
  expect_match(html, "Mummichog ORA pathway analysis")
  # the optional tabs are simply absent, and the methodology says why
  expect_false(grepl('id="gsea-results-table"', html, fixed = TRUE))
  expect_false(grepl('id="supporting-evidence"', html, fixed = TRUE))
  expect_match(html, "was not available for this run")
})

test_that("the mummichog report section renders with no results at all", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("DT")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  skip_if_not(file.exists(METAB_TEMPLATE), "report template not found")

  html <- render_mummichog_section(list(), withr::local_tempdir())
  expect_true(nzchar(html))
  expect_false(grepl('id="supporting-evidence"', html, fixed = TRUE))
  expect_false(grepl('id="pathway-plot"', html, fixed = TRUE))
})
