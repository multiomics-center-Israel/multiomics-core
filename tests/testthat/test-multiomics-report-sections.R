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

test_that("the cross-omics heatmap legends stay within 100 words", {
    # One rule for every legend: at most 100 words, counted with the optional
    # multi-contrast sentence included, since that is the longest it renders.
    src <- template_lines()
    for (label in c("enrichment-heatmap-legend", "enrichment-ora-heatmap-legend")) {
        txt <- chunk_text(src, label)
        n_words <- length(strsplit(trimws(txt), "\\s+")[[1]])
        expect_lte(n_words, 100, label = label)
    }
})

test_that("the shortened legends keep the limits a reader must not read past", {
    src <- template_lines()
    nes <- chunk_text(src, "enrichment-heatmap-legend")
    for (claim in c("not p near 1",
                    "A missing layer was not necessarily tested and found negative",
                    "is not a significance threshold",
                    "scored only with its other test is left out for that layer",
                    "the run log names it")) {
        expect_true(grepl(claim, gsub("\\s+", " ", nes), fixed = TRUE), info = claim)
    }
    ora <- gsub("\\s+", " ", chunk_text(src, "enrichment-ora-heatmap-legend"))
    for (claim in c("nothing is re-adjusted across contrasts",
                    "GSEA results do not enter",
                    "not p near 1",
                    "a display choice, not an error rate controlled across contrasts")) {
        expect_true(grepl(claim, ora, fixed = TRUE), info = claim)
    }
    # External ORA tables need not be BH, so the legend claims only "adjusted".
    expect_false(grepl("BH", ora, fixed = TRUE))
})

test_that("GO figures say why metabolomics is absent from them", {
    src <- paste(template_lines(), collapse = "\n")
    expect_true(grepl("collection_note <- function(coll)", src, fixed = TRUE))
    expect_true(grepl('grepl("^GO", coll, ignore.case = TRUE)', src, fixed = TRUE))
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

test_that("maps without a tab of their own leave no near-empty tab or second note", {
    src <- template_lines()
    headings <- static_headings(src)
    # Overview is the one tab that always renders; Multi-Omics Maps only with maps.
    expect_false(any(grepl("Multi-Omics Maps", headings, fixed = TRUE)))
    expect_true(any(grepl('`r if (has_pathview) "## Multi-Omics Maps {.unnumbered}"`',
                          src, fixed = TRUE)))
    expect_false(any(grepl("pathview-note|pathview-none", src)))
    # Names and contrasts come from the shared index, not parsed again in a tab.
    joined <- paste(src, collapse = "\n")
    expect_true(grepl("pv_index <- pathview_map_index(", joined, fixed = TRUE))
    expect_true(grepl("pathview_missing_note(multi_cfg)", joined, fixed = TRUE))
    expect_false(grepl("pathview_map_title <- function", joined, fixed = TRUE))
    expect_false(grepl("contrast_part", joined, fixed = TRUE))
})

# Writes empty PNG files, plus a KGML carrying a title for the ids in `titled`.
pathview_fixture <- function(files, titled = list()) {
    d <- withr::local_tempdir(.local_envir = parent.frame())
    for (f in files) file.create(file.path(d, f))
    for (stem in names(titled)) {
        writeLines(c('<?xml version="1.0"?>',
                     sprintf('<pathway name="path:%s" title="%s">', stem, titled[[stem]])),
                   file.path(d, paste0(stem, ".xml")))
    }
    d
}

test_that("a map's contrast is read from its filename, with the recorded spelling first", {
    expect_identical(pathview_map_contrast("ko00010.multi_ora_A.vs.B.multi.png"), "A vs B")
    expect_identical(pathview_map_contrast("ko00010.multi_ora_A.vs.B.png"), "A vs B")
    # The renderer's record wins over the dotted key.
    expect_identical(pathview_map_contrast("ko00010.multi_ora_X1.5.vs.X0.multi.png",
                                           list(X1.5.vs.X0 = "1.5 - 0")), "1.5 - 0")
    # A single-contrast map and the per-layer maps name no contrast.
    expect_identical(pathview_map_contrast(c("ko00010.multi_ora.multi.png",
                                             "hsa00010.metab_top.png",
                                             "hsa00010.prot_top.png")),
                     c("", "", ""))
    expect_identical(pathview_map_contrast(character(0)), character(0))
})

test_that("the map index gives each contrast of a pathway its own row", {
    d <- pathview_fixture(
        c("ko00010.multi_ora_A.vs.B.multi.png", "ko00010.multi_ora_C.vs.D.multi.png",
          "ko00010.png", "hsa00020.metab_top.png", "hsa00030.prot_top.png"),
        titled = list(ko00010 = "Glycolysis"))
    idx <- pathview_map_index(d)
    expect_identical(idx$set, c("multi_ora", "multi_ora", "metab_top", "prot_top"))
    expect_identical(idx$contrast, c("A vs B", "C vs D", "", ""))
    expect_identical(idx$kegg_id, c("00010", "00010", "00020", "00030"))
    expect_identical(idx$pathway, c("Glycolysis", "Glycolysis", "", ""))
    expect_identical(idx$heading, c("Glycolysis (00010) - A vs B",
                                    "Glycolysis (00010) - C vs D",
                                    "Pathway 00020", "Pathway 00030"))
    # Without the renderer's record, the set claims no number of layers.
    expect_identical(idx$map_set[1], "Cross-omics pathways")
    # The blank KEGG template beside the overlays is not a map.
    expect_false(any(basename(idx$png) == "ko00010.png"))
    # No two rows describe the same thing.
    expect_false(anyDuplicated(idx[, c("set", "contrast", "kegg_id")]) > 0)
})

test_that("the map index includes only the sets the report shows", {
    d <- pathview_fixture(c("ko00010.multi_ora_A.vs.B.multi.png",
                            "hsa00020.metab_top.png", "hsa00030.prot_top.png"))
    idx <- pathview_map_index(d, multi_ora = FALSE, prot_top = FALSE)
    expect_identical(idx$set, "metab_top")
    union <- pathview_map_index(d, metab_top = FALSE, prot_top = FALSE,
                                multi_ora_label = "Enriched in a gene layer")
    expect_identical(union$map_set, "Enriched in a gene layer")

    none <- pathview_map_index(file.path(d, "absent"))
    expect_identical(nrow(none), 0L)
    expect_identical(names(none), c("set", "map_set", "contrast", "kegg_id",
                                    "pathway", "heading", "png"))
    expect_identical(nrow(pathview_map_index(d, FALSE, FALSE, FALSE)), 0L)
})

test_that("with no maps the section says whether they were switched off", {
    off <- pathview_missing_note(list(enrichment = list(pathview = list(run_pathview = FALSE))))
    expect_match(off, "switched off", fixed = TRUE)
    expect_match(off, "enrichment.pathview.run_pathview", fixed = TRUE)
    # Same default as run_multi_ora(): an absent switch means maps were asked for.
    for (cfg in list(NULL, list(), list(enrichment = list(pathview = list(run_pathview = TRUE))))) {
        on <- pathview_missing_note(cfg)
        expect_match(on, "No pathway maps are available", fixed = TRUE)
        expect_false(grepl("switched off", on, fixed = TRUE))
    }
})

test_that("per-layer maps pathview wrote as .multi.png are indexed too", {
    # The per-layer renderer keeps whichever of the two names pathview wrote;
    # the combined PDF already took both, so the table and tab must as well.
    d <- pathview_fixture(c("hsa00020.metab_top.multi.png", "hsa00030.prot_top.multi.png",
                            "hsa00040.metab_top.png"),
                          titled = list(hsa00020 = "Citrate cycle"))
    idx <- pathview_map_index(d, multi_ora = FALSE)
    expect_identical(idx$set, c("metab_top", "metab_top", "prot_top"))
    expect_identical(idx$kegg_id, c("00020", "00040", "00030"))
    expect_identical(idx$heading[1], "Citrate cycle (00020)")
    expect_identical(idx$contrast, c("", "", ""))
})

test_that("the map set is named from the tier the renderer recorded", {
    two <- pathview_multi_ora_support(FALSE, 2L)
    expect_identical(two$label, "Enriched in >= 2 layers")
    expect_match(two$intro, ">= 2 omics layers", fixed = TRUE)

    # The single-layer fallback must not read as a two-layer result.
    one <- pathview_multi_ora_support(FALSE, 1L)
    expect_identical(one$label, "Enriched in one layer or more")
    expect_false(grepl(">= 2", paste(one$label, one$intro), fixed = TRUE))
    expect_match(one$intro, "one layer alone", fixed = TRUE)

    # No record (a run from before it existed): no claim about layers at all.
    for (none in list(NULL, NA, "unknown")) {
        res <- pathview_multi_ora_support(FALSE, none)
        expect_identical(res$label, "Cross-omics pathways")
        expect_false(grepl("layer", paste(res$label, res$intro), fixed = TRUE))
    }

    # The union renderer is its own case whatever the supported record says.
    expect_identical(pathview_multi_ora_support(TRUE, 2L)$label, "Enriched in a gene layer")

    # The report takes both from the resolver, never from fixed text.
    src <- paste(template_lines(), collapse = "\n")
    expect_true(grepl("multi_ora_label = pathview_support$label", src, fixed = TRUE))
    expect_true(grepl("cat(pathview_support$intro,", src, fixed = TRUE))
    expect_false(grepl("**>= 2 omics layers**", src, fixed = TRUE))
})

test_that("contrast spellings are read from the renderer that drew the maps", {
    src <- paste(template_lines(), collapse = "\n")
    expect_true(grepl(paste0("pathview_active_meta <- if (pathview_is_union) ",
                             "pathview_union_meta else pathview_supported_meta"),
                      src, fixed = TRUE))
    expect_true(grepl("pathview_active_meta$contrast_labels %||% list()", src, fixed = TRUE))
    # With the recorded spelling, a decimal contrast survives make.names().
    key <- make.names("1.5 - 0")
    png <- sprintf("hsa00010.multi_ora_%s.multi.png", key)
    labels <- stats::setNames(list("1.5 - 0"), key)
    expect_identical(pathview_map_contrast(png, labels), "1.5 - 0")
    # Without it (an older run), the existing fallback still applies.
    expect_identical(pathview_map_contrast(png), gsub("\\.", " ", key))
})

test_that("the supported renderer records its tier and contrast spellings", {
    # Reaching the sidecar needs a pathview call, so it is pinned at the source,
    # as the compound_nodes record already is.
    body_src <- paste(deparse(body(generate_multi_ora_pathview)), collapse = " ")
    expect_true(grepl("support_layers = selection$support_layers", body_src, fixed = TRUE))
    expect_true(grepl("contrast_labels = contrast_labels", body_src, fixed = TRUE))
    expect_true(grepl("make.names(names(all_generated_pngs))", body_src, fixed = TRUE))
})

test_that("the pathway selection says whether it fell back to one layer", {
    tbl <- function(support, support_pval = support) {
        data.frame(ID = paste0("P", seq_along(support)), pathway = "x",
                   n_omics_support = support, n_omics_support_pval = support_pval,
                   pooled_pvalue = seq_along(support) / 100)
    }
    # A pathway in two layers: the two-layer rule holds.
    two <- select_multi_ora_pathview_pathways(tbl(c(2, 1)), min_support = 2)
    expect_identical(two$support_layers, 2L)
    expect_identical(two$supported$ID, "P1")
    # Two layers only on raw p: still the two-layer rule.
    raw <- select_multi_ora_pathview_pathways(tbl(c(1, 0), c(2, 1)), min_support = 2)
    expect_identical(raw$support_layers, 2L)
    expect_identical(raw$supported$ID, "P1")
    # No pathway in two layers: the fallback, recorded as one layer.
    one <- suppressMessages(select_multi_ora_pathview_pathways(tbl(c(1, 0)), min_support = 2))
    expect_identical(one$support_layers, 1L)
    expect_identical(one$supported$ID, "P1")
    # Both metabolomics and proteomics under FDR: the paired tier, two layers.
    paired <- cbind(tbl(c(0, 0)), metabolomics_padj = c(0.01, 0.5),
                    proteomics_padj = c(0.02, 0.5))
    res <- suppressMessages(select_multi_ora_pathview_pathways(paired, min_support = 2))
    expect_identical(res$support_layers, 2L)
    expect_identical(res$supported$ID, "P1")
})

# ---- every figure legend within 100 words ----------------------------------

# String constants inside each figure_legend() call. Read off the parser's
# tokens rather than by walking the syntax tree: an argument left empty, as in
# df[i, ], is R's missing value, and passing it to any function -- including the
# walker itself -- raises "argument is missing". Tokens have no such trap, and
# paste()-built legends are still counted whole. Both branches of an if () are
# counted together, which is the worst case. Text a legend gets at render time
# -- a variable or a function call -- is not a constant: its names are recorded
# in `dynamic` so each one can be checked on its own below.
legend_calls <- function(src) {
    starts <- grep("^```\\{r ", src)
    out <- data.frame(text = character(0), dynamic = character(0),
                      stringsAsFactors = FALSE)
    for (s in starts) {
        end <- s + which(grepl("^```\\s*$", src[(s + 1):length(src)]))[1]
        if (is.na(end) || end <= s + 1) next
        pd <- tryCatch(utils::getParseData(parse(text = src[(s + 1):(end - 1)],
                                                 keep.source = TRUE)),
                       error = function(e) NULL)
        if (is.null(pd) || nrow(pd) == 0) next
        pd <- pd[pd$terminal, , drop = FALSE]
        pd <- pd[order(pd$line1, pd$col1), , drop = FALSE]

        i <- 1
        while (i <= nrow(pd)) {
            if (pd$token[i] == "SYMBOL_FUNCTION_CALL" &&
                pd$text[i] == "figure_legend") {
                depth <- 0
                strs <- character(0)
                syms <- character(0)
                j <- i + 1
                while (j <= nrow(pd)) {
                    tk <- pd$token[j]
                    if (tk == "'('") {
                        depth <- depth + 1
                    } else if (tk == "')'") {
                        depth <- depth - 1
                        if (depth == 0) break
                    } else if (tk == "STR_CONST") {
                        strs <- c(strs, pd$text[j])
                    } else if (tk %in% c("SYMBOL", "SYMBOL_FUNCTION_CALL") &&
                               !pd$text[j] %in% c("paste", "paste0", "sprintf")) {
                        syms <- c(syms, pd$text[j])
                    }
                    j <- j + 1
                }
                out <- rbind(out, data.frame(
                    text = paste(gsub('^["\']|["\']$', "", strs), collapse = " "),
                    dynamic = paste(sort(unique(syms)), collapse = ","),
                    stringsAsFactors = FALSE))
                i <- j
            }
            i <- i + 1
        }
    }
    out
}

legend_texts <- function(src) legend_calls(src)$text

n_words <- function(x) lengths(regmatches(x, gregexpr("[^[:space:]]+", x)))

test_that("every figure legend in the report stays within 100 words", {
    texts <- legend_texts(template_lines())
    # The legends are there to be counted; finding none means the parse failed.
    expect_gt(length(texts), 40)
    over <- texts[n_words(texts) > 100]
    expect_identical(length(over), 0L,
                     info = paste(substr(over, 1, 60), collapse = " | "))
})

test_that("the pathway-map legend stays within 100 words with its threshold text", {
    # Read from the template, compound sentence included, so an edit there is
    # counted rather than a copy of the text kept here.
    calls <- legend_calls(template_lines())
    pv <- calls[grepl("KEGG maps drawn with pathview", calls$text, fixed = TRUE), ]
    expect_identical(nrow(pv), 1L)
    expect_match(pv$text, "where available", fixed = TRUE)
    expect_lte(n_words(paste(pv$text, pathview_significance_caption())), 100)
})

test_that("every legend with text from outside its call is checked", {
    # A legend the parse cannot count whole has to be named here with how its
    # length is bounded, so a new one fails until it is covered.
    calls <- legend_calls(template_lines())
    dyn <- calls[nzchar(calls$dynamic), ]
    known <- c(
        # conditions: both branches are already counted as constants
        "single_contrast", "pathview_has_compounds",
        # a MOFA2 view name, a few words at most; see the next test
        "view_name",
        # counted with the pathway-map legend above
        "pathview_significance_caption",
        # worst case checked in test-enzyme-metabolite-pairs.R
        "describe_enzyme_metabolite_table", "enzyme_metabolite_config", "em", "cfg",
        # evaluated from the template below
        "loadings_gsea_legend", "diablo_loadings_caveat")
    used <- unique(unlist(strsplit(dyn$dynamic, ",", fixed = TRUE)))
    expect_identical(setdiff(used, known), character(0))
    # A name-filled legend leaves room for a name of three words.
    named <- dyn[grepl("view_name", dyn$dynamic, fixed = TRUE), ]
    expect_gt(nrow(named), 0)
    expect_true(all(n_words(named$text) <= 97), info = paste(named$text, collapse = " | "))
})

test_that("the loadings legends stay within 100 words with their caveat", {
    # Built from two variables in a setup chunk; evaluated from the template.
    src <- template_lines()
    joined <- paste(src, collapse = "\n")
    grab <- function(name) {
        start <- grep(sprintf("^%s <- paste\\(", name), src)
        expect_length(start, 1)
        end <- start + which(grepl("\\)\\s*$", src[start:length(src)]))[1] - 1
        eval(parse(text = src[start:end]), envir = new.env())
    }
    gsea <- grab("loadings_gsea_legend")
    caveat <- grab("diablo_loadings_caveat")
    expect_lte(n_words(gsea), 100)
    expect_true(grepl("figure_legend(paste(loadings_gsea_legend, diablo_loadings_caveat))",
                      joined, fixed = TRUE))
    expect_lte(n_words(paste(gsea, caveat)), 100)
    calls <- legend_calls(src)
    ora <- calls[grepl("diablo_loadings_caveat", calls$dynamic, fixed = TRUE) &
                 nzchar(calls$text), ]
    expect_identical(nrow(ora), 1L)
    expect_lte(n_words(paste(ora$text, caveat)), 100)
})

test_that("the DIABLO variable loadings plot stays in the report", {
    src <- paste(template_lines(), collapse = "\n")
    expect_true(grepl("diablo_variable_plot.png", src, fixed = TRUE))
    expect_true(grepl("```{r diablo-variable-legend", src, fixed = TRUE))
})

# ---- figures dropped on request --------------------------------------------

test_that("the combined MultiGSEA figure and the DIABLO log2FC loadings are gone", {
    src <- paste(template_lines(), collapse = "\n")
    for (gone in c("multigsea-plots", "multigsea-note",
                   "diablo-loadings-log2fc", "diablo_loadings_log2fc",
                   "Loadings Colored by log2FC",
                   "mofa-loadings-log2fc", "mofa_loadings_log2fc",
                   "Weights Colored by log2FC")) {
        expect_false(grepl(gone, src, fixed = TRUE), info = gone)
    }
    # The per-contrast MultiGSEA tab and the loadings enrichment stay.
    expect_true(grepl("multigsea-per-contrast", src, fixed = TRUE))
    expect_true(grepl("Loadings-Based Pathway Enrichment", src, fixed = TRUE))
})

test_that("the section is named for what is left in it", {
    headings <- static_headings(template_lines())
    expect_true("# MultiGSEA, Multi-ORA and Compound Enrichment {.tabset}" %in% headings)
    expect_false(any(grepl("MultiGSEA Plots|Multi-Omics Enrichment", headings)))
})

test_that("the pairwise MultiGSEA figures show whenever there is no per-contrast set", {
    src <- template_lines()
    joined <- paste(src, collapse = "\n")
    # The view follows what MultiGSEA wrote: it writes per_contrast/ only when
    # its own inputs hold more than one contrast, which a config listing several
    # can fall short of. Without it, the top-level pairs are shown -- under the
    # lone contrast's name in a single-contrast run, as #249 does elsewhere.
    # A figure, not a directory: MultiGSEA creates a contrast's directory before
    # finding it has no pair to draw.
    expect_true(grepl('has_mg_per_contrast <- length(list.files(file.path(multigsea_dir, "per_contrast"),',
                      joined, fixed = TRUE))
    expect_true(grepl('pattern = "^multigsea_.*\\\\.png$",\n                                         recursive = TRUE)) > 0',
                      joined, fixed = TRUE))
    expect_false(grepl('list.dirs(file.path(multigsea_dir, "per_contrast")', joined, fixed = TRUE))
    expect_true(grepl("show_mg_per_contrast <- has_multigsea && has_mg_per_contrast && !single_contrast",
                      joined, fixed = TRUE))
    expect_true(grepl("show_mg_top_pairs <- has_multigsea && !show_mg_per_contrast",
                      joined, fixed = TRUE))
    expect_true(grepl('`r if (show_mg_top_pairs) combined_tab(2, "All Contrasts (Pooled)")`',
                      joined, fixed = TRUE))
    expect_true(grepl('`r if (show_mg_per_contrast) "## Per-Contrast MultiGSEA {.tabset .unnumbered}"`',
                      joined, fixed = TRUE))
    top <- grep("^```\\{r multigsea-top-pairs,", src, value = TRUE)
    expect_length(top, 1)
    expect_true(grepl("eval=show_mg_top_pairs", top, fixed = TRUE))
    pc <- grep("^```\\{r multigsea-per-contrast,", src, value = TRUE)
    expect_true(grepl("eval=show_mg_per_contrast", pc, fixed = TRUE))
    # The same figures the per-contrast tab draws: the pairwise ones only.
    expect_true(grepl('pattern = "^multigsea_.*_vs_.*\\\\.png$"', joined, fixed = TRUE))
    # Drawn and captioned by one helper in both places (calls, not definitions).
    expect_identical(lengths(regmatches(joined, gregexpr("show_multigsea_pair(",
                                                         joined, fixed = TRUE))), 2L)
    expect_identical(lengths(regmatches(joined, gregexpr("multigsea_pair_legend()",
                                                         joined, fixed = TRUE))), 2L)
    expect_false(grepl("Per-contrast MultiGSEA results not available", joined, fixed = TRUE))
})

test_that("the pairwise MultiGSEA legend says what the figure can and cannot show", {
    calls <- legend_calls(template_lines())
    leg <- calls$text[grepl("Two layers compared pathway by pathway", calls$text, fixed = TRUE)]
    expect_length(leg, 1)
    expect_lte(n_words(leg), 100)
    for (phrase in c("Each point is a pathway", "-log10", "adjusted p-value",
                     "sits at zero", "size or colour key", "unscored layer counting as zero",
                     "union of pathways, zeros included",
                     "not a measure of biological agreement")) {
        expect_true(grepl(phrase, leg, fixed = TRUE), info = phrase)
    }
    # The correlation is over the union, not the pathways both layers scored.
    for (phrase in c("overlapping pathways", "intersection")) {
        expect_false(grepl(phrase, leg, fixed = TRUE), info = phrase)
    }
})

test_that("paths pasted into generated chunks use forward slashes", {
    # Child chunks are built with sprintf(), so every image path lands in R
    # source as a string literal; a Windows backslash there is an escape.
    src <- paste(template_lines(), collapse = "\n")
    expect_true(grepl('run_dir <- normalizePath(dirname(knitr::current_input()), winslash = "/", mustWork = TRUE)',
                      src, fixed = TRUE))
})
