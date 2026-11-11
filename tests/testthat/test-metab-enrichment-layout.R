# tests/testthat/test-metab-enrichment-layout.R
#
# The metabolomics enrichment stage writes each result unit to its own
# directory under Enrichment/ (R/domain/metabolomics/06p_enrichment_layout.R).
# These tests cover the layout only: where files go, that their content is the
# table the method returned, that the returned object is unchanged, and that
# nothing outside the stage's own files is removed.
#
# Everything is synthetic. The four enrichment methods are stubbed, because what
# is under test is what the module does with their results, not what they
# compute.

# ---- helpers ----------------------------------------------------------------

# The pipeline functions are sourced into an environment, not exported from a
# package, so with_mocked_bindings() has no namespace to rebind. Assign into the
# module's own environment and restore on exit instead.
local_metab_enrichment_stubs <- function(stubs, env = parent.frame()) {
    target <- environment(mod_metabolomics_enrichment)
    nms <- names(stubs)
    had <- vapply(nms, exists, logical(1), envir = target, inherits = FALSE)
    old <- lapply(nms[had], get, envir = target, inherits = FALSE)
    names(old) <- nms[had]
    withr::defer({
        for (nm in nms) {
            if (nm %in% names(old)) {
                assign(nm, old[[nm]], envir = target)
            } else if (exists(nm, envir = target, inherits = FALSE)) {
                rm(list = nm, envir = target)
            }
        }
    }, envir = env)
    for (nm in nms) assign(nm, stubs[[nm]], envir = target)
    invisible(NULL)
}

layout_fixture <- function(contrasts = c("A vs B", "C/D")) {
    samples <- paste0("S", 1:8)
    groups  <- rep(c("A", "B", "C", "D"), each = 2)
    mat <- matrix(seq_len(3 * 8) / 10, nrow = 3,
                  dimnames = list(paste0("F", 1:3), samples))
    num <- c("A", "C")[seq_along(contrasts)]
    den <- c("B", "D")[seq_along(contrasts)]
    list(
        pre = list(expr_raw = mat, expr_filt = mat, expr_work = mat,
                   meta = data.frame(sample_id = samples, sample_type = groups,
                                     stringsAsFactors = FALSE),
                   row_data = data.frame(feature_id = rownames(mat),
                                         stringsAsFactors = FALSE)),
        de_res = list(de_tables = stats::setNames(
            lapply(contrasts, function(x) data.frame(feature_id = rownames(mat))),
            contrasts)),
        inputs = list(contrasts = data.frame(Contrast_name = contrasts,
                                             Numerator = num, Denominator = den,
                                             stringsAsFactors = FALSE)),
        config = list(modes = list(metabolomics = list(
            enrichment = list(run_enrichment = TRUE))))
    )
}

# One synthetic table per call, tagged so a table written to the wrong unit is
# visible. FDR values straddle 0.05 on purpose.
stub_table <- function(tag, fdr = c(0.01, 0.2, NA)) {
    data.frame(pathway = paste0("P", seq_along(fdr)), tag = tag, FDR = fdr,
               stringsAsFactors = FALSE)
}

method_stubs <- function(plot = function(...) NULL) {
    list(
        run_metabolomics_qea = function(pre, config)
            list(table = stub_table(paste0("qea_n", ncol(pre$expr_raw)))),
        run_metabolomics_ssgsea = function(pre, config, contrast_groups = NULL,
                                           contrast_column = NULL)
            list(table = stub_table(paste0("ssgsea_n", ncol(pre$expr_raw)),
                                    fdr = c(0.5, 0.6)),
                 scores = NULL),
        run_metabolomics_ora = function(pre, de_res, config, contrast = NULL)
            list(table = stub_table(paste0("ora_", contrast))),
        run_metabolomics_gsea = function(pre, de_res, config, contrast = NULL)
            list(table = stub_table(paste0("gsea_", contrast), fdr = c(0.3, 0.4))),
        plot_enrichment_barplot = plot, plot_qea_lollipop = plot,
        plot_ssgsea_lollipop = plot, plot_ssgsea_boxplots = plot,
        plot_ora_dotplot = plot, plot_ora_lollipop = plot,
        plot_gsea_nes_barplot = plot, plot_gsea_lollipop = plot
    )
}

rel_files <- function(root) sort(list.files(root, recursive = TRUE))

expect_tsv_holds <- function(path, table) {
    ref <- tempfile(fileext = ".tsv")
    readr::write_tsv(table, ref)
    expect_identical(readLines(path), readLines(ref))
}

# ---- path builders ----------------------------------------------------------

test_that("each method has its unit directory, with the context in the path", {
    expect_identical(metab_enrichment_unit_rel_dir("ora", "A_vs_B"),
                     "ORA/combined/contrasts/A_vs_B/de")
    expect_identical(metab_enrichment_unit_rel_dir("gsea", "A_vs_B"),
                     "GSEA/combined/ranking_by_fc/A_vs_B")
    expect_identical(metab_enrichment_unit_rel_dir("qea"),
                     "QEA/combined/all_samples")
    expect_identical(metab_enrichment_unit_rel_dir("qea", "A_vs_B"),
                     "QEA/combined/contrasts/A_vs_B")
    expect_identical(metab_enrichment_unit_rel_dir("ssgsea"),
                     "ssGSEA/combined/all_samples")
    expect_identical(metab_enrichment_unit_rel_dir("ssgsea", "A_vs_B"),
                     "ssGSEA/combined/contrasts/A_vs_B")
    expect_identical(metab_enrichment_unit_dir("out/Enrichment", "qea"),
                     file.path("out/Enrichment", "QEA/combined/all_samples"))
})

test_that("a unit directory cannot be built from an unsafe or missing name", {
    expect_error(metab_enrichment_unit_rel_dir("mummichog"), "unknown method")
    expect_error(metab_enrichment_unit_rel_dir("ora"), "contrast_dir")
    expect_error(metab_enrichment_unit_rel_dir("gsea"), "contrast_dir")
    for (bad in c("", ".", "..", "a/b", "a\\b", NA_character_)) {
        expect_error(metab_enrichment_unit_rel_dir("qea", bad), "not a usable")
    }
})

test_that("contrast directory names are distinct after sanitisation", {
    labels <- c("A vs B", "A/vs/B", "a_VS_b", "C vs D")
    tokens <- c("A_vs_B", "A_vs_B", "a_VS_b", "C_vs_D")
    expect_warning(dirs <- metab_enrichment_contrast_dirs(labels, tokens),
                   "collide")

    expect_identical(names(dirs), labels)
    expect_identical(unname(dirs[c("A vs B", "C vs D")]), c("A_vs_B", "C_vs_D"))
    # Distinct even on a case-insensitive file system.
    expect_identical(anyDuplicated(tolower(dirs)), 0L)
    # Every name is accepted by the path builder.
    for (d in dirs) expect_silent(metab_enrichment_unit_rel_dir("ora", d))
})

test_that("a suffix given to one contrast cannot collide with another", {
    # "A_2" is taken by a real contrast, so the second "A" must not become "A_2".
    expect_warning(
        dirs <- metab_enrichment_contrast_dirs(c("A", "A_2", "A "),
                                               c("A", "A_2", "A")),
        "collide")
    expect_identical(anyDuplicated(tolower(dirs)), 0L)
})

test_that("names that are no directory of their own are replaced", {
    dirs <- suppressWarnings(
        metab_enrichment_contrast_dirs(c("x", "y", "z", "w"),
                                       c("", "..", "A.", NA)))
    expect_false(any(dirs %in% c("", ".", "..")))
    expect_false(any(grepl("\\.$", dirs)))
    expect_identical(anyDuplicated(tolower(dirs)), 0L)
})

test_that("a contrast named like a Windows device gets a usable directory name", {
    labels <- c("CON", "nul", "Com1", "LPT9.txt", "aux.", "CONTROL", "COM10", "A_CON")
    expect_silent(dirs <- metab_enrichment_contrast_dirs(labels, labels))

    expect_identical(unname(dirs[c("CON", "nul", "Com1", "LPT9.txt", "aux.")]),
                     c("contrast_CON", "contrast_nul", "contrast_Com1",
                       "contrast_LPT9.txt", "contrast_aux"))
    # Names that merely contain or start with one are left alone.
    expect_identical(unname(dirs[c("CONTROL", "COM10", "A_CON")]),
                     c("CONTROL", "COM10", "A_CON"))
    reserved <- "^(con|prn|aux|nul|com[1-9]|lpt[1-9])(\\..*)?$"
    expect_false(any(grepl(reserved, dirs, ignore.case = TRUE)))
    for (d in dirs) expect_silent(metab_enrichment_unit_rel_dir("gsea", d))
})

test_that("clean contrast names pass through silently and unchanged", {
    expect_silent(dirs <- metab_enrichment_contrast_dirs(c("A vs B", "C vs D"),
                                                         c("A_vs_B", "C_vs_D")))
    expect_identical(unname(dirs), c("A_vs_B", "C_vs_D"))
    expect_length(metab_enrichment_contrast_dirs(character(0), character(0)), 0L)
    expect_error(metab_enrichment_contrast_dirs("A", c("A", "B")), "same length")
})

# ---- flat (pre-layout) names -------------------------------------------------

test_that("every unit plot maps back to the flat name it had before the layout", {
    root <- file.path("res", "metabolomics", "Enrichment")
    unit <- function(method, contrast = NULL)
        metab_enrichment_unit_dir(root, method, contrast)

    new <- c(
        file.path(unit("qea"), c("barplot.png", "lollipop.png", "library_kegg.png")),
        file.path(unit("ssgsea"), c("barplot.png", "boxplots.png", "lollipop.png")),
        file.path(unit("ora", "A_vs_B"), c("dotplot.png", "lollipop.png", "barplot.png")),
        file.path(unit("gsea", "A_vs_B"), c("nes_barplot.png", "lollipop.png", "barplot.png")),
        file.path(unit("qea", "A_vs_B"), c("lollipop.png", "barplot.png")),
        file.path(unit("ssgsea", "A_vs_B"), c("lollipop.png", "barplot.png")),
        # the first-contrast copies stay flat and keep their own names
        file.path(root, c("enrichment_ora_dotplot.png", "enrichment_gsea_barplot.png"))
    )
    old <- c(
        "enrichment_qea_barplot.png", "enrichment_qea_lollipop.png",
        "enrichment_qea_kegg.png",
        "enrichment_ssgsea_barplot.png", "enrichment_ssgsea_boxplots.png",
        "enrichment_ssgsea_lollipop.png",
        "enrichment_ora_dotplot_A_vs_B.png", "enrichment_ora_lollipop_A_vs_B.png",
        "enrichment_ora_barplot_A_vs_B.png",
        "enrichment_gsea_nes_barplot_A_vs_B.png", "enrichment_gsea_lollipop_A_vs_B.png",
        "enrichment_gsea_barplot_A_vs_B.png",
        "enrichment_qea_lollipop_A_vs_B.png", "enrichment_qea_barplot_A_vs_B.png",
        "enrichment_ssgsea_lollipop_A_vs_B.png", "enrichment_ssgsea_barplot_A_vs_B.png",
        "enrichment_ora_dotplot.png", "enrichment_gsea_barplot.png"
    )
    expect_identical(metab_enrichment_legacy_plot_name(new, root), old)
    # Windows separators give the same answer.
    expect_identical(
        metab_enrichment_legacy_plot_name(gsub("/", "\\\\", new), root), old)
    # The flat names are distinct, so ordering by them is a total order.
    expect_identical(anyDuplicated(old), 0L)
})

test_that("a file outside Enrichment/ keeps its own name", {
    expect_identical(
        metab_enrichment_legacy_plot_name("res/metabolomics/Network/net.png",
                                          "res/metabolomics/Enrichment"),
        "net.png")
    expect_identical(
        metab_enrichment_legacy_plot_name(character(0), "res/Enrichment"),
        character(0))
})

# ---- manifest ----------------------------------------------------------------

test_that("the manifest has the RNA-seq manifest columns and one row per unit", {
    units <- list(
        list(method = "qea", contrast = NA_character_,
             rel_dir = "QEA/combined/all_samples",
             table = data.frame(FDR = c(0.01, 0.04, 0.05, NA))),
        list(method = "ora", contrast = "A vs B",
             rel_dir = "ORA/combined/contrasts/A_vs_B/de",
             table = data.frame(FDR = c(0.2, 0.9))),
        list(method = "gsea", contrast = "A vs B",
             rel_dir = "GSEA/combined/ranking_by_fc/A_vs_B",
             table = data.frame(FDR = 0.001)),
        list(method = "ssgsea", contrast = "A vs B",
             rel_dir = "ssGSEA/combined/contrasts/A_vs_B",
             table = data.frame(FDR = c(NA_real_, NA_real_)))
    )
    m <- build_metab_enrichment_manifest(units)

    expect_identical(names(m), names(.empty_enrichment_manifest()))
    expect_identical(m$analysis, c("QEA", "ORA", "GSEA", "ssGSEA"))
    expect_true(all(m$database == "combined"))
    expect_identical(m$group, c("all_samples", "contrasts", "fc", "contrasts"))
    # The label, not the directory name, identifies the contrast.
    expect_identical(m$item, c("all_samples", "A vs B", "A vs B", "A vs B"))
    expect_identical(m$subgroup, c(NA, "all", NA, NA))
    expect_identical(m$group_type, c(NA, "direction", NA, NA))
    # Strictly below 0.05; a missing FDR is not significant.
    expect_identical(m$n_significant, c(2L, 0L, 1L, 0L))
    expect_identical(m$status, c("significant", "empty", "significant", "empty"))
    expect_true(all(m$evaluated))
    expect_false(any(m$has_simplify))
    expect_identical(m$storage_key, vapply(units, `[[`, character(1), "rel_dir"))
})

test_that("a run with no units gives the empty manifest", {
    m <- build_metab_enrichment_manifest(list())
    expect_identical(m, .empty_enrichment_manifest())
})

# ---- the module: what it writes, and what it returns ---------------------------

test_that("the module writes one unit directory per result, plus the manifest", {
    fx  <- layout_fixture()
    out <- file.path(withr::local_tempdir(), "metabolomics")
    local_metab_enrichment_stubs(method_stubs())

    res <- expect_silent(mod_metabolomics_enrichment(
        fx$pre, fx$de_res, fx$config, out, inputs = fx$inputs))

    expect_identical(rel_files(file.path(out, "Enrichment")), sort(c(
        "enrichment_manifest.csv",
        "QEA/combined/all_samples/results.tsv",
        "QEA/combined/contrasts/A_vs_B/results.tsv",
        "QEA/combined/contrasts/C_D/results.tsv",
        "ssGSEA/combined/all_samples/results.tsv",
        "ssGSEA/combined/contrasts/A_vs_B/results.tsv",
        "ssGSEA/combined/contrasts/C_D/results.tsv",
        "ORA/combined/contrasts/A_vs_B/de/results.tsv",
        "ORA/combined/contrasts/C_D/de/results.tsv",
        "GSEA/combined/ranking_by_fc/A_vs_B/results.tsv",
        "GSEA/combined/ranking_by_fc/C_D/results.tsv")))
    # Only the first-contrast compatibility copies remain in Datasets/.
    expect_identical(rel_files(file.path(out, "Datasets")),
                     c("enrichment_gsea_results.tsv", "enrichment_ora_results.tsv"))

    # Every file the module reports exists, and it reports every file it wrote.
    expect_true(all(file.exists(res$files)))
    expect_setequal(
        normalizePath(res$files),
        normalizePath(c(
            list.files(file.path(out, "Enrichment"), recursive = TRUE, full.names = TRUE),
            list.files(file.path(out, "Datasets"), recursive = TRUE, full.names = TRUE))))
})

test_that("each results.tsv holds exactly the table its method returned", {
    fx  <- layout_fixture()
    out <- file.path(withr::local_tempdir(), "metabolomics")
    local_metab_enrichment_stubs(method_stubs())
    res <- mod_metabolomics_enrichment(fx$pre, fx$de_res, fx$config, out,
                                       inputs = fx$inputs)
    enr <- file.path(out, "Enrichment")
    at  <- function(method, contrast = NULL)
        file.path(metab_enrichment_unit_dir(enr, method, contrast), "results.tsv")

    expect_tsv_holds(at("qea"), res$qea$table)
    expect_tsv_holds(at("ssgsea"), res$ssgsea$table)
    dirs <- c("A vs B" = "A_vs_B", "C/D" = "C_D")
    for (ct in names(dirs)) {
        expect_tsv_holds(at("ora", dirs[[ct]]), res$ora_by_contrast[[ct]]$res$table)
        expect_tsv_holds(at("gsea", dirs[[ct]]), res$gsea_by_contrast[[ct]]$res$table)
        expect_tsv_holds(at("qea", dirs[[ct]]), res$qea_by_contrast[[ct]]$res$table)
        expect_tsv_holds(at("ssgsea", dirs[[ct]]), res$ssgsea_by_contrast[[ct]]$res$table)
    }
    # The compatibility copies are still the FIRST contrast, as before.
    expect_tsv_holds(file.path(out, "Datasets", "enrichment_ora_results.tsv"),
                     res$ora_by_contrast[["A vs B"]]$res$table)
    expect_tsv_holds(file.path(out, "Datasets", "enrichment_gsea_results.tsv"),
                     res$gsea_by_contrast[["A vs B"]]$res$table)
})

test_that("the returned object keeps its fields and routes each contrast's result", {
    fx  <- layout_fixture()
    out <- file.path(withr::local_tempdir(), "metabolomics")
    local_metab_enrichment_stubs(method_stubs())
    res <- mod_metabolomics_enrichment(fx$pre, fx$de_res, fx$config, out,
                                       inputs = fx$inputs)

    expect_identical(names(res),
                     c("qea", "ssgsea", "ora", "gsea", "qea_by_contrast",
                       "ssgsea_by_contrast", "ora_by_contrast", "gsea_by_contrast",
                       "plots", "files"))
    for (by in c("qea_by_contrast", "ssgsea_by_contrast", "ora_by_contrast",
                 "gsea_by_contrast")) {
        expect_identical(names(res[[by]]), c("A vs B", "C/D"))
        for (ct in names(res[[by]]))
            expect_identical(names(res[[by]][[ct]]), c("res", "plots", "files"))
    }
    # The global runs see all samples, the per-contrast runs their two groups.
    expect_identical(unique(res$qea$table$tag), "qea_n8")
    expect_identical(unique(res$qea_by_contrast[["C/D"]]$res$table$tag), "qea_n4")
    expect_identical(unique(res$ora_by_contrast[["C/D"]]$res$table$tag), "ora_C/D")
    # `ora` / `gsea` are still the first contrast.
    expect_identical(res$ora, res$ora_by_contrast[["A vs B"]]$res)
    expect_identical(res$gsea, res$gsea_by_contrast[["A vs B"]]$res)
})

test_that("the manifest on disk lists the units of this run", {
    fx  <- layout_fixture()
    out <- file.path(withr::local_tempdir(), "metabolomics")
    local_metab_enrichment_stubs(method_stubs())
    mod_metabolomics_enrichment(fx$pre, fx$de_res, fx$config, out,
                                inputs = fx$inputs)
    enr <- file.path(out, "Enrichment")
    m <- utils::read.csv(file.path(enr, "enrichment_manifest.csv"),
                         stringsAsFactors = FALSE)

    expect_identical(names(m), names(.empty_enrichment_manifest()))
    expect_identical(nrow(m), 10L)
    # Each listed unit is a directory holding its table.
    expect_true(all(file.exists(file.path(enr, m$storage_key, "results.tsv"))))
    expect_setequal(m$item, c("all_samples", "A vs B", "C/D"))
    # The stub tables: ORA and QEA have one pathway under 0.05, the others none.
    expect_true(all(m$n_significant[m$analysis %in% c("ORA", "QEA")] == 1L))
    expect_true(all(m$n_significant[m$analysis %in% c("GSEA", "ssGSEA")] == 0L))
})

test_that("contrasts that sanitise to one name keep separate outputs", {
    fx  <- layout_fixture(contrasts = c("A vs B", "A/vs/B"))
    out <- file.path(withr::local_tempdir(), "metabolomics")
    local_metab_enrichment_stubs(method_stubs())

    expect_warning(
        res <- mod_metabolomics_enrichment(fx$pre, fx$de_res, fx$config, out,
                                           inputs = fx$inputs),
        "collide")
    enr <- file.path(out, "Enrichment")
    first  <- file.path(enr, "ORA/combined/contrasts/A_vs_B/de/results.tsv")
    second <- file.path(enr, "ORA/combined/contrasts/A_vs_B_2/de/results.tsv")
    expect_tsv_holds(first,  res$ora_by_contrast[["A vs B"]]$res$table)
    expect_tsv_holds(second, res$ora_by_contrast[["A/vs/B"]]$res$table)
})

# ---- what the module leaves alone -----------------------------------------------

# An output folder as an earlier run could have left it. `stale` are files this
# stage owns (by path shape AND name) that the stubbed run below does not
# write; `foreign` are files it does not own. Most of the foreign ones are
# deliberate look-alikes: a PNG inside a method folder, a file with the stage's
# prefix but not one of its plot names, a stage file name in the wrong place.
seed_earlier_run <- function(out) {
    enr <- file.path(out, "Enrichment")
    old_unit <- file.path(enr, "ORA", "combined", "contrasts", "Old", "de")
    cur_unit <- file.path(enr, "ORA", "combined", "contrasts", "A_vs_B", "de")
    qea_all  <- file.path(enr, "QEA", "combined", "all_samples")
    for (d in c(old_unit, cur_unit, qea_all,
                file.path(enr, "GSEA", "combined", "ranking_by_fc", "Old", "extra"),
                file.path(out, "Datasets"),
                file.path(out, "mummichog_pinned", "Old"))) {
        dir.create(d, recursive = TRUE)
    }

    stale <- c(
        # flat files from before the unit layout
        file.path(enr, "enrichment_ora_dotplot_Old.png"),
        file.path(enr, "enrichment_qea_barplot.png"),
        file.path(enr, "enrichment_qea_kegg.png"),
        file.path(out, "Datasets", "enrichment_ora_Old_results.tsv"),
        file.path(out, "Datasets", "enrichment_qea_results.tsv"),
        # a first-contrast copy this run does not redraw
        file.path(enr, "enrichment_ora_dotplot.png"),
        # a unit of a contrast that is gone
        file.path(old_unit, "results.tsv"),
        file.path(old_unit, "dotplot.png"),
        # plots in units this run rewrites, but does not redraw
        file.path(cur_unit, "dotplot.png"),
        file.path(qea_all, "library_old.png"))
    foreign <- c(
        file.path(out, "Datasets", "final_results.tsv"),
        file.path(out, "Datasets", "enrichment_summary.tsv"),
        file.path(enr, "notes.txt"),
        file.path(enr, "other.png"),
        file.path(enr, "mummichog_pathway_bubble_Old.png"),
        # the stage's prefix, but not one of its plot names
        file.path(enr, "enrichment_ora_myfigure.png"),
        # PNGs someone put inside the method folders
        file.path(enr, "ORA", "overview.png"),
        file.path(enr, "ORA", "combined", "summary.png"),
        file.path(old_unit, "my_figure.png"),
        file.path(old_unit, "notes.txt"),
        file.path(cur_unit, "custom.png"),
        # a GSEA plot name in an ORA unit, and a stage name one level too deep
        file.path(cur_unit, "nes_barplot.png"),
        file.path(enr, "GSEA", "combined", "ranking_by_fc", "Old", "extra", "barplot.png"),
        # a library plot without the prefix, and a results table that is not .tsv
        file.path(qea_all, "kegg.png"),
        file.path(qea_all, "results.csv"),
        file.path(out, "mummichog_pinned", "Old", "input.tsv"),
        file.path(out, "mummichog_pinned", "enrichment_ora_dotplot_Old.png"))
    file.create(c(stale, foreign))
    manifest <- file.path(enr, "enrichment_manifest.csv")
    writeLines("earlier run", manifest)
    list(stale = stale, foreign = foreign, manifest = manifest)
}

test_that("only the stage's own files are listed as owned", {
    out <- file.path(withr::local_tempdir(), "metabolomics")
    old <- seed_earlier_run(out)
    owned <- list_metab_enrichment_owned_files(file.path(out, "Enrichment"),
                                               file.path(out, "Datasets"))
    expect_setequal(normalizePath(owned), normalizePath(old$stale))
    # The manifest is handled separately and is never in the list.
    expect_false(any(basename(owned) == "enrichment_manifest.csv"))
    expect_length(list_metab_enrichment_owned_files(
        file.path(out, "no_such_dir"), file.path(out, "nor_this")), 0L)
})

test_that("a disabled run deletes nothing, not even the earlier manifest", {
    fx  <- layout_fixture()
    out <- file.path(withr::local_tempdir(), "metabolomics")
    old <- seed_earlier_run(out)
    fx$config$modes$metabolomics$enrichment$run_enrichment <- FALSE
    before <- rel_files(out)

    expect_null(suppressMessages(mod_metabolomics_enrichment(
        fx$pre, fx$de_res, fx$config, out, inputs = fx$inputs)))

    expect_identical(rel_files(out), before)
    expect_identical(readLines(old$manifest), "earlier run")
})

test_that("a completed run removes its own earlier files and nothing else", {
    fx  <- layout_fixture()
    out <- file.path(withr::local_tempdir(), "metabolomics")
    old <- seed_earlier_run(out)
    local_metab_enrichment_stubs(method_stubs())

    expect_message(
        res <- mod_metabolomics_enrichment(fx$pre, fx$de_res, fx$config, out,
                                           inputs = fx$inputs),
        "removed 10 file")

    expect_false(any(file.exists(old$stale)))
    expect_true(all(file.exists(old$foreign)))
    # What this run wrote is all still there.
    expect_true(all(file.exists(res$files)))

    m <- utils::read.csv(old$manifest, stringsAsFactors = FALSE)
    # It is this run's manifest, not the seeded one-line file.
    expect_identical(names(m), names(.empty_enrichment_manifest()))
    expect_identical(nrow(m), 10L)
    expect_false("Old" %in% m$item)
    # Every enrichment table left on disk is one the manifest accounts for,
    # apart from the two first-contrast copies.
    enr <- file.path(out, "Enrichment")
    expect_setequal(
        list.files(enr, pattern = "^results\\.tsv$", recursive = TRUE),
        file.path(m$storage_key, "results.tsv"))
})

test_that("a clean run has nothing to remove and says nothing", {
    fx  <- layout_fixture()
    out <- file.path(withr::local_tempdir(), "metabolomics")
    local_metab_enrichment_stubs(method_stubs())
    expect_silent(mod_metabolomics_enrichment(fx$pre, fx$de_res, fx$config, out,
                                              inputs = fx$inputs))
    # Running again over its own output removes nothing either.
    expect_silent(res <- mod_metabolomics_enrichment(fx$pre, fx$de_res, fx$config,
                                                     out, inputs = fx$inputs))
    expect_true(all(file.exists(res$files)))
})

test_that("a PNG that is not the stage's own survives inside a method folder", {
    fx  <- layout_fixture()
    out <- file.path(withr::local_tempdir(), "metabolomics")
    local_metab_enrichment_stubs(method_stubs())
    # A first run creates the real unit folders ...
    mod_metabolomics_enrichment(fx$pre, fx$de_res, fx$config, out,
                                inputs = fx$inputs)
    enr <- file.path(out, "Enrichment")
    mine <- c(
        file.path(enr, "ORA", "combined", "contrasts", "A_vs_B", "de", "figure_for_paper.png"),
        file.path(enr, "GSEA", "combined", "ranking_by_fc", "C_D", "annotated.png"),
        file.path(enr, "QEA", "combined", "all_samples", "kegg.png"),
        file.path(enr, "ssGSEA", "combined", "contrasts", "A_vs_B", "scores_by_hand.png"),
        file.path(enr, "ssGSEA", "panel.png"))
    file.create(mine)

    # ... and a second run over them leaves every such file in place.
    expect_silent(mod_metabolomics_enrichment(fx$pre, fx$de_res, fx$config, out,
                                              inputs = fx$inputs))
    expect_true(all(file.exists(mine)))
    expect_length(intersect(
        normalizePath(mine),
        normalizePath(list_metab_enrichment_owned_files(
            enr, file.path(out, "Datasets")))), 0L)
})

test_that("nothing is removed when the manifest cannot be written", {
    fx  <- layout_fixture()
    out <- file.path(withr::local_tempdir(), "metabolomics")
    old <- seed_earlier_run(out)
    stubs <- method_stubs()
    stubs$.write_metab_enrichment_manifest <- function(...) stop("manifest failed")
    local_metab_enrichment_stubs(stubs)

    expect_error(mod_metabolomics_enrichment(fx$pre, fx$de_res, fx$config, out,
                                             inputs = fx$inputs),
                 "manifest failed")

    expect_true(all(file.exists(c(old$stale, old$foreign))))
    expect_false(file.exists(old$manifest))
})

test_that("a file that could not be deleted is reported as left, not as removed", {
    fx  <- layout_fixture()
    out <- file.path(withr::local_tempdir(), "metabolomics")
    old <- seed_earlier_run(out)
    stubs <- method_stubs()
    # Deletion that silently does nothing, as unlink() does on a locked file.
    # The earlier manifest is not locked: that case stops the run (see below).
    stubs$.delete_metab_enrichment_files <- function(paths) {
        is_manifest <- grepl("enrichment_manifest\\.csv", paths)
        unlink(paths[is_manifest])
        is_manifest
    }
    local_metab_enrichment_stubs(stubs)

    msgs <- testthat::capture_messages(
        expect_warning(
            res <- mod_metabolomics_enrichment(fx$pre, fx$de_res, fx$config, out,
                                               inputs = fx$inputs),
            "10 file\\(s\\) left by an earlier run could not be removed"))

    expect_false(any(grepl("removed", msgs)))
    expect_true(all(file.exists(old$stale)))
    # The run itself is complete: its results and its manifest are there.
    expect_true(all(file.exists(res$files)))
    m <- utils::read.csv(old$manifest, stringsAsFactors = FALSE)
    expect_identical(nrow(m), 10L)
})

test_that("a run stops before writing anything when the earlier manifest is stuck", {
    fx  <- layout_fixture()
    out <- file.path(withr::local_tempdir(), "metabolomics")
    old <- seed_earlier_run(out)
    before <- rel_files(out)
    stubs <- method_stubs()
    # The earlier manifest is open elsewhere: unlink() leaves it and says nothing.
    stubs$.delete_metab_enrichment_files <- function(paths) rep(FALSE, length(paths))
    local_metab_enrichment_stubs(stubs)

    expect_error(mod_metabolomics_enrichment(fx$pre, fx$de_res, fx$config, out,
                                             inputs = fx$inputs),
                 "manifest of an earlier run could not be removed")

    # No result of the new run exists, and the earlier folder is as it was, so
    # the earlier manifest still describes what is on disk.
    expect_identical(rel_files(out), before)
    expect_identical(readLines(old$manifest), "earlier run")
})

test_that("a hidden unit directory is still part of the stage's files", {
    out <- file.path(withr::local_tempdir(), "metabolomics")
    enr <- file.path(out, "Enrichment")
    unit <- file.path(enr, "ORA", "combined", "contrasts", ".treated", "de")
    dir.create(unit, recursive = TRUE)
    mine <- file.path(unit, c("results.tsv", "dotplot.png"))
    file.create(c(mine, file.path(unit, ".notes"), file.path(enr, ".DS_Store")))

    owned <- list_metab_enrichment_owned_files(enr, file.path(out, "Datasets"))
    expect_setequal(normalizePath(owned), normalizePath(mine))
})

test_that("only the files that are really gone are returned as removed", {
    d <- withr::local_tempdir()
    enr <- file.path(d, "Enrichment"); ds <- file.path(d, "Datasets")
    dir.create(enr); dir.create(ds)
    f <- file.path(enr, c("enrichment_ora_dotplot_A.png", "enrichment_ora_dotplot_B.png"))
    file.create(f)
    local_metab_enrichment_stubs(list(
        .delete_metab_enrichment_files = function(paths) {
            gone <- grepl("_A\\.png$", paths)
            unlink(paths[gone])
            gone
        }))

    expect_message(
        expect_warning(
            removed <- .remove_stale_metab_enrichment_files(enr, ds, written = character(0)),
            "1 file\\(s\\) left by an earlier run could not be removed"),
        "removed 1 file")
    expect_identical(basename(removed), "enrichment_ora_dotplot_A.png")
})

test_that("deleting files reports what is actually gone", {
    d <- withr::local_tempdir()
    f <- file.path(d, c("a.png", "b.png"))
    file.create(f)
    expect_identical(.delete_metab_enrichment_files(f), c(TRUE, TRUE))
    expect_false(any(file.exists(f)))
    expect_identical(.delete_metab_enrichment_files(character(0)), logical(0))
    # unlink() does not delete a directory unless asked to, and says nothing.
    sub <- file.path(d, "kept_dir")
    dir.create(sub)
    expect_identical(.delete_metab_enrichment_files(sub), FALSE)
    expect_true(dir.exists(sub))
})

test_that("a run that produces no result removes nothing", {
    fx  <- layout_fixture()
    out <- file.path(withr::local_tempdir(), "metabolomics")
    old <- seed_earlier_run(out)
    stubs <- method_stubs()
    for (m in c("run_metabolomics_qea", "run_metabolomics_ssgsea",
                "run_metabolomics_ora", "run_metabolomics_gsea")) {
        stubs[[m]] <- function(...) NULL
    }
    local_metab_enrichment_stubs(stubs)

    expect_message(
        res <- mod_metabolomics_enrichment(fx$pre, fx$de_res, fx$config, out,
                                           inputs = fx$inputs),
        "produced no result")

    expect_true(all(file.exists(c(old$stale, old$foreign))))
    # The manifest is this run's: complete, and empty.
    m <- utils::read.csv(old$manifest, stringsAsFactors = FALSE)
    expect_identical(names(m), names(.empty_enrichment_manifest()))
    expect_identical(nrow(m), 0L)
    expect_identical(basename(res$files), "enrichment_manifest.csv")
})

test_that("a run that stops part-way removes no results and leaves no manifest", {
    fx  <- layout_fixture()
    out <- file.path(withr::local_tempdir(), "metabolomics")
    old <- seed_earlier_run(out)
    stubs <- method_stubs()
    stubs$.save_legacy_enrichment <- function(...) stop("stopped part-way")
    local_metab_enrichment_stubs(stubs)

    expect_error(mod_metabolomics_enrichment(fx$pre, fx$de_res, fx$config, out,
                                             inputs = fx$inputs),
                 "stopped part-way")

    expect_false(file.exists(old$manifest))
    expect_false(file.exists(paste0(old$manifest, ".tmp")))
    # The earlier results it did not get to replace are still there.
    expect_true(all(file.exists(old$stale[!grepl("A_vs_B", old$stale)])))
    expect_true(all(file.exists(old$foreign)))
})

# ---- manifest write ---------------------------------------------------------------

test_that("the manifest is written whole, with no temporary file left", {
    path <- file.path(withr::local_tempdir(), "enrichment_manifest.csv")
    m <- build_metab_enrichment_manifest(list(
        list(method = "qea", contrast = NA_character_,
             rel_dir = "QEA/combined/all_samples",
             table = data.frame(FDR = 0.01))))

    expect_identical(.write_metab_enrichment_manifest(m, path), path)
    expect_false(file.exists(paste0(path, ".tmp")))
    back <- utils::read.csv(path, stringsAsFactors = FALSE)
    expect_identical(names(back), names(m))
    expect_identical(back$storage_key, m$storage_key)
})

test_that("a failed manifest write leaves the target as it was and no partial file", {
    path <- file.path(withr::local_tempdir(), "enrichment_manifest.csv")
    writeLines("complete earlier manifest", path)
    # write.csv() cannot encode a list column, so the write fails after the
    # temporary file has been opened.
    bad <- data.frame(analysis = c("ORA", "GSEA"))
    bad$broken <- list(1, 2)

    expect_error(.write_metab_enrichment_manifest(bad, path))
    expect_identical(readLines(path), "complete earlier manifest")
    expect_false(file.exists(paste0(path, ".tmp")))
})

# ---- plots ----------------------------------------------------------------------

# ---- a plot counts as written only when the file is really there ------------------
#
# These tests replace the one call that draws to a graphics device, so they do
# not depend on how a device reports a failure (an error on some platforms, a
# warning or nothing at all on others) or on a device being available.

# Stand-ins for .draw_metab_enrichment_png(plot, path, width, height).
draw_ok      <- function(plot, path, width, height) writeLines("new plot", path)
draw_nothing <- function(plot, path, width, height) invisible(NULL)
draw_warns   <- function(plot, path, width, height) {
    warning("could not open file")
    invisible(NULL)
}
draw_empty   <- function(plot, path, width, height) file.create(path)

temp_plot_files <- function(dir) {
    list.files(dir, pattern = "^\\.saving_", all.files = TRUE, recursive = TRUE)
}

test_that("a save that produced a file reports it and leaves no temporary file", {
    d <- withr::local_tempdir()
    target <- file.path(d, "dotplot.png")
    local_metab_enrichment_stubs(list(.draw_metab_enrichment_png = draw_ok))

    expect_identical(.save_metab_enrichment_plot("a plot", target), target)
    expect_identical(readLines(target), "new plot")
    expect_length(temp_plot_files(d), 0L)

    # Saving again replaces the file in place.
    local_metab_enrichment_stubs(list(.draw_metab_enrichment_png =
        function(plot, path, width, height) writeLines("newer plot", path)))
    .save_metab_enrichment_plot("a plot", target)
    expect_identical(readLines(target), "newer plot")
})

test_that("a device that writes nothing is a failed save, whatever it signals", {
    d <- withr::local_tempdir()
    target <- file.path(d, "dotplot.png")
    writeLines("earlier run's plot", target)

    # Silent: no error, no warning, no file.
    local_metab_enrichment_stubs(list(.draw_metab_enrichment_png = draw_nothing))
    expect_error(.save_metab_enrichment_plot("a plot", target), "was not written")

    # A warning instead of an error.
    local_metab_enrichment_stubs(list(.draw_metab_enrichment_png = draw_warns))
    expect_warning(
        expect_error(.save_metab_enrichment_plot("a plot", target), "was not written"),
        "could not open file")

    # A file that was opened but never filled.
    local_metab_enrichment_stubs(list(.draw_metab_enrichment_png = draw_empty))
    expect_error(.save_metab_enrichment_plot("a plot", target), "was not written")

    # In every case the earlier file is untouched and nothing is left behind.
    expect_identical(readLines(target), "earlier run's plot")
    expect_length(temp_plot_files(d), 0L)
})

test_that("a plot that cannot be put in place is a failed save", {
    d <- withr::local_tempdir()
    # Something that cannot be replaced by a file stands in for a PNG that is
    # open in another program. The draw itself succeeds.
    target <- file.path(d, "dotplot.png")
    dir.create(target)
    local_metab_enrichment_stubs(list(.draw_metab_enrichment_png = draw_ok))

    expect_error(.save_metab_enrichment_plot("a plot", target),
                 "could not be put in place")
    expect_true(dir.exists(target))
    expect_length(temp_plot_files(d), 0L)
})

test_that("a first-contrast copy whose save failed is not reported as written", {
    d <- withr::local_tempdir()
    enr <- file.path(d, "Enrichment"); ds <- file.path(d, "Datasets")
    dir.create(enr); dir.create(ds)
    earlier <- file.path(enr, "enrichment_ora_dotplot.png")
    writeLines("earlier run's plot", earlier)
    res <- list(table = data.frame(pathway = "P1", FDR = 0.01))
    plots <- list(dotplot = "a plot", lollipop = "a plot")

    # The device fails for the dotplot only, without raising an error. The file
    # it is handed is a temporary one that does not carry the plot name, so the
    # two saves are told apart by their order: the dotplot is saved first.
    calls <- 0L
    local_metab_enrichment_stubs(list(.draw_metab_enrichment_png =
        function(plot, path, width, height) {
            calls <<- calls + 1L
            if (calls > 1L) writeLines("new plot", path)
            invisible(NULL)
        }))

    files <- .save_legacy_enrichment("ora", res, plots, out_enr = enr, out_ds = ds)

    expect_identical(calls, 2L)
    expect_false(earlier %in% files)
    expect_setequal(basename(files),
                    c("enrichment_ora_results.tsv", "enrichment_ora_lollipop.png"))
    expect_true(all(file.exists(files)))
    # The earlier plot is still on disk, and is not claimed by this run.
    expect_identical(readLines(earlier), "earlier run's plot")
})

test_that("a run whose plots were not written reports none and clears the old ones", {
    fx  <- layout_fixture(contrasts = "A vs B")
    out <- file.path(withr::local_tempdir(), "metabolomics")
    enr <- file.path(out, "Enrichment")
    dir.create(enr, recursive = TRUE)
    earlier <- file.path(enr, c("enrichment_ora_dotplot.png",
                                "enrichment_gsea_nes_barplot.png"))
    file.create(earlier)
    stubs <- method_stubs(plot = function(...) "a plot")
    stubs$.draw_metab_enrichment_png <- draw_nothing
    local_metab_enrichment_stubs(stubs)

    res <- suppressMessages(suppressWarnings(mod_metabolomics_enrichment(
        fx$pre, fx$de_res, fx$config, out, inputs = fx$inputs)))

    # No plot is reported, so none can reach the PowerPoint ...
    expect_false(any(grepl("\\.png$", res$files)))
    expect_true(all(file.exists(res$files)))
    # ... and the earlier first-contrast copies were removed, not kept as current.
    expect_false(any(file.exists(earlier)))
    expect_length(temp_plot_files(out), 0L)
    # The tables are unaffected.
    expect_true(file.exists(file.path(enr, "ORA/combined/contrasts/A_vs_B/de/results.tsv")))
})

test_that("plots are written beside their table under fixed names", {
    skip_if_not_installed("ggplot2")
    skip_if_not(capabilities("png") || requireNamespace("ragg", quietly = TRUE),
                "no PNG device")
    fx  <- layout_fixture(contrasts = "A vs B")
    out <- file.path(withr::local_tempdir(), "metabolomics")
    local_metab_enrichment_stubs(method_stubs(plot = function(...) ggplot2::ggplot()))

    res <- mod_metabolomics_enrichment(fx$pre, fx$de_res, fx$config, out,
                                       inputs = fx$inputs)
    enr  <- file.path(out, "Enrichment")
    pngs <- res$files[grepl("\\.png$", res$files)]

    expect_setequal(
        sub(paste0(enr, "/"), "", pngs, fixed = TRUE),
        c("QEA/combined/all_samples/barplot.png",
          "QEA/combined/all_samples/lollipop.png",
          "ssGSEA/combined/all_samples/barplot.png",
          "ssGSEA/combined/all_samples/lollipop.png",
          "ORA/combined/contrasts/A_vs_B/de/dotplot.png",
          "ORA/combined/contrasts/A_vs_B/de/lollipop.png",
          "ORA/combined/contrasts/A_vs_B/de/barplot.png",
          "GSEA/combined/ranking_by_fc/A_vs_B/nes_barplot.png",
          "GSEA/combined/ranking_by_fc/A_vs_B/lollipop.png",
          "GSEA/combined/ranking_by_fc/A_vs_B/barplot.png",
          "QEA/combined/contrasts/A_vs_B/lollipop.png",
          "QEA/combined/contrasts/A_vs_B/barplot.png",
          "ssGSEA/combined/contrasts/A_vs_B/lollipop.png",
          "ssGSEA/combined/contrasts/A_vs_B/barplot.png",
          "enrichment_ora_dotplot.png", "enrichment_ora_lollipop.png",
          "enrichment_ora_barplot.png", "enrichment_gsea_nes_barplot.png",
          "enrichment_gsea_lollipop.png", "enrichment_gsea_barplot.png"))
    expect_true(all(file.exists(pngs)))

    # The deck labels and orders slides by the flat names these plots had before
    # the layout: one distinct name per plot, the same set as before.
    flat <- metab_enrichment_legacy_plot_name(pngs, enr)
    expect_identical(anyDuplicated(flat), 0L)
    expect_setequal(flat, c(
        "enrichment_qea_barplot.png", "enrichment_qea_lollipop.png",
        "enrichment_ssgsea_barplot.png", "enrichment_ssgsea_lollipop.png",
        "enrichment_ora_dotplot_A_vs_B.png", "enrichment_ora_lollipop_A_vs_B.png",
        "enrichment_ora_barplot_A_vs_B.png",
        "enrichment_gsea_nes_barplot_A_vs_B.png",
        "enrichment_gsea_lollipop_A_vs_B.png", "enrichment_gsea_barplot_A_vs_B.png",
        "enrichment_qea_lollipop_A_vs_B.png", "enrichment_qea_barplot_A_vs_B.png",
        "enrichment_ssgsea_lollipop_A_vs_B.png",
        "enrichment_ssgsea_barplot_A_vs_B.png",
        "enrichment_ora_dotplot.png", "enrichment_ora_lollipop.png",
        "enrichment_ora_barplot.png", "enrichment_gsea_nes_barplot.png",
        "enrichment_gsea_lollipop.png", "enrichment_gsea_barplot.png"))
})
