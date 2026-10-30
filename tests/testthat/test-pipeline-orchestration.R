# tests/testthat/test-pipeline-orchestration.R
# Tests for pipeline orchestration functions (pipe_rnaseq, pipe_proteomics, etc.)
# Validates that each pipe function returns correctly structured lists of tar_target objects
# without evaluating the underlying data (testing DAG creation only).

skip_if_not_installed("targets")
library(targets)

# =============================================================================
# Helper: Extract target names from a list of tar_target objects
# =============================================================================
get_target_names <- function(target_list) {
    vapply(target_list, function(t) t$settings$name, character(1))
}

# =============================================================================
# Helper: Validate that every element in a target list is a proper target object
# =============================================================================
assert_valid_targets <- function(target_list) {
    expect_type(target_list, "list")
    expect_true(length(target_list) > 0, info = "Target list should not be empty")

    for (i in seq_along(target_list)) {
        tgt <- target_list[[i]]
        expect_true(
            !is.null(tgt$settings$name),
            info = sprintf("Element %d should have $settings$name", i)
        )
    }
}


# =============================================================================
# Tests: pipe_rnaseq()
# =============================================================================

test_that("pipe_rnaseq returns a valid target list", {
    result <- pipe_rnaseq()
    assert_valid_targets(result)
})

test_that("pipe_rnaseq contains expected target names", {
    result <- pipe_rnaseq()
    names <- get_target_names(result)
    expect_true("rna_input_files" %in% names)
    expect_true("rna_inputs" %in% names)
    expect_true("rna_pre" %in% names)
    expect_true("rna_de_res" %in% names)
    expect_true("rna_out_dir" %in% names)
})

test_that("pipe_rnaseq returns expected number of targets", {
    result <- pipe_rnaseq()
    expect_gte(length(result), 15)
})

test_that("everything that reads an RNA input file depends on rna_input_files", {
    # rna_input_files is format = "file", so editing a tracked file invalidates
    # whatever depends on it. rna_annot and rna_trinotate_main read
    # config$modes$rna$files themselves rather than going through rna_inputs --
    # so wiring only rna_inputs would leave an edited annotation file rebuilding
    # the inputs while these two kept serving objects built from the old file.
    result <- pipe_rnaseq()
    by_name <- setNames(result, get_target_names(result))

    expect_identical(by_name$rna_input_files$settings$format, "file")

    # $command$deps is what targets builds the DAG from, but it is internal;
    # falling back to the symbols in the command expression keeps this asserting
    # something real rather than passing vacuously if that field ever moves.
    target_deps <- function(tgt) {
        deps <- tgt$command$deps
        if (length(deps) > 0) return(deps)
        unlist(lapply(as.list(tgt$command$expr), all.vars))
    }

    for (nm in c("rna_inputs", "rna_annot", "rna_trinotate_main")) {
        deps <- target_deps(by_name[[nm]])
        expect_true(length(deps) > 0,
                    info = sprintf("%s reports no dependencies at all", nm))
        expect_true("rna_input_files" %in% deps,
                    info = sprintf("%s does not depend on rna_input_files", nm))
    }
})

test_that("every configured RNA input file is tracked, de_table lists included", {
    # load_omics_inputs() skips multi-value keys, but load_precomputed_rna_de()
    # reads files$de_table directly and supports one table per contrast. Left
    # out of the file target, editing one of them would leave the RNA branch
    # looking up to date -- exactly the staleness this target exists to stop.
    #
    # Evaluated rather than pattern-matched on the source: what matters is the
    # set of paths the target actually yields.
    cfg <- list(
        project = list(dir = "/proj"),
        paths = list(raw = "data"),
        modes = list(rna = list(files = list(
            metadata          = "rna/meta.csv",
            contrasts         = "rna/contrasts.csv",
            de_table          = list("rna/de_A.tsv", "rna/de_B.tsv"),
            annotation        = "",           # configured but empty: skipped
            is_logtransformed = FALSE         # a flag, not a path
        )))
    )

    result <- pipe_rnaseq()
    tgt <- setNames(result, get_target_names(result))$rna_input_files
    expr <- tgt$command$expr
    if (is.list(expr) || is.expression(expr)) expr <- expr[[1]]

    env <- new.env(parent = globalenv())
    env$config <- cfg
    paths <- eval(expr, envir = env)

    expect_true(all(c("/proj/data/rna/de_A.tsv", "/proj/data/rna/de_B.tsv")
                    %in% paths))
    expect_true(all(c("/proj/data/rna/meta.csv", "/proj/data/rna/contrasts.csv")
                    %in% paths))
    # An empty path and a logical flag are not files; a duplicate would make
    # targets track the same file twice.
    expect_length(paths, 4L)
    expect_false(any(duplicated(paths)))
})


# =============================================================================
# Tests: pipe_proteomics()
# =============================================================================

test_that("pipe_proteomics returns a valid target list", {
    result <- pipe_proteomics()
    assert_valid_targets(result)
})

test_that("pipe_proteomics contains expected target names", {
    result <- pipe_proteomics()
    names <- get_target_names(result)
    expect_true("prot_out_dir" %in% names)
    expect_true("prot_pre" %in% names)
    expect_true("prot_de_res" %in% names)
})

test_that("pipe_proteomics adds the pooled-ranking target beside the legacy one", {
    result <- pipe_proteomics()
    by_name <- setNames(result, get_target_names(result))
    expect_true("prot_pathway_pooled_res" %in% names(by_name))
    expect_true("prot_pathway_res" %in% names(by_name))

    target_expr <- function(tgt) {
        expr <- tgt$command$expr
        if (is.list(expr) || is.expression(expr)) expr <- expr[[1]]
        expr
    }
    # The legacy target still calls the legacy module, unchanged.
    expect_identical(target_expr(by_name$prot_pathway_res)[[1]], quote(mod_proteomics_pathway))
    expect_identical(target_expr(by_name$prot_pathway_pooled_res)[[1]],
                     quote(mod_proteomics_pathway_pooled))

    # The report waits for it.
    report_vars <- all.vars(target_expr(by_name$prot_report))
    expect_true("prot_pathway_pooled_res" %in% report_vars)
})

test_that("the pooled-ranking target is a single-omics output only", {
    # multiomics consumes prot_pathway_res; the pooled target must not be built
    # when only the core proteomics targets are requested.
    names_core <- get_target_names(pipe_proteomics(skip_outputs = TRUE))
    expect_false("prot_pathway_pooled_res" %in% names_core)
    expect_true("prot_pathway_res" %in% names_core)
})

test_that("pipe_proteomics returns expected number of targets", {
    result <- pipe_proteomics()
    expect_gte(length(result), 10)
})


# =============================================================================
# Tests: pipe_metabolomics()
# =============================================================================

test_that("pipe_metabolomics returns a valid target list", {
    result <- pipe_metabolomics("pqn")
    assert_valid_targets(result)
})

test_that("pipe_metabolomics contains expected target names", {
    result <- pipe_metabolomics("pqn")
    names <- get_target_names(result)
    expect_true("metab_out_dir" %in% names)
    expect_true("metab_pre" %in% names)
    expect_true("metab_de_res" %in% names)
})

test_that("pipe_metabolomics returns expected number of targets", {
    result <- pipe_metabolomics("pqn")
    expect_gte(length(result), 8)
})

test_that("the metab shiny payload target passes only arguments the builder takes", {
    # save_shiny_payload_metabolomics() takes `...` and forwards it straight
    # into build_shiny_payload_metabolomics(), so an argument name the builder
    # has no formal for is an "unused argument" error R raises only while the
    # target builds. The unit tests call the builder directly and never see it,
    # and so does lint -- which leaves this call site uncovered. Check it here.
    result <- pipe_metabolomics("pqn")
    tgt <- setNames(result, get_target_names(result))$metab_shiny_payload
    expect_false(is.null(tgt))

    expr <- tgt$command$expr
    if (is.list(expr) || is.expression(expr)) expr <- expr[[1]]
    expect_identical(expr[[1]], quote(save_shiny_payload_metabolomics))

    # out_file is the wrapper's own formal; out_dir it derives and supplies.
    passed  <- setdiff(names(as.list(expr)), c("", "out_file"))
    expect_true(length(passed) > 0)
    expect_equal(setdiff(passed, names(formals(build_shiny_payload_metabolomics))),
                 character(0))
})

test_that("the report's mummichog file list tracks content, not just paths", {
    # metab_mummichog_report_files only re-exposes the pinned stage's tracked
    # paths. As an ordinary target it would hash that character vector, so a
    # rerun rewriting an EmpiricalCompound table in place would leave its hash
    # unchanged and metab_report up to date on last run's evidence. format =
    # "file" makes {targets} hash the files themselves, as it does upstream.
    # mummichog is opt-in, and the alias only exists in the enabled branch --
    # pipe_metabolomics() defaults to FALSE, which builds the NULL stand-ins.
    on <- pipe_metabolomics("pqn", mummichog_enabled = TRUE)
    on <- setNames(on, get_target_names(on))

    expect_identical(on$metab_mummichog_pinned_files$settings$format, "file")
    expect_identical(on$metab_mummichog_report_files$settings$format, "file")

    # Disabled branch unchanged: the same two symbols the report depends on,
    # as ordinary targets, and no pinned stage.
    off <- pipe_metabolomics("pqn")
    off <- setNames(off, get_target_names(off))

    expect_true(all(c("metab_mummichog_report_pathways",
                      "metab_mummichog_report_files") %in% names(off)))
    expect_false("metab_mummichog_pinned_files" %in% names(off))
    expect_false(identical(off$metab_mummichog_report_files$settings$format, "file"))
})


# -----------------------------------------------------------------------------
# Issue #306: the mummichog exports must not hang off HTML rendering.
#
# These assert the DAG, not the filenames: {targets}' rebuild behaviour follows
# from which targets exist, in which block, with which format and which edges.
# Checking that here costs nothing and needs no pipeline run.
# -----------------------------------------------------------------------------

# $command$deps is what targets builds the DAG from, but it is internal; falling
# back to the symbols in the command expression keeps these asserting something
# real rather than passing vacuously if that field ever moves.
mummi_target_deps <- function(tgt) {
    deps <- tgt$command$deps
    if (length(deps) > 0) return(deps)
    unlist(lapply(as.list(tgt$command$expr), all.vars))
}

test_that("a deleted mummichog export rebuilds only the export target", {
    on <- pipe_metabolomics("pqn", mummichog_enabled = TRUE)
    on <- setNames(on, get_target_names(on))

    # format = "file" is the whole self-healing mechanism: targets hashes the
    # written files, so deleting one makes THIS target outdated.
    expect_identical(on$metab_mummichog_exports$settings$format, "file")

    # ... and nothing else, because the edge runs sections -> exports and never
    # back. If the sections target depended on the exports, a deleted PNG would
    # re-run the GSEA permutations that produced it.
    sec_deps <- mummi_target_deps(on$metab_mummichog_sections)
    expect_true(length(sec_deps) > 0)
    expect_false("metab_mummichog_exports" %in% sec_deps)

    expect_true("metab_mummichog_sections" %in%
                    mummi_target_deps(on$metab_mummichog_exports))
})

test_that("the HTML report does not depend on the mummichog exports", {
    on <- pipe_metabolomics("pqn", mummichog_enabled = TRUE)
    on <- setNames(on, get_target_names(on))

    report_deps <- mummi_target_deps(on$metab_report)
    expect_true(length(report_deps) > 0)

    # The report consumes the SECTIONS, so it renders what the exports wrote
    # without being invalidated when one of those files is touched.
    expect_true("metab_mummichog_sections" %in% report_deps)
    expect_false("metab_mummichog_exports" %in% report_deps)

    # Nothing at all consumes the exports -- it is a leaf. A consumer would drag
    # the export writing back into some other target's rebuild.
    consumers <- names(Filter(function(tgt)
        "metab_mummichog_exports" %in% mummi_target_deps(tgt), on))
    expect_equal(consumers, character(0))
})

test_that("multiomics mode still builds the mummichog sections and exports", {
    # skip_outputs drops the whole outputs block -- report, Shiny payload, Excel
    # -- but the mummichog stage runs in core, and #306 is that its analytical
    # outputs used to be dropped with the report.
    multi <- get_target_names(
        pipe_metabolomics("pqn", skip_outputs = TRUE, mummichog_enabled = TRUE))

    expect_true("metab_mummichog_pinned_files" %in% multi)
    expect_true("metab_mummichog_sections"     %in% multi)
    expect_true("metab_mummichog_exports"      %in% multi)

    # The outputs block really is gone, so this is not passing vacuously.
    expect_false("metab_report"         %in% multi)
    expect_false("metab_shiny_payload"  %in% multi)

    # Every input the two targets need survives the truncation too.
    expect_true(all(c("metab_out_dir", "metab_pre", "metab_de_res",
                      "metab_mummichog_report_pathways",
                      "metab_mummichog_report_files") %in% multi))
})

test_that("the mummichog exports are not reachable through the report module", {
    on <- pipe_metabolomics("pqn", mummichog_enabled = TRUE)
    on <- setNames(on, get_target_names(on))

    # The exports target calls the export wrapper directly. If it went through
    # mod_metabolomics_report(), a missing `rmarkdown` would return early and
    # silently take the evidence and GSEA exports with it -- which is the second
    # half of #306.
    expr <- on$metab_mummichog_exports$command$expr
    if (is.list(expr) || is.expression(expr)) expr <- expr[[1]]
    expect_identical(expr[[1]], quote(mod_mummichog_exports))

    # And the report no longer takes the raw mummichog inputs, so it cannot
    # rebuild sections or write exports of its own.
    report_formals <- names(formals(mod_metabolomics_report))
    expect_true("mummichog_sections" %in% report_formals)
    expect_false("mummichog_files"    %in% report_formals)
    expect_false("mummichog_pathways" %in% report_formals)
})

test_that("a disabled mummichog config adds no sections or exports targets", {
    off <- get_target_names(pipe_metabolomics("pqn"))

    # The report still needs a stable symbol to depend on in both branches ...
    expect_true("metab_mummichog_sections" %in% off)
    # ... but a disabled stage has nothing to write.
    expect_false("metab_mummichog_exports" %in% off)
})


# =============================================================================
# Tests: pipe_multiomics() — removed
# =============================================================================
# Multi-omics integration layer was split into a separate review branch.
# These tests will be restored when pipe_multiomics() is merged back.


# =============================================================================
# Tests: Target name uniqueness across pipes
# =============================================================================

test_that("target names are unique within each pipe", {
    pipes <- list(
        rna = pipe_rnaseq(),
        prot = pipe_proteomics(),
        metab = pipe_metabolomics()
    )

    for (pipe_name in names(pipes)) {
        target_names <- get_target_names(pipes[[pipe_name]])
        expect_equal(
            length(target_names), length(unique(target_names)),
            info = sprintf("Duplicate target names in pipe_%s", pipe_name)
        )
    }
})
