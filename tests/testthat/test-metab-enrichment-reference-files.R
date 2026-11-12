# tests/testthat/test-metab-enrichment-reference-files.R
#
# The metabolomics enrichment reads two reference files that sit outside the
# data folder: the GMT pathway libraries and the HMDB -> KEGG table. No target
# depended on them, so an edited library left the cached enrichment in place.
#
# They are now fingerprinted while the plan is defined, and the fingerprint is
# the command of a target (R/domain/metabolomics/06q_enrichment_reference_files.R).
# These tests cover:
#   1. which files are found, and that a missing one is left out, not an error;
#   2. the DAG: which targets depend on which fingerprint;
#   3. {targets} itself, in a throwaway project: what is rerun when a library
#      is edited, when a configured library that was missing appears, when the
#      mapping table changes, and when nothing changes;
#   4. the same for the mummichog report sections, which read both files
#      themselves and so must not rely on an upstream value changing.
#
# Everything runs on temporary files and a private store. No project input is
# read or written.

ref_config <- function(proj, gmt = NULL, mapping = NULL) {
    list(project = list(dir = proj), paths = list(raw = "data"),
         modes = list(metabolomics = list(enrichment = list(
             gmt_file = gmt, mapping_file = mapping))))
}

# A project folder with data/refs/{a.gmt, map.tsv}; b.gmt is deliberately absent.
make_ref_project <- function(env = parent.frame()) {
    proj <- withr::local_tempdir(.local_envir = env)
    refs <- file.path(proj, "data", "refs")
    dir.create(refs, recursive = TRUE)
    writeLines("P1\tfirst pathway\tC00001\tC00002", file.path(refs, "a.gmt"))
    writeLines(c("HMDB\tKEGG", "HMDB0000001\tC00001"), file.path(refs, "map.tsv"))
    list(proj = proj, refs = refs,
         a = file.path(refs, "a.gmt"), b = file.path(refs, "b.gmt"),
         map = file.path(refs, "map.tsv"))
}

norm <- function(x) normalizePath(x, winslash = "/", mustWork = TRUE)

# ---- 1. which files are found ------------------------------------------------

test_that("configured GMT files are found where the enrichment code reads them", {
    p <- make_ref_project()
    outside <- file.path(withr::local_tempdir(), "abs.gmt")
    writeLines("P9\tx\tC00009\tC00010", outside)

    cfg <- ref_config(p$proj, gmt = list("refs/a.gmt", outside))
    found <- find_metab_enrichment_gmt_files(cfg)

    # A relative path is located under the raw data folder, an absolute one kept:
    # the same resolution run_metabolomics_ora() and the others apply.
    expect_identical(found, norm(c(p$a, outside)))
    expect_identical(found, norm(resolve_input_path(cfg, cfg$modes$metabolomics$enrichment$gmt_file)))
})

test_that("a configured GMT file that does not exist is left out, not an error", {
    p <- make_ref_project()
    cfg <- ref_config(p$proj, gmt = list("refs/a.gmt", "refs/b.gmt"))

    expect_silent(found <- find_metab_enrichment_gmt_files(cfg))
    expect_identical(found, norm(p$a))

    # Once it exists it is found, after the one that was already there.
    writeLines("P2\tsecond pathway\tC00003\tC00004", p$b)
    expect_identical(find_metab_enrichment_gmt_files(cfg), norm(c(p$a, p$b)))
})

test_that("nothing configured gives no files, for every way of saying nothing", {
    p <- make_ref_project()
    empty <- list(
        ref_config(p$proj),
        ref_config(p$proj, gmt = "", mapping = ""),
        ref_config(p$proj, gmt = list(), mapping = NULL),
        list(project = list(dir = p$proj), modes = list(metabolomics = list())),
        list(project = list(dir = p$proj), modes = list())
    )
    for (cfg in empty) {
        expect_identical(find_metab_enrichment_gmt_files(cfg), character(0))
        expect_identical(find_metab_kegg_mapping_files(cfg), character(0))
        expect_identical(metab_reference_fingerprints(cfg),
                         list(gmt = character(0), mapping = character(0)))
    }
})

test_that("a file listed twice, or a directory, is not tracked twice or at all", {
    p <- make_ref_project()
    cfg <- ref_config(p$proj, gmt = list("refs/a.gmt", p$a, "refs"))
    expect_identical(find_metab_enrichment_gmt_files(cfg), norm(p$a))
})

test_that("the mapping file is found at each path a reader uses", {
    p <- make_ref_project()

    # Absolute: both readers use the same file.
    expect_identical(find_metab_kegg_mapping_files(ref_config(p$proj, mapping = p$map)),
                     norm(p$map))

    # Relative: the enrichment stage reads it under the raw data folder ...
    cfg <- ref_config(p$proj, mapping = "refs/map.tsv")
    expect_identical(find_metab_kegg_mapping_files(cfg), norm(p$map))

    # ... and mod_met_raw() reads the path as written, from the working
    # directory. Where that is a second file, both are tracked.
    wd <- withr::local_tempdir()
    dir.create(file.path(wd, "refs"))
    writeLines(c("HMDB\tKEGG", "HMDB0000002\tC00002"), file.path(wd, "refs", "map.tsv"))
    withr::local_dir(wd)
    expect_setequal(find_metab_kegg_mapping_files(cfg),
                    norm(c(p$map, file.path(wd, "refs", "map.tsv"))))
})

test_that("a fingerprint is the md5 of each file, by path", {
    p <- make_ref_project()
    fp <- fingerprint_reference_files(c(p$a, p$map))

    expect_identical(names(fp), c(p$a, p$map))
    expect_identical(unname(fp), unname(tools::md5sum(c(p$a, p$map))))
    expect_identical(fingerprint_reference_files(character(0)), character(0))

    # Content, not time: rewriting the same bytes keeps it, a new byte changes it.
    before <- fingerprint_reference_files(p$a)
    writeLines(readLines(p$a), p$a)
    expect_identical(fingerprint_reference_files(p$a), before)
    cat("P2\tsecond pathway\tC00003\tC00004\n", file = p$a, append = TRUE)
    expect_false(identical(fingerprint_reference_files(p$a), before))
})

# ---- 2. the DAG ----------------------------------------------------------------

skip_if_not_installed("targets")

ref_target_names <- function(target_list) {
    vapply(target_list, function(t) t$settings$name, character(1))
}

# $command$deps is what {targets} builds the DAG from, but it is internal;
# falling back to the symbols in the command keeps this asserting something real.
ref_target_deps <- function(tgt) {
    deps <- tgt$command$deps
    if (length(deps) > 0) return(deps)
    unlist(lapply(as.list(tgt$command$expr), all.vars))
}

ref_target_command <- function(tgt) {
    expr <- tgt$command$expr
    if (is.list(expr) || is.expression(expr)) expr <- expr[[1]]
    expr
}

test_that("a reference target's command is its fingerprint, and so is its value", {
    fp <- c("/proj/data/refs/a.gmt" = "0123abcd", "/proj/data/refs/b.gmt" = "4567ef00")
    tgt <- metab_reference_target("metab_enrichment_gmt_ref", fp)

    expect_identical(tgt$settings$name, "metab_enrichment_gmt_ref")
    expect_identical(eval(ref_target_command(tgt)), fp)
    # The md5 is in the command text: that is what makes a changed file a
    # changed command.
    expect_true(grepl("4567ef00", paste(deparse(ref_target_command(tgt)), collapse = ""),
                      fixed = TRUE))

    none <- metab_reference_target("metab_enrichment_gmt_ref", character(0))
    expect_identical(eval(ref_target_command(none)), character(0))
    expect_identical(eval(ref_target_command(
        metab_reference_target("metab_enrichment_gmt_ref", NULL))), character(0))

    # A different fingerprint is a different command; the same one is the same.
    changed <- fp; changed[[1]] <- "ffffffff"
    expect_false(identical(ref_target_command(tgt),
                           ref_target_command(metab_reference_target("x", changed))))
    expect_identical(ref_target_command(tgt),
                     ref_target_command(metab_reference_target("x", fp)))

    expect_error(metab_reference_target("x", c("0123abcd")), "named character")
    expect_error(metab_reference_target("x", list(a = "b")), "named character")
})

test_that("met_raw depends on the mapping file only, enrichment on both", {
    fp <- list(gmt     = c("/proj/data/refs/a.gmt" = "0123abcd"),
               mapping = c("/proj/data/refs/map.tsv" = "89abcdef"))
    result  <- pipe_metabolomics("pqn", reference_files = fp)
    by_name <- setNames(result, ref_target_names(result))

    expect_identical(eval(ref_target_command(by_name$metab_enrichment_gmt_ref)), fp$gmt)
    expect_identical(eval(ref_target_command(by_name$metab_kegg_mapping_ref)), fp$mapping)

    raw_deps <- ref_target_deps(by_name$met_raw)
    expect_true("metab_kegg_mapping_ref" %in% raw_deps)
    # A changed library must not rebuild the raw matrix and everything after it.
    expect_false("metab_enrichment_gmt_ref" %in% raw_deps)

    enr_deps <- ref_target_deps(by_name$metab_enrichment_res)
    expect_true(all(c("metab_kegg_mapping_ref", "metab_enrichment_gmt_ref") %in% enr_deps))

    # Both targets still call what they called before.
    expect_true("mod_met_raw" %in% all.names(ref_target_command(by_name$met_raw)))
    expect_true("mod_metabolomics_enrichment" %in%
                    all.names(ref_target_command(by_name$metab_enrichment_res)))

    # No other target reads these files, so no other target depends on them.
    others <- setdiff(names(by_name), c("met_raw", "metab_enrichment_res"))
    for (nm in others) {
        deps <- ref_target_deps(by_name[[nm]])
        expect_false(any(c("metab_kegg_mapping_ref", "metab_enrichment_gmt_ref") %in% deps),
                     info = sprintf("%s depends on a reference fingerprint", nm))
    }
})

test_that("the mummichog sections depend on both fingerprints, the mummichog run on neither", {
    fp <- list(gmt     = c("/proj/data/refs/a.gmt" = "0123abcd"),
               mapping = c("/proj/data/refs/map.tsv" = "89abcdef"))
    refs <- c("metab_kegg_mapping_ref", "metab_enrichment_gmt_ref")

    for (engine in c("metaboanalyst", "python")) {
        result  <- pipe_metabolomics("pqn", mummichog_enabled = TRUE,
                                     mummichog_engine = engine, reference_files = fp)
        by_name <- setNames(result, ref_target_names(result))

        sec <- by_name$metab_mummichog_sections
        expect_true(all(refs %in% ref_target_deps(sec)), info = engine)
        # It still calls what it called before, on the same inputs.
        expect_true("mod_mummichog_sections" %in% all.names(ref_target_command(sec)),
                    info = engine)
        expect_true(all(c("metab_mummichog_report_pathways", "metab_mummichog_report_files",
                          "metab_pre", "metab_de_res", "metab_enrichment_res") %in%
                            ref_target_deps(sec)), info = engine)

        # The run and the targets that only pass its files on read neither
        # file, so a changed library must not reach them directly.
        for (nm in c("metab_mummichog_pinned_files", "metab_mummichog_report_pathways",
                     "metab_mummichog_report_files", "metab_mummichog_exports")) {
            expect_false(any(refs %in% ref_target_deps(by_name[[nm]])),
                         info = sprintf("%s (%s) depends on a reference fingerprint",
                                        nm, engine))
        }

        # With mummichog enabled the readers are exactly these three.
        readers <- names(by_name)[vapply(by_name, function(t) {
            any(refs %in% ref_target_deps(t))
        }, logical(1))]
        expect_setequal(readers, c("met_raw", "metab_enrichment_res",
                                   "metab_mummichog_sections"))
    }

    # Under multi-omics the sections are still built, and still tracked.
    core <- pipe_metabolomics("pqn", skip_outputs = TRUE, mummichog_enabled = TRUE,
                              reference_files = fp)
    core <- setNames(core, ref_target_names(core))
    expect_true(all(refs %in% names(core)))
    expect_true(all(refs %in% ref_target_deps(core$metab_mummichog_sections)))

    # Disabled: the sections are a constant NULL and depend on nothing.
    off <- pipe_metabolomics("pqn", reference_files = fp)
    off <- setNames(off, ref_target_names(off))
    expect_false(any(refs %in% ref_target_deps(off$metab_mummichog_sections)))
})

test_that("the reference targets exist in every mode that has their readers", {
    # No argument: the two targets are there with nothing to track, so callers
    # that do not pass fingerprints build the same plan as before plus two
    # empty targets.
    plain <- pipe_metabolomics("pqn")
    names_plain <- ref_target_names(plain)
    expect_true(all(c("metab_kegg_mapping_ref", "metab_enrichment_gmt_ref") %in% names_plain))
    by_name <- setNames(plain, names_plain)
    expect_identical(eval(ref_target_command(by_name$metab_kegg_mapping_ref)), character(0))
    expect_identical(eval(ref_target_command(by_name$metab_enrichment_gmt_ref)), character(0))

    # Review mode has met_raw but no enrichment stage.
    review <- ref_target_names(pipe_metabolomics(NULL))
    expect_true("metab_kegg_mapping_ref" %in% review)
    expect_false("metab_enrichment_res" %in% review)
    expect_false("metab_enrichment_gmt_ref" %in% review)

    # Under multi-omics the enrichment stage still runs, and is still tracked.
    core <- ref_target_names(pipe_metabolomics("pqn", skip_outputs = TRUE))
    expect_true(all(c("metab_kegg_mapping_ref", "metab_enrichment_gmt_ref",
                      "metab_enrichment_res") %in% core))
})

# ---- 3. {targets}, in a throwaway project ---------------------------------------
#
# The readers are stand-ins that only record that they ran; the fingerprints and
# the reference targets are the real ones. Section 2 shows the real readers
# depend on the same two targets in the same way.

test_that("{targets} reruns the readers of a changed or newly present reference file, and nothing else", {
    p <- make_ref_project()
    log_file <- file.path(withr::local_tempdir(), "ran.log")
    r_file <- function(...) norm(file.path(root_dir, "R", ...))

    script <- c(
        "library(targets)",
        sprintf("source(%s)", deparse(r_file("core", "00_paths.R"))),
        sprintf("source(%s)", deparse(r_file("domain", "metabolomics",
                                             "06q_enrichment_reference_files.R"))),
        sprintf("source(%s)", deparse(r_file("pipeline", "metabolomics",
                                             "00_pipe_metabolomics.R"))),
        sprintf("ran <- function(name) cat(name, '\\n', sep = '', file = %s, append = TRUE)",
                deparse(log_file)),
        sprintf(paste0("cfg <- list(project = list(dir = %s), paths = list(raw = 'data'), ",
                       "modes = list(metabolomics = list(enrichment = list(",
                       "gmt_file = list('refs/a.gmt', 'refs/b.gmt'), ",
                       "mapping_file = 'refs/map.tsv'))))"), deparse(p$proj)),
        "refs <- metab_reference_fingerprints(cfg)",
        "list(",
        "  metab_reference_target('metab_kegg_mapping_ref', refs$mapping),",
        "  metab_reference_target('metab_enrichment_gmt_ref', refs$gmt),",
        "  tar_target(toy_met_raw, { metab_kegg_mapping_ref; ran('toy_met_raw'); 1 }),",
        "  tar_target(toy_enrichment, { metab_kegg_mapping_ref; metab_enrichment_gmt_ref;",
        "                               ran('toy_enrichment'); toy_met_raw + 1 }),",
        "  tar_target(toy_report, { ran('toy_report'); toy_enrichment + 1 })",
        ")"
    )

    targets::tar_dir({
        writeLines(script, "_targets.R")

        outdated <- function() {
            targets::tar_outdated(callr_function = NULL, reporter = "silent",
                                  envir = new.env(parent = globalenv()))
        }
        # Run the plan and return the stand-in readers that actually ran.
        make <- function() {
            unlink(log_file)
            suppressMessages(targets::tar_make(callr_function = NULL, reporter = "silent",
                                               envir = new.env(parent = globalenv())))
            if (file.exists(log_file)) sort(readLines(log_file)) else character(0)
        }
        refs_only <- function(x) intersect(
            x, c("metab_kegg_mapping_ref", "metab_enrichment_gmt_ref"))

        # First run: everything is built once.
        expect_identical(make(), c("toy_enrichment", "toy_met_raw", "toy_report"))

        # Nothing changed: nothing is outdated and nothing runs.
        expect_identical(outdated(), character(0))
        expect_identical(make(), character(0))

        # An existing GMT is edited in place.
        cat("P2\tsecond pathway\tC00003\tC00004\n", file = p$a, append = TRUE)
        out <- outdated()
        expect_identical(refs_only(out), "metab_enrichment_gmt_ref")
        expect_true("toy_enrichment" %in% out)
        expect_false("toy_met_raw" %in% out)
        # Enrichment reruns; the raw stage does not. The stand-in returns the
        # same value, so what follows it does not rerun either.
        expect_identical(make(), "toy_enrichment")
        expect_identical(outdated(), character(0))

        # A configured GMT that was missing appears.
        writeLines("P3\tthird pathway\tC00005\tC00006", p$b)
        out <- outdated()
        expect_identical(refs_only(out), "metab_enrichment_gmt_ref")
        expect_true("toy_enrichment" %in% out)
        expect_false("toy_met_raw" %in% out)
        expect_identical(make(), "toy_enrichment")
        expect_identical(outdated(), character(0))

        # The mapping table changes: both of its readers rerun.
        cat("HMDB0000002\tC00002\n", file = p$map, append = TRUE)
        out <- outdated()
        expect_identical(refs_only(out), "metab_kegg_mapping_ref")
        expect_true(all(c("toy_met_raw", "toy_enrichment") %in% out))
        expect_identical(make(), c("toy_enrichment", "toy_met_raw"))
        expect_identical(outdated(), character(0))

        # A tracked GMT is removed: enrichment reruns without it.
        unlink(p$b)
        out <- outdated()
        expect_identical(refs_only(out), "metab_enrichment_gmt_ref")
        expect_false("toy_met_raw" %in% out)
        expect_identical(make(), "toy_enrichment")

        # And again nothing changed: nothing is outdated and nothing runs.
        expect_identical(outdated(), character(0))
        expect_identical(make(), character(0))
    })
})

# ---- 4. the mummichog report sections ---------------------------------------------
#
# The sections read the GMT libraries and the mapping table themselves. Before
# they depended on the fingerprints, they reran only if the enrichment result
# or the preprocessed data came back different -- and neither has to. Here the
# stand-ins are built so that neither does:
#   toy_enrichment  returns NULL, as the enrichment stage does when
#                   `run_enrichment` is off;
#   toy_pre         follows toy_met_raw, whose value does not change, so it is
#                   never rebuilt after the first run.
# toy_sections_indirect is the sections target as it was: the same inputs
# without the two fingerprints. It is there to show the gap the dependency
# closes, and that these scenarios would have caught it.

test_that("{targets} reruns the mummichog sections on a changed reference file when nothing upstream changed", {
    p <- make_ref_project()
    log_file <- file.path(withr::local_tempdir(), "ran.log")
    r_file <- function(...) norm(file.path(root_dir, "R", ...))

    script <- c(
        "library(targets)",
        sprintf("source(%s)", deparse(r_file("core", "00_paths.R"))),
        sprintf("source(%s)", deparse(r_file("domain", "metabolomics",
                                             "06q_enrichment_reference_files.R"))),
        sprintf("source(%s)", deparse(r_file("pipeline", "metabolomics",
                                             "00_pipe_metabolomics.R"))),
        sprintf("ran <- function(name) cat(name, '\\n', sep = '', file = %s, append = TRUE)",
                deparse(log_file)),
        sprintf(paste0("cfg <- list(project = list(dir = %s), paths = list(raw = 'data'), ",
                       "modes = list(metabolomics = list(enrichment = list(",
                       "run_enrichment = FALSE, gmt_file = 'refs/a.gmt', ",
                       "mapping_file = 'refs/map.tsv'))))"), deparse(p$proj)),
        "refs <- metab_reference_fingerprints(cfg)",
        "list(",
        "  metab_reference_target('metab_kegg_mapping_ref', refs$mapping),",
        "  metab_reference_target('metab_enrichment_gmt_ref', refs$gmt),",
        "  tar_target(toy_met_raw, { metab_kegg_mapping_ref; ran('toy_met_raw'); 1 }),",
        "  tar_target(toy_pre, { ran('toy_pre'); toy_met_raw + 1 }),",
        "  tar_target(toy_enrichment, { metab_kegg_mapping_ref; metab_enrichment_gmt_ref;",
        "                               toy_pre; ran('toy_enrichment'); NULL }),",
        "  tar_target(toy_mummichog_run, { ran('toy_mummichog_run'); toy_pre + 1 }),",
        "  tar_target(toy_sections, { metab_kegg_mapping_ref; metab_enrichment_gmt_ref;",
        "                             ran('toy_sections');",
        "                             list(toy_mummichog_run, toy_pre, toy_enrichment) }),",
        "  tar_target(toy_sections_indirect, { ran('toy_sections_indirect');",
        "                             list(toy_mummichog_run, toy_pre, toy_enrichment) }),",
        "  tar_target(toy_exports, { ran('toy_exports'); length(toy_sections) })",
        ")"
    )

    targets::tar_dir({
        writeLines(script, "_targets.R")

        outdated <- function() {
            targets::tar_outdated(callr_function = NULL, reporter = "silent",
                                  envir = new.env(parent = globalenv()))
        }
        make <- function() {
            unlink(log_file)
            suppressMessages(targets::tar_make(callr_function = NULL, reporter = "silent",
                                               envir = new.env(parent = globalenv())))
            if (file.exists(log_file)) sort(readLines(log_file)) else character(0)
        }

        # First run: everything is built once.
        expect_identical(make(), sort(c("toy_met_raw", "toy_pre", "toy_enrichment",
                                        "toy_mummichog_run", "toy_sections",
                                        "toy_sections_indirect", "toy_exports")))
        expect_null(targets::tar_read(toy_enrichment))

        # No reference change: nothing is outdated, the sections do not rerun.
        expect_identical(outdated(), character(0))
        expect_identical(make(), character(0))

        # A GMT changes while identified enrichment is off.
        cat("P2\tsecond pathway\tC00003\tC00004\n", file = p$a, append = TRUE)
        out <- outdated()
        expect_true("toy_sections" %in% out)
        # Upstream isolation: not the raw stage, not the data, not the run.
        expect_false(any(c("toy_met_raw", "toy_pre", "toy_mummichog_run") %in% out))
        ran_now <- make()
        expect_identical(ran_now, c("toy_enrichment", "toy_sections"))
        # Enrichment reran and is still NULL: its value gave the sections no
        # reason to rerun, the fingerprint did. Without it they stay stale.
        expect_null(targets::tar_read(toy_enrichment))
        expect_false("toy_sections_indirect" %in% ran_now)
        expect_identical(outdated(), character(0))

        # The mapping table changes. The raw stage reruns and returns what it
        # returned before, so the preprocessed data is not rebuilt.
        cat("HMDB0000002\tC00002\n", file = p$map, append = TRUE)
        out <- outdated()
        expect_true(all(c("toy_met_raw", "toy_sections") %in% out))
        ran_now <- make()
        expect_identical(ran_now, c("toy_enrichment", "toy_met_raw", "toy_sections"))
        expect_false(any(c("toy_pre", "toy_mummichog_run", "toy_sections_indirect") %in%
                             ran_now))
        expect_identical(outdated(), character(0))

        # And again no reference change: nothing is outdated and nothing runs.
        expect_identical(outdated(), character(0))
        expect_identical(make(), character(0))
    })
})
