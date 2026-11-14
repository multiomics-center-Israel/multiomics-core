# tests/testthat/test-metab-shiny-export-enrichment.R
#
# The enrichment block of the metabolomics Shiny payload. It is added whenever
# the run has enrichment results; there is no option to set.
#
#   1. with results: the nine-component `enrichment` block and
#      `enrichment_provenance` are added, and no other key changes;
#   2. without results: nothing is added, or the block says there are none;
#   3. results that cannot be encoded: the payload is still written, with a
#      warning and the reasons on record;
#   4. mummichog is built on its own, with or without enrichment results;
#   5. the saved RDS carries the block and can be read back;
#   6. the {targets} wiring of the payload target.
#
# Everything is synthetic and lives in a temporary folder.

# ---- the synthetic run -----------------------------------------------------------
#
# Features g01..g30 are compounds C00001..C00030; m1a, m1b, m1c are three
# features of C00101. One contrast, "A_vs_B", in which g01..g08, m1a and m1c are
# significant, so that PW1 (C00001..10 + C00101) is over-represented.
make_export_enrichment_fixture <- function(env = parent.frame()) {
    dir <- withr::local_tempdir(.local_envir = env)
    samples <- paste0("s", 1:6)
    feats   <- c(sprintf("g%02d", 1:30), "m1a", "m1b", "m1c")
    kegg    <- c(sprintf("C%05d", 1:30), "C00101", "C00101", "C00101")

    mat <- withr::with_seed(1, matrix(
        stats::rnorm(length(feats) * 6, mean = 10), nrow = length(feats),
        dimnames = list(feats, samples)))
    meta <- data.frame(sample_id = samples, grp = rep(c("A", "B"), 3),
                       row.names = samples, stringsAsFactors = FALSE)
    row_data <- data.frame(feature_id = feats, original_id = paste0("orig_", feats),
                           name = paste0("met_", feats), KEGG = kegg,
                           row.names = feats, stringsAsFactors = FALSE)
    pre <- list(expr_raw = mat, expr_filt = mat, expr_work = mat,
                meta = meta, row_data = row_data)

    i <- 1:30
    logfc <- c(ifelse(i <= 10, 2 + 0.01 * i, (-1)^i * (0.2 + 0.01 * i)), 2.5, -3.0, 2.6)
    significant <- feats %in% c(sprintf("g%02d", 1:8), "m1a", "m1c")
    pval <- ifelse(significant, 0.001, 0.5)
    de_table <- data.frame(feature_id = feats, logFC = logfc, P.Value = pval,
                           adj.P.Val = pval, stringsAsFactors = FALSE)
    summary_df <- data.frame(
        feature_id        = feats,
        linearFC.A_vs_B   = sign(logfc) * 2^abs(logfc),
        pvalue.A_vs_B     = pval,
        padj.A_vs_B       = pval,
        pass.A_vs_B       = as.integer(significant),
        pass_any_contrast = as.integer(significant),
        stringsAsFactors  = FALSE)
    de_res <- list(summary_df = summary_df, de_tables = list(A_vs_B = de_table))
    inputs <- list(contrasts = data.frame(Contrast_name = "A_vs_B",
                                          stringsAsFactors = FALSE))

    gmt_line <- function(id, desc, members) paste(c(id, desc, members), collapse = "\t")
    lib <- file.path(dir, "libA.gmt")
    writeLines(c(gmt_line("PW1", "First",  c(kegg[1:10], "C00101")),
                 gmt_line("PW2", "Second", kegg[11:20]),
                 gmt_line("PW3", "Third",  kegg[21:30])), lib)

    config <- list(
        project = list(dir = dir), paths = list(raw = "data"),
        modes = list(metabolomics = list(
            de            = list(padj_cutoff = 0.05, p_cutoff = 0.05,
                                 linear_fc_cutoff = 1.5),
            preprocessing = list(chosen_norm = "pqn", transform = "log2", scaling = "none"),
            effects       = list(samples = "sample_id", color = "grp"),
            excel         = list(group_cv = FALSE),
            enrichment    = list(run_enrichment = TRUE, gmt_file = lib))))

    quiet <- function(expr) suppressWarnings(suppressMessages(expr))
    ora  <- quiet(run_metabolomics_ora(pre, de_res, config, contrast = "A_vs_B"))
    gsea <- quiet(withr::with_seed(11, run_metabolomics_gsea(pre, de_res, config,
                                                            contrast = "A_vs_B")))
    enrichment_res <- list(
        ora = ora, gsea = gsea,
        ora_by_contrast  = list(A_vs_B = list(res = ora)),
        gsea_by_contrast = list(A_vs_B = list(res = gsea)))

    list(pre = pre, de_res = de_res, inputs = inputs, config = config,
         enrichment_res = enrichment_res, feats = feats)
}

build_export <- function(fx, enrichment_res = fx$enrichment_res, include_legacy = TRUE,
                         ...) {
    suppressMessages(suppressWarnings(
        build_shiny_payload_metabolomics(
            pre = fx$pre, de_res = fx$de_res, inputs = fx$inputs, config = fx$config,
            enrichment_res = enrichment_res, include_legacy = include_legacy, ...)))
}

# The payload as it was before the enrichment block existed: the same build,
# with the block's keys left out.
build_export_without_block <- function(fx, ...) {
    local_export_stub("metab_shiny_enrichment_keys", function(...) list())
    build_export(fx, ...)
}

# Build the payload and return it together with every warning it raised.
build_export_with_warnings <- function(fx, enrichment_res) {
    warnings <- character(0)
    payload <- withCallingHandlers(
        suppressMessages(build_shiny_payload_metabolomics(
            pre = fx$pre, de_res = fx$de_res, inputs = fx$inputs, config = fx$config,
            enrichment_res = enrichment_res, include_legacy = TRUE)),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    list(payload = payload, warnings = warnings)
}

# Replace a function where the pipeline functions live, for one test.
local_export_stub <- function(name, fn, env = parent.frame()) {
    target <- environment(build_shiny_payload_metabolomics)
    old <- get(name, envir = target, inherits = FALSE)
    withr::defer(assign(name, old, envir = target), envir = env)
    assign(name, fn, envir = target)
    invisible(NULL)
}

new_keys <- c("enrichment", "enrichment_provenance")

# Two builds differ in the time they were made, and in nothing else.
comparable <- function(payload) {
    payload[setdiff(names(payload), "payload_created_at")]
}

# ---- 1. with results ---------------------------------------------------------------

test_that("the block is added without any option, and nothing else changes", {
    skip_if_not_installed("fgsea")
    fx <- make_export_enrichment_fixture()
    # The config has run_enrichment and the library, and no Shiny setting.
    expect_identical(names(fx$config$modes$metabolomics$enrichment),
                     c("run_enrichment", "gmt_file"))

    on <- build_export(fx)
    before <- build_export_without_block(fx)

    # Same keys in the same order, plus the provenance; every other key equal.
    expect_identical(setdiff(names(on), names(before)), "enrichment_provenance")
    expect_identical(names(on)[names(on) != "enrichment_provenance"], names(before))
    others <- setdiff(names(before), c(new_keys, "payload_created_at"))
    expect_identical(on[others], before[others])
    # ... the legacy enrichment tables among them.
    expect_true(all(c("enrichment_ora", "enrichment_gsea") %in% others))
    expect_null(before$enrichment)

    enr <- on$enrichment
    expect_identical(names(enr),
                     c("available", "gene_index", "config", "manifest", "ora", "gsea",
                       "gsea_leading_edge", "gsea_rankings", "pathway_membership"))
    expect_true(enr$available)
    expect_true(validate_enrichment_payload(enr, strict = TRUE))
    expect_identical(enr$gene_index, rownames(on$expr_norm))
    expect_identical(enr$ora$combined$contrasts_without_direction$A_vs_B$ID, "PW1")
    expect_identical(enr$manifest$item, c("A_vs_B", "A_vs_B"))

    prov <- on$enrichment_provenance
    expect_identical(names(prov$ora), "A_vs_B")
    expect_identical(names(prov$gsea), "A_vs_B")
    expect_identical(prov$status$state[prov$status$analysis == "ORA"], "significant")
    expect_true(all(prov$feature_index$in_expr_norm))

    # The legacy first-contrast tables and the new block describe one result.
    # Compared column by column: the row names of the two differ, the rows do not.
    cols <- setdiff(names(on$enrichment_ora), "overlap_genes")
    expect_identical(as.list(prov$ora$A_vs_B$all_tested[cols]),
                     as.list(on$enrichment_ora[cols]))
    expect_identical(as.list(prov$gsea$A_vs_B$all_tested[names(on$enrichment_gsea)]),
                     as.list(on$enrichment_gsea))
})

test_that("the block is the same with or without the legacy keys", {
    skip_if_not_installed("fgsea")
    fx <- make_export_enrichment_fixture()
    with_legacy    <- build_export(fx, include_legacy = TRUE)
    without_legacy <- build_export(fx, include_legacy = FALSE)

    expect_identical(without_legacy[new_keys], with_legacy[new_keys])
    expect_false("enrichment_ora" %in% names(without_legacy))
})

# ---- 2. without results ------------------------------------------------------------------

test_that("a run without enrichment gets no block and no provenance", {
    skip_if_not_installed("fgsea")
    fx <- make_export_enrichment_fixture()
    none   <- build_export(fx, enrichment_res = NULL)
    before <- build_export_without_block(fx, enrichment_res = NULL)

    expect_identical(comparable(none), comparable(before))
    expect_null(none$enrichment)
    expect_true("enrichment" %in% names(none))
    expect_false("enrichment_provenance" %in% names(none))
    expect_identical(metab_shiny_enrichment_keys(NULL, rownames(fx$pre$expr_work)), list())
})

test_that("an enrichment run with no ORA or GSEA result gives an unavailable block", {
    skip_if_not_installed("fgsea")
    fx <- make_export_enrichment_fixture()
    none <- list(ora = NULL, gsea = NULL, ora_by_contrast = list(),
                 gsea_by_contrast = list())
    on <- build_export(fx, enrichment_res = none)

    expect_false(on$enrichment$available)
    expect_identical(nrow(on$enrichment$manifest), 0L)
    # Nothing is made up: no table, no ranking, no membership.
    expect_identical(on$enrichment$ora, list())
    expect_identical(on$enrichment$gsea, list())
    expect_identical(on$enrichment$gsea_rankings, list())
    expect_identical(on$enrichment$pathway_membership, list())
    expect_true(validate_enrichment_payload(on$enrichment, strict = TRUE))
    expect_identical(on$enrichment$gene_index, rownames(on$expr_norm))
    # Every contrast is on record as having no result, for both analyses.
    st <- on$enrichment_provenance$status
    expect_identical(st$contrast, c("A_vs_B", "A_vs_B"))
    expect_identical(st$state, c("no_result", "no_result"))
})

test_that("ORA without GSEA gives a block with the ORA part only", {
    skip_if_not_installed("fgsea")
    fx <- make_export_enrichment_fixture()
    only_ora <- fx$enrichment_res
    only_ora$gsea <- NULL
    only_ora$gsea_by_contrast <- list()
    on <- build_export(fx, enrichment_res = only_ora)

    expect_true(on$enrichment$available)
    expect_identical(names(on$enrichment$ora), "combined")
    expect_identical(on$enrichment$gsea, list())
    expect_identical(on$enrichment$gsea_rankings, list())
    expect_identical(on$enrichment$manifest$analysis, "ORA")
    st <- on$enrichment_provenance$status
    expect_identical(st$state[st$analysis == "GSEA"], "no_result")
    expect_true(validate_enrichment_payload(on$enrichment, strict = TRUE))
})

# ---- 3. results that cannot be encoded ------------------------------------------------------

test_that("results that cannot be encoded leave the payload whole, with the reason on record", {
    skip_if_not_installed("fgsea")
    fx <- make_export_enrichment_fixture()
    broken <- fx$enrichment_res
    broken$ora_by_contrast$A_vs_B$res$provenance$status <- "failed"
    broken$ora_by_contrast$A_vs_B$res$provenance$error <- list(step = "build",
                                                               message = "boom")

    built <- build_export_with_warnings(fx, broken)
    p <- built$payload
    # One warning says the block was left out, and why.
    about_block <- grep("enrichment block was not added", built$warnings, value = TRUE)
    expect_length(about_block, 1L)
    expect_match(about_block, "provenance:ORA:A_vs_B: it could not be built: boom",
                 fixed = TRUE)
    expect_match(about_block, "Every other key of the payload is written as usual",
                 fixed = TRUE)
    before <- build_export_without_block(fx, enrichment_res = broken)

    # No block, and every other key exactly as without it.
    expect_null(p$enrichment)
    others <- setdiff(names(before), "payload_created_at")
    expect_identical(p[others], before[others])
    expect_identical(p$enrichment_ora, fx$enrichment_res$ora$table)

    # Not dropped quietly: the provenance key says what happened.
    prov <- p$enrichment_provenance
    expect_identical(names(prov), c("schema_version", "available", "reason", "preflight"))
    expect_false(prov$available)
    expect_identical(prov$reason, "blocked")
    expect_identical(prov$preflight$check[prov$preflight$level == "blocker"],
                     "provenance:ORA:A_vs_B")
    expect_error(metab_enrichment_pathway_detail(p, "ORA", "A_vs_B", "PW1"),
                 "no enrichment block \\(reason: blocked\\)")
})

test_that("results cached before provenance existed are reported, not guessed at", {
    skip_if_not_installed("fgsea")
    fx <- make_export_enrichment_fixture()
    legacy <- fx$enrichment_res
    legacy$ora_by_contrast$A_vs_B$res$provenance  <- NULL
    legacy$gsea_by_contrast$A_vs_B$res$provenance <- NULL

    built <- build_export_with_warnings(fx, legacy)
    p <- built$payload
    expect_length(grep("enrichment block was not added", built$warnings), 1L)
    expect_null(p$enrichment)
    expect_false(p$enrichment_provenance$available)
    blocked <- p$enrichment_provenance$preflight
    expect_setequal(blocked$check[blocked$level == "blocker"],
                    c("provenance:ORA:A_vs_B", "provenance:GSEA:A_vs_B"))
    # The legacy tables are still exported as before.
    expect_identical(p$enrichment_ora, fx$enrichment_res$ora$table)
    expect_identical(p$enrichment_gsea, fx$enrichment_res$gsea$table)
})

# ---- 4. mummichog --------------------------------------------------------------------------

test_that("mummichog is exported on its own, with or without enrichment results", {
    skip_if_not_installed("fgsea")
    fx <- make_export_enrichment_fixture()
    pathways <- list(A_vs_B = data.frame(
        check.names = FALSE, stringsAsFactors = FALSE,
        pathway = c("Pathway one", "Pathway two", "Pathway three"),
        overlap_size = c(5, 3, 8), pathway_size = c(10, 6, 12),
        "p-value" = c(0.02, 0.20, 0.60)))

    with_enrichment    <- build_export(fx, mummichog_pathways = pathways)
    without_enrichment <- build_export(fx, enrichment_res = NULL,
                                       mummichog_pathways = pathways)
    no_mummichog       <- build_export(fx)

    # Mummichog without identified-metabolite enrichment.
    expect_named(without_enrichment$mummichog, "A_vs_B")
    expect_null(without_enrichment$enrichment)
    # Enrichment without mummichog.
    expect_null(no_mummichog$mummichog)
    expect_true(no_mummichog$enrichment$available)
    # Both: each is what it is alone.
    expect_named(with_enrichment$mummichog, "A_vs_B")
    expect_identical(names(with_enrichment$mummichog$A_vs_B),
                     names(without_enrichment$mummichog$A_vs_B))
    expect_identical(with_enrichment$mummichog$A_vs_B$table,
                     without_enrichment$mummichog$A_vs_B$table)
    expect_identical(with_enrichment[new_keys], no_mummichog[new_keys])
})

# ---- 5. the saved RDS ---------------------------------------------------------------------

test_that("the saved RDS carries the block and can be inspected from the file alone", {
    skip_if_not_installed("fgsea")
    fx <- make_export_enrichment_fixture()
    out_file <- file.path(withr::local_tempdir(), "shiny_payload.rds")

    saved <- suppressMessages(suppressWarnings(save_shiny_payload_metabolomics(
        pre = fx$pre, de_res = fx$de_res, inputs = fx$inputs, config = fx$config,
        enrichment_res = fx$enrichment_res, include_legacy = TRUE,
        out_file = out_file)))
    expect_identical(saved, out_file)

    p <- readRDS(out_file)
    expect_true(p$enrichment$available)
    expect_true(validate_enrichment_payload(p$enrichment, strict = TRUE))

    d <- metab_enrichment_pathway_detail(p, "ORA", "A_vs_B", "PW1")
    expect_identical(d$result$overlap, 9L)
    expect_identical(d$features$feature_id[d$features$compound_id == "C00101"],
                     c("m1a", "m1b", "m1c"))
    g <- metab_enrichment_pathway_detail(p, "GSEA", "A_vs_B", "PW1")
    expect_identical(g$compounds$representative_feature[g$compounds$compound_id == "C00101"],
                     "m1b")
    # A pathway that is not significant in ORA is there too.
    expect_false(metab_enrichment_pathway_detail(p, "ORA", "A_vs_B", "PW2")$result$significant)
})

# ---- 6. {targets} wiring ---------------------------------------------------------------------

skip_if_not_installed("targets")
# pipe_metabolomics() calls tar_target() unqualified, as _targets.R does after
# attaching {targets}; attach it here so the file runs in a fresh session.
suppressPackageStartupMessages(library(targets))

test_that("the payload target depends on what the enrichment block is built from", {
    plan <- pipe_metabolomics("pqn")
    names(plan) <- vapply(plan, function(t) t$settings$name, character(1))
    target <- plan$metab_shiny_payload
    deps <- target$command$deps
    if (length(deps) == 0) deps <- unlist(lapply(as.list(target$command$expr), all.vars))

    # The results, the DE contrasts, the matrix it exports and the config that
    # carries the switch: a change in any of them rebuilds the payload.
    expect_true(all(c("metab_enrichment_res", "metab_de_res", "metab_pre", "config") %in% deps))
    command <- target$command$expr
    if (is.list(command) || is.expression(command)) command <- command[[1]]
    expect_true("save_shiny_payload_metabolomics" %in% all.names(command))
    expect_identical(target$settings$format, "file")

    # No new target was needed, and the enrichment stage does not depend on the
    # payload: the block is assembled inside the payload target.
    expect_false(any(grepl("enrichment_payload", names(plan))))
    enr_deps <- plan$metab_enrichment_res$command$deps
    expect_false("metab_shiny_payload" %in% enr_deps)

    # Under multi-omics the per-omic payload is not built, with or without the option.
    core <- vapply(pipe_metabolomics("pqn", skip_outputs = TRUE),
                   function(t) t$settings$name, character(1))
    expect_false("metab_shiny_payload" %in% core)
})

test_that("the payload builder reaches the block builder, so {targets} tracks it", {
    # {targets} follows global functions called by a target's command. The
    # chain from the payload target to the block builder must be unbroken.
    body_names <- function(fn) all.names(body(fn))
    expect_true("build_shiny_payload_metabolomics" %in%
                    body_names(save_shiny_payload_metabolomics))
    expect_true("metab_shiny_enrichment_keys" %in%
                    body_names(build_shiny_payload_metabolomics))
    expect_true(all(c("preflight_metab_enrichment_payload",
                      "build_metab_enrichment_payload") %in%
                        body_names(metab_shiny_enrichment_keys)))
    # The option and its reader are gone, from the code and from the template.
    expect_false(exists("metab_shiny_enrichment_enabled"))
    skip_if_not_installed("yaml")
    template <- yaml::read_yaml(file.path(root_dir, "config", "templates",
                                          "metabolomics_template.yaml"))
    expect_false("shiny_payload" %in% names(template$modes$metabolomics$enrichment))
})
