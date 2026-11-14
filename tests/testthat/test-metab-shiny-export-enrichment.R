# tests/testthat/test-metab-shiny-export-enrichment.R
#
# The opt-in enrichment block of the metabolomics Shiny payload
# (modes.metabolomics.enrichment.shiny_payload).
#
#   1. off (the default): the payload is what it was, and the block's builder
#      is not even called;
#   2. on: the nine-component `enrichment` block and `enrichment_provenance`
#      are added, and no other key changes;
#   3. on without enrichment results: nothing is added, or the block says so;
#   4. on with a result that cannot be encoded: the payload stops;
#   5. the saved RDS carries the block and can be read back;
#   6. the {targets} wiring of the payload target.
#
# Everything is synthetic and lives in a temporary folder.

# ---- the synthetic run -----------------------------------------------------------
#
# Features g01..g30 are compounds C00001..C00030; m1a, m1b, m1c are three
# features of C00101. One contrast, "A_vs_B", in which g01..g08, m1a and m1c are
# significant, so that PW1 (C00001..10 + C00101) is over-represented.
make_export_enrichment_fixture <- function(shiny_payload = NULL, env = parent.frame()) {
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

    enrichment <- list(run_enrichment = TRUE, gmt_file = lib)
    if (!is.null(shiny_payload)) enrichment$shiny_payload <- shiny_payload
    config <- list(
        project = list(dir = dir), paths = list(raw = "data"),
        modes = list(metabolomics = list(
            de            = list(padj_cutoff = 0.05, p_cutoff = 0.05,
                                 linear_fc_cutoff = 1.5),
            preprocessing = list(chosen_norm = "pqn", transform = "log2", scaling = "none"),
            effects       = list(samples = "sample_id", color = "grp"),
            excel         = list(group_cv = FALSE),
            enrichment    = enrichment)))

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

with_option <- function(fx, value) {
    fx$config$modes$metabolomics$enrichment$shiny_payload <- value
    fx
}

build_export <- function(fx, enrichment_res = fx$enrichment_res, include_legacy = TRUE) {
    suppressMessages(suppressWarnings(
        build_shiny_payload_metabolomics(
            pre = fx$pre, de_res = fx$de_res, inputs = fx$inputs, config = fx$config,
            enrichment_res = enrichment_res, include_legacy = include_legacy)))
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

# ---- the switch -------------------------------------------------------------------

test_that("only a literal TRUE turns the enrichment block on", {
    cfg <- function(value) {
        list(modes = list(metabolomics = list(enrichment = list(shiny_payload = value))))
    }
    expect_true(metab_shiny_enrichment_enabled(cfg(TRUE)))
    for (value in list(FALSE, NULL, NA, "true", "yes", 1, c(TRUE, TRUE))) {
        expect_false(metab_shiny_enrichment_enabled(cfg(value)))
    }
    expect_false(metab_shiny_enrichment_enabled(list()))
    expect_false(metab_shiny_enrichment_enabled(
        list(modes = list(metabolomics = list(enrichment = list(run_enrichment = TRUE))))))

    # The shipped template leaves it off.
    template <- file.path(root_dir, "config", "templates", "metabolomics_template.yaml")
    skip_if_not_installed("yaml")
    expect_false(metab_shiny_enrichment_enabled(yaml::read_yaml(template)))
    expect_identical(
        yaml::read_yaml(template)$modes$metabolomics$enrichment$shiny_payload, FALSE)
})

# ---- 1. off ------------------------------------------------------------------------

test_that("with the option absent or false the payload is what it was", {
    skip_if_not_installed("fgsea")
    fx <- make_export_enrichment_fixture()
    # If the builder were reached at all, this would fail the test.
    local_export_stub("build_metab_enrichment_payload", function(...) {
        stop("the enrichment block must not be built when the option is off")
    })

    absent <- build_export(fx)
    off    <- build_export(with_option(fx, FALSE))

    expect_identical(comparable(off), comparable(absent))
    expect_null(absent$enrichment)
    expect_false("enrichment_provenance" %in% names(absent))
    # The canonical key is still there, as it always was, and still NULL.
    expect_true("enrichment" %in% names(absent))
    # The legacy tables are the run's own.
    expect_identical(absent$enrichment_ora, fx$enrichment_res$ora$table)
    expect_identical(absent$enrichment_gsea, fx$enrichment_res$gsea$table)
})

# ---- 2. on -------------------------------------------------------------------------

test_that("with the option on the block is added and nothing else changes", {
    skip_if_not_installed("fgsea")
    fx  <- make_export_enrichment_fixture()
    off <- build_export(fx)
    on  <- build_export(with_option(fx, TRUE))

    # Same keys in the same order, plus the provenance; every other key equal.
    expect_identical(setdiff(names(on), names(off)), "enrichment_provenance")
    expect_identical(names(on)[names(on) != "enrichment_provenance"], names(off))
    others <- setdiff(names(off), c(new_keys, "payload_created_at"))
    expect_identical(on[others], off[others])
    # ... the legacy enrichment tables among them.
    expect_true(all(c("enrichment_ora", "enrichment_gsea") %in% others))

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
    fx <- with_option(make_export_enrichment_fixture(), TRUE)
    with_legacy    <- build_export(fx, include_legacy = TRUE)
    without_legacy <- build_export(fx, include_legacy = FALSE)

    expect_identical(without_legacy[new_keys], with_legacy[new_keys])
    expect_false("enrichment_ora" %in% names(without_legacy))
})

# ---- 3. on, without results ------------------------------------------------------------

test_that("with enrichment disabled the option adds nothing", {
    skip_if_not_installed("fgsea")
    fx  <- make_export_enrichment_fixture()
    off <- build_export(fx, enrichment_res = NULL)
    on  <- build_export(with_option(fx, TRUE), enrichment_res = NULL)

    expect_identical(comparable(on), comparable(off))
    expect_null(on$enrichment)
    expect_false("enrichment_provenance" %in% names(on))
})

test_that("an enrichment run with no ORA or GSEA result gives an unavailable block", {
    skip_if_not_installed("fgsea")
    fx <- with_option(make_export_enrichment_fixture(), TRUE)
    none <- list(ora = NULL, gsea = NULL, ora_by_contrast = list(),
                 gsea_by_contrast = list())
    on <- build_export(fx, enrichment_res = none)

    expect_false(on$enrichment$available)
    expect_identical(nrow(on$enrichment$manifest), 0L)
    expect_true(validate_enrichment_payload(on$enrichment, strict = TRUE))
    expect_identical(on$enrichment$gene_index, rownames(on$expr_norm))
    # Every contrast is on record as having no result, for both analyses.
    st <- on$enrichment_provenance$status
    expect_identical(st$contrast, c("A_vs_B", "A_vs_B"))
    expect_identical(st$state, c("no_result", "no_result"))
})

test_that("ORA without GSEA gives a block with the ORA part only", {
    skip_if_not_installed("fgsea")
    fx <- with_option(make_export_enrichment_fixture(), TRUE)
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

# ---- 4. on, with a result that cannot be encoded ----------------------------------------

test_that("a result without usable provenance stops the payload when the option is on", {
    skip_if_not_installed("fgsea")
    fx <- make_export_enrichment_fixture()
    broken <- fx$enrichment_res
    broken$ora_by_contrast$A_vs_B$res$provenance$status <- "failed"
    broken$ora_by_contrast$A_vs_B$res$provenance$error <- list(step = "build",
                                                               message = "boom")

    expect_error(build_export(with_option(fx, TRUE), enrichment_res = broken),
                 "provenance:ORA:A_vs_B: it could not be built: boom")
    # With the option off the same run builds its payload as before.
    off <- build_export(fx, enrichment_res = broken)
    expect_null(off$enrichment)
    expect_identical(off$enrichment_ora, fx$enrichment_res$ora$table)
})

# ---- 5. the saved RDS ---------------------------------------------------------------------

test_that("the saved RDS carries the block and can be inspected from the file alone", {
    skip_if_not_installed("fgsea")
    fx <- with_option(make_export_enrichment_fixture(), TRUE)
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
    expect_true(all(c("metab_shiny_enrichment_enabled", "build_metab_enrichment_payload") %in%
                        body_names(build_shiny_payload_metabolomics)))
})
