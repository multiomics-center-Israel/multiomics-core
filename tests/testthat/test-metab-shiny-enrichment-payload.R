# tests/testthat/test-metab-shiny-enrichment-payload.R
#
# The metabolomics ORA and GSEA results as the Shiny app's enrichment block
# (R/domain/metabolomics/07b_shiny_enrichment_payload.R).
#
# The block is a re-encoding of results and provenance. These tests cover:
#   1. the nine components, validated by the RNA-seq validator;
#   2. gene_index: the exported matrix first, other features flagged;
#   3. ORA: compound counts and feature lists kept apart;
#   4. GSEA: rankings and leading edges in representative features only;
#   5. significant-only display tables, every tested pathway in provenance;
#   6. ran-and-empty versus no result;
#   7. the pathway detail reader, for significant and non-significant pathways;
#   8. the preflight and its blockers.
#
# Everything is synthetic and lives in a temporary folder.

# ---- the synthetic case ---------------------------------------------------------
#
# Features g01..g30 are compounds C00001..C00030. Besides them:
#   m1a m1b m1c  three features of C00101; m1a and m1c are significant, m1b has
#                the largest |logFC| and is not.
#   d1           C00108, measured, absent from every DE table.
#   n1           name only; dropped because the data carries KEGG ids.
# One library: PW1 = C00001..10 + C00101, PW2 = C00011..20 + C00108,
# PW3 = C00021..30.
# Contrasts:
#   "A vs B"  g01..g08, m1a, m1c significant: PW1 is strongly over-represented.
#   "C vs D"  nothing significant: ORA has no result; GSEA still runs.
#   "E vs F"  g01, g11, g21 significant: ORA runs and finds nothing.
# The exported matrix lacks g30, which the enrichment still refers to.
make_payload_case <- function(second_library = character(0), rename = NULL,
                              env = parent.frame()) {
    dir <- withr::local_tempdir(.local_envir = env)
    base_feat <- sprintf("g%02d", 1:30)
    base_kegg <- sprintf("C%05d", 1:30)
    feature_id <- c(base_feat, "m1a", "m1b", "m1c", "d1", "n1")
    row_data <- data.frame(
        feature_id = feature_id,
        Name = c(paste0("cmpd", 1:30), "multiA", "multiB", "multiC", "NoDE", "Citrate"),
        KEGG = c(base_kegg, "C00101", "C00101", "C00101", "C00108", NA),
        stringsAsFactors = FALSE)
    if (!is.null(rename)) {
        row_data$feature_id[row_data$feature_id == names(rename)] <- unname(rename)
    }
    ids <- row_data$feature_id
    expr <- withr::with_seed(1, matrix(
        stats::rnorm(nrow(row_data) * 6, mean = 10), nrow = nrow(row_data),
        dimnames = list(ids, paste0("S", 1:6))))

    i <- 1:30
    de_ids <- ids[ids != "d1"]
    small <- (-1)^i * (0.2 + 0.01 * i)
    fc_ab <- c(ifelse(i <= 10, 2 + 0.01 * i, small), 2.5, -3.0, 2.6, 0.9)
    fc_ef <- c((-1)^i * (1 + 0.01 * i), 2.5, -3.0, 2.6, 0.9)
    de_table <- function(logfc, significant) {
        data.frame(feature_id = de_ids, logFC = logfc,
                   P.Value = ifelse(de_ids %in% significant, 0.001, 0.5),
                   adj.P.Val = ifelse(de_ids %in% significant, 0.001, 0.5),
                   stringsAsFactors = FALSE)
    }
    g <- function(n) ids[n]
    de_res <- list(de_tables = list(
        "A vs B" = de_table(fc_ab, c(g(1:8), g(31), g(33))),
        "C vs D" = de_table(fc_ab, character(0)),
        "E vs F" = de_table(fc_ef, c(g(1), g(11), g(21)))))

    gmt_line <- function(id, desc, members) paste(c(id, desc, members), collapse = "\t")
    lib <- file.path(dir, "libA.gmt")
    writeLines(c(
        gmt_line("PW1", "First",  c(base_kegg[1:10], "C00101")),
        gmt_line("PW2", "Second", c(base_kegg[11:20], "C00108")),
        gmt_line("PW3", "Third",  base_kegg[21:30])), lib)
    gmt_files <- list(lib)
    if (length(second_library) > 0) {
        lib_b <- file.path(dir, "libB.gmt")
        writeLines(second_library, lib_b)
        gmt_files <- list(lib, lib_b)
    }

    config <- list(
        project = list(dir = dir), paths = list(raw = "data"),
        modes = list(metabolomics = list(
            de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5),
            enrichment = list(run_enrichment = TRUE, gmt_file = gmt_files))))
    pre <- list(row_data = row_data, expr_raw = expr)

    quiet <- function(expr) suppressWarnings(suppressMessages(expr))
    one <- function(fun, contrast, seed = 11) {
        res <- quiet(withr::with_seed(seed, fun(pre, de_res, config, contrast = contrast)))
        if (is.null(res)) NULL else list(res = res)
    }
    contrasts <- names(de_res$de_tables)
    by <- function(fun) {
        out <- lapply(stats::setNames(contrasts, contrasts), function(ct) one(fun, ct))
        out[!vapply(out, is.null, logical(1))]
    }
    list(
        enrichment_res = list(ora_by_contrast  = by(run_metabolomics_ora),
                              gsea_by_contrast = by(run_metabolomics_gsea)),
        expr_features = ids[ids != ids[30]],
        contrasts = contrasts, ids = ids)
}

build_case_payload <- function(case) {
    build_metab_enrichment_payload(case$enrichment_res, case$expr_features,
                                   contrasts = case$contrasts)
}

decode <- function(payload, idx) payload$enrichment$gene_index[idx]

# ---- 1. structure ------------------------------------------------------------------

test_that("the block has exactly the nine RNA-seq components and validates", {
    skip_if_not_installed("fgsea")
    out <- build_case_payload(make_payload_case())

    expect_identical(names(out), c("enrichment", "enrichment_provenance"))
    expect_identical(names(out$enrichment),
                     c("available", "gene_index", "config", "manifest", "ora", "gsea",
                       "gsea_leading_edge", "gsea_rankings", "pathway_membership"))
    expect_true(out$enrichment$available)
    expect_true(validate_enrichment_payload(out$enrichment, strict = TRUE))

    expect_identical(names(out$enrichment$ora), "combined")
    expect_identical(names(out$enrichment$ora$combined), "contrasts_without_direction")
    expect_identical(names(out$enrichment$gsea$combined), "fc")
    expect_identical(names(out$enrichment$gsea_rankings), "fc")
    expect_identical(out$enrichment$config$databases_run, "combined")
    expect_identical(out$enrichment$config$rankings, "fc")
    expect_identical(out$enrichment$config$pvalue_cutoff, .METAB_ENRICHMENT_SIG_FDR)
    expect_identical(out$enrichment$config$libraries, "libA")

    expect_identical(names(out$enrichment_provenance),
                     c("schema_version", "id_type", "database", "ranking",
                       "significance", "status", "feature_index", "feature_compounds",
                       "reference_files", "ora", "gsea", "preflight"))
    expect_identical(out$enrichment_provenance$significance,
                     list(column = "FDR", cutoff = .METAB_ENRICHMENT_SIG_FDR))
})

test_that("no enrichment gives no block, and no result gives an unavailable one", {
    expect_identical(build_metab_enrichment_payload(NULL, c("f1", "f2")),
                     list(enrichment = NULL, enrichment_provenance = NULL))

    out <- build_metab_enrichment_payload(
        list(ora_by_contrast = list(), gsea_by_contrast = list()),
        c("f1", "f2"), contrasts = "A vs B")
    expect_false(out$enrichment$available)
    expect_identical(nrow(out$enrichment$manifest), 0L)
    expect_true(validate_enrichment_payload(out$enrichment, strict = TRUE))
    expect_identical(out$enrichment$gene_index, c("f1", "f2"))
    expect_identical(out$enrichment_provenance$status$state, c("no_result", "no_result"))
})

# ---- 2. gene_index -------------------------------------------------------------------

test_that("gene_index starts with the exported matrix, and flags everything else", {
    skip_if_not_installed("fgsea")
    case <- make_payload_case()
    out <- build_case_payload(case)
    gi <- out$enrichment$gene_index
    n <- length(case$expr_features)

    expect_identical(gi[seq_len(n)], case$expr_features)
    expect_false(anyDuplicated(gi) > 0)
    # g30 is not a row of the exported matrix, and the enrichment refers to it.
    expect_identical(gi[-seq_len(n)], "g30")

    fi <- out$enrichment_provenance$feature_index
    expect_identical(fi$feature_id, gi)
    expect_identical(fi$feature_id[!fi$in_expr_norm], "g30")
    expect_identical(sum(fi$in_expr_norm), n)
})

# ---- 3. ORA ----------------------------------------------------------------------------

test_that("ORA keeps compound counts and feature lists apart", {
    skip_if_not_installed("fgsea")
    case <- make_payload_case()
    out <- build_case_payload(case)
    tab <- out$enrichment$ora$combined$contrasts_without_direction[["A vs B"]]

    expect_identical(names(tab), c("ID", "Description", "GeneRatio", "BgRatio", "pvalue",
                                   "p.adjust", "geneID", "Count", "FoldEnrichment",
                                   "n_features", "library"))
    expect_identical(tab$ID, "PW1")
    expect_identical(tab$Description, "First")
    # 9 foreground compounds, all in PW1; PW1 has 11 of the 32 universe compounds.
    expect_identical(tab$Count, 9L)
    expect_identical(tab$GeneRatio, "9/9")
    expect_identical(tab$BgRatio, "11/32")
    # Ten features: C00101 contributes two significant features to one compound.
    expect_identical(tab$n_features, 10L)
    expect_setequal(strsplit(tab$geneID, "/", fixed = TRUE)[[1]],
                    c(sprintf("g%02d", 1:8), "m1a", "m1c"))
    # The non-significant sibling is not in the list.
    expect_false("m1b" %in% strsplit(tab$geneID, "/", fixed = TRUE)[[1]])

    # The statistics are the result table's own.
    res <- case$enrichment_res$ora_by_contrast[["A vs B"]]$res$table
    row <- res[res$pathway == "PW1 - First", ]
    expect_identical(tab$pvalue, row$raw_p)
    expect_identical(tab$p.adjust, row$FDR)
    expect_identical(tab$FoldEnrichment, row$fold_enrichment)
})

# ---- 4. GSEA ----------------------------------------------------------------------------

test_that("the GSEA ranking round-trips through gene_index in representative features", {
    skip_if_not_installed("fgsea")
    case <- make_payload_case()
    out <- build_case_payload(case)

    for (contrast in names(case$enrichment_res$gsea_by_contrast)) {
        res <- case$enrichment_res$gsea_by_contrast[[contrast]]$res
        rk  <- out$enrichment$gsea_rankings$fc[[contrast]]
        ranked <- res$provenance$ranked

        expect_identical(names(rk), c("idx", "score"))
        expect_identical(decode(out, rk$idx), ranked$representative_feature)
        expect_identical(rk$score, unname(res$ranks))
        # One feature per compound: the siblings are not in the ranking.
        expect_true("m1b" %in% decode(out, rk$idx))
        expect_false(any(c("m1a", "m1c") %in% decode(out, rk$idx)))
    }
})

test_that("GSEA tables and leading edges describe the significant pathways only", {
    skip_if_not_installed("fgsea")
    case <- make_payload_case()
    out <- build_case_payload(case)

    for (contrast in names(case$enrichment_res$gsea_by_contrast)) {
        res <- case$enrichment_res$gsea_by_contrast[[contrast]]$res
        tab <- out$enrichment$gsea$combined$fc[[contrast]]
        le  <- out$enrichment$gsea_leading_edge$combined$fc[[contrast]]
        sig <- res$table[!is.na(res$table$FDR) & res$table$FDR < .METAB_ENRICHMENT_SIG_FDR, ]

        expect_identical(names(tab), c("ID", "Description", "setSize", "enrichmentScore",
                                       "NES", "pvalue", "p.adjust", "library"))
        expect_false("core_enrichment" %in% names(tab))
        expect_identical(nrow(tab), nrow(sig))
        expect_identical(tab$NES, sig$NES)
        expect_identical(tab$enrichmentScore, sig$ES)
        expect_identical(tab$p.adjust, sig$FDR)
        expect_identical(names(le), tab$ID)

        ranked <- res$provenance$ranked
        for (r in seq_len(nrow(sig))) {
            compounds <- strsplit(sig$leading_edge[r], ";", fixed = TRUE)[[1]]
            expect_identical(
                decode(out, le[[r]]),
                ranked$representative_feature[match(compounds, ranked$compound_id)])
        }
        # A leading edge never holds a feature that was not ranked.
        expect_false(any(c("m1a", "m1c", "d1") %in%
                             decode(out, unlist(le, use.names = FALSE))))
    }
})

test_that("pathway membership lists every measured feature of the pathway compounds", {
    skip_if_not_installed("fgsea")
    out <- build_case_payload(make_payload_case())
    pm <- out$enrichment$pathway_membership$combined

    expect_setequal(names(pm), c("PW1", "PW2", "PW3"))
    # All three features of C00101, representative or not.
    expect_setequal(decode(out, pm$PW1), c(sprintf("g%02d", 1:10), "m1a", "m1b", "m1c"))
    # d1 is measured and in no DE table: it is a member all the same.
    expect_setequal(decode(out, pm$PW2), c(sprintf("g%02d", 11:20), "d1"))
    expect_setequal(decode(out, pm$PW3), sprintf("g%02d", 21:30))
})

# ---- 5. display tables versus every tested pathway ----------------------------------------

test_that("every tested pathway is kept in provenance with its original statistics", {
    skip_if_not_installed("fgsea")
    case <- make_payload_case()
    out <- build_case_payload(case)
    prov <- out$enrichment_provenance

    for (method in c("ora", "gsea")) {
        by_contrast <- case$enrichment_res[[paste0(method, "_by_contrast")]]
        expect_identical(names(prov[[method]]), names(by_contrast))
        for (contrast in names(by_contrast)) {
            res <- by_contrast[[contrast]]$res$table
            all_tested <- prov[[method]][[contrast]]$all_tested

            expect_identical(nrow(all_tested), nrow(res))
            for (col in setdiff(names(res), "overlap_genes")) {
                expect_identical(all_tested[[col]], res[[col]])
            }
            expect_identical(all_tested$pathway_id, sub(" - .*$", "", res$pathway))
            expect_identical(all_tested$significant,
                             !is.na(res$FDR) & res$FDR < .METAB_ENRICHMENT_SIG_FDR)

            # The display table is the significant subset, in the same order.
            display <- if (method == "ora") {
                out$enrichment$ora$combined$contrasts_without_direction[[contrast]]
            } else {
                out$enrichment$gsea$combined$fc[[contrast]]
            }
            expect_identical(display$ID, all_tested$pathway_id[all_tested$significant])
        }
    }

    # Non-significant ORA pathways are still explorable.
    ab <- prov$ora[["A vs B"]]$all_tested
    expect_identical(nrow(ab), 3L)
    expect_identical(sum(ab$significant), 1L)
    expect_true(all(c("raw_p", "FDR", "fold_enrichment", "overlap", "pathway_size",
                      "fg_size", "bg_size", "library", "feature_ids", "n_features") %in%
                        names(ab)))
})

# ---- 6. ran-and-empty versus no result ------------------------------------------------------

test_that("an analysis that found nothing is told apart from one with no result", {
    skip_if_not_installed("fgsea")
    case <- make_payload_case()
    out <- build_case_payload(case)
    m <- out$enrichment$manifest
    st <- out$enrichment_provenance$status
    ora_m <- m[m$analysis == "ORA", ]

    expect_setequal(m$analysis, c("ORA", "GSEA"))
    expect_identical(unique(m$database), "combined")
    expect_identical(unique(ora_m$group), "contrasts")
    expect_identical(unique(ora_m$subgroup), "all")
    expect_identical(unique(m$group[m$analysis == "GSEA"]), "fc")
    # storage_key names the unit's place in this block, not a directory on disk.
    expect_identical(ora_m$storage_key[ora_m$item == "A vs B"],
                     "ora/combined/contrasts_without_direction/A vs B")

    # Ran and significant.
    ab <- ora_m[ora_m$item == "A vs B", ]
    expect_identical(c(ab$status, as.character(ab$n_significant)), c("significant", "1"))
    # Ran and found nothing: a manifest row, and a table with no rows.
    ef <- ora_m[ora_m$item == "E vs F", ]
    expect_identical(c(ef$status, as.character(ef$n_significant)), c("empty", "0"))
    expect_true(ef$evaluated)
    empty <- out$enrichment$ora$combined$contrasts_without_direction[["E vs F"]]
    expect_identical(nrow(empty), 0L)
    expect_identical(names(empty)[1:2], c("ID", "Description"))
    # No result: no manifest row and no table.
    expect_false("C vs D" %in% ora_m$item)
    expect_null(out$enrichment$ora$combined$contrasts_without_direction[["C vs D"]])

    state <- function(analysis, contrast) {
        st$state[st$analysis == analysis & st$contrast == contrast]
    }
    expect_identical(state("ORA", "A vs B"), "significant")
    expect_identical(state("ORA", "E vs F"), "no_significant_pathways")
    expect_identical(state("ORA", "C vs D"), "no_result")
    expect_true(state("GSEA", "C vs D") %in% c("significant", "no_significant_pathways"))
    expect_identical(nrow(st), 6L)
    expect_identical(st$n_tested[st$analysis == "ORA" & st$contrast == "E vs F"], 3L)
    expect_true(is.na(st$n_tested[st$analysis == "ORA" & st$contrast == "C vs D"]))
})

# ---- 7. the detail reader ---------------------------------------------------------------------

test_that("the detail reader shows compounds and features behind an ORA pathway", {
    skip_if_not_installed("fgsea")
    out <- build_case_payload(make_payload_case())

    d <- metab_enrichment_pathway_detail(out, "ORA", "A vs B", "PW1")
    expect_identical(names(d), c("analysis", "contrast", "pathway_id", "result",
                                 "compounds", "features"))
    expect_identical(nrow(d$result), 1L)
    expect_true(d$result$significant)
    expect_identical(d$result$overlap, 9L)

    expect_identical(nrow(d$compounds), 11L)
    multi <- d$compounds[d$compounds$compound_id == "C00101", ]
    expect_identical(c(multi$n_features, multi$n_significant_features), c(3L, 2L))
    expect_true(multi$in_foreground)

    f <- d$features
    expect_identical(nrow(f), 13L)
    m <- f[f$compound_id == "C00101", ]
    expect_identical(m$feature_id, c("m1a", "m1b", "m1c"))
    expect_identical(m$feature_significant, c(TRUE, FALSE, TRUE))
    expect_identical(m$compound_in_foreground, rep(TRUE, 3))
    expect_identical(unique(f$id_source), "kegg_column")
    # In the universe and not in the foreground.
    g9 <- f[f$feature_id == "g09", ]
    expect_identical(c(g9$compound_in_universe, g9$feature_significant,
                       g9$compound_in_foreground, g9$compound_in_overlap),
                     c(TRUE, FALSE, FALSE, FALSE))

    # A pathway that is NOT significant is just as inspectable.
    ns <- metab_enrichment_pathway_detail(out, "ora", "A vs B", "PW2")
    expect_false(ns$result$significant)
    expect_identical(ns$result$overlap, 0L)
    expect_identical(nrow(ns$compounds), 11L)
    expect_true("d1" %in% ns$features$feature_id)
    expect_false(any(ns$features$compound_in_overlap))
    expect_false(ns$features$in_de_table[ns$features$feature_id == "d1"])
})

test_that("the detail reader shows representatives and the leading edge of a GSEA pathway", {
    skip_if_not_installed("fgsea")
    out <- build_case_payload(make_payload_case())

    d <- metab_enrichment_pathway_detail(out, "GSEA", "A vs B", "PW1")
    expect_identical(nrow(d$result), 1L)
    expect_true(all(c("NES", "ES", "pval", "FDR", "pathway_size", "leading_edge",
                      "leading_edge_features", "significant") %in% names(d$result)))

    multi <- d$compounds[d$compounds$compound_id == "C00101", ]
    expect_identical(multi$representative_feature, "m1b")
    expect_identical(multi$logFC, -3.0)
    expect_true(multi$in_effective_statistic)

    f <- d$features
    m <- f[f$compound_id == "C00101", ]
    expect_identical(m$feature_id, c("m1a", "m1b", "m1c"))
    expect_identical(m$role, c("not_representative", "representative", "not_representative"))
    expect_identical(m$feature_is_representative, c(FALSE, TRUE, FALSE))
    expect_identical(m$contributes_to_score, c(FALSE, TRUE, FALSE))
    expect_identical(m$in_leading_edge[c(1, 3)], c(FALSE, FALSE))
    expect_identical(m$logFC, c(2.5, -3.0, 2.6))
    # The leading edge on record is the one in the result row.
    expect_setequal(f$feature_id[f$in_leading_edge],
                    strsplit(d$result$leading_edge_features, "/", fixed = TRUE)[[1]])

    # A measured compound that was not ranked: shown, credited nothing.
    pw2 <- metab_enrichment_pathway_detail(out, "GSEA", "C vs D", "PW2")
    unranked <- pw2$features[pw2$features$feature_id == "d1", ]
    expect_identical(unranked$role, "not_in_de_table")
    expect_identical(c(unranked$compound_ranked, unranked$contributes_to_score,
                       unranked$in_leading_edge), c(FALSE, FALSE, FALSE))
    c108 <- pw2$compounds[pw2$compounds$member_as_listed == "C00108", ]
    expect_identical(c(c108$in_input_set, c108$measured), c(FALSE, TRUE))
    expect_true(is.na(c108$rank))
})

test_that("the detail reader says what is missing", {
    skip_if_not_installed("fgsea")
    out <- build_case_payload(make_payload_case())

    expect_error(metab_enrichment_pathway_detail(out, "ORA", "C vs D", "PW1"),
                 "no ORA result for contrast 'C vs D' \\(state: no_result\\)")
    expect_error(metab_enrichment_pathway_detail(out, "ORA", "A vs B", "PW9"),
                 "pathway 'PW9' was not tested")
    expect_error(metab_enrichment_pathway_detail(out, "QEA", "A vs B", "PW1"),
                 "analysis must be")
    expect_error(metab_enrichment_pathway_detail(list(), "ORA", "A vs B", "PW1"),
                 "no enrichment_provenance")
})

# ---- 8. preflight ---------------------------------------------------------------------------

test_that("the preflight passes a clean run and describes it", {
    skip_if_not_installed("fgsea")
    case <- make_payload_case()
    pf <- preflight_metab_enrichment_payload(case$enrichment_res, case$expr_features,
                                             contrasts = case$contrasts)

    expect_identical(names(pf), c("check", "level", "detail"))
    expect_true(all(pf$level %in% c("ok", "info", "blocker")))
    expect_false(any(pf$level == "blocker"))
    expect_match(pf$detail[pf$check == "features_outside_expr_norm"], "^1 feature")
    expect_match(pf$detail[pf$check == "units_without_result"], "ORA:C vs D")
    # Significant and without a usable compound: the name-only feature is not
    # significant here, so none.
    expect_match(pf$detail[pf$check == "significant_features_without_compound:ORA:A vs B"],
                 "^0 significant")
})

test_that("a feature id containing '/' blocks the payload", {
    skip_if_not_installed("fgsea")
    case <- make_payload_case(rename = c(g05 = "g05/alt"))
    pf <- preflight_metab_enrichment_payload(case$enrichment_res, case$expr_features)
    expect_identical(pf$level[pf$check == "feature_ids_without_slash"], "blocker")
    expect_error(build_case_payload(case), "feature_ids_without_slash.*g05/alt")
})

test_that("two tested pathways with one id block the payload", {
    skip_if_not_installed("fgsea")
    # The same id in a second library, with another description, is another
    # label: both are tested, and both would be "PW1" to the app.
    twin <- paste(c("PW1", "Other", sprintf("C%05d", 21:25)), collapse = "\t")
    case <- make_payload_case(second_library = twin)
    pf <- preflight_metab_enrichment_payload(case$enrichment_res, case$expr_features)
    expect_true(any(pf$level == "blocker" & grepl("^pathway_ids:", pf$check)))
    expect_error(build_case_payload(case), "two tested pathways share an id")
})

test_that("a result without usable provenance blocks the payload", {
    skip_if_not_installed("fgsea")
    case <- make_payload_case()
    case$enrichment_res$ora_by_contrast[["A vs B"]]$res$provenance$status <- "inconsistent"
    expect_error(build_case_payload(case), "provenance:ORA:A vs B")

    case <- make_payload_case()
    case$enrichment_res$gsea_by_contrast[["A vs B"]]$res$provenance <- NULL
    expect_error(build_case_payload(case),
                 "provenance:GSEA:A vs B: there is no provenance record")
})

test_that("duplicated feature ids in the exported matrix block the payload", {
    skip_if_not_installed("fgsea")
    case <- make_payload_case()
    case$expr_features <- c(case$expr_features, case$expr_features[1])
    expect_error(build_case_payload(case), "expr_feature_ids: duplicated feature ids")
})
