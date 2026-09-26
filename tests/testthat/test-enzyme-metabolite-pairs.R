# What the enzyme-metabolite pair table does with its annotation.
#
# Everything here is synthetic: invented protein ids, invented metabolite names
# and a handful of KEGG-shaped identifiers. No KEGG call is made — the four
# link tables, the compound-pathway table and the reaction equations are all
# stubbed, because what is under test is the joining and filtering, not what
# KEGG returns.
#
# The functions are sourced into the global environment by tar_source(), not
# exported from a package, so with_mocked_bindings() has no namespace to
# rebind. Assign into the function's own environment and restore on exit.

local_stubs <- function(stubs, env = parent.frame()) {
    target <- environment(build_enzyme_metabolite_pairs)
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

# Two enzymes and three metabolites, wired so that every filter has something
# to act on: P1 (EC 1.1.1.1) acts on C00100 and on the currency compound
# C00002; P2 (EC 2.2.2.2) acts on C00200, which sits in no pathway P2 is in;
# P3 carries no EC at all. C00300 is a compound no measured metabolite has.
prot_de <- function(hits = c(TRUE, TRUE, TRUE)) {
    list(tables = NULL, summary_df = NULL,
         .std = data.frame(
             feature_id = c("P1", "P2", "P3"),
             log2fc = c(1.5, -2, 0.1),
             pvalue = c(0.001, 0.002, 0.9),
             padj = ifelse(hits, c(0.01, 0.02, 0.9), c(0.5, 0.6, 0.9)),
             stringsAsFactors = FALSE))
}

metab_de <- function() {
    list(de_tables = NULL,
         .std = data.frame(
             feature_id = c("alpha", "beta"),
             log2fc = c(2, -1),
             pvalue = c(0.004, 0.2),
             padj = c(0.04, 0.3),
             stringsAsFactors = FALSE))
}

# extract_de_tables() is stubbed to hand back the fixtures above under one
# contrast per layer, spelled differently on purpose: the two layers are paired
# on the canonical contrast key, not on the literal name.
stub_extract <- function() {
    function(de_data, omics_type, harmonization_res = NULL) {
        if (identical(omics_type, "proteomics")) {
            list(A_vs_B = de_data$.std)
        } else {
            list(`A vs. B` = de_data$.std)
        }
    }
}

stub_annotation <- function(link_overrides = list()) {
    links <- list(
        # gene -> EC; P3's gene is deliberately absent
        enzyme = data.frame(from = c("tst:1", "tst:2"),
                            to = c("1.1.1.1", "2.2.2.2"),
                            stringsAsFactors = FALSE),
        # EC -> reaction
        reaction = data.frame(from = c("1.1.1.1", "1.1.1.1", "2.2.2.2"),
                              to = c("R00001", "R00002", "R00003"),
                              stringsAsFactors = FALSE),
        # reaction -> compound
        compound = data.frame(
            from = c("R00001", "R00001", "R00002", "R00003", "R00003"),
            to = c("C00100", "C00002", "C00100", "C00200", "C00300"),
            stringsAsFactors = FALSE),
        # gene -> pathway
        pathway = data.frame(from = c("tst:1", "tst:2"),
                             to = c("tst00010", "tst00020"),
                             stringsAsFactors = FALSE)
    )
    links <- utils::modifyList(links, link_overrides)

    list(
        extract_de_tables = stub_extract(),
        resolve_kegg_org_code = function(...) "tst",
        get_organism_db = function(...) "OrgDb.stub",
        map_feature_ids_to_entrez = function(...) data.frame(
            feature_id = c("P1", "P2", "P3"), ENTREZID = c("1", "2", "3"),
            stringsAsFactors = FALSE),
        convert_entrez_to_kegg = function(...) c("1" = "tst:1", "2" = "tst:2",
                                                 "3" = "tst:3"),
        kegg_link_table = function(target, source, cache_dir = NULL) links[[target]],
        map_metabolite_ids_to_kegg = function(...) data.frame(
            feature_id = c("alpha", "beta"), KEGG_CPD = c("C00100", "C00200"),
            stringsAsFactors = FALSE),
        get_kegg_compound_pathways = function(...) data.frame(
            compound = c("C00100", "C00200", "C00002"),
            pathway = c("map00010", "map00030", "map00010"),
            name = c("Alpha pathway", "Beta pathway", "Alpha pathway"),
            stringsAsFactors = FALSE),
        fetch_reaction_roles = function(reaction_ids, cache_dir = NULL) data.frame(
            reaction = c("R00001", "R00002"), compound = c("C00100", "C00100"),
            role = c("substrate", "product"), stringsAsFactors = FALSE),
        protein_group_gene_symbols = function(...) c(P1 = "Aaa", P2 = "Bbb",
                                                     P3 = "Ccc"),
        kegg_compound_names = function(cache_dir = NULL) c(C00100 = "Alphanoate")
    )
}

run_pairs <- function(overrides = list(), config = list(), quiet = TRUE,
                      env = parent.frame()) {
    stubs <- utils::modifyList(stub_annotation(), overrides)
    local_stubs(stubs, env = env)
    # The table is opt-in; these tests exercise it switched on unless a test
    # says otherwise.
    cfg <- utils::modifyList(
        list(global = list(organism = "Test organism"),
             modes = list(multiomics = list(enrichment = list(
                 enzyme_metabolite = list(enabled = TRUE))))),
        config)
    call <- function() build_enzyme_metabolite_pairs(
        de_results = list(proteomics = prot_de(), metabolomics = metab_de()),
        harmonization_res = list(), config = cfg,
        out_dir = withr::local_tempdir(.local_envir = env))
    # The message tests need the messages through, the rest do not want them.
    if (quiet) suppressMessages(call()) else call()
}


# Rows with a compound; the table also lists enzymes nothing could be paired
# with, which the tests below cover on their own.
only_pairs <- function(x) x[!is.na(x$compound), , drop = FALSE]

test_that("an enzyme is paired with every measured metabolite it acts on", {
    pairs <- only_pairs(run_pairs())

    expect_s3_class(pairs, "data.frame")
    # P1 -> C00100 (alpha) survives; P1 -> C00002 is a currency compound,
    # P2 -> C00200 shares no pathway, P2 -> C00300 was never measured, and P3
    # has no EC.
    expect_identical(pairs$protein, "P1")
    expect_identical(pairs$metabolite, "alpha")
    expect_identical(pairs$ec, "1.1.1.1")
    expect_identical(pairs$gene_symbol, "Aaa")
    expect_identical(pairs$contrast, "A_vs_B")
})

test_that("the reactions behind one pair are kept in one row", {
    pairs <- only_pairs(run_pairs())

    # R00001 and R00002 both link EC 1.1.1.1 to C00100: one pair, not two rows.
    expect_identical(nrow(pairs), 1L)
    expect_identical(pairs$reaction, "R00001;R00002")
    # The two reactions disagree on the role, and the row says so.
    expect_identical(pairs$role, "product;substrate")
})

test_that("a metabolite the enzyme does not share a pathway with is dropped", {
    with_pathway <- run_pairs(list(get_kegg_compound_pathways = function(...) {
        data.frame(compound = c("C00100", "C00200"),
                   pathway = c("map00010", "map00020"),
                   name = c("Alpha pathway", "Beta pathway"),
                   stringsAsFactors = FALSE)
    }))
    # C00200 now shares tst00020 with P2, so beta joins alpha.
    expect_setequal(with_pathway$metabolite, c("alpha", "beta"))

    relaxed <- run_pairs(config = list(modes = list(multiomics = list(
        enrichment = list(enzyme_metabolite = list(require_shared_pathway = FALSE))))))
    # Without the filter, P2 -> C00200 comes back even with no shared pathway.
    expect_setequal(relaxed$metabolite, c("alpha", "beta"))
})

test_that("currency metabolites are dropped unless the config keeps them", {
    # ATP has to be a measured metabolite as well as a mapped one: a pair is
    # kept only where the metabolomics layer has a DE row for it.
    measured_atp <- list(
        map_metabolite_ids_to_kegg = function(...) data.frame(
            feature_id = c("alpha", "atp"), KEGG_CPD = c("C00100", "C00002"),
            stringsAsFactors = FALSE),
        extract_de_tables = function(de_data, omics_type, harmonization_res = NULL) {
            if (identical(omics_type, "proteomics")) return(list(A_vs_B = de_data$.std))
            list(`A vs. B` = data.frame(
                feature_id = c("alpha", "atp"), log2fc = c(2, -1),
                pvalue = c(0.004, 0.2), padj = c(0.04, 0.3),
                stringsAsFactors = FALSE))
        })

    kept <- run_pairs(measured_atp,
        config = list(modes = list(multiomics = list(enrichment = list(
            enzyme_metabolite = list(drop_currency_metabolites = FALSE))))))
    expect_true("atp" %in% kept$metabolite)

    dropped <- run_pairs(measured_atp)
    expect_false("atp" %in% dropped$metabolite)
})

test_that("only enzymes that changed are paired, unless the config says otherwise", {
    # P2's own pair needs no shared-pathway filter to survive, so relax it and
    # then read which enzymes appear.
    relax <- list(modes = list(multiomics = list(enrichment = list(
        enzyme_metabolite = list(require_shared_pathway = FALSE)))))
    expect_setequal(run_pairs(config = relax)$protein, c("P1", "P2"))

    stubs <- stub_annotation()
    stubs$extract_de_tables <- function(de_data, omics_type, harmonization_res = NULL) {
        std <- de_data$.std
        if (identical(omics_type, "proteomics")) {
            std$padj <- c(0.01, 0.9, 0.9)   # only P1 is a hit
            return(list(A_vs_B = std))
        }
        list(`A vs. B` = std)
    }
    expect_identical(only_pairs(run_pairs(stubs, config = relax))$protein, c("P1"))

    relax_all <- relax
    relax_all$modes$multiomics$enrichment$enzyme_metabolite$enzyme_hits_only <- FALSE
    # Paired proteins only: P3 is now listed too, as a row with no metabolite.
    expect_setequal(only_pairs(run_pairs(stubs, config = relax_all))$protein,
                    c("P1", "P2"))
})

test_that("a metabolite with no KEGG compound and a protein with no EC are dropped", {
    pairs <- run_pairs(list(map_metabolite_ids_to_kegg = function(...) data.frame(
        feature_id = "alpha", KEGG_CPD = "C00100", stringsAsFactors = FALSE)))
    expect_false(any(only_pairs(pairs)$metabolite == "beta"))
    # P3 is in the DE table and maps to a KEGG gene, but that gene has no EC.
    expect_false(any(pairs$protein == "P3"))
})

test_that("the hit flags carry the rule that produced them", {
    stubs <- stub_annotation()
    stubs$extract_de_tables <- function(de_data, omics_type, harmonization_res = NULL) {
        std <- de_data$.std
        std$padj <- rep(0.9, nrow(std))   # nothing clears adjustment
        if (identical(omics_type, "proteomics")) return(list(A_vs_B = std))
        list(`A vs. B` = std)
    }
    pairs <- run_pairs(stubs)
    expect_identical(unique(pairs$enzyme_hit_rule), "raw p")
    expect_true(all(pairs$enzyme_hit))

    expect_identical(de_hit_flags(data.frame(pvalue = c(0.001, 0.4),
                                             padj = c(0.01, 0.9)))$rule, "padj")
    expect_identical(de_hit_flags(data.frame(pvalue = c(0.001, 0.4),
                                             padj = c(0.2, 0.9)))$rule, "raw p")
})

test_that("KEGG being unreachable gives NULL and a message, not an error", {
    stubs <- stub_annotation()
    stubs$kegg_link_table <- function(target, source, cache_dir = NULL) NULL
    expect_message(pairs <- run_pairs(stubs, quiet = FALSE),
                   "KEGG enzyme annotation unavailable")
    expect_null(pairs)

    stubs2 <- stub_annotation()
    stubs2$resolve_kegg_org_code <- function(...) NULL
    expect_message(pairs2 <- run_pairs(stubs2, quiet = FALSE), "KEGG organism code")
    expect_null(pairs2)
})

test_that("a run with one layer, or no shared contrast, yields NULL", {
    local_stubs(stub_annotation())
    # Opted in, so the run reaches the layer and contrast checks rather than
    # stopping at the (default-off) switch.
    cfg <- list(global = list(organism = "Test organism"),
                modes = list(multiomics = list(enrichment = list(
                    enzyme_metabolite = list(enabled = TRUE)))))
    out <- withr::local_tempdir()

    expect_message(
        one_layer <- build_enzyme_metabolite_pairs(
            list(proteomics = prot_de()), list(), cfg, out),
        "need both a proteomics and a metabolomics layer")
    expect_null(one_layer)

    local_stubs(list(extract_de_tables = function(de_data, omics_type, ...) {
        if (identical(omics_type, "proteomics")) list(A_vs_B = de_data$.std)
        else list(X_vs_Y = de_data$.std)
    }))
    expect_message(
        no_contrast <- build_enzyme_metabolite_pairs(
            list(proteomics = prot_de(), metabolomics = metab_de()), list(), cfg, out),
        "No contrast is shared")
    expect_null(no_contrast)
})

test_that("the config switch turns the whole table off", {
    pairs <- run_pairs(config = list(modes = list(multiomics = list(
        enrichment = list(enzyme_metabolite = list(enabled = FALSE))))))
    expect_null(pairs)
})


# ---- equation parsing -------------------------------------------------------

test_that("an equation splits into substrates and products", {
    roles <- parse_kegg_equation("C00001 + 2 C00002 <=> C00003 + C00004")
    expect_identical(roles$compound, c("C00001", "C00002", "C00003", "C00004"))
    expect_identical(roles$role, c("substrate", "substrate", "product", "product"))
})

test_that("a compound on both sides gets both roles, and junk gives no rows", {
    both <- parse_kegg_equation("C00001 <=> C00001")
    expect_identical(both$role, c("substrate", "product"))

    expect_identical(nrow(parse_kegg_equation("no arrow here")), 0L)
    expect_identical(nrow(parse_kegg_equation(NA_character_)), 0L)
    expect_identical(nrow(parse_kegg_equation("G00001 <=> G00002")), 0L)
})

test_that("a KEGG flat-file response is parsed into reaction roles", {
    lines <- c(
        "ENTRY       R00001                      Reaction",
        "NAME        test reaction",
        "EQUATION    C00100 + C00002 <=> C00200",
        "///",
        "ENTRY       R00002                      Reaction",
        "EQUATION    C00300 <=>",
        "            C00400",
        "///")
    roles <- .parse_reaction_records(lines)

    expect_identical(unique(roles$reaction), c("R00001", "R00002"))
    expect_identical(roles$role[roles$compound == "C00200"], "product")
    # The wrapped continuation line belongs to the same equation.
    expect_identical(roles$role[roles$compound == "C00400"], "product")
})


# ---- writing ----------------------------------------------------------------

test_that("nothing to write leaves no file, and clears a previous one", {
    out <- withr::local_tempdir()
    path <- file.path(out, "enzyme_metabolite_pairs.tsv")
    writeLines("left over from an earlier run", path)

    expect_null(write_enzyme_metabolite_pairs(NULL, out))
    expect_false(file.exists(path))

    written <- write_enzyme_metabolite_pairs(
        data.frame(contrast = "A_vs_B", protein = "P1", stringsAsFactors = FALSE), out)
    expect_identical(written, path)
    expect_identical(nrow(read.delim(path, stringsAsFactors = FALSE)), 1L)
})

test_that("a bare KEGG gene id still joins the prefixed one in the link table", {
    # convert_entrez_to_kegg() returns "29740" where /link/ returns "tst:29740";
    # merging them as they come matched nothing and the table came back empty.
    prefixed <- run_pairs(list(
        convert_entrez_to_kegg = function(...) c("1" = "1", "2" = "2", "3" = "3"),
        kegg_link_table = function(target, source, cache_dir = NULL) {
            links <- list(
                enzyme = data.frame(from = c("tst:1", "tst:2"),
                                    to = c("1.1.1.1", "2.2.2.2"),
                                    stringsAsFactors = FALSE),
                reaction = data.frame(from = c("1.1.1.1", "1.1.1.1", "2.2.2.2"),
                                      to = c("R00001", "R00002", "R00003"),
                                      stringsAsFactors = FALSE),
                compound = data.frame(
                    from = c("R00001", "R00001", "R00002", "R00003", "R00003"),
                    to = c("C00100", "C00002", "C00100", "C00200", "C00300"),
                    stringsAsFactors = FALSE),
                pathway = data.frame(from = c("tst:1", "tst:2"),
                                     to = c("tst00010", "tst00020"),
                                     stringsAsFactors = FALSE))
            links[[target]]
        }))
    expect_true("alpha" %in% prefixed$metabolite)
})

test_that("kegg_gene_key strips only the organism prefix", {
    expect_identical(kegg_gene_key(c("rno:29740", "29740", "hsa:1"), "rno"),
                     c("29740", "29740", "hsa:1"))
    expect_identical(kegg_gene_key("29740", NULL), "29740")
})

test_that("every changed enzyme is listed, with a note when nothing pairs", {
    # P1 pairs with alpha; P2 carries an EC but its only compound shares no
    # pathway. P3 has no EC at all, so it is not an enzyme and is left out even
    # with the hits-only filter relaxed.
    all_rows <- run_pairs(config = list(modes = list(multiomics = list(
        enrichment = list(enzyme_metabolite = list(enzyme_hits_only = FALSE))))))
    expect_setequal(all_rows$protein, c("P1", "P2"))
    unpaired <- all_rows[is.na(all_rows$compound), ]
    expect_setequal(unpaired$protein, "P2")
    expect_match(unpaired$note[unpaired$protein == "P2"], "shared pathway")
    # An unpaired row carries the enzyme's own numbers and nothing invented.
    expect_true(all(is.na(unpaired$metabolite_log2fc)))
    expect_false(any(is.na(unpaired$enzyme_log2fc)))
    expect_false(any(is.na(all_rows$ec)))
})

test_that("the unpaired enzymes can be switched off", {
    only_pairs <- run_pairs(config = list(modes = list(multiomics = list(
        enrichment = list(enzyme_metabolite = list(list_unpaired_enzymes = FALSE))))))
    expect_true(all(!is.na(only_pairs$compound)))
    expect_setequal(only_pairs$protein, "P1")
})

test_that("pairs whose metabolite cleared FDR come first, unpaired enzymes last", {
    rows <- run_pairs(config = list(modes = list(multiomics = list(
        enrichment = list(enzyme_metabolite = list(require_shared_pathway = FALSE))))))
    # alpha is the only metabolite under FDR 0.05, so its pair leads.
    expect_identical(rows$metabolite[1], "alpha")
    # Any row without a compound sits below every row with one.
    if (any(is.na(rows$compound)) && any(!is.na(rows$compound))) {
        expect_gt(min(which(is.na(rows$compound))), max(which(!is.na(rows$compound))))
    }
})

test_that("a mapped metabolite carries KEGG's name for it where KEGG has one", {
    rows <- run_pairs()
    paired <- only_pairs(rows)
    expect_identical(paired$compound_name[paired$compound == "C00100"], "Alphanoate")
    # A compound KEGG does not name, and an enzyme with no pair, stay empty
    # rather than inventing a label.
    expect_true(all(is.na(rows$compound_name[is.na(rows$compound)])))

    unnamed <- run_pairs(list(kegg_compound_names = function(cache_dir = NULL) NULL))
    expect_true(all(is.na(only_pairs(unnamed)$compound_name)))
})

test_that("a failed KEGG request is retried before the table gives up", {
    calls <- 0L
    fake_read <- function(url, warn = FALSE) {
        calls <<- calls + 1L
        if (calls < 3L) stop("cannot open the connection")
        c("cpd:C00100\trn:R00001")
    }
    # The reader is called through readLines(); swap it in this function's env.
    target <- environment(.kegg_rest_lines)
    old <- get("readLines", envir = target, inherits = TRUE)
    assign("readLines", fake_read, envir = target)
    withr::defer(assign("readLines", old, envir = target))

    lines <- suppressMessages(.kegg_rest_lines("https://example.invalid", "test",
                                               attempts = 3L, pause = 0))
    expect_identical(calls, 3L)
    expect_length(lines, 1)

    calls <- 0L
    fake_read <- function(url, warn = FALSE) {
        calls <<- calls + 1L
        stop("cannot open the connection")
    }
    assign("readLines", fake_read, envir = target)
    expect_null(suppressMessages(.kegg_rest_lines("https://example.invalid", "test",
                                                  attempts = 2L, pause = 0)))
    expect_identical(calls, 2L)
})


# ---- review of #255 against main --------------------------------------------

relaxed_pathways <- list(modes = list(multiomics = list(enrichment = list(
    enzyme_metabolite = list(require_shared_pathway = FALSE)))))

test_that("two measured metabolites on one KEGG compound are two rows", {
    # alpha and alpha_neg are the same compound measured twice (say, in both
    # ionisation modes), with different statistics.
    stubs <- stub_annotation()
    stubs$extract_de_tables <- function(de_data, omics_type, harmonization_res = NULL) {
        std <- de_data$.std
        if (identical(omics_type, "proteomics")) return(list(A_vs_B = std))
        std <- rbind(std, data.frame(feature_id = "alpha_neg", log2fc = -0.5,
                                     pvalue = 0.3, padj = 0.6,
                                     stringsAsFactors = FALSE))
        list(`A vs. B` = std)
    }
    stubs$map_metabolite_ids_to_kegg <- function(...) data.frame(
        feature_id = c("alpha", "alpha_neg", "beta"),
        KEGG_CPD = c("C00100", "C00100", "C00200"), stringsAsFactors = FALSE)

    pairs <- only_pairs(run_pairs(stubs))
    p1 <- pairs[pairs$protein == "P1" & pairs$compound == "C00100", ]
    expect_setequal(p1$metabolite, c("alpha", "alpha_neg"))
    # Each row keeps its own metabolite's numbers, not the other one's.
    expect_equal(p1$metabolite_log2fc[p1$metabolite == "alpha"], 2)
    expect_equal(p1$metabolite_log2fc[p1$metabolite == "alpha_neg"], -0.5)
    # The reactions behind the pair are still joined within each row.
    expect_true(all(p1$reaction == "R00001;R00002"))
})

test_that("with the shared-pathway filter on, missing pathway annotation skips the table", {
    no_gene_paths <- stub_annotation(link_overrides = list(pathway = NULL))
    expect_message(pairs <- run_pairs(no_gene_paths, quiet = FALSE),
                   "pathway annotation unavailable.*skipping")
    expect_null(pairs)

    no_cpd_paths <- stub_annotation()
    no_cpd_paths$get_kegg_compound_pathways <- function(...) NULL
    expect_null(run_pairs(no_cpd_paths))
})

test_that("with the filter off, missing pathway annotation is said, not implied", {
    stubs <- stub_annotation(link_overrides = list(pathway = NULL))
    expect_message(pairs <- run_pairs(stubs, config = relaxed_pathways, quiet = FALSE),
                   "listed without pathways")
    paired <- only_pairs(pairs)
    expect_true(nrow(paired) > 0)
    expect_true(all(is.na(paired$pathway_id)))
    expect_true(all(paired$note == "KEGG pathway annotation unavailable"))
    # No unpaired enzyme is told it lacks a shared pathway: there was no filter.
    unpaired <- pairs[is.na(pairs$compound), , drop = FALSE]
    expect_false(any(grepl("shared pathway", unpaired$note)))
})

test_that("gene symbols are joined by feature id, not by row position", {
    stubs <- stub_annotation()
    stubs$protein_group_gene_symbols <- NULL
    local_stubs(stubs)
    # row_data is in a different order from the expression rows; a positional
    # pairing would give P1 the symbol of P3.
    harm <- list(inputs = list(proteomics = list(
        expr_work = matrix(0, 3, 1, dimnames = list(c("P1", "P2", "P3"), "S1")),
        row_data = data.frame(Protein.Group = c("P3", "P1", "P2"),
                              Genes = c("Ccc", "Aaa", "Bbb"),
                              stringsAsFactors = FALSE))))
    cfg <- list(global = list(organism = "Test organism"),
                modes = list(multiomics = list(enrichment = list(
                    enzyme_metabolite = list(enabled = TRUE)))))
    pairs <- suppressMessages(build_enzyme_metabolite_pairs(
        de_results = list(proteomics = prot_de(), metabolomics = metab_de()),
        harmonization_res = harm, config = cfg, out_dir = withr::local_tempdir()))
    expect_identical(unique(pairs$gene_symbol[pairs$protein == "P1"]), "Aaa")
})

test_that("the legend says when 'changed' rests on raw p, and where", {
    cfg <- enzyme_metabolite_config(list())
    fdr_only <- data.frame(contrast = c("A_vs_B", "C_vs_D"),
                           enzyme_hit_rule = "padj", metabolite_hit_rule = "padj",
                           note = NA, stringsAsFactors = FALSE)
    leg <- describe_enzyme_metabolite_table(fdr_only, cfg)
    expect_match(leg, "A changed enzyme means FDR < 0.05.", fixed = TRUE)
    expect_false(grepl("raw p", leg))

    mixed <- fdr_only
    mixed$enzyme_hit_rule[2] <- "raw p"
    leg <- describe_enzyme_metabolite_table(mixed, cfg)
    expect_match(leg, "except in C_vs_D, where no enzyme cleared FDR and it means raw p < 0.05, uncorrected (enzyme_hit_rule)",
                 fixed = TRUE)
    expect_match(leg, "A changed metabolite means FDR < 0.05.", fixed = TRUE)
})

test_that("the legend follows the settings the table was built with", {
    tab <- data.frame(contrast = "A_vs_B", enzyme_hit_rule = "padj",
                      metabolite_hit_rule = "padj", note = NA,
                      stringsAsFactors = FALSE)
    on <- describe_enzyme_metabolite_table(tab, enzyme_metabolite_config(list()))
    expect_match(on, "Every enzyme that changed", fixed = TRUE)
    expect_match(on, "a pathway both are in", fixed = TRUE)
    expect_match(on, "listed too, and `note` says why", fixed = TRUE)

    off <- enzyme_metabolite_config(list(modes = list(multiomics = list(enrichment = list(
        enzyme_metabolite = list(enzyme_hits_only = FALSE, require_shared_pathway = FALSE,
                                 list_unpaired_enzymes = FALSE))))))
    leg <- describe_enzyme_metabolite_table(tab, off)
    expect_match(leg, "Every measured enzyme, changed or not", fixed = TRUE)
    expect_match(leg, "no shared-pathway requirement", fixed = TRUE)
    expect_match(leg, "Only enzymes with at least one pair are listed", fixed = TRUE)
    expect_false(grepl("a pathway both are in", leg, fixed = TRUE))

    no_paths <- tab
    no_paths$note <- "KEGG pathway annotation unavailable"
    expect_match(describe_enzyme_metabolite_table(no_paths, off),
                 "could not be fetched for this run", fixed = TRUE)
})

test_that("the report shows both hit rules and builds its legend from the table", {
    rmd <- paste(readLines(testthat::test_path(
        "..", "..", "R", "domain", "multiomics", "report_template_multiomics.Rmd")),
        collapse = "\n")
    cols <- regmatches(rmd, regexpr("em_cols <- c\\([^)]*\\)", rmd))
    expect_length(cols, 1)
    expect_match(cols, '"enzyme_hit_rule"', fixed = TRUE)
    expect_match(cols, '"metabolite_hit_rule"', fixed = TRUE)
    expect_match(rmd, "describe_enzyme_metabolite_table(em, enzyme_metabolite_config(cfg))",
                 fixed = TRUE)
})

test_that("changed proteins without an EC are counted once, whatever the config", {
    tabs <- list(
        A = data.frame(feature_id = c("P1", "P3", "P4", "P5"), pvalue = 0.001,
                       padj = c(0.01, 0.01, 0.9, 0.01), stringsAsFactors = FALSE),
        B = data.frame(feature_id = c("P3", "P4"), pvalue = 0.001,
                       padj = c(0.01, 0.01), stringsAsFactors = FALSE))
    ec <- data.frame(feature_id = "P1", ec = "1.1.1.1", stringsAsFactors = FALSE)
    # P3 changed in both contrasts and counts once; P4 changed only in B; P5
    # never reached a KEGG gene, which is a mapping gap, not a missing EC.
    out <- .changed_not_listed(tabs, ec, mapped = c("P1", "P3", "P4"), cutoff = 0.05)
    expect_identical(out$no_ec, c("P3", "P4"))
    expect_identical(out$unmapped, "P5")

    # Through the builder: P3 is a hit with no EC. The count is the same whether
    # unchanged enzymes and unpaired enzymes are listed or not.
    stubs <- stub_annotation()
    stubs$extract_de_tables <- function(de_data, omics_type, harmonization_res = NULL) {
        std <- de_data$.std
        if (identical(omics_type, "proteomics")) {
            std$padj <- c(0.01, 0.9, 0.01)
            return(list(A_vs_B = std))
        }
        list(`A vs. B` = std)
    }
    for (hits_only in c(TRUE, FALSE)) for (unpaired in c(TRUE, FALSE)) {
        cfg <- list(modes = list(multiomics = list(enrichment = list(
            enzyme_metabolite = list(enzyme_hits_only = hits_only,
                                     list_unpaired_enzymes = unpaired)))))
        expect_message(run_pairs(stubs, config = cfg, quiet = FALSE),
                       "Not listed: 1 changed protein\\(s\\) with no EC number, and 0 whose id")
    }

    # The same P3, now failing its KEGG gene conversion: counted as unmapped.
    unmapped <- stubs
    unmapped$convert_entrez_to_kegg <- function(...) c("1" = "tst:1", "2" = "tst:2")
    expect_message(run_pairs(unmapped, quiet = FALSE),
                   "Not listed: 0 changed protein\\(s\\) with no EC number, and 1 whose id")
})


# ---- Codex review of b979590 -------------------------------------------------

test_that("the shared-pathway filter uses each protein's own pathways", {
    # P4 carries the same EC as P1 but its gene sits only in tst00099, a
    # pathway C00100 is not in. P1's pathway must not qualify P4's pair.
    # Replaced whole rather than through link_overrides: modifyList() merges
    # data frames column by column and cannot grow them.
    stubs <- stub_annotation()
    base_links <- stubs$kegg_link_table
    stubs$kegg_link_table <- function(target, source, cache_dir = NULL) {
        switch(target,
               enzyme = data.frame(from = c("tst:1", "tst:2", "tst:4"),
                                   to = c("1.1.1.1", "2.2.2.2", "1.1.1.1"),
                                   stringsAsFactors = FALSE),
               pathway = data.frame(from = c("tst:1", "tst:2", "tst:4"),
                                    to = c("tst00010", "tst00020", "tst00099"),
                                    stringsAsFactors = FALSE),
               base_links(target, source, cache_dir))
    }
    stubs$extract_de_tables <- function(de_data, omics_type, harmonization_res = NULL) {
        std <- de_data$.std
        if (identical(omics_type, "proteomics")) {
            std <- rbind(std, data.frame(feature_id = "P4", log2fc = 1, pvalue = 0.001,
                                         padj = 0.01, stringsAsFactors = FALSE))
            return(list(A_vs_B = std))
        }
        list(`A vs. B` = std)
    }
    stubs$map_feature_ids_to_entrez <- function(...) data.frame(
        feature_id = c("P1", "P2", "P3", "P4"), ENTREZID = c("1", "2", "3", "4"),
        stringsAsFactors = FALSE)
    stubs$convert_entrez_to_kegg <- function(...) c("1" = "tst:1", "2" = "tst:2",
                                                    "3" = "tst:3", "4" = "tst:4")

    pairs <- run_pairs(stubs)
    expect_true("P1" %in% only_pairs(pairs)$protein)
    expect_false("P4" %in% only_pairs(pairs)$protein)
    expect_identical(pairs$note[pairs$protein == "P4"],
                     "no measured metabolite in a shared pathway")
})

test_that("an enzyme whose only measured compounds are currency ones says so", {
    # P1 links C00100 (not measured here) and ATP (measured, and a currency
    # compound); nothing is left to pair once ATP is dropped.
    stubs <- stub_annotation()
    stubs$extract_de_tables <- function(de_data, omics_type, harmonization_res = NULL) {
        std <- de_data$.std
        if (identical(omics_type, "proteomics")) return(list(A_vs_B = std))
        std$feature_id <- c("atp", "beta")
        list(`A vs. B` = std)
    }
    stubs$map_metabolite_ids_to_kegg <- function(...) data.frame(
        feature_id = c("atp", "beta"), KEGG_CPD = c("C00002", "C00200"),
        stringsAsFactors = FALSE)

    pairs <- run_pairs(stubs)
    expect_identical(pairs$note[pairs$protein == "P1"],
                     "only currency metabolites measured, and those are dropped")
    # Keeping currency compounds pairs it with ATP, which shares its pathway.
    keep <- list(modes = list(multiomics = list(enrichment = list(
        enzyme_metabolite = list(drop_currency_metabolites = FALSE)))))
    expect_true("atp" %in% only_pairs(run_pairs(stubs, config = keep))$metabolite)

    # And an enzyme with nothing measured at all is told apart from both.
    cfg <- enzyme_metabolite_config(list())
    expect_identical(.unpaired_note(FALSE, FALSE, cfg), "no measured metabolite")
    expect_identical(.unpaired_note(TRUE, TRUE, cfg),
                     "no measured metabolite in a shared pathway")
})

test_that("the report reads the TSV's blanks back as missing values", {
    dir <- withr::local_tempdir()
    tab <- data.frame(contrast = c("A_vs_B", "A_vs_B"), protein = c("P1", "P2"),
                      compound = c("C00100", NA), note = c(NA, "no measured metabolite"),
                      stringsAsFactors = FALSE)
    path <- write_enzyme_metabolite_pairs(tab, dir)

    rmd <- paste(readLines(testthat::test_path(
        "..", "..", "R", "domain", "multiomics", "report_template_multiomics.Rmd")),
        collapse = "\n")
    call <- regmatches(rmd, regexpr("em <- read\\.delim\\([^)]*\\)", rmd))
    expect_length(call, 1)
    expect_match(call, 'na.strings = ""', fixed = TRUE)

    # Read the way the report reads it: one pair and one unpaired enzyme.
    em <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE,
                     na.strings = "")
    expect_identical(sum(!is.na(em$compound)), 1L)
    expect_identical(sum(is.na(em$compound)), 1L)
})


test_that("the table is opt-in, and its flags must be real booleans", {
    expect_false(enzyme_metabolite_config(list())$enabled)
    expect_true(enzyme_metabolite_config(list())$require_shared_pathway)
    # A string never switches anything on.
    yes <- enzyme_metabolite_config(list(modes = list(multiomics = list(enrichment = list(
        enzyme_metabolite = list(enabled = "yes"))))))
    expect_false(yes$enabled)

    validate <- function(em) suppressMessages(suppressWarnings(
        validate_multiomics_config(list(integration = list(methods = "SNF"),
                                        enrichment = list(enzyme_metabolite = em)))))
    v <- validate(list())
    expect_false(v$enrichment$enzyme_metabolite$enabled)
    expect_true(v$enrichment$enzyme_metabolite$drop_currency_metabolites)
    expect_error(validate(list(enabled = "no")), "enabled must be true or false")
    expect_error(validate(list(require_shared_pathway = NA)), "must be true or false")
    expect_error(validate(list(list_unpaired_enzymes = c(TRUE, FALSE))),
                 "must be true or false")
    expect_true(validate(list(enabled = TRUE))$enrichment$enzyme_metabolite$enabled)
})


# ---- Codex review of 9cf232a -------------------------------------------------

test_that("an unpaired enzyme has no role, and a pair without an equation is unknown", {
    pairs <- run_pairs()
    # P2 has nothing to pair with: no reaction, no compound, so no role.
    unpaired <- pairs[is.na(pairs$compound), , drop = FALSE]
    expect_true("P2" %in% unpaired$protein)
    expect_true(all(is.na(unpaired$role)))
    expect_identical(only_pairs(pairs)$role, "product;substrate")

    # A real pair whose equation could not be fetched is still "unknown".
    no_eq <- stub_annotation()
    no_eq$fetch_reaction_roles <- function(...) NULL
    fetched <- run_pairs(no_eq)
    expect_identical(only_pairs(fetched)$role, "unknown")
    expect_true(all(is.na(fetched$role[is.na(fetched$compound)])))
})


# ---- Copilot review of 4fc7ca1 -----------------------------------------------

test_that("the unlisted-protein count is logged even when no pair survives", {
    # Only P2 changes, and its one measured compound shares no pathway with it,
    # so no pair survives; P3 changes too, with no EC.
    stubs <- stub_annotation()
    stubs$extract_de_tables <- function(de_data, omics_type, harmonization_res = NULL) {
        std <- de_data$.std
        if (identical(omics_type, "proteomics")) {
            std$padj <- c(0.9, 0.01, 0.01)
            return(list(A_vs_B = std))
        }
        list(`A vs. B` = std)
    }
    no_list <- list(modes = list(multiomics = list(enrichment = list(
        enzyme_metabolite = list(list_unpaired_enzymes = FALSE)))))
    msgs <- character(0)
    res <- withCallingHandlers(run_pairs(stubs, config = no_list, quiet = FALSE),
                               message = function(m) {
                                   msgs <<- c(msgs, conditionMessage(m))
                                   invokeRestart("muffleMessage")
                               })
    expect_null(res)
    expect_true(any(grepl("No enzyme-metabolite pair survived", msgs)))
    expect_true(any(grepl("Not listed: 1 changed protein\\(s\\) with no EC number", msgs)))

    # And on the pathway-skip return.
    skip <- stub_annotation(link_overrides = list(pathway = NULL))
    expect_message(run_pairs(skip, quiet = FALSE), "Not listed: 0 changed protein")
})

test_that("an unknown enzyme_metabolite key is rejected, not ignored", {
    validate <- function(em) suppressMessages(suppressWarnings(
        validate_multiomics_config(list(integration = list(methods = "SNF"),
                                        enrichment = list(enzyme_metabolite = em)))))
    expect_error(validate(list(enabled = TRUE, require_shared_pathways = FALSE)),
                 "Unknown key\\(s\\).*require_shared_pathways")
})
