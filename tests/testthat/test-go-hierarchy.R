# tests/testthat/test-go-hierarchy.R
#
# GO hierarchy figure of the layers' existing GO fgsea results. The GO graph
# here is hand-made, so nothing depends on GO.db:
#
#   GO:0000001                      (root; its GO.db parent "all" is dropped)
#    |- is_a -    GO:0000002
#    |- part_of - GO:0000003
#   GO:0000004  is_a GO:0000002, is_a GO:0000003
#   GO:0000005  is_a GO:0000002
#   GO:0000006  regulates GO:0000004   (not an inherited relation: dropped)

toy_parents_raw <- function() {
    data.frame(
        child    = c("GO:0000001", "GO:0000002", "GO:0000003", "GO:0000004",
                     "GO:0000004", "GO:0000005", "GO:0000006"),
        parent   = c("all", "GO:0000001", "GO:0000001", "GO:0000002",
                     "GO:0000003", "GO:0000002", "GO:0000004"),
        relation = c("isa", "isa", "part_of", "isa", "isa", "isa", "regulates"),
        stringsAsFactors = FALSE)
}

toy_terms <- function() {
    data.frame(go_id = sprintf("GO:%07d", 1:6),
               term = paste("term", 1:6),
               ontology = "BP", stringsAsFactors = FALSE)
}

# Columns exactly as run_pathway_analysis() writes an fgsea table.
fgsea_table <- function(ids, nes, padj, contrast = "Perfused_vs_Control",
                        database = "GO_BP") {
    data.frame(pathway = ids, pathway_name = paste("name", ids),
               pval = padj / 2, padj = padj, log2err = 0.5, ES = nes / 2,
               NES = nes, size = 20L, leadingEdge = "G1,G2",
               contrast = contrast, database = database, method = "fgsea",
               stringsAsFactors = FALSE)
}

gh_config <- function(...) {
    list(modes = list(multiomics = list(go_hierarchy = utils::modifyList(
        list(enabled = TRUE, ontologies = list("GO_BP"), fdr = 0.05,
             max_seed_terms = 30L, max_nodes = 120L), list(...)))))
}

# ---- GO graph -----------------------------------------------------------------

test_that("relations are normalised and only is_a / part_of are kept", {
    p <- normalize_go_relations(toy_parents_raw())
    expect_setequal(unique(p$relation), c("is_a", "part_of"))
    expect_false("all" %in% p$parent)
    expect_false("GO:0000006" %in% p$child)
})

test_that("the closure holds the seeds and every ancestor, nothing else", {
    p <- normalize_go_relations(toy_parents_raw())
    expect_setequal(go_dag_closure("GO:0000004", p),
                    c("GO:0000004", "GO:0000002", "GO:0000003", "GO:0000001"))
    expect_setequal(go_dag_closure("GO:0000005", p),
                    c("GO:0000005", "GO:0000002", "GO:0000001"))
    expect_identical(go_dag_closure("GO:0000001", p), "GO:0000001")
})

test_that("edges are the parent links with both ends in the node set", {
    p <- normalize_go_relations(toy_parents_raw())
    dag <- build_go_dag(data.frame(go_id = "GO:0000004", padj = 0.01), p,
                        term_info = toy_terms())
    expect_true(all(dag$edges$child %in% dag$nodes$go_id))
    expect_true(all(dag$edges$parent %in% dag$nodes$go_id))
    expect_equal(nrow(dag$edges), 4)
    expect_false("GO:0000005" %in% c(dag$edges$child, dag$edges$parent))
    expect_identical(dag$edges$relation[dag$edges$child == "GO:0000003"], "part_of")
    expect_identical(dag$nodes$is_seed, dag$nodes$go_id == "GO:0000004")
    expect_true(all(is.na(dag$nodes$padj[!dag$nodes$is_seed])))
    expect_identical(dag$nodes$term[dag$nodes$go_id == "GO:0000001"], "term 1")
})

test_that("the node limit drops the last seeds first and says how many", {
    p <- normalize_go_relations(toy_parents_raw())
    seeds <- data.frame(go_id = c("GO:0000004", "GO:0000005"), padj = c(0.001, 0.01))
    dag <- build_go_dag(seeds, p, term_info = toy_terms(), max_nodes = 4)
    expect_equal(dag$n_seeds_requested, 2)
    expect_equal(dag$n_seeds_shown, 1)
    expect_identical(dag$nodes$go_id[dag$nodes$is_seed], "GO:0000004")
    expect_lte(nrow(dag$nodes), 4)

    full <- build_go_dag(seeds, p, term_info = toy_terms(), max_nodes = 120)
    expect_equal(full$n_seeds_shown, 2)
    expect_equal(nrow(full$nodes), 5)
})

test_that("adjusted p bins stop at the cutoff and order by strength", {
    b <- go_padj_bins(c(1e-10, 0.03, NA), fdr = 0.05)
    expect_true(is.ordered(b))
    expect_true(b[1] > b[2])
    expect_true(is.na(b[3]))
    expect_length(levels(go_padj_bins(0.0001, fdr = 0.001)), 4)
})

# ---- reading the existing fgsea results ----------------------------------------

test_that("GO fgsea tables are read in both result shapes, other tables skipped", {
    contrast_res <- list(
        GO_BP_fgsea  = fgsea_table(c("GO:0000004", "GO:0000005"), c(2, -1.5), c(0.01, 0.02)),
        KEGG_fgsea   = fgsea_table("hsa04110", 2, 0.01, database = "KEGG"),
        GO_BP_ora_up = transform(fgsea_table("GO:0000002", 2, 0.01), method = "ora"))
    prot_shape <- list(Perfused_vs_Control = contrast_res)
    rna_shape  <- list(annotation = NULL,
                       pathway_results = list(Perfused_vs_Control = contrast_res),
                       plot_files = list())
    for (res in list(prot_shape, rna_shape)) {
        tbl <- collect_go_fgsea_tables(res)
        expect_setequal(tbl$go_id, c("GO:0000004", "GO:0000005"))
        expect_true(all(tbl$ontology == "BP"))
        expect_true(all(tbl$contrast == "Perfused_vs_Control"))
        expect_identical(tbl$term[tbl$go_id == "GO:0000004"], "name GO:0000004")
        expect_equal(tbl$NES[tbl$go_id == "GO:0000005"], -1.5)
    }
    expect_equal(nrow(collect_go_fgsea_tables(NULL)), 0)
})

test_that("the mixed GO collection is split by each term's ontology", {
    res <- list(A_vs_B = list(GO_fgsea = fgsea_table(
        c("GO:0000004", "GO:0000099"), c(2, 2), c(0.01, 0.01), database = "GO")))
    tbl <- collect_go_fgsea_tables(res, ontology_of = function(ids)
        ifelse(ids == "GO:0000004", "BP", "MF"))
    expect_identical(tbl$ontology[order(tbl$go_id)], c("BP", "MF"))
})

test_that("seeds pass padj < fdr, split by NES sign, ordered and capped", {
    tbl <- data.frame(
        go_id = sprintf("GO:%07d", 1:6),
        term  = paste("t", 1:6),
        NES   = c(2, 1.8, 1.5, -2, 0, 2),
        padj  = c(0.02, 0.001, 0.01, 0.03, 0.001, 0.05),
        pval  = 0.001, size = 15L, stringsAsFactors = FALSE)
    s <- select_go_dag_seeds(tbl, fdr = 0.05, max_seed_terms = 2)
    expect_identical(s$positive$seeds$go_id, c("GO:0000002", "GO:0000003"))
    expect_equal(s$positive$n_significant, 3)       # padj = 0.05 is not < 0.05
    expect_identical(s$negative$seeds$go_id, "GO:0000004")
    expect_false("GO:0000005" %in% c(s$positive$seeds$go_id, s$negative$seeds$go_id))
})

test_that("directions are named after the contrast arms when they can be read", {
    expect_identical(go_direction_label("Perfused_vs_Control", "positive"), "higher in Perfused")
    expect_identical(go_direction_label("Perfused_vs_Control", "negative"), "higher in Control")
    expect_identical(go_direction_label("Treatedvs.Control", "positive"), "higher in Treated")
    expect_identical(go_direction_label("ConditionEffect", "negative"), "negative NES")
})

# ---- building and writing -------------------------------------------------------

toy_results <- function() {
    list(proteomics = list(Perfused_vs_Control = list(
        GO_BP_fgsea = fgsea_table(c("GO:0000004", "GO:0000005"), c(2, 1.2),
                                  c(0.001, 0.02)))))
}

test_that("disabled returns NULL and builds nothing", {
    expect_null(build_go_hierarchy_dags(toy_results(), gh_config(enabled = FALSE)))
})

test_that("one unit per contrast, ontology and direction, empty ones noted", {
    skip_if_not_installed("ggraph")
    skip_if_not_installed("igraph")
    built <- build_go_hierarchy_dags(
        toy_results(), gh_config(),
        parent_tables = list(BP = normalize_go_relations(toy_parents_raw())),
        term_info = toy_terms())
    dirs <- vapply(built$units, `[[`, "", "direction")
    expect_setequal(dirs, c("positive", "negative"))
    pos <- built$units[[which(dirs == "positive")]]
    neg <- built$units[[which(dirs == "negative")]]
    expect_equal(pos$n_significant, 2)
    expect_s3_class(pos$plot, "ggplot")
    expect_identical(pos$direction_label, "higher in Perfused")
    expect_null(neg$dag)
    expect_match(neg$note, "No GO BP term")
})

test_that("a layer with no GO fgsea result gets one note, not an error", {
    res <- list(transcriptomics = list(pathway_results = list(A_vs_B = list(
        KEGG_fgsea = fgsea_table("hsa04110", 2, 0.01, database = "KEGG")))))
    built <- build_go_hierarchy_dags(res, gh_config())
    expect_length(built$units, 1)
    expect_true(is.na(built$units[[1]]$contrast))
    expect_match(built$units[[1]]$note, "No fgsea result on GO")
})

test_that("outputs, index and the portal edge table are written and read back", {
    skip_if_not_installed("ggraph")
    skip_if_not_installed("igraph")
    run_dir <- withr::local_tempdir()
    out_dir <- file.path(run_dir, "multiomics", "go_hierarchy")
    built <- build_go_hierarchy_dags(
        toy_results(), gh_config(),
        parent_tables = list(BP = normalize_go_relations(toy_parents_raw())),
        term_info = toy_terms())
    files <- write_go_hierarchy_outputs(built, out_dir, run_dir, term_info = toy_terms())
    expect_true(all(file.exists(files)))

    idx <- read_go_hierarchy_index(out_dir)
    expect_equal(nrow(idx), 2)
    drawn <- idx[!is.na(idx$png), ]
    expect_equal(nrow(drawn), 1)
    expect_true(file.exists(drawn$png))
    expect_true(file.exists(drawn$pdf))
    expect_match(drawn$png, "proteomics/Perfused_vs_Control/BP_positive.png", fixed = TRUE)

    tsv <- read.delim(file.path(run_dir, "proteomics", "Datasets", "go_dag_edges.tsv"),
                      stringsAsFactors = FALSE)
    expect_identical(colnames(tsv),
                     c("child", "parent", "relation", "parent_name", "parent_ontology"))
    expect_false(any(duplicated(tsv[, c("child", "parent", "relation")])))

    # A later run with the figure switched off leaves nothing behind.
    write_go_hierarchy_outputs(NULL, out_dir, run_dir)
    expect_false(dir.exists(out_dir))
    expect_false(file.exists(file.path(run_dir, "proteomics", "Datasets", "go_dag_edges.tsv")))
})

# ---- config -------------------------------------------------------------------

test_that("the resolver defaults to off, GO BP, fdr 0.05, 30 seeds", {
    gh <- go_hierarchy_config(list(modes = list(multiomics = list())))
    expect_false(gh$enabled)
    expect_identical(gh$ontologies, "BP")
    expect_equal(gh$fdr, 0.05)
    expect_equal(gh$max_seed_terms, 30L)
})

test_that("the validator fills defaults and rejects bad values", {
    base <- list(integration = list(methods = "SNF", snf = list(K = 3, n_clusters = 2)),
                 feature_selection = list(method = "variance", top_n = 500))
    v <- validate_multiomics_config(base)
    expect_false(v$go_hierarchy$enabled)
    expect_equal(unlist(v$go_hierarchy$ontologies), "GO_BP")
    expect_equal(v$go_hierarchy$fdr, 0.05)
    expect_equal(v$go_hierarchy$max_seed_terms, 30L)

    bad <- function(...) { b <- base; b$go_hierarchy <- list(...); b }
    expect_error(validate_multiomics_config(bad(ontologies = "BP")), "GO_BP")
    expect_error(validate_multiomics_config(bad(fdr = 2)), "fdr")
    expect_error(validate_multiomics_config(bad(enabled = "yes")), "true or false")
    expect_error(validate_multiomics_config(bad(max_seed_terms = 0)), "max_seed_terms")
    expect_error(validate_multiomics_config(bad(seeds = 10)), "Unknown key")
})

test_that("the shipped template leaves the figure off", {
    cfg <- yaml::read_yaml(file.path(root_dir, "config", "templates",
                                     "multiomics_config.yaml"))
    gh <- cfg$modes$multiomics$go_hierarchy
    expect_false(isTRUE(gh$enabled))
    expect_identical(unlist(gh$ontologies), "GO_BP")
})

# ---- report gating ------------------------------------------------------------

test_that("the report section is gated on the resolver and reads the index", {
    f <- file.path(root_dir, "R", "domain", "multiomics", "report_template_multiomics.Rmd")
    skip_if_not(file.exists(f), "multiomics report template not found")
    src <- paste(readLines(f, warn = FALSE), collapse = "\n")
    expect_true(grepl('`r if (show_go_hierarchy) "# GO hierarchy of GSEA results {.tabset}"`',
                      src, fixed = TRUE))
    expect_true(grepl("show_go_hierarchy <- go_hier_cfg$enabled", src, fixed = TRUE))
    expect_true(grepl("go_hier_cfg <- go_hierarchy_config(cfg)", src, fixed = TRUE))
    expect_true(grepl("read_go_hierarchy_index(go_hier_dir)", src, fixed = TRUE))
    for (chunk in c("go-hierarchy-methods", "go-hierarchy-body")) {
        line <- regmatches(src, regexpr(sprintf("```\\{r %s,[^}]*\\}", chunk), src))
        expect_length(line, 1)
        expect_true(grepl("eval=show_go_hierarchy", line, fixed = TRUE), info = chunk)
    }
})

test_that("the GO hierarchy legend stays within 100 words in its worst case", {
    n_words <- function(x) lengths(regmatches(x, gregexpr("[^[:space:]]+", x)))
    # Longest note: both caps applied; a three-word arm name.
    note <- go_hierarchy_note(n_significant = 500, n_requested = 30, n_shown = 12,
                              max_nodes = 120)
    txt <- go_hierarchy_legend("BP", 0.05, "positive",
                               "higher in long arm name", note)
    expect_lte(n_words(txt), 100)
    tpl <- paste(readLines(file.path(root_dir, "R/domain/multiomics/report_template_multiomics.Rmd")),
                 collapse = "\n")
    expect_true(grepl("figure_legend(go_legend)", tpl, fixed = TRUE))
})
