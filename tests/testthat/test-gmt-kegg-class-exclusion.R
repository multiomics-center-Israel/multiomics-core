# Excluded KEGG classes on the GMT fallback path.
#
# The class exclusion is a reporting decision a project makes once, in its
# config. #205 wired it into the KEGG branch of the gene-based multi-ORA and
# left the GMT branch unfiltered. That branch is not an exotic corner: it is
# what runs when an organism has no KEGG code or no OrgDb, and a GMT assembled
# for such an organism routinely carries KEGG sets alongside GO, Pfam and
# InterPro ones. So the same config was honoured on a model organism and
# ignored on exactly the runs the exclusion exists for.
#
# Two properties this holds to:
#
#   The exclusion is applied to a finished table. Which branch the wrapper
#   takes -- adjusted hits, or the raw-p fallback -- is decided on the
#   unfiltered result, so a project's reporting choice can never enlarge its
#   own result set by emptying the adjusted branch and falling through to the
#   wider one.
#
#   Only rows is_kegg_pathway_accession() recognises have a class at all. A GO
#   term, a Pfam accession or a custom gene-set name in the same GMT is not a
#   KEGG pathway of unknown class -- it is not a KEGG pathway -- and comes back
#   untouched.
#
# Everything here is synthetic and offline. clusterProfiler::enricher runs for
# real against a hand-built TERM2GENE, which reaches no network; the BRITE
# classification is stubbed, so both the excluded classes and the fail-open
# path are asserted on every machine rather than only on one with network.

# enricher is namespaced and cannot be rebound this way, but the classification
# fetch is a sourced global -- and it is the one that would otherwise reach out.
# tar_source() puts these in the global environment rather than a namespace, so
# with_mocked_bindings() has nothing to rebind; assign into the function's own
# environment and restore on exit.
local_gmt_stubs <- function(stubs, env = parent.frame()) {
    target <- environment(run_multi_ora_enricher)
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

# A miniature br08901 in the shape kegg_pathway_categories() returns: two
# classes, so one can be excluded while the other stays.
brite_classes <- function() {
    data.frame(
        pathway_id   = c("00010", "00020", "04260", "04620"),
        category     = c("Metabolism", "Metabolism",
                         "Organismal Systems", "Organismal Systems"),
        subcategory  = c("Carbohydrate metabolism", "Carbohydrate metabolism",
                         "Circulatory system", "Immune system"),
        pathway_name = c("Glycolysis / Gluconeogenesis", "Citrate cycle",
                         "Cardiac muscle contraction",
                         "Toll-like receptor signaling pathway"),
        stringsAsFactors = FALSE
    )
}

use_fixture_classification <- function(env = parent.frame()) {
    local_gmt_stubs(list(kegg_pathway_categories = function(...) brite_classes()),
                    env = env)
}


# ---- the fixture ------------------------------------------------------------
#
# A GMT of the kind this path exists for: two KEGG maps under the run's own
# organism prefix, and three sets from collections KEGG knows nothing about.
#
# Membership is chosen for the p-values it produces, not for biology. Five
# tested sets of twenty genes partition g001-g100; the significant genes are
# g001-g042, split so the two KEGG maps take twelve each and the three
# non-KEGG sets six each. Six hundred further genes sit in filler sets holding
# no significant gene at all -- enricher never tests a set the query does not
# touch, so they stay out of the BH family while still making the background
# large enough for both overlap levels to be significant, the KEGG ones by
# several orders of magnitude more. That gap is what lets the tests below say
# which branch ran, and it holds whichever way a clusterProfiler version
# reconciles the universe with the annotated pool, because here they are equal.
KEGG_SETS     <- c("hsa00010", "hsa04260")
NON_KEGG_SETS <- c("GO:0006915", "PF00069", "custom_set_A")

gmt_sig <- function() sprintf("g%03d", 1:42)

gmt_term2gene <- function() {
    set <- function(term, idx) {
        data.frame(term = term, gene = sprintf("g%03d", idx),
                   stringsAsFactors = FALSE)
    }
    tested <- rbind(
        set("hsa00010",     c(1:12,  43:50)),   # Metabolism -- kept below
        set("hsa04260",     c(13:24, 51:58)),   # Organismal Systems -- excluded
        set("GO:0006915",   c(25:30, 59:72)),
        set("PF00069",      c(31:36, 73:86)),
        set("custom_set_A", c(37:42, 87:100))
    )
    filler <- do.call(rbind, lapply(1:30, function(i) {
        set(sprintf("filler_%02d", i), 100 + (i - 1) * 20 + 1:20)
    }))
    rbind(tested, filler)
}

gmt_universe <- function() sort(unique(gmt_term2gene()$gene))

gmt_term2name <- function() {
    data.frame(
        term = c(KEGG_SETS, NON_KEGG_SETS),
        name = c("Glycolysis / Gluconeogenesis", "Cardiac muscle contraction",
                 "Apoptosis", "Protein kinase domain", "Curated set A"),
        stringsAsFactors = FALSE
    )
}

# One entry point, so every test drives identical inputs and differs only in
# what it asks to be excluded.
run_gmt_ora <- function(exclude_classes = NULL, pval_cutoff = 0.1,
                        kegg_org = "hsa") {
    suppressMessages(run_multi_ora_enricher(
        gmt_sig(), gmt_universe(), gmt_term2gene(), gmt_term2name(),
        label = "pooled", pval_cutoff = pval_cutoff,
        exclude_classes = exclude_classes, kegg_org = kegg_org))
}


# ---- what the exclusion removes, and what it leaves alone -------------------

test_that("an excluded KEGG class disappears from the GMT table", {
    skip_if_not_installed("clusterProfiler")
    use_fixture_classification()

    unfiltered <- run_gmt_ora()
    # The fixture's premise, asserted rather than assumed: every set reaches
    # the finished table, so the exclusion has something to remove and
    # something to leave.
    expect_setequal(unfiltered$ID, c(KEGG_SETS, NON_KEGG_SETS))

    filtered <- run_gmt_ora(exclude_classes = "Organismal Systems")

    expect_false("hsa04260" %in% filtered$ID)
    expect_true("hsa00010" %in% filtered$ID)
})

test_that("a subcategory can be excluded on this path too", {
    skip_if_not_installed("clusterProfiler")
    use_fixture_classification()

    filtered <- run_gmt_ora(exclude_classes = "Circulatory system")

    expect_false("hsa04260" %in% filtered$ID)
    expect_true("hsa00010" %in% filtered$ID)
})

test_that("a retained KEGG row carries exactly the statistics it arrived with", {
    skip_if_not_installed("clusterProfiler")
    use_fixture_classification()

    unfiltered <- run_gmt_ora()
    filtered   <- run_gmt_ora(exclude_classes = "Organismal Systems")

    # The exclusion is a reporting layer. The tested universe, the Fisher test
    # and the BH family behind these numbers were all settled before it ran, so
    # dropping a row must not move a number on the rows that stay.
    expect_identical(filtered[filtered$ID == "hsa00010", ],
                     unfiltered[unfiltered$ID == "hsa00010", ])
})

test_that("a mixed GMT loses only its excluded KEGG rows", {
    skip_if_not_installed("clusterProfiler")
    use_fixture_classification()

    unfiltered <- run_gmt_ora()
    filtered   <- run_gmt_ora(exclude_classes = c("Organismal Systems", "Metabolism"))

    # Excluding every class the hierarchy lists still leaves the GO, Pfam and
    # custom rows whole: they have no class to be excluded by, so they come
    # back with the same values, the same column shape and the same order.
    expect_identical(filtered[filtered$ID %in% NON_KEGG_SETS, ],
                     unfiltered[unfiltered$ID %in% NON_KEGG_SETS, ])
    expect_setequal(filtered$ID, NON_KEGG_SETS)
})

test_that("a KEGG accession from another organism is not this run's to exclude", {
    skip_if_not_installed("clusterProfiler")
    use_fixture_classification()

    # The same table read with no organism code: "hsa04260" is then not a
    # recognised accession, so the #201 identity contract keeps it rather than
    # guessing that any three letters before five digits must be KEGG.
    filtered <- run_gmt_ora(exclude_classes = "Organismal Systems", kegg_org = NULL)

    expect_true("hsa04260" %in% filtered$ID)
})


# ---- the branch is chosen before the exclusion runs -------------------------

test_that("emptying the adjusted branch does not fall through to the raw-p fallback", {
    skip_if_not_installed("clusterProfiler")
    use_fixture_classification()

    # A cutoff tight enough that only the two KEGG maps clear it, while the
    # non-KEGG sets still sit under raw p < 0.05. Filtering before the branch
    # decision would empty the adjusted hits, drop into the fallback, and hand
    # back the very rows the adjusted branch had already rejected.
    strict <- run_gmt_ora(pval_cutoff = 1e-5)
    expect_setequal(strict$ID, KEGG_SETS)

    open <- run_gmt_ora()
    expect_true(all(open$pvalue[open$ID %in% NON_KEGG_SETS] < 0.05))

    filtered <- run_gmt_ora(exclude_classes = c("Organismal Systems", "Metabolism"),
                            pval_cutoff = 1e-5)

    # Nothing survived the branch that ran, which is the honest answer. The
    # fallback's wider row set is not it.
    expect_null(filtered)
})

test_that("the fallback branch is filtered too when it is the one that ran", {
    skip_if_not_installed("clusterProfiler")
    use_fixture_classification()

    # No adjusted hit at a cutoff this small, so the wrapper announces the
    # raw-p fallback -- and that table is filtered on its way out as well,
    # rather than the exclusion applying to only one of the two exits.
    expect_message(
        filtered <- run_multi_ora_enricher(
            gmt_sig(), gmt_universe(), gmt_term2gene(), gmt_term2name(),
            label = "pooled", pval_cutoff = 1e-300,
            exclude_classes = "Organismal Systems", kegg_org = "hsa"),
        "padj too strict")

    expect_false("hsa04260" %in% filtered$ID)
    expect_true("hsa00010" %in% filtered$ID)
})


# ---- the paths that must not reach the network ------------------------------

test_that("no exclusion list means the classification is never consulted", {
    skip_if_not_installed("clusterProfiler")
    # The default is lazy on purpose: almost every project excludes nothing,
    # and those runs must not pay for -- or fail on -- a KEGG fetch.
    local_gmt_stubs(list(kegg_pathway_categories = function(...)
        stop("classification fetched with nothing to exclude")))

    unfiltered <- run_gmt_ora()
    # The wrapper catches its own errors, so a fetch that did happen would
    # surface as a silent NULL rather than as a failure. Pin the rows first.
    expect_setequal(unfiltered$ID, c(KEGG_SETS, NON_KEGG_SETS))

    expect_identical(run_gmt_ora(exclude_classes = character(0)), unfiltered)
    expect_identical(run_gmt_ora(exclude_classes = list()), unfiltered)
})

test_that("a classification handed in is used as it is, and nothing is fetched", {
    skip_if_not_installed("clusterProfiler")
    fetches <- 0L
    local_gmt_stubs(list(kegg_pathway_categories = function(...) {
        fetches <<- fetches + 1L
        brite_classes()
    }))

    supplied <- suppressMessages(run_multi_ora_enricher(
        gmt_sig(), gmt_universe(), gmt_term2gene(), gmt_term2name(),
        label = "pooled", exclude_classes = "Organismal Systems",
        kegg_org = "hsa", classification = brite_classes()))

    expect_equal(fetches, 0L)
    expect_false("hsa04260" %in% supplied$ID)

    # ...and the lazy default is still the default: a caller that passes
    # nothing resolves it itself, exactly once for that call.
    defaulted <- run_gmt_ora(exclude_classes = "Organismal Systems")

    expect_equal(fetches, 1L)
    expect_identical(defaulted, supplied)
})

test_that("an unavailable classification keeps everything and says so", {
    skip_if_not_installed("clusterProfiler")
    local_gmt_stubs(list(kegg_pathway_categories = function(...) NULL))

    expect_message(
        filtered <- run_multi_ora_enricher(
            gmt_sig(), gmt_universe(), gmt_term2gene(), gmt_term2name(),
            label = "pooled", exclude_classes = "Organismal Systems",
            kegg_org = "hsa"),
        "keeping all")

    # Fail open. An unreachable hierarchy must not read as "nothing enriched".
    expect_true("hsa04260" %in% filtered$ID)
})


# ---- one filter, at the producer --------------------------------------------

test_that("the summary consumes what it is given and applies no exclusion of its own", {
    # Filtering lives at the producer so the summary, the CSV and the pooled
    # figure cannot disagree about which pathways exist. Handed an unfiltered
    # table the summary keeps the excluded class, which is exactly why the
    # wiring test below matters.
    ora <- data.frame(
        pathway = c("Glycolysis / Gluconeogenesis", "Cardiac muscle contraction"),
        ID      = c("hsa00010", "hsa04260"),
        pvalue  = c(0.001, 0.002),
        padj    = c(0.01, 0.02),
        Count   = c(5L, 4L),
        stringsAsFactors = FALSE
    )

    summarise <- function(df) {
        build_multi_ora_summary(df, list(transcriptomics = df), NULL)
    }

    expect_setequal(summarise(ora)$ID, ora$ID)

    # Given the already-filtered table, the excluded class is simply absent:
    # the summary neither restores it nor filters a second time.
    expect_identical(summarise(ora[ora$ID != "hsa04260", , drop = FALSE])$ID,
                     "hsa00010")
})

# ---- one classification, resolved once for the whole run --------------------
#
# kegg_pathway_categories() caches a successful fetch and nothing else: a
# failed lookup returns NULL without writing, so it is paid again every time it
# is asked for. Letting each producer fall through to .exclude_kegg_classes()'s
# lazy default therefore cost one endpoint timeout for the pooled table and one
# more per omics layer -- and this path refuses to run with fewer than two --
# before the run gave up and kept everything. Resolved once in
# run_multi_ora_gmt() and handed to every producer instead.

# A GMT in the shape read_gmt() wants: term, description, then members.
write_gmt_fixture <- function(terms, env = parent.frame()) {
    path <- withr::local_tempfile(fileext = ".gmt", .local_envir = env)
    writeLines(vapply(names(terms),
                      function(t) paste(c(t, t, terms[[t]]), collapse = "\t"),
                      character(1)),
               path)
    path
}

# Only the keys run_multi_ora_gmt() reads. Absolute temp paths pass through
# resolve_input_path() untouched, so no project layout is needed.
gmt_run_config <- function(rna_gmt, prot_gmt, exclude) {
    list(
        global = list(organism = "human"),
        modes = list(
            rna        = list(pathway = list(gmt_file = rna_gmt)),
            proteomics = list(pathway = list(gmt_file = prot_gmt)),
            multiomics = list(enrichment =
                                  list(exclude_pathway_classes = exclude))
        )
    )
}

# Drives a two-layer GMT run and reports what the producers were handed.
# extract_de_tables() and run_multi_ora_enricher() are both stubbed: what is
# under test is how many times the classification is resolved and which object
# each producer receives, not DE extraction or enrichment. With every producer
# returning NULL the run stops at the empty summary, so nothing is written.
capture_gmt_producers <- function(exclude, env = parent.frame()) {
    state <- new.env(parent = emptyenv())
    state$fetches <- 0L
    state$seen <- list()

    local_gmt_stubs(list(
        kegg_pathway_categories = function(...) {
            state$fetches <- state$fetches + 1L
            brite_classes()
        },
        extract_de_tables = function(de_data, omics_type,
                                     harmonization_res = NULL) {
            list(contrast = data.frame(
                feature_id = sprintf("g%03d", 1:6),
                pvalue     = rep(0.001, 6),
                padj       = rep(0.001, 6),
                stringsAsFactors = FALSE))
        },
        run_multi_ora_enricher = function(..., classification = NULL) {
            state$seen <- c(state$seen, list(classification))
            NULL
        }
    ), env = env)

    config <- gmt_run_config(
        write_gmt_fixture(list(hsa00010 = sprintf("g%03d", 1:6)), env = env),
        write_gmt_fixture(list(hsa04260 = sprintf("g%03d", 1:6)), env = env),
        exclude)

    suppressMessages(run_multi_ora_gmt(
        de_results = list(transcriptomics = list(), proteomics = list()),
        harmonization_res = NULL, config = config,
        out_dir = withr::local_tempdir(.local_envir = env)))

    list(fetches = state$fetches, seen = state$seen)
}

test_that("one GMT run resolves the classification once for all its producers", {
    got <- capture_gmt_producers("Organismal Systems")

    # Pooled plus one per omics layer -- three producers, one resolution.
    expect_equal(length(got$seen), 3L)
    expect_equal(got$fetches, 1L)

    # And all three were handed the same table, so the pooled and per-omics
    # results cannot be filtered against different hierarchies.
    expect_identical(got$seen[[1]], got$seen[[2]])
    expect_identical(got$seen[[2]], got$seen[[3]])
    expect_identical(got$seen[[1]], brite_classes())
})

test_that("a run with nothing to exclude resolves no classification at all", {
    got <- capture_gmt_producers(character(0))

    expect_equal(length(got$seen), 3L)
    expect_equal(got$fetches, 0L)
    # Every producer is told there is nothing, rather than left to find out.
    expect_true(all(vapply(got$seen, is.null, logical(1))))
})

test_that("run_multi_ora_gmt resolves the exclusion once and passes it to every producer", {
    src <- paste(deparse(body(run_multi_ora_gmt)), collapse = " ")

    expect_true(grepl(".excluded_pathway_classes(config)", src, fixed = TRUE))
    # resolve_kegg_org_code() rather than get_kegg_organism(): this branch runs
    # when either the KEGG code or the OrgDb is missing, so an organism KEGG
    # knows but Bioconductor has no OrgDb for arrives here with a usable code.
    # The narrower registry would leave its accessions unclassified, and so
    # unexcluded, on precisely the runs this path exists to serve.
    expect_true(grepl("resolve_kegg_org_code(config$global$organism)",
                      src, fixed = TRUE))

    # Pooled and per-omics both, or one layer's KEGG rows outlive the other's.
    count_in_src <- function(needle) {
        lengths(regmatches(src, gregexpr(needle, src, fixed = TRUE)))[[1]]
    }
    expect_equal(count_in_src("run_multi_ora_enricher"), 2L)
    expect_equal(count_in_src("exclude_classes = exclude_classes"), 2L)
    expect_equal(count_in_src("classification = pathway_classes"), 2L)

    expect_true(all(c("exclude_classes", "kegg_org", "classification") %in%
                    names(formals(run_multi_ora_enricher))))
})
