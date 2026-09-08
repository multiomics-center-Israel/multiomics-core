# tests/testthat/test-mummichog-gsea-parity.R
#
# Parity test: our MetaboAnalyst-style GSEA engine must reproduce the PINNED
# MetaboAnalystR implementation exactly on a stored reference.
#
#   repo   xia-lab/MetaboAnalystR
#   commit 398476ae2a0c996e390925fa62adc46b5ecde334
#   path   PerformPSEA -> .init.RT.Permutations -> .compute.mummichog.RT.fgsea
#            -> fgsea2 -> my.fgsea -> .run_fgsea_inner
#
# The reference was produced by running that very function (see
# fixtures/mummichog_gsea_parity/generate_reference.R and REFERENCE.md). It is
# NEVER regenerated here: CI compares against the committed fixture.
#
# Deterministic quantities (tested pathway set, sizes, ES, leading edge) are
# compared for EXACT equality. ES/NES/p/padj are floating point and go through a
# tight tolerance only; the permutation tallies themselves are deterministic
# because upstream seeds each batch explicitly (set.seed(123) -> sample.int).
#
# fgsea VERSION MATTERS: the engine calls fgsea internals
# (calcGseaStatCumulativeBatch), so the reference is generated under the fgsea
# series renv.lock pins, and these tests refuse to compare across series rather
# than reporting a cross-version match as parity.

fixture_dir <- file.path(
  normalizePath(if (dir.exists("R")) "." else "../..", mustWork = FALSE),
  "tests", "testthat", "fixtures", "mummichog_gsea_parity")
ref_file <- file.path(fixture_dir, "reference.rds")

# Tight but non-zero: same primitives, same seeds, so only accumulation order
# can differ at all.
PARITY_TOL <- 1e-8

# x.y series of a package version ("1.36.2" -> "1.36").
parity_series <- function(v) {
  paste(strsplit(as.character(v), ".", fixed = TRUE)[[1]][1:2], collapse = ".")
}

# Skip unless the running fgsea is in the same series the fixture was generated
# under. A match found on a different series would not be evidence of parity.
skip_unless_parity_fgsea <- function(ref) {
  skip_if_not_installed("fgsea")
  got <- parity_series(utils::packageVersion("fgsea"))
  skip_if_not(identical(got, ref$fgsea_series),
              sprintf(paste0("fgsea %s is in series %s; the parity reference ",
                             "was generated under series %s (renv.lock pins %s)"),
                      utils::packageVersion("fgsea"), got, ref$fgsea_series,
                      ref$fgsea_locked))
}


test_that("the parity fixture is present and records its provenance", {
  skip_if_not(file.exists(ref_file), "parity reference fixture not found")
  ref <- readRDS(ref_file)

  expect_identical(ref$metaboanalystr_commit,
                   "398476ae2a0c996e390925fa62adc46b5ecde334")
  # the pinned commit our code claims to follow must be the one that generated
  # the reference
  expect_identical(ref$metaboanalystr_commit, .MMC_GSEA_REF_COMMIT)
  expect_identical(ref$nperm, 100L)
  expect_equal(ref$gsea_param, 1)
  expect_equal(ref$min_size, 1)
  expect_equal(ref$max_size, Inf)
  expect_identical(ref$set_seed, 123L)

  # The reference must have been generated under the fgsea series this project
  # locks, because the engine depends on fgsea internals.
  root <- normalizePath(if (dir.exists("R")) "." else "../..", mustWork = FALSE)
  lock <- jsonlite::fromJSON(file.path(root, "renv.lock"))
  locked <- lock$Packages$fgsea$Version
  expect_identical(ref$fgsea_locked, locked)
  expect_identical(ref$fgsea_series, parity_series(locked))
  expect_match(ref$fgsea_branch, "post-1\\.24\\.0")
  # and the exact build is recorded, not just the series
  expect_true(is.character(ref$fgsea_version) && nzchar(ref$fgsea_version))
})


test_that("our engine reproduces the pinned MetaboAnalystR result exactly", {
  skip_if_not(file.exists(ref_file), "parity reference fixture not found")
  skip_if_not_installed("fastmatch")
  skip_if_not_installed("data.table")
  skip_if_not_installed("BiocParallel")
  ref <- readRDS(ref_file)
  skip_unless_parity_fgsea(ref)
  # the exact build actually running, for the record
  expect_true(nzchar(as.character(utils::packageVersion("fgsea"))))

  got <- mmc_gsea_metaboanalyst(
    pathways   = ref$pathways,
    stats      = ref$stats,
    ranks      = ref$ranks,
    n_perm     = ref$nperm,
    min_size   = ref$min_size,
    max_size   = ref$max_size,
    gsea_param = ref$gsea_param,
    seed       = ref$set_seed
  )
  exp <- ref$reference

  # --- tested pathway universe: exact, same order --------------------------
  expect_identical(got$pathway, exp$pathway)
  # PW_absent has no detected EC -> size 0 -> excluded by minSize = 1
  expect_false("PW_absent" %in% got$pathway)
  # PW_single (size 1) IS tested, because minSize = 1
  expect_true("PW_single" %in% got$pathway)

  # --- deterministic quantities: exact ------------------------------------
  expect_identical(got$size, exp$size)
  # upstream writes an empty leading edge as "", we carry NA; compare on the
  # same representation, contents and order otherwise identical
  got_lead <- ifelse(is.na(got$leadingEdge), "", got$leadingEdge)
  expect_identical(got_lead, exp$leadingEdgeMatched)

  # --- numeric quantities: tight tolerance --------------------------------
  expect_equal(got$ES,  exp$ES,  tolerance = PARITY_TOL)
  expect_equal(got$NES, exp$NES, tolerance = PARITY_TOL)
  expect_equal(got$pval, exp$pval, tolerance = PARITY_TOL)
  expect_equal(got$padj, exp$padj, tolerance = PARITY_TOL)
  expect_equal(got$nMoreExtreme, exp$nMoreExtreme, tolerance = PARITY_TOL)
})


test_that("BH adjustment is taken over exactly the tested pathway universe", {
  skip_if_not(file.exists(ref_file), "parity reference fixture not found")
  ref <- readRDS(ref_file)
  skip_unless_parity_fgsea(ref)

  got <- mmc_gsea_metaboanalyst(ref$pathways, ref$stats, ref$ranks,
                                n_perm = ref$nperm, min_size = ref$min_size,
                                max_size = ref$max_size,
                                gsea_param = ref$gsea_param,
                                seed = ref$set_seed)

  # padj is BH over the kept pathways only — the dropped PW_absent must not
  # inflate the universe.
  expect_equal(got$padj, stats::p.adjust(got$pval, method = "fdr"),
               tolerance = PARITY_TOL)
  expect_equal(nrow(got), length(ref$pathways) - 1L)
  # and it matches the reference's own adjustment
  expect_equal(got$padj, ref$reference$padj, tolerance = PARITY_TOL)
})


test_that("the tied-score ES = 0 fallback is reproduced and flagged", {
  skip_if_not(file.exists(ref_file), "parity reference fixture not found")
  ref <- readRDS(ref_file)
  skip_unless_parity_fgsea(ref)

  got <- mmc_gsea_metaboanalyst(ref$pathways, ref$stats, ref$ranks,
                                n_perm = ref$nperm, min_size = ref$min_size,
                                max_size = ref$max_size,
                                gsea_param = ref$gsea_param,
                                seed = ref$set_seed)

  # E02 and E03 share a score, so PW_top and PW_overlap_b pass a
  # non-strictly-increasing selectedStats; upstream's tryCatch substitutes
  # ES = 0 and an empty leading edge, and the reference carries exactly that.
  tied <- c("PW_top", "PW_overlap_b")
  expect_true(all(got$es_defaulted[got$pathway %in% tied]))
  expect_false(any(got$es_defaulted[!got$pathway %in% tied]))
  expect_equal(got$ES[got$pathway %in% tied], c(0, 0))
  expect_equal(ref$reference$ES[ref$reference$pathway %in% tied], c(0, 0))
  # upstream emitted one warning per defaulted pathway
  expect_length(ref$upstream_warnings, 2L)

  # the flag is generic and the CAUSE is recorded beside it, classified from
  # the actual condition (duplicate ranked positions), not assumed
  expect_true(all(c("es_defaulted", "es_reason") %in% names(got)))
  expect_true(all(grepl("^tied EC scores \\(1 duplicate ranked position\\)$",
                        got$es_reason[got$pathway %in% tied])))
  expect_true(all(is.na(got$es_reason[!got$pathway %in% tied])))
})


test_that("our ranked-input construction is the one the reference was built on", {
  skip_if_not(file.exists(ref_file), "parity reference fixture not found")
  ref <- readRDS(ref_file)

  ours <- mmc_gsea_ranked_inputs(ref$ec_scores)

  # `stats`: one entry per UNIQUE score, descending by SIGNED score, named by
  # the "; "-joined ids that share it.
  expect_identical(names(ours$stats), names(ref$stats))
  expect_equal(unname(ours$stats), unname(ref$stats), tolerance = PARITY_TOL)
  expect_true(all(diff(unname(ours$stats)) < 0))          # strictly descending
  expect_identical(unname(ours$stats[["E02; E03"]]), 3.2) # the tie group

  # `ranks`: one entry per EC, valued with its tie-group index.
  expect_identical(names(ours$ranks), names(ref$ranks))
  expect_equal(unname(ours$ranks), unname(ref$ranks), tolerance = PARITY_TOL)
  expect_equal(unname(ours$ranks[["E02"]]), unname(ours$ranks[["E03"]]))
  expect_length(ours$ranks, length(ref$ec_scores))        # every EC present
  expect_length(ours$stats, length(ref$ec_scores) - 1L)   # one tie collapsed
})


test_that("the MetaboAnalyst-equivalent defaults are the shipped defaults", {
  # These are the pinned reference's values; changing them silently would break
  # the semantics the parity fixture locks down.
  expect_identical(.MMC_GSEA_DEFAULTS$n_perm, 100L)
  expect_equal(.MMC_GSEA_DEFAULTS$gsea_param, 1)
  expect_equal(.MMC_GSEA_DEFAULTS$min_size, 1)
  expect_equal(.MMC_GSEA_DEFAULTS$max_size, Inf)
  expect_identical(.MMC_GSEA_DEFAULTS$seed, 123L)

  # and they are what the engine and the runner actually default to
  eng <- formals(mmc_gsea_metaboanalyst)
  expect_identical(eng$n_perm,     quote(.MMC_GSEA_DEFAULTS$n_perm))
  expect_identical(eng$min_size,   quote(.MMC_GSEA_DEFAULTS$min_size))
  expect_identical(eng$max_size,   quote(.MMC_GSEA_DEFAULTS$max_size))
  expect_identical(eng$gsea_param, quote(.MMC_GSEA_DEFAULTS$gsea_param))
  expect_identical(eng$seed,       quote(.MMC_GSEA_DEFAULTS$seed))

  run <- formals(run_mummichog_gsea)
  expect_identical(run$n_perm,     quote(.MMC_GSEA_DEFAULTS$n_perm))
  expect_identical(run$min_size,   quote(.MMC_GSEA_DEFAULTS$min_size))
  expect_identical(run$max_size,   quote(.MMC_GSEA_DEFAULTS$max_size))
  expect_identical(run$gsea_param, quote(.MMC_GSEA_DEFAULTS$gsea_param))
  expect_identical(run$seed,       quote(.MMC_GSEA_DEFAULTS$seed))
})


test_that("minSize = 1 keeps size-1 pathways and maxSize = Inf keeps large ones", {
  skip_if_not(file.exists(ref_file), "parity reference fixture not found")
  ref <- readRDS(ref_file)
  skip_unless_parity_fgsea(ref)

  # default (1 / Inf): every pathway with at least one detected EC is tested
  kept <- mmc_gsea_metaboanalyst(ref$pathways, ref$stats, ref$ranks,
                                 n_perm = 100L, seed = 123L)$pathway
  expect_true("PW_single" %in% kept)
  expect_setequal(kept, setdiff(names(ref$pathways), "PW_absent"))

  # a deliberate override drops it again — proving the filter is live, and that
  # 2/500 (the values this branch used to default to) is a real change
  narrowed <- mmc_gsea_metaboanalyst(ref$pathways, ref$stats, ref$ranks,
                                     n_perm = 100L, seed = 123L,
                                     min_size = 2, max_size = 4)$pathway
  expect_false("PW_single" %in% narrowed)       # size 1 < min_size
  expect_false("PW_overlap_a" %in% narrowed)    # size 5 > max_size
})


test_that("an override of the reference defaults is reported as a deviation", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)

  clean <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "HL_vs_LL")
  expect_length(clean$deviations, 0L)
  expect_identical(clean$params$n_perm, 100L)
  expect_equal(clean$params$min_size, 1)
  expect_equal(clean$params$max_size, Inf)
  expect_identical(clean$reference$commit, .MMC_GSEA_REF_COMMIT)

  dev <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "HL_vs_LL",
                            min_size = 3, n_perm = 200L)
  expect_true(any(grepl("min_size = 3", dev$deviations)))
  expect_true(any(grepl("n_perm = 200", dev$deviations)))
})


# ---------------------------------------------------------------------------
# ES fallback classification: a tie must be called a tie, and nothing else may
# ---------------------------------------------------------------------------

test_that("a non-tie fallback is classified as such, not as a tie", {
  skip_if_not(file.exists(ref_file), "parity reference fixture not found")
  ref <- readRDS(ref_file)
  skip_unless_parity_fgsea(ref)

  # A pathway holding one EC per tie group covers EVERY ranked position with no
  # duplicate index, so calcGseaStat refuses it for a completely different
  # reason ("GSEA statistic is not defined when all genes are selected").
  # Upstream still substitutes ES = 0 — but this is NOT a tie.
  one_per_group <- names(ref$ranks)[!duplicated(ref$ranks)]
  expect_length(one_per_group, length(ref$stats))          # == N positions
  expect_false(any(duplicated(ref$ranks[one_per_group])))  # ...and no ties

  got <- mmc_gsea_metaboanalyst(
    c(ref$pathways, list(PW_all = one_per_group)),
    ref$stats, ref$ranks, n_perm = ref$nperm, min_size = ref$min_size,
    max_size = ref$max_size, gsea_param = ref$gsea_param, seed = ref$set_seed)

  row <- got[got$pathway == "PW_all", ]
  expect_equal(nrow(row), 1L)
  expect_true(row$es_defaulted)
  expect_equal(row$ES, 0)                                  # upstream's fallback
  expect_identical(row$es_reason, "pathway selects every ranked position")
  expect_false(grepl("tie", row$es_reason, ignore.case = TRUE))

  # the genuinely tied pathways are still labelled as tied, side by side
  expect_true(all(grepl("^tied EC scores",
                        got$es_reason[got$pathway %in% c("PW_top",
                                                         "PW_overlap_b")])))
})

test_that("the fallback reason classifier covers each condition it claims", {
  reason <- .mmc_gsea_fallback_reason

  # duplicate ranked positions -> tie, counted
  expect_identical(reason(c(1, 1, 3), "all(head(S, -1) < tail(S, -1)) is not TRUE"),
                   "tied EC scores (1 duplicate ranked position)")
  expect_identical(reason(c(2, 2, 2, 5), "anything"),
                   "tied EC scores (2 duplicate ranked positions)")

  # whole-ranking selection -> its own reason, never a tie
  expect_identical(
    reason(1:5, "GSEA statistic is not defined when all genes are selected"),
    "pathway selects every ranked position")

  # upstream's own pre-check branch (no error was raised)
  expect_identical(reason(integer(0), NULL), "pathway has no ranked members")
  expect_identical(reason(c(1, NA), NULL),
                   "pathway members include a missing rank")

  # anything else is surfaced verbatim as UNEXPECTED, not explained away
  unexpected <- reason(c(1, 2, 3), "some future fgsea assertion blew up")
  expect_true(startsWith(unexpected, .MMC_GSEA_FALLBACK_UNEXPECTED))
  expect_match(unexpected, "some future fgsea assertion blew up", fixed = TRUE)
  expect_false(grepl("tie", unexpected, ignore.case = TRUE))
})

test_that("an unexpected fallback raises a warning rather than passing quietly", {
  skip_if_not(file.exists(ref_file), "parity reference fixture not found")
  ref <- readRDS(ref_file)
  skip_unless_parity_fgsea(ref)

  # Force calcGseaStat to fail for a reason we do not recognise, so the
  # unexpected-classification path is exercised end to end.
  local_mocked_bindings(
    calcGseaStat = function(...) stop("synthetic non-tie failure"),
    .package = "fgsea"
  )
  expect_warning(
    got <- mmc_gsea_metaboanalyst(ref$pathways, ref$stats, ref$ranks,
                                  n_perm = ref$nperm, min_size = ref$min_size,
                                  max_size = ref$max_size,
                                  gsea_param = ref$gsea_param,
                                  seed = ref$set_seed),
    "UNEXPECTED calcGseaStat failure")

  # tied pathways are still classified by their real condition...
  tied <- got$pathway %in% c("PW_top", "PW_overlap_b")
  expect_true(all(grepl("^tied EC scores", got$es_reason[tied])))
  # ...and the rest report the synthetic cause verbatim, as unexpected
  expect_true(all(startsWith(got$es_reason[!tied],
                             .MMC_GSEA_FALLBACK_UNEXPECTED)))
  expect_true(all(grepl("synthetic non-tie failure", got$es_reason[!tied],
                        fixed = TRUE)))
  expect_true(all(got$es_defaulted))
})
