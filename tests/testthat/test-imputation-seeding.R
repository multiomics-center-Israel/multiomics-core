# tests/testthat/test-imputation-seeding.R
#
# Reproducibility of the proteomics imputation sequence.
#
# params$seed is the single source. The QC draw in preprocess_proteomics()
# takes params$seed with no offset, and DE imputation run i takes
# params$seed + i, so no two draws share a seed and the whole sequence is
# reproducible from one number.
#
# The bug this pins: impute_proteomics_qrilc() and impute_proteomics_dep2()
# used to call set.seed() with a configured constant INSIDE each call, which
# overwrote the per-run seed. Every repetition of those stochastic methods then
# returned the same matrix, so no_repetitions > 1 looked like multiple
# imputation while repeating one draw.
#
# These tests assert behaviour, not implementation: they compare matrices
# across runs and across executions rather than looking for set.seed() calls.
#
# All fixtures are synthetic (f1..f12, S1..S6).

seed_cfg <- function(method = "perseus_like", n_reps = 3, seed = 1,
                     dep2_method = NULL) {
    imp <- list(method = method, multi_imputation = TRUE,
                no_repetitions = n_reps, min_no_passed = 1)
    if (!is.null(dep2_method)) imp$dep2_method <- dep2_method
    list(
        params = list(seed = seed),
        modes = list(proteomics = list(
            imputation = imp,
            id_columns = list(protein_id = "Protein.Group", sample_col = "SampleID"),
            normalization = list(method = "none")
        ))
    )
}

# A matrix with enough spread for a stochastic draw to be visibly different
# between runs, and a fixed block of missing values to fill.
seed_expr <- function(n_feat = 12, n_samp = 6, seed = 99) {
    withr::with_seed(seed, {
        m <- matrix(stats::rnorm(n_feat * n_samp, mean = 20, sd = 2),
                    nrow = n_feat,
                    dimnames = list(paste0("f", seq_len(n_feat)),
                                    paste0("S", seq_len(n_samp))))
    })
    # Missing values in a fixed pattern, leaving every column with plenty of
    # observed values so the imputation methods have a distribution to work from.
    m[1:3, 1] <- NA
    m[4:5, 3] <- NA
    m[6, 5]   <- NA
    m
}

observed_cells <- function(m) !is.na(m)

# The methods that draw independently per repetition, and the ones that do not.
stochastic_cases <- list(
    perseus_like = seed_cfg("perseus_like"),
    qrilc        = seed_cfg("qrilc"),
    dep2_minprob = seed_cfg("dep2", dep2_method = "MinProb")
)
deterministic_cases <- list(
    minval      = seed_cfg("minval"),
    dep2_mindet = seed_cfg("dep2", dep2_method = "MinDet"),
    none        = seed_cfg("none")
)


# =============================================================================
# The same seed reproduces the same ordered sequence
# =============================================================================

test_that("the same seed and config reproduce the whole run sequence", {
    expr <- seed_expr()
    for (nm in names(stochastic_cases)) {
        cfg <- stochastic_cases[[nm]]
        # Draw something in between, so a run that depended on ambient RNG
        # state rather than on its own seed would come back different.
        first  <- make_imputations_proteomics(expr, cfg)
        invisible(stats::runif(10))
        second <- make_imputations_proteomics(expr, cfg)

        expect_equal(first, second, info = nm)
        expect_length(first, 3L)
    }
})

test_that("deterministic methods reproduce their sequence too", {
    expr <- seed_expr()
    for (nm in names(deterministic_cases)) {
        cfg <- deterministic_cases[[nm]]
        expect_equal(make_imputations_proteomics(expr, cfg),
                     make_imputations_proteomics(expr, cfg), info = nm)
    }
})

test_that("a different seed gives different stochastic draws", {
    expr <- seed_expr()
    for (nm in names(stochastic_cases)) {
        cfg <- stochastic_cases[[nm]]
        alt <- cfg
        alt$params$seed <- cfg$params$seed + 1000L

        a <- make_imputations_proteomics(expr, cfg)[[1]]
        b <- make_imputations_proteomics(expr, alt)[[1]]
        miss <- is.na(expr)
        expect_false(isTRUE(all.equal(a[miss], b[miss])), info = nm)
    }
})


# =============================================================================
# Repetitions of a stochastic method are genuinely different draws
# =============================================================================

test_that("a stochastic method draws differently in every repetition", {
    # The regression test for #229. Before the fix, qrilc and dep2 MinProb
    # returned byte-identical matrices here because they reseeded internally
    # from a configured constant; perseus_like was already correct.
    expr <- seed_expr()
    miss <- is.na(expr)

    for (nm in names(stochastic_cases)) {
        runs <- make_imputations_proteomics(expr, stochastic_cases[[nm]])
        expect_length(runs, 3L)

        # At least one imputed cell differs between the first two runs...
        expect_false(isTRUE(all.equal(runs[[1]][miss], runs[[2]][miss])),
                     info = paste(nm, "runs 1 and 2"))
        # ...and no two runs are wholly identical.
        for (i in seq_along(runs)) {
            for (j in seq_along(runs)) {
                if (i >= j) next
                expect_false(identical(runs[[i]], runs[[j]]),
                             info = sprintf("%s runs %d and %d", nm, i, j))
            }
        }
    }
})

test_that("a deterministic method repeats the same matrix, by design", {
    expr <- seed_expr()
    for (nm in names(deterministic_cases)) {
        runs <- make_imputations_proteomics(expr, deterministic_cases[[nm]])
        expect_length(runs, 3L)
        for (i in seq_along(runs)) {
            expect_identical(runs[[i]], runs[[1]], info = paste(nm, "run", i))
        }
    }
})


# =============================================================================
# Imputation only ever fills the gaps
# =============================================================================

test_that("measured values are untouched in every run of every method", {
    expr <- seed_expr()
    keep <- observed_cells(expr)

    for (cases in list(stochastic_cases, deterministic_cases)) {
        for (nm in names(cases)) {
            runs <- make_imputations_proteomics(expr, cases[[nm]])
            for (i in seq_along(runs)) {
                expect_equal(runs[[i]][keep], expr[keep],
                             info = sprintf("%s run %d", nm, i))
                expect_identical(dimnames(runs[[i]]), dimnames(expr), info = nm)
            }
        }
    }
})


# =============================================================================
# The seed contract: no collision between the QC draw and the DE runs
# =============================================================================

test_that("the QC draw uses params$seed and the DE runs use params$seed + i", {
    # Pins the documented offsets from the outside: reproducing run i by hand
    # from params$seed + i must give exactly what the wrapper produced, and the
    # QC draw's seed (offset 0) must not reproduce any of them.
    expr <- seed_expr()
    cfg <- seed_cfg("perseus_like", n_reps = 3, seed = 7)
    runs <- make_imputations_proteomics(expr, cfg)
    miss <- is.na(expr)

    for (i in seq_along(runs)) {
        set.seed(7L + i)
        by_hand <- impute_proteomics(expr, cfg = cfg$modes$proteomics)
        expect_equal(runs[[i]], by_hand, info = paste("run", i))
    }

    # Offset 0 is reserved for the QC draw, so it matches no DE run.
    set.seed(7L)
    qc_draw <- impute_proteomics(expr, cfg = cfg$modes$proteomics)
    for (i in seq_along(runs)) {
        expect_false(isTRUE(all.equal(qc_draw[miss], runs[[i]][miss])),
                     info = paste("QC draw vs run", i))
    }
})

test_that("a config without params$seed still runs", {
    # seed_base had no default, so set.seed(integer(0)) errored out.
    cfg <- seed_cfg("perseus_like", n_reps = 2)
    cfg$params <- NULL
    expect_no_error(runs <- make_imputations_proteomics(seed_expr(), cfg))
    expect_length(runs, 2L)
})


# =============================================================================
# The preprocessing draw
# =============================================================================

test_that("the QC draw is seeded from params$seed at its call site", {
    # Previously unseeded for perseus_like: the QC draw came from whatever RNG
    # state preceded it, so the default method's expr_imp_single could not be
    # reproduced at all.
    #
    # Asserted against the source rather than by calling preprocess_proteomics()
    # end to end. Driving that function needs a loader-shaped `inputs` and takes
    # the draw through filtering and normalization first, none of which this
    # change touches; what would silently regress is the seed line being
    # dropped or moved after the call, and that is what this pins. The seeding
    # CONTRACT -- offset 0 for the QC draw, offset i for DE run i -- is pinned
    # behaviourally in the seed-collision test above. Same idiom as
    # test-pca-protein-subsets.R, which pins the ordering of this same call.
    f <- testthat::test_path("..", "..", "R", "domain", "proteomics", "04_preprocess.R")
    skip_if(!file.exists(f), "preprocessing source not found from the test working directory")
    src <- readLines(f, warn = FALSE)

    line_of <- function(pattern) {
        hits <- grep(pattern, src, fixed = TRUE)
        expect_length(hits, 1)
        hits
    }
    seed_line <- line_of("set.seed(as.integer(config$params$seed %||% 1L))")
    imp_line  <- line_of("imp_res <- impute_proteomics(")

    # Immediately before, so nothing can consume the stream in between.
    expect_equal(imp_line, seed_line + 1L)

    # Read from the full config: `cfg` in that function is mode-level and has
    # no params, so cfg$params$seed would silently be NULL and every project
    # would land on the same fallback.
    expect_match(src[seed_line], "config$params$seed", fixed = TRUE)
    expect_false(grepl("cfg$params$seed", src[seed_line], fixed = TRUE))
})


# =============================================================================
# The deprecated per-method seeds
# =============================================================================

test_that("a deprecated seed key warns and is not required", {
    base <- seed_cfg("qrilc")$modes$proteomics
    base$imputation$qrilc_random_seed <- 7

    expect_warning(validate_proteomics_config(base), "deprecated and ignored")
    expect_warning(validate_proteomics_config(base), "params\\$seed")

    # dep2 no longer requires its seed key: this used to abort validation.
    d <- seed_cfg("dep2", dep2_method = "MinDet")$modes$proteomics
    expect_silent(validate_proteomics_config(d))

    d$imputation$dep2_random_seed <- 1
    expect_warning(validate_proteomics_config(d), "dep2_random_seed")
})

test_that("a deprecated seed key no longer changes any result", {
    expr <- seed_expr()
    cfg <- seed_cfg("qrilc")
    with_key <- cfg
    with_key$modes$proteomics$imputation$qrilc_random_seed <- 12345

    expect_equal(make_imputations_proteomics(expr, cfg),
                 make_imputations_proteomics(expr, with_key))
})
