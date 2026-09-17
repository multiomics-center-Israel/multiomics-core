# Regression tests: info$normalization$transform must name the transform that
# was ACTUALLY APPLIED to the matrix, not the one that was configured.
#
# mod_met_normalize_linear() (tss/pqn) and both EigenMS targets normalise on the
# LINEAR matrix and then call transform_metab(method = "log2") themselves,
# ignoring preprocessing$transform. Only the median path runs on mod_met_log()'s
# output, so only there do configured and applied coincide. mod_met_corrected()
# used to record norm_cfg$transform unconditionally, so a lane configured
# "PQN + glog10" reported glog10 while holding a log2 matrix.
#
# That label is not decoration. fc_to_log2_units()
# (R/domain/metabolomics/03_differential.R) reads it to put the fold change into
# log2 units, and multiplies by log2(10) when it says log10/glog10 -- so the
# mislabel inflated every reported log2FC on such a lane by 3.32x.
#
# A correlation gate cannot see this: scaling every logFC by one constant leaves
# r = 1.0. The check that sees it is a SLOPE -- median(logFC / implied) -- and
# that is what the end-to-end test below asserts.

met_cfg <- function(chosen_norm, transform, scaling = "auto") {
    list(modes = list(metabolomics = list(
        effects       = list(samples = "sample_id", color = "Group"),
        preprocessing = list(transform = transform, pseudocount = 1,
                             chosen_norm = chosen_norm, scaling = scaling,
                             biological_factor_col = "total_protein",
                             drift_correction = list(enabled = FALSE)),
        de            = list(method = "t_test_equal", p_cutoff = 0.05,
                             linear_fc_cutoff = 1)
    )))
}

met_fixture <- function(n_feat = 30, seed = 3) {
    set.seed(seed)
    samples <- c(sprintf("A%d", 1:5), sprintf("B%d", 1:5))
    meta <- data.frame(sample_id = samples,
                       Group = rep(c("A", "B"), each = 5),
                       # the per-sample covariate the bio_factor lane divides by
                       total_protein = seq(1, 2, length.out = 10),
                       stringsAsFactors = FALSE)
    mat <- matrix(stats::rlnorm(n_feat * length(samples), meanlog = 12, sdlog = 0.2),
                  nrow = n_feat,
                  dimnames = list(sprintf("F%03d", seq_len(n_feat)), samples))
    mat[1:10, meta$Group == "B"] <- mat[1:10, meta$Group == "B"] * 2
    list(mat = mat, meta = meta, row_data = NULL)
}

# Build the mod_met_corrected() output the targets pipeline would build.
corrected_for <- function(f, cfg) {
    logged <- mod_met_log(f, cfg)
    chosen <- cfg$modes$metabolomics$preprocessing$chosen_norm
    suppressWarnings(suppressMessages(mod_met_corrected(
        norm_tss        = if (chosen == "tss") mod_met_normalize_linear(f, "tss", cfg),
        norm_median     = if (chosen == "median") mod_met_normalize_log(logged, cfg),
        norm_pqn        = if (chosen == "pqn") mod_met_normalize_linear(f, "pqn", cfg),
        norm_bio_factor = if (chosen == "bio_factor") mod_met_normalize_bio_factor(f, cfg),
        logged          = logged,
        meta            = f$meta,
        out_dir         = withr::local_tempdir(),
        config          = cfg)))
}

# ---- 1. the recorded transform is the applied one ---------------------------

for (nm in c("tss", "pqn", "bio_factor")) {
    test_that(paste0("mod_met_corrected [", nm,
                     "]: records log2, the transform it actually applied"), {
        f   <- met_fixture()
        cfg <- met_cfg(nm, transform = "glog10")
        out <- corrected_for(f, cfg)

        ni <- out$info$normalization
        expect_equal(ni$transform, "log2")
        # the configured value is not lost, it is just not passed off as applied
        expect_equal(ni$transform_configured, "glog10")

        # and the matrix really is log2, not glog10: it reproduces
        # transform_metab(..., "log2") of the normalised linear matrix, and it
        # does NOT reproduce the glog10 of it.
        norm_lin <- switch(nm,
            tss        = norm_total_sum(f$mat),
            pqn        = norm_pqn(f$mat),
            bio_factor = norm_biological_factor(f$mat, f$meta,
                                                factor_col = "total_protein",
                                                sample_col = "sample_id"))
        expect_equal(unname(out$mat_pre_scale),
                     unname(transform_metab(norm_lin, method = "log2",
                                            pseudocount = 1)))
        expect_gt(max(abs(out$mat_pre_scale -
                          transform_metab(norm_lin, method = "glog10",
                                          pseudocount = 1))), 1)
    })
}

test_that("mod_met_corrected [median]: records the CONFIGURED transform, because that is the one applied", {
    # The median path is a log-shift on mod_met_log()'s output, so it carries
    # whatever transform was configured. This branch must not move.
    f   <- met_fixture()
    cfg <- met_cfg("median", transform = "glog10")
    out <- corrected_for(f, cfg)

    ni <- out$info$normalization
    expect_equal(ni$transform, "glog10")
    expect_equal(ni$transform_configured, "glog10")
    expect_equal(unname(out$mat_pre_scale),
                 unname(mod_met_normalize_log(mod_met_log(f, cfg), cfg)$mat))
})

test_that("mod_met_corrected [none]: records the CONFIGURED transform, because that is the one applied", {
    # chosen_norm = "none" means "the table was normalised upstream": no sample
    # normalisation runs and mod_met_corrected() reads mod_met_log()'s target
    # directly. That target applies the CONFIGURED transform, so here -- as on
    # the median path -- configured and applied coincide, and they must be
    # recorded as coinciding rather than by accident.
    f   <- met_fixture()
    cfg <- met_cfg("none", transform = "glog10")
    out <- corrected_for(f, cfg)

    ni <- out$info$normalization
    expect_equal(ni$transform, "glog10")
    expect_equal(ni$transform_configured, "glog10")

    # and the matrix really is the glog10 of the filtered matrix: no sample
    # normalisation, no log2 substituted for the configured transform.
    expect_equal(unname(out$mat_pre_scale), unname(mod_met_log(f, cfg)$mat))
    expect_equal(unname(out$mat_pre_scale),
                 unname(transform_metab(f$mat, method = "glog10", pseudocount = 1)))
    expect_gt(max(abs(out$mat_pre_scale -
                      transform_metab(f$mat, method = "log2", pseudocount = 1))), 1)
})

test_that("mod_met_corrected: a target that declares no transform warns instead of guessing", {
    f   <- met_fixture()
    cfg <- met_cfg("pqn", transform = "glog10")
    # (PQN's own "no is_QC column" warning is not the one under test)
    pqn <- suppressWarnings(mod_met_normalize_linear(f, "pqn", cfg))
    pqn$transform <- NULL   # a normalisation written before this contract
    expect_warning(
        suppressMessages(mod_met_corrected(
            norm_tss = NULL, norm_median = NULL, norm_pqn = pqn,
            logged = mod_met_log(f, cfg), meta = f$meta,
            out_dir = withr::local_tempdir(), config = cfg)),
        "does not declare the transform it applied")
})

# ---- 2. the contract, so a NEW normalisation cannot silently mislabel -------

# chosen_norm -> the target mod_met_corrected() will read for it. Enumerated
# here so that a normalisation added to the switch() in mod_met_corrected()
# without a line in this list fails the next test, instead of shipping a lane
# whose recorded transform is the configured one by default. (That is how
# `bio_factor` and `none` arrived: both were added to the switch long after
# this contract was written.)
norm_targets <- list(
    none           = function(f, cfg) mod_met_log(f, cfg),
    tss            = function(f, cfg) mod_met_normalize_linear(f, "tss", cfg),
    median         = function(f, cfg) mod_met_normalize_log(mod_met_log(f, cfg), cfg),
    pqn            = function(f, cfg) mod_met_normalize_linear(f, "pqn", cfg),
    # the two EigenMS targets need the EigenMS package; the body scan below
    # covers them without it.
    eigenms        = NULL,
    eigenms_forced = NULL,
    bio_factor     = function(f, cfg) mod_met_normalize_bio_factor(f, cfg)
)

# The named branches of `switch(chosen_norm, ...)` inside mod_met_corrected(),
# read off the parsed function body so the list above cannot drift from it.
switch_branches <- function(fn, var) {
    found <- NULL
    walk <- function(e) {
        if (!is.null(found) || !is.call(e)) return(invisible(NULL))
        if (identical(e[[1]], as.name("switch")) &&
            length(e) >= 2L && identical(e[[2]], as.name(var))) {
            found <<- e
            return(invisible(NULL))
        }
        parts <- as.list(e)
        for (i in seq_along(parts)) {
            # a call can hold the empty symbol (an omitted argument); touching it
            # errors, so ask tryCatch rather than testing it directly.
            is_sub <- tryCatch(is.call(parts[[i]]), error = function(err) FALSE)
            if (isTRUE(is_sub)) walk(parts[[i]])
        }
        invisible(NULL)
    }
    walk(body(fn))
    if (is.null(found)) stop("no switch(", var, ") found in the function body")
    nms <- names(found)[-(1:2)]
    nms[nzchar(nms)]
}

test_that("every chosen_norm branch of mod_met_corrected is covered by this file", {
    expect_setequal(switch_branches(mod_met_corrected, "chosen_norm"),
                    names(norm_targets))
})

test_that("every buildable normalisation target declares the transform it applied", {
    f <- met_fixture()
    for (nm in names(Filter(Negate(is.null), norm_targets))) {
        tgt <- suppressWarnings(suppressMessages(
            norm_targets[[nm]](f, met_cfg(nm, transform = "glog10"))))
        expect_false(is.null(tgt$transform),
                     info = paste0("chosen_norm = '", nm,
                                   "': its target must return a `transform` field"))
    }
})

test_that("every mod_met_normalize_* target declares the transform it applied", {
    # mod_met_corrected() reads the transform off the chosen target rather than
    # from a list of normalisation names, so this is the one place a new
    # normalisation has to be correct. Adding one without a `transform` field
    # fails here rather than shipping mislabelled fold changes.
    fns <- ls(envir = globalenv(), pattern = "^mod_met_(normalize_|log$)")
    expect_gt(length(fns), 0)
    for (fn in fns) {
        src <- paste(deparse(body(get(fn, envir = globalenv()))), collapse = "\n")
        expect_match(src, "transform(_in)?\\s*=", all = FALSE,
                     info = paste0(fn, "() must return the transform its `mat` carries"))
    }
})

# ---- 3. the consequence: the fold-change unit conversion -------------------

test_that("PQN + configured glog10: logFC is NOT inflated by log2(10)", {
    # The slope check. A correlation gate is blind here -- multiplying every
    # logFC by one constant keeps r = 1.0 -- so assert the slope itself:
    # median(reported logFC / the group difference on the tested matrix) must be
    # 1, and specifically not log2(10) = 3.3219.
    f   <- met_fixture()
    cfg <- met_cfg("pqn", transform = "glog10", scaling = "auto")
    out <- corrected_for(f, cfg)

    pre <- list(expr_raw = f$mat, expr_filt = f$mat,
                expr_log = mod_met_log(f, cfg)$mat,
                expr_pre_scale = out$mat_pre_scale,
                expr_work = out$mat,
                meta = out$meta, row_data = NULL, info = out$info)

    ctr <- data.frame(Contrast_name = "B_vs_A", Numerator = "B",
                      Denominator = "A", stringsAsFactors = FALSE)
    de <- suppressMessages(
        run_metabolomics_de(pre, cfg, ctr))$de_tables[["B_vs_A"]]

    tested <- pre$expr_pre_scale
    idx_A <- which(out$meta$Group == "A")
    idx_B <- which(out$meta$Group == "B")
    implied <- rowMeans(tested[de$feature_id, idx_B, drop = FALSE]) -
               rowMeans(tested[de$feature_id, idx_A, drop = FALSE])

    slope <- stats::median(de$logFC / implied)
    expect_equal(slope, 1, tolerance = 1e-8)
    expect_false(isTRUE(all.equal(slope, log2(10), tolerance = 1e-3)))

    # the correlation gate that cannot see this, stated so nobody reaches for it
    expect_equal(unname(stats::cor(de$logFC, implied)), 1, tolerance = 1e-8)
})
