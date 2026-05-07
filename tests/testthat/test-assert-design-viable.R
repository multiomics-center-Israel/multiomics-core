test_that("resolve_group_col walks the canonical chain", {
    expect_equal(
        resolve_group_col(list(filtering = list(group_col = "treatment"),
                               de_table  = list(group_col = "ignored"),
                               effects   = list(color = "ignored2"))),
        "treatment"
    )
    expect_equal(
        resolve_group_col(list(de_table = list(group_col = "Condition"),
                               effects  = list(color = "ignored"))),
        "Condition"
    )
    expect_equal(
        resolve_group_col(list(effects = list(color = c("Group", "batch")))),
        "Group"
    )
    expect_equal(resolve_group_col(list()), "Condition")
})

test_that("resolve_group_col is strict — empty/NA fall through", {
    cfg <- list(filtering = list(group_col = ""),
                de_table  = list(group_col = NA_character_),
                effects   = list(color = "treatment"))
    expect_equal(resolve_group_col(cfg), "treatment")
})

test_that("resolve_group_col uses only the first element of effects$color", {
    cfg <- list(effects = list(color = c("primary", "secondary", "tertiary")))
    expect_equal(resolve_group_col(cfg), "primary")
})

test_that("assert_design_viable stops on group with n=1", {
    meta <- data.frame(
        sample    = paste0("S", 1:5),
        Condition = c("ctrl", "ctrl", "ctrl", "trt", "trt")
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    meta1 <- data.frame(sample = "S1", Condition = "ctrl")
    expect_error(assert_design_viable(meta1, cfg), "n >= 2 per group")
})

test_that("assert_design_viable stops naming all underpowered groups when multiple", {
    meta <- data.frame(
        sample    = paste0("S", 1:2),
        Condition = c("ctrl", "trt")
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_error(assert_design_viable(meta, cfg), "Groups with n < 2")
})

test_that("assert_design_viable warns severe at n=2 (single composite warning)", {
    meta <- data.frame(
        sample    = paste0("S", 1:4),
        Condition = c("ctrl", "ctrl", "trt", "trt")
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_warning(assert_design_viable(meta, cfg), "n=2: severe")
})

test_that("assert_design_viable warns low-power at 3 <= n < 5", {
    meta <- data.frame(
        sample    = paste0("S", 1:6),
        Condition = rep(c("ctrl", "trt"), each = 3)
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_warning(assert_design_viable(meta, cfg), "low statistical power")
})

test_that("assert_design_viable warns on global imbalance separately", {
    meta <- data.frame(
        sample    = paste0("S", 1:23),
        Condition = c(rep("ctrl", 20), rep("trt", 3))
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_warning(assert_design_viable(meta, cfg), "imbalanced")
})

test_that("assert_design_viable is silent at n=5 per group (boundary)", {
    meta <- data.frame(
        sample    = paste0("S", 1:10),
        Condition = rep(c("ctrl", "trt"), each = 5)
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_silent(assert_design_viable(meta, cfg))
})

test_that("assert_design_viable is silent at n=5/5/5 three-group design", {
    meta <- data.frame(
        sample    = paste0("S", 1:15),
        Condition = rep(c("a", "b", "c"), each = 5)
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_silent(assert_design_viable(meta, cfg))
})

test_that("assert_design_viable stops when group column is missing from meta", {
    meta <- data.frame(sample = paste0("S", 1:6),
                       treatment = rep(c("ctrl", "trt"), each = 3))
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_error(
        assert_design_viable(meta, cfg),
        "not found in metadata"
    )
})

test_that("assert_design_viable stops on NA values in the group column", {
    meta <- data.frame(
        sample    = paste0("S", 1:6),
        Condition = c("ctrl", "ctrl", NA, "trt", "trt", "trt"),
        stringsAsFactors = FALSE
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_error(assert_design_viable(meta, cfg), "Assign a group label")
})

test_that("assert_design_viable stops on empty-string values in the group column", {
    meta <- data.frame(
        sample    = paste0("S", 1:6),
        Condition = c("ctrl", "ctrl", "", "trt", "trt", "trt"),
        stringsAsFactors = FALSE
    )
    cfg <- list(filtering = list(group_col = "Condition"),
                effects   = list(samples = "sample"))
    expect_error(assert_design_viable(meta, cfg), "Assign a group label")
})

# ---------------------------------------------------------------------------
# Integration tests: assertion fires from preprocess_rna at the post-
# sample_filter call site. These exercise the actual call site, not the
# pure function — a unit-level call would have passed even with the buggy
# pre-filter placement that codex caught.
# ---------------------------------------------------------------------------

# Fixture builder. The "include_singleton" group ('outlier' n=1) is the
# trigger for the design-viability gate; whether it survives or not depends
# on the sample_filter config passed to preprocess_rna.
make_rna_fixture <- function(include_singleton = TRUE) {
    base_samples <- c(paste0("C", 1:5), paste0("T", 1:5))
    base_cond <- c(rep("ctrl", 5), rep("trt", 5))
    if (include_singleton) {
        samples <- c(base_samples, "X1")
        conditions <- c(base_cond, "outlier")
    } else {
        samples <- base_samples
        conditions <- base_cond
    }
    n_genes <- 20
    set.seed(42)
    counts_mat <- matrix(rpois(n_genes * length(samples), lambda = 50),
                         nrow = n_genes, ncol = length(samples))
    counts_df <- data.frame(gene_id = paste0("G", seq_len(n_genes)),
                            counts_mat, check.names = FALSE,
                            stringsAsFactors = FALSE)
    colnames(counts_df)[-1] <- samples

    metadata <- data.frame(
        SampleID = samples, Condition = conditions,
        stringsAsFactors = FALSE
    )
    list(counts = counts_df, metadata = metadata)
}

# Config builder. filter_mode="none" short-circuits gene filtering so the
# test reaches normalization quickly; we don't care whether normalization
# completes — we only care whether [design] errors fire.
make_rna_config <- function(sample_filter_enabled = FALSE,
                            sample_filter_rules = NULL) {
    list(
        modes = list(
            rna = list(
                id_columns    = list(sample_col = "SampleID", gene_id = "gene_id"),
                effects       = list(samples = "SampleID"),
                filtering     = list(group_col = "Condition", mode = "none"),
                sample_filter = list(
                    enabled = sample_filter_enabled,
                    rules   = sample_filter_rules
                )
            )
        )
    )
}

# Run preprocess_rna and return either NULL (ran cleanly through assertion)
# or the error message string (so we can pattern-match on [design]).
capture_preprocess_error <- function(inputs, config) {
    tryCatch({
        suppressWarnings(suppressMessages(
            preprocess_rna(inputs, config, verbose = FALSE)
        ))
        NULL
    }, error = function(e) conditionMessage(e))
}

test_that("preprocess_rna fires [design] error on a singleton group with no sample_filter", {
    inputs <- make_rna_fixture(include_singleton = TRUE)
    config <- make_rna_config(sample_filter_enabled = FALSE)

    err <- capture_preprocess_error(inputs, config)
    expect_true(!is.null(err) && grepl("\\[design\\]", err),
                info = sprintf("expected a [design] error, got: %s", err %||% "(no error)"))
})

test_that("preprocess_rna does NOT fire [design] errors when sample_filter removes the singleton (regression for pre-filter placement bug)", {
    # Fixture has an n=1 'outlier' group that would trigger the design gate
    # if the gate ran on raw input metadata. sample_filter intentionally
    # whitelists only ctrl/trt, removing the singleton before the gate runs.
    inputs <- make_rna_fixture(include_singleton = TRUE)
    config <- make_rna_config(
        sample_filter_enabled = TRUE,
        sample_filter_rules   = list(Condition = c("ctrl", "trt"))
    )

    err <- capture_preprocess_error(inputs, config)
    # preprocess_rna may still fail later for reasons unrelated to design
    # viability (e.g. normalization on a tiny fixture); we only assert that
    # if it DID fail, the failure was not from the [design] gate.
    if (!is.null(err)) {
        expect_false(grepl("\\[design\\]", err),
                     info = sprintf("unexpected [design] error: %s", err))
    } else {
        succeed("preprocess_rna ran past the design-viability gate")
    }
})

test_that("preprocess_rna fires [design] error on the same fixture when sample_filter is disabled (negative twin — proves the filter is doing the work, not silent skipping)", {
    inputs <- make_rna_fixture(include_singleton = TRUE)
    # sample_filter rules present but enabled=FALSE → get_sample_filter_rules
    # returns NULL → apply_sample_filter is not called → singleton survives
    # to the gate.
    config <- make_rna_config(
        sample_filter_enabled = FALSE,
        sample_filter_rules   = list(Condition = c("ctrl", "trt"))
    )

    err <- capture_preprocess_error(inputs, config)
    expect_true(!is.null(err) && grepl("\\[design\\]", err),
                info = sprintf("expected [design] error when filter is disabled, got: %s",
                               err %||% "(no error)"))
})
