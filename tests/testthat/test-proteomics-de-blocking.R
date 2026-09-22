# Tests for de$block_col / duplicateCorrelation blocking in the limma path
# test-proteomics-de-blocking.R

make_block_meta <- function() {
    data.frame(
        SampleName = c("L1_ctl", "L1_trt", "L2_ctl", "L2_trt", "L3_ctl", "L3_trt"),
        Group      = rep(c("ctl", "trt"), times = 3),
        Line       = rep(c("L1", "L2", "L3"), each = 2),
        stringsAsFactors = FALSE
    )
}

block_cfg <- function(block_col = "Line") {
    list(de = list(block_col = block_col))
}

test_that("resolve_de_block returns NULL when no block_col is configured", {
    expect_null(resolve_de_block(make_block_meta(), list(de = list()), "SampleName"))
    expect_null(resolve_de_block(make_block_meta(), block_cfg(""), "SampleName"))
})

test_that("resolve_de_block returns a factor aligned to sample order", {
    blk <- resolve_de_block(make_block_meta(), block_cfg(), "SampleName")
    expect_s3_class(blk, "factor")
    expect_length(blk, 6)
    expect_equal(nlevels(blk), 3)
    # order must follow the metadata rows, not the sorted level order
    expect_equal(as.character(blk), rep(c("L1", "L2", "L3"), each = 2))
})

test_that("resolve_de_block rejects a missing column", {
    expect_error(
        resolve_de_block(make_block_meta(), block_cfg("Donor"), "SampleName"),
        "not in the sample metadata"
    )
})

test_that("resolve_de_block rejects NA block labels", {
    meta <- make_block_meta()
    meta$Line[2] <- NA
    expect_error(
        resolve_de_block(meta, block_cfg(), "SampleName"),
        "missing values"
    )
})

test_that("resolve_de_block rejects a single-level block", {
    meta <- make_block_meta()
    meta$Line <- "only_one"
    expect_error(
        resolve_de_block(meta, block_cfg(), "SampleName"),
        "single level"
    )
})

test_that("resolve_de_block rejects a block with one sample per level", {
    expect_error(
        resolve_de_block(make_block_meta(), block_cfg("SampleName"), "SampleName"),
        "no two samples share a block"
    )
})

test_that("blocking recovers a paired effect that the unpaired fit misses", {
    skip_if_not_installed("limma")

    # Six samples, three blocks. A large per-block offset swamps a small but
    # perfectly consistent treatment effect, which is what blocking is for.
    withr::with_seed(42, {
        n_feat <- 300
        meta <- make_block_meta()
        block_offset <- rep(c(0, 4, 8), each = 2)
        treat_effect <- rep(c(0, 1), times = 3)

        expr <- matrix(rnorm(n_feat * 6, sd = 0.2), nrow = n_feat)
        expr <- sweep(expr, 2, block_offset, "+")
        expr[1:50, ] <- expr[1:50, ] + matrix(rep(treat_effect, each = 50), nrow = 50)
        rownames(expr) <- paste0("PROT", seq_len(n_feat))
        colnames(expr) <- meta$SampleName
    })

    design <- stats::model.matrix(~ 0 + factor(meta$Group))
    colnames(design) <- c("ctl", "trt")
    cm <- limma::makeContrasts(trt - ctl, levels = design)

    unpaired <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr, design), cm))
    blk <- resolve_de_block(meta, block_cfg(), "SampleName")
    dup <- limma::duplicateCorrelation(expr, design, block = blk)
    blocked <- limma::eBayes(limma::contrasts.fit(
        limma::lmFit(expr, design, block = blk, correlation = dup$consensus), cm))

    n_sig <- function(fit) sum(limma::topTable(fit, number = Inf)$adj.P.Val <= 0.05)

    expect_gt(dup$consensus, 0.5)
    expect_gt(n_sig(blocked), n_sig(unpaired))
    expect_equal(unique(unpaired$df.residual), unique(blocked$df.residual))
})
