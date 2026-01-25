test_that("normalization returns correct dimensions and attributes", {
    # Helper to load project source if not a package
    # We assume tests are run from project root or we find it
    # For now, let's just source the specific file needed, or rely on a helper (see next step)

    # Logic: If we are running checks, we want to ensure definitions are available.
    # Let's try sourcing explicitly for this unit test first.
    if (!exists("normalize_counts")) {
        source(file.path("../../R/preprocess/normalize.R"))
    }

    # 1. Setup Data
    counts <- matrix(rpois(30, lambda = 50), nrow = 5, ncol = 6)
    meta <- data.frame(group = rep(c("A", "B"), each = 3))
    rownames(counts) <- paste0("gene", 1:5)
    colnames(counts) <- paste0("S", 1:6)

    # Ensure meta has matching rownames if your function requires alignment by name
    rownames(meta) <- colnames(counts)

    # 2. Test TMM (mocking or ensuring method exists)
    # verify method argument support
    # Assuming "TMMlogCPM" is a valid method in your codebase

    # We wrap in try because we haven't seen the full dependency chain of normalize_counts
    # (e.g. does it depend on 'edgeR' attached?)
    suppressPackageStartupMessages(library(edgeR))
    suppressPackageStartupMessages(library(limma))

    res_tmm <- normalize_counts(counts, method = "TMMlogCPM")

    expect_equal(dim(res_tmm), c(5, 6))
    expect_equal(attr(res_tmm, "method"), "TMMlogCPM")
    expect_true(is.numeric(res_tmm))

    # 3. Test VST
    # check dependency for VST (DESeq2 usually)
    if (requireNamespace("DESeq2", quietly = TRUE)) {
        res_vst <- normalize_counts(counts, meta = meta, method = "VST")
        expect_equal(dim(res_vst), c(5, 6))
        expect_equal(attr(res_vst, "method"), "VST")
    }
})
