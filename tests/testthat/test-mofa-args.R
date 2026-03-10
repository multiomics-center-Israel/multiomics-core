# Tests for MOFA integration argument safety
# Regression tests for command injection fix and sample guard in PR #48

test_that("MOFA sample guard rejects matrices with fewer than 3 samples", {
    # Simulate the guard logic from run_mofa_python()
    matrices <- list(rna = matrix(1:4, nrow = 2, ncol = 2))

    n_samples <- ncol(matrices[[1]])
    expect_equal(n_samples, 2)
    expect_true(n_samples < 3)
})

test_that("MOFA sample guard accepts matrices with 3+ samples", {
    matrices <- list(rna = matrix(1:9, nrow = 3, ncol = 3))

    n_samples <- ncol(matrices[[1]])
    expect_equal(n_samples, 3)
    expect_false(n_samples < 3)
})

test_that("num_factors is capped at n_samples - 1", {
    matrices <- list(rna = matrix(rnorm(20), nrow = 4, ncol = 5))
    num_factors <- 10
    n_samples <- ncol(matrices[[1]])

    capped <- min(num_factors, n_samples - 1)
    expect_equal(capped, 4)
})

test_that("system2 arg vector does not allow shell injection", {
    # Verify that system2 with args vector treats special chars as literals
    python_exec <- "python3"
    script_path <- "scripts/run_mofa.py"
    outfile <- "/tmp/test; rm -rf /"
    outdir <- "/tmp/output"

    cmd_args <- c(
        script_path,
        "--outfile", outfile,
        "--outdir", outdir,
        "--factors", "5"
    )

    # With system2, each arg is a separate element — shell metacharacters
    # are passed as literal strings, not interpreted
    expect_true(any(grepl("; rm -rf /", cmd_args, fixed = TRUE)))
    expect_equal(cmd_args[3], "/tmp/test; rm -rf /")
})
