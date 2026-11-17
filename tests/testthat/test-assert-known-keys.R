# assert_known_keys(): the one check every config validator uses to refuse a key
# no code reads. Synthetic blocks only.

known <- c("p_cutoff", "method")

test_that("an absent block and a block of known keys pass", {
    expect_true(assert_known_keys(NULL, known, "de"))
    expect_true(assert_known_keys(list(), known, "de"))
    expect_true(assert_known_keys(list(p_cutoff = 0.05), known, "de"))
    expect_true(assert_known_keys(list(method = "limma", p_cutoff = 0.05), known, "de"))
})

test_that("an unknown key stops, naming the block, the key and the known keys", {
    expect_error(assert_known_keys(list(p_cuttoff = 0.05), known, "modes$x$de"),
                 "modes$x$de has unknown key(s): p_cuttoff. Known: p_cutoff, method.",
                 fixed = TRUE)
})

test_that("every unknown key is named, in the order written", {
    expect_error(assert_known_keys(list(foo = 1, method = "limma", bar = 2), known, "de"),
                 "unknown key(s): foo, bar.", fixed = TRUE)
})

test_that("keys are matched exactly, not by prefix or case", {
    # `$` on a list matches by prefix; a validator must not.
    expect_error(assert_known_keys(list(p_cut = 0.05), known, "de"),
                 "unknown key(s): p_cut.", fixed = TRUE)
    expect_error(assert_known_keys(list(Method = "limma"), known, "de"),
                 "unknown key(s): Method.", fixed = TRUE)
})

test_that("a block that is not a map of named settings stops", {
    msg <- "de must be a map of settings (key: value)."
    expect_error(assert_known_keys("limma", known, "de"), msg, fixed = TRUE)
    expect_error(assert_known_keys(list("limma", 0.05), known, "de"), msg, fixed = TRUE)
    # One named and one unnamed entry: the unnamed one has an empty name.
    expect_error(assert_known_keys(list(method = "limma", 0.05), known, "de"), msg,
                 fixed = TRUE)
    expect_error(assert_known_keys(stats::setNames(list(1), " "), known, "de"), msg,
                 fixed = TRUE)
})

test_that("a hint is appended to the unknown-key message only", {
    expect_error(assert_known_keys(list(foo = 1), known, "de", hint = "Remove the key."),
                 "Known: p_cutoff, method.\n  Remove the key.", fixed = TRUE)
    # Without a hint the message ends at the list of known keys.
    err <- tryCatch(assert_known_keys(list(foo = 1), known, "de"),
                    error = function(e) conditionMessage(e))
    expect_identical(err, "de has unknown key(s): foo. Known: p_cutoff, method.")
})
