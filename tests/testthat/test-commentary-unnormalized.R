# tests/testthat/test-commentary-unnormalized.R
#
# Commentary for values the pipeline did not normalize (RNA preprocessed input,
# or norm_method "none" in the other omics) must not claim that normalization
# removed technical variation.

test_that("commentary_values_unnormalized() recognises the 'none' methods", {
    expect_true(commentary_values_unnormalized(list(norm_method = "none")))
    expect_true(commentary_values_unnormalized(list(norm_method = "none (preprocessed input)")))
    expect_false(commentary_values_unnormalized(list(norm_method = "VST")))
    expect_false(commentary_values_unnormalized(list(norm_method = "TMMlogCPM")))
    expect_false(commentary_values_unnormalized(list()))
    expect_false(commentary_values_unnormalized(list(norm_method = NA_character_)))
})

fallback_text <- function(plot_type, figure_id, norm_method) {
    fig <- data.frame(figure_id = figure_id, title = "T", plot_type = plot_type,
                      stringsAsFactors = FALSE)
    cm <- generate_fallback_commentary(fig, list(norm_method = norm_method))
    paste(unlist(cm[c("what_is_this", "biological_interpretation",
                      "observations", "issues_checks")]), collapse = " ")
}

test_that("fallback commentary does not claim normalization for unnormalized values", {
    for (pt in c("boxplot", "density")) {
        fid <- if (pt == "boxplot") "norm_boxplot" else "density_plot"
        txt <- fallback_text(pt, fid, "none (preprocessed input)")
        expect_match(txt, "did not normalize", fixed = TRUE)
        expect_no_match(txt, "After normalization", fixed = TRUE)
        expect_no_match(txt, "adequately removed", fixed = TRUE)
    }
    expect_match(fallback_text("boxplot", "norm_boxplot", "VST"),
                 "After normalization (VST)", fixed = TRUE)
})

test_that("AI guidance tells the model not to claim normalization for unnormalized values", {
    ctx <- list(plot_type = "boxplot", figure_id = "norm_boxplot", norm_method = "none")
    expect_match(build_figure_guidance(ctx), "did NOT normalize", fixed = TRUE)
    ctx$norm_method <- "VST"
    expect_match(build_figure_guidance(ctx), "NORMALIZED EXPRESSION BOXPLOT", fixed = TRUE)
})
