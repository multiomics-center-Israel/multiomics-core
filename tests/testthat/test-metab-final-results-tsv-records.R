# Regression: embedded record separators in chemical names must not turn one
# feature into multiple rows in Final_results_metabolomics.tsv.

test_that("Final_results TSV retains one physical record per feature", {
    df <- data.frame(
        feature_id = c("F1", "F2"),
        original_id = c("chemical\ncontinuation", "simple"),
        Name = c("name\r\ncontinuation", NA_character_),
        "Unique ID" = c("part1\tpart2", "other"),
        intensity = c(25.125, NA_real_),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    output_dir <- tempfile("metab_final_tsv_")
    on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
    path <- write_metabolomics_final_results_tsv(df, output_dir)

    # Original annotations passed to the Excel writer remain unchanged.
    expect_identical(df$original_id[1], "chemical\ncontinuation")
    expect_identical(df$Name[1], "name\r\ncontinuation")
    expect_identical(df[["Unique ID"]][1], "part1\tpart2")

    # Header + one physical line per feature, with no parsing failures.
    expect_length(readLines(path, warn = FALSE), nrow(df) + 1L)
    result <- readr::read_tsv(path, show_col_types = FALSE)
    expect_equal(nrow(readr::problems(result)), 0L)
    expect_equal(nrow(result), nrow(df))
    expect_identical(as.character(result$feature_id), df$feature_id)
    expect_identical(result$original_id[1], "chemical continuation")
    expect_identical(result$Name[1], "name continuation")
    expect_identical(result[["Unique ID"]][1], "part1 part2")
    expect_equal(result$intensity, df$intensity)
})
