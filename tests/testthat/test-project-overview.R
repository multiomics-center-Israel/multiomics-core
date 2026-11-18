# tests/testthat/test-project-overview.R
#
# project.user_notes is documented in the config templates as "appears in
# report as Project Overview". project_notes_html() is the report-side
# reading; the RNA report shows the section only when there are notes.
#
# All text here is synthetic.

test_that("no notes means no section", {
    expect_null(project_notes_html(NULL))
    expect_null(project_notes_html(""))
    expect_null(project_notes_html("   \n  "))
    expect_null(project_notes_html(NA_character_))
    expect_null(project_notes_html(character(0)))
    expect_null(project_notes_html(list("not a string")))
})

test_that("a multi-line note keeps its line breaks", {
    html <- project_notes_html("First line.\nSecond line.\n\nAfter a blank line.\n")
    expect_match(html, "First line.<br>\nSecond line.<br>\n<br>\nAfter a blank line.",
                 fixed = TRUE)
    # Trailing newline from a YAML `|` block is not shown as an empty line.
    expect_false(grepl("blank line.<br>", html, fixed = TRUE))
    expect_match(html, "^<div class=\"project-overview-box\"")
})

test_that("Windows line endings do not leave stray carriage returns", {
    html <- project_notes_html("One\r\nTwo")
    expect_false(grepl("\r", html, fixed = TRUE))
    expect_match(html, "One<br>\nTwo", fixed = TRUE)
})

test_that("notes are shown as text, not as HTML", {
    html <- project_notes_html("Compare <A> & <B>")
    expect_match(html, "Compare &lt;A&gt; &amp; &lt;B&gt;", fixed = TRUE)
})

test_that("a YAML list of lines is one note", {
    html <- project_notes_html(c("Line one", NA, "Line two"))
    expect_match(html, "Line one<br>\nLine two", fixed = TRUE)
})

test_that("the RNA report renders the overview before the Technical Summary", {
    rmd <- file.path(root_dir, "R", "domain", "rnaseq", "report_template.Rmd")
    skip_if_not(file.exists(rmd), "report_template.Rmd not found")
    src <- readLines(rmd, warn = FALSE)

    lines <- find_rmd_chunk_lines(rmd, "project-overview")
    expect_false(is.null(lines))
    chunk <- src[seq(lines[1], lines[2])]
    expect_true(any(grepl("project_notes_html(cfg$project$user_notes)", chunk, fixed = TRUE)))
    # The heading is written only inside the NULL check, never unconditionally.
    expect_false("# Project Overview" %in% src)
    expect_lt(lines[2], which(src == "# Technical Summary"))
})
