# tests/testthat/test-annotate-de-stats.R
#
# Unit tests for annotate_de_stats(): left-join of row_data annotation
# columns into a payload$de_stats data.frame.

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

make_de <- function(ids   = c("f1", "f2", "f3"),
                    extra = list(linearFC.A       = c(1.5, -2.0, 0.5),
                                 pass_any_contrast = c(1, NA, 1))) {
    df <- data.frame(feature_id = ids, stringsAsFactors = FALSE)
    for (nm in names(extra)) df[[nm]] <- extra[[nm]]
    df
}

make_rd <- function(ids   = c("f1", "f2", "f3"),
                    extra = list(Name    = c("alpha", "beta", "gamma"),
                                 Pathway = c("p1",    "p2",   "p1"))) {
    df <- data.frame(feature_id = ids, stringsAsFactors = FALSE)
    for (nm in names(extra)) df[[nm]] <- extra[[nm]]
    df
}

# ---------------------------------------------------------------------------
# Basic merge
# ---------------------------------------------------------------------------

test_that("merge preserves row count", {
    de  <- make_de()
    rd  <- make_rd()
    out <- annotate_de_stats(de, rd)
    expect_equal(nrow(out), nrow(de))
})

test_that("annotation columns are merged with correct values", {
    de  <- make_de()
    rd  <- make_rd()
    out <- annotate_de_stats(de, rd)
    expect_true(all(c("Name", "Pathway") %in% colnames(out)))
    expect_equal(out$Name,    c("alpha", "beta", "gamma"))
    expect_equal(out$Pathway, c("p1",    "p2",   "p1"))
})

test_that("column order is id, annotations, stats", {
    de  <- make_de()
    rd  <- make_rd()
    out <- annotate_de_stats(de, rd)
    expect_equal(
        colnames(out),
        c("feature_id", "Name", "Pathway", "linearFC.A", "pass_any_contrast")
    )
})

test_that("original DE statistic values are preserved", {
    de  <- make_de()
    rd  <- make_rd()
    out <- annotate_de_stats(de, rd)
    expect_equal(out$linearFC.A,        de$linearFC.A)
    expect_equal(out$pass_any_contrast, de$pass_any_contrast)
    expect_equal(out$feature_id,        de$feature_id)
})

test_that("merge works when join keys differ in order between de_stats and row_data", {
    de <- make_de(ids   = c("f3", "f1", "f2"),
                  extra = list(linearFC.A = c(0.5, 1.5, -2.0)))
    rd <- make_rd()  # ids in order f1, f2, f3
    out <- annotate_de_stats(de, rd)
    expect_equal(out$Name, c("gamma", "alpha", "beta"))
})

# ---------------------------------------------------------------------------
# No-op paths
# ---------------------------------------------------------------------------

test_that("NULL row_data returns de_stats unchanged", {
    de <- make_de()
    expect_message(out <- annotate_de_stats(de, NULL), "NULL or empty")
    expect_identical(out, de)
})

test_that("zero-row row_data returns de_stats unchanged", {
    de <- make_de()
    rd <- make_rd()[0, , drop = FALSE]
    expect_message(out <- annotate_de_stats(de, rd), "NULL or empty")
    expect_identical(out, de)
})

test_that("row_data with only the join key is a no-op", {
    de <- make_de()
    rd <- data.frame(feature_id = c("f1", "f2", "f3"),
                     stringsAsFactors = FALSE)
    expect_message(out <- annotate_de_stats(de, rd), "no annotation columns")
    expect_identical(out, de)
})

test_that("NULL de_stats returns NULL", {
    expect_null(annotate_de_stats(NULL, make_rd()))
})

# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

test_that("missing id_col_de errors", {
    expect_error(
        annotate_de_stats(make_de(), make_rd(), id_col_de = "not_a_col"),
        "id_col_de"
    )
})

test_that("missing id_col_row errors", {
    expect_error(
        annotate_de_stats(make_de(), make_rd(), id_col_row = "not_a_col"),
        "id_col_row"
    )
})

test_that("non-data.frame de_stats errors", {
    expect_error(annotate_de_stats(list(a = 1), make_rd()), "data.frame")
})

test_that("non-data.frame row_data errors", {
    expect_error(annotate_de_stats(make_de(), list(a = 1)), "data.frame")
})

# ---------------------------------------------------------------------------
# Collisions
# ---------------------------------------------------------------------------

test_that("collision skip (default) preserves the existing column", {
    de       <- make_de()
    de$Name  <- c("DE-alpha", "DE-beta", "DE-gamma")  # collides with row_data
    rd       <- make_rd()
    expect_message(out <- annotate_de_stats(de, rd), "skipping")
    expect_equal(out$Name,    c("DE-alpha", "DE-beta", "DE-gamma"))
    expect_equal(out$Pathway, c("p1", "p2", "p1"))  # non-colliding still merged
})

test_that("collision overwrite replaces the existing column from row_data", {
    de      <- make_de()
    de$Name <- c("DE-alpha", "DE-beta", "DE-gamma")
    rd      <- make_rd()
    out     <- annotate_de_stats(de, rd, on_collision = "overwrite")
    expect_equal(out$Name, c("alpha", "beta", "gamma"))
})

test_that("all-colliding skip returns de_stats unchanged", {
    de         <- make_de()
    de$Name    <- c("a", "b", "c")
    de$Pathway <- c("x", "y", "z")
    rd         <- make_rd()
    expect_message(out <- annotate_de_stats(de, rd), "skipping")
    expect_identical(out, de)
})

# ---------------------------------------------------------------------------
# Unmatched features / NA in join key
# ---------------------------------------------------------------------------

test_that("unmatched features become NA in annotation columns", {
    de <- make_de(
        ids   = c("f1", "f2", "f3", "f4"),
        extra = list(linearFC.A        = c(1, 2, 3, 4),
                     pass_any_contrast = c(1, NA, 1, NA))
    )
    rd  <- make_rd()  # only f1, f2, f3
    suppressWarnings(out <- annotate_de_stats(de, rd))
    expect_equal(out$Name, c("alpha", "beta", "gamma", NA))
    expect_equal(nrow(out), 4L)
})

test_that("low match rate emits a warning", {
    de <- make_de(ids   = paste0("f", 1:10),
                  extra = list(pass_any_contrast = rep(1, 10)))
    rd <- make_rd()  # matches 3 of 10
    expect_warning(annotate_de_stats(de, rd), "low match rate")
})

test_that("NA in de_stats join key yields NA in annotations", {
    de  <- make_de(ids = c("f1", NA, "f3"))
    rd  <- make_rd()
    out <- annotate_de_stats(de, rd)
    expect_true(is.na(out$Name[2]))
    expect_equal(out$Name[c(1, 3)], c("alpha", "gamma"))
})

# ---------------------------------------------------------------------------
# Duplicate keys in row_data
# ---------------------------------------------------------------------------

test_that("duplicate keys in row_data warn and first match wins", {
    de <- make_de()
    rd <- rbind(
        make_rd(),
        data.frame(feature_id = "f1", Name = "DUPLICATE", Pathway = "px",
                   stringsAsFactors = FALSE)
    )
    expect_warning(out <- annotate_de_stats(de, rd), "duplicate key")
    expect_equal(out$Name[1], "alpha")
})

# ---------------------------------------------------------------------------
# Custom id column names (RNA-seq / proteomics scenarios)
# ---------------------------------------------------------------------------

test_that("merge works with different id column names on each side", {
    de <- data.frame(
        feature_id        = c("g1", "g2"),
        linearFC.A        = c(0.5, -0.5),
        pass_any_contrast = c(1, NA),
        stringsAsFactors  = FALSE
    )
    rd <- data.frame(
        gene_id          = c("g1", "g2"),
        Symbol           = c("GENE1", "GENE2"),
        stringsAsFactors = FALSE
    )
    out <- annotate_de_stats(de, rd,
                             id_col_de  = "feature_id",
                             id_col_row = "gene_id")
    expect_equal(colnames(out),
                 c("feature_id", "Symbol", "linearFC.A", "pass_any_contrast"))
    expect_equal(out$Symbol, c("GENE1", "GENE2"))
})
