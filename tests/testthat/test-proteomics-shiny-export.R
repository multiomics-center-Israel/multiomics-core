# tests/testthat/test-proteomics-shiny-export.R
#
# build_shiny_payload_proteomics(): de_final_table is a clean DE data.frame built
# from the final-results data.frame (same pattern as the RNA-seq and
# metabolomics payloads), not read back from the Excel workbook, whose Results
# sheet is a presentation layout. Fixtures are synthetic and deterministic.

make_prot_export_fixture <- function(pass = c(1, 1, 0, 0, 1)) {
    samples <- paste0("s", 1:4)
    feats   <- paste0("P", 1:5)

    mat <- matrix(as.numeric(1:20), nrow = 5, ncol = 4,
                  dimnames = list(feats, samples))
    meta <- data.frame(sample_id = samples, grp = c("A", "A", "B", "B"),
                       row.names = samples, stringsAsFactors = FALSE)
    pre <- list(expr_filt = mat, expr_imp_single = mat, meta = meta,
                row_data = data.frame(protein_id = feats,
                                      stringsAsFactors = FALSE))

    # The ID column is deliberately NOT the first column and is not the default
    # ("FeatureID"), so the test fails if the builder guesses by position or
    # ignores modes$proteomics$de_table$id_col.
    final_results <- data.frame(
        Genes = paste0("G", 1:5),
        Protein = feats,
        s1 = 1:5 + 0.5, s2 = 1:5 + 1.5, s3 = 1:5 + 2.5, s4 = 1:5 + 3.5,
        linearFC.imputs.A_vs_B = c(2, -2, 1.1, -1.2, 3),
        pvalue.imputs.A_vs_B   = c(0.001, 0.001, 0.5, 0.6, 0.002),
        padj.imputs.A_vs_B     = c(0.01, 0.01, 0.6, 0.7, 0.02),
        upDown.imputs.A_vs_B   = c("up", "down", "", "", "up"),
        manual_cutoffs.imputs.A_vs_B = NA,
        pass_any_contrast = pass,
        stringsAsFactors = FALSE
    )

    config <- list(modes = list(proteomics = list(
        de            = list(padj_cutoff = 0.05, linear_fc_cutoff = 1.5),
        normalization = list(method = "none"),
        effects       = list(color = list("grp")),
        id_columns    = list(protein_id = "protein_id"),
        de_table      = list(id_col = "Protein")
    )))

    zscore_mat <- matrix(0, nrow = 3, ncol = 4,
                         dimnames = list(c("P1", "P2", "P5"),
                                         paste0(samples, ".zscore")))
    clustering_res <- list(
        excel_order = list(ordered_ids = c("P5", "P1", "P2"),
                           zscore_mat  = zscore_mat)
    )

    list(pre = pre, inputs = list(contrasts = NULL), config = config,
         final_results = final_results, clustering_res = clustering_res)
}

# de_res is an empty (non-NULL) list on purpose: the de_final_table block only
# needs de_res to be present, and this keeps the DE-statistics code out of scope.
build_prot_payload <- function(fx, final_results = fx$final_results,
                               clustering_res = fx$clustering_res, ...) {
    suppressMessages(suppressWarnings(
        build_shiny_payload_proteomics(
            pre = fx$pre, de_res = list(), inputs = fx$inputs,
            config = fx$config, clustering_res = clustering_res,
            final_results = final_results, ...
        )
    ))
}

test_that("de_final_table holds only DE rows, without the cutoff/pass helper columns", {
    de <- build_prot_payload(make_prot_export_fixture())$de_final_table

    expect_s3_class(de, "data.frame")
    expect_setequal(de$Protein, c("P1", "P2", "P5"))
    expect_false("pass_any_contrast" %in% names(de))
    expect_false(any(startsWith(names(de), "manual_cutoffs")))
    expect_true(all(c("linearFC.imputs.A_vs_B", "pvalue.imputs.A_vs_B",
                      "padj.imputs.A_vs_B", "upDown.imputs.A_vs_B") %in% names(de)))
})

test_that("de_final_table is a clean table: numeric columns stay numeric, no layout artefacts", {
    de <- build_prot_payload(make_prot_export_fixture())$de_final_table

    expect_true(is.numeric(de$s1))
    expect_true(is.numeric(de$linearFC.imputs.A_vs_B))
    expect_false(any(grepl("^\\.\\.\\.[0-9]+$", names(de))))
    expect_false("Sample_ID" %in% de$Protein)     # no sheet-layout rows above the data
})

test_that("de_final_table adds clustering order and z-scores, matched by the configured id column", {
    de <- build_prot_payload(make_prot_export_fixture())$de_final_table
    ord <- setNames(de$order, de$Protein)

    expect_identical(ord[["P5"]], 1L)
    expect_identical(ord[["P1"]], 2L)
    expect_identical(ord[["P2"]], 3L)
    expect_true(all(paste0("s", 1:4, ".zscore") %in% names(de)))
})

test_that("de_final_table has no order or z-scores without a clustering order", {
    de <- build_prot_payload(make_prot_export_fixture(), clustering_res = NULL)$de_final_table

    expect_setequal(de$Protein, c("P1", "P2", "P5"))
    expect_false("order" %in% names(de))
    expect_false(any(grepl("\\.zscore$", names(de))))
})

test_that("de_final_table is a zero-row data.frame when no feature passes", {
    de <- build_prot_payload(make_prot_export_fixture(pass = rep(0, 5)))$de_final_table

    expect_s3_class(de, "data.frame")
    expect_identical(nrow(de), 0L)
})

test_that("de_final_table is NULL, with the key present, when no final_results is passed", {
    p <- build_prot_payload(make_prot_export_fixture(), final_results = NULL)

    expect_true("de_final_table" %in% names(p))
    expect_null(p$de_final_table)
})

test_that("a final_results without the configured id column warns and leaves de_final_table NULL", {
    fx <- make_prot_export_fixture()
    fx$config$modes$proteomics$de_table$id_col <- "NoSuchColumn"

    expect_warning(
        p <- suppressMessages(build_shiny_payload_proteomics(
            pre = fx$pre, de_res = list(), inputs = fx$inputs, config = fx$config,
            clustering_res = fx$clustering_res, final_results = fx$final_results
        )),
        "de_final_table"
    )
    expect_null(p$de_final_table)
})

test_that("de_final_xlsx is still attached from xlsx_files when no final_results is passed", {
    f <- file.path(tempdir(), "Final_results_DE_P_0.05.xlsx")
    bytes <- as.raw(1:10)          # content is irrelevant: the bytes are stored, not parsed
    writeBin(bytes, f)
    on.exit(unlink(f), add = TRUE)

    p <- build_prot_payload(make_prot_export_fixture(), final_results = NULL,
                            xlsx_files = f)

    expect_identical(p$de_final_xlsx, bytes)
    expect_null(p$de_final_table)
})
