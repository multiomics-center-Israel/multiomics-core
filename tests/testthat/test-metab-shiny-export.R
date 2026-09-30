# tests/testthat/test-metab-shiny-export.R
#
# build_shiny_payload_metabolomics(): de_final_table is a clean DE data.frame
# built from build_final_results_metabolomics() (same pattern as the RNA-seq
# payload), and clust_patterns_list is carried over from the clustering objects.
# All fixtures are synthetic and deterministic; nothing here is random.

make_metab_export_fixture <- function(pass = c(1L, 1L, 0L, 0L, 1L),
                                      contrast_label = "A_vs_B") {
    samples <- paste0("s", 1:6)
    feats   <- paste0("F", 1:5)

    mat <- matrix(as.numeric(1:30), nrow = 5, ncol = 6,
                  dimnames = list(feats, samples))
    meta <- data.frame(sample_id = samples, grp = rep(c("A", "B"), 3),
                       row.names = samples, stringsAsFactors = FALSE)
    # rownames(row_data) must be the feature ids: the final-results builder
    # matches annotation columns by rowname, as in the real pipeline.
    row_data <- data.frame(feature_id = feats,
                           original_id = paste0("orig", 1:5),
                           name = paste0("met", 1:5),
                           row.names = feats, stringsAsFactors = FALSE)
    pre <- list(expr_raw = mat, expr_filt = mat, expr_work = mat,
                meta = meta, row_data = row_data)

    # Contrast column names follow get_contrast_cols(mode = "metabolomics").
    summary_df <- data.frame(
        feature_id          = feats,
        linearFC.A_vs_B     = c(2, -2, 1.1, -1.2, 3),
        pvalue.A_vs_B       = c(0.001, 0.001, 0.5, 0.6, 0.002),
        padj.A_vs_B         = c(0.01, 0.01, 0.6, 0.7, 0.02),
        pass.A_vs_B         = pass,
        pass_any_contrast   = pass,
        stringsAsFactors    = FALSE
    )
    de_res <- list(
        summary_df = summary_df,
        de_tables  = setNames(list(data.frame(feature_id = feats)), contrast_label)
    )
    inputs <- list(contrasts = data.frame(Contrast_name = "A_vs_B",
                                          stringsAsFactors = FALSE))
    config <- list(modes = list(metabolomics = list(
        de            = list(padj_cutoff = 0.05, linear_fc_cutoff = 1.5),
        preprocessing = list(chosen_norm = "pqn", transform = "log2", scaling = "none"),
        effects       = list(samples = "sample_id", color = "grp"),
        # Per-group CV has its own log2-whitelist guards; it is not under test here.
        excel         = list(group_cv = FALSE)
    )))

    zscore_mat <- matrix(0, nrow = 3, ncol = 6,
                         dimnames = list(c("F1", "F2", "F5"),
                                         paste0(samples, ".zscore")))
    clustering_res <- list(
        excel_order = list(ordered_ids = c("F5", "F1", "F2"),
                           zscore_mat  = zscore_mat),
        objects = list(patterns      = data.frame(pattern = "p1"),
                       patterns_list = c("010", "101"))
    )

    list(pre = pre, de_res = de_res, inputs = inputs, config = config,
         clustering_res = clustering_res)
}

build_metab_payload <- function(fx, clustering_res = fx$clustering_res,
                                include_legacy = FALSE, ...) {
    suppressMessages(suppressWarnings(
        build_shiny_payload_metabolomics(
            pre = fx$pre, de_res = fx$de_res, inputs = fx$inputs,
            config = fx$config, clustering_res = clustering_res,
            include_legacy = include_legacy, ...
        )
    ))
}

test_that("de_final_table holds only DE rows, without the cutoff/pass helper columns", {
    de <- build_metab_payload(make_metab_export_fixture())$de_final_table

    expect_s3_class(de, "data.frame")
    expect_setequal(de$feature_id, c("F1", "F2", "F5"))
    expect_false("pass_any_contrast" %in% names(de))
    expect_false(any(startsWith(names(de), "manual_cutoffs")))
    expect_true(all(c("linearFC.A_vs_B", "pvalue.A_vs_B", "padj.A_vs_B",
                      "upDown.A_vs_B") %in% names(de)))
})

test_that("de_final_table is richer than de_sig_stats and keeps original_id after feature_id", {
    p  <- build_metab_payload(make_metab_export_fixture())
    de <- p$de_final_table

    # Per-sample expression comes from the final-results builder; de_sig_stats has none.
    expect_true(all(paste0("s", 1:6) %in% names(de)))
    expect_false(any(paste0("s", 1:6) %in% names(p$de_sig_stats)))
    expect_identical(names(de)[1:2], c("feature_id", "original_id"))
})

test_that("de_final_table adds clustering order and z-scores from excel_order", {
    de <- build_metab_payload(make_metab_export_fixture())$de_final_table
    ord <- setNames(de$order, de$feature_id)

    expect_identical(ord[["F5"]], 1L)
    expect_identical(ord[["F1"]], 2L)
    expect_identical(ord[["F2"]], 3L)
    expect_true(all(paste0("s", 1:6, ".zscore") %in% names(de)))
})

test_that("de_final_table omits order and z-scores when there is no clustering order", {
    fx <- make_metab_export_fixture()
    de <- build_metab_payload(fx, clustering_res = NULL)$de_final_table

    expect_setequal(de$feature_id, c("F1", "F2", "F5"))
    expect_false("order" %in% names(de))
    expect_false(any(grepl("\\.zscore$", names(de))))
})

test_that("de_final_table is a zero-row data.frame when no feature passes", {
    fx <- make_metab_export_fixture(pass = rep(0L, 5))
    de <- build_metab_payload(fx)$de_final_table

    expect_s3_class(de, "data.frame")
    expect_identical(nrow(de), 0L)
})

test_that("a failing final-results build warns and leaves de_final_table NULL", {
    # The contrast label is absent from summary_df, so the builder cannot find
    # its stat columns.
    fx <- make_metab_export_fixture(contrast_label = "A_vs_C")

    expect_warning(
        p <- suppressMessages(build_shiny_payload_metabolomics(
            pre = fx$pre, de_res = fx$de_res, inputs = fx$inputs,
            config = fx$config, clustering_res = fx$clustering_res,
            include_legacy = FALSE
        )),
        "de_final_table"
    )
    expect_null(p$de_final_table)
})

test_that("clust_patterns_list is carried over from the clustering objects", {
    p <- build_metab_payload(make_metab_export_fixture())
    expect_identical(p$clust_patterns_list, c("010", "101"))
})

test_that("clust_patterns_list stays present as NULL when patterns have no list", {
    fx <- make_metab_export_fixture()
    fx$clustering_res$objects$patterns_list <- NULL
    p <- build_metab_payload(fx)

    expect_true("clust_patterns_list" %in% names(p))
    expect_null(p$clust_patterns_list)
})

# --- mummichog extension --------------------------------------------------

# Columns as mummichog delivers them (literal `p-value`), see test-mummichog-plots.R.
make_mummichog_pathways <- function() {
    data.frame(
        check.names      = FALSE,
        stringsAsFactors = FALSE,
        pathway      = c("Pathway one", "Pathway two", "Pathway three"),
        overlap_size = c(5, 3, 8),
        pathway_size = c(10, 6, 12),
        "p-value"    = c(0.02, 0.20, 0.60)
    )
}

test_that("mummichog is NULL, with the key present, when it did not run or has no result", {
    fx <- make_metab_export_fixture()
    for (input in list(NULL, list(), list(A_vs_B = NULL))) {
        p <- build_metab_payload(fx, mummichog_pathways = input)
        expect_true("mummichog" %in% names(p))
        expect_null(p$mummichog)
    }
    expect_null(build_metab_payload(fx)$mummichog)   # argument omitted
})

test_that("mummichog sections match build_mummichog_report_sections() output", {
    fx <- make_metab_export_fixture()
    by_contrast <- list(A_vs_B = make_mummichog_pathways())

    p      <- build_metab_payload(fx, mummichog_pathways = by_contrast)
    direct <- build_mummichog_report_sections(by_contrast, fx$config)

    expect_named(p$mummichog, names(direct))
    sec <- p$mummichog[["A_vs_B"]]
    dir <- direct[["A_vs_B"]]

    expect_named(sec, c("title", "subtitle", "plot", "table", "slug"),
                 ignore.order = TRUE)
    expect_identical(sec$title, dir$title)
    expect_identical(sec$subtitle, dir$subtitle)
    expect_identical(sec$slug, dir$slug)
    expect_identical(sec$table, dir$table)
    # Two separately built ggplots never compare identical() (each holds its own
    # plot_env), so compare what defines the plot: class, labels and data.
    expect_s3_class(sec$plot, "ggplot")
    expect_identical(sec$plot$labels$title, dir$plot$labels$title)
    expect_identical(sec$plot$data, dir$plot$data)
    expect_length(sec$plot$layers, length(dir$plot$layers))
})

test_that("mummichog keeps only contrasts that have a result", {
    fx <- make_metab_export_fixture()
    p <- build_metab_payload(
        fx, mummichog_pathways = list(A_vs_B = make_mummichog_pathways(),
                                      no_result = NULL))
    expect_named(p$mummichog, "A_vs_B")
})

test_that("mummichog does not depend on include_legacy", {
    fx <- make_metab_export_fixture()
    by_contrast <- list(A_vs_B = make_mummichog_pathways())

    for (legacy in c(FALSE, TRUE)) {
        p <- build_metab_payload(fx, include_legacy = legacy,
                                 mummichog_pathways = by_contrast)
        expect_named(p$mummichog, "A_vs_B", info = paste("include_legacy =", legacy))
    }
})
