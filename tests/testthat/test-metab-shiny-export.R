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

    expect_named(sec, c("title", "subtitle", "plot", "table", "slug",
                        "summary_table", "ora_plot", "ora_table",
                        "gsea_plot", "gsea_table",
                        "overview_table", "ora_detail_table", "input_counts",
                        "evidence_pathway_table", "evidence_ec_table",
                        "evidence_feature_table", "evidence_note",
                        "evidence_log2fc_note", "p_cutoff"))
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


# --- metabolite names -----------------------------------------------------------
# Name is an annotation: it comes from row_data$Name, joined by feature_id, in
# every DE table. original_id is the researcher's own id and may be a generic
# "feature_N"; it must never stand in for the name.

METAB_NAMES <- c(F1 = "alanine", F2 = NA, F3 = "citrate", F4 = "glucose",
                 F5 = "lactate")

make_metab_name_fixture <- function() {
    fx <- make_metab_export_fixture()
    feats <- names(METAB_NAMES)
    row_data <- data.frame(feature_id  = feats,
                           original_id = paste0("feature_", seq_along(feats)),
                           Name        = unname(METAB_NAMES),
                           row.names   = feats, stringsAsFactors = FALSE)
    # Annotation rows in a different order from the DE table.
    fx$pre$row_data <- row_data[c(3, 1, 5, 2, 4), , drop = FALSE]
    fx
}

test_that("de_stats and de_sig_stats take Name from row_data by feature_id", {
    fx <- make_metab_name_fixture()
    p  <- build_metab_payload(fx)

    for (key in c("de_stats", "de_sig_stats")) {
        tab <- p[[key]]
        expect_identical(tab$Name, unname(METAB_NAMES[tab$feature_id]), info = key)
        expect_identical(tab$original_id,
                         sub("^F", "feature_", tab$feature_id), info = key)
    }
    # The ids themselves are untouched, in the DE table's own order.
    expect_identical(p$de_stats$feature_id, fx$de_res$summary_df$feature_id)
})

test_that("de_final_table carries the same names, beside feature_id and original_id", {
    p  <- build_metab_payload(make_metab_name_fixture())
    de <- p$de_final_table

    expect_identical(names(de)[1:3], c("feature_id", "original_id", "Name"))
    expect_identical(de$Name, unname(METAB_NAMES[de$feature_id]))
    expect_identical(de$original_id, sub("^F", "feature_", de$feature_id))
})

test_that("Name agrees across de_stats, de_final_table and feature_annot", {
    p <- build_metab_payload(make_metab_name_fixture())

    from_annot <- function(ids) p$feature_annot[ids, "Name"]
    expect_identical(p$de_stats$Name, from_annot(p$de_stats$feature_id))
    expect_identical(p$de_final_table$Name, from_annot(p$de_final_table$feature_id))
})

test_that("a feature without an annotation row gets NA, not another feature's name", {
    fx <- make_metab_name_fixture()
    fx$pre$row_data <- fx$pre$row_data[rownames(fx$pre$row_data) != "F1", , drop = FALSE]
    p  <- build_metab_payload(fx)

    expected <- unname(METAB_NAMES[p$de_stats$feature_id])
    expected[p$de_stats$feature_id == "F1"] <- NA
    expect_identical(p$de_stats$Name, expected)
    expect_true(is.na(p$de_final_table$Name[p$de_final_table$feature_id == "F1"]))
    expect_identical(p$de_final_table$Name[p$de_final_table$feature_id == "F5"], "lactate")
})

test_that("without a Name column in row_data, the DE tables have none", {
    fx <- make_metab_name_fixture()
    fx$pre$row_data$Name <- NULL
    p  <- build_metab_payload(fx)

    expect_false("Name" %in% names(p$de_stats))
    expect_false("Name" %in% names(p$de_final_table))
    expect_true("original_id" %in% names(p$de_stats))
})

test_that("the final-results table written to TSV and Excel names features from row_data", {
    fx <- make_metab_name_fixture()
    final <- suppressMessages(suppressWarnings(build_final_results_metabolomics(
        pre             = fx$pre,
        summary_df      = fx$de_res$summary_df,
        contrast_labels = names(fx$de_res$de_tables),
        row_data        = fx$pre$row_data,
        feature_id_col  = "feature_id",
        cv_contrasts_df = fx$inputs$contrasts,
        config          = fx$config
    )))

    expect_identical(names(final)[1:3], c("feature_id", "original_id", "Name"))
    expect_setequal(final$feature_id, fx$de_res$summary_df$feature_id)
    expect_identical(final$Name, unname(METAB_NAMES[final$feature_id]))
    expect_false(any(final$Name %in% final$original_id))
})


# --- mummichog: explicitly named fields from the report sections --------------
# The pipeline hands the payload builder the sections the HTML report renders
# (metab_mummichog_sections). Each payload entry keeps the five original fields
# and adds summary_table, ora_plot, ora_table, gsea_plot and gsea_table, which
# are objects of the section, not recomputed.

MUMMICHOG_LEGACY_KEYS <- c("title", "subtitle", "plot", "table", "slug")
MUMMICHOG_ADDED_KEYS  <- c("summary_table", "ora_plot", "ora_table",
                           "gsea_plot", "gsea_table")
MUMMICHOG_REPORT_KEYS <- c("overview_table", "ora_detail_table", "input_counts",
                           "evidence_pathway_table", "evidence_ec_table",
                           "evidence_feature_table", "evidence_note",
                           "evidence_log2fc_note", "p_cutoff")
MUMMICHOG_ALL_KEYS    <- c(MUMMICHOG_LEGACY_KEYS, MUMMICHOG_ADDED_KEYS,
                           MUMMICHOG_REPORT_KEYS)

# More pathways than the Summary keeps, so a truncated table would show.
make_many_mummichog_pathways <- function(n = 30) {
    data.frame(
        check.names      = FALSE,
        stringsAsFactors = FALSE,
        pathway      = sprintf("Pathway %02d", seq_len(n)),
        overlap_size = rep(c(5, 3, 8), length.out = n),
        pathway_size = rep(c(10, 6, 12), length.out = n),
        "p-value"    = seq(0.001, 0.9, length.out = n)
    )
}

# Shaped like run_mummichog_gsea()$table: one row per tested pathway.
make_mummichog_gsea_table <- function(n = 30) {
    data.frame(
        check.names      = FALSE,
        stringsAsFactors = FALSE,
        "Pathway"      = sprintf("Pathway %02d", seq_len(n)),
        "Tested size"  = rep(c(4L, 7L, 9L), length.out = n),
        "ES"           = seq(-0.8, 0.8, length.out = n),
        "NES"          = seq(-2, 2, length.out = n),
        "P.Value"      = seq(0.002, 0.95, length.out = n),
        "padj"         = seq(0.02, 0.99, length.out = n),
        "Leading-edge EmpiricalCompounds" = rep(c("E1; E2", "E3", ""), length.out = n)
    )
}

# Two contrasts as the sections target holds them: the first with GSEA, the
# second without (GSEA switched off or unavailable for it).
make_mummichog_sections <- function(fx) {
    sections <- suppressMessages(build_mummichog_report_sections(
        list(A_vs_B = make_many_mummichog_pathways(),
             C_vs_B = make_mummichog_pathways()),
        fx$config))
    gsea_tbl <- make_mummichog_gsea_table()
    sections$A_vs_B$gsea <- list(table = gsea_tbl, metric_label = "t-statistic")
    # A stand-in for the scatter: the payload passes the section's plot object
    # through, whatever it is.
    sections$A_vs_B$gsea_plot <- ggplot2::ggplot(
        gsea_tbl, ggplot2::aes(x = NES, y = -log10(P.Value))) + ggplot2::geom_point()
    # A stand-in for the evidence a real section carries.
    sections$A_vs_B$evidence <- list(pathway_summary = data.frame(Pathway = "Pathway 01"))
    sections
}

test_that("every mummichog entry keeps its ten fields first, then the report objects", {
    fx <- make_metab_export_fixture()
    p  <- build_metab_payload(fx, mummichog_sections = make_mummichog_sections(fx))

    expect_named(p$mummichog, c("A_vs_B", "C_vs_B"))
    for (contrast in names(p$mummichog)) {
        expect_named(p$mummichog[[contrast]], MUMMICHOG_ALL_KEYS, info = contrast)
    }
    expect_identical(.METAB_SHINY_MUMMICHOG_ADDED_KEYS, MUMMICHOG_ADDED_KEYS)
    expect_identical(.METAB_SHINY_MUMMICHOG_REPORT_KEYS, MUMMICHOG_REPORT_KEYS)
})

test_that("the original mummichog fields are the section's own, unchanged", {
    fx       <- make_metab_export_fixture()
    sections <- make_mummichog_sections(fx)
    p        <- build_metab_payload(fx, mummichog_sections = sections)

    for (contrast in names(sections)) {
        for (key in MUMMICHOG_LEGACY_KEYS) {
            expect_identical(p$mummichog[[contrast]][[key]], sections[[contrast]][[key]],
                             info = paste(contrast, key))
        }
    }
})

test_that("ora_plot and ora_table are the section's ORA plot and table", {
    fx       <- make_metab_export_fixture()
    sections <- make_mummichog_sections(fx)
    entry    <- build_metab_payload(fx, mummichog_sections = sections)$mummichog$A_vs_B

    expect_identical(entry$ora_plot, sections$A_vs_B$plot)
    expect_identical(entry$ora_plot, entry$plot)
    expect_identical(entry$ora_table, sections$A_vs_B$table)
    expect_identical(entry$ora_table, entry$table)
})

test_that("summary_table is the report's summary object, with all its columns", {
    fx       <- make_metab_export_fixture()
    sections <- make_mummichog_sections(fx)
    p        <- build_metab_payload(fx, mummichog_sections = sections)

    for (contrast in names(sections)) {
        expect_identical(p$mummichog[[contrast]]$summary_table,
                         sections[[contrast]]$summary, info = contrast)
    }
    smry <- p$mummichog$A_vs_B$summary_table
    # The report shows the display columns; the values behind them follow.
    expect_identical(names(smry)[seq_along(.MMC_SUMMARY_DISPLAY_COLS)],
                     .MMC_SUMMARY_DISPLAY_COLS)
    expect_true(all(c("Pathway ID", "Nominal Mummichog signal", "ORA p.value",
                      "GSEA P.Value", "GSEA padj", "Supporting ECs",
                      "Supporting features", "Identified methods tested") %in%
                      names(smry)))
    expect_length(names(smry), 15L)
})

test_that("the ORA and GSEA tables keep every row; only the Summary is a top list", {
    fx       <- make_metab_export_fixture()
    sections <- make_mummichog_sections(fx)
    entry    <- build_metab_payload(fx, mummichog_sections = sections)$mummichog$A_vs_B

    expect_identical(nrow(entry$ora_table), 30L)
    expect_identical(nrow(entry$gsea_table), 30L)
    expect_setequal(entry$ora_table$Pathway, make_many_mummichog_pathways()$pathway)
    expect_identical(entry$gsea_table, sections$A_vs_B$gsea$table)
    expect_lte(nrow(entry$summary_table), .MMC_SUMMARY_TOP_N)
    expect_lt(nrow(entry$summary_table), nrow(entry$ora_table))
})

test_that("gsea_plot is the section's GSEA plot", {
    fx       <- make_metab_export_fixture()
    sections <- make_mummichog_sections(fx)
    entry    <- build_metab_payload(fx, mummichog_sections = sections)$mummichog$A_vs_B

    expect_s3_class(entry$gsea_plot, "ggplot")
    expect_identical(entry$gsea_plot, sections$A_vs_B$gsea_plot)
})

test_that("a contrast without GSEA has NULL gsea fields, with the keys present", {
    fx    <- make_metab_export_fixture()
    entry <- build_metab_payload(
        fx, mummichog_sections = make_mummichog_sections(fx))$mummichog$C_vs_B

    expect_true(all(c("gsea_plot", "gsea_table") %in% names(entry)))
    expect_null(entry$gsea_plot)
    expect_null(entry$gsea_table)
    # The ORA side of that contrast is unaffected.
    expect_s3_class(entry$ora_plot, "ggplot")
    expect_identical(nrow(entry$ora_table), 3L)
    expect_s3_class(entry$summary_table, "data.frame")
})

test_that("p-value and statistic columns keep their numeric types", {
    fx       <- make_metab_export_fixture()
    sections <- make_mummichog_sections(fx)
    entry    <- build_metab_payload(fx, mummichog_sections = sections)$mummichog$A_vs_B

    expect_type(entry$ora_table[["p.value"]], "double")
    for (col in c("ES", "NES", "P.Value", "padj")) {
        expect_type(entry$gsea_table[[col]], "double")
        expect_identical(entry$gsea_table[[col]], sections$A_vs_B$gsea$table[[col]],
                         info = col)
    }
    expect_type(entry$gsea_table[["Tested size"]], "integer")
    for (col in c("ORA p.value", "GSEA P.Value", "GSEA padj")) {
        expect_type(entry$summary_table[[col]], "double")
        expect_identical(entry$summary_table[[col]], sections$A_vs_B$summary[[col]],
                         info = col)
    }
})

test_that("every mummichog table can be written as CSV", {
    fx <- make_metab_export_fixture()
    p  <- build_metab_payload(fx, mummichog_sections = make_mummichog_sections(fx))

    for (contrast in names(p$mummichog)) {
        for (key in c("summary_table", "ora_table", "gsea_table")) {
            tbl <- p$mummichog[[contrast]][[key]]
            if (is.null(tbl)) next
            expect_false(any(vapply(tbl, is.list, logical(1))),
                         info = paste(contrast, key))
            out <- withr::local_tempfile(fileext = ".csv")
            utils::write.csv(tbl, out, row.names = FALSE)
            expect_identical(nrow(utils::read.csv(out, check.names = FALSE)), nrow(tbl),
                             info = paste(contrast, key))
        }
    }
})

test_that("a list column is flattened in the payload only, and nothing else changes", {
    fx       <- make_metab_export_fixture()
    sections <- make_mummichog_sections(fx)
    members  <- rep(list(c("E1", "E2"), "E3", character(0)), length.out = 30)
    sections$A_vs_B$gsea$table[["Leading-edge EmpiricalCompounds"]] <- members
    before   <- sections$A_vs_B$gsea$table

    gsea <- build_metab_payload(fx, mummichog_sections = sections)$mummichog$A_vs_B$gsea_table

    expect_identical(gsea[["Leading-edge EmpiricalCompounds"]][1:3], c("E1; E2", "E3", ""))
    expect_type(gsea[["Leading-edge EmpiricalCompounds"]], "character")
    other <- setdiff(names(before), "Leading-edge EmpiricalCompounds")
    expect_identical(gsea[other], before[other])
    expect_identical(names(gsea), names(before))
    # The report section itself still holds the list column.
    expect_true(is.list(sections$A_vs_B$gsea$table[["Leading-edge EmpiricalCompounds"]]))
})

test_that("flatten_list_columns() returns a table without list columns untouched", {
    tbl <- make_mummichog_gsea_table()
    expect_identical(flatten_list_columns(tbl), tbl)
    expect_null(flatten_list_columns(NULL))
})

test_that("evidence and the GSEA internals stay out of the payload", {
    fx    <- make_metab_export_fixture()
    entry <- build_metab_payload(
        fx, mummichog_sections = make_mummichog_sections(fx))$mummichog$A_vs_B

    expect_false(any(c("evidence", "gsea", "summary", "polarity", "engine") %in%
                       names(entry)))
})

test_that("no sections means a NULL mummichog key, still present", {
    fx <- make_metab_export_fixture()
    for (sections in list(NULL, list())) {
        p <- build_metab_payload(fx, mummichog_sections = sections)
        expect_true("mummichog" %in% names(p))
        expect_null(p$mummichog)
    }
})

test_that("sections take precedence over pathway tables; pathway tables alone still work", {
    fx       <- make_metab_export_fixture()
    sections <- make_mummichog_sections(fx)

    both <- build_metab_payload(fx, mummichog_sections = sections,
                                mummichog_pathways = list(other = make_mummichog_pathways()))
    expect_named(both$mummichog, names(sections))

    only_pathways <- build_metab_payload(
        fx, mummichog_pathways = list(A_vs_B = make_mummichog_pathways()))$mummichog$A_vs_B
    expect_named(only_pathways, MUMMICHOG_ALL_KEYS)
    expect_identical(only_pathways$ora_table, only_pathways$table)
    expect_null(only_pathways$gsea_plot)
    expect_null(only_pathways$gsea_table)
    expect_s3_class(only_pathways$overview_table, "data.frame")
    expect_null(only_pathways$evidence_feature_table)
})

# --- mummichog: the report objects, on the MetaboAnalyst-compatible fixture ------
# ma_two_contrast_sections() (helper-mummichog-ma.R), SYNTHETIC:
#   HL_vs_LL  Glycolysis  ORA adjusted;  Pyruvate metabolism  in the ORA table,
#             nominal in neither ORA nor GSEA;  TCA cycle  GSEA only, no
#             significant EC, so not in the ORA table
#   LL_vs_DD  Glycolysis  in the ORA table, not nominal in ORA; GSEA adjusted

test_that("the report objects are the section's own, after the unchanged ten fields", {
    secs <- ma_two_contrast_sections()$sections
    for (contrast in names(secs)) {
        s <- secs[[contrast]]
        e <- metab_shiny_mummichog_entry(s)
        expect_named(e, MUMMICHOG_ALL_KEYS)
        expect_identical(e$summary_table, flatten_list_columns(s$summary), info = contrast)
        expect_identical(e$ora_table, flatten_list_columns(s$table), info = contrast)
        expect_identical(e$gsea_table, flatten_list_columns(s$gsea$table), info = contrast)
        expect_identical(e$overview_table, s$overview, info = contrast)
        expect_identical(e$ora_detail_table, s$ora_detail, info = contrast)
        expect_identical(e$input_counts, s$input_counts, info = contrast)
        expect_identical(e$evidence_pathway_table, s$evidence$pathway_summary,
                         info = contrast)
        expect_identical(e$evidence_ec_table, s$evidence$ec_table, info = contrast)
        expect_identical(e$p_cutoff, s$p_cutoff, info = contrast)
        expect_null(e$evidence_note)
        expect_null(e$evidence_log2fc_note)
        # Complete ORA results: every pathway of the engine's table.
        expect_identical(nrow(e$ora_detail_table), nrow(s$table), info = contrast)
    }
})

test_that("evidence keeps every pathway that has it, whatever its p-value", {
    secs <- ma_two_contrast_sections()$sections
    hl <- metab_shiny_mummichog_entry(secs$HL_vs_LL)
    ll <- metab_shiny_mummichog_entry(secs$LL_vs_DD)
    ov <- hl$overview_table
    expect_gte(ov[["ORA raw p"]][ov$Pathway == "Pyruvate metabolism"], 0.05)
    expect_true("Pyruvate metabolism" %in% hl$evidence_ec_table$Pathway)
    expect_true("Glycolysis" %in% ll$evidence_ec_table$Pathway)
    # GSEA only, with no ORA overlap: in the overview, with no evidence rows.
    expect_true("TCA cycle" %in% ov$Pathway)
    expect_false("TCA cycle" %in% hl$evidence_ec_table$Pathway)
    # The evidence tables name pathways exactly as the overview does.
    for (e in list(hl, ll)) {
        expect_false(anyDuplicated(e$overview_table$Pathway) > 0)
        expect_setequal(e$evidence_pathway_table$Pathway,
                        unique(e$evidence_ec_table$Pathway))
        expect_true(all(e$evidence_ec_table$Pathway %in% e$overview_table$Pathway))
        expect_true(all(e$evidence_feature_table$Pathway %in%
                          e$evidence_pathway_table$Pathway))
    }
})

test_that("the feature evidence gains only the workbook's candidate names", {
    s  <- ma_two_contrast_sections()$sections$HL_vs_LL
    ft <- metab_shiny_mummichog_entry(s)$evidence_feature_table
    orig <- s$evidence$feature_table
    # Every original column, row and value as the section holds it.
    expect_identical(as.list(ft[names(orig)]), as.list(orig))
    expect_identical(setdiff(names(ft), names(orig)), .MMC_FEATURE_CAND_NAME_COL)
    expect_identical(match(.MMC_FEATURE_CAND_NAME_COL, names(ft)),
                     match(.MMC_FEATURE_CAND_COL, names(ft)) + 1L)
    # The k-th name belongs to the k-th candidate id, from the section's map.
    ids <- lapply(ft[[.MMC_FEATURE_CAND_COL]], .mmc_split_cell, sep = ";")
    nms <- lapply(ft[[.MMC_FEATURE_CAND_NAME_COL]], .mmc_split_cell, sep = ";")
    expect_identical(lengths(nms), lengths(ids))
    expected <- unname(s$candidate_names[unlist(ids)])
    expected[is.na(expected)] <- .MMC_NO_MODEL_NAME
    expect_identical(unlist(nms), expected)
    # A feature in a supporting EC that matched no pathway candidate itself.
    f04 <- ft[ft$Feature == "F04", , drop = FALSE]
    expect_true(nrow(f04) > 0)
    expect_true(all(is.na(f04[[.MMC_FEATURE_CAND_NAME_COL]])))
    expect_false(any(f04[["Matched via pathway candidate"]]))
    expect_true(all(f04$Agreement == "Not assessed"))
    # Identifiers stay text.
    expect_true("0012" %in% ft$Feature)
})

test_that("a contrast without evidence has NULL evidence tables and says why", {
    fx  <- ma_evidence_fixture()
    bad <- fx$model
    attr(bad, "sha256") <- strrep("a", 64)
    s <- suppressWarnings(suppressMessages(build_mummichog_report_sections(
        list(HL_vs_LL = fx$pathways), fx$config, files = fx$files,
        row_data = fx$row_data, de_res = fx$de_res, model = bad)))$HL_vs_LL
    e <- metab_shiny_mummichog_entry(s)
    expect_named(e, MUMMICHOG_ALL_KEYS)
    for (key in c("evidence_pathway_table", "evidence_ec_table",
                  "evidence_feature_table")) {
        expect_null(e[[key]])
    }
    expect_identical(e$evidence_note, s$evidence_note)
    expect_match(e$evidence_note, "not the model this contrast was analysed with")
    expect_identical(e$overview_table, s$overview)
})

test_that("the pipeline hands the payload the report's own mummichog sections", {
    skip_if_not_installed("targets")
    # The sections target exists whether or not mummichog is enabled (NULL when
    # it is not), so the payload target names it in both cases.
    for (enabled in c(FALSE, TRUE)) {
        plan <- pipe_metabolomics("pqn", mummichog_enabled = enabled)
        tgt  <- Filter(function(t) identical(t$settings$name, "metab_shiny_payload"),
                       plan)[[1]]
        deps <- unlist(lapply(as.list(tgt$command$expr), all.vars))

        expect_true("metab_mummichog_sections" %in% deps, info = enabled)
        expect_false("metab_mummichog_report_pathways" %in% deps, info = enabled)
    }
})
