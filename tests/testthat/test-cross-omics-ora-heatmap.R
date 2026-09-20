# What reaches the cross-omics ORA figure, and what a blank cell means.
#
# The figure shows the adjusted p-value each omics layer's own ORA already
# reported. It does not re-adjust anything, it does not accept a GSEA row, and
# it does not substitute a raw p-value when the adjusted one is absent -- the
# caption says "adjusted", so a cell that cannot honour that stays empty.
#
# A blank cell means no adjusted ORA p-value for that pathway and layer reached
# the table behind the figure. Some enrichment tables arrive already filtered,
# so that is not the same as "the pathway was not tested there", and nothing
# here can tell the two apart.
#
# All fixtures are synthetic.

ora_rows <- function(id, padj, pvalue = padj / 100, method = "ora",
                     contrast = "A_vs_B", pathway = paste("Pathway", id)) {
    data.frame(
        pathway  = pathway,
        ID       = id,
        pvalue   = pvalue,
        padj     = padj,
        method   = method,
        contrast = contrast,
        stringsAsFactors = FALSE
    )
}


# ---- ORA membership --------------------------------------------------------

test_that("a GSEA row never supplies a cell, even beside an ORA row for the same pathway", {
    tabs <- list(transcriptomics = rbind(
        ora_rows("map00010", padj = 1e-2),
        ora_rows("map00010", padj = 1e-8, method = "fgsea")
    ))

    m <- build_ora_adjusted_p_matrix(tabs, "00010", "transcriptomics")

    # The minimum over the layer would be the fgsea value if membership were not
    # enforced; it is the ORA one.
    expect_equal(unname(m["00010", "transcriptomics"]), 1e-2)
})

test_that("a layer holding only GSEA results contributes nothing, and that is not an error", {
    tabs <- list(
        transcriptomics = ora_rows("map00010", padj = 1e-3),
        proteomics      = ora_rows("map00010", padj = 1e-9, method = "fgsea")
    )

    # A valid table that simply has no ORA row is a normal outcome, not a schema
    # problem, so it must not warn.
    expect_no_warning(
        m <- build_ora_adjusted_p_matrix(tabs, "00010",
                                         c("transcriptomics", "proteomics"))
    )
    expect_equal(unname(m["00010", "transcriptomics"]), 1e-3)
    expect_true(is.na(m["00010", "proteomics"]))
})

test_that("a layer with no method column is refused rather than assumed to be ORA", {
    no_method <- ora_rows("map00010", padj = 1e-6)
    no_method$method <- NULL
    tabs <- list(transcriptomics = ora_rows("map00010", padj = 1e-3),
                 proteomics      = no_method)

    expect_warning(
        m <- build_ora_adjusted_p_matrix(tabs, "00010",
                                         c("transcriptomics", "proteomics")),
        "proteomics"
    )
    expect_true(is.na(m["00010", "proteomics"]))
    # The other layer is unaffected: one malformed table does not blank the rest.
    expect_equal(unname(m["00010", "transcriptomics"]), 1e-3)
})

test_that("rows whose method is missing are dropped without dropping the layer", {
    # bind_rows() NA-fills `method` when a layer stacks tables that do not all
    # carry it, so a present column can still hold unknown rows.
    mixed <- rbind(ora_rows("map00010", padj = 1e-3),
                   ora_rows("map00020", padj = 1e-9))
    mixed$method[2] <- NA_character_

    m <- build_ora_adjusted_p_matrix(list(transcriptomics = mixed),
                                     c("00010", "00020"), "transcriptomics")

    expect_equal(unname(m["00010", "transcriptomics"]), 1e-3)
    expect_true(is.na(m["00020", "transcriptomics"]))
})

test_that("method matching ignores case", {
    tabs <- list(transcriptomics = ora_rows("map00010", padj = 1e-3,
                                            method = "ORA"))

    m <- build_ora_adjusted_p_matrix(tabs, "00010", "transcriptomics")

    expect_equal(unname(m["00010", "transcriptomics"]), 1e-3)
})


# ---- which statistic the cell holds ----------------------------------------

test_that("a layer with no adjusted p-value column never falls back to the raw one", {
    raw_only <- ora_rows("map00010", padj = 1e-6)
    raw_only$padj <- NULL
    raw_only$pvalue <- 1e-12   # very small, and still not what this figure shows
    tabs <- list(transcriptomics = ora_rows("map00010", padj = 1e-3),
                 proteomics      = raw_only)

    expect_warning(
        m <- build_ora_adjusted_p_matrix(tabs, "00010",
                                         c("transcriptomics", "proteomics")),
        "proteomics"
    )
    expect_true(is.na(m["00010", "proteomics"]))
})

test_that("the adjusted p-value is resolved per row, not per table", {
    # extract_enrichment_df() binds a layer's sub-results with bind_rows(), which
    # NA-fills, so one table can carry fgsea rows keyed on padj beside ORA rows
    # carrying only clusterProfiler's p.adjust. Picking one column for the whole
    # table reads NA for every row that used the other -- silently, because the
    # column check passed.
    bound <- data.frame(
        pathway  = c("Glycolysis", "Citrate cycle"),
        ID       = c("map00010", "map00020"),
        pvalue   = c(1e-5, 1e-6),
        padj     = c(1e-4, NA_real_),
        p.adjust = c(NA_real_, 1e-3),
        method   = c("fgsea", "ora"),
        stringsAsFactors = FALSE
    )

    vals <- .ora_adjusted_p_values(bound)
    expect_equal(vals, c(1e-4, 1e-3))

    # And end to end: the ORA row's value reaches the matrix, the fgsea row's
    # does not, and the layer is not written off as having nothing.
    m <- build_ora_adjusted_p_matrix(list(transcriptomics = bound),
                                     c("00010", "00020"), "transcriptomics")
    expect_true(is.na(m["00010", "transcriptomics"]))
    expect_equal(unname(m["00020", "transcriptomics"]), 1e-3)
})

test_that("resolving per row prefers padj and never invents a value", {
    both <- data.frame(padj = 0.01, p.adjust = 0.99, stringsAsFactors = FALSE)
    expect_equal(.ora_adjusted_p_values(both), 0.01)

    neither <- data.frame(pvalue = 1e-9, qvalue = 1e-9, stringsAsFactors = FALSE)
    expect_true(is.na(.ora_adjusted_p_values(neither)))

    empty_row <- data.frame(padj = NA_real_, p.adjust = NA_real_)
    expect_true(is.na(.ora_adjusted_p_values(empty_row)))
})

test_that("clusterProfiler's p.adjust is read, and qvalue is not", {
    cp <- ora_rows("map00010", padj = 1e-4)
    names(cp)[names(cp) == "padj"] <- "p.adjust"

    m <- build_ora_adjusted_p_matrix(list(proteomics = cp), "00010", "proteomics")
    expect_equal(unname(m["00010", "proteomics"]), 1e-4)

    # qvalue is a different adjustment, so it is not an alias for this one.
    qv <- ora_rows("map00010", padj = 1e-4)
    names(qv)[names(qv) == "padj"] <- "qvalue"
    expect_warning(
        m2 <- build_ora_adjusted_p_matrix(list(proteomics = qv), "00010",
                                          "proteomics"),
        "adjusted"
    )
    expect_true(is.na(m2["00010", "proteomics"]))
})

test_that("the cell carries the layer's value unchanged, with no adjustment of its own", {
    # Two rows would give a BH-adjusted pair different from either input; the
    # cell is the input.
    tabs <- list(transcriptomics = rbind(ora_rows("map00010", padj = 0.04),
                                         ora_rows("map00020", padj = 0.04)))

    m <- build_ora_adjusted_p_matrix(tabs, c("00010", "00020"), "transcriptomics")

    expect_equal(unname(m[, "transcriptomics"]), c(0.04, 0.04))
})


# ---- several contrasts within one layer ------------------------------------

test_that("contrast-level results for one pathway reduce to the smallest adjusted p-value", {
    tabs <- list(transcriptomics = rbind(
        ora_rows("map00010", padj = 1e-2, contrast = "A_vs_B"),
        ora_rows("map00010", padj = 1e-5, contrast = "C_vs_D"),
        ora_rows("map00010", padj = 1e-3, contrast = "E_vs_F")
    ))

    m <- build_ora_adjusted_p_matrix(tabs, "00010", "transcriptomics")

    expect_equal(unname(m["00010", "transcriptomics"]), 1e-5)
})

test_that("a single contrast is left exactly as it arrived", {
    tabs <- list(transcriptomics = ora_rows("map00010", padj = 0.031))

    m <- build_ora_adjusted_p_matrix(tabs, "00010", "transcriptomics")

    expect_equal(unname(m["00010", "transcriptomics"]), 0.031)
})


# ---- pathway identity ------------------------------------------------------

test_that("the layers' spellings of one KEGG pathway land on one row", {
    tabs <- list(
        transcriptomics = ora_rows("hsa00010", padj = 1e-3),
        proteomics      = ora_rows("ko00010",  padj = 1e-4),
        metabolomics    = ora_rows("map00010", padj = 1e-5)
    )

    m <- build_ora_adjusted_p_matrix(
        tabs, "00010", c("transcriptomics", "proteomics", "metabolomics"),
        kegg_org = "hsa")

    expect_equal(nrow(m), 1L)
    expect_equal(unname(m["00010", ]), c(1e-3, 1e-4, 1e-5))
})

test_that("identity is resolved per row, so a mixed table keys each row on what it has", {
    mixed <- rbind(ora_rows("map00010", padj = 1e-3),
                   ora_rows("map00020", padj = 1e-4))
    # A custom collection row carries a gene-set name and no accession.
    mixed$ID[2] <- NA_character_
    mixed$pathway[2] <- "Custom response module"

    m <- build_ora_adjusted_p_matrix(list(transcriptomics = mixed),
                                     c("00010", "Custom response module"),
                                     "transcriptomics")

    expect_equal(unname(m[, "transcriptomics"]), c(1e-3, 1e-4))
})


# ---- shape, and what a blank cell means ------------------------------------

test_that("a pathway no layer enriched is all NA rather than absent or zero", {
    tabs <- list(transcriptomics = ora_rows("map00010", padj = 1e-3),
                 proteomics      = ora_rows("map00010", padj = 1e-4))

    m <- build_ora_adjusted_p_matrix(tabs, c("00010", "00020"),
                                     c("transcriptomics", "proteomics"))

    expect_equal(rownames(m), c("00010", "00020"))
    expect_true(all(is.na(m["00020", ])))
})

test_that("complementary missingness survives the builder", {
    # Two pathways sharing no observed layer -- the shape a union candidate
    # universe makes routine, and the one that breaks row clustering.
    tabs <- list(
        transcriptomics = ora_rows("map00010", padj = 1e-3),
        proteomics      = ora_rows("map00020", padj = 1e-4)
    )

    m <- build_ora_adjusted_p_matrix(tabs, c("00010", "00020"),
                                     c("transcriptomics", "proteomics"))

    expect_true(is.na(m["00010", "proteomics"]))
    expect_true(is.na(m["00020", "transcriptomics"]))
    expect_equal(unname(m["00010", "transcriptomics"]), 1e-3)
    expect_equal(unname(m["00020", "proteomics"]), 1e-4)
})

test_that("rows outside the target set do not appear", {
    tabs <- list(transcriptomics = rbind(ora_rows("map00010", padj = 1e-3),
                                         ora_rows("map00999", padj = 1e-9)))

    m <- build_ora_adjusted_p_matrix(tabs, "00010", "transcriptomics")

    expect_equal(rownames(m), "00010")
})

test_that("an absent or empty layer is a column of NA, not an error or a warning", {
    tabs <- list(transcriptomics = ora_rows("map00010", padj = 1e-3),
                 proteomics      = NULL,
                 metabolomics    = ora_rows("map00010", padj = 1e-3)[0, ])

    expect_no_warning(
        m <- build_ora_adjusted_p_matrix(
            tabs, "00010", c("transcriptomics", "proteomics", "metabolomics"))
    )
    expect_true(all(is.na(m["00010", c("proteomics", "metabolomics")])))
})

test_that("no pathways and no layers both give an empty matrix rather than an error", {
    tabs <- list(transcriptomics = ora_rows("map00010", padj = 1e-3))

    expect_equal(dim(build_ora_adjusted_p_matrix(tabs, character(0),
                                                 "transcriptomics")),
                 c(0L, 1L))
    expect_equal(dim(build_ora_adjusted_p_matrix(tabs, "00010", character(0))),
                 c(1L, 0L))
})


# ---- which rows the figure leads with --------------------------------------

test_that("the selector ranks on the columns it is given, not on fixed ones", {
    # A table carrying both pairs: the meta-analysis pair would lead with 00010,
    # the ORA pair with 00020. The ORA figure must not be ranked on raw-p
    # meta-analysis columns, so the caller names its own.
    both <- data.frame(
        norm_id       = c("00010", "00020"),
        n_omics       = c(2L, 1L),
        combined_pval = c(1e-9, 1e-2),
        n_ora_layers  = c(1L, 2L),
        best_ora_padj = c(1e-9, 1e-2),
        stringsAsFactors = FALSE
    )

    expect_identical(select_multi_omics_pathways(both, top_n = 2)$norm_id,
                     c("00010", "00020"))
    expect_identical(
        select_multi_omics_pathways(both, top_n = 2,
                                    count_col = "n_ora_layers",
                                    score_col = "best_ora_padj")$norm_id,
        c("00020", "00010"))
})

test_that("the ORA ranking keeps #207's shape: more layers first, then the best value", {
    ranking <- data.frame(
        norm_id       = c("00010", "00020", "00030"),
        n_ora_layers  = c(1L, 2L, 2L),
        best_ora_padj = c(1e-12, 1e-5, 1e-3),
        stringsAsFactors = FALSE
    )

    expect_identical(
        select_multi_omics_pathways(ranking, top_n = 3,
                                    count_col = "n_ora_layers",
                                    score_col = "best_ora_padj")$norm_id,
        c("00020", "00030", "00010"))
})


# ---- drawing ---------------------------------------------------------------

test_that("a matrix with nothing to show says so instead of drawing an empty grid", {
    m <- matrix(NA_real_, nrow = 2, ncol = 2,
                dimnames = list(c("00010", "00020"),
                                c("transcriptomics", "proteomics")))

    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    expect_null(plot_cross_omics_ora_heatmap(m))
})

test_that("pathways with no adjusted ORA p-value anywhere do not take a slot", {
    # Two empty rows sit ahead of the two with evidence. With room for two, it
    # must be the two with something to show.
    m <- matrix(NA_real_, nrow = 4, ncol = 2,
                dimnames = list(c("00010", "00020", "00030", "00040"),
                                c("transcriptomics", "proteomics")))
    m["00030", ] <- c(1e-3, 1e-4)
    m["00040", "proteomics"] <- 1e-9

    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    drawn <- plot_cross_omics_ora_heatmap(m, top_n = 2)

    # 00030 leads despite the weaker value: two layers contribute to it and one
    # to 00040. Labels are the keys here because no pathway_tables were given.
    expect_identical(rownames(drawn), c("00030", "00040"))
})

test_that("a missing cell reaches the drawn matrix as NA, not as zero", {
    # Complementary missingness on purpose: the cap assigns through a logical
    # subscript, and a logical subscript carrying NA is an error in `[<-`, so
    # this is also the drawing path a union candidate universe makes routine.
    m <- matrix(c(1e-3, NA, NA, 1e-4), nrow = 2, byrow = TRUE,
                dimnames = list(c("00010", "00020"),
                                c("transcriptomics", "proteomics")))

    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    drawn <- plot_cross_omics_ora_heatmap(m)

    expect_true(is.na(drawn["00010", "proteomics"]))
    expect_true(is.na(drawn["00020", "transcriptomics"]))
    # And the values that are there are the -log10 of what arrived, capped at 10.
    expect_equal(unname(drawn["00010", "transcriptomics"]), 3)
    expect_equal(unname(drawn["00020", "proteomics"]), 4)
})

test_that("very small adjusted p-values are capped rather than erroring on the NA beside them", {
    m <- matrix(c(1e-30, NA, NA, 1e-4), nrow = 2, byrow = TRUE,
                dimnames = list(c("00010", "00020"),
                                c("transcriptomics", "proteomics")))

    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    drawn <- plot_cross_omics_ora_heatmap(m)

    expect_equal(unname(drawn["00010", "transcriptomics"]), 10)
    expect_true(is.na(drawn["00010", "proteomics"]))
})

test_that("a matrix with a single row draws instead of failing", {
    # heatmap() refuses fewer than two rows or columns even with both dendrograms
    # off, and one surviving pathway is ordinary in a per-contrast view. Drawn
    # through either branch, this must produce a figure rather than the outer
    # tryCatch's "failed" placeholder.
    m <- matrix(c(1e-3, 1e-4), nrow = 1,
                dimnames = list("00010", c("transcriptomics", "proteomics")))

    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    drawn <- NULL
    expect_no_error(drawn <- plot_cross_omics_ora_heatmap(m))
    expect_identical(rownames(drawn), "00010")
})

test_that("a single cell of evidence draws instead of failing", {
    # One row and one column at once, and a zero-width value range with it.
    m <- matrix(1e-3, nrow = 1, dimnames = list("00010", "metabolomics"))

    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    expect_no_error(plot_cross_omics_ora_heatmap(m))
})

test_that("a matrix whose values all agree draws instead of failing", {
    # pheatmap derives its breaks from the minimum and maximum; identical values
    # make every break the same and cut() refuses them.
    m <- matrix(1e-3, nrow = 2, ncol = 2,
                dimnames = list(c("00010", "00020"),
                                c("transcriptomics", "proteomics")))

    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    expect_no_error(plot_cross_omics_ora_heatmap(m))
})

test_that("the meta-analysis heatmap survives the same two shapes", {
    # The same defects were in plot_cross_omics_pathway_heatmap(), which sits one
    # section away in the report. Swept with this change rather than left.
    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    one_row <- data.frame(pathway = "Glycolysis", norm_id = "00010",
                          pval_transcriptomics = 1e-3, pval_proteomics = 1e-4,
                          n_omics = 2L, combined_pval = 1e-4,
                          stringsAsFactors = FALSE)
    expect_no_error(plot_cross_omics_pathway_heatmap(
        one_row, c("transcriptomics", "proteomics")))

    all_equal <- data.frame(
        pathway = c("Glycolysis", "Citrate cycle"),
        norm_id = c("00010", "00020"),
        pval_transcriptomics = c(1e-3, 1e-3),
        pval_proteomics = c(1e-3, 1e-3),
        n_omics = c(2L, 2L), combined_pval = c(1e-3, 1e-3),
        stringsAsFactors = FALSE
    )
    expect_no_error(plot_cross_omics_pathway_heatmap(
        all_equal, c("transcriptomics", "proteomics")))
})

test_that("a widened range is only ever used where the real one has no width", {
    expect_equal(.nondegenerate_range(matrix(c(2, 5))), c(2, 5))
    expect_equal(.nondegenerate_range(matrix(c(3, NA, 3))), c(2.5, 3.5))
    expect_equal(.nondegenerate_range(matrix(c(4, NA))), c(3.5, 4.5))
    # Nothing finite to draw from: a usable range rather than c(Inf, -Inf).
    expect_equal(.nondegenerate_range(matrix(NA_real_)), c(0, 1))
    # Whatever it returns, the two ends differ -- that is the whole contract.
    expect_true(.nondegenerate_range(matrix(NA_real_))[1] <
                .nondegenerate_range(matrix(NA_real_))[2])
})

test_that("the base fallback keeps the NA handling the primary path has", {
    # Not reachable behaviourally: which of the two branches draws depends on
    # whether pheatmap is installed on the machine running the tests, so the
    # fallback's own NA handling is pinned at the source instead.
    body_src <- paste(deparse(body(plot_cross_omics_ora_heatmap)), collapse = " ")

    # Both dendrograms off -- rows whose missing layers do not overlap share no
    # observed cell, so dist() is NA between them and hclust() stops on it.
    expect_true(grepl("Rowv = NA", body_src, fixed = TRUE))
    expect_true(grepl("Colv = NA", body_src, fixed = TRUE))
    # image() does not paint an NA cell, so the device background shows through.
    expect_true(grepl("with_par", body_src, fixed = TRUE))
    expect_true(grepl('bg = "grey90"', body_src, fixed = TRUE))
    # And nothing fills a missing value to make either of those unnecessary.
    expect_false(grepl("is.na(log_padj_matrix)] <- ", body_src, fixed = TRUE))
    expect_true(grepl("!is.na(log_padj_matrix)", body_src, fixed = TRUE))
})

test_that("drawing leaves the graphics state as it found it", {
    # Only the fallback branch sets a background, so on a machine with pheatmap
    # this asserts that the primary branch leaves par alone; the mechanism the
    # fallback uses to restore it is pinned separately below.
    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)
    before <- graphics::par("bg")

    m <- matrix(c(1e-3, NA, NA, 1e-4), nrow = 2, byrow = TRUE,
                dimnames = list(c("00010", "00020"),
                                c("transcriptomics", "proteomics")))
    plot_cross_omics_ora_heatmap(m)

    expect_identical(graphics::par("bg"), before)
})

test_that("a background set for the fallback does not leak into the next plot", {
    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)
    before <- graphics::par("bg")

    m <- matrix(c(2, NA, NA, 3), nrow = 2, byrow = TRUE,
                dimnames = list(c("00010", "00020"),
                                c("transcriptomics", "proteomics")))
    withr::with_par(list(bg = "grey90"), {
        stats::heatmap(m, scale = "none", Rowv = NA, Colv = NA,
                       col = grDevices::colorRampPalette(c("white", "red"))(5))
    })

    expect_identical(graphics::par("bg"), before)
})


# ---- the figure a run declines to draw --------------------------------------

.ora_run_tables <- function(method = "ora") {
    layer <- function(p1, p2) {
        df <- data.frame(
            pathway = c("Glycolysis", "Citrate cycle"),
            ID      = c("map00010", "map00020"),
            pvalue  = c(p1, p2) / 10,
            padj    = c(p1, p2),
            stringsAsFactors = FALSE
        )
        if (!is.null(method)) df$method <- method
        df
    }
    list(transcriptomics = layer(1e-2, 1e-3),
         proteomics      = layer(1e-4, 1e-5))
}

.ora_run_config <- list(
    global = list(organism = "Homo sapiens"),
    modes  = list(multiomics = list(enrichment = list()))
)

test_that("a run that can draw the ORA figure writes it", {
    # One file per gene-set collection. These fixtures are KEGG map accessions
    # throughout, so KEGG is the only collection with anything to draw.
    out_dir <- withr::local_tempdir()

    suppressWarnings(suppressMessages(analyze_cross_omics_enrichment(
        .ora_run_tables("ora"), .ora_run_config, out_dir = out_dir)))

    expect_true(file.exists(file.path(out_dir,
                                      "cross_omics_ora_heatmap_KEGG.png")))
})

test_that("a run that skips the ORA figure deletes the previous run's copy", {
    # targets reuses an output directory, and the report finds these figures by
    # name. Left in place, last run's ORA evidence would be read as this run's.
    # Both spellings must go: the fixed name these figures used before they were
    # split per collection, and a collection this run does not produce.
    out_dir <- withr::local_tempdir()
    stale <- c(file.path(out_dir, "cross_omics_ora_heatmap.png"),
               file.path(out_dir, "cross_omics_ora_heatmap_KEGG.png"))
    for (f in stale) writeLines("left over from an earlier run", f)

    # No method column anywhere, so no layer can be confirmed as ORA.
    suppressWarnings(suppressMessages(analyze_cross_omics_enrichment(
        .ora_run_tables(NULL), .ora_run_config, out_dir = out_dir)))

    expect_false(any(file.exists(stale)))
})

test_that("skipping the figure leaves the rest of the run's output alone", {
    out_dir <- withr::local_tempdir()
    bystander <- file.path(out_dir, "cross_omics_pathway_heatmap_KEGG.png")

    res <- suppressWarnings(suppressMessages(analyze_cross_omics_enrichment(
        .ora_run_tables(NULL), .ora_run_config, out_dir = out_dir)))

    # Only the ORA figures are absent; the meta-analysis heatmap is written on
    # every run and is still here, and the result carries no ORA figure path
    # under any collection.
    expect_true(file.exists(bystander))
    expect_length(grep("^ora_heatmap", names(res$plots), value = TRUE), 0L)
})


# ---- the drawing path pheatmap hides ----------------------------------------

# Which heatmap branch runs depends on whether pheatmap is installed, so the
# tests above exercise the fallback only on a machine without it -- which is not
# the machine CI runs on. These call it directly, so the shapes it exists for
# are covered either way.

test_that("the small-matrix path draws every shape heatmap() refuses", {
    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    cols <- grDevices::colorRampPalette(c("white", "red"))(50)
    one_cell <- matrix(3, nrow = 1, dimnames = list("00010", "metabolomics"))
    one_row  <- matrix(c(3, 4), nrow = 1,
                       dimnames = list("00010", c("transcriptomics", "proteomics")))
    one_col  <- matrix(c(3, 4), ncol = 1,
                       dimnames = list(c("00010", "00020"), "transcriptomics"))

    for (m in list(one_cell, one_row, one_col)) {
        expect_no_error(.draw_small_heatmap(m, main = "t", col = cols,
                                            zlim = .nondegenerate_range(m)))
    }
})

test_that("the small-matrix path keeps a missing cell missing", {
    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    # Nothing is filled in to make the shape drawable; the cell stays NA and the
    # grey device background shows through it, as in the two-dimensional path.
    m <- matrix(c(3, NA), nrow = 1,
                dimnames = list("00010", c("transcriptomics", "proteomics")))
    cols <- grDevices::colorRampPalette(c("white", "red"))(50)

    expect_no_error(.draw_small_heatmap(m, main = "t", col = cols,
                                        zlim = .nondegenerate_range(m)))
    expect_true(is.na(m["00010", "proteomics"]))
})
