# tests/testthat/test-de-integration-inputs.R
#
# Reading DE-integration layers into one shape. Synthetic tables only, written
# to temporary files in the shapes the pipeline's own exports use.

write_table_tmp <- function(df, ext = "tsv") {
    path <- tempfile(fileext = paste0(".", ext))
    utils::write.table(df, path, sep = if (ext == "csv") "," else "\t",
                       quote = FALSE, row.names = FALSE, na = "NA")
    path
}

dei_config <- function(layers, ...) {
    list(project = list(dir = tempdir()), paths = list(raw = "data"),
         params = list(seed = 1),
         modes = list(de_integration = validate_de_integration_config(
             c(list(layers = layers), list(...)))))
}

hits_default <- list(use_table_flag = TRUE, p_cutoff = 0.05, use_adjusted = TRUE,
                     linear_fc_cutoff = 1.5)

# A proteomics summary for contrast A_vs_B: P1 up, P2 down, P3 flat.
prot_summary <- function(contrast = "A_vs_B") {
    df <- data.frame(
        FeatureID = c("P1", "P2;P9", "P3"),
        Genes = c("G1", "G2;G9", ""),
        First.Protein.Description = c("d1", "d2", "d3"),
        stringsAsFactors = FALSE)
    df[[paste0("log2FC.imputs.", contrast)]]      <- c(1.2, -2.0, 0.1)
    df[[paste0("linearRatio.imputs.", contrast)]] <- 2^c(1.2, -2.0, 0.1)
    df[[paste0("linearFC.imputs.", contrast)]]    <- c(2.3, -4, 1.07)
    df[[paste0("pvalue.imputs.", contrast)]]      <- c(0.001, 0.002, 0.8)
    df[[paste0("padj.imputs.", contrast)]]        <- c(0.01, 0.02, 0.9)
    df[[paste0("pass.imputs.", contrast)]]        <- c(1, 1, NA)
    df$pass_any_contrast <- c(1, 1, NA)
    df
}

layer_cfg <- function(path, format = "proteomics_summary", omics = "proteomics",
                      name = "cells", ...) {
    c(list(name = name, label = name, omics_type = omics, format = format, path = path),
      list(...))
}

# An observed block over four samples, two per group, for the features of
# prot_summary(); `contr` is its contrasts file.
obs_block <- function(contr = data.frame(Contrast_name = "A_vs_B", Factor = "Group",
                                         Numerator = "A", Denominator = "B"),
                      ids = c("P1", "P2;P9", "P3"), groups = c("A", "A", "B", "B")) {
    mat <- data.frame(Protein.Group = ids,
                      s1 = c(20, NA, 15), s2 = c(21, NA, 15),
                      s3 = c(19, 22, NA), s4 = c(20, 23, 16),
                      stringsAsFactors = FALSE)
    sheet <- data.frame(SampleName = c("s1", "s2", "s3", "s4"), Group = groups)
    list(matrix = write_table_tmp(mat), samplesheet = write_table_tmp(sheet, "csv"),
         sample_col = "SampleName", contrasts_file = write_table_tmp(contr, "csv"))
}

# A proteomics final table carrying the export's own N.observed.<group> columns.
prot_final <- function(contrast = "A_vs_B", groups = c("A", "B")) {
    df <- prot_summary(contrast)
    names(df)[names(df) == paste0("pass.imputs.", contrast)] <- paste0("upDown.imputs.", contrast)
    df[[paste0("upDown.imputs.", contrast)]] <- c("Up", "Down", NA)
    df[[paste0("N.observed.", groups[1])]] <- c(3, 1, 2)
    df[[paste0("N.observed.", groups[2])]] <- c(3, 3, 2)
    df
}

test_that("a proteomics summary is read with its stored log2FC and its own hit flag", {
    path <- write_table_tmp(prot_summary())
    ly <- read_de_layer(layer_cfg(path), "A_vs_B", dei_config(list(
        layer_cfg(path), layer_cfg(path, name = "media"))), hits_default)
    tab <- ly$tables$A_vs_B
    expect_identical(tab$feature_id, c("P1", "P2;P9", "P3"))
    expect_equal(tab$log2fc, c(1.2, -2.0, 0.1))
    expect_identical(tab$hit, c(TRUE, TRUE, FALSE))
    expect_identical(ly$provenance$log2fc_source, "log2FC.imputs.A_vs_B")
    expect_identical(ly$provenance$padj_source, "padj.imputs.A_vs_B")
    expect_identical(ly$provenance$hit_source, "pass.imputs.A_vs_B")
    # First gene of a multi-gene group; no symbol for P3.
    expect_identical(tab$symbol, c("G1", "G2", NA))
    expect_identical(tab$symbol_source, c("symbol", "symbol", "none"))
})

test_that("without a stored log2FC the sign comes from the signed linear FC", {
    df <- prot_summary()
    df$log2FC.imputs.A_vs_B <- NULL
    # The pre-computed loader writes an unsigned ratio, 2^|log2FC|.
    df$linearRatio.imputs.A_vs_B <- 2^abs(c(1.2, -2.0, 0.1))
    path <- write_table_tmp(df)
    ly <- read_de_layer(layer_cfg(path), "A_vs_B", NULL, hits_default)
    expect_equal(ly$tables$A_vs_B$log2fc, c(1.2, -2.0, 0.1))
    expect_match(ly$provenance$log2fc_source, "^sign\\(linearFC")
})

test_that("a table with only the rounded linear FC is converted, and says so", {
    df <- prot_summary()
    df$log2FC.imputs.A_vs_B <- NULL
    df$linearRatio.imputs.A_vs_B <- NULL
    path <- write_table_tmp(df)
    ly <- read_de_layer(layer_cfg(path), "A_vs_B", NULL, hits_default)
    expect_equal(ly$tables$A_vs_B$log2fc, c(log2(2.3), -2, log2(1.07)))
    expect_match(ly$provenance$log2fc_source, "significant digits")
})

test_that("the hit flag is per contrast, never pass_any_contrast", {
    df <- merge(prot_summary("A_vs_B"), prot_summary("C_vs_D")[, c(
        "FeatureID", grep("C_vs_D", names(prot_summary("C_vs_D")), value = TRUE))],
        by = "FeatureID")
    df$pass.imputs.C_vs_D <- c(NA, NA, NA)
    path <- write_table_tmp(df)
    ly <- read_de_layer(layer_cfg(path), c("A_vs_B", "C_vs_D"), NULL, hits_default)
    expect_identical(sum(ly$tables$A_vs_B$hit), 2L)
    expect_identical(sum(ly$tables$C_vs_D$hit), 0L)
})

test_that("cutoffs decide when the table flag is switched off", {
    path <- write_table_tmp(prot_summary())
    h <- modifyList(hits_default, list(use_table_flag = FALSE, linear_fc_cutoff = 3))
    ly <- read_de_layer(layer_cfg(path), "A_vs_B", NULL, h)
    # |log2FC| >= log2(3) = 1.58: only P2 (-2.0) clears it.
    expect_identical(ly$tables$A_vs_B$hit, c(FALSE, TRUE, FALSE))
    expect_match(ly$provenance$hit_source, "padj <= 0.05")
})

test_that("a differently spelled contrast resolves by key", {
    path <- write_table_tmp(prot_summary("A_vs_B"))
    ly <- read_de_layer(layer_cfg(path), "A vs. B", NULL, hits_default)
    expect_identical(ly$provenance$contrast, "A_vs_B")
    expect_identical(ly$provenance$requested_contrast, "A vs. B")
})

test_that("a configured contrast the table lacks is an error, even when it holds one", {
    path <- write_table_tmp(prot_summary("A_vs_B"))
    expect_error(read_de_layer(layer_cfg(path), "Treated_vs_Control", NULL, hits_default),
                 "'Treated_vs_Control' cannot be resolved.*The table holds: A_vs_B")

    # Configured end to end: never swapped for the table's only contrast.
    a <- write_table_tmp(prot_summary("A_vs_B"))
    b <- write_table_tmp(prot_summary("A_vs_B"))
    config <- dei_config(list(layer_cfg(a), layer_cfg(b, name = "media")),
                         comparisons = list(list(name = "c1", members = list(
                             cells = "A_vs_B", media = "Treated_vs_Control"))))
    expect_error(mod_dei_load_layers(config), "Layer 'media'.*Treated_vs_Control")
})

test_that("an unknown contrast in a multi-contrast table names what is there", {
    df <- merge(prot_summary("A_vs_B"), prot_summary("C_vs_D")[, c(
        "FeatureID", grep("C_vs_D", names(prot_summary("C_vs_D")), value = TRUE))],
        by = "FeatureID")
    path <- write_table_tmp(df)
    expect_error(read_de_layer(layer_cfg(path), "E_vs_F", NULL, hits_default),
                 "A_vs_B, C_vs_D")
})

test_that("a final table gives observed counts from its own N.observed columns", {
    path <- write_table_tmp(prot_final())
    ly <- read_de_layer(layer_cfg(path, format = "proteomics_final"), "A_vs_B", NULL,
                        hits_default)
    tab <- ly$tables$A_vs_B
    expect_identical(tab$hit, c(TRUE, TRUE, FALSE))
    expect_identical(tab$n_obs_num, c(3L, 1L, 2L))
    expect_identical(tab$n_obs_den, c(3L, 3L, 2L))
    expect_identical(tab$well_observed, c(TRUE, FALSE, TRUE))
    expect_identical(ly$provenance$n_obs_source, "N.observed.A, N.observed.B")
})

test_that("the groups are read off the columns for the export's space-stripped labels", {
    # "S vs NS" is written as SvsNS: nothing to split on, so the groups come
    # from the N.observed.<group> columns the table carries.
    path <- write_table_tmp(prot_final("SvsNS", groups = c("S", "NS")))
    ly <- read_de_layer(layer_cfg(path, format = "proteomics_final"), "SvsNS", NULL,
                        hits_default)
    tab <- ly$tables$SvsNS
    expect_identical(tab$n_obs_num, c(3L, 1L, 2L))
    expect_identical(tab$n_obs_den, c(3L, 3L, 2L))
    expect_identical(ly$provenance$n_obs_source, "N.observed.S, N.observed.NS")
    expect_identical(contrast_groups("SvsNS", groups = c("NS", "S"))[1:2],
                     list(numerator = "S", denominator = "NS"))
})

test_that("an old n_obs. prefix is not part of the native contract", {
    df <- prot_summary()
    df$n_obs.A <- c(3, 1, 2)
    df$n_obs.B <- c(3, 3, 2)
    ly <- read_de_layer(layer_cfg(write_table_tmp(df)), "A_vs_B", NULL, hits_default)
    expect_true(all(is.na(ly$tables$A_vs_B$n_obs_num)))
    expect_identical(ly$provenance$n_obs_source, "not available")
})

test_that("groups that fit a contrast more than one way are not guessed", {
    df <- prot_final("A_vs_B", groups = c("A", "B"))
    df$N.observed.a <- c(1, 1, 1)
    df$N.observed.b <- c(1, 1, 1)
    path <- write_table_tmp(df)
    expect_warning(
        ly <- read_de_layer(layer_cfg(path, format = "proteomics_final"), "A_vs_B", NULL,
                            hits_default),
        "ambiguous")
    expect_true(all(is.na(ly$tables$A_vs_B$n_obs_num)))
    expect_true(all(is.na(ly$tables$A_vs_B$well_observed)))
})

test_that("the contrasts file names the groups before the columns do", {
    # The contrasts file says the contrast runs NS over S; the columns alone
    # would have read it the other way.
    contr <- data.frame(Contrast_name = "S vs NS", Factor = "Group",
                        Numerator = "NS", Denominator = "S")
    path <- write_table_tmp(prot_final("SvsNS", groups = c("S", "NS")))
    ly <- read_de_layer(layer_cfg(path, format = "proteomics_final",
                                  observed = obs_block(contr, groups = c("S", "S", "NS", "NS"))),
                        "SvsNS", NULL, hits_default)
    expect_identical(ly$tables$SvsNS$n_obs_num, c(3L, 3L, 2L))
    expect_identical(ly$provenance$n_obs_source, "N.observed.NS, N.observed.S")
})

test_that("observed counts come from the unimputed matrix when the table has none", {
    path <- write_table_tmp(prot_summary())
    ly <- read_de_layer(layer_cfg(path, observed = obs_block()), "A_vs_B", NULL, hits_default)
    tab <- ly$tables$A_vs_B
    expect_identical(tab$n_obs_num, c(2L, 0L, 2L))
    expect_identical(tab$n_obs_den, c(2L, 2L, 1L))
    expect_identical(tab$well_observed, c(TRUE, FALSE, FALSE))
    expect_match(ly$provenance$n_obs_source, "unimputed matrix")
})

test_that("a generic table without mapped counts falls through to its observed block", {
    ext <- data.frame(prot = c("P1", "P2;P9", "P3"), logFC = c(2, -0.1, -3),
                      P = c(0.001, 0.5, 0.004), stringsAsFactors = FALSE)
    path <- write_table_tmp(ext, "csv")
    ly <- read_de_layer(layer_cfg(path, format = "generic", contrast = "A_vs_B",
        columns = list(id = "prot", log2fc = "logFC", pvalue = "P"),
        observed = obs_block()), "A_vs_B", NULL, hits_default)
    tab <- ly$tables$A_vs_B
    expect_identical(tab$n_obs_num, c(2L, 0L, 2L))
    expect_identical(tab$n_obs_den, c(2L, 2L, 1L))
    expect_match(ly$provenance$n_obs_source, "unimputed matrix")
})

test_that("a generic table's mapped counts are used, and must exist", {
    ext <- data.frame(prot = c("X1", "X2"), logFC = c(1, -1), P = c(0.01, 0.2),
                      nA = c(3, 1), nB = c(2, 2))
    path <- write_table_tmp(ext, "csv")
    cols <- list(id = "prot", log2fc = "logFC", pvalue = "P", n_obs_num = "nA",
                 n_obs_den = "nB")
    ly <- read_de_layer(layer_cfg(path, format = "generic", contrast = "T_vs_C",
                                  columns = cols), "T_vs_C", NULL, hits_default)
    expect_identical(ly$tables$T_vs_C$n_obs_num, c(3L, 1L))
    expect_identical(ly$provenance$n_obs_source, "nA, nB")

    cols$n_obs_den <- "nC"
    expect_error(read_de_layer(layer_cfg(path, format = "generic", contrast = "T_vs_C",
                                         columns = cols), "T_vs_C", NULL, hits_default),
                 "Layer 'cells': columns.n_obs_den names 'nC', not found")
})

test_that("an optional column the map names must exist, not quietly fall back", {
    ext <- data.frame(prot = c("X1", "X2"), logFC = c(1, -1), P = c(0.01, 0.2))
    path <- write_table_tmp(ext, "csv")
    ly <- layer_cfg(path, format = "generic", contrast = "T_vs_C",
                    columns = list(id = "prot", log2fc = "logFC", pvalue = "P",
                                   padj = "adj.P", hit = "sig"))
    expect_error(read_de_layer(ly, "T_vs_C", NULL, hits_default),
                 "columns.padj names 'adj.P', columns.hit names 'sig', not found")
})

test_that("an observed block that cannot serve its purpose stops the run", {
    path <- write_table_tmp(prot_summary())
    read_with <- function(obs) {
        read_de_layer(layer_cfg(path, observed = obs), "A_vs_B", NULL, hits_default)
    }

    obs <- obs_block(); obs$sample_col <- "Sample"
    expect_error(read_with(obs), "Layer 'cells': observed.samplesheet .*no column 'Sample'")

    obs <- obs_block(data.frame(Contrast_name = "A_vs_B", Factor = "Group", Numerator = "A"))
    expect_error(read_with(obs), "observed.contrasts_file .*lacks column\\(s\\): Denominator")

    obs <- obs_block(data.frame(Contrast_name = "A_vs_B", Factor = "Treatment",
                                Numerator = "A", Denominator = "B"))
    expect_error(read_with(obs), "observed.contrasts_file .*groups by Treatment")

    obs <- obs_block(); obs$samplesheet <- file.path(tempdir(), "no_such_sheet.csv")
    expect_error(read_with(obs), "observed.samplesheet .*cannot be read")

    obs <- obs_block(); obs$id_col <- "FeatureID"
    expect_error(read_with(obs), "observed.matrix .*no column 'FeatureID'")

    sheet <- data.frame(SampleName = c("x1", "x2"), Group = c("A", "B"))
    obs <- obs_block(); obs$samplesheet <- write_table_tmp(sheet, "csv")
    expect_error(read_with(obs), "observed.matrix .*no column named after a sample")
})

test_that("a contrasts-file group that no counted sample carries stops the run", {
    path <- write_table_tmp(prot_summary())
    contr <- data.frame(Contrast_name = "A_vs_B", Factor = "Group",
                        Numerator = "Treatd", Denominator = "B")
    expect_error(read_de_layer(layer_cfg(path, observed = obs_block(contr)), "A_vs_B", NULL,
                               hits_default),
                 "Layer 'cells': observed.contrasts_file .*group\\(s\\) Treatd for contrast 'A_vs_B'.*'Group' column")
})

test_that("p-values outside [0, 1], or not numbers, stop the run", {
    read_p <- function(p, padj = NULL) {
        ext <- data.frame(prot = c("X1", "X2"), logFC = c(2, -2), P = p,
                          stringsAsFactors = FALSE)
        cols <- list(id = "prot", log2fc = "logFC", pvalue = "P")
        if (!is.null(padj)) { ext$Q <- padj; cols$padj <- "Q" }
        read_de_layer(layer_cfg(write_table_tmp(ext, "csv"), format = "generic",
                                contrast = "T_vs_C", columns = cols),
                      "T_vs_C", NULL, hits_default)
    }
    expect_error(read_p(c(-0.01, 0.5)), "Layer 'cells' .*column 'P' must hold p-values")
    expect_error(read_p(c(0.01, 1.5)), "column 'P' must hold p-values")
    expect_error(.dei_numeric_col(data.frame(P = c(0.01, Inf)), "P", "Layer 'x'", 0, 1,
                                  finite = TRUE, what = "p-values"),
                 "Layer 'x': column 'P' must hold p-values; it has values that are not finite")
    expect_error(read_p(c("0.01", "low")), "column 'P' must hold p-values but is not numeric")
    expect_error(read_p(c(0.01, 0.5), padj = c(0.02, 2)),
                 "column 'Q' must hold adjusted p-values")
    # Missing values stay missing rather than stopping the run.
    expect_identical(read_p(c(0.01, NA))$tables$T_vs_C$hit, c(TRUE, FALSE))
})

test_that("a fold-change column that is not numeric stops the run", {
    ext <- data.frame(prot = c("X1", "X2"), logFC = c("up", "down"), P = c(0.01, 0.2))
    expect_error(read_de_layer(layer_cfg(write_table_tmp(ext, "csv"), format = "generic",
                                         contrast = "T_vs_C",
                                         columns = list(id = "prot", log2fc = "logFC",
                                                        pvalue = "P")),
                               "T_vs_C", NULL, hits_default),
                 "column 'logFC' must hold fold changes but is not numeric")
})

test_that("mapped observed counts must be whole, non-negative numbers", {
    counts <- function(nA) {
        ext <- data.frame(prot = c("X1", "X2"), logFC = c(1, -1), P = c(0.01, 0.2),
                          nA = nA, nB = c(2, 2))
        read_de_layer(layer_cfg(write_table_tmp(ext, "csv"), format = "generic",
                                contrast = "T_vs_C",
                                columns = list(id = "prot", log2fc = "logFC", pvalue = "P",
                                               n_obs_num = "nA", n_obs_den = "nB")),
                      "T_vs_C", NULL, hits_default)
    }
    expect_error(counts(c(-1, 2)), "column 'nA' must hold counts of observed values")
    expect_error(counts(c(1.5, 2)), "column 'nA' must hold whole counts")
})

test_that("observed inputs with repeated or blank ids stop the run", {
    path <- write_table_tmp(prot_summary())
    read_with <- function(obs) {
        read_de_layer(layer_cfg(path, observed = obs), "A_vs_B", NULL, hits_default)
    }

    obs <- obs_block()
    sheet <- data.frame(SampleName = c("s1", "s1", "s3", "s4"), Group = c("A", "B", "B", "B"))
    obs$samplesheet <- write_table_tmp(sheet, "csv")
    expect_error(read_with(obs), "observed.samplesheet .*repeats sample id\\(s\\) s1")

    obs <- obs_block(ids = c("P1", "P1", "P3"))
    expect_error(read_with(obs), "observed.matrix .*repeats feature id\\(s\\) P1")

    obs <- obs_block(data.frame(Contrast_name = c("A_vs_B", "A vs B"), Factor = "Group",
                                Numerator = c("A", "B"), Denominator = c("B", "A")))
    expect_error(read_with(obs), "observed.contrasts_file .*names one contrast more than once")
})

test_that("an annotation giving one gene id two symbols stops the run", {
    rna <- data.frame(Gene = c("ENSG1", "ENSG2"),
                      log2FC.A_vs_B = c(1, -1), pvalue.A_vs_B = c(0.01, 0.02),
                      padj.A_vs_B = c(0.04, 0.04), stringsAsFactors = FALSE)
    path <- write_table_tmp(rna)
    ann <- write_table_tmp(data.frame(gene_id = c("ENSG1", "ENSG1", "ENSG2"),
                                      symbol = c("G1", "G1b", "G2")), "csv")
    expect_error(read_de_layer(layer_cfg(path, format = "rnaseq_summary", omics = "rnaseq",
                                         annotation_file = ann), "A_vs_B", NULL, hits_default),
                 "more than one symbol for gene id\\(s\\) ENSG1")
    # The same pairing written twice is not a conflict.
    ann <- write_table_tmp(data.frame(gene_id = c("ENSG1", "ENSG1"), symbol = c("G1", "G1")),
                           "csv")
    ly <- read_de_layer(layer_cfg(path, format = "rnaseq_summary", omics = "rnaseq",
                                  annotation_file = ann), "A_vs_B", NULL, hits_default)
    expect_identical(ly$tables$A_vs_B$symbol, c("G1", "ENSG2"))
})

test_that("BH runs over the rows the table keeps, not the ones it drops", {
    ext <- data.frame(prot = c("X1", "X2", "X2", NA), logFC = c(2, 2, 2, 2),
                      P = c(0.01, 0.02, 0.5, 0.03), stringsAsFactors = FALSE)
    expect_warning(
        ly <- read_de_layer(layer_cfg(write_table_tmp(ext, "csv"), format = "generic",
                                      contrast = "T_vs_C",
                                      columns = list(id = "prot", log2fc = "logFC", pvalue = "P")),
                            "T_vs_C", NULL, hits_default),
        "1 duplicated feature id")
    tab <- ly$tables$T_vs_C
    expect_identical(tab$feature_id, c("X1", "X2"))
    expect_equal(tab$pvalue, c(0.01, 0.02))
    expect_equal(tab$padj, p.adjust(c(0.01, 0.02), "BH"))
})

test_that("a signed linear fold change between -1 and 1 stops the run", {
    ext <- data.frame(prot = c("X1", "X2"), FC = c(2, -0.5), P = c(0.01, 0.02))
    expect_error(read_de_layer(layer_cfg(write_table_tmp(ext, "csv"), format = "generic",
                                         contrast = "T_vs_C",
                                         columns = list(id = "prot", linear_fc = "FC",
                                                        pvalue = "P")),
                               "T_vs_C", NULL, hits_default),
                 "column 'FC' must hold signed linear fold changes")
})

test_that("fold changes must be finite, and ratios positive", {
    read_prot <- function(df) {
        read_de_layer(layer_cfg(write_table_tmp(df)), "A_vs_B", NULL, hits_default)
    }
    df <- prot_summary(); df$log2FC.imputs.A_vs_B <- NULL
    df$linearRatio.imputs.A_vs_B <- c(0, 4, 1.07)
    expect_error(read_prot(df), "column 'linearRatio.imputs.A_vs_B' must hold positive")
    df$linearRatio.imputs.A_vs_B <- c(-2, 4, 1.07)
    expect_error(read_prot(df), "must hold positive linear ratios")
    # The sign source is held to the signed convention too.
    df <- prot_summary(); df$log2FC.imputs.A_vs_B <- NULL
    df$linearFC.imputs.A_vs_B <- c(2.3, -0.25, 1.07)
    expect_error(read_prot(df), "linearFC.imputs.A_vs_B' must hold signed linear fold changes")
    expect_error(.dei_numeric_col(data.frame(L = c(1, -Inf)), "L", "Layer 'x'", finite = TRUE,
                                  what = "fold changes"),
                 "must hold fold changes; it has values that are not finite")
})

test_that("an annotation row without a symbol neither clashes nor hides one", {
    rna <- data.frame(Gene = c("ENSG1", "ENSG2"),
                      log2FC.A_vs_B = c(1, -1), pvalue.A_vs_B = c(0.01, 0.02),
                      padj.A_vs_B = c(0.04, 0.04), stringsAsFactors = FALSE)
    path <- write_table_tmp(rna)
    # The blank row comes first, so a first-match lookup would find no symbol.
    ann <- write_table_tmp(data.frame(gene_id = c("ENSG1", "ENSG1"), symbol = c(NA, "G1")),
                           "csv")
    ly <- read_de_layer(layer_cfg(path, format = "rnaseq_summary", omics = "rnaseq",
                                  annotation_file = ann), "A_vs_B", NULL, hits_default)
    expect_identical(ly$tables$A_vs_B$symbol, c("G1", "ENSG2"))
    expect_identical(ly$tables$A_vs_B$symbol_source, c("symbol", "gene_id"))
})

test_that("a repeated column name in any input stops the run", {
    write_raw <- function(lines, ext) {
        f <- tempfile(fileext = paste0(".", ext))
        writeLines(lines, f)
        f
    }
    path <- write_table_tmp(prot_summary())
    obs <- obs_block()
    obs$matrix <- write_raw(c("Protein.Group\ts1\ts2\ts1\ts4",
                              "P1\t20\t21\t19\t20", "P2;P9\tNA\tNA\t22\t23",
                              "P3\t15\t15\tNA\t16"), "tsv")
    expect_error(read_de_layer(layer_cfg(path, observed = obs), "A_vs_B", NULL, hits_default),
                 "Layer 'cells': observed.matrix .*repeated column name \\(s1\\)")

    # The DE table itself: two p-value columns for one contrast.
    df <- prot_summary()
    lines <- c(paste(c(names(df), "pvalue.imputs.A_vs_B"), collapse = "\t"),
               apply(cbind(df, 0.5), 1, paste, collapse = "\t"))
    expect_error(read_de_layer(layer_cfg(write_raw(lines, "tsv")), "A_vs_B", NULL, hits_default),
                 "Layer 'cells': path .*repeated column name \\(pvalue.imputs.A_vs_B\\)")

    # A blank header cell.
    obs <- obs_block()
    obs$matrix <- write_raw(c("Protein.Group\ts1\t\ts3\ts4",
                              "P1\t20\t21\t19\t20", "P2;P9\tNA\tNA\t22\t23",
                              "P3\t15\t15\tNA\t16"), "tsv")
    expect_error(read_de_layer(layer_cfg(path, observed = obs), "A_vs_B", NULL, hits_default),
                 "observed.matrix .*a blank column name")
})

test_that("a unique column name that merely ends in ...<n> is accepted", {
    # Checked on the header as written, not on the names readr hands back.
    f <- tempfile(fileext = ".csv")
    writeLines(c("sample...1,sample...2", "1,2"), f)
    df <- data.frame(a = 1, b = 2)
    expect_silent(.dei_check_table(df, f, "Layer 'x'"))
    writeLines(c("s1,s1", "1,2"), f)
    expect_error(.dei_check_table(df, f, "Layer 'x'"), "repeated column name \\(s1\\)")
})

test_that("a malformed value past readr's type guess stops the run", {
    # readr guesses the column type from its first 1000 rows; a bad value
    # later becomes NA with only a warning.
    n <- 1100
    ext <- data.frame(prot = paste0("X", seq_len(n)), logFC = rep(1, n),
                      P = as.character(rep(0.01, n)), stringsAsFactors = FALSE)
    ext$P[n] <- "oops"
    path <- write_table_tmp(ext, "csv")
    suppressWarnings(expect_error(
        read_de_layer(layer_cfg(path, format = "generic", contrast = "T_vs_C",
                                columns = list(id = "prot", log2fc = "logFC", pvalue = "P")),
                      "T_vs_C", NULL, hits_default),
        "Layer 'cells': path .*column 'P' has values that are not numbers \\(e.g. 'oops'\\)"))
})

test_that("an annotation or observed matrix that matches no feature id stops the run", {
    rna <- data.frame(Gene = c("ENSG1", "ENSG2"),
                      log2FC.A_vs_B = c(1, -1), pvalue.A_vs_B = c(0.01, 0.02),
                      padj.A_vs_B = c(0.04, 0.04), stringsAsFactors = FALSE)
    ann <- write_table_tmp(data.frame(gene_id = c("ENSG1.5", "ENSG2.5"),
                                      symbol = c("G1", "G2")), "csv")
    expect_error(read_de_layer(layer_cfg(write_table_tmp(rna), format = "rnaseq_summary",
                                         omics = "rnaseq", annotation_file = ann),
                               "A_vs_B", NULL, hits_default),
                 "annotation_file matches none of the layer's feature ids")

    obs <- obs_block(ids = c("Q1", "Q2", "Q3"))
    expect_error(read_de_layer(layer_cfg(write_table_tmp(prot_summary()), observed = obs),
                               "A_vs_B", NULL, hits_default),
                 "observed.matrix .*matches none of the DE table's feature ids")
})

test_that("an observed matrix with infinite values stops the run", {
    path <- write_table_tmp(prot_summary())
    obs <- obs_block()
    mat <- data.frame(Protein.Group = c("P1", "P2;P9", "P3"), s1 = c(20, -Inf, 15),
                      s2 = c(21, NA, 15), s3 = c(19, 22, NA), s4 = c(20, 23, 16))
    # readr reads "-Inf" back as a number, as a log of zero would be written.
    obs$matrix <- write_table_tmp(mat)
    expect_error(read_de_layer(layer_cfg(path, observed = obs), "A_vs_B", NULL, hits_default),
                 "observed.matrix .*infinite values in sample column\\(s\\) s1")
})

test_that("a contrasts row with one group on both sides, or a blank one, stops the run", {
    path <- write_table_tmp(prot_summary())
    read_with <- function(contr) {
        read_de_layer(layer_cfg(path, observed = obs_block(contr)), "A_vs_B", NULL,
                      hits_default)
    }
    expect_error(read_with(data.frame(Contrast_name = "A_vs_B", Factor = "Group",
                                      Numerator = "A", Denominator = "A")),
                 "same group as Numerator and Denominator for A_vs_B")
    expect_error(read_with(data.frame(Contrast_name = "A_vs_B", Factor = "Group",
                                      Numerator = "A", Denominator = NA)),
                 "blank Numerator or Denominator")
})

test_that("a contrasts file with a blank Contrast_name stops the run", {
    path <- write_table_tmp(prot_summary())
    obs <- obs_block(data.frame(Contrast_name = c("A_vs_B", NA), Factor = "Group",
                                Numerator = c("A", "B"), Denominator = c("B", "A")))
    expect_error(read_de_layer(layer_cfg(path, observed = obs), "A_vs_B", NULL, hits_default),
                 "observed.contrasts_file .*blank or missing Contrast_name")
})

test_that("a contrast missing from the contrasts file leaves counts unavailable, with a warning", {
    path <- write_table_tmp(prot_summary())
    contr <- data.frame(Contrast_name = "C_vs_D", Factor = "Group",
                        Numerator = "A", Denominator = "B")
    expect_warning(
        ly <- read_de_layer(layer_cfg(path, observed = obs_block(contr)), "A_vs_B", NULL,
                            hits_default),
        "'A_vs_B' is not in observed.contrasts_file")
    expect_identical(ly$provenance$n_obs_source, "not available")
})

test_that("a native export can name its own id column", {
    df <- prot_summary()
    names(df)[names(df) == "FeatureID"] <- "Protein.Group"
    path <- write_table_tmp(df)
    expect_error(read_de_layer(layer_cfg(path), "A_vs_B", NULL, hits_default),
                 "no feature id column \\(looked for FeatureID\\)")
    ly <- read_de_layer(layer_cfg(path, id_col = "Protein.Group"), "A_vs_B", NULL,
                        hits_default)
    expect_identical(ly$tables$A_vs_B$feature_id, c("P1", "P2;P9", "P3"))
    expect_identical(ly$tables$A_vs_B$symbol, c("G1", "G2", NA))
})

test_that("features are ordered the same way under any locale", {
    df <- prot_summary()
    df$FeatureID <- c("b", "B", "a")
    df <- rbind(df, df[3, ])
    df$FeatureID[4] <- "A"
    ly <- read_de_layer(layer_cfg(write_table_tmp(df)), "A_vs_B", NULL, hits_default)
    # Bytewise: upper case before lower case, whatever the collation.
    expect_identical(ly$tables$A_vs_B$feature_id, c("A", "B", "a", "b"))
})

test_that("an RNA layer uses the gene id as its key unless an annotation names symbols", {
    rna <- data.frame(Gene = c("ENSG1", "ENSG2"),
                      log2FC.A_vs_B = c(1, -1), linearFC.A_vs_B = c(2, -2),
                      pvalue.A_vs_B = c(0.01, 0.02), padj.A_vs_B = c(0.04, 0.04),
                      A_vs_B_pass = c(1, 0), stringsAsFactors = FALSE)
    path <- write_table_tmp(rna)
    ly <- read_de_layer(layer_cfg(path, format = "rnaseq_summary", omics = "rnaseq"),
                        "A_vs_B", NULL, hits_default)
    tab <- ly$tables$A_vs_B
    expect_identical(tab$symbol, c("ENSG1", "ENSG2"))
    expect_identical(tab$symbol_source, c("gene_id", "gene_id"))
    expect_identical(tab$hit, c(TRUE, FALSE))

    ann <- write_table_tmp(data.frame(gene_id = c("ENSG1"), symbol = c("G1")), "csv")
    ly <- read_de_layer(layer_cfg(path, format = "rnaseq_summary", omics = "rnaseq",
                                  annotation_file = ann), "A_vs_B", NULL, hits_default)
    tab <- ly$tables$A_vs_B
    expect_identical(tab$symbol, c("G1", "ENSG2"))
    expect_identical(tab$symbol_source, c("symbol", "gene_id"))

    bad <- write_table_tmp(data.frame(id = c("ENSG1"), name = c("G1")), "csv")
    expect_error(read_de_layer(layer_cfg(path, format = "rnaseq_summary", omics = "rnaseq",
                                         annotation_file = bad), "A_vs_B", NULL, hits_default),
                 "annotation_file .*lacks column\\(s\\): gene_id, symbol")
})

test_that("a generic table is read through its column map", {
    ext <- data.frame(prot = c("X1", "X2", "X3"), logFC = c(2, -0.1, -3),
                      P = c(0.001, 0.5, 0.004), gid = c("g1", "g2", "g3"))
    path <- write_table_tmp(ext, "csv")
    ly <- read_de_layer(layer_cfg(path, format = "generic", contrast = "T_vs_C",
        columns = list(id = "prot", log2fc = "logFC", pvalue = "P", gene_id = "gid")),
        "T_vs_C", NULL, hits_default)
    tab <- ly$tables$T_vs_C
    expect_equal(tab$padj, p.adjust(c(0.001, 0.5, 0.004), "BH"))
    expect_identical(ly$provenance$padj_source, "BH of P")
    expect_identical(tab$hit, c(TRUE, FALSE, TRUE))
    expect_identical(tab$symbol_source, rep("gene_id", 3))
})

test_that("duplicated ids keep the row with the smallest p-value", {
    df <- prot_summary()
    df <- rbind(df, df[1, ])
    df$pvalue.imputs.A_vs_B[4] <- 1e-6
    df$log2FC.imputs.A_vs_B[4] <- 9
    path <- write_table_tmp(df)
    expect_warning(ly <- read_de_layer(layer_cfg(path), "A_vs_B", NULL, hits_default),
                   "1 duplicated feature id")
    tab <- ly$tables$A_vs_B
    expect_identical(nrow(tab), 3L)
    expect_equal(tab$log2fc[tab$feature_id == "P1"], 9)
    expect_identical(ly$provenance$n_duplicates_dropped, 1L)
})

test_that("a missing column is named in the error", {
    ext <- data.frame(prot = c("X1", "X2"), logFC = c(1, -1), Q = c(0.01, 0.2))
    path <- write_table_tmp(ext, "csv")
    ly <- layer_cfg(path, format = "generic", contrast = "T_vs_C",
                    columns = list(id = "prot", log2fc = "logFC", pvalue = "P"))
    expect_error(read_de_layer(ly, "T_vs_C", NULL, hits_default), "column 'P' not found")

    ly$columns <- list(id = "prot", log2fc = "FC", pvalue = "Q")
    expect_error(read_de_layer(ly, "T_vs_C", NULL, hits_default),
                 "No fold-change column found; looked for: FC")
})

test_that("contrasts are paired across layers by key, or as lone contrasts", {
    cfg <- list(comparisons = list())
    comps <- resolve_dei_comparisons(cfg, list(cells = c("A_vs_B", "C_vs_D"),
                                               media = c("A vs. B")))
    expect_length(comps, 1)
    expect_identical(unname(comps[[1]]$members), c("A_vs_B", "A vs. B"))
    expect_identical(comps[[1]]$name, "A_vs_B")

    expect_message(comps <- resolve_dei_comparisons(cfg, list(cells = "X_vs_Y",
                                                              media = "Treated_vs_Control")),
                   "pairing each layer's only contrast")
    expect_identical(names(comps[[1]]$members), c("cells", "media"))

    expect_error(resolve_dei_comparisons(cfg, list(cells = c("A_vs_B", "C_vs_D"),
                                                   media = "E_vs_F")),
                 "comparisons")
})

test_that("a layer with two contrasts on one key is not quietly left out", {
    cfg <- list(comparisons = list())
    expect_error(resolve_dei_comparisons(cfg, list(cells = c("A_vs_B", "A vs. B"),
                                                   media = "A_vs_B", tissue = "A - B")),
                 "Layer 'cells' holds contrasts that name the same comparison \\(A_vs_B, A vs. B\\)")
})

test_that("configured comparisons are taken as written", {
    cfg <- list(comparisons = list(list(name = "mine",
        members = list(cells = "A_vs_B", media = "Z_vs_W"), flip = "media")))
    comps <- resolve_dei_comparisons(cfg, list(cells = "A_vs_B", media = "Z_vs_W"))
    expect_identical(comps[[1]]$name, "mine")
    expect_identical(comps[[1]]$flip, "media")
})

test_that("every input file is tracked, and a missing one stops the run", {
    a <- write_table_tmp(prot_summary())
    b <- write_table_tmp(prot_summary())
    config <- dei_config(list(layer_cfg(a), layer_cfg(b, name = "media")))
    expect_setequal(de_integration_input_files(config), c(a, b))

    config$modes$de_integration$layers[[2]]$annotation_file <- file.path(tempdir(), "nope.csv")
    expect_error(de_integration_input_files(config), "annotation_file not found")
})

test_that("two layers load end to end and the summary is written", {
    a <- write_table_tmp(prot_summary())
    # Proteomics exports strip spaces from contrast names; a dotted spelling of
    # the same contrast still pairs by key.
    b <- write_table_tmp(prot_summary("A.vs.B"))
    config <- dei_config(list(layer_cfg(a), layer_cfg(b, name = "media")))
    res <- mod_dei_load_layers(config)
    expect_named(res$layers, c("cells", "media"))
    expect_length(res$comparisons, 1)
    expect_identical(nrow(res$summary), 2L)

    out <- withr::local_tempdir()
    files <- write_dei_layer_summary(res, out)
    expect_true(all(file.exists(files)))
    expect_true(all(c("layer_summary.tsv", "comparisons.tsv") %in% basename(files)))
    comps <- utils::read.delim(files[basename(files) == "comparisons.tsv"])
    expect_true("flip_requested" %in% names(comps))
    expect_false("flipped" %in% names(comps))
})

test_that("the pipeline tracks the inputs and loads layers through them", {
    skip_if_not_installed("targets")
    suppressPackageStartupMessages(library(targets))
    tl <- pipe_de_integration()
    nm <- vapply(tl, function(t) t$settings$name, character(1))
    expect_identical(nm, c("dei_out_dir", "dei_input_files", "dei_layers",
                           "dei_layer_summary_file"))
    by <- stats::setNames(tl, nm)
    expect_identical(by$dei_input_files$settings$format, "file")
    expect_identical(by$dei_layer_summary_file$settings$format, "file")
    deps <- by$dei_layers$command$deps
    if (length(deps) == 0) deps <- unlist(lapply(as.list(by$dei_layers$command$expr), all.vars))
    expect_true("dei_input_files" %in% deps)
})
