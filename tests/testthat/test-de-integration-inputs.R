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

test_that("a proteomics summary is read with its stored log2FC and its own hit flag", {
    path <- write_table_tmp(prot_summary())
    ly <- read_de_layer(layer_cfg(path), "A_vs_B", dei_config(list(
        layer_cfg(path), layer_cfg(path, name = "media"))), hits_default)
    tab <- ly$tables$A_vs_B
    expect_identical(tab$feature_id, c("P1", "P2;P9", "P3"))
    expect_equal(tab$log2fc, c(1.2, -2.0, 0.1))
    expect_identical(tab$hit, c(TRUE, TRUE, FALSE))
    expect_identical(ly$provenance$log2fc_source, "log2FC.imputs.A_vs_B")
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

test_that("a differently spelled contrast resolves by key, and a lone one by default", {
    path <- write_table_tmp(prot_summary("A_vs_B"))
    ly <- read_de_layer(layer_cfg(path), "A vs. B", NULL, hits_default)
    expect_identical(ly$provenance$contrast, "A_vs_B")
    expect_message(
        ly <- read_de_layer(layer_cfg(path), "Treated_vs_Control", NULL, hits_default),
        "only contrast")
    expect_identical(ly$provenance$contrast, "A_vs_B")
})

test_that("an unknown contrast in a multi-contrast table names what is there", {
    df <- merge(prot_summary("A_vs_B"), prot_summary("C_vs_D")[, c(
        "FeatureID", grep("C_vs_D", names(prot_summary("C_vs_D")), value = TRUE))],
        by = "FeatureID")
    path <- write_table_tmp(df)
    expect_error(read_de_layer(layer_cfg(path), "E_vs_F", NULL, hits_default),
                 "A_vs_B, C_vs_D")
})

test_that("a final table gives observed counts from its own columns", {
    df <- prot_summary()
    names(df)[names(df) == "pass.imputs.A_vs_B"] <- "upDown.imputs.A_vs_B"
    df$upDown.imputs.A_vs_B <- c("Up", "Down", NA)
    df$n_obs.A <- c(3, 1, 2)
    df$n_obs.B <- c(3, 3, 2)
    path <- write_table_tmp(df)
    ly <- read_de_layer(layer_cfg(path, format = "proteomics_final"), "A_vs_B", NULL,
                        hits_default)
    tab <- ly$tables$A_vs_B
    expect_identical(tab$hit, c(TRUE, TRUE, FALSE))
    expect_identical(tab$n_obs_num, c(3L, 1L, 2L))
    expect_identical(tab$well_observed, c(TRUE, FALSE, TRUE))
})

test_that("observed counts come from the unimputed matrix when the table has none", {
    mat <- data.frame(Protein.Group = c("P1", "P2;P9", "P3"),
                      s1 = c(20, NA, 15), s2 = c(21, NA, 15),
                      s3 = c(19, 22, NA), s4 = c(20, 23, 16),
                      stringsAsFactors = FALSE)
    sheet <- data.frame(SampleName = c("s1", "s2", "s3", "s4"),
                        Group = c("A", "A", "B", "B"))
    contr <- data.frame(Contrast_name = "A_vs_B", Factor = "Group",
                        Numerator = "A", Denominator = "B")
    obs <- list(matrix = write_table_tmp(mat), samplesheet = write_table_tmp(sheet, "csv"),
                sample_col = "SampleName", contrasts_file = write_table_tmp(contr, "csv"))
    path <- write_table_tmp(prot_summary())
    ly <- read_de_layer(layer_cfg(path, observed = obs), "A_vs_B", NULL, hits_default)
    tab <- ly$tables$A_vs_B
    expect_identical(tab$n_obs_num, c(2L, 0L, 2L))
    expect_identical(tab$n_obs_den, c(2L, 2L, 1L))
    expect_identical(tab$well_observed, c(TRUE, FALSE, FALSE))
    expect_match(ly$provenance$n_obs_source, "unimputed matrix")
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
