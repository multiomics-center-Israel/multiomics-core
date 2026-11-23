# tests/testthat/test-rna-input-source.R
#
# The RNA report's gene counts and its Input Data section describe the input
# the run actually used. With both `files$counts` and `files$txi` configured,
# the run builds on one of them; the report used to re-read `files$counts` and
# so could show the gene count of a file the run never used.
#
# All data here is synthetic.

# ---- Config-key helpers -----------------------------------------------------

test_that("rna_configured_expression_inputs lists the keys holding a path", {
    expect_identical(
        rna_configured_expression_inputs(list(counts = "rna/c.csv", txi = "rna/txi.rds",
                                              metadata = "rna/meta.csv")),
        c("counts", "txi"))
    expect_identical(rna_configured_expression_inputs(list(counts = "", txi = "rna/txi.rds")),
                     "txi")
    expect_identical(rna_configured_expression_inputs(list(counts = NULL, txi = NA_character_,
                                                            preprocessed_counts = "p.tsv")),
                     "preprocessed_counts")
    expect_identical(rna_configured_expression_inputs(list()), character(0))
})

test_that("rna_rds_input_key names the key whose .rds file was kept as txi", {
    expect_identical(rna_rds_input_key(list(counts = "rna/c.csv", txi = "rna/txi.rds")), "txi")
    # load_omics_inputs() stores any .rds under `txi`, whatever its key.
    expect_identical(rna_rds_input_key(list(counts = "rna/c.RDS")), "counts")
    expect_identical(rna_rds_input_key(list(counts = "rna/c.csv")), NA_character_)
})

# ---- describe_rna_input_source() --------------------------------------------

both_files <- list(counts = "rna/counts.csv", txi = "rna/txi.rds")

test_that("with counts and txi both set, the description names the one used", {
    src <- describe_rna_input_source(
        list(configured = c("counts", "txi"), txi = "txi", counts = "counts"),
        "tximport", both_files)

    expect_true(src$conflict)
    expect_identical(src$used_key, "txi")
    expect_identical(src$used_file, "rna/txi.rds")
    expect_identical(src$configured, c("counts", "txi"))
    expect_match(src$note, "`files$counts` (`rna/counts.csv`) and `files$txi` (`rna/txi.rds`)",
                 fixed = TRUE)
    expect_match(src$note, "This run used `files$txi` (`rna/txi.rds`), read as a tximport object",
                 fixed = TRUE)
    expect_match(src$note, "every gene count reported for it refers to that input only",
                 fixed = TRUE)
})

test_that("the input used follows the branch preprocess_rna took, not a key priority", {
    keys <- list(configured = c("counts", "txi"), txi = "txi", counts = "counts")
    # Preprocessed input reads the counts slot even when a txi is loaded.
    src <- describe_rna_input_source(keys, "preprocessed", both_files)
    expect_identical(src$used_key, "counts")
    expect_match(src$note, "This run used `files$counts`", fixed = TRUE)
    expect_match(src$description, "a preprocessed expression matrix", fixed = TRUE)
})

test_that("a single configured input is described without a note", {
    src <- describe_rna_input_source(
        list(configured = "counts", txi = NA_character_, counts = "counts"),
        "matrix", list(counts = "rna/counts.csv"))
    expect_false(src$conflict)
    expect_null(src$note)
    expect_identical(src$description, "`files$counts` (`rna/counts.csv`), read as a raw count matrix")
})

test_that("the preprocessed_counts fallback is reported as such", {
    src <- describe_rna_input_source(
        list(configured = "preprocessed_counts", txi = NA_character_,
             counts = "preprocessed_counts"),
        "matrix", list(preprocessed_counts = "rna/pre.tsv"))
    expect_identical(src$used_key, "preprocessed_counts")
    expect_false(src$conflict)
})

test_that("inputs assembled without load_rna_inputs() say the key is not recorded", {
    src <- describe_rna_input_source(NULL, "matrix", list(counts = "rna/counts.csv"))
    expect_identical(src$used_key, NA_character_)
    expect_false(src$conflict)
    expect_match(src$description, "not recorded", fixed = TRUE)
})

# ---- rna_recorded_gene_counts() ---------------------------------------------

test_that("gene counts come from the run's info, NA when not recorded", {
    expect_identical(rna_recorded_gene_counts(list(n_genes_raw = 120L, n_genes_filt = 80L)),
                     list(before = 120L, after = 80L))
    expect_identical(rna_recorded_gene_counts(NULL),
                     list(before = NA_integer_, after = NA_integer_))
    expect_identical(rna_recorded_gene_counts(list(n_genes_raw = 50))$after, NA_integer_)
})

# ---- load_rna_inputs() + preprocess_rna(): counts and txi both set ----------

# A project whose counts file has more genes than its tximport object, so the
# two inputs are told apart by their gene count.
stage_two_inputs <- function(n_counts = 300L, n_txi = 100L, seed = 7L) {
    proj <- withr::local_tempdir(.local_envir = parent.frame())
    raw <- file.path(proj, "data", "rna")
    generate_synthetic_rna(raw, seed = seed, n_genes = n_counts)

    samples <- c("CTRL_01", "CTRL_02", "CTRL_03", "TRT_01", "TRT_02", "TRT_03")
    genes <- sprintf("TX_GENE_%04d", seq_len(n_txi))
    counts <- withr::with_seed(seed, matrix(
        stats::rnbinom(n_txi * 6, mu = 200, size = 8), nrow = n_txi,
        dimnames = list(genes, samples)))
    storage.mode(counts) <- "numeric"
    txi <- list(
        counts    = counts,
        abundance = counts / 10,
        length    = matrix(1000, nrow = n_txi, ncol = 6, dimnames = dimnames(counts)),
        countsFromAbundance = "no"
    )
    saveRDS(txi, file.path(raw, "txi.rds"))
    proj
}

two_input_config <- function(proj, files) {
    list(
        project = list(dir = proj, name = "input_source_test"),
        paths   = list(raw = "data", out = file.path(proj, "out")),
        modes   = list(rna = list(
            files         = c(files, list(metadata = "rna/metadata.csv",
                                          contrasts = "rna/contrasts.csv")),
            id_columns    = list(sample_col = "SampleID", gene_id = "gene_id"),
            effects       = list(samples = "SampleID", color = "group"),
            filtering     = list(mode = "none", group_col = "group"),
            normalization = list(method = "VST")
        ))
    )
}

test_that("load_rna_inputs records which key filled each expression slot", {
    proj <- stage_two_inputs()
    cfg <- two_input_config(proj, list(counts = "rna/counts.csv", txi = "rna/txi.rds"))
    inputs <- suppressMessages(load_rna_inputs(cfg))

    expect_identical(inputs$input_keys$configured, c("counts", "txi"))
    expect_identical(inputs$input_keys$txi, "txi")
    expect_identical(inputs$input_keys$counts, "counts")
    expect_identical(inputs$source_type, "tximport")
})

test_that("with counts and txi both set, the log warns and the gene counts are the txi's", {
    skip_if_not_installed("DESeq2")
    proj <- stage_two_inputs(n_counts = 300L, n_txi = 100L)
    cfg <- two_input_config(proj, list(counts = "rna/counts.csv", txi = "rna/txi.rds"))
    inputs <- suppressMessages(load_rna_inputs(cfg))

    warnings <- character(0)
    pre <- withCallingHandlers(
        suppressMessages(preprocess_rna(inputs, cfg, verbose = FALSE)),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })

    input_warnings <- grep("More than one expression input", warnings, value = TRUE)
    expect_length(input_warnings, 1)
    expect_match(input_warnings, "This run used `files$txi`", fixed = TRUE)

    src <- pre$info$input_source
    expect_identical(src$source_type, "tximport")
    expect_identical(src$used_key, "txi")
    expect_true(src$conflict)

    # The counts are those of the tximport object the run used, not of the
    # 300-gene counts file that was also configured.
    expect_identical(rna_recorded_gene_counts(pre$info), list(before = 100L, after = 100L))
    expect_identical(nrow(pre$expr_raw), 100L)
})

test_that("a single counts input raises no input warning and is recorded", {
    skip_if_not_installed("edgeR")
    proj <- stage_two_inputs()
    cfg <- two_input_config(proj, list(counts = "rna/counts.csv"))
    cfg$modes$rna$normalization$method <- "TMMlogCPM"
    inputs <- suppressMessages(load_rna_inputs(cfg))

    warnings <- character(0)
    pre <- withCallingHandlers(
        suppressMessages(preprocess_rna(inputs, cfg, verbose = FALSE)),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })

    expect_false(any(grepl("More than one expression input", warnings)))
    expect_identical(pre$info$input_source$used_key, "counts")
    expect_false(pre$info$input_source$conflict)
    expect_identical(pre$info$n_genes_raw, 300L)
})

test_that("recording the input source leaves the analysis matrices unchanged", {
    skip_if_not_installed("edgeR")
    proj <- stage_two_inputs()
    cfg <- two_input_config(proj, list(counts = "rna/counts.csv"))
    cfg$modes$rna$normalization$method <- "TMMlogCPM"
    inputs <- suppressMessages(load_rna_inputs(cfg))
    # As the inputs were before load_rna_inputs() recorded the keys.
    legacy <- inputs
    legacy$input_keys <- NULL

    run <- function(x) suppressWarnings(suppressMessages(preprocess_rna(x, cfg, verbose = FALSE)))
    now <- run(inputs)
    before <- run(legacy)

    for (slot in c("expr_raw", "expr_filt", "expr_work", "de_input", "row_data", "meta")) {
        expect_identical(now[[slot]], before[[slot]], info = slot)
    }
    expect_identical(now$info[names(now$info) != "input_source"],
                     before$info[names(before$info) != "input_source"])
})

# ---- Report template --------------------------------------------------------

rna_report_src <- function() {
    rmd <- file.path(root_dir, "R", "domain", "rnaseq", "report_template.Rmd")
    skip_if_not(file.exists(rmd), "report_template.Rmd not found")
    readLines(rmd, warn = FALSE)
}

test_that("the report takes gene counts from the run, not from files$counts", {
    src <- rna_report_src()
    code <- src[!grepl("^\\s*#", src)]

    # The only remaining read of the raw input is the "before" histogram,
    # gated on the run having used files$counts as a raw count matrix.
    expect_false(any(grepl("rna_cfg$files$counts", code, fixed = TRUE)))
    expect_true(any(grepl("rna_gene_n <- rna_recorded_gene_counts(params$rna_info)", code,
                          fixed = TRUE)))
    expect_true(any(grepl('identical(rna_input_src$used_key, "counts")', code, fixed = TRUE)))
    expect_false(any(grepl("nrow(read.csv(", code, fixed = TRUE)))
})

test_that("the Input Data section shows the input used and the conflict note", {
    src <- rna_report_src()
    heading <- which(src == "## Input Data")
    expect_length(heading, 1)
    # Inside the Analysis Parameters tabset, after the Project tab.
    expect_gt(heading, which(src == "## Project"))
    expect_lt(heading, which(src == "## Experimental Design"))

    lines <- find_rmd_chunk_lines(file.path(root_dir, "R", "domain", "rnaseq",
                                            "report_template.Rmd"), "params-input")
    expect_false(is.null(lines))
    chunk <- src[seq(lines[1], lines[2])]
    expect_true(any(grepl("rna_input_src$description", chunk, fixed = TRUE)))
    expect_true(any(grepl("rna_input_src$note", chunk, fixed = TRUE)))
})
