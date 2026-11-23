# Which expression input an RNA run used
#
# load_rna_inputs() records which config key filled each expression slot,
# preprocess_rna() records which slot it used, and the log warning and the
# report's Input Data section and gene counts all read that record through
# the helpers here.

#' Expression-input keys set in the RNA config
#'
#' `counts`, `txi` and `preprocessed_counts` can each supply the expression
#' matrix. Setting more than one is not an error -- the run builds on one of
#' them -- but the others are still read from disk, and a report that looked
#' at the wrong one would describe data the run never used.
#'
#' @param files The `config$modes$rna$files` list.
#' @return The keys among `counts`, `txi` and `preprocessed_counts` that hold
#'   a path, in that order. The test is the one `load_omics_inputs()` applies
#'   before reading a file.
#' @examples
#' rna_configured_expression_inputs(list(counts = "rna/c.csv", txi = "rna/txi.rds"))
#' rna_configured_expression_inputs(list(counts = "", txi = "rna/txi.rds"))  # "txi"
rna_configured_expression_inputs <- function(files) {
    keys <- c("counts", "txi", "preprocessed_counts")
    keys[vapply(keys, function(k) rna_is_path_value(files[[k]]), logical(1))]
}

#' Config key whose file `load_omics_inputs()` loaded as the tximport object
#'
#' `load_omics_inputs()` stores every `.rds` file under `inputs$txi`, whatever
#' its key, so a later one replaces an earlier one. This applies that same
#' rule to the keys, so the key reported is the one whose object was kept.
#'
#' @param files The `config$modes$rna$files` list.
#' @return The key, or `NA_character_` when no `.rds` path is configured.
#' @examples
#' rna_rds_input_key(list(counts = "rna/c.csv", txi = "rna/txi.rds"))  # "txi"
rna_rds_input_key <- function(files) {
    is_rds <- vapply(names(files), function(k) {
        rna_is_path_value(files[[k]]) && tolower(tools::file_ext(files[[k]])) == "rds"
    }, logical(1))
    keys <- names(files)[is_rds]
    if (length(keys) == 0) NA_character_ else keys[length(keys)]
}

#' Does a config value name a file path?
#'
#' @param value One entry of `config$modes$rna$files`.
#' @return `TRUE` for a single non-empty string, the only values
#'   `load_omics_inputs()` reads.
rna_is_path_value <- function(value) {
    is.character(value) && length(value) == 1 && !is.na(value) && nzchar(value)
}

#' Describe the expression input a run was built on
#'
#' The one reading of "which input did this run use", shared by the log
#' warning in [preprocess_rna()] and the report's Input Data section. The
#' input is the slot [preprocess_rna()] took (`source_type`), traced back to
#' the config key [load_rna_inputs()] filled it from -- not a guess from the
#' keys that happen to be set.
#'
#' @param input_keys `inputs$input_keys` from [load_rna_inputs()]; `NULL` for
#'   inputs assembled without it, in which case the key is not known.
#' @param source_type The input type [preprocess_rna()] used: `"tximport"`,
#'   `"matrix"` or `"preprocessed"`.
#' @param files The `config$modes$rna$files` list.
#' @return A list with `source_type`, `used_key` and `used_file` (`NA` when
#'   unknown), `configured` (the keys set), `conflict` (`TRUE` when more than
#'   one is set), `description` (one phrase naming the input used) and `note`
#'   (the warning text when `conflict`, otherwise `NULL`).
#' @examples
#' describe_rna_input_source(
#'     list(configured = c("counts", "txi"), txi = "txi", counts = "counts"),
#'     "tximport", list(counts = "rna/c.csv", txi = "rna/txi.rds"))$note
describe_rna_input_source <- function(input_keys, source_type, files) {
    configured <- as.character(input_keys$configured %||% character(0))
    used_key <- if (identical(source_type, "tximport")) input_keys$txi else input_keys$counts
    if (!rna_is_path_value(used_key)) used_key <- NA_character_
    used_file <- if (is.na(used_key) || !rna_is_path_value(files[[used_key]])) {
        NA_character_
    } else {
        files[[used_key]]
    }

    read_as <- switch(source_type,
        tximport     = "a tximport object",
        matrix       = "a raw count matrix",
        preprocessed = "a preprocessed expression matrix",
        sprintf("input of type '%s'", source_type))
    description <- if (is.na(used_key)) {
        sprintf("%s (config key not recorded)", read_as)
    } else if (is.na(used_file)) {
        sprintf("`files$%s`, read as %s", used_key, read_as)
    } else {
        sprintf("`files$%s` (`%s`), read as %s", used_key, used_file, read_as)
    }

    conflict <- length(configured) > 1
    note <- if (conflict) {
        set <- vapply(configured, function(k) sprintf("`files$%s` (`%s`)", k, files[[k]]),
                      character(1))
        sprintf(paste0(
            "More than one expression input is set in the config: %s. This run used %s, ",
            "and every gene count reported for it refers to that input only. ",
            "Remove the unused key from modes.rna.files to make the choice explicit."),
            paste(set, collapse = " and "), description)
    } else {
        NULL
    }

    list(source_type = source_type, used_key = used_key, used_file = used_file,
         configured = configured, conflict = conflict, description = description,
         note = note)
}

#' Gene counts before and after the gene filter, as the run recorded them
#'
#' Read from [preprocess_rna()]'s `info`, i.e. the matrix the run actually
#' filtered, never from the input files on disk: with more than one input
#' configured, the file named by `files$counts` need not be the one used.
#'
#' @param info `pre$info` from [preprocess_rna()] (the report's
#'   `params$rna_info`), or `NULL`.
#' @return A list with integer `before` and `after`; `NA` when not recorded.
#' @examples
#' rna_recorded_gene_counts(list(n_genes_raw = 120L, n_genes_filt = 80L))
#' rna_recorded_gene_counts(NULL)  # both NA
rna_recorded_gene_counts <- function(info) {
    count <- function(x) {
        if (is.numeric(x) && length(x) == 1 && !is.na(x)) as.integer(x) else NA_integer_
    }
    list(before = count(info$n_genes_raw), after = count(info$n_genes_filt))
}
