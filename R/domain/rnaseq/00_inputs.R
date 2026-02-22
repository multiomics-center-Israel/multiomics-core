#' Load RNA inputs
#'
#' Supports both raw count matrices and tximport objects.
#' - For raw counts: specify files$counts pointing to CSV/TSV
#' - For tximport: specify files$txi pointing to RDS file containing tximport object
#'
#' @param config list as returned by load_config()
#' @return list with one of: (counts, ...) for raw counts OR (txi, ...) for tximport
load_rna_inputs <- function(config) {
    inputs <- load_omics_inputs(config, mode = "rna")

    # Check if txi was loaded
    if (!is.null(inputs$txi)) {
        # Validate tximport structure
        if (!is_valid_tximport_structure(inputs$txi, validate_only = TRUE)) {
            stop(
                "[load_rna_inputs] File loaded as 'txi' is not a valid tximport object. ",
                "Expected list with 'counts', 'abundance', 'length' matrices.",
                call. = FALSE
            )
        }
        message("[load_rna_inputs] Loaded tximport object from RDS")
        inputs$source_type <- "tximport"
    } else if (!is.null(inputs$counts)) {
        message("[load_rna_inputs] Loaded raw count matrix")
        inputs$source_type <- "matrix"
    }

    inputs
}

#' Validate RNA inputs (for raw count matrix path only)
#'
#' @param inputs List returned by \code{load_rna_inputs()}
#' @param cfg    The \code{config$modes$rna} sub-list
validate_rna_inputs <- function(inputs, cfg) {
    # This validation is only for raw count matrices
    # tximport validation is handled separately by validate_tximport()
    if (identical(inputs$source_type, "tximport") || !is.null(inputs$txi)) {
        # For tximport, metadata validation only
        sample_col <- cfg$id_columns$sample_col %||% "SampleID"
        meta <- inputs$metadata
        if (!sample_col %in% names(meta)) {
            stop("metadata missing sample column: ", sample_col)
        }
        return(invisible(TRUE))
    }

    gene_id_col <- cfg$id_columns$gene_id
    sample_col <- cfg$id_columns$sample_col %||% "SampleID"

    if (!gene_id_col %in% names(inputs$counts)) {
        stop("RNA counts missing gene id column: ", gene_id_col)
    }

    meta <- inputs$metadata
    if (!sample_col %in% names(meta)) {
        stop("metadata missing sample column: ", sample_col)
    }

    # Basic checks for counts (numeric columns only — exclude annotation columns)
    non_id_cols <- setdiff(names(inputs$counts), gene_id_col)
    is_numeric <- vapply(inputs$counts[, non_id_cols, drop = FALSE],
                         function(x) is.numeric(x) || is.integer(x),
                         logical(1))
    sample_cols <- non_id_cols[is_numeric]
    if (length(sample_cols) == 0) stop("Counts table has no numeric sample columns.")

    # Check alignment consistency (informational, strict check happens later)
    meta_samples <- unique(as.character(meta[[sample_col]]))
    count_samples <- sample_cols

    missing_in_counts <- setdiff(meta_samples, count_samples)
    if (length(missing_in_counts) > 0) {
        stop(
            "metadata contains samples not present in counts: ",
            paste(missing_in_counts, collapse = ", ")
        )
    }

    invisible(TRUE)
}

#' Generic loader helper
#'
#' Loads files specified in config. Supports:
#' - CSV/TSV files (loaded via read_table_auto)
#' - RDS files (loaded via readRDS, used for tximport objects)
#'
#' @param config Configuration list
#' @param mode One of "proteomics" or "rna"
#' @return List of loaded objects
load_omics_inputs <- function(config, mode = c("proteomics", "rna")) {
    mode <- match.arg(mode)
    cfg <- config$modes[[mode]]
    if (is.null(cfg)) stop("No config for mode ", mode)

    files <- cfg$files

    # Validate required files are specified in config
    # Note: RNA 'counts' is not strictly required here because tximport uses 'txi' key;
    # counts-vs-txi validation is handled separately in load_rna_inputs / preprocess_rna
    required_files <- switch(mode,
        proteomics = c("protein", "sample_map", "metadata", "contrasts"),
        rna = c("metadata", "contrasts"),
        character(0)
    )

    # Check 1: key completely missing from config
    missing_keys <- setdiff(required_files, names(files))
    if (length(missing_keys) > 0) {
        stop(
            sprintf(
                "[%s] Missing required file key(s) in config$modes$%s$files: %s",
                mode, mode, paste(missing_keys, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    # Check 2: key present but set to null or empty string
    for (nm in required_files) {
        if (is.null(files[[nm]]) || !nzchar(files[[nm]])) {
            stop(
                sprintf(
                    "[%s] File '%s' is required but has null or empty path in config$modes$%s$files",
                    mode, nm, mode
                ),
                call. = FALSE
            )
        }
    }

    inputs <- list()

    for (nm in names(files)) {
        rel <- files[[nm]]
        # Skip NULL or empty file paths (only for optional files — required already validated above)
        if (is.null(rel) || !nzchar(rel)) {
            next
        }
        abs <- resolve_raw_path(config, rel)
        # Check 3: path present but file doesn't exist
        if (!file.exists(abs)) stop("File not found: ", abs)

        # Detect file type and load appropriately
        ext <- tolower(tools::file_ext(abs))

        if (ext == "rds") {
            # RDS file - use readRDS (for tximport objects, etc.)
            message(sprintf("[load_omics_inputs] Loading RDS file: %s", nm))
            nm <- "txi"
            inputs[[nm]] <- readRDS(abs)
        } else {
            # Tabular file - use standard reader
            inputs[[nm]] <- read_table_auto(abs)
        }
    }

    if (!is.null(cfg$engine)) inputs$engine <- cfg$engine

    # Check 4: validate contrasts file content (at least 1 row + expected columns)
    if ("contrasts" %in% required_files && !is.null(inputs$contrasts)) {
        validate_contrasts_content(inputs$contrasts, mode)
    }

    inputs
}

#' Validate contrasts file content
#'
#' Ensures the loaded contrasts data frame has at least one row and contains
#' all required columns: Contrast_name, Factor, Numerator, Denominator.
#'
#' @param contrasts_df Data frame loaded from the contrasts file
#' @param mode Character string identifying the omics mode (for error messages)
#' @return invisible(TRUE) on success, stops with error on failure
validate_contrasts_content <- function(contrasts_df, mode = "omics") {
    if (!is.data.frame(contrasts_df)) {
        stop(
            sprintf("[%s] Contrasts file did not load as a data frame.", mode),
            call. = FALSE
        )
    }

    if (nrow(contrasts_df) == 0) {
        stop(
            sprintf(
                "[%s] Contrasts file is empty (0 rows). At least one contrast is required.",
                mode
            ),
            call. = FALSE
        )
    }

    required_cols <- c("Contrast_name", "Factor", "Numerator", "Denominator")
    missing_cols <- setdiff(required_cols, colnames(contrasts_df))
    if (length(missing_cols) > 0) {
        stop(
            sprintf(
                "[%s] Contrasts file is missing required column(s): %s. Expected: %s",
                mode,
                paste(missing_cols, collapse = ", "),
                paste(required_cols, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}
