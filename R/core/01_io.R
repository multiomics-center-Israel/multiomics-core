# =============================================================================
# Shared omics input loader
# =============================================================================

#' Normalize a contrast name for safe use in column names
#'
#' Strips spaces from contrast names so that column names like
#' "linearFC.imputs.<contrast>" are consistent between producers
#' (DE summary builders) and consumers (export, pathway, reports).
#'
#' @param x Character scalar: contrast name (e.g., "1.56ppm vs. 0ppm")
#' @return Character scalar with spaces removed (e.g., "1.56ppmvs.0ppm")
normalize_contrast_name <- function(x) {
    gsub(" ", "", x)
}

#' Load omics input files from config
#'
#' Generic loader for any omics mode. Validates required files, loads CSV/TSV
#' and RDS files, and validates contrasts content.
#'
#' @param config Configuration list
#' @param mode One of "proteomics" or "rna"
#' @return List of loaded objects
load_omics_inputs <- function(config, mode = c("proteomics", "rna")) {
    mode <- match.arg(mode)
    cfg <- config$modes[[mode]]
    if (is.null(cfg)) stop("No config for mode ", mode)

    files <- cfg$files

    # Determine required files based on mode and input format
    is_preprocessed <- identical(cfg$input$format, "preprocessed")
    required_files <- switch(mode,
        proteomics = if (is_preprocessed) c("preprocessed_protein", "metadata", "contrasts") else c("protein", "metadata", "contrasts"),
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
        # Skip NULL, empty, non-character (e.g. is_logtransformed: false),
        # or multi-value entries (e.g. de_table: [...]) — only load scalar paths
        if (is.null(rel) || !is.character(rel) || length(rel) != 1 || !nzchar(rel)) {
            next
        }
        abs <- resolve_raw_path(config, rel)
        if (!file.exists(abs)) stop("File not found: ", abs)

        # Detect file type and load appropriately
        ext <- tolower(tools::file_ext(abs))

        if (ext == "rds") {
            # RDS file - use readRDS (for tximport objects, etc.)
            message(sprintf("[load_omics_inputs] Loading RDS file: %s", nm))
            inputs[["txi"]] <- readRDS(abs)
        } else {
            inputs[[nm]] <- read_table_auto(abs)
        }
    }

    if (!is.null(cfg$engine)) inputs$engine <- cfg$engine

    # Validate contrasts file content (at least 1 row + expected columns)
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

# =============================================================================
# TSV/CSV utilities
# =============================================================================

#' Save a data frame as TSV (creating parent dir if needed)
#'
#' @param x Data frame.
#' @param dir Directory path.
#' @param filename Filename.
#' @return The full path.
save_tsv <- function(x, dir, filename) {
    ensure_dir(dir)
    path <- file.path(dir, filename)
    readr::write_tsv(x, path)
    path
}

#' Save a data frame as TSV to a full path
save_tsv_path <- function(x, path) {
    ensure_dir(dirname(path))
    readr::write_tsv(x, path)
    path
}

#' Read a table automatically detecting TSV vs CSV by extension
read_table_auto <- function(path) {
    ext <- tolower(tools::file_ext(path))
    df <- if (ext %in% c("tsv", "txt")) {
        readr::read_tsv(path, show_col_types = FALSE)
    } else {
        readr::read_csv(path, show_col_types = FALSE)
    }
    # Convert tibble to data.frame to support rownames and proper subsetting
    as.data.frame(df)
}
