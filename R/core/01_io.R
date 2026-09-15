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

#' Read a sample sheet, guarding both ways read.csv() mangles one
#'
#' \code{read.csv()} fails on a sample sheet in two independent ways, and a
#' reader that handles one still falls to the other:
#'
#' \enumerate{
#'   \item \strong{Wrong delimiter.} Sample sheets are frequently tab-separated
#'     (.txt/.tsv). Read as CSV, every column collapses into one, and a consumer
#'     then finds none of the columns it expects.
#'   \item \strong{Ragged rows.} \code{read.csv()} treats column 1 as row names
#'     whenever a data row has MORE fields than the header. One unquoted comma in
#'     a free-text column is enough. The failure is silent: the frame keeps the
#'     right column NAMES while every value sits one column to the left. On a
#'     real run that put the raw file path into \code{SampleName}, no expression
#'     column matched a sample id, and the report's explorer rendered empty with
#'     no error. \code{row.names = NULL} does not fix it -- only a reader that
#'     respects the header's column count does.
#' }
#'
#' \code{read_table_auto()} already reads with readr, so it does not shift; what
#' it cannot do is spot a sheet whose extension lies about its delimiter, since
#' it decides from the extension alone. This wrapper reads the separator off the
#' header line and hands it down, and returns NULL instead of erroring so a
#' report chunk can carry on without the sheet. Everything else -- the Latin1
#' retry, the character sanitization, the data.frame conversion -- is
#' \code{read_table_auto()}'s and is not duplicated here.
#'
#' @param path Path to the sample sheet (CSV or TSV).
#' @return A data.frame, or \code{NULL} when \code{path} is missing, empty or
#'   unreadable.
read_samplesheet <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)

  l1  <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) character(0))
  sep <- if (length(l1) > 0 && grepl("\t", l1, fixed = TRUE)) "\t" else ","

  tryCatch(read_table_auto(path, sep = sep), error = function(e) NULL)
}

#' Load omics input files from config
#'
#' Generic loader for any omics mode. Validates required files, loads CSV/TSV
#' and RDS files, and validates contrasts content.
#'
#' @param config Configuration list
#' @param mode One of "proteomics", "rna", or "metabolomics"
#' @return List of loaded objects
load_omics_inputs <- function(config, mode = c("proteomics", "rna", "metabolomics")) {
  mode <- match.arg(mode)
  cfg <- config$modes[[mode]]
  if (is.null(cfg)) stop("No config for mode ", mode)
  
  files <- cfg$files
  
  # Determine required files based on mode and input format
  is_preprocessed <- identical(cfg$input$format, "preprocessed")
  required_files <- switch(mode,
                           proteomics = if (is_preprocessed) c("preprocessed_protein", "metadata", "contrasts") else c("protein", "metadata", "contrasts"),
                           rna = c("metadata", "contrasts"),
                           metabolomics = c("metadata", "contrasts"),
                           character(0)
  )

  # The proteomics "limma_percontrast" method auto-generates control-referenced
  # contrasts from de$control_condition, so a contrasts file is not required in
  # that configuration. Any explicitly supplied file is still loaded and validated
  # below; every other mode/method continues to require contrasts.
  if (mode == "proteomics" &&
      identical(cfg$de$method, "limma_percontrast") &&
      !is.null(cfg$de$control_condition) && nzchar(cfg$de$control_condition)) {
    required_files <- setdiff(required_files, "contrasts")
  }
  
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
    if (dir.exists(abs)) next    # skip directory paths (e.g., data_dir)
    if (!file.exists(abs)) stop(sprintf("[%s] File '%s' not found at: %s", mode, nm, abs), call. = FALSE)
    
    # Detect file type and load appropriately
    ext <- tolower(tools::file_ext(abs))
    
    if (ext == "rds") {
      # RDS file - use readRDS (for tximport objects, etc.)
      message(sprintf("[load_omics_inputs] Loading RDS file: %s", nm))
      inputs[["txi"]] <- readRDS(abs)
    } else {
      inputs[[nm]] <- read_metab_file(abs)
    }
  }
  
  if (!is.null(cfg$engine)) inputs$engine <- cfg$engine
  
  # Validate contrasts file content (at least 1 row + expected columns) whenever a
  # contrasts file was loaded - including the optional case where limma_percontrast
  # could auto-generate them but the user still supplied a file.
  if (!is.null(inputs$contrasts)) {
    validate_contrasts_content(inputs$contrasts, mode)
  }

  # Warn early about commas in a tab-separated sample sheet — they parse fine
  # here but break the CSV-assuming report readers at render time.
  if (!is.null(inputs$metadata) && is.character(files$metadata) && nzchar(files$metadata)) {
    check_metadata_delimiter_safety(
      as.data.frame(inputs$metadata),
      resolve_raw_path(config, files$metadata),
      mode
    )
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

#' Sanitize character columns: convert to UTF-8, replace NBSP, trim whitespace.
#' Emits a warning listing affected columns and modification counts.
sanitize_character_columns <- function(df, source = "input") {
  chr_cols <- which(vapply(df, is.character, logical(1)))
  if (length(chr_cols) == 0L) return(df)
  
  modified_summary <- character(0)
  
  for (j in chr_cols) {
    orig <- df[[j]]
    cleaned <- iconv(orig, from = "", to = "UTF-8", sub = "")
    cleaned <- gsub("\u00A0", " ", cleaned, fixed = TRUE)
    cleaned <- trimws(cleaned)
    
    n_changed <- sum(orig != cleaned, na.rm = TRUE)
    if (n_changed > 0L) {
      modified_summary <- c(modified_summary,
                            sprintf("  - %s: %d value(s)", names(df)[j], n_changed))
    }
    df[[j]] <- cleaned
  }
  
  if (length(modified_summary) > 0L) {
    warning(
      sprintf("Encoding/whitespace issues sanitized in %s:\n", source),
      paste(modified_summary, collapse = "\n"),
      call. = FALSE
    )
  }
  
  df
}

#' Read a table automatically detecting TSV vs CSV by extension
#'
#' @param path Path to the file.
#' @param sep Optional separator, \code{"\t"} or \code{","}. Overrides the
#'   extension, for callers that have determined the delimiter another way (see
#'   \code{\link{read_samplesheet}}, which reads it off the header because a
#'   sample sheet's extension often lies). \code{NULL} keeps the extension rule,
#'   so existing callers are unaffected.
#' @return A data.frame.
read_table_auto <- function(path, sep = NULL) {
  ext <- tolower(tools::file_ext(path))
  use_tsv <- if (is.null(sep)) ext %in% c("tsv", "txt") else identical(sep, "\t")
  read_fn <- if (use_tsv) readr::read_tsv else readr::read_csv
  df <- tryCatch(
    read_fn(path, show_col_types = FALSE),
    error = function(e) {
      if (grepl("invalid.*UTF-8|invalid.*utf8", conditionMessage(e),
                ignore.case = TRUE)) {
        warning(
          sprintf("Non-UTF-8 encoding detected in %s; re-reading as Latin1.",
                  basename(path)),
          call. = FALSE
        )
        read_fn(path, show_col_types = FALSE,
                locale = readr::locale(encoding = "Latin1"))
      } else {
        stop(e)
      }
    }
  )
  # Convert tibble to data.frame to support rownames and proper subsetting
  df <- as.data.frame(df)
  # Sanitize character columns: NBSP, whitespace
  df <- sanitize_character_columns(df, source = basename(path))
  df
}
