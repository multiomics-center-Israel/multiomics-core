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
read_table_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  read_fn <- if (ext %in% c("tsv", "txt")) readr::read_tsv else readr::read_csv
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


# =============================================================================
# Pre-computed DE summary tables
# =============================================================================

#' Resolve a column in a per-contrast DE summary table
#'
#' Our own \code{Datasets/*_summary_p0.05.tsv} exports hold every contrast in one
#' table and suffix the statistic columns with the contrast name, e.g.
#' \code{log2FC.S_vs_NS} or \code{padj.imputs.SP_vs_NSP}. The pre-computed DE
#' loaders originally matched bare names only, so pointing them at one of our own
#' exports silently produced all-NA statistics. This resolves both shapes.
#'
#' Matching order: an exact bare name first, then a suffixed column. When several
#' contrasts are present the label decides; when only one is present it is used
#' regardless of what it is called, since the multiomics contrast label often
#' differs from the short code the single-omics run used.
#'
#' @param cn Character vector of column names in the table.
#' @param bare Candidate bare column names, in preference order.
#' @param prefixes Candidate prefixes for the suffixed form (e.g. "log2FC",
#'   "linearFC.imputs"), in preference order.
#' @param contrast_label Contrast name to prefer when the table holds several.
#' @return The resolved column name, or NA_character_ when nothing matches.
resolve_de_summary_col <- function(cn, bare, prefixes, contrast_label = NULL) {
    hit <- cn[cn %in% bare]
    if (length(hit) > 0) return(hit[1])

    for (pre in prefixes) {
        stem <- paste0(pre, ".")
        suffixed <- cn[startsWith(cn, stem)]
        if (length(suffixed) == 0) next
        if (!is.null(contrast_label)) {
            exact <- suffixed[suffixed == paste0(stem, contrast_label)]
            if (length(exact) > 0) return(exact[1])
        }
        # One contrast in the file: its short code need not match the label.
        if (length(suffixed) == 1) return(suffixed[1])
        if (!is.null(contrast_label)) {
            stop("Cannot resolve column '", pre, "' for contrast '", contrast_label,
                 "': the table holds several contrasts (",
                 paste(substring(suffixed, nchar(stem) + 1L), collapse = ", "),
                 "). Rename the contrast or split the table.")
        }
    }
    NA_character_
}


#' Convert a signed linear fold change to log2
#'
#' Proteomics summaries store linearFC as a signed linear ratio: 2^log2FC when
#' the change is non-negative and -1 / 2^log2FC when it is negative (see
#' build_provenance_notes() in R/core/05_export_excel.R). A plain log2() returns
#' NaN for every down-regulated feature, silently dropping about half the
#' proteome wherever the value is reused.
#'
#' @param x Numeric vector of signed linear fold changes.
#' @return Numeric vector of log2 fold changes; NA where x is NA or zero.
signed_linear_fc_to_log2 <- function(x) {
    x <- as.numeric(x)
    out <- rep(NA_real_, length(x))
    # Index rather than ifelse(): ifelse evaluates both branches, so log2() of
    # the negative values emits a NaN warning even though those are discarded.
    up   <- !is.na(x) & x > 0
    down <- !is.na(x) & x < 0
    out[up]   <- log2(x[up])
    out[down] <- -log2(abs(x[down]))
    out
}
