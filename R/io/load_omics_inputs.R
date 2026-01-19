#' Load omics inputs (generic)
#'
#' Generic loader for any omics mode (proteomics / rna / etc.).
#' Expects:
#'   config$modes[[mode]]$files – a named list of paths
#'   config$modes[[mode]]$engine (optional)
#'
#' @param config list as returned by load_config()
#' @param mode character, e.g. "proteomics", "rna"
#' @return named list of tibbles + optional engine element
load_omics_inputs <- function(config, mode = c("proteomics", "rna")) {
  mode <- match.arg(mode)
  
  # ---- get mode config ----
  cfg <- config$modes[[mode]]
  if (is.null(cfg)) {
    stop("No configuration found for mode '", mode, "' in config$modes$")
  }
  if (is.null(cfg$files) || !length(cfg$files)) {
    stop("config$modes$", mode, "$files is missing or empty.")
  }
  
  files <- cfg$files
  
  inputs <- list()
  
  for (nm in names(files)) {
    f <- files[[nm]]
    
    # ---- OPTIONAL file (NULL / NA / empty) ----
    if (is.null(f) || is.na(f) || !nzchar(f)) {
      inputs[[nm]] <- NULL
      message(
        sprintf(
          "[%s] Optional input '%s' not provided (NULL) – skipping.",
          mode, nm
        )
      )
      next
    }
    
    # ---- path provided but file missing ----
    if (!file.exists(f)) {
      stop(
        "In mode '", mode, "': file not found for '", nm, "': ", f
      )
    }
    
    # ---- read file ----
    inputs[[nm]] <- read_table_auto(f)
  }
  
  # ---- attach engine if present ----
  if (!is.null(cfg$engine)) {
    inputs$engine <- cfg$engine
  }
  
  inputs
}

read_table_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else {
    readr::read_csv(path, show_col_types = FALSE)
  }
}


#' Load proteomics inputs 
#'
#' @param config list as returned by load_config()
#' @return list (protein, sample_map, meta, contrasts, engine, ...)
load_proteomics_inputs <- function(config) {
  load_omics_inputs(config, mode = "proteomics")
}

#' Load rna inputs
#'
#' @param config list as returned by load_config()
#' @return list (genes, sample_map, meta, contrasts, engine, ...)
load_rna_inputs <- function(config) {
  load_omics_inputs(config, mode = "rna")
}

