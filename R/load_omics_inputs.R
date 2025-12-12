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
    stop("No configuration found for mode '", mode, "' in config$modes$", mode)
  }
  if (is.null(cfg$files) || !length(cfg$files)) {
    stop("config$modes$", mode, "$files is missing or empty.")
  }
  
  # ---- check files exist ----
  files <- cfg$files  # named list: protein / sample_map / metadata / contrasts / ...
  for (nm in names(files)) {
    f <- files[[nm]]
    if (is.null(f) || !nzchar(f)) {
      stop("In mode '", mode, "': path for '", nm, "' is missing/empty in cfg$files.")
    }
    if (!file.exists(f)) {
      stop("In mode '", mode, "': file not found for '", nm, "': ", f)
    }
  }
  
  # ---- read all CSVs in a generic way ----
  inputs <- lapply(
    files,
    function(path) readr::read_csv(path, show_col_types = FALSE)
  )
  names(inputs) <- names(files)
  
  # ---- attach engine if present ----
  if (!is.null(cfg$engine)) {
    inputs$engine <- cfg$engine
  }
  
  inputs
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

