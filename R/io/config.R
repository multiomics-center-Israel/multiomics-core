# R/config.R

#' Load pipeline configuration from YAML.
#'
#' @param path Path to a YAML config file (default: "config/config.yml").
#' @return A named list with the parsed configuration.
load_config <- function(path = "config/config.yml") {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Install it with install.packages('yaml').")
  }
  
  cfg <- yaml::read_yaml(path)
  
  cfg$.config_path  <- normalizePath(path, mustWork = FALSE)
  cfg$.config_mtime <- file.info(path)$mtime
  
  return(cfg)
}


# Get effects (color / shape / samples) for a given modality (e.g. "rna", "proteomics")
# R/config.R
get_effects <- function(cfg, mode) {
  mode_cfg <- cfg$modes[[mode]]
  if (is.null(mode_cfg)) {
    stop(sprintf("Mode '%s' not found under config$modes.", mode))
  }
  
  eff <- mode_cfg$effects
  if (is.null(eff)) eff <- list()
  
  list(
    color   = eff$color   %||% NULL,
    shape   = eff$shape   %||% NULL,
    samples = eff$samples %||% "SampleID"
  )
}


