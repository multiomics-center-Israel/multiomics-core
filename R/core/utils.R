
# simple helper for "x %||% default"
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

get_sample_filter_rules <- function(config, mode) {
  cfg <- config$modes[[mode]]
  if (is.null(cfg)) return(NULL)
  sf <- cfg$sample_filter
  if (is.null(sf) || !isTRUE(sf$enabled)) return(NULL)
  sf$rules %||% NULL
}
