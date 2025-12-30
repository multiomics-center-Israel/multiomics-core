
# simple helper for "x %||% default"
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}