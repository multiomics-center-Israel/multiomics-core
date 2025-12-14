#' Extract a matrix from an imputation result
#'
#' This helper standardizes access to imputed expression matrices.
#' It supports both:
#'   - a matrix returned directly, or
#'   - a list containing the matrix under a known key.
#'
#' This avoids repetitive if/else logic in downstream code.
#'
#' @param x Imputation result (matrix or list).
#' @param keys Character vector of candidate list element names
#'        that may contain the matrix.
#'
#' @return A numeric matrix.
#'
extract_imputed_matrix <- function(x,
                                   keys = c("expr_imp", "expr", "imputed")) {
  if (is.matrix(x)) return(x)
  
  stopifnot(is.list(x))
  found <- intersect(keys, names(x))
  stopifnot(length(found) >= 1)
  
  mat <- x[[found[1]]]
  stopifnot(is.matrix(mat))
  
  mat
}
