#' Determine which features pass the filter based on min_count per condition
#'
#' Returns a logical vector where TRUE means the feature passed the filter.
#'
#' @param expr_mat Numeric matrix/data.frame (features x samples)
#' @param group Condition/factor for each sample
#' @param min_per_group Named numeric vector of min_count thresholds
#'
#' @return Logical vector (TRUE/FALSE per feature)
pass_filter <- function(expr_mat, group, min_per_group) {
  expr_mat <- as.matrix(expr_mat)
  
  if (length(group) != ncol(expr_mat)) {
    stop("Length of `group` must equal number of columns in `expr_mat`.")
  }
  
  group  <- as.character(group)
  groups <- unique(group)
  
  # Allow single threshold without names
  if (length(min_per_group) == 1 && is.null(names(min_per_group))) {
    min_per_group <- setNames(rep(min_per_group, length(groups)), groups)
  }
  
  if (!all(groups %in% names(min_per_group))) {
    stop(
      "min_per_group must have names for all groups. Missing: ",
      paste(setdiff(groups, names(min_per_group)), collapse = ", ")
    )
  }
  
  passes_per_group <- sapply(groups, function(g) {
    cols <- which(group == g)
    if (length(cols) == 0) return(rep(FALSE, nrow(expr_mat)))
    
    sums <- rowSums(!is.na(expr_mat[, cols, drop = FALSE]))
    sums >= min_per_group[[g]]
  })
  
  apply(passes_per_group, 1, any)
}

#' Extract min_count per condition from config
#'
#' Supports all of these config shapes:
#' - min_count: 3
#' - min_count: { default: 3 }
#' - min_count: { default: 3, C: 2 }
#' - min_count: { C: 3, S: 3, SH: 2 }
#'
#' @param min_cfg config$modes$proteomics$filtering$min_count
#' @param groups  character vector of condition values (one per sample)
#'
#' @return Named numeric vector of length unique(groups)
extract_min_count <- function(min_cfg, groups) {
  groups <- unique(as.character(groups))
  
  # Case 0: numeric scalar (e.g. min_count: 3)
  if (is.numeric(min_cfg) && length(min_cfg) == 1 && is.null(names(min_cfg))) {
    return(setNames(rep(min_cfg, length(groups)), groups))
  }
  
  # Case 1: list-like with $default and optional overrides
  if (!is.null(min_cfg$default)) {
    out <- setNames(rep(min_cfg$default, length(groups)), groups)
    
    overrides <- min_cfg[setdiff(names(min_cfg), "default")]
    if (length(overrides) > 0) {
      for (nm in names(overrides)) {
        if (nm %in% names(out)) {
          out[nm] <- overrides[[nm]]
        }
      }
    }
    return(out)
  }
  
  # Case 2: fully named thresholds (no default)
  out <- unlist(min_cfg)
  names(out) <- names(min_cfg)
  out
}