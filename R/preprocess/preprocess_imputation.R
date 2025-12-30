#' Perseus-style imputation (downshifted & narrowed normal distribution)
#'
#' For each sample (column), missing values are imputed from a normal
#' distribution with:
#'   mean = mean(sample, na.rm = TRUE) - downshift * sd(sample, na.rm = TRUE)
#'   sd   = width * sd(sample, na.rm = TRUE)
#'
#' @param expr_mat  Numeric matrix/data.frame (features x samples), usually log2.
#' @param width  Fraction of the sample sd used as sd of the imputation dist.
#' @param downshift How many sd to subtract from the mean for the imputation mean.
#'
#' @return A list with:
#'   - imputed: numeric matrix with imputed values
#'   - imputed_flag: logical matrix, TRUE where a value was imputed
perseus_impute_with_flags <- function(expr_mat, width = 0.3, downshift = 1.8) {
  # Coerce to matrix for safety
  expr_mat <- as.matrix(expr_mat)
  sample_cols <- colnames(expr_mat)
  
  # Logical matrix of missing values (these will be imputed)
  imputed_flag <- is.na(expr_mat)
  
  # Copy for imputation
  imputed <- expr_mat
  
  for (j in seq_len(ncol(imputed))) {
    x <- imputed[, j]
    
    # stop if all values are NA
    if (all(is.na(x))) {
      stop("Imputation failed: sample '", colnames(imputed)[j], "' is all-NA after filtering.")
    }

    
    # Compute sd and mean of observed values
    obs <- x[!is.na(x)]
    s <- stats::sd(obs)
    if (!is.finite(s) || s == 0) s <- 1e-8
    m <- mean(obs)
    
    imp_sd   <- width * s
    imp_mean <- m - downshift * s
    
    n_missing <- sum(is.na(x))
    if (n_missing > 0) {
      x[is.na(x)] <- stats::rnorm(n_missing, mean = imp_mean, sd = imp_sd)
    }
    
    imputed[, j] <- x
  }
  
  colnames(imputed)      <- sample_cols
  rownames(imputed_flag) <- rownames(imputed)
  colnames(imputed_flag) <- sample_cols
  
  list(
    imputed      = imputed,
    imputed_flag = imputed_flag
  )
}

make_imputations_proteomics <- function(expr_mat,
                                        cfg,
                                        seed_base = NULL,
                                        verbose = FALSE) {
  
  imp_cfg <- cfg$modes$proteomics$imputation
  n_imputations <- as.integer(imp_cfg$no_repetitions)
  
  stopifnot(is.matrix(expr_mat))
  stopifnot(!is.na(n_imputations), n_imputations >= 1)
  
  imps <- vector("list", n_imputations)
  
  for (i in seq_len(n_imputations)) {
    if (isTRUE(verbose)) {
      message(sprintf("Imputation: %d / %d", i, n_imputations))
    }
    
    if (!is.null(seed_base)) {
      set.seed(as.integer(seed_base) + i)
    }
    
    expr_imp_i <- impute_proteomics_perseus(
      expr_mat,
      cfg = cfg$modes$proteomics,
      return_flags = FALSE
    )
    
    stopifnot(is.matrix(expr_imp_i))
    stopifnot(identical(dim(expr_imp_i), dim(expr_mat)))
    stopifnot(identical(colnames(expr_imp_i), colnames(expr_mat)))
    stopifnot(identical(rownames(expr_imp_i), rownames(expr_mat)))
    
    imps[[i]] <- expr_imp_i
  }
  
  imps
}

