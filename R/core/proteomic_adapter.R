# R/core/proteomics_adapter.R

as_legacy_pre <- function(pre) {
  pre$expr_filt <- pre$expr_filt_mat %||% pre$expr_filt
  pre$expr_imp  <- pre$expr_imp_single %||% pre$expr_imp
  
  if (is.null(pre$imputation)) pre$imputation <- list()
  if (is.null(pre$imputation$imputed_flag) && !is.null(pre$imputation_qc$imputed_flag)) {
    pre$imputation$imputed_flag <- pre$imputation_qc$imputed_flag
  }
  
  pre
}

#' Build a standardized proteomics expression object (engine-aware)
#'
#' For now this supports DIANN only, but keeps a stable contract so that
#' adding MaxQuant later does not require changing downstream code.
#'
#' Contract:
#'   - assay_log2   : matrix (features × samples), log2 scale (downstream contract)
#'   - assay_linear : optional linear-scale matrix (NULL unless available/needed)
#'   - row_data     : feature annotations aligned to assay rows
#'   - col_data     : sample metadata aligned to assay columns
#'   - info         : engine/scale metadata
#'
#' @param inputs List from load_proteomics_inputs().
#' @param config Full config list.
#'
#' @return List with fields: assay_log2, assay_linear, row_data, col_data, info.
get_proteomics_expression_matrix <- function(inputs, config) {
  
  cfg <- config$modes$proteomics
  
  # Build matrix depending on engine
  if (cfg$engine == "DIANN") {
    
    # DIA-NN → expression matrix (your current builder already returns log2)
    assay_log2 <- get_measurements_per_sample_diann(
      protein    = inputs$protein,
      sample_map = inputs$sample_map,
      meta       = inputs$metadata,
      cfg        = cfg
    )
    
    assay_linear <- NULL
    
    # Feature annotations (row_data)
    row_data <- inputs$protein[ , c(cfg$id_columns$protein_id, unlist(cfg$id_columns$protein_annot)), drop = FALSE]
    
  } else {
    stop(sprintf("Unsupported proteomics engine: %s", cfg$engine))
  }
  
  # Align col_data to assay columns 
  col_data <- inputs$metadata
  # sample_col <- cfg$effects$samples %||% cfg$id_columns$sample_col
  col_data <- align_meta_to_expr(assay_log2, inputs$metadata, cfg)

  list(
    assay_log2   = assay_log2,
    assay_linear = assay_linear,
    row_data     = row_data,
    col_data     = col_data,
    info = list(
      mode         = "proteomics",
      engine       = cfg$engine,
      scale_in     = cfg$scale_in %||% "log2",
      target_scale = cfg$transform$target_scale %||% "log2"
    )
  )
}

