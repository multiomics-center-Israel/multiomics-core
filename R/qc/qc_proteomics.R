#------- HELPERS ------- 

#' Prepare data for QC plotting
#'
#' Handles common tasks: matrix conversion, metadata matching, and annotation creation.
#' Includes strict validation to fail fast on data inconsistencies.
#'
#' @param expr Numeric matrix or data.frame.
#' @param meta Metadata data.frame.
#' @param cfg Config object containing effects$samples and effects$color.
#' @return A list containing aligned expr, meta, annot, and column names.
#' Prepare data for QC plotting
#'
#' Handles common tasks: matrix conversion, metadata matching, and annotation creation.
#' Includes strict validation to fail fast on data inconsistencies (duplicates, missing samples).
#'
#' @param expr Numeric matrix or data.frame.
#' @param meta Metadata data.frame.
#' @param cfg Config object containing effects$samples and effects$color.
#' @return A list containing aligned expr, meta, annot, and column names.
prepare_qc_data <- function(expr, meta, cfg) {
  
  # 1. Ensure Expression is a Matrix with Names
  expr <- as.matrix(expr)
  sample_ids <- colnames(expr)
  
  # FIX: Handle empty matrix or NULL colnames
  if (is.null(sample_ids) || length(sample_ids) == 0) {
    stop("Expression matrix must have non-empty column names (sample IDs).")
  }
  
  if (anyDuplicated(sample_ids)) {
    stop("Expression matrix contains duplicate column names/sample IDs.")
  }
  
  # 2. Extract Config & Validate Keys
  eff <- cfg$effects
  if (is.null(eff$samples) || is.null(eff$color)) {
    stop("Config must contain cfg$effects$samples and cfg$effects$color")
  }
  sample_col <- eff$samples
  color_col  <- eff$color
  
  # 3. Ensure Metadata is a base data.frame (safe against tibbles)
  meta <- as.data.frame(meta)
  
  # 4. Validate Columns Existence in Metadata
  if (!sample_col %in% colnames(meta)) {
    stop(sprintf("Sample column '%s' not found in metadata.", sample_col))
  }
  if (!color_col %in% colnames(meta)) {
    stop(sprintf("Color/Condition column '%s' not found in metadata.", color_col))
  }
  
  # FIX: Critical check for duplicates in metadata ID column
  # match() returns the first hit, so duplicates would be silently ignored without this.
  if (anyDuplicated(meta[[sample_col]])) {
    stop(sprintf("Metadata column '%s' contains duplicate sample IDs.", sample_col))
  }
  
  # 5. Match Metadata to Expression Columns
  # This guarantees that meta_sub rows are 1:1 with expr columns
  meta_sub <- meta[match(sample_ids, meta[[sample_col]]), , drop = FALSE]
  
  # 6. Verify Integrity (Missing samples)
  # FIX: Detailed error message with count + examples
  if (any(is.na(meta_sub[[sample_col]]))) {
    missing_mask <- is.na(meta_sub[[sample_col]])
    missing_samples <- sample_ids[missing_mask]
    
    stop(sprintf(
      "Found %d samples in expression matrix missing from metadata column '%s'. Examples: %s", 
      length(missing_samples), 
      sample_col, 
      paste(head(missing_samples, 3), collapse = ", ")
    ))
  }
  
  # FIX: Extra sanity check for 1:1 alignment
  if (nrow(meta_sub) != ncol(expr)) {
    stop("Internal error: metadata alignment failed (meta_sub rows != expr columns).")
  }
  
  # 7. Create Generic Annotation (for pheatmap)
  annot <- data.frame(
    Condition = meta_sub[[color_col]],
    row.names = sample_ids,
    stringsAsFactors = FALSE
  )
  
  list(
    expr = expr,
    meta = meta_sub, # Perfectly aligned with expr columns
    annot = annot,
    sample_col = sample_col,
    color_col = color_col,
    sample_ids = sample_ids
  )
}
#' Convert Expression Matrix to Long Format
#' 
#' Used for ggplot2. Optimized to use 'rep' since data is already aligned.
to_long_format <- function(prep_data) {
  
  n_features <- nrow(prep_data$expr)
  
  # Vectorized creation of long vectors
  df_long <- data.frame(
    sample = rep(prep_data$sample_ids, each = n_features),
    value  = as.vector(prep_data$expr),
    stringsAsFactors = FALSE
  )
  
  # OPTIMIZATION: Use rep() instead of match()
  # Since prep_data$meta is guaranteed to be sorted by sample_ids (from prepare_qc_data),
  # we can simply repeat each condition n_features times.
  cond_values <- prep_data$meta[[prep_data$color_col]]
  df_long[[prep_data$color_col]] <- rep(cond_values, each = n_features)
  
  # Filter non-finite values (NA/NaN/Inf) to keep plots clean
  df_long[is.finite(df_long$value), , drop = FALSE]
}



#------- MAIN PLOTTING FUNCTIONS ------- 


#' Boxplots of normalized expression per sample
norm_boxplot <- function(expr_norm, meta, cfg, out_file = NULL) {
  
  d <- prepare_qc_data(expr_norm, meta, cfg)
  df_long <- to_long_format(d)
  
  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = sample, y = value, colour = .data[[d$color_col]])
  ) +
    ggplot2::geom_boxplot(outlier.size = 0.4) +
    ggplot2::labs(
      title  = "Normalized expression boxplots",
      x      = "Sample",
      y      = "log2(normalized intensity)",
      colour = d$color_col
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
  
  if (!is.null(out_file)) {
    ggplot2::ggsave(out_file, plot = p, width = 8, height = 5)
  }
  
  invisible(p)
}


#' Histogram summary of normalized expression by condition
norm_histogram_summary <- function(expr_norm, meta, cfg, out_file = NULL) {
  
  d <- prepare_qc_data(expr_norm, meta, cfg)
  df_long <- to_long_format(d)
  
  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = value, fill = .data[[d$color_col]])
  ) +
    ggplot2::geom_histogram(alpha = 0.6, bins = 60, position = "identity") +
    ggplot2::facet_wrap(as.formula(paste("~", d$color_col)), nrow = 1, scales = "free_y") +
    ggplot2::labs(
      title = "Normalized expression histograms by condition",
      x     = "log2(normalized intensity)",
      y     = "Frequency",
      fill  = d$color_col
    ) +
    ggplot2::theme_minimal()
  
  if (!is.null(out_file)) {
    ggplot2::ggsave(out_file, plot = p, width = 8, height = 4)
  }
  
  invisible(p)
}

wrap_qc_heatmap <- function(expr_mat, meta, cfg, out_file = NULL) {
  
  # 1. Prepare Data
  # (Assuming prepare_qc_data logic is simple annotation creation)
  annot <- data.frame(
    Condition = meta[[cfg$effects$color]],
    row.names = meta[[cfg$effects$samples]]
  )
  
  # 2. Plot
  ph <- plot_heatmap_core(
    expr_mat = expr_mat,
    annotation_col = annot,
    title = "QC: Sample Protein Expression",
    max_rows = 2000,       # QC usually needs subsampling
    cluster_rows = TRUE,
    cluster_cols = TRUE
  )
  
  # 3. Save & Return
  if (!is.null(out_file)) save_heatmap_to_file(ph, out_file)
  return(ph)
}

#' Sample–sample correlation heatmap (proteomics QC)
qc_sample_correlation_heatmap <- function(expr_mat, meta, cfg, out_file, 
                                          method = "pearson", fontsize = 12) {
  d <- prepare_qc_data(expr_mat, meta, cfg)
  
  ph <- plot_sample_correlation_heatmap(
    expr_mat = d$expr,
    method = method,
    annotation_col = d$annot,
    fontsize = fontsize
  )
  
  save_heatmap_to_file(ph, out_file, width = 1600, height = 1200, res = 150)
  invisible(ph)
}


#' Sample–sample distance heatmap (QC)
qc_sample_distance_heatmap <- function(expr_mat, meta, cfg, out_file, 
                                       with_na = FALSE, fontsize = 12) {
  
  # FIX: Ensure matrix conversion happens BEFORE complete.cases logic
  # This prevents data.frame type coercion issues
  expr_mat <- as.matrix(expr_mat)
  
  if (!with_na) {
    keep <- stats::complete.cases(expr_mat)
    expr_mat <- expr_mat[keep, , drop = FALSE]
  }
  
  # Validate & Prepare (after filtering NAs)
  d <- prepare_qc_data(expr_mat, meta, cfg)
  
  ph <- plot_sample_distance_heatmap(
    expr_mat = d$expr,
    annotation_col = d$annot,
    fontsize = fontsize
  )
  
  save_heatmap_to_file(ph, out_file, width = 1600, height = 1200, res = 150)
  invisible(ph)
}


#' Density overlay of normalized expression
qc_proteomics_density <- function(expr_mat, meta, cfg, out_file = NULL,
                                  alpha = 0.3, show_legend = TRUE,
                                  title = "Density plot of normalized intensities") {
  
  # Keeping prepare_qc_data here intentionally.
  # Even though density doesn't strictly need metadata to run, 
  # we want to ensure the input data is consistent with the pipeline standards.
  d <- prepare_qc_data(expr_mat, meta, cfg)
  
  p <- plot_density_overlay(
    expr_mat = d$expr,
    title    = title,
    alpha    = alpha
  )
  
  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  if (!is.null(out_file)) {
    if (!grepl("\\.png$", out_file, ignore.case = TRUE)) {
      out_file <- paste0(out_file, ".png")
    }
    ggplot2::ggsave(out_file, p, width = 10, height = 6, dpi = 150)
  }
  
  invisible(p)
}

# (Existing function, kept for completeness)
imputed_histograms_summary <- function(imputed, imputed_flag, cfg = NULL, out_file = NULL) {
  p <- plot_imputation_summary(imputed, imputed_flag)
  if (!is.null(out_file)) {
    ggplot2::ggsave(out_file, plot = p, width = 12, height = 5, dpi = 150)
  }
  invisible(p)
}

#' Write per-sample imputation histograms (Batch writer)
#' Refactored: Uses ggsave and returns objects
write_imputation_histograms_per_sample <- function(expr_mat, imputed_flag, out_dir) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  sample_ids <- colnames(expr_mat)
  files <- character(0)
  plots <- list() 
  for (s in sample_ids) {
    p <- plot_imputation_histogram_one_sample(expr_mat, imputed_flag, s)
    f <- file.path(out_dir, paste0("X", s, ".png"))
    ggplot2::ggsave(filename = f, plot = p, width = 3.33, height = 6, dpi = 150)
    files <- c(files, f)
    plots[[paste0("sample_", s)]] <- p
  }
  
  return(list(files = files, plots = plots))
  
  # res_hist <- write_imputation_histograms_per_sample(    
  #                    pre$expr_imp_single, pre$imputation$imputed_flag, out_dir = out_qc)
  # files <- c(files, res_hist$files)
  # plots <- c(plots, res_hist$plots)
}

