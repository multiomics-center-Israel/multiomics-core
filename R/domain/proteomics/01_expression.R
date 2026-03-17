#' Apply custom annotation mapping to populate Genes and description columns
#'
#' @param row_data Data frame of feature annotations
#' @param cfg Proteomics mode config
#' @return Updated row_data with gene names and descriptions filled in
apply_custom_annotation <- function(row_data, cfg) {
  ann_cfg <- cfg$annotation %||% list()
  mapping_file <- ann_cfg$custom_mapping_file
  if (is.null(mapping_file) || !nzchar(mapping_file) || !file.exists(mapping_file)) {
    return(row_data)
  }
  
  protein_id_col <- cfg$id_columns$protein_id
  mapping <- read.csv(mapping_file, stringsAsFactors = FALSE)
  
  # Strip version suffix from protein IDs for matching (e.g. KAE8301209.1 -> KAE8301209)
  ids_stripped <- sub("\\.[0-9]+$", "", as.character(row_data[[protein_id_col]]))
  
  # Match to mapping$protein_id
  idx <- match(ids_stripped, mapping$protein_id)
  
  mapped_genes <- ifelse(!is.na(idx), mapping$gene_name[idx], NA_character_)
  mapped_desc  <- ifelse(!is.na(idx), mapping$description[idx], NA_character_)
  
  # Build display label: "gene_name | description" for readable display
  mapped_display <- ifelse(
    !is.na(mapped_genes) & !is.na(mapped_desc) & nzchar(mapped_desc),
    paste0(mapped_genes, " | ", mapped_desc),
    mapped_genes
  )
  
  # Only overwrite if existing values are empty or just "|"
  if ("Genes" %in% colnames(row_data)) {
    empty_genes <- is.na(row_data$Genes) | !nzchar(trimws(gsub("\\|", "", row_data$Genes)))
    row_data$Genes[empty_genes & !is.na(mapped_display)] <- mapped_display[empty_genes & !is.na(mapped_display)]
  } else {
    row_data$Genes <- mapped_display
  }
  
  if ("First.Protein.Description" %in% colnames(row_data)) {
    empty_desc <- is.na(row_data$First.Protein.Description) | !nzchar(trimws(row_data$First.Protein.Description))
    row_data$First.Protein.Description[empty_desc & !is.na(mapped_desc)] <- mapped_desc[empty_desc & !is.na(mapped_desc)]
  } else {
    row_data$First.Protein.Description <- mapped_desc
  }
  
  n_mapped <- sum(!is.na(mapped_genes))
  message(sprintf("Custom annotation: mapped %d / %d protein IDs to gene names.", n_mapped, nrow(row_data)))
  row_data
}

#' Build a standardized proteomics expression object (engine-aware)
#'
#' @param inputs List from load_proteomics_inputs().
#' @param config Full config list.
#' @return List with fields: assay_log2, assay_linear, row_data, col_data, info.
get_proteomics_expression_matrix <- function(inputs, config) {
  cfg <- config$modes$proteomics
  is_preprocessed <- identical(inputs$source_type, "preprocessed") ||
    identical(cfg$input$format, "preprocessed")
  
  # Resolved input scale for reporting / transform decisions
  scale_in <- cfg$scale_in %||% NULL
  
  if (is_preprocessed) {
    # ---- Preprocessed path: protein table is already an expression matrix ----
    protein_id_col <- cfg$id_columns$protein_id
    protein <- inputs$protein
    
    if (!protein_id_col %in% colnames(protein)) {
      stop(sprintf(
        "[proteomics preprocessed] Protein ID column '%s' not found. Available: %s",
        protein_id_col, paste(colnames(protein), collapse = ", ")
      ))
    }
    
    feat_ids <- as.character(protein[[protein_id_col]])
    
    # Identify numeric (sample) columns
    non_id_cols <- setdiff(colnames(protein), protein_id_col)
    is_numeric <- vapply(
      protein[, non_id_cols, drop = FALSE],
      function(x) is.numeric(x) || is.integer(x),
      logical(1)
    )
    anno_cols <- non_id_cols[!is_numeric]
    sample_cols <- non_id_cols[is_numeric]
    
    if (length(sample_cols) == 0) {
      stop("[proteomics preprocessed] No numeric sample columns found in protein table.")
    }
    
    # Build row_data from ID + annotation columns
    row_data <- protein[, c(protein_id_col, anno_cols), drop = FALSE]
    
    # Build expression matrix
    assay_log2 <- as.matrix(protein[, sample_cols, drop = FALSE])
    storage.mode(assay_log2) <- "numeric"
    rownames(assay_log2) <- feat_ids
    
    # Resolve input scale
    is_log_flag <- cfg$files$is_logtransformed %||% NULL
    
    if (!is.null(scale_in) && !is.null(is_log_flag)) {
      if ((identical(scale_in, "log2") && !isTRUE(is_log_flag)) ||
          (identical(scale_in, "linear") && isTRUE(is_log_flag))) {
        stop("Inconsistent proteomics config: scale_in conflicts with files$is_logtransformed")
      }
    }
    
    if (is.null(scale_in)) {
      scale_in <- if (isTRUE(is_log_flag)) "log2" else "linear"
    }
    
    if (identical(scale_in, "linear")) {
      message("[proteomics preprocessed] Applying log2(x + 1) transform")
      assay_log2 <- log2(assay_log2 + 1)
    } else if (!identical(scale_in, "log2")) {
      stop(sprintf(
        "[proteomics preprocessed] Unsupported scale_in '%s'. Expected 'linear' or 'log2'.",
        as.character(scale_in)
      ))
    }
    
    assay_linear <- NULL
    message(sprintf(
      "[proteomics preprocessed] %d features x %d samples (scale_in=%s)",
      nrow(assay_log2), ncol(assay_log2), scale_in
    ))
    
  } else if (cfg$engine == "DIANN") {
    # ---- DIANN engine path ----
    if (is.null(scale_in)) {
      scale_in <- "linear"
    }
    
    assay_log2 <- get_measurements_per_sample_diann(
      protein    = inputs$protein,
      sample_map = inputs$sample_map,
      meta       = inputs$metadata,
      cfg        = cfg,
      scale_in   = scale_in
    )
    assay_linear <- NULL
    
    # Feature annotations (row_data)
    row_data <- inputs$protein[, c(cfg$id_columns$protein_id, unlist(cfg$id_columns$protein_annot)), drop = FALSE]
    row_data <- apply_custom_annotation(row_data, cfg)
    
  } else {
    stop(sprintf("Unsupported proteomics engine: %s", cfg$engine))
  }
  
  # Align col_data to assay columns
  col_data <- inputs$metadata
  col_data <- align_meta_to_expr(assay_log2, inputs$metadata, cfg)
  
  list(
    assay_log2 = assay_log2,
    assay_linear = assay_linear,
    row_data = row_data,
    col_data = col_data,
    info = list(
      mode         = "proteomics",
      engine       = cfg$engine,
      scale_in     = scale_in %||% "log2",
      target_scale = cfg$transform$target_scale %||% "log2"
    )
  )
}

#' Extract DIA-NN measurements per sample (features × samples) and ensure log2 output
get_measurements_per_sample_diann <- function(protein, sample_map, meta, cfg, scale_in = "linear") {
  id_cols <- cfg$id_columns
  eff_cols <- cfg$effects
  
  # 0) Protein ID for rownames
  protein_id_col <- id_cols$protein_id
  check_has_cols(protein, protein_id_col, df_name = "protein")
  feat_ids <- as.character(protein[[protein_id_col]])
  
  if (anyNA(feat_ids) || any(feat_ids == "")) {
    stop(sprintf("Protein ID column '%s' contains NA/empty values.", protein_id_col))
  }
  if (anyDuplicated(feat_ids) > 0) {
    dups <- unique(feat_ids[duplicated(feat_ids)])
    stop(sprintf(
      "Protein ID column '%s' contains duplicated IDs (e.g. %s).",
      protein_id_col, paste(head(dups, 3), collapse = ", ")
    ))
  }
  
  # 1) Identify non-sample columns in the protein table
  annot_cols <- id_cols$protein_annot
  annot_cols <- if (is.null(annot_cols)) character(0) else unlist(annot_cols)
  non_sample_cols <- unique(c(protein_id_col, annot_cols))
  
  # Sample columns = everything else
  sample_cols <- setdiff(colnames(protein), non_sample_cols)
  if (length(sample_cols) == 0) {
    stop("No sample columns detected in protein table after excluding ID/annotation columns.")
  }
  
  df_m <- protein[, sample_cols, drop = FALSE]
  
  # 2) Coerce to numeric (DIA-NN often has blanks)
  df_m <- as.data.frame(
    suppressWarnings(sapply(df_m, function(x) as.numeric(as.character(x)))),
    check.names = FALSE
  )
  
  # 3) Transform only if needed
  if (identical(scale_in, "linear")) {
    if (any(df_m <= 0, na.rm = TRUE)) {
      warning("Non-positive values found before log2 on DIA-NN input; will produce -Inf/NaN")
    }
    message("[proteomics DIANN] Applying log2 transform")
    df_m <- as.data.frame(
      log2(df_m),
      check.names = FALSE
    )
  } else if (identical(scale_in, "log2")) {
    message("[proteomics DIANN] Input declared as already log2; skipping log2 transform")
  } else {
    stop(sprintf(
      "[proteomics DIANN] Unsupported scale_in '%s'. Expected 'linear' or 'log2'.",
      as.character(scale_in)
    ))
  }
  
  # 4) Rename columns: raw sample names -> SampleID using sample_map
  if (!is.null(sample_map)) {
    df_m <- apply_sample_map_to_colnames(df_m, sample_map, id_cols$map_from, id_cols$map_to)
  }
  # else: no sample_map — column names are already the sample IDs
  
  # 5) Reorder columns to match metadata sample order
  meta_sample_col <- eff_cols$samples
  check_has_cols(meta, meta_sample_col, df_name = "metadata")
  meta_sample_ids <- meta[[meta_sample_col]]
  
  missing_in_df <- setdiff(meta_sample_ids, colnames(df_m))
  if (length(missing_in_df) > 0) {
    warning(
      "These samples from metadata were not found in DIA-NN columns after renaming: ",
      paste(head(missing_in_df, 10), collapse = ", ")
    )
  }
  
  ordered_cols <- intersect(meta_sample_ids, colnames(df_m))
  df_m_ordered <- df_m[, ordered_cols, drop = FALSE]
  
  # 6) Return as numeric matrix with proper rownames
  coerce_df_to_numeric_matrix(df_m_ordered, rownames_vec = feat_ids, name = "diann_measurements")
}