write_rnaseq_outputs_legacy <- function(pre, de_res, inputs, config, out_dir) {
  files <- character(0)
  dirs  <- create_legacy_output_dirs(out_dir)
  
  # 1) datasets
  files <- c(files, write_rnaseq_datasets_legacy(pre, config, out_dir))
  
  # 2) summary_df (built here, not inside de_res)
  summary_df <- build_rnaseq_summary_df(de_res, config)
  files <- c(files, save_tsv(summary_df, dirs$datasets, sprintf("deseq2_summary_p%s.tsv", p_tag(config))))
  
  # 3) final results TSV
  if (!is.null(inputs$contrasts)) {
    final_results <- build_final_results_rnaseq(
      pre          = pre,
      summary_df   = summary_df,
      contrasts_df = inputs$contrasts,
      row_data     = pre$row_data
    )
    files <- c(files, save_tsv(final_results, dirs$datasets, "final_results.tsv"))
    
    # 4) Excel outputs (optional; reuse legacy writer with small adapter)
    files <- c(files, write_final_results_excels_legacy_rna(final_results, pre, config, out_dir))
  }
  
  unique(files)
}


#' Writes intermediate RNA matrices (legacy-like)
write_rnaseq_datasets_legacy <- function(pre, config, out_dir) {
  dirs <- create_legacy_output_dirs(out_dir)
  files <- character(0)
  
  # raw filtered counts
  files <- c(files, save_tsv(
    as.data.frame(pre$expr_filt, check.names = FALSE),
    dirs$datasets,
    "rna_counts_filtered.tsv"
  ))
  
  # expr_work (for QC/DE; usually logCPM/VST)
  files <- c(files, save_tsv(
    as.data.frame(pre$expr_work, check.names = FALSE),
    dirs$datasets,
    sprintf("rna_expr_work_%s.tsv", attr(pre$expr_work, "method") %||% "norm")
  ))
  
  unique(files)
}

build_rnaseq_summary_df <- function(de_res, config) {
  stopifnot(is.list(de_res), !is.null(de_res$tables), !is.null(de_res$contrasts))
  
  contrasts_df <- de_res$contrasts
  stopifnot(all(c("Contrast_name","Factor","Numerator","Denominator") %in% colnames(contrasts_df)))
  
  # cutoffs (אם אין בconfig, ברירות מחדל)
  de_cfg <- config$modes$rna$de %||% list()
  padj_cutoff <- as.numeric(de_cfg$padj_cutoff %||% 0.05)
  linear_fc_cutoff <- as.numeric(de_cfg$linear_fc_cutoff %||% 1.5)
  
  # Anchor FeatureID from first contrast table
  cn0 <- contrasts_df$Contrast_name[[1]]
  tab0 <- de_res$tables[[cn0]]
  if (is.null(tab0)) stop("RNA DE tables missing first contrast: ", cn0)
  if (!"FeatureID" %in% colnames(tab0)) stop("RNA DE table missing FeatureID column for: ", cn0)
  
  base_ids <- tab0$FeatureID
  out <- data.frame(FeatureID = base_ids, stringsAsFactors = FALSE, check.names = FALSE)
  
  pass_counts <- integer(length(base_ids))
  
  for (i in seq_len(nrow(contrasts_df))) {
    cn <- contrasts_df$Contrast_name[i]
    tab <- de_res$tables[[cn]]
    if (is.null(tab)) stop("RNA DE tables missing contrast: ", cn)
    
    # align rows by FeatureID
    tab <- tab[match(base_ids, tab$FeatureID), , drop = FALSE]
    
    log2fc <- tab$log2FoldChange
    pval   <- tab$pvalue
    padj   <- tab$padj
    
    # legacy-style signed linearFC:
    #   if log2fc>0 => 2^log2fc
    #   if log2fc<0 => -1/(2^log2fc)
    linearFC <- ifelse(
      is.na(log2fc), NA_real_,
      ifelse(log2fc > 0, 2^log2fc, -1/(2^log2fc))
    )
    
    # magnitude ratio (always positive)
    linearRatio <- ifelse(is.na(log2fc), NA_real_, 2^log2fc)
    
    pass_dir <- ifelse(
      !is.na(padj) &
        (padj <= padj_cutoff) &
        (abs(linearFC) >= linear_fc_cutoff),
      ifelse(linearFC > 0, "up", "down"),
      NA_character_
    )
    
    out[[paste0("sum.pass.", cn)]] <- ifelse(is.na(pass_dir), 0, 1)
    
    # keep proteomics naming for compatibility
    out[[paste0("pass.imputs.", cn)]] <- pass_dir
    out[[paste0("linearRatio.imputs.", cn)]] <- signif(linearRatio, 6)
    out[[paste0("linearFC.imputs.", cn)]] <- signif(linearFC, 6)
    out[[paste0("pvalue.imputs.", cn)]] <- signif(pval, 6)
    out[[paste0("padj.imputs.", cn)]] <- signif(padj, 6)
    
    pass_counts <- pass_counts + ifelse(is.na(pass_dir), 0L, 1L)
  }
  
  out$n_pass_contrasts <- pass_counts
  out$pass_any_contrast <- ifelse(pass_counts > 0, 1L, NA_integer_)
  
  out
}

build_final_results_rnaseq <- function(pre, summary_df, contrasts_df, row_data = NULL) {
  expr_df <- as.data.frame(pre$expr_work, check.names = FALSE)
  
  # base ids
  gene_ids <- summary_df$FeatureID
  
  base <- data.frame(
    Gene = gene_ids,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Optional row_data annotations (אם בעתיד תוסיפי)
  if (!is.null(row_data)) {
    for (col in intersect(c("gene_name","symbol","description","GeneSymbol","Description"), colnames(row_data))) {
      # try match by rownames or a gene_id column
      if (!is.null(rownames(row_data)) && all(gene_ids %in% rownames(row_data))) {
        base[[col]] <- row_data[gene_ids, col]
      }
    }
  }
  
  # Add expression (aligned by rownames)
  if (!is.null(rownames(pre$expr_work))) {
    expr_df2 <- expr_df[match(gene_ids, rownames(pre$expr_work)), , drop = FALSE]
  } else {
    stop("pre$expr_work must have rownames for final results.")
  }
  base <- cbind(base, expr_df2)
  
  # Add contrast stats (same naming as proteomics helper)
  contrast_names <- contrasts_df$Contrast_name
  
  for (cn in contrast_names) {
    cols <- get_contrast_cols(cn)
    needed <- c(cols$fc, cols$p, cols$padj)
    missing_cols <- setdiff(needed, colnames(summary_df))
    if (length(missing_cols) > 0) {
      stop(sprintf("Summary DF missing columns for contrast '%s': %s", cn, paste(missing_cols, collapse = ", ")))
    }
  }
  
  m <- match(base$Gene, summary_df$FeatureID)
  
  for (cn in contrast_names) {
    cols <- get_contrast_cols(cn)
    
    fc_vals <- summary_df[[cols$fc]][m]
    pass_vals <- if (cols$pass %in% colnames(summary_df)) summary_df[[cols$pass]][m] else rep(NA, length(m))
    
    base[[cols$fc]]   <- fc_vals
    base[[cols$p]]    <- summary_df[[cols$p]][m]
    base[[cols$padj]] <- summary_df[[cols$padj]][m]
    
    base[[cols$updown]] <- ifelse(!is.na(pass_vals),
                                  ifelse(as.numeric(fc_vals) >= 0, "up", "down"), "")
    base[[cols$manual]] <- NA
  }
  
  # pass_any_contrast
  pass_cols <- paste0("pass.imputs.", contrast_names)
  existing_pass_cols <- intersect(pass_cols, colnames(summary_df))
  
  if (length(existing_pass_cols) == 0) {
    warning("No 'pass.imputs' columns found in summary_df. pass_any_contrast will be NA.")
    base$pass_any_contrast <- NA
  } else {
    pass_mat <- summary_df[m, existing_pass_cols, drop = FALSE]
    base$pass_any_contrast <- ifelse(rowSums(!is.na(pass_mat)) > 0, 1, NA)
  }
  
  base
}

