p_tag <- function(config, default = "NA") {
  p <- config$modes$proteomics$de$p_cutoff
  if (is.null(p) || is.na(p)) return(default)
  format(p, trim = TRUE, scientific = FALSE)
}

#' Run limma differential analysis for proteomics with contrast support
#'
#' Fits a limma linear model on an imputed proteomics expression matrix and
#' returns per-contrast result tables with feature annotations.
#'
#' @return A list with aligned metadata, design matrix, contrasts, fitted model and per-contrast DE tables.
#' @export
run_limma_proteomics <- function(expr_imp, meta, contrasts_df, prot_tbl, cfg) {
  p_cfg <- cfg$modes$proteomics
  
  sample_col <- p_cfg$effects$samples %||% "SampleID"
  group_col  <- p_cfg$effects$color   %||% "Condition"
  
  protein_id_col <- p_cfg$id_columns$protein_id %||% "Protein.Group"
  default_annot <- c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description")
  annot_cols <- unique(c(protein_id_col, p_cfg$id_columns$protein_annot %||% default_annot))
  
  assert_numeric_matrix(expr_imp, "expr_imp")
  meta_aligned <- align_meta_to_expr(expr_imp, meta, p_cfg)
  
  meta_aligned[[group_col]] <- factor(meta_aligned[[group_col]])
  design <- model.matrix(stats::as.formula(paste0("~ 0 + ", group_col)), data = meta_aligned)
  colnames(design) <- levels(meta_aligned[[group_col]])
  
  stopifnot(all(contrasts_df$Factor == group_col))
  contrast_formulas <- setNames(
    paste(contrasts_df$Numerator, contrasts_df$Denominator, sep = " - "),
    contrasts_df$Contrast_name
  )
  
  contrast_matrix <- limma::makeContrasts(contrasts = contrast_formulas, levels = design)
  colnames(contrast_matrix) <- names(contrast_formulas)
  
  fit2 <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr_imp, design), contrast_matrix))
  
  ann <- align_annotations_to_expr(expr_imp, prot_tbl, protein_id_col, annot_cols)
  feature_id <- ann[[protein_id_col]]
  annot_out <- setdiff(annot_cols, protein_id_col)
  
  de_tables <- lapply(colnames(contrast_matrix), function(cn) {
    de <- limma::topTable(fit2, coef = cn, adjust.method="BH", sort.by="none", number=Inf)
    de <- align_de_to_expr(de, expr_imp, contrast_name = cn)
    data.frame(
      FeatureID = feature_id,
      Contrast  = cn,
      ann[, annot_out, drop = FALSE],
      de[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")],
      check.names = FALSE
    )
  })
  
  # IMPORTANT: names(de_tables) are used by summarize_limma_mult_imputation()
  names(de_tables) <- colnames(contrast_matrix)
  
  list(
    expr_imp = expr_imp,
    meta_aligned = meta_aligned,
    design = design,
    contrast_formulas = contrast_formulas,
    contrast_matrix = contrast_matrix,
    fit2 = fit2,
    de_tables = de_tables
  )
}


#' Run limma proteomics DE on a precomputed list of imputed datasets
#'
#' Use this function when imputations were already generated and you want to:
#' - avoid recomputing imputations (e.g., when trying a different DE method), or
#' - decouple imputation from DE in a `{targets}` pipeline.
#'
#' @return A list of length `length(imputations)`. Each element is the output of
#'         `run_limma_proteomics()` for the corresponding imputed dataset.

run_limma_multimp <- function(imputations, meta, contrasts_df, prot_tbl, cfg,
                              verbose = FALSE) {
  
  validate_proteomics_imputations(imputations = imputations,
                                  meta = meta, cfg = cfg)
  
  runs <- vector("list", length(imputations))
  
  for (i in seq_along(imputations)) {
    if (isTRUE(verbose)) {
      message(sprintf("Limma on imputations: %d / %d", i, length(imputations)))
    }
    
    runs[[i]] <- run_limma_proteomics(
      expr_imp     = imputations[[i]],
      meta         = meta,
      contrasts_df = contrasts_df,
      prot_tbl     = prot_tbl,
      cfg          = cfg
    )
  }
  
  runs
}

.build_contrasts_matrix <- function(contrasts_df, design, factor_col) {
  # contrasts_df: Contrast_name, Factor, Numerator, Denominator
  stopifnot(all(c("Contrast_name", "Factor", "Numerator", "Denominator") %in% names(contrasts_df)))
  
  contrasts_df <- contrasts_df[contrasts_df$Factor == factor_col, , drop = FALSE]
  if (nrow(contrasts_df) == 0) stop("No contrasts match Factor='", factor_col, "'.")
  
  # design columns look like: factor_colLevel (because ~0 + factor)
  levs <- sub(paste0("^", factor_col), "", colnames(design))
  
  make_one <- function(num, den) {
    num_col <- paste0(factor_col, num)
    den_col <- paste0(factor_col, den)
    if (!num_col %in% colnames(design)) stop("Numerator level not in design: ", num)
    if (!den_col %in% colnames(design)) stop("Denominator level not in design: ", den)
    v <- rep(0, ncol(design)); names(v) <- colnames(design)
    v[num_col] <- 1
    v[den_col] <- -1
    v
  }
  
  cm <- vapply(seq_len(nrow(contrasts_df)), function(i) {
    make_one(contrasts_df$Numerator[i], contrasts_df$Denominator[i])
  }, FUN.VALUE = numeric(ncol(design)))
  
  cm <- t(cm)
  rownames(cm) <- contrasts_df$Contrast_name
  cm
}

#' Extract DE feature IDs from a DE result object
#'
#' Omics-agnostic helper that returns feature identifiers passing
#' differential expression thresholds as defined in config.
#'
#' @param de_res DE result object containing a summary table
#' @param config Global config
#'
#' @return character vector of feature IDs
get_de_features <- function(de_res, cfg) {
  stopifnot(is.list(de_res), !is.null(de_res$summary_df))
  df <- de_res$summary_df
  stopifnot(is.data.frame(df))
  
  # Accept legacy + new column names
  id_col <- cfg$de_table$id_col
  pass_any_col <- cfg$de_table$pass_any_col
  if (is.na(id_col) || is.null(id_col)) {
    stop("DE summary missing id column. Expected : ", id_col)
  }
  
  # Standardize internally
  df$feature_id <- df[[id_col]]
  
  feats <- df$feature_id[!is.na(df[[pass_any_col]]) & df[[pass_any_col]] == 1]
  unique(as.character(feats))
}

#' Standardized per-contrast column names in the DE summary table
#'
#' This project uses a "wide" summary format where each contrast creates
#' multiple columns (fc/p/padj/pass/upDown/manual_cutoffs).
#'
#' @param contrast Character scalar. Contrast name (e.g., "T_vs_C").
#' @return Named list with standardized column names.
get_contrast_cols <- function(contrast) {
  stopifnot(is.character(contrast), length(contrast) == 1, nzchar(contrast))
  
  list(
    fc     = paste0("linearFC.imputs.", contrast),
    p      = paste0("pvalue.imputs.",   contrast),
    padj   = paste0("padj.imputs.",     contrast),
    pass   = paste0("pass.imputs.",     contrast),
    updown = paste0("upDown.imputs.",   contrast),
    manual = paste0("manual_cutoffs.",  contrast)
  )
}

signed_fc_to_log2 <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x) & x != 0
  out[ok] <- sign(x[ok]) * log2(abs(x[ok]))
  out
}
assert_de_contract <- function(de_res, stage = "proteomics") {
  stopifnot(is.list(de_res))
  
  required <- c("summary_df")
  missing <- setdiff(required, names(de_res))
  if (length(missing) > 0) {
    stop(sprintf(
      "DE contract failed for %s. Missing fields: %s",
      stage, paste(missing, collapse = ", ")
    ))
  }
  
  invisible(TRUE)
}


assert_de_contract <- function(de_res, stage = "proteomics") {
  stopifnot(is.list(de_res))
  
  required <- c("summary_df")
  missing <- setdiff(required, names(de_res))
  if (length(missing) > 0) {
    stop(sprintf(
      "DE contract failed for %s. Missing fields: %s",
      stage, paste(missing, collapse = ", ")
    ))
  }
  
  invisible(TRUE)
}

#' Mark differential expression pass for a single result table
#'
#' Applies the legacy "pass1" rule: a feature passes if it meets BOTH
#' the p-value cutoff (raw or adjusted) and the absolute log2FC cutoff.
#'
#' Output encoding matches the old script:
#'   - 1  : pass
#'   - NA : fail
#'
#' @param de_table A data.frame from topTable()
#' @param p_cutoff Numeric cutoff for p-value or adjusted p-value.
#' @param lfc_cutoff Numeric cutoff for |logFC| (log2 scale).
#' @param use_adj Logical; TRUE uses 'adj.P.Val', FALSE uses 'P.Value'.
#'
#' @return Numeric vector of length nrow(de_table) with values {1, NA}.
mark_pass1 <- function(de_table, p_cutoff, lfc_cutoff, use_adj = TRUE) {
  pcol <- if (isTRUE(use_adj)) "adj.P.Val" else "P.Value"
  stopifnot(all(c("logFC", pcol) %in% colnames(de_table)))
  
  pass <- (de_table[[pcol]] <= p_cutoff) &
    (abs(de_table[["logFC"]]) >= lfc_cutoff)
  
  ifelse(pass, 1, NA)
}

#' Helper to add pass_any_contrast column
#' 
#' Scans for columns starting with 'pass.imputs.' and creates a summary flag.
#'
#' @param summary_df Data frame containing pass columns.
#' @param pass_prefix Regex prefix for pass columns.
#' @param out_col Name of the output binary column.
#' @param out_n_col Name of the output count column.
#' 
#' @return The data frame with added columns.
add_pass_any_contrast <- function(summary_df,
                                  pass_prefix = "^pass\\.imputs\\.",
                                  out_col = "pass_any_contrast",
                                  out_n_col = "n_pass_contrasts") {
  
  pass_cols <- grep(pass_prefix, names(summary_df), value = TRUE)
  
  if (length(pass_cols) == 0) {
    # If no pass columns found, initialize with 0/NA to avoid downstream errors
    summary_df[[out_n_col]] <- 0L
    summary_df[[out_col]]   <- NA
    return(summary_df)
  }
  
  pass_mat <- as.matrix(summary_df[, pass_cols, drop = FALSE])
  
  # Count how many contrasts passed (ignoring NAs)
  n_pass <- rowSums(!is.na(pass_mat) & pass_mat == 1, na.rm = TRUE)
  
  summary_df[[out_n_col]] <- as.integer(n_pass)
  summary_df[[out_col]]   <- ifelse(n_pass > 0, 1, NA)
  
  summary_df
}

#' Summarize differential expression across multiple imputations (legacy-compatible)
#'
#' Reproduces the old script logic:
#'   - "stable" if passed in >= min_no_passed imputations
#'   - average linear ratio across imputations and convert to signed linearFC
#'   - summarize p-value and padj using quantile
#'
#' @param runs_de_tables List of length cfg$no_repetitions.
#' @param config Full config list.
#'
#' @return A data.frame with one row per feature.
summarize_limma_mult_imputation <- function(runs_de_tables, config) {
  de_cfg <- config$modes$proteomics$de
  imp_cfg   <- config$modes$proteomics$imputation
  
  NO_REPETITIONS <- as.integer(imp_cfg$no_repetitions)
  MIN_NO_PASSED  <- as.integer(imp_cfg$min_no_passed)
  
  use_adj_for_pass1 <- isTRUE(de_cfg$use_adj_for_pass1)
  p_cutoff <- as.numeric(de_cfg$p_cutoff)
  
  linear_fc_cutoff <- as.numeric(de_cfg$linear_fc_cutoff)
  lfc_cutoff <- log2(linear_fc_cutoff)
  
  stopifnot(length(runs_de_tables) == NO_REPETITIONS)
  stopifnot(MIN_NO_PASSED >= 1, MIN_NO_PASSED <= NO_REPETITIONS)
  
  q <- MIN_NO_PASSED / NO_REPETITIONS
  
  contrasts <- names(runs_de_tables[[1]])
  stopifnot(length(contrasts) > 0)
  
  id_cols <- c("FeatureID", "Protein.Names", "Genes", "First.Protein.Description")
  ref_df  <- runs_de_tables[[1]][[contrasts[1]]]
  stopifnot(all(id_cols %in% colnames(ref_df)))
  
  out <- ref_df[, id_cols, drop = FALSE]
  ref_ids <- ref_df$FeatureID
  
  # Validate alignment
  for (n in seq_len(NO_REPETITIONS)) {
    for (cn in contrasts) {
      cur <- runs_de_tables[[n]][[cn]]
      stopifnot(nrow(cur) == length(ref_ids))
      stopifnot(all(cur$FeatureID == ref_ids))
    }
  }
  
  # Loop over contrasts to calculate stability
  for (cn in contrasts) {
    contrast_print <- gsub(" ", "", cn)
    
    logfc_mat <- sapply(seq_len(NO_REPETITIONS),
                        function(n) runs_de_tables[[n]][[cn]][["logFC"]])
    p_mat     <- sapply(seq_len(NO_REPETITIONS),
                        function(n) runs_de_tables[[n]][[cn]][["P.Value"]])
    padj_mat  <- sapply(seq_len(NO_REPETITIONS),
                        function(n) runs_de_tables[[n]][[cn]][["adj.P.Val"]])
    
    pass1_mat <- sapply(seq_len(NO_REPETITIONS), function(n)
      mark_pass1(runs_de_tables[[n]][[cn]],
                 p_cutoff    = p_cutoff,
                 lfc_cutoff  = lfc_cutoff,
                 use_adj     = use_adj_for_pass1)
    )
    
    sum_pass    <- rowSums(pass1_mat, na.rm = TRUE)
    pass_imputs <- ifelse(sum_pass >= MIN_NO_PASSED, 1, NA)
    
    linearRatio_imputs <- rowMeans(2^logfc_mat, na.rm = TRUE)
    linearFC_imputs <- ifelse(linearRatio_imputs >= 1,
                              linearRatio_imputs,
                              -1 / linearRatio_imputs)
    
    pvalue_imputs <- apply(p_mat,   1, quantile, probs = q, na.rm = TRUE)
    padj_imputs   <- apply(padj_mat,1, quantile, probs = q, na.rm = TRUE)
    
    out[[paste0("sum.pass.", contrast_print)]]           <- sum_pass
    out[[paste0("pass.imputs.", contrast_print)]]        <- pass_imputs
    out[[paste0("linearRatio.imputs.", contrast_print)]] <- linearRatio_imputs
    out[[paste0("linearFC.imputs.", contrast_print)]]    <- signif(linearFC_imputs, 4)
    out[[paste0("pvalue.imputs.", contrast_print)]]      <- pvalue_imputs
    out[[paste0("padj.imputs.", contrast_print)]]        <- padj_imputs
  }
  
  # Use the helper to add the global pass flags (DRY principle)
  out <- add_pass_any_contrast(out)
  
  return(out)
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


#' Writes intermediate proteomics matrices
write_proteomics_datasets_legacy <- function(pre, runs = NULL, config, out_dir) {
  dirs <- create_legacy_output_dirs(out_dir)
  cfg  <- config$modes$proteomics
  files <- character(0)
  
  files <- c(files, save_tsv(as.data.frame(pre$expr_filt, check.names = FALSE),
                             dirs$datasets, "protein_log2_filtered_unimputed.tsv"))
  
  fname_imp <- sprintf(
    "protein_log2_filtered_imputed_once_width_%s_shift_%s.tsv",
    cfg$imputation$width, cfg$imputation$downshift
  )
  
  files <- c(files, save_tsv(as.data.frame(pre$expr_imp_single,  check.names = FALSE),
                             dirs$datasets, fname_imp))
  
  
  if (!is.null(runs) && length(runs) > 0) {
    rep_dir <- file.path(dirs$datasets, "imputed_repetitions")
    ensure_dir(rep_dir)
    
    for (i in seq_along(runs)) {
      expr_i <- runs[[i]]$expr_imp
      if (is.null(expr_i)) next
      
      # repetitions
      files <- c(files, save_tsv(as.data.frame(expr_i, check.names = FALSE),
                                 rep_dir, sprintf("protein_log2_filtered_imputed_%02d.tsv", i)))
      
    }
  }
  
  unique(files)
}

#' Write legacy-style  multi-imputation summary
write_limma_multimp_summary_legacy <- function(summary_df, config, out_dir) {
  dirs <- create_legacy_output_dirs(out_dir)
  save_tsv(summary_df, dirs$datasets, sprintf("limma_multimp_summary_p%s.tsv", p_tag(config)))
}

#' Build legacy-style wide limma table across imputations
build_limma_results_multimp_wide <- function(runs_de_tables, contrast_name, 
                                             stats_cols = c("logFC", "P.Value", "adj.P.Val")) {
  stopifnot(length(runs_de_tables) >= 1)
  
  # Anchor with first run
  base <- runs_de_tables[[1]][[contrast_name]]
  
  # FIX: Fail-fast if first run is broken
  if (is.null(base)) stop(sprintf("Contrast '%s' missing in imputation run 1", contrast_name))
  stopifnot("FeatureID" %in% colnames(base))
  
  id_cols <- intersect(c("FeatureID", "Protein.Names", "Genes", "First.Protein.Description", "Contrast"), colnames(base))
  out <- base[, id_cols, drop = FALSE]
  
  # Append stats per imputation
  for (i in seq_along(runs_de_tables)) {
    tab <- runs_de_tables[[i]][[contrast_name]]
    
    # FIX: Strict check - missing imputation data is critical
    if (is.null(tab)) {
      stop(sprintf("Critical: Contrast '%s' is missing in imputation run %d (but existed in run 1)", contrast_name, i))
    }
    
    # Ensure alignment
    tab <- align_de_table_by_feature_id(
      tab = tab,
      ref_ids = out$FeatureID,
      run_i = i,
      contrast_name = contrast_name,
      id_col = "FeatureID"
    )
    
    
    stat_block <- tab[, intersect(stats_cols, colnames(tab)), drop = FALSE]
    colnames(stat_block) <- paste0(colnames(stat_block), ".", i)
    out <- cbind(out, stat_block)
  }
  
  out
}

#' Write legacy-style limma results across multiple imputations
write_limma_results_multimp_legacy <- function(de_res, contrast_name, config, out_dir) {
  dirs <- create_legacy_output_dirs(out_dir)
  
  wide_df <- build_limma_results_multimp_wide(
    runs_de_tables = de_res$runs_de_tables,
    contrast_name  = contrast_name
  )
  
  fname <- sprintf("limma_results_multimp_p%s.tsv", p_tag(config))
  save_tsv(wide_df, dirs$datasets, fname)
}


#' Build legacy-style Final_results table
build_final_results_proteomics <- function(pre, summary_df, contrasts_df, row_data = NULL) {
  
  expr_df <- as.data.frame(pre$expr_filt, check.names = FALSE)
  if (is.null(row_data)) row_data <- pre$row_data
  
  stopifnot(!is.null(row_data), "Protein.Group" %in% colnames(row_data))
  
  # Initialize Base
  base <- data.frame(
    Protein = row_data[["Protein.Group"]],
    stringsAsFactors = FALSE, check.names = FALSE
  )
  
  # Match to Summary
  m <- match(base$Protein, summary_df$FeatureID)
  
  # Add Annotations
  for (col in c("Protein.Names", "Genes", "First.Protein.Description")) {
    val <- if (col %in% colnames(summary_df)) summary_df[[col]][m] else NA
    if (col %in% colnames(row_data)) {
      fallback <- row_data[[col]]
      val <- ifelse(is.na(val), fallback, val)
    }
    base[[col]] <- val
  }
  
  # Add Expression
  base <- cbind(base, expr_df)
  
  # Add Contrast Stats
  contrast_names <- contrasts_df$Contrast_name
  
  # FIX: Pre-validate that all required columns exist in summary_df
  # This prevents partial failure inside the loop
  for (cn in contrast_names) {
    cols <- get_contrast_cols(cn)
    needed <- c(cols$fc, cols$p, cols$padj) # pass is optional-ish but usually needed
    missing_cols <- setdiff(needed, colnames(summary_df))
    if (length(missing_cols) > 0) {
      stop(sprintf("Summary DF missing columns for contrast '%s': %s", cn, paste(missing_cols, collapse=", ")))
    }
  }
  
  for (cn in contrast_names) {
    cols <- get_contrast_cols(cn) 
    
    fc_vals   <- summary_df[[cols$fc]][m]
    # FIX: Robust check for 'pass' column existence
    pass_vals <- if (cols$pass %in% colnames(summary_df)) summary_df[[cols$pass]][m] else rep(NA, length(m))
    
    base[[cols$fc]]   <- fc_vals
    base[[cols$p]]    <- summary_df[[cols$p]][m]
    base[[cols$padj]] <- summary_df[[cols$padj]][m]
    
    # Calculate Up/Down logic
    base[[cols$updown]] <- ifelse(!is.na(pass_vals),
                                  ifelse(as.numeric(fc_vals) >= 0, "up", "down"), "")
    base[[cols$manual]] <- NA 
  }
  
  # FIX: Robust pass_any_contrast logic
  pass_cols <- paste0("pass.imputs.", contrast_names)
  existing_pass_cols <- intersect(pass_cols, colnames(summary_df))
  
  if (length(existing_pass_cols) == 0) {
    warning("No 'pass.imputs' columns found in summary_df. pass_any_contrast will be NA.")
    base$pass_any_contrast <- NA
  } else {
    # Create matrix, handle potential NAs in 'm' automatically (returns NA row)
    pass_mat <- summary_df[m, existing_pass_cols, drop = FALSE]
    # rowSums ignoring NAs, check if > 0
    # Note: !is.na(NA) is FALSE, so NA entries don't contribute to the sum
    base$pass_any_contrast <- ifelse(rowSums(!is.na(pass_mat)) > 0, 1, NA)
  }
  
  base
}

#' Orchestrator: write all proteomics multi-imputation outputs (legacy-compatible)
#'
#' @param pre Preprocessed object
#' @param de_res List with $runs_de_tables, $runs, $summary_df
#' @param inputs proteomics inputs (expects $contrasts at least)
#' @param config global config
#' @param out_dir run directory
#' @param write_runs logical; if TRUE writes per-imputation matrices (if available in de_res$runs)
#'
#' @return character vector of written file paths
write_proteomics_multimpute_outputs <- function(pre, de_res,
                                                inputs, config, out_dir,excel_order = NULL,
                                                write_runs = FALSE) {
  
  files <- character(0)
  dirs  <- create_legacy_output_dirs(out_dir)
  
  # 1) datasets
  runs_for_datasets <- if (isTRUE(write_runs)) de_res$runs else NULL
  files <- c(files, write_proteomics_datasets_legacy(pre, runs_for_datasets, config, out_dir))
  
  # 2) summary
  if (!is.null(de_res$summary_df)) {
    files <- c(files, write_limma_multimp_summary_legacy(de_res$summary_df, config, out_dir))
  }
  
  # 3) wide limma per contrast (optional)
  if (!is.null(de_res$runs_de_tables) && length(de_res$runs_de_tables) > 0) {
    contrast_names <- names(de_res$runs_de_tables[[1]])
    for (cn in contrast_names) {
      files <- c(files, write_limma_results_multimp_legacy(de_res = de_res, contrast_name = cn, config = config, out_dir = out_dir))
    }
  }
  
  # 4) final results TSV (optional but useful)
  if (!is.null(inputs$contrasts) && !is.null(de_res$summary_df)) {
    final_results <- build_final_results_proteomics(
      pre          = pre,
      summary_df   = de_res$summary_df,
      contrasts_df = inputs$contrasts,
      row_data     = pre$row_data
    )
    files <- c(files, save_tsv(final_results, dirs$datasets, "final_results.tsv"))
    
    # 5) Excel outputs (delegated to dedicated file)
    # files <- c(files, write_final_results_excels_proteomics(final_results, pre, config, out_dir))
    files <- c(files, write_final_results_excels_proteomics(
      final_results = final_results,
      pre          = pre,
      config       = config,
      out_dir      = out_dir,
      excel_order  = excel_order   # NEW
    ))
    
  }
  
  unique(files)
}

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
    files <- c(files, write_final_results_excels_legacy_rna(final_results,pre, config, out_dir))
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

#' Add legacy Cutoffs sheet + named regions (PVAL_CO/FDR_CO/LFC_CO)
#'
#' This mirrors the old pipeline exactly:
#' - labels at Cutoffs!B4:B6
#' - values at Cutoffs!C4:C6
#' - named regions point to C4/C5/C6
#' - green/red styles mark which cutoff is used (p-value vs FDR)
#' Add legacy Cutoffs sheet + named regions (PVAL_CO/FDR_CO/LFC_CO)
#'
#' Generic version for any omics mode (proteomics / rna / metabolomics ...)
#'
#' @param wb openxlsx workbook
#' @param config full config
#' @param mode string; e.g. "proteomics", "rna"
#' @param sheet sheet name (default "Cutoffs")
add_cutoffs_sheet_legacy <- function(wb,
                                     config,
                                     mode = "proteomics",
                                     sheet = "Cutoffs") {
  stopifnot(requireNamespace("openxlsx", quietly = TRUE))
  
  de_cfg <- config$modes[[mode]]$de
  if (is.null(de_cfg)) stop("No de config found for mode: ", mode)
  
  FDR_ADJ <- isTRUE(de_cfg$use_adj_for_pass1)
  
  P_CUTOFF         <- de_cfg$p_cutoff %||% 0.05
  LINEAR_FC_CUTOFF <- de_cfg$linear_fc_cutoff %||% 1.5
  
  # recreate sheet to avoid stale named regions
  if (sheet %in% names(wb)) openxlsx::removeWorksheet(wb, sheet)
  openxlsx::addWorksheet(wb, sheetName = sheet, gridLines = TRUE)
  
  # labels (B4:B6)
  openxlsx::writeData(
    wb, sheet,
    x = c("p-value", "Adjusted pvalue (FDR)", "linear Fold Change (linearFC)"),
    startCol = 2, startRow = 4,
    colNames = FALSE, rowNames = FALSE
  )
  
  # values (C4:C6)
  openxlsx::writeData(
    wb, sheet,
    x = c(
      ifelse(FDR_ADJ, "", P_CUTOFF),
      ifelse(FDR_ADJ, P_CUTOFF, ""),
      LINEAR_FC_CUTOFF
    ),
    startCol = 3, startRow = 4,
    colNames = FALSE, rowNames = FALSE
  )
  
  # named regions (MUST match legacy)
  openxlsx::createNamedRegion(wb, sheet = sheet, cols = 3, rows = 4, name = "PVAL_CO")
  openxlsx::createNamedRegion(wb, sheet = sheet, cols = 3, rows = 5, name = "FDR_CO")
  openxlsx::createNamedRegion(wb, sheet = sheet, cols = 3, rows = 6, name = "LFC_CO")
  
  # styles
  style_used <- openxlsx::createStyle(
    border = "TopBottomLeftRight",
    borderStyle = "thick",
    fgFill = "green",
    halign = "center"
  )
  style_not_used <- openxlsx::createStyle(
    border = "TopBottomLeftRight",
    borderStyle = "thick",
    fgFill = "red",
    halign = "center"
  )
  
  if (FDR_ADJ) {
    openxlsx::addStyle(wb, sheet, style_not_used, rows = 4, cols = 3, stack = FALSE)
    openxlsx::addStyle(wb, sheet, style_used,     rows = 5:6, cols = 3, stack = FALSE)
  } else {
    openxlsx::addStyle(wb, sheet, style_not_used, rows = 5, cols = 3, stack = FALSE)
    openxlsx::addStyle(wb, sheet, style_used,     rows = c(4, 6), cols = 3, stack = FALSE)
  }
  
  openxlsx::setColWidths(wb, sheet, cols = 2, widths = "auto")
  
  invisible(TRUE)
}

#' Write legacy-style Final_results excels (generic for any mode)
#'
#' @param final_results data.frame (must include pass_any_contrast + manual_cutoffs.* optionally)
#' @param config config
#' @param out_dir output dir (run_dir/mode)
#' @param mode string, e.g. "proteomics" or "rna"
#' @param id_col column in final_results used as feature id (e.g., "Protein" / "Gene")
#' @param expr_for_de matrix used for zscore/order in DE-only excel (rows are feature IDs)
#' @param with_cutoffs logical; if TRUE add cutoffs sheet + fill manual formulas (if manual_cutoffs.* exists)
write_final_results_excels_legacy_generic <- function(final_results,
                                                      config,
                                                      out_dir,
                                                      mode,
                                                      id_col,
                                                      expr_for_de,
                                                      with_cutoffs = TRUE,
                                                      excel_order = NULL) {
  stopifnot(requireNamespace("openxlsx", quietly = TRUE))
  
  ptag <- p_tag(config)
  
  f_all <- file.path(out_dir, sprintf("Final_results_P_%s.xlsx", ptag))
  f_de  <- file.path(out_dir, sprintf("Final_results_DE_P_%s.xlsx", ptag))
  
  save_wb_results <- function(df, path, with_cutoffs = FALSE) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Results")
    openxlsx::writeData(wb, "Results", df)
    
    if (with_cutoffs) {
      add_cutoffs_sheet_legacy(wb, config, mode = mode)
      fill_manual_cutoffs_formulas_legacy(wb, "Results", df, config, mode = mode)
    }
    
    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  }
  
  # 1) ALL
  save_wb_results(final_results, f_all, with_cutoffs = isTRUE(with_cutoffs))
  
  # 2) DE-only
  is_de <- !is.na(final_results$pass_any_contrast) & final_results$pass_any_contrast == 1
  de_df <- final_results[is_de, , drop = FALSE]
  de_df <- de_df[, !startsWith(names(de_df), "manual_cutoffs"), drop = FALSE]
  
  # optional: add order + zscore columns if expr_for_de is provided
  if (!is.null(expr_for_de)) {
    expr_for_de <- as.matrix(expr_for_de)
    stopifnot(!is.null(rownames(expr_for_de)))
    
    de_ids <- intersect(de_df[[id_col]], rownames(expr_for_de))
    mat_de <- expr_for_de[de_ids, , drop = FALSE]
    
    if (nrow(mat_de) == 0) {
      save_wb_results(mat_de, f_de, with_cutoffs = FALSE)
      return(c(f_all, f_de))
    }
    
    if (!is.null(excel_order)) {
      ordered_ids <- excel_order$ordered_ids
      z <- excel_order$zscore_mat
      
      ordered_ids <- intersect(ordered_ids, de_df[[id_col]])
      order_tbl <- data.frame(
        tmp_id = ordered_ids,
        order = seq_along(ordered_ids),
        stringsAsFactors = FALSE
      )
      names(order_tbl)[1] <- id_col
      
      z_tbl <- as.data.frame(z, check.names = FALSE)
      z_tbl[[id_col]] <- rownames(z_tbl)
      
      de_df <- dplyr::as_tibble(de_df) |>
        dplyr::left_join(order_tbl, by = id_col) |>
        dplyr::left_join(z_tbl, by = id_col) |>
        dplyr::distinct()
    } else {
      z <- zscore_rows(mat_de)
      colnames(z) <- paste0(colnames(z), ".zscore")
      
      ordered_ids <- get_hclust_row_order(z, row_distance = "correlation", hclust_method = "complete")
      
      order_tbl <- data.frame(
        feature_id = ordered_ids,
        order = seq_along(ordered_ids),
        stringsAsFactors = FALSE
      )
      names(order_tbl)[1] <- id_col
      
      z_tbl <- as.data.frame(z, check.names = FALSE)
      z_tbl[[id_col]] <- rownames(z_tbl)
      
      de_df <- dplyr::as_tibble(de_df) |>
        dplyr::left_join(order_tbl, by = id_col) |>
        dplyr::left_join(z_tbl, by = id_col) |>
        dplyr::distinct()
    }
    
    
  }
  
  
  save_wb_results(de_df, f_de, with_cutoffs = FALSE)
  c(f_all, f_de)
}

#' Fill manual_cutoffs.* columns with Excel formulas (legacy)
#'
#' Writes one formula for the entire manual_cutoffs.* column using relative
#' row references and named cutoff cells: PVAL_CO / FDR_CO / LFC_CO.
#'
#' @param wb openxlsx workbook
#' @param sheet sheet name (usually "Results")
#' @param final_results data.frame written to the sheet
#' @param config full config (uses modes$proteomics$de$use_adj_for_pass1)

fill_manual_cutoffs_formulas_legacy <- function(wb, sheet, final_results, config,
                                                mode = "proteomics") {
  de_cfg <- config$modes[[mode]]$de
  use_fdr <- isTRUE(de_cfg$use_adj_for_pass1)
  
  manual_cols <- grep("^manual_cutoffs\\.", names(final_results), value = TRUE)
  if (length(manual_cols) == 0) return(invisible(NULL))
  
  start_row <- 2
  n_rows <- nrow(final_results)
  
  if (n_rows < 1) return(invisible(NULL))
  
  # Generate a sequence of row numbers (e.g., 2, 3, 4...)
  rows_seq <- start_row:(start_row + n_rows - 1)
  
  for (mcol in manual_cols) {
    contrast <- sub("^manual_cutoffs\\.", "", mcol)
    cols <- get_contrast_cols(contrast)
    
    # Require the stats columns for this contrast (skip quietly if missing)
    if (!all(c(cols$fc, cols$p, cols$padj) %in% names(final_results))) next
    
    # Excel column letters
    fc_L   <- openxlsx::int2col(match(cols$fc,   names(final_results)))
    p_L    <- openxlsx::int2col(match(cols$p,    names(final_results)))
    padj_L <- openxlsx::int2col(match(cols$padj, names(final_results)))
    m_i    <- match(mcol, names(final_results))
    
    fc_refs   <- paste0(fc_L,   rows_seq)
    p_refs    <- paste0(p_L,    rows_seq)
    padj_refs <- paste0(padj_L, rows_seq)
    
    #  Use vectorized sprintf to create a unique formula for each row
    cond <- if (use_fdr) {
      sprintf('AND(ISNUMBER(%s),%s<=FDR_CO,ABS(%s)>=LFC_CO)', 
              padj_refs, padj_refs, fc_refs)
    } else {
      sprintf('AND(ISNUMBER(%s),%s<=PVAL_CO,ABS(%s)>=LFC_CO)', 
              p_refs, p_refs, fc_refs)
    }
    
    formulas <- sprintf('IF(%s,IF(%s>0,"up","down"),"")', cond, fc_refs)
    
    # Write the vector of formulas
    openxlsx::writeFormula(
      wb, sheet,
      x        = formulas,  # <--- Now passing a vector of different formulas
      startCol = m_i,
      startRow = start_row
    )
  }
  
  invisible(TRUE)
}
# ---- Proteomics final_file ----
write_final_results_excels_proteomics <- function(final_results, pre, config, out_dir, excel_order = NULL) {
  write_final_results_excels_legacy_generic(
    final_results = final_results,
    config       = config,
    out_dir      = out_dir,
    mode         = "proteomics",
    id_col       = "Protein",
    expr_for_de  = pre$expr_imp_single,  
    with_cutoffs = TRUE,
    excel_order  = excel_order            
  )
}


# ---- RNASeq final_file ----

write_final_results_excels_legacy_rna <- function(final_results, pre, config, out_dir) {
  write_final_results_excels_legacy_generic(
    final_results = final_results,
    config       = config,
    out_dir      = out_dir,
    mode         = "rna",
    id_col       = "Gene",
    expr_for_de  = pre$expr_work,
    with_cutoffs = TRUE 
  )
}


