#' Run DESeq2 differential expression analysis
#'
#' Accepts either raw integer count matrices or tximport objects. Input type is
#' detected automatically and the appropriate DESeq2 constructor is used.
#'
#' @param counts Expression input: integer count matrix (genes x samples) OR tximport object.
#'   For tximport, must contain 'counts', 'abundance', and 'length' matrices.
#' @param meta Data frame with sample metadata
#' @param contrasts_df Data frame with columns: Contrast_name, Factor, Numerator, Denominator
#' @param de_cfg List with DE configuration (p_cutoff, linear_fc_cutoff, sample_col, sample_alignment, etc.)
#' @return List with 'dds' (DESeqDataSet) and 'tables' (list of DE result data frames)
#' @export
run_deseq2_de <- function(counts, meta, contrasts_df, de_cfg) {
    # Validate contrasts
    if (nrow(contrasts_df) == 0) {
        stop("contrasts_df is empty")
    }

    # Get factor column
    factor_col <- unique(contrasts_df$Factor)
    if (length(factor_col) != 1) {
        stop("All contrasts must use the same factor column")
    }
    factor_col <- factor_col[[1]]

    if (!factor_col %in% colnames(meta)) {
        stop(sprintf("Factor column '%s' not found in metadata", factor_col))
    }

    # Ensure factor
    meta[[factor_col]] <- as.factor(meta[[factor_col]])

    # Store levels for validation
    valid_levels <- levels(meta[[factor_col]])

    # Extract sample column from config (default to rownames-based alignment)
    sample_col <- de_cfg$sample_col %||% de_cfg$id_columns$sample_col %||% "SampleID"

    # Ensure sample_col exists in meta; if not, try to use rownames
    if (!sample_col %in% colnames(meta)) {
        if (!is.null(rownames(meta)) && !any(rownames(meta) == "")) {
            meta[[sample_col]] <- rownames(meta)
            message(sprintf("[run_deseq2_de] Using rownames as sample_col '%s'", sample_col))
        } else {
            stop(sprintf("Sample column '%s' not found in metadata and rownames are not usable", sample_col))
        }
    }

    # Determine alignment mode from config (strict by default)
    lenient_alignment <- identical(de_cfg$sample_alignment, "lenient")

    # Create DESeqDataSet using factory (handles both matrix and tximport)
    design_formula <- stats::as.formula(paste0("~ ", factor_col))
    dds0 <- create_deseq_dataset(
        expr = counts,
        meta = meta,
        design = design_formula,
        sample_col = sample_col,
        lenient_alignment = lenient_alignment
    )

    # DESeq2 mode selection from config
    deseq_mode <- de_cfg$deseq_mode %||% "default"

    if (deseq_mode == "legacy") {
        # Legacy mode: betaPrior=TRUE, modelMatrixType="expanded"
        # This matches the old "Neat RNA-Seq" pipeline behavior
        message("Using DESeq2 LEGACY mode (betaPrior=TRUE, expanded model)")
        dds <- DESeq2::DESeq(
            dds0,
            betaPrior = TRUE,
            modelMatrixType = "expanded"
        )
    } else {
        # Default mode: uses package defaults (betaPrior=FALSE in modern DESeq2)
        message("Using DESeq2 DEFAULT mode (package defaults)")
        dds <- DESeq2::DESeq(dds0)
    }

    # Extract results per contrast
    tables <- list()
    for (i in 1:nrow(contrasts_df)) {
        cn <- contrasts_df$Contrast_name[i]
        num <- contrasts_df$Numerator[i]
        den <- contrasts_df$Denominator[i]
        # FIX: Ensure contrast is a character vector, not a data.frame row
        # Using $ accessor guarantees vector extraction
        contrast <- c(as.character(contrasts_df$Factor[i]), as.character(num), as.character(den))
        # Level Validation
        if (!num %in% valid_levels) {
            warning(sprintf("Numerator '%s' not in factor levels for contrast '%s'. Skipping.", num, cn))
            next
        }
        if (!den %in% valid_levels) {
            warning(sprintf("Denominator '%s' not in factor levels for contrast '%s'. Skipping.", den, cn))
            next
        }

        # Compute alpha cutoff
        alpha_cut <- de_cfg$p_cutoff %||% de_cfg$padj_cutoff %||% 0.05

        # Extract results
        res <- DESeq2::results(
            dds,
            contrast = contrast,
            alpha = alpha_cut
        )

        # Convert to data frame with explicit column enforcement
        tab <- as.data.frame(res)

        # CRITICAL: Ensure required columns exist (even if NA)
        required_cols <- c("log2FoldChange", "pvalue", "padj")
        for (col in required_cols) {
            if (!col %in% colnames(tab)) {
                warning(sprintf("Column '%s' missing in results for contrast '%s'. Adding as NA.", col, cn))
                tab[[col]] <- NA_real_
            }
        }

        # Add gene IDs
        tab$FeatureID <- rownames(tab)

        tables[[cn]] <- tab
    }

    # Return both dds and tables for Shiny export
    return(list(
        dds = dds,
        tables = tables
    ))
}


#' Build RNA-seq summary data frame from DESeq2 results
#'
#' @param de_tables List of DESeq2 result tables (output from run_deseq2_de)
#' @param de_cfg DE configuration list
#' @return Data frame with summary statistics and pass flags
#' @export
build_rnaseq_summary_df <- function(de_tables, de_cfg) {

  if (length(de_tables) == 0) {
    warning("No DE tables provided to build_rnaseq_summary_df")
    return(data.frame())
  }
  
  padj_cutoff <- de_cfg$p_cutoff %||% de_cfg$padj_cutoff %||% 0.05
  linear_fc_cutoff <- de_cfg$linear_fc_cutoff %||% 1.5
  
  all_genes <- unique(unlist(lapply(de_tables, function(x) x$FeatureID)))
  
  summary_df <- data.frame(
    FeatureID = all_genes,
    stringsAsFactors = FALSE
  )
  
  for (cn in names(de_tables)) {
    tab <- de_tables[[cn]]
    idx <- match(summary_df$FeatureID, tab$FeatureID)
    
    
    lfc <- tab$log2FoldChange[idx]
    raw_fc <- ifelse(lfc >= 0, 2^lfc, -1 * (2^-lfc))
    rounded_fc <- signif(raw_fc, 3) 
    
    fc_col <- paste0("linearFC.", cn)
    summary_df[[fc_col]] <- rounded_fc
    
    summary_df[[paste0("pvalue.", cn)]] <- tab$pvalue[idx]
    summary_df[[paste0("padj.", cn)]] <- tab$padj[idx]
    
    pass_col <- paste0(cn, "_pass")
    
    is_sig <- !is.na(tab$padj[idx]) & 
      tab$padj[idx] <= padj_cutoff & 
      abs(as.numeric(rounded_fc)) >= linear_fc_cutoff
    
    summary_df[[pass_col]] <- ifelse(is_sig, 1, NA)
                                    
  }
  
  # 4. Pass Any 
  pass_cols <- grep("_pass$", colnames(summary_df), value = TRUE)
  if (length(pass_cols) > 0) {
    summary_df$pass_any_contrast <- apply(summary_df[, pass_cols, drop = FALSE], 1, function(x) {
      val <- any(x == 1 & !is.na(x))
      if(val) return(1) else return(NA)
    })
  } else {
    summary_df$pass_any_contrast <- NA
  }
  
  return(summary_df)
}

#' Build DE summary counts table
#'
#' Creates a summary table showing counts of upregulated, downregulated,
#' and total DE genes for each contrast, plus an overall pass_any row.
#'
#' @param summary_df Summary data frame from build_rnaseq_summary_df()
#' @param contrasts_df Contrasts data frame with Contrast_name column
#' @return Data frame with columns: Name, up, down, any
#' @export
build_de_summary_counts <- function(summary_df, contrasts_df) {
    if (is.null(summary_df) || nrow(summary_df) == 0) {
        warning("Summary data frame is empty, cannot build summary counts")
        return(data.frame(Name = character(0), up = integer(0), down = integer(0), any = integer(0)))
    }

    # Get contrast names from contrasts_df
    contrast_names <- contrasts_df$Contrast_name

    # Initialize result data frame
    result <- data.frame(
        Name = character(),
        up = integer(),
        down = integer(),
        any = integer(),
        stringsAsFactors = FALSE
    )

    # Process each contrast
    for (cn in contrast_names) {
        # Get pass column for this contrast
        pass_col <- paste0(cn, "_pass")
        fc_col <- paste0("linearFC.", cn)

        # Check if columns exist
        if (!pass_col %in% colnames(summary_df) || !fc_col %in% colnames(summary_df)) {
            warning(sprintf("Missing columns for contrast '%s', skipping", cn))
            next
        }

        # Get pass mask (only genes passing significance)
        pass_mask <- summary_df[[pass_col]] %in% TRUE

        # Count up (positive linearFC and pass)
        up_mask <- pass_mask & (summary_df[[fc_col]] > 0)
        n_up <- sum(up_mask, na.rm = TRUE)

        # Count down (negative linearFC and pass)
        down_mask <- pass_mask & (summary_df[[fc_col]] < 0)
        n_down <- sum(down_mask, na.rm = TRUE)

        # Any = up + down
        n_any <- n_up + n_down

        # Add row
        result <- rbind(result, data.frame(
            Name = cn,
            up = n_up,
            down = n_down,
            any = n_any,
            stringsAsFactors = FALSE
        ))
    }

    # Add pass_any row
    if ("pass_any_contrast" %in% colnames(summary_df)) {
        n_pass_any <- sum(summary_df$pass_any_contrast %in% TRUE, na.rm = TRUE)
        result <- rbind(result, data.frame(
            Name = "pass_any",
            up = 0L,
            down = 0L,
            any = n_pass_any,
            stringsAsFactors = FALSE
        ))
    }

    return(result)
}


#' Run limma differential expression for preprocessed RNA data
#'
#' Used when input is already normalized (preprocessed). Expects log-scale
#' expression values. Returns the same structure as run_deseq2_de() so
#' downstream code (build_rnaseq_summary_df, Shiny export, etc.) works
#' without changes.
#'
#' @param expr Numeric matrix (genes x samples), already normalized/log-transformed.
#' @param meta Data frame with sample metadata.
#' @param contrasts_df Data frame with columns: Contrast_name, Factor, Numerator, Denominator.
#' @param de_cfg DE configuration list.
#' @return List with 'dds' (NULL) and 'tables' (named list of DE result data frames).
run_limma_rna_de <- function(expr, meta, contrasts_df, de_cfg) {
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("Package 'limma' is required for DE on preprocessed data.", call. = FALSE)
    }

    if (nrow(contrasts_df) == 0) {
        stop("contrasts_df is empty")
    }

    factor_col <- unique(contrasts_df$Factor)
    if (length(factor_col) != 1) {
        stop("All contrasts must use the same factor column")
    }
    factor_col <- factor_col[[1]]

    if (!factor_col %in% colnames(meta)) {
        stop(sprintf("Factor column '%s' not found in metadata", factor_col))
    }

    # Align metadata to expression matrix columns
    sample_col <- de_cfg$sample_col %||% de_cfg$id_columns$sample_col %||% "SampleID"
    if (!sample_col %in% colnames(meta)) {
        if (!is.null(rownames(meta)) && !any(rownames(meta) == "")) {
            meta[[sample_col]] <- rownames(meta)
        } else {
            stop(sprintf("Sample column '%s' not found in metadata", sample_col))
        }
    }

    # Match metadata rows to expression columns
    meta <- meta[match(colnames(expr), meta[[sample_col]]), , drop = FALSE]

    meta[[factor_col]] <- factor(meta[[factor_col]])
    orig_levels <- levels(meta[[factor_col]])
    safe_levels <- make.names(orig_levels)
    levels(meta[[factor_col]]) <- safe_levels

    message("[run_limma_rna_de] Using limma for preprocessed data")

    design <- stats::model.matrix(stats::as.formula(paste0("~ 0 + ", factor_col)), data = meta)
    colnames(design) <- safe_levels

    # Build contrast formulas
    safe_num <- make.names(contrasts_df$Numerator)
    safe_den <- make.names(contrasts_df$Denominator)
    contrast_formulas <- setNames(
        paste(safe_num, safe_den, sep = " - "),
        contrasts_df$Contrast_name
    )

    contrast_matrix <- limma::makeContrasts(contrasts = contrast_formulas, levels = design)
    colnames(contrast_matrix) <- names(contrast_formulas)

    fit <- limma::lmFit(expr, design)
    fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix))

    # Extract per-contrast tables with DESeq2-compatible column names
    tables <- list()
    for (cn in colnames(contrast_matrix)) {
        de <- limma::topTable(fit2, coef = cn, adjust.method = "BH",
                              sort.by = "none", number = Inf)

        tab <- data.frame(
            baseMean       = de$AveExpr,
            log2FoldChange = de$logFC,
            lfcSE          = NA_real_,
            stat           = de$t,
            pvalue         = de$P.Value,
            padj           = de$adj.P.Val,
            FeatureID      = rownames(de),
            stringsAsFactors = FALSE
        )
        rownames(tab) <- rownames(de)

        tables[[cn]] <- tab
    }

    list(
        dds    = NULL,
        tables = tables
    )
}


#' Load pre-computed RNA-seq DE tables from config$modes$rna$files$de_table
#'
#' Reads CSV files and maps columns to the standard RNA DE contract
#' (same structure as run_deseq2_de / run_limma_rna_de).
#'
#' @param config Full pipeline config.
#' @param contrasts_df Optional contrasts data frame (with Contrast_name column).
#'   When provided, contrast names from this table are used instead of filenames.
#' @return List with 'dds' (NULL) and 'tables' (named list of DE data frames).
load_precomputed_rna_de <- function(config, contrasts_df = NULL) {
    cfg <- config$modes$rna

    de_files <- cfg$files$de_table
    if (is.list(de_files)) de_files <- unlist(de_files)

    # Use contrast names from contrasts_df when available (must match file count)
    if (!is.null(contrasts_df) && "Contrast_name" %in% colnames(contrasts_df) &&
        nrow(contrasts_df) == length(de_files)) {
        contrast_labels <- as.character(contrasts_df$Contrast_name)
    } else {
        # Fallback: derive from file names
        contrast_labels <- vapply(de_files, function(f) {
            bn <- tools::file_path_sans_ext(basename(f))
            sub("^de_", "", bn)
        }, character(1), USE.NAMES = FALSE)
    }

    tables <- list()
    for (i in seq_along(de_files)) {
        abs_path <- resolve_raw_path(config, de_files[i])
        if (!file.exists(abs_path)) {
            stop("Pre-computed RNA DE table not found: ", abs_path)
        }

        raw <- read_table_auto(abs_path)
        cn <- colnames(raw)

        # Feature IDs: try named columns first, then unnamed first column
        id_col <- cn[cn %in% c("FeatureID", "gene_id", "feature_id", "GeneID")][1]
        if (is.na(id_col)) {
            # Unnamed first column (readr: "...1", base R: "X", "V1")
            id_col_idx <- match(TRUE, cn %in% c("...1", "", "X", "V1"))
            feat_ids <- if (!is.na(id_col_idx)) {
                as.character(raw[[id_col_idx]])
            } else {
                rownames(raw)
            }
        } else {
            feat_ids <- as.character(raw[[id_col]])
        }

        # log2FoldChange
        lfc_col <- cn[cn %in% c("log2FoldChange", "logFC", "log2FC",
                                 "log2(FC)", "log2.FC.")][1]
        lfc_vals <- if (!is.na(lfc_col)) as.numeric(raw[[lfc_col]]) else NA_real_

        # pvalue
        pval_col <- cn[cn %in% c("pvalue", "P.Value", "PValue", "p.value",
                                  "raw.pval")][1]
        pval_vals <- if (!is.na(pval_col)) as.numeric(raw[[pval_col]]) else NA_real_

        # padj
        padj_col <- cn[cn %in% c("padj", "adj.P.Val", "FDR", "q.value",
                                  "p.adjust", "qvalue")][1]
        padj_vals <- if (!is.na(padj_col)) {
            as.numeric(raw[[padj_col]])
        } else {
            stats::p.adjust(pval_vals, method = "BH")
        }

        # baseMean
        base_col <- cn[cn %in% c("baseMean", "AveExpr", "logCPM")][1]
        base_vals <- if (!is.na(base_col)) as.numeric(raw[[base_col]]) else NA_real_

        # stat
        stat_col <- cn[cn %in% c("stat", "t", "statistic")][1]
        stat_vals <- if (!is.na(stat_col)) as.numeric(raw[[stat_col]]) else NA_real_

        tab <- data.frame(
            baseMean       = base_vals,
            log2FoldChange = lfc_vals,
            lfcSE          = NA_real_,
            stat           = stat_vals,
            pvalue         = pval_vals,
            padj           = padj_vals,
            FeatureID      = feat_ids,
            stringsAsFactors = FALSE
        )
        rownames(tab) <- feat_ids

        tables[[contrast_labels[i]]] <- tab
        message("  Loaded ", nrow(tab), " features from ", basename(de_files[i]),
                " (label: ", contrast_labels[i], ")")
    }

    message("RNA DE: loaded ", length(tables), " pre-computed contrast(s)")

    list(
        dds    = NULL,
        tables = tables
    )
}
