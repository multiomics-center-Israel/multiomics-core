#' Add legacy Cutoffs sheet + named regions (PVAL_CO/FDR_CO/LFC_CO)
#'
#' Generic version for any omics mode (proteomics / rna / metabolomics ...)
#'
#' @param wb openxlsx workbook
#' @param config full config
#' @param mode string; e.g. "proteomics", "rna"
#' @param sheet sheet name (default "Cutoffs")
add_cutoffs_sheet_legacy <- function(wb, config, mode = "proteomics", sheet = "Cutoffs") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Package 'openxlsx' is required.")

    de_cfg <- config$modes[[mode]]$de
    if (is.null(de_cfg)) stop("No de config found for mode: ", mode)

    FDR_ADJ <- isTRUE(de_cfg$use_adj_for_pass1)
    P_CUTOFF <- de_cfg$p_cutoff %||% 0.05
    LINEAR_FC_CUTOFF <- de_cfg$linear_fc_cutoff %||% 1.5

    if (sheet %in% names(wb)) openxlsx::removeWorksheet(wb, sheet)
    openxlsx::addWorksheet(wb, sheetName = sheet, gridLines = TRUE)

    openxlsx::writeData(wb, sheet, x = c("p-value", "Adjusted pvalue (FDR)", "linear Fold Change (linearFC)"), startCol = 2, startRow = 4, colNames = FALSE, rowNames = FALSE)
    openxlsx::writeData(wb, sheet, x = c(ifelse(FDR_ADJ, "", P_CUTOFF), ifelse(FDR_ADJ, P_CUTOFF, ""), LINEAR_FC_CUTOFF), startCol = 3, startRow = 4, colNames = FALSE, rowNames = FALSE)

    openxlsx::createNamedRegion(wb, sheet = sheet, cols = 3, rows = 4, name = "PVAL_CO")
    openxlsx::createNamedRegion(wb, sheet = sheet, cols = 3, rows = 5, name = "FDR_CO")
    openxlsx::createNamedRegion(wb, sheet = sheet, cols = 3, rows = 6, name = "LFC_CO")

    style_used <- openxlsx::createStyle(border = "TopBottomLeftRight", borderStyle = "thick", fgFill = "green", halign = "center")
    style_not_used <- openxlsx::createStyle(border = "TopBottomLeftRight", borderStyle = "thick", fgFill = "red", halign = "center")

    if (FDR_ADJ) {
        openxlsx::addStyle(wb, sheet, style_not_used, rows = 4, cols = 3, stack = FALSE)
        openxlsx::addStyle(wb, sheet, style_used, rows = 5:6, cols = 3, stack = FALSE)
    } else {
        openxlsx::addStyle(wb, sheet, style_not_used, rows = 5, cols = 3, stack = FALSE)
        openxlsx::addStyle(wb, sheet, style_used, rows = c(4, 6), cols = 3, stack = FALSE)
    }
    openxlsx::setColWidths(wb, sheet, cols = 2, widths = "auto")
    invisible(TRUE)
}

#' Get p-value cutoff tag for filename (generic for any mode)
p_tag_generic <- function(config, mode, default = "NA") {
    p <- config$modes[[mode]]$de$p_cutoff
    if (is.null(p) || is.na(p)) {
        return(default)
    }
    format(p, trim = TRUE, scientific = FALSE)
}

#' Write legacy-style Final_results excels (generic for any mode)
write_final_results_excels_legacy_generic <- function(final_results, config, out_dir, mode, id_col, expr_for_de, with_cutoffs = TRUE, clustering_res = NULL) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Package 'openxlsx' is required.")
    # Validate inputs
    if (is.null(final_results)) {
        stop("final_results is NULL. Cannot export to Excel. Check that DE analysis produced results.")
    }
    if (!is.data.frame(final_results)) {
        stop("final_results must be a data.frame, got: ", class(final_results)[1])
    }

    # Get p-value cutoff tag for filename
    p_tag_val <- p_tag_generic(config, mode)

    # Renamed: Final_results_P -> Final_results_ALL
    f_all <- file.path(out_dir, sprintf("Final_results_ALL_P_%s.xlsx", p_tag_val))
    f_de <- file.path(out_dir, sprintf("Final_results_DE_P_%s.xlsx", p_tag_val))

    save_wb_results <- function(df, path, with_cutoffs = FALSE) {
        wb <- openxlsx::createWorkbook()
        openxlsx::addWorksheet(wb, "Results")
        openxlsx::writeData(wb, "Results", df)

        if (with_cutoffs) {
            add_cutoffs_sheet_legacy(wb, config, mode = mode)
            fill_manual_cutoffs_formulas_legacy(wb, "Results", df, config, mode = mode)
        }
        openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
        path # Return the file path for targets tracking
    }

    # Save full results (ALL table) and capture path
    f_all_created <- save_wb_results(final_results, f_all, with_cutoffs = isTRUE(with_cutoffs))

    # Build DE table with proper column ordering
    is_de <- !is.na(final_results$pass_any_contrast) & final_results$pass_any_contrast == 1
    de_df <- final_results[is_de, , drop = FALSE]

    # Remove manual_cutoffs and pass_any_contrast columns
    de_df <- de_df[, !startsWith(names(de_df), "manual_cutoffs") & names(de_df) != "pass_any_contrast", drop = FALSE]

    if (!is.null(expr_for_de)) {
        expr_for_de <- as.matrix(expr_for_de)
        de_ids <- intersect(de_df[[id_col]], rownames(expr_for_de))
        mat_de <- expr_for_de[de_ids, , drop = FALSE]

        if (nrow(mat_de) == 0) {
            f_de_created <- save_wb_results(de_df, f_de, with_cutoffs = FALSE)
            return(c(f_all_created, f_de_created))
        }

        # Integrate clustering order and z-scores if available
        cl_obj <- clustering_res %||% list()
        excel_ord <- cl_obj$excel_order %||% NULL
        
        
        if (!is.null(excel_ord) && !is.null(excel_ord$ordered_ids)) {
            ordered_ids <- excel_ord$ordered_ids

            # 1. Add 'order' column (hierarchical clustering rank)
            ranks <- match(de_df[[id_col]], ordered_ids)
            de_df$order <- ranks

            # 2. Add Z-score columns
            if (!is.null(excel_ord$zscore_mat)) {
                zmat <- excel_ord$zscore_mat
                idx <- match(de_df[[id_col]], rownames(zmat))
                z_sub <- zmat[idx, , drop = FALSE]
                de_df <- cbind(de_df, z_sub)
            }
        } else {
            # No clustering: add order column with NA
            de_df$order <- NA_integer_
        }
        
        de_df <- add_default_order_if_missing(de_df, mat_de, id_col)
        de_df <- add_zscores_if_missing(de_df, mat_de, id_col)

        # Reorder columns to: ID, raw_counts, DE_stats, order, z-scores
        # Identify column groups
        id_cols <- id_col
        expr_cols <- colnames(mat_de)
        de_stat_cols <- grep("^(linearFC|pvalue|padj|upDown)\\.", names(de_df), value = TRUE)
        order_col <- "order"
        zscore_cols <- grep("\\.zscore$", names(de_df), value = TRUE)

        # Check which expression columns are already present (from build_final_results_generic)
        expr_cols_present <- intersect(expr_cols, names(de_df))

        # Build desired column order
        desired_order <- c(
            id_cols,
            expr_cols_present,  # Raw counts (already in de_df)
            de_stat_cols,
            order_col,
            zscore_cols
        )

        # Reorder columns
        de_df <- de_df[, desired_order, drop = FALSE]
    }

    # Save DE results and capture path
    f_de_created <- save_wb_results(de_df, f_de, with_cutoffs = FALSE)
    c(f_all_created, f_de_created)
}

#' Get standard column names for a contrast
#' @param contrast Contrast name
#' @param mode "proteomics" (uses .imputs.) or "rna" (no .imputs.)
get_contrast_cols <- function(contrast, mode = "proteomics") {
    stopifnot(is.character(contrast), length(contrast) == 1, nzchar(contrast))

    # Proteomics DE summary strips spaces from contrast names;
    # RNA-seq keeps them as-is.  Only normalize for proteomics.
    if (mode != "rna") {
        contrast <- normalize_contrast_name(contrast)
    }

    # RNA doesn't use ".imputs." in column names
    if (mode == "rna") {
        list(
            fc     = paste0("linearFC.", contrast),
            p      = paste0("pvalue.", contrast),
            padj   = paste0("padj.", contrast),
            pass   = paste0(contrast, "_pass"),
            updown = paste0("upDown.", contrast),
            manual = paste0("manual_cutoffs.", contrast)
        )
    } else {
        # Proteomics (uses imputation naming)
        list(
            fc     = paste0("linearFC.imputs.", contrast_safe),
            p      = paste0("pvalue.imputs.", contrast_safe),
            padj   = paste0("padj.imputs.", contrast_safe),
            pass   = paste0("pass.imputs.", contrast_safe),
            updown = paste0("upDown.imputs.", contrast_safe),
            manual = paste0("manual_cutoffs.", contrast_safe)
        )
    }
}

fill_manual_cutoffs_formulas_legacy <- function(wb, sheet, final_results, config, mode = "proteomics") {
    de_cfg <- config$modes[[mode]]$de
    use_fdr <- isTRUE(de_cfg$use_adj_for_pass1)

    manual_cols <- grep("^manual_cutoffs\\.", names(final_results), value = TRUE)
    if (length(manual_cols) == 0) {
        return(invisible(NULL))
    }

    start_row <- 2
    n_rows <- nrow(final_results)
    if (n_rows < 1) {
        return(invisible(NULL))
    }

    rows_seq <- start_row:(start_row + n_rows - 1)

    for (mcol in manual_cols) {
        contrast <- sub("^manual_cutoffs\\.", "", mcol)
        cols <- get_contrast_cols(contrast, mode = mode)
        if (!all(c(cols$fc, cols$p, cols$padj) %in% names(final_results))) next

        fc_L <- openxlsx::int2col(match(cols$fc, names(final_results)))
        p_L <- openxlsx::int2col(match(cols$p, names(final_results)))
        padj_L <- openxlsx::int2col(match(cols$padj, names(final_results)))
        m_i <- match(mcol, names(final_results))

        fc_refs <- paste0(fc_L, rows_seq)
        p_refs <- paste0(p_L, rows_seq)
        padj_refs <- paste0(padj_L, rows_seq)

        cond <- if (use_fdr) {
            sprintf("AND(ISNUMBER(%s),%s<=FDR_CO,ABS(%s)>=LFC_CO)", padj_refs, padj_refs, fc_refs)
        } else {
            sprintf("AND(ISNUMBER(%s),%s<=PVAL_CO,ABS(%s)>=LFC_CO)", p_refs, p_refs, fc_refs)
        }
        formulas <- sprintf('IF(%s,IF(%s>0,"up","down"),"")', cond, fc_refs)

        openxlsx::writeFormula(wb, sheet, x = formulas, startCol = m_i, startRow = start_row)
    }
    invisible(TRUE)
}

#' Build final results table (generic for any mode)
#'
#' SEMANTICS (CRITICAL):
#' - pass.imputs.<contrast> is {NA, 1}: NA = did not pass, 1 = passed
#' - pass_any_contrast already exists in summary_df (computed by add_pass_any_contrast)
#' - upDown.imputs.<contrast> should ONLY be populated when pass == 1
#' - FC columns are SIGNED (linearFC or logFC), not ratios
#'
#' @param summary_df DE summary with FeatureID, contrast columns, pass_any_contrast
#' @param expr_df Expression matrix (features x samples)
#' @param contrasts_df Contrasts table with Contrast_name column
#' @param feature_id_col Name of ID column in output ("Protein", "Gene", etc.)
#' @param annot_cols Named vector: c(output_col = "source_col_in_summary_or_rowdata")
#' @param row_data Optional annotation data.frame
#' @param fc_is_signed Logical; if TRUE (default), FC is signed (linearFC/logFC).
#'                     If FALSE, must provide fc_direction_col.
#' @param fc_direction_col Optional; column name for direction if FC is unsigned ratio
#'
#' @return data.frame with ID, annotations, expression, DE stats, pass_any_contrast
build_final_results_generic <- function(
  summary_df,
  expr_df,
  contrasts_df,
  feature_id_col = "FeatureID",
  annot_cols = NULL,
  row_data = NULL,
  fc_is_signed = TRUE,
  fc_direction_col = NULL,
  mode = "proteomics"  # FIX 2: Add mode parameter for column naming
) {
    # ============================================================
    # VALIDATION (explicit errors, not stopifnot)
    # ============================================================

    if (!is.data.frame(summary_df)) {
        stop("summary_df must be a data.frame, got: ", class(summary_df)[1])
    }

    if (!(feature_id_col %in% colnames(summary_df))) {
        stop(
            "summary_df must have a '", feature_id_col, "' column. Found columns: ",
            paste(head(colnames(summary_df), 10), collapse = ", ")
        )
    }

    if (!is.data.frame(contrasts_df)) {
        stop("contrasts_df must be a data.frame, got: ", class(contrasts_df)[1])
    }

    if (!("Contrast_name" %in% colnames(contrasts_df))) {
        stop(
            "contrasts_df must have a 'Contrast_name' column. Found columns: ",
            paste(colnames(contrasts_df), collapse = ", ")
        )
    }

    if (nrow(contrasts_df) == 0) {
        stop("contrasts_df is empty. At least one contrast is required.")
    }

    # FC direction validation
    if (!isTRUE(fc_is_signed) && is.null(fc_direction_col)) {
        stop("If fc_is_signed = FALSE, must provide fc_direction_col for determining up/down direction")
    }

    # ============================================================
    # BUILD BASE WITH ID
    # ============================================================

    base <- data.frame(
        ID = summary_df[[feature_id_col]],
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    names(base)[1] <- feature_id_col

    # ============================================================
    # ADD ANNOTATIONS WITH MATCH-RATE CHECK
    # ============================================================

    if (!is.null(annot_cols) && length(annot_cols) > 0) {
        for (out_col in names(annot_cols)) {
            src_col <- annot_cols[out_col]

            # Try summary_df first
            if (src_col %in% colnames(summary_df)) {
                base[[out_col]] <- summary_df[[src_col]]
            }
            # Then try row_data
            else if (!is.null(row_data) && src_col %in% colnames(row_data)) {
                if (is.null(rownames(row_data))) {
                    warning(sprintf(
                        "row_data has no rownames. Cannot match annotation column '%s'. Setting to NA.",
                        src_col
                    ))
                    base[[out_col]] <- NA
                } else {
                    matched_vals <- row_data[[src_col]][match(base[[feature_id_col]], rownames(row_data))]

                    # Check match rate
                    match_rate <- sum(!is.na(matched_vals)) / length(matched_vals)
                    if (match_rate < 0.95) {
                        warning(sprintf(
                            "Low match rate for annotation '%s': %.1f%% (%d/%d features matched). Check that row_data rownames match feature IDs.",
                            out_col, match_rate * 100, sum(!is.na(matched_vals)), length(matched_vals)
                        ))
                    }

                    base[[out_col]] <- matched_vals
                }
            }
            # Not found anywhere
            else {
                warning(sprintf(
                    "Annotation column '%s' not found in summary_df or row_data. Setting to NA.",
                    src_col
                ))
                base[[out_col]] <- NA
            }
        }
    }

    # ============================================================
    # ADD EXPRESSION VALUES
    # ============================================================

    expr_df <- as.data.frame(expr_df, check.names = FALSE)

    if (is.null(rownames(expr_df))) {
        warning("expr_df has no rownames. Cannot add expression values.")
    } else {
        expr_matched <- expr_df[match(base[[feature_id_col]], rownames(expr_df)), , drop = FALSE]

        # Check match rate
        match_rate <- sum(base[[feature_id_col]] %in% rownames(expr_df)) / nrow(base)
        if (match_rate < 0.95) {
            warning(sprintf(
                "Low match rate for expression values: %.1f%% (%d/%d features matched). Check that expr_df rownames match feature IDs.",
                match_rate * 100, sum(base[[feature_id_col]] %in% rownames(expr_df)), nrow(base)
            ))
        }

        base <- cbind(base, expr_matched)
    }

    # ============================================================
    # ADD CONTRAST STATISTICS
    # ============================================================

    contrast_names <- contrasts_df$Contrast_name
    m <- match(base[[feature_id_col]], summary_df[[feature_id_col]])

    for (cn in contrast_names) {
        cols <- get_contrast_cols(cn, mode = mode)  # FIX 2: Pass mode parameter

        # Validate required columns exist
        needed <- c(cols$fc, cols$p, cols$padj)
        missing <- setdiff(needed, colnames(summary_df))
        if (length(missing) > 0) {
            stop(sprintf(
                "Missing required columns for contrast '%s': %s\nAvailable columns: %s",
                cn,
                paste(missing, collapse = ", "),
                paste(head(colnames(summary_df), 20), collapse = ", ")
            ))
        }

        # Extract values
        fc_vals <- summary_df[[cols$fc]][m]
        pass_vals <- if (cols$pass %in% colnames(summary_df)) {
            summary_df[[cols$pass]][m]
        } else {
            rep(NA, length(m))
        }

        # Populate FC, p-value, adjusted p-value
        base[[cols$fc]] <- fc_vals
        base[[cols$p]] <- summary_df[[cols$p]][m]
        base[[cols$padj]] <- summary_df[[cols$padj]][m]

        # ============================================================
        # DETERMINE UP/DOWN DIRECTION
        # ============================================================

        if (isTRUE(fc_is_signed)) {
            # FC is signed (linearFC or logFC) - use sign directly
            # Validate that FC looks signed (has negative values)
            if (all(fc_vals >= 0, na.rm = TRUE) && any(!is.na(fc_vals))) {
                warning(sprintf(
                    "Contrast '%s': FC column '%s' has no negative values. Are you sure fc_is_signed = TRUE?",
                    cn, cols$fc
                ))
            }

            direction <- ifelse(
                !is.na(pass_vals) & pass_vals == 1,
                ifelse(as.numeric(fc_vals) >= 0, "up", "down"),
                ""
            )
        } else {
            # FC is unsigned ratio - use direction column
            if (is.null(fc_direction_col) || !(fc_direction_col %in% colnames(summary_df))) {
                stop(sprintf(
                    "fc_direction_col '%s' not found in summary_df for contrast '%s'",
                    fc_direction_col, cn
                ))
            }

            dir_vals <- summary_df[[fc_direction_col]][m]
            direction <- ifelse(
                !is.na(pass_vals) & pass_vals == 1,
                dir_vals,
                ""
            )
        }

        base[[cols$updown]] <- direction
        base[[cols$manual]] <- NA
    }

    # ============================================================
    # USE EXISTING pass_any_contrast
    # ============================================================

    if ("pass_any_contrast" %in% colnames(summary_df)) {
        base$pass_any_contrast <- summary_df$pass_any_contrast[m]
    } else {
        # Fallback: compute if missing (shouldn't happen in normal pipeline)
        warning("pass_any_contrast not found in summary_df. Computing from pass columns.")
        pass_cols <- paste0("pass.imputs.", contrast_names)
        existing_pass_cols <- intersect(pass_cols, colnames(summary_df))

        if (length(existing_pass_cols) > 0) {
            pass_mat <- summary_df[m, existing_pass_cols, drop = FALSE]
            # Count features where pass == 1 (not just non-NA)
            n_pass <- rowSums(!is.na(pass_mat) & pass_mat == 1, na.rm = TRUE)
            base$pass_any_contrast <- ifelse(n_pass > 0, 1, NA)
        } else {
            base$pass_any_contrast <- NA
        }
    }

    base
}


#' Safe Z-score calculation for DE tables
add_zscores_if_missing <- function(df, expr_mat, id_col) {
  # If Z-scores already exist, just return the data frame
  if (any(grepl("\\.zscore$", names(df)))) return(df)
  
  # Otherwise, calculate and bind
  row_means <- rowMeans(expr_mat, na.rm = TRUE)
  row_sds   <- apply(expr_mat, 1, sd, na.rm = TRUE)
  z_mat     <- (expr_mat - row_means) / (row_sds + 1e-6) # Avoid div by zero
  colnames(z_mat) <- paste0(colnames(expr_mat), ".zscore")
  
  # Match and join
  idx <- match(df[[id_col]], rownames(z_mat))
  return(cbind(df, z_mat[idx, , drop = FALSE]))
}

#' Fallback hierarchical clustering for row order
add_default_order_if_missing <- function(df, expr_mat, id_col) {
  # If order is already populated (not all NA), keep it
  if ("order" %in% names(df) && !all(is.na(df$order))) {
    return(df)
  }
  
  # Calculate a simple hierarchical clustering order
  # We use a basic Euclidean/Complete linkage
  message("Calculating fallback hierarchical clustering order...")
  
  # Handle NAs by replacing with row mean (safe for proteomics)
  mat_clean <- expr_mat
  if (any(is.na(mat_clean))) {
    row_means <- rowMeans(mat_clean, na.rm = TRUE)
    idx_na <- which(is.na(mat_clean), arr.ind = TRUE)
    mat_clean[idx_na] <- row_means[idx_na[,1]]
  }
  
  # Simple clustering
  dists <- dist(mat_clean)
  hc <- hclust(dists, method = "complete")
  
  # Create an order mapping
  ordered_ids <- rownames(mat_clean)[hc$order]
  ranks <- match(df[[id_col]], ordered_ids)
  
  df$order <- ranks
  return(df)
}
