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

#' Column glossary + fold-change reconciliation notes for a mode
#'
#' Returns the text shown on the "How to read" sheet: what each column block
#' holds, and how to get from the per-sample values to the reported fold change.
#' Split out from \code{\link{add_provenance_sheet}} so the wording is testable
#' without openxlsx.
#'
#' @param mode "rna", "proteomics", or any other mode (generic fallback).
#' @return List with \code{glossary} (data.frame: Column, Meaning) and
#'   \code{notes} (character vector, one paragraph per element).
build_provenance_notes <- function(mode = "rna") {
    fc_rule <- paste(
        "linearFC = 2^log2FC when log2FC is at least 0, and -1 / 2^log2FC when it is negative.",
        "It is a signed linear fold change, not a log: linearFC = -1.61 means 1.61-fold lower in the numerator group."
    )

    if (identical(mode, "rna")) {
        glossary <- data.frame(
            Column = c(
                "<sample>",
                "<sample>.norm",
                "Mean.<group>",
                "CV.<group>",
                "log2FC.<contrast>",
                "linearFC.<contrast>",
                "pvalue.<contrast>",
                "padj.<contrast>",
                "upDown.<contrast>"
            ),
            Meaning = c(
                "Raw counts, after gene filtering. Not corrected for sequencing depth.",
                "DESeq2 normalized counts: the raw count divided by that sample's size factor. This is the matrix the DE model was fitted on.",
                "Arithmetic mean of the .norm values across the replicates of that group.",
                "Coefficient of variation (%) within the group, computed on linear CPM.",
                "DESeq2 log2 fold change for the contrast (numerator relative to denominator).",
                fc_rule,
                "Wald test p-value.",
                "Benjamini-Hochberg adjusted p-value.",
                "Direction, filled only for features that pass both cutoffs (see the Cutoffs sheet)."
            ),
            stringsAsFactors = FALSE
        )
        notes <- c(
            "Reconciling the fold change from this row:",
            "1. log2( Mean.<numerator> / Mean.<denominator> ) should land close to log2FC.",
            "2. Apply the linearFC rule above to log2FC. This step is exact.",
            "The same arithmetic on the raw <sample> columns will not reproduce log2FC. Raw counts are not corrected for library size, so any difference in sequencing depth between the groups shifts the ratio.",
            "Step 1 is an approximation, not an identity. log2FC is a negative-binomial GLM coefficient, not a ratio of arithmetic means; the two agree closely for well-expressed features and can differ for low-count ones."
        )
    } else if (identical(mode, "proteomics")) {
        glossary <- data.frame(
            Column = c(
                "<sample>",
                "<sample>.norm",
                "Mean.<group>",
                "CV.<group>",
                "log2FC.imputs.<contrast>",
                "linearFC.imputs.<contrast>",
                "pvalue.imputs.<contrast>",
                "padj.imputs.<contrast>",
                "upDown.imputs.<contrast>"
            ),
            Meaning = c(
                "log2 intensity after filtering and normalization, before imputation. Blank cells were not measured.",
                "The same matrix after imputation, for one representative imputation run. This is the kind of matrix limma was fitted on.",
                "Arithmetic mean of the .norm log2 values across the replicates of that group.",
                "Coefficient of variation (%) within the group, on linear intensities back-transformed from the unimputed values, so measured values only.",
                "log2 of the mean linear ratio across the imputation runs: log2( mean of 2^logFC over runs ).",
                fc_rule,
                "Quantile across imputation runs of the moderated t-test p-value.",
                "Quantile across imputation runs of the Benjamini-Hochberg adjusted p-value.",
                "Direction, filled only for features that pass the consensus rule (see the Cutoffs sheet)."
            ),
            stringsAsFactors = FALSE
        )
        notes <- c(
            "Reconciling the fold change from this row:",
            "1. Mean.<numerator> minus Mean.<denominator> (a difference, because the values are log2) should land close to log2FC.imputs.",
            "2. Apply the linearFC rule above to log2FC.imputs. This step is exact.",
            "Step 1 is an approximation for two reasons: the reported statistic averages over all imputation runs while the .norm block shows a single run, and limma reports a moderated model coefficient rather than a difference of group means.",
            "The unimputed <sample> columns will differ again, because features with missing values contribute to the model only after imputation."
        )
    } else {
        glossary <- data.frame(
            Column = c("<sample>", "Mean.<group>", "CV.<group>",
                       "log2FC.<contrast>", "linearFC.<contrast>"),
            Meaning = c(
                "Per-sample values as exported by this mode.",
                "Arithmetic mean across the replicates of that group.",
                "Coefficient of variation (%) within the group.",
                "log2 fold change for the contrast, as reported by the DE model.",
                fc_rule
            ),
            stringsAsFactors = FALSE
        )
        notes <- "Group means are descriptive summaries. The reported fold change is a model estimate, so the two agree closely but are not identical."
    }

    list(glossary = glossary, notes = notes)
}

#' Add a "How to read" sheet explaining the column blocks and fold changes
#'
#' Written to every Final_results workbook so the numbers in the Results sheet
#' can be rechecked without reading the pipeline source.
#'
#' @param wb openxlsx workbook.
#' @param mode Mode string, e.g. "rna" or "proteomics".
#' @param sheet Sheet name (default "How to read").
#' @return TRUE, invisibly.
add_provenance_sheet <- function(wb, mode = "rna", sheet = "How to read") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Package 'openxlsx' is required.")

    info <- build_provenance_notes(mode)

    if (sheet %in% names(wb)) openxlsx::removeWorksheet(wb, sheet)
    openxlsx::addWorksheet(wb, sheetName = sheet, gridLines = FALSE)

    title_style   <- openxlsx::createStyle(textDecoration = "bold", fontSize = 12)
    header_style  <- openxlsx::createStyle(textDecoration = "bold",
                                           border = "bottom", borderStyle = "thin")
    wrap_style    <- openxlsx::createStyle(wrapText = TRUE, valign = "top")

    openxlsx::writeData(wb, sheet, x = "Columns in the Results sheet",
                        startCol = 1, startRow = 1, colNames = FALSE)
    openxlsx::addStyle(wb, sheet, title_style, rows = 1, cols = 1, stack = TRUE)

    glossary_start <- 3
    openxlsx::writeData(wb, sheet, info$glossary, startRow = glossary_start)
    openxlsx::addStyle(wb, sheet, header_style, rows = glossary_start, cols = 1:2, stack = TRUE)
    openxlsx::addStyle(wb, sheet, wrap_style,
                       rows = (glossary_start + 1):(glossary_start + nrow(info$glossary)),
                       cols = 2, stack = TRUE)

    notes_start <- glossary_start + nrow(info$glossary) + 2
    openxlsx::writeData(wb, sheet, x = "Checking the numbers",
                        startCol = 1, startRow = notes_start, colNames = FALSE)
    openxlsx::addStyle(wb, sheet, title_style, rows = notes_start, cols = 1, stack = TRUE)
    openxlsx::writeData(wb, sheet, x = info$notes,
                        startCol = 1, startRow = notes_start + 1, colNames = FALSE)
    openxlsx::addStyle(wb, sheet, wrap_style,
                       rows = (notes_start + 1):(notes_start + length(info$notes)),
                       cols = 1, stack = TRUE)

    openxlsx::setColWidths(wb, sheet, cols = 1, widths = 34)
    openxlsx::setColWidths(wb, sheet, cols = 2, widths = 95)
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
#'
#' @param sample_meta Optional sample metadata data.frame for annotation rows.
#' @param sample_id_col Column name in sample_meta that matches expression column names.
#' @param annotation_rows Character vector of metadata column names to show as
#'   annotation rows above the data. NULL = all non-ID columns.
#' @param sample_label_cols Character vector of metadata columns to concatenate
#'   (with "_") for informative sample headers. NULL = keep original IDs.
#' @param provenance_sheet Logical; add a "How to read" sheet documenting the
#'   column blocks and how to reconcile log2FC/linearFC with the per-sample
#'   values. Opt-in so modes that do not export the normalized block are
#'   unaffected.
write_final_results_excels_legacy_generic <- function(final_results, config, out_dir, mode, id_col,
                                                       expr_for_de, with_cutoffs = TRUE,
                                                       clustering_res = NULL,
                                                       sample_meta = NULL,
                                                       sample_id_col = NULL,
                                                       annotation_rows = NULL,
                                                       sample_label_cols = NULL,
                                                       provenance_sheet = FALSE) {
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

    # ---- Build sample label mapping and annotation data ----
    sample_label_map <- NULL
    annot_row_data   <- NULL
    has_sample_annots <- !is.null(sample_meta) && !is.null(sample_id_col) &&
        sample_id_col %in% colnames(sample_meta)

    if (has_sample_annots) {
        meta_ids <- as.character(sample_meta[[sample_id_col]])

        # Build informative sample labels if configured
        if (!is.null(sample_label_cols) && length(sample_label_cols) > 0) {
            valid_cols <- intersect(sample_label_cols, colnames(sample_meta))
            if (length(valid_cols) > 0) {
                labels <- apply(sample_meta[, valid_cols, drop = FALSE], 1, function(row) {
                    paste(as.character(row), collapse = "_")
                })
                # Ensure uniqueness by appending index for duplicates
                labels <- make.unique(labels, sep = "_")
                sample_label_map <- stats::setNames(labels, meta_ids)
            }
        }

        # Determine which metadata columns to use as annotation rows
        if (is.null(annotation_rows)) {
            annotation_rows <- setdiff(colnames(sample_meta), sample_id_col)
        } else {
            annotation_rows <- intersect(annotation_rows, colnames(sample_meta))
        }

        if (length(annotation_rows) > 0) {
            # Build a named list: annotation_label -> named vector (sample_id -> value)
            annot_row_data <- lapply(annotation_rows, function(col) {
                stats::setNames(as.character(sample_meta[[col]]), meta_ids)
            })
            names(annot_row_data) <- annotation_rows
        }
    }

    save_wb_results <- function(df, path, with_cutoffs = FALSE) {
        wb <- openxlsx::createWorkbook()
        sheet <- "Results"
        openxlsx::addWorksheet(wb, sheet, gridLines = TRUE)

        # Identify sample columns in df (expression + zscore columns that match meta IDs)
        df_cols <- colnames(df)
        sample_cols_in_df <- character(0)
        if (has_sample_annots) {
            meta_ids <- as.character(sample_meta[[sample_id_col]])
            sample_cols_in_df <- intersect(df_cols, meta_ids)
        }

        # ---- Apply sample label renaming ----
        df_out <- df
        if (!is.null(sample_label_map) && length(sample_cols_in_df) > 0) {
            for (sc in sample_cols_in_df) {
                if (sc %in% names(sample_label_map)) {
                    new_name <- sample_label_map[[sc]]
                    colnames(df_out)[colnames(df_out) == sc] <- new_name
                }
            }
            # Also rename the derived per-sample blocks (normalized values,
            # z-scores) so a renamed header stays readable across the whole row
            for (sfx in c(".norm", ".zscore")) {
                sfx_re <- paste0(gsub(".", "\\.", sfx, fixed = TRUE), "$")
                derived_cols <- grep(sfx_re, colnames(df_out), value = TRUE)
                for (dc in derived_cols) {
                    base_id <- sub(sfx_re, "", dc)
                    if (base_id %in% names(sample_label_map)) {
                        new_name <- paste0(sample_label_map[[base_id]], sfx)
                        colnames(df_out)[colnames(df_out) == dc] <- new_name
                    }
                }
            }
        }

        # ---- Detect DE stat columns and group by contrast ----
        de_col_pattern <- "^(log2FC|linearFC|pvalue|padj|upDown)\\."
        de_col_indices <- grep(de_col_pattern, colnames(df_out))
        contrast_groups <- list()
        if (length(de_col_indices) > 0) {
            for (ci in de_col_indices) {
                col_name <- colnames(df_out)[ci]
                contrast <- sub("^[^.]+\\.", "", col_name)
                contrast <- sub("^imputs\\.", "", contrast)
                contrast_groups[[contrast]] <- c(contrast_groups[[contrast]], ci)
            }
        }
        has_contrast_groups <- length(contrast_groups) > 0
        n_contrast_header_rows <- if (has_contrast_groups) 1L else 0L

        # ---- Number of annotation rows above the data ----
        n_annot_rows <- 0
        if (!is.null(annot_row_data) && length(annot_row_data) > 0 &&
            length(sample_cols_in_df) > 0) {
            # Always include "Sample_ID" row for traceability
            n_annot_rows <- length(annot_row_data) + 1  # +1 for Sample_ID row
        }

        # ---- Write annotation rows (rows 1..n_annot_rows) ----
        if (n_annot_rows > 0) {
            # Map every per-sample column to its sample ID and its column index
            # in df_out (after potential renaming). Both the measured block and
            # the normalized block belong to the same sample, so the group
            # labels above must repeat over both — otherwise the second block
            # sits under a blank header and reads as unrelated.
            sample_col_map <- do.call(rbind, lapply(sample_cols_in_df, function(sc) {
                display <- if (!is.null(sample_label_map) && sc %in% names(sample_label_map)) {
                    sample_label_map[[sc]]
                } else {
                    sc
                }
                idx <- match(c(display, paste0(display, ".norm")), colnames(df_out))
                data.frame(sample_id = sc, col_idx = idx,
                           stringsAsFactors = FALSE)
            }))
            sample_col_map <- sample_col_map[!is.na(sample_col_map$col_idx), , drop = FALSE]

            row_idx <- 1
            # Write Sample_ID row first
            openxlsx::writeData(wb, sheet, x = "Sample_ID", startCol = 1, startRow = row_idx,
                                colNames = FALSE, rowNames = FALSE)
            for (j in seq_len(nrow(sample_col_map))) {
                openxlsx::writeData(wb, sheet, x = sample_col_map$sample_id[j],
                                    startCol = sample_col_map$col_idx[j],
                                    startRow = row_idx,
                                    colNames = FALSE, rowNames = FALSE)
            }
            row_idx <- row_idx + 1

            # Write metadata annotation rows
            for (annot_name in names(annot_row_data)) {
                vals <- annot_row_data[[annot_name]]
                openxlsx::writeData(wb, sheet, x = annot_name, startCol = 1,
                                    startRow = row_idx, colNames = FALSE, rowNames = FALSE)
                for (j in seq_len(nrow(sample_col_map))) {
                    sc <- sample_col_map$sample_id[j]
                    val <- if (sc %in% names(vals)) vals[[sc]] else NA
                    openxlsx::writeData(wb, sheet, x = val,
                                        startCol = sample_col_map$col_idx[j],
                                        startRow = row_idx,
                                        colNames = FALSE, rowNames = FALSE)
                }
                row_idx <- row_idx + 1
            }

            # ---- Style annotation rows ----
            annot_fill <- openxlsx::createStyle(fgFill = "#F0F0F0", textDecoration = "bold")
            label_style <- openxlsx::createStyle(fgFill = "#D9D9D9", textDecoration = "bold")
            separator_style <- openxlsx::createStyle(
                border = "bottom", borderStyle = "thick", borderColour = "#333333"
            )

            for (r in seq_len(n_annot_rows)) {
                openxlsx::addStyle(wb, sheet, annot_fill, rows = r,
                                   cols = seq_along(df_cols), stack = TRUE)
                # Bold label in column 1
                openxlsx::addStyle(wb, sheet, label_style, rows = r, cols = 1, stack = TRUE)
            }
            # Thick border on last annotation row
            openxlsx::addStyle(wb, sheet, separator_style, rows = n_annot_rows,
                               cols = seq_along(df_cols), stack = TRUE)
        }

        # ---- Grouped contrast header row ----
        if (has_contrast_groups) {
            contrast_row <- n_annot_rows + 1
            contrast_style <- openxlsx::createStyle(
                textDecoration = "bold", halign = "center",
                border = "bottom", borderStyle = "thin"
            )
            for (contrast_name in names(contrast_groups)) {
                cols <- contrast_groups[[contrast_name]]
                openxlsx::writeData(wb, sheet, x = contrast_name,
                                    startCol = cols[1], startRow = contrast_row,
                                    colNames = FALSE, rowNames = FALSE)
                if (length(cols) > 1) {
                    openxlsx::mergeCells(wb, sheet,
                                         cols = cols, rows = contrast_row)
                }
                openxlsx::addStyle(wb, sheet, contrast_style,
                                   rows = contrast_row, cols = cols, stack = TRUE)
            }
        }

        # ---- Write data table (header + data) ----
        data_start_row <- n_annot_rows + n_contrast_header_rows + 1
        openxlsx::writeData(wb, sheet, df_out, startRow = data_start_row)

        # ---- Overwrite DE header cells with short stat names ----
        if (has_contrast_groups) {
            for (ci in de_col_indices) {
                stat_name <- sub("\\..*", "", colnames(df_out)[ci])
                openxlsx::writeData(wb, sheet, x = stat_name,
                                    startCol = ci, startRow = data_start_row,
                                    colNames = FALSE, rowNames = FALSE)
            }
        }

        # ---- Header styling ----
        header_style <- openxlsx::createStyle(textDecoration = "bold",
                                               border = "bottom", borderStyle = "thin")
        openxlsx::addStyle(wb, sheet, header_style, rows = data_start_row,
                           cols = seq_along(colnames(df_out)), stack = TRUE)

        # ---- Freeze panes ----
        # Freeze after annotation rows + header row, and after the ID column
        openxlsx::freezePane(wb, sheet,
                             firstActiveRow = data_start_row + 1,
                             firstActiveCol = 2)

        # ---- Auto-width for first column (ID) ----
        openxlsx::setColWidths(wb, sheet, cols = 1, widths = "auto")

        if (with_cutoffs) {
            add_cutoffs_sheet_legacy(wb, config, mode = mode)
            fill_manual_cutoffs_formulas_legacy(wb, sheet, df_out, config,
                                                 mode = mode,
                                                 start_row = data_start_row + 1)
        }
        if (isTRUE(provenance_sheet)) {
            add_provenance_sheet(wb, mode = mode)
        }
        openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
        path
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

            # 1. Hierarchical_Order (rank in hierarchical dendrogram)
            ranks <- match(de_df[[id_col]], ordered_ids)
            de_df$Hierarchical_Order <- ranks

            # 2. Add Z-score columns
            if (!is.null(excel_ord$zscore_mat)) {
                zmat <- excel_ord$zscore_mat
                idx <- match(de_df[[id_col]], rownames(zmat))
                z_sub <- zmat[idx, , drop = FALSE]
                de_df <- cbind(de_df, z_sub)
            }

            # 3. Partition clustering columns
            if (!is.null(excel_ord$partition_clusters)) {
                pc <- excel_ord$partition_clusters
                de_df$Partition_Cluster_ID <- pc[match(de_df[[id_col]], names(pc))]
                # Partition_Order: rank features by cluster ID, then Hierarchical_Order within
                sort_key <- order(
                    ifelse(is.na(de_df$Partition_Cluster_ID), Inf, de_df$Partition_Cluster_ID),
                    ifelse(is.na(de_df$Hierarchical_Order), Inf, de_df$Hierarchical_Order)
                )
                de_df$Partition_Order <- NA_integer_
                de_df$Partition_Order[sort_key] <- seq_len(nrow(de_df))
            }

            # 4. Binary pattern columns
            if (!is.null(excel_ord$binary_best)) {
                bp <- excel_ord$binary_best
                bp_idx <- match(de_df[[id_col]], bp$feature_id)
                de_df$Binary_Pattern <- bp$best_pattern[bp_idx]
                de_df$Binary_Corr <- bp$best_corr[bp_idx]
            }
        } else {
            # No clustering: add Hierarchical_Order with NA
            de_df$Hierarchical_Order <- NA_integer_
        }

        de_df <- add_default_order_if_missing(de_df, mat_de, id_col)
        de_df <- add_zscores_if_missing(de_df, mat_de, id_col)

        # Sort DE table by Hierarchical_Order to match heatmap row ordering
        if ("Hierarchical_Order" %in% names(de_df) && !all(is.na(de_df$Hierarchical_Order))) {
            de_df <- de_df[order(de_df$Hierarchical_Order, na.last = TRUE), , drop = FALSE]
        }

        # Reorder columns to: ID, annotations, expression, CV, DE_stats, clustering, z-scores
        id_cols <- id_col
        expr_cols <- colnames(mat_de)
        norm_cols_present <- intersect(paste0(expr_cols, ".norm"), names(de_df))
        mean_cols_present <- grep("^Mean\\.", names(de_df), value = TRUE)
        cv_cols_present <- grep("^CV\\.", names(de_df), value = TRUE)
        de_stat_cols <- grep("^(log2FC|linearFC|pvalue|padj|upDown)\\.", names(de_df), value = TRUE)
        clustering_cols <- intersect(
            c("Hierarchical_Order", "Partition_Cluster_ID", "Partition_Order",
              "Binary_Pattern", "Binary_Corr"),
            names(de_df)
        )
        zscore_cols <- grep("\\.zscore$", names(de_df), value = TRUE)

        # Annotation columns = everything not in ID, expression, normalized
        # expression, group means, CV, DE stats, clustering, z-scores, or 'order'
        all_known <- c(id_cols, expr_cols, norm_cols_present, mean_cols_present,
                       cv_cols_present, de_stat_cols, clustering_cols, zscore_cols, "order")
        annot_cols_present <- setdiff(names(de_df), all_known)

        # Check which expression columns are already present (from build_final_results_generic)
        expr_cols_present <- intersect(expr_cols, names(de_df))

        # Build desired column order: the measured block, then the block the
        # model saw, then the summaries derived from it, then the DE stats
        desired_order <- c(
            id_cols,
            annot_cols_present,
            expr_cols_present,
            norm_cols_present,
            mean_cols_present,
            cv_cols_present,
            de_stat_cols,
            clustering_cols,
            zscore_cols
        )

        # Reorder columns (only keep columns that exist)
        desired_order <- intersect(desired_order, names(de_df))
        de_df <- de_df[, desired_order, drop = FALSE]
    }

    # Save DE results and capture path
    f_de_created <- save_wb_results(de_df, f_de, with_cutoffs = FALSE)
    c(f_all_created, f_de_created)
}

#' Get standard column names for a contrast
#' @param contrast Contrast name
#' @param mode "proteomics" (uses .imputs.), "rna" or "metabolomics" (no .imputs.)
get_contrast_cols <- function(contrast, mode = "proteomics") {
    stopifnot(is.character(contrast), length(contrast) == 1, nzchar(contrast))

    # Proteomics DE summary strips spaces from contrast names; rna, metabolomics
    # and lipidomics keep them as-is.  Only normalize for proteomics.
    if (!mode %in% c("rna", "metabolomics", "lipidomics")) {
        contrast <- normalize_contrast_name(contrast)
    }

    # RNA-seq: no .imputs. infix; pass column is "<contrast>_pass"
    # (see build_rnaseq_summary_df() in R/domain/rnaseq/04_de_summary.R).
    if (mode == "rna") {
        list(
            log2fc = paste0("log2FC.", contrast),
            fc     = paste0("linearFC.", contrast),
            p      = paste0("pvalue.", contrast),
            padj   = paste0("padj.", contrast),
            pass   = paste0(contrast, "_pass"),
            updown = paste0("upDown.", contrast),
            manual = paste0("manual_cutoffs.", contrast)
        )
    } else if (mode %in% c("metabolomics", "lipidomics")) {
        # Metabolomics summary_df columns use pvalue./padj. prefixes
        # (see load_precomputed_metabolomics_de and build_de_summary in
        # R/domain/metabolomics/03_differential.R). Pass column is pass.<ctr>.
        list(
            fc     = paste0("linearFC.", contrast),
            p      = paste0("pvalue.", contrast),
            padj   = paste0("padj.", contrast),
            pass   = paste0("pass.", contrast),
            updown = paste0("upDown.", contrast),
            manual = paste0("manual_cutoffs.", contrast)
        )
    } else {
        # Proteomics (uses imputation naming)
        list(
            log2fc = paste0("log2FC.imputs.", contrast),
            fc     = paste0("linearFC.imputs.", contrast),
            p      = paste0("pvalue.imputs.", contrast),
            padj   = paste0("padj.imputs.", contrast),
            pass   = paste0("pass.imputs.", contrast),
            updown = paste0("upDown.imputs.", contrast),
            manual = paste0("manual_cutoffs.", contrast)
        )
    }
}

fill_manual_cutoffs_formulas_legacy <- function(wb, sheet, final_results, config,
                                                mode = "proteomics", start_row = 2) {
    de_cfg <- config$modes[[mode]]$de
    use_fdr <- isTRUE(de_cfg$use_adj_for_pass1)

    manual_cols <- grep("^manual_cutoffs\\.", names(final_results), value = TRUE)
    if (length(manual_cols) == 0) {
        return(invisible(NULL))
    }

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

#' Per-row coefficient of variation (CV%) on a linear matrix
#'
#' Arithmetic CV expressed as a percentage: \code{100 * sd(x) / mean(x)} per
#' row, using the sample (n-1) standard deviation (R's \code{stats::sd}).
#' Intended for LINEAR, positive data only — CV is not meaningful on log-scale
#' values (callers back-transform first).
#'
#' Edge cases (all return \code{NA_real_}, never 0 or Inf):
#' - rows with fewer than 2 non-missing values (CV/sample-SD undefined),
#' - rows whose mean is 0 or non-finite (division undefined).
#'
#' @param mat Numeric matrix (features x samples) on a linear scale. NAs are
#'   dropped per row (\code{na.rm = TRUE}).
#' @return Named numeric vector of CV% per row (names = \code{rownames(mat)}).
cv_percent <- function(mat) {
    mat <- as.matrix(mat)
    if (nrow(mat) == 0L) return(stats::setNames(numeric(0), rownames(mat)))

    n_obs <- rowSums(!is.na(mat))
    means <- rowMeans(mat, na.rm = TRUE)
    sds   <- apply(mat, 1, stats::sd, na.rm = TRUE)

    cv <- 100 * sds / means
    # Undefined cases -> NA (single/empty observations, or zero/non-finite mean)
    cv[n_obs < 2L] <- NA_real_
    cv[!is.finite(means) | means == 0] <- NA_real_
    cv[!is.finite(cv)] <- NA_real_

    if (!is.null(rownames(mat))) names(cv) <- rownames(mat)
    cv
}

#' Per-group coefficient of variation (CV%) columns for contrast groups
#'
#' Computes \code{CV.<group>} columns (percent CV across the biological
#' replicates of each group) for every group that appears in at least one
#' contrast, i.e. \code{union(contrasts_df$Numerator, contrasts_df$Denominator)}.
#' Group membership is resolved from \code{sample_meta}: samples sharing the
#' same value in the \code{Factor} column are biological replicates.
#'
#' CV is computed on \code{expr_linear}, which MUST already be on a linear,
#' positive scale (the caller is responsible for any back-transform — e.g.
#' linear CPM for RNA-seq, \code{2^x} for proteomics/metabolomics). Math and
#' edge-case handling are delegated to \code{\link{cv_percent}}.
#'
#' @param expr_linear Numeric matrix (features x samples), linear scale.
#'   Column names must match \code{sample_meta[[sample_id_col]]}.
#' @param sample_meta Sample metadata data.frame (one row per sample).
#' @param sample_id_col Column in \code{sample_meta} holding sample IDs that
#'   match \code{colnames(expr_linear)}.
#' @param contrasts_df Contrasts table with \code{Factor}, \code{Numerator},
#'   \code{Denominator} columns.
#' @param group_col Optional override for the grouping column. Defaults to the
#'   (single) value of \code{contrasts_df$Factor}.
#' @return A feature-indexed data.frame of \code{CV.<group>} columns
#'   (rownames = \code{rownames(expr_linear)}), or \code{NULL} if CV cannot be
#'   computed (missing inputs, no contrast groups, or ambiguous Factor).
compute_group_cv_columns <- function(expr_linear, sample_meta, sample_id_col,
                                     contrasts_df, group_col = NULL) {
    compute_group_stat_columns(
        expr          = expr_linear,
        sample_meta   = sample_meta,
        sample_id_col = sample_id_col,
        contrasts_df  = contrasts_df,
        stat_fn       = cv_percent,
        prefix        = "CV.",
        group_col     = group_col
    )
}

#' Per-group mean columns for contrast groups
#'
#' Companion to \code{\link{compute_group_cv_columns}}: emits
#' \code{Mean.<group>} columns, the arithmetic per-feature mean across the
#' replicates of each group that appears in at least one contrast.
#'
#' The mean is taken on whatever scale \code{expr} is on — the caller decides.
#' Pass the matrix the DE model actually saw (DESeq2-normalized counts for
#' RNA-seq, the imputed log2 matrix for proteomics) so that a reader can walk
#' from the per-sample values to the reported fold change.
#'
#' @param expr Numeric matrix (features x samples). Column names must match
#'   \code{sample_meta[[sample_id_col]]}.
#' @param sample_meta Sample metadata data.frame (one row per sample).
#' @param sample_id_col Column in \code{sample_meta} holding the sample IDs.
#' @param contrasts_df Contrasts table with \code{Factor}, \code{Numerator},
#'   \code{Denominator} columns.
#' @param group_col Optional override for the grouping column.
#' @return Feature-indexed data.frame of \code{Mean.<group>} columns, or NULL.
compute_group_mean_columns <- function(expr, sample_meta, sample_id_col,
                                       contrasts_df, group_col = NULL) {
    compute_group_stat_columns(
        expr          = expr,
        sample_meta   = sample_meta,
        sample_id_col = sample_id_col,
        contrasts_df  = contrasts_df,
        stat_fn       = function(m) rowMeans(m, na.rm = TRUE),
        prefix        = "Mean.",
        group_col     = group_col
    )
}

#' Per-group summary columns for contrast groups (shared machinery)
#'
#' Resolves the grouping column from the contrasts table, collects the samples
#' belonging to each group that appears in at least one contrast, and applies
#' \code{stat_fn} to that sub-matrix. Backs both
#' \code{\link{compute_group_cv_columns}} and
#' \code{\link{compute_group_mean_columns}} so the two can never disagree about
#' which samples make up a group.
#'
#' @param expr Numeric matrix (features x samples). Column names must match
#'   \code{sample_meta[[sample_id_col]]}.
#' @param sample_meta Sample metadata data.frame (one row per sample).
#' @param sample_id_col Column in \code{sample_meta} holding the sample IDs.
#' @param contrasts_df Contrasts table with \code{Factor}, \code{Numerator},
#'   \code{Denominator} columns.
#' @param stat_fn Function taking a features x replicates matrix and returning
#'   one value per feature (row).
#' @param prefix Column-name prefix for the emitted columns, e.g. \code{"CV."}.
#' @param group_col Optional override for the grouping column. Defaults to the
#'   (single) value of \code{contrasts_df$Factor}.
#' @return A feature-indexed data.frame of \code{<prefix><group>} columns
#'   (rownames = \code{rownames(expr)}), or \code{NULL} if the statistic cannot
#'   be computed (missing inputs, no contrast groups, or ambiguous Factor).
compute_group_stat_columns <- function(expr, sample_meta, sample_id_col,
                                       contrasts_df, stat_fn, prefix,
                                       group_col = NULL) {
    if (is.null(expr) || is.null(sample_meta) || is.null(sample_id_col)) {
        return(NULL)
    }
    if (!is.data.frame(contrasts_df) ||
        !all(c("Numerator", "Denominator") %in% colnames(contrasts_df))) {
        return(NULL)
    }
    if (!sample_id_col %in% colnames(sample_meta)) {
        return(NULL)
    }

    # Resolve the grouping column: explicit override, else the contrasts' Factor
    factor_col <- group_col %||% (if ("Factor" %in% colnames(contrasts_df)) {
        unique(as.character(contrasts_df$Factor))
    } else NULL)
    factor_col <- factor_col[!is.na(factor_col) & factor_col %in% colnames(sample_meta)]
    if (length(factor_col) != 1L) {
        warning(sprintf(
            "compute_group_stat_columns('%s'): could not resolve a single grouping column; skipping.",
            prefix
        ))
        return(NULL)
    }

    # Groups that appear in at least one contrast (Q5: union of Num/Den)
    groups <- unique(c(as.character(contrasts_df$Numerator),
                       as.character(contrasts_df$Denominator)))
    groups <- groups[!is.na(groups) & nzchar(groups)]
    if (length(groups) == 0L) return(NULL)

    expr <- as.matrix(expr)
    meta_ids <- as.character(sample_meta[[sample_id_col]])
    meta_grp <- as.character(sample_meta[[factor_col]])
    col_grp  <- meta_grp[match(colnames(expr), meta_ids)]

    stat_list <- lapply(groups, function(g) {
        cols <- which(col_grp == g)
        if (length(cols) == 0L) {
            return(rep(NA_real_, nrow(expr)))
        }
        unname(stat_fn(expr[, cols, drop = FALSE]))
    })
    names(stat_list) <- paste0(prefix, groups)

    stat_df <- as.data.frame(stat_list, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(stat_df) <- rownames(expr)
    stat_df
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
#' @param cv_cols Optional feature-indexed data.frame of per-group CV% columns
#'   (e.g. from \code{\link{compute_group_cv_columns}}). When supplied, these
#'   columns are inserted immediately after the per-sample expression block and
#'   before the per-contrast statistics, matched by feature ID.
#' @param norm_expr Optional second per-sample matrix (features x samples) — the
#'   values the DE model actually saw, as opposed to the raw/observed values in
#'   \code{expr_df}. Written as a block right after \code{expr_df}, with
#'   \code{norm_suffix} appended to each column name so the two blocks stay
#'   distinguishable. This is what lets a reader reproduce the reported fold
#'   change from the table itself.
#' @param norm_suffix Suffix appended to the \code{norm_expr} column names.
#' @param mean_cols Optional feature-indexed data.frame of per-group mean columns
#'   (e.g. from \code{\link{compute_group_mean_columns}}), inserted after the
#'   \code{norm_expr} block and before \code{cv_cols}.
#'
#' @return data.frame with ID, annotations, expression, [normalized expression],
#'   [Mean.<group>], [CV.<group>], DE stats, pass_any_contrast
build_final_results_generic <- function(
  summary_df,
  expr_df,
  contrasts_df,
  feature_id_col = "FeatureID",
  annot_cols = NULL,
  row_data = NULL,
  fc_is_signed = TRUE,
  fc_direction_col = NULL,
  mode = "proteomics",  # FIX 2: Add mode parameter for column naming
  cv_cols = NULL,
  norm_expr = NULL,
  norm_suffix = ".norm",
  mean_cols = NULL
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
    # ADD MODEL-INPUT ("NORMALIZED") EXPRESSION BLOCK
    # ============================================================
    # expr_df above is what was measured (raw counts / observed log2
    # intensities); this block is what the DE model was fitted on. Showing both
    # is the whole point: a fold change computed on normalized values cannot be
    # rederived from raw values alone, which is exactly the mismatch readers hit
    # when they eyeball the per-sample columns.

    if (!is.null(norm_expr)) {
        norm_df <- as.data.frame(norm_expr, check.names = FALSE)
        if (is.null(rownames(norm_df))) {
            warning("norm_expr has no rownames. Cannot add normalized expression values.")
        } else {
            norm_matched <- norm_df[match(base[[feature_id_col]], rownames(norm_df)), , drop = FALSE]
            colnames(norm_matched) <- paste0(colnames(norm_df), norm_suffix)
            rownames(norm_matched) <- NULL

            match_rate <- sum(base[[feature_id_col]] %in% rownames(norm_df)) / nrow(base)
            if (match_rate < 0.95) {
                warning(sprintf(
                    "Low match rate for normalized expression values: %.1f%% (%d/%d features matched). Check that norm_expr rownames match feature IDs.",
                    match_rate * 100, sum(base[[feature_id_col]] %in% rownames(norm_df)), nrow(base)
                ))
            }

            base <- cbind(base, norm_matched)
        }
    }

    # ============================================================
    # ADD PER-GROUP MEAN COLUMNS (on the model-input scale)
    # ============================================================

    if (!is.null(mean_cols) && is.data.frame(mean_cols) && ncol(mean_cols) > 0) {
        if (is.null(rownames(mean_cols))) {
            warning("mean_cols has no rownames. Cannot add per-group mean columns.")
        } else {
            mean_matched <- mean_cols[match(base[[feature_id_col]], rownames(mean_cols)), , drop = FALSE]
            rownames(mean_matched) <- NULL
            base <- cbind(base, mean_matched)
        }
    }

    # ============================================================
    # ADD PER-GROUP CV COLUMNS (linear-scale CV%, after expression)
    # ============================================================

    if (!is.null(cv_cols) && is.data.frame(cv_cols) && ncol(cv_cols) > 0) {
        if (is.null(rownames(cv_cols))) {
            warning("cv_cols has no rownames. Cannot add per-group CV columns.")
        } else {
            cv_matched <- cv_cols[match(base[[feature_id_col]], rownames(cv_cols)), , drop = FALSE]
            rownames(cv_matched) <- NULL
            base <- cbind(base, cv_matched)
        }
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
            # Loud signal instead of a silent all-NA fallback: when the pass
            # column the builder expects is absent, upDown would be blank for
            # every row even when fc/p/padj are present (the see-saw bug, where
            # get_contrast_cols() and the DE step disagree on the pass affix).
            warning(sprintf(
                "Contrast '%s' (mode '%s'): expected pass column '%s' not found in summary_df; upDown will be blank for all rows. Check get_contrast_cols() naming vs the DE step's pass column.",
                cn, mode, cols$pass
            ))
            rep(NA, length(m))
        }

        # Populate log2FC (when the mode's DE step supplies it), FC, p, padj.
        # log2FC comes first because it is the primitive the model reports;
        # linearFC is a presentation of it (signed reciprocal for FC < 1).
        if (!is.null(cols$log2fc) && cols$log2fc %in% colnames(summary_df)) {
            base[[cols$log2fc]] <- summary_df[[cols$log2fc]][m]
        }
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
  # If Hierarchical_Order is already populated (not all NA), keep it
  if ("Hierarchical_Order" %in% names(df) && !all(is.na(df$Hierarchical_Order))) {
    return(df)
  }
  # Legacy compat: also check old 'order' column
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


  # hclust requires >= 2 rows; assign rank 1 for a single feature
  if (nrow(mat_clean) < 2) {
    ranks <- if (nrow(mat_clean) == 1) match(df[[id_col]], rownames(mat_clean)) else rep(NA_integer_, nrow(df))
    if ("Hierarchical_Order" %in% names(df)) df$Hierarchical_Order <- ranks else df$order <- ranks
    return(df)
  }

  # Simple clustering
  dists <- dist(mat_clean)
  hc <- hclust(dists, method = "complete")

  # Create an order mapping
  ordered_ids <- rownames(mat_clean)[hc$order]
  ranks <- match(df[[id_col]], ordered_ids)


  # Use Hierarchical_Order if column exists, otherwise fall back to 'order'
  if ("Hierarchical_Order" %in% names(df)) {
    df$Hierarchical_Order <- ranks
  } else {
    df$order <- ranks
  }
  return(df)
}

# ============================================================
# Shiny-payload helpers
# ============================================================

# Keep these in sync with the sprintf() patterns in
# write_final_results_excels_legacy_generic() (lines 80-81).
.final_xlsx_all_pattern <- "Final_results_ALL_P_.*\\.xlsx$"
.final_xlsx_de_pattern  <- "Final_results_DE_P_.*\\.xlsx$"

#' Attach raw .xlsx bytes for Final_results_{ALL,DE} to a Shiny payload
#'
#' Adds two optional keys to a canonical Shiny payload:
#'   - \code{all_final_xlsx}: raw bytes of \code{Final_results_ALL_P_<p>.xlsx}
#'   - \code{de_final_xlsx}:  raw bytes of \code{Final_results_DE_P_<p>.xlsx}
#'
#' Both files are produced by
#' \code{\link{write_final_results_excels_legacy_generic}}. The Shiny app can
#' \code{writeBin()} either field to a temp file and re-open it as a workbook.
#'
#' Existing \code{de_final_table} (a data.frame) is left untouched. NULL-safe:
#' if a matching path is absent from \code{xlsx_files} or the file does not
#' exist on disk, the corresponding key stays NULL and a \code{message()} is
#' emitted.
#'
#' @param payload Named list (canonical Shiny payload).
#' @param xlsx_files Character vector of file paths returned by the legacy
#'   exports target (e.g. \code{rna_outputs_legacy}, \code{metab_final_results},
#'   or \code{excel_files} inside \code{mod_proteomics_exports}).
#' @return \code{payload} with \code{all_final_xlsx} / \code{de_final_xlsx}
#'   populated as \code{raw} vectors (or left NULL).
#' @export
attach_final_results_xlsx_bytes <- function(payload, xlsx_files) {
    if (is.null(xlsx_files) || length(xlsx_files) == 0) {
        message("[shiny payload] no xlsx_files provided; ",
                "all_final_xlsx / de_final_xlsx left NULL")
        return(payload)
    }
    xlsx_files <- as.character(xlsx_files)

    read_one <- function(pattern, key) {
        hits <- xlsx_files[grepl(pattern, basename(xlsx_files))]
        if (length(hits) == 0) {
            message(sprintf("[shiny payload] no match for %s; %s left NULL",
                            pattern, key))
            return(NULL)
        }
        path <- hits[[1]]
        if (!file.exists(path)) {
            message(sprintf("[shiny payload] file missing: %s; %s left NULL",
                            path, key))
            return(NULL)
        }
        readBin(path, what = "raw", n = file.info(path)$size)
    }

    payload$all_final_xlsx <- read_one(.final_xlsx_all_pattern, "all_final_xlsx")
    payload$de_final_xlsx  <- read_one(.final_xlsx_de_pattern,  "de_final_xlsx")
    payload
}
