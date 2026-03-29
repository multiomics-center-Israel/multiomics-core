#' Cell Type Deconvolution for RNA-seq
#'
#' Performs cell type deconvolution using xCell2 (or xCell fallback).
#' Only supported for Homo sapiens and Mus musculus.
#'
#' Adapted from multiomics framework for multiomics-core architecture.

# ==============================================================================
# ORGANISM GUARD
# ==============================================================================

#' Check if deconvolution is supported for organism
#'
#' @param organism Character, organism name (e.g., "Homo sapiens", "Giardia lamblia")
#' @return Logical TRUE only for human and mouse
#' @export
is_deconvolution_supported <- function(organism) {
    if (is.null(organism) || is.na(organism)) return(FALSE)
    organism_lower <- tolower(organism)
    organism_lower %in% c(
        "homo sapiens", "human", "hsapiens", "hsa",
        "mus musculus", "mouse", "mmusculus", "mmu"
    )
}

# ==============================================================================
# MAIN ENTRY POINT
# ==============================================================================

#' Run deconvolution analysis
#'
#' @param pre Preprocessed data list (expr_filt, meta, etc.)
#' @param de_res DE results list (summary_df, etc.)
#' @param cfg RNA-seq mode config (config$modes$rna)
#' @param out_dir Output directory
#' @return List with deconvolution results, or NULL if skipped
#' @export
run_deconvolution <- function(pre, de_res, cfg, out_dir) {
    # Get organism
    organism <- cfg$annotation$organism %||% "auto"

    # If auto, try to detect from data
    if (organism == "auto") {
        gene_ids <- rownames(pre$expr_filt)
        if (!is.null(gene_ids)) {
            org_detect <- tryCatch(detect_organism_from_ids(gene_ids), error = function(e) NULL)
            if (!is.null(org_detect) && org_detect$organism != "unknown") {
                organism <- org_detect$organism
            }
        }
    }

    organism_norm <- tryCatch(normalize_organism_name(organism), error = function(e) organism)

    if (!is_deconvolution_supported(organism_norm)) {
        message("Deconvolution not supported for organism: ", organism_norm,
                " (only human/mouse). Skipping.")
        return(list(
            status = "skipped",
            reason = paste0("Organism '", organism_norm, "' not supported (human/mouse only)")
        ))
    }

    message("=== Running Cell Type Deconvolution ===")
    message("Organism: ", organism_norm)

    deconv_dir <- file.path(out_dir, "Deconvolution")
    if (!dir.exists(deconv_dir)) dir.create(deconv_dir, recursive = TRUE)

    group_col <-  cfg$filtering$group_col %||% cfg$de_table$group_col %||% cfg$effects$color %||% "Condition"
    
    # Prepare expression for xCell2 (needs non-log CPM)
    expr_for_xcell <- prepare_expression_for_xcell(pre$expr_filt)

    # Determine if human or mouse
    is_human <- tolower(organism_norm) %in% c("homo sapiens", "human")

    # Run xCell2 or fallback
    xcell_results <- if (is_human) {
        run_xcell2_human(expr_for_xcell, deconv_dir)
    } else {
        run_xcell2_mouse(expr_for_xcell, deconv_dir)
    }

    if (is.null(xcell_results)) {
        message("Deconvolution failed. Skipping downstream analysis.")
        return(list(status = "failed", reason = "xCell2/xCell not available or failed"))
    }

    # Differential composition analysis
    diff_comp <- tryCatch(
        analyze_differential_composition(
            xcell_results$scores, pre$meta, group_col, deconv_dir
        ),
        error = function(e) {
            message("Differential composition failed: ", e$message)
            NULL
        }
    )

    # Generate plots
    figures <- tryCatch(
        generate_deconvolution_plots(
            xcell_results, diff_comp, pre$meta, group_col, deconv_dir
        ),
        error = function(e) {
            message("Deconvolution plot generation failed: ", e$message)
            list()
        }
    )

    # Save scores
    utils::write.csv(xcell_results$scores,
                     file.path(deconv_dir, "xcell2_scores.csv"),
                     row.names = FALSE)

    saveRDS(
        list(
            status = "completed",
            organism = organism_norm,
            scores = xcell_results$scores,
            differential = diff_comp,
            figures = figures
        ),
        file.path(deconv_dir, "deconvolution_results.rds")
    )

    message("=== Deconvolution Complete ===")

    list(
        status = "completed",
        organism = organism_norm,
        scores = xcell_results$scores,
        differential = diff_comp,
        figures = figures
    )
}

# ==============================================================================
# EXPRESSION PREPARATION
# ==============================================================================

#' Prepare expression matrix for xCell2 (non-log CPM)
#' @noRd
prepare_expression_for_xcell <- function(expr_filt) {
    message("Preparing expression matrix for xCell2...")

    # expr_filt is typically raw filtered counts — convert to CPM
    if (requireNamespace("edgeR", quietly = TRUE)) {
        cpm_mat <- edgeR::cpm(as.matrix(expr_filt), log = FALSE)
    } else {
        # Manual CPM
        lib_sizes <- colSums(expr_filt)
        cpm_mat <- sweep(as.matrix(expr_filt), 2, lib_sizes, "/") * 1e6
    }

    message("  Expression matrix: ", nrow(cpm_mat), " genes x ", ncol(cpm_mat), " samples")
    cpm_mat
}

# ==============================================================================
# xCell2 WRAPPERS
# ==============================================================================

#' Run xCell2 for human samples
#' @noRd
run_xcell2_human <- function(expr_mat, deconv_dir) {
    if (requireNamespace("xCell2", quietly = TRUE)) {
        message("Running xCell 2.0 (human)...")
        tryCatch({
            ref <- xCell2::BlueprintEncode
            scores <- xCell2::xCell2Analysis(expr_mat, ref = ref)
            scores_df <- as.data.frame(t(scores))
            scores_df$sample_id <- rownames(scores_df)
            return(list(scores = scores_df, reference = "BlueprintEncode"))
        }, error = function(e) {
            message("xCell2 failed: ", e$message, ". Trying legacy xCell...")
        })
    }

    # Fallback to legacy xCell
    if (requireNamespace("xCell", quietly = TRUE)) {
        message("Running legacy xCell (human)...")
        tryCatch({
            scores <- xCell::xCellAnalysis(expr_mat)
            scores_df <- as.data.frame(t(scores))
            scores_df$sample_id <- rownames(scores_df)
            return(list(scores = scores_df, reference = "xCell_legacy"))
        }, error = function(e) {
            message("xCell failed: ", e$message)
        })
    }

    message("Neither xCell2 nor xCell available.")
    NULL
}

#' Run xCell2 for mouse samples
#' @noRd
run_xcell2_mouse <- function(expr_mat, deconv_dir) {
    if (requireNamespace("xCell2", quietly = TRUE)) {
        message("Running xCell 2.0 (mouse)...")
        tryCatch({
            ref <- xCell2::ImmGen
            scores <- xCell2::xCell2Analysis(expr_mat, ref = ref)
            scores_df <- as.data.frame(t(scores))
            scores_df$sample_id <- rownames(scores_df)
            return(list(scores = scores_df, reference = "ImmGen"))
        }, error = function(e) {
            message("xCell2 mouse failed: ", e$message)
        })
    }

    message("xCell2 not available for mouse.")
    NULL
}

# ==============================================================================
# DIFFERENTIAL COMPOSITION
# ==============================================================================

#' Analyze differential cell type composition between conditions
#' @noRd
analyze_differential_composition <- function(scores, metadata, group_col, deconv_dir) {
    if (!(group_col %in% colnames(metadata))) {
        message("Condition column '", group_col, "' not found. Skipping differential analysis.")
        return(NULL)
    }

    message("Analyzing differential cell type composition...")

    # Merge scores with metadata
    sample_col <- "sample_id"
    if (!(sample_col %in% colnames(scores))) {
        scores$sample_id <- rownames(scores)
    }

    scores_meta <- merge(scores, metadata, by.x = "sample_id", by.y = "SampleID", all.x = TRUE)
    if (!(group_col %in% colnames(scores_meta))) return(NULL)

    celltype_cols <- setdiff(colnames(scores), "sample_id")
    conditions <- scores_meta[[group_col]]
    cond_levels <- unique(conditions)

    if (length(cond_levels) < 2) return(NULL)

    results <- lapply(celltype_cols, function(ct) {
        vals <- scores_meta[[ct]]
        if (all(is.na(vals))) return(NULL)

        tryCatch({
            if (length(cond_levels) == 2) {
                g1 <- vals[conditions == cond_levels[1]]
                g2 <- vals[conditions == cond_levels[2]]
                test <- wilcox.test(g2, g1)
                data.frame(
                    cell_type = ct,
                    comparison = paste0(cond_levels[2], "_vs_", cond_levels[1]),
                    pvalue = test$p.value,
                    stringsAsFactors = FALSE
                )
            } else {
                test <- kruskal.test(vals ~ factor(conditions))
                data.frame(
                    cell_type = ct,
                    comparison = "multi_group",
                    pvalue = test$p.value,
                    stringsAsFactors = FALSE
                )
            }
        }, error = function(e) NULL)
    })

    diff_df <- do.call(rbind, results[!vapply(results, is.null, logical(1))])
    if (is.null(diff_df) || nrow(diff_df) == 0) return(NULL)

    diff_df$padj <- p.adjust(diff_df$pvalue, method = "BH")
    diff_df <- diff_df[order(diff_df$padj), ]

    utils::write.csv(diff_df, file.path(deconv_dir, "differential_composition.csv"),
                     row.names = FALSE)

    message("  Found ", sum(diff_df$padj < 0.05, na.rm = TRUE),
            " significantly different cell types (FDR < 0.05)")
    diff_df
}

# ==============================================================================
# VISUALIZATION
# ==============================================================================

#' Generate deconvolution plots
#' @noRd
generate_deconvolution_plots <- function(xcell_results, diff_comp, metadata, group_col, deconv_dir) {
    message("Generating deconvolution plots...")
    plots_dir <- file.path(deconv_dir, "plots")
    if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

    figures <- list()
    scores <- xcell_results$scores
    celltype_cols <- setdiff(colnames(scores), "sample_id")

    # Heatmap (if ComplexHeatmap available)
    if (requireNamespace("ComplexHeatmap", quietly = TRUE) && length(celltype_cols) > 0) {
        tryCatch({
            score_mat <- as.matrix(scores[, celltype_cols, drop = FALSE])
            rownames(score_mat) <- scores$sample_id
            score_mat_scaled <- t(scale(t(score_mat)))

            fig_path <- file.path(plots_dir, "xcell2_heatmap.png")
            grDevices::png(fig_path, width = 12, height = 10, units = "in", res = 150)
            ht <- ComplexHeatmap::Heatmap(
                t(score_mat_scaled),
                name = "z-score",
                show_column_names = ncol(score_mat_scaled) <= 50,
                column_title = "Cell Type Enrichment Scores"
            )
            ComplexHeatmap::draw(ht)
            grDevices::dev.off()
            figures$heatmap <- fig_path
        }, error = function(e) {
            message("  Heatmap generation failed: ", e$message)
        })
    }

    # Stacked bar plot
    tryCatch({
        scores_long <- tidyr::pivot_longer(
            scores, cols = tidyselect::all_of(celltype_cols),
            names_to = "cell_type", values_to = "score"
        )
        scores_long <- merge(
            scores_long, metadata[, c("SampleID", group_col), drop = FALSE],
            by.x = "sample_id", by.y = "SampleID", all.x = TRUE
        )

        if (group_col %in% colnames(scores_long)) {
            summary_df <- stats::aggregate(
                score ~ cell_type + get(group_col),
                data = scores_long, FUN = mean, na.rm = TRUE
            )
            colnames(summary_df)[2] <- group_col

            p <- ggplot2::ggplot(summary_df,
                ggplot2::aes(x = .data[[group_col]], y = score, fill = cell_type)) +
                ggplot2::geom_bar(stat = "identity", position = "stack") +
                ggplot2::labs(title = "Cell Type Composition by Condition",
                              x = "Condition", y = "Score", fill = "Cell Type") +
                ggplot2::theme_minimal() +
                ggplot2::theme(legend.position = "right")

            fig_path <- file.path(plots_dir, "celltype_stacked_bar.png")
            ggplot2::ggsave(fig_path, p, width = 10, height = 8)
            figures$stacked_bar <- fig_path
        }
    }, error = function(e) {
        message("  Stacked bar plot failed: ", e$message)
    })

    message("Deconvolution plots saved to: ", plots_dir)
    figures
}
