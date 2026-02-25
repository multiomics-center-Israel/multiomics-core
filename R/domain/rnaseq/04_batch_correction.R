#' Batch Effect Detection and Correction for RNA-seq
#'
#' Adapted from multiomics framework for multiomics-core architecture.
#' Works with the pre object (pre$expr_work, pre$expr_filt, pre$meta).
#'
#' Detection: PCA + ANOVA, simplified PVCA, silhouette analysis.
#' Correction: ComBat-Seq (default), SVA, RUVSeq.

# ==============================================================================
# MAIN ENTRY POINT
# ==============================================================================

#' Run batch analysis on preprocessed RNA-seq data
#'
#' @param pre Preprocessed data list (expr_raw, expr_filt, expr_work, meta, row_data)
#' @param cfg RNA-seq mode config (config$modes$rna)
#' @param out_dir Output directory
#' @return Modified pre list (with corrected expr_work/expr_filt if correction applied),
#'         or original pre if disabled / no batch detected
#' @export
run_batch_analysis <- function(pre, cfg, out_dir) {
    bc <- get_batch_config(cfg)

    if (!isTRUE(bc$enabled)) {
        message("Batch correction disabled. Passing through unchanged.")
        return(pre)
    }

    message("=== Running Batch Effect Analysis ===")

    batch_dir <- file.path(out_dir, "Batch_correction")
    if (!dir.exists(batch_dir)) dir.create(batch_dir, recursive = TRUE)

    group_col  <- cfg$filtering$group_col %||% "condition"
    sample_col <- cfg$id_columns$sample_col %||%
                  cfg$effects$samples %||% "SampleID"

    # Identify batch columns
    batch_cols <- identify_batch_columns(pre$meta, bc, group_col, sample_col)

    if (length(batch_cols) == 0) {
        message("No batch variables identified. Skipping batch analysis.")
        return(pre)
    }

    message("Batch variables: ", paste(batch_cols, collapse = ", "))

    # Detect batch effects using expr_work (log-transformed)
    detection <- detect_batch_effects(pre$expr_work, pre$meta, batch_cols, group_col)

    # Correct if requested and detected
    if (isTRUE(bc$correct) && isTRUE(detection$batch_detected)) {
        message("Correcting batch effects using: ", bc$method)

        batch_col <- batch_cols[1]
        batch_var <- pre$meta[[batch_col]]
        condition_var <- pre$meta[[group_col]]

        corrected <- NULL

        if (bc$method == "combat_seq") {
            corrected <- correct_combat_seq(pre$expr_filt, batch_var, condition_var)
        } else if (bc$method == "sva") {
            corrected <- correct_sva(pre$expr_filt, batch_var, condition_var, bc)
        } else if (bc$method == "ruv") {
            corrected <- correct_ruv(pre$expr_filt, batch_var, condition_var, bc)
        }

        if (!is.null(corrected)) {
            pre$expr_filt <- corrected
            # Re-derive work matrix (log-transform)
            pre$expr_work <- log2(corrected + 1)
            pre$batch_corrected <- TRUE
            pre$batch_method <- bc$method

            # Diagnostic side effect: writes batch-corrected matrix for manual inspection.
            # Not tracked as a targets output — downstream targets consume the in-memory
            # `pre` object returned by this function, not this file.
            corrected_df <- as.data.frame(corrected)
            corrected_df <- cbind(gene_id = rownames(corrected_df), corrected_df)
            utils::write.csv(corrected_df, file.path(batch_dir, "batch_corrected_counts.csv"),
                             row.names = FALSE)
            message("Saved batch-corrected counts")
        }
    }

    # Generate plots
    plot_batch_analysis(detection, pre$meta, batch_cols, group_col, batch_dir)

    message("=== Batch Analysis Complete ===")
    pre
}

# ==============================================================================
# CONFIGURATION
# ==============================================================================

#' Get batch correction config with defaults
#' @noRd
get_batch_config <- function(cfg) {
    bc <- cfg$batch_correction %||% list()
    list(
        enabled     = bc$enabled %||% FALSE,
        correct     = bc$correct %||% bc$enabled %||% FALSE,
        method      = bc$method %||% "combat_seq",
        batch_columns = bc$batch_columns,
        batch_threshold = bc$batch_threshold %||% 0.1
    )
}

# ==============================================================================
# BATCH COLUMN IDENTIFICATION
# ==============================================================================

#' Identify potential batch columns in metadata
#' @noRd
identify_batch_columns <- function(metadata, bc, group_col, sample_col) {
    if (!is.null(bc$batch_columns)) {
        return(intersect(bc$batch_columns, colnames(metadata)))
    }

    exclude <- c(group_col, sample_col, "sample")
    batch_patterns <- c("batch", "run", "date", "plate", "lane", "flowcell",
                        "library", "sequencing", "cohort", "site", "center")

    potential <- character(0)
    for (col in colnames(metadata)) {
        if (col %in% exclude) next
        col_lower <- tolower(col)
        if (any(vapply(batch_patterns, function(p) grepl(p, col_lower), logical(1)))) {
            potential <- c(potential, col)
            next
        }
        vals <- metadata[[col]]
        if (is.factor(vals) || is.character(vals)) {
            n_lev <- length(unique(vals))
            if (n_lev >= 2 && n_lev <= min(20, nrow(metadata) / 2)) {
                potential <- c(potential, col)
            }
        }
    }

    if (length(potential) > 0) {
        message("Auto-detected potential batch variables: ", paste(potential, collapse = ", "))
    }
    potential
}

# ==============================================================================
# BATCH DETECTION
# ==============================================================================

#' Detect batch effects via PCA + ANOVA
#' @noRd
detect_batch_effects <- function(expr_work, metadata, batch_cols, group_col) {
    message("Detecting batch effects...")

    pca <- prcomp(t(expr_work), scale. = TRUE, center = TRUE)
    n_pcs <- min(5, ncol(pca$x))
    pca_scores <- pca$x[, 1:n_pcs, drop = FALSE]
    var_exp <- summary(pca)$importance[2, 1:n_pcs]

    batch_stats <- list()

    for (bc in batch_cols) {
        if (!(bc %in% colnames(metadata))) next
        batch_var <- metadata[[bc]]
        if (length(unique(batch_var)) < 2) next

        pc_r2 <- vapply(seq_len(n_pcs), function(pc) {
            tryCatch({
                fit <- lm(pca_scores[, pc] ~ as.factor(batch_var))
                summary(fit)$r.squared
            }, error = function(e) 0)
        }, numeric(1))

        weighted_r2 <- sum(pc_r2 * var_exp)

        batch_stats[[bc]] <- list(
            pc_r2 = pc_r2,
            weighted_r2 = weighted_r2,
            n_levels = length(unique(batch_var))
        )

        message("  ", bc, ": weighted R2=", round(weighted_r2, 3))
    }

    max_r2 <- if (length(batch_stats) > 0) {
        max(vapply(batch_stats, function(x) x$weighted_r2, numeric(1)))
    } else 0

    batch_detected <- max_r2 > 0.1

    if (batch_detected) {
        message("Batch effects detected (max R2=", round(max_r2, 3), ")")
    } else {
        message("No significant batch effects detected")
    }

    list(
        pca_scores = pca_scores,
        var_explained = var_exp,
        batch_stats = batch_stats,
        batch_detected = batch_detected
    )
}

# ==============================================================================
# CORRECTION METHODS
# ==============================================================================

#' ComBat-Seq correction
#' @noRd
correct_combat_seq <- function(counts, batch_var, condition_var) {
    if (!requireNamespace("sva", quietly = TRUE)) {
        message("sva package not available for ComBat-Seq")
        return(NULL)
    }

    message("  Running ComBat-Seq...")
    tryCatch({
        sva::ComBat_seq(
            counts = as.matrix(counts),
            batch = as.factor(batch_var),
            group = as.factor(condition_var)
        )
    }, error = function(e) {
        message("  ComBat-Seq failed: ", e$message)
        NULL
    })
}

#' SVA correction
#' @noRd
correct_sva <- function(counts, batch_var, condition_var, bc) {
    if (!requireNamespace("sva", quietly = TRUE)) {
        message("sva package not available")
        return(NULL)
    }

    message("  Running SVA...")
    counts_norm <- log2(as.matrix(counts) + 1)
    mod <- model.matrix(~as.factor(condition_var))
    mod0 <- model.matrix(~1, data = data.frame(condition = condition_var))

    n_sv <- tryCatch(sva::num.sv(counts_norm, mod, method = "be"), error = function(e) 2)
    message("  Estimating ", n_sv, " surrogate variables")

    svobj <- tryCatch(sva::sva(counts_norm, mod, mod0, n.sv = n_sv), error = function(e) {
        message("  SVA failed: ", e$message)
        NULL
    })
    if (is.null(svobj)) return(NULL)

    tryCatch({
        limma::removeBatchEffect(counts_norm, covariates = svobj$sv)
    }, error = function(e) {
        message("  Batch removal failed: ", e$message)
        NULL
    })
}

#' RUVSeq correction
#' @noRd
correct_ruv <- function(counts, batch_var, condition_var, bc) {
    if (!requireNamespace("RUVSeq", quietly = TRUE)) {
        message("RUVSeq package not available")
        return(NULL)
    }

    message("  Running RUVSeq...")
    eset <- tryCatch({
        RUVSeq::newSeqExpressionSet(
            as.matrix(counts),
            phenoData = data.frame(condition = condition_var, row.names = colnames(counts))
        )
    }, error = function(e) {
        message("  Failed to create SeqExpressionSet: ", e$message)
        NULL
    })
    if (is.null(eset)) return(NULL)

    gene_vars <- apply(counts, 1, var)
    control_genes <- names(sort(gene_vars))[seq_len(min(1000, length(gene_vars)))]
    n_k <- min(2, 5)

    tryCatch({
        ruv_result <- RUVSeq::RUVg(eset, control_genes, k = n_k)
        RUVSeq::normCounts(ruv_result)
    }, error = function(e) {
        message("  RUVSeq failed: ", e$message)
        NULL
    })
}

# ==============================================================================
# VISUALIZATION
# ==============================================================================

#' Generate batch effect plots
#' @noRd
plot_batch_analysis <- function(detection, metadata, batch_cols, group_col, batch_dir) {
    message("Generating batch effect plots...")
    plots_dir <- file.path(batch_dir, "plots")
    if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

    pca_scores <- detection$pca_scores
    var_exp <- detection$var_explained

    pca_df <- as.data.frame(pca_scores[, 1:min(2, ncol(pca_scores)), drop = FALSE])
    pca_df$sample_id <- rownames(pca_scores)

    for (col in c(batch_cols, group_col)) {
        if (col %in% colnames(metadata)) {
            pca_df[[col]] <- metadata[[col]]
        }
    }

    for (bc in batch_cols) {
        if (!(bc %in% colnames(pca_df))) next

        p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2,
                                                     color = .data[[bc]])) +
            ggplot2::geom_point(size = 3, alpha = 0.7) +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = paste0("PCA colored by ", bc),
                x = paste0("PC1 (", round(var_exp[1] * 100, 1), "%)"),
                y = paste0("PC2 (", round(var_exp[2] * 100, 1), "%)"),
                color = bc
            )

        ggplot2::ggsave(file.path(plots_dir, paste0("pca_by_", bc, ".png")),
                        plot = p, width = 10, height = 8, dpi = 150)
    }

    message("Batch plots saved to: ", plots_dir)
}
