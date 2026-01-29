#' Run DESeq2 differential expression analysis
#'
#' @param counts Integer count matrix (genes x samples)
#' @param meta Data frame with sample metadata
#' @param contrasts_df Data frame with columns: Contrast_name, Factor, Numerator, Denominator
#' @param de_cfg List with DE configuration (p_cutoff, linear_fc_cutoff, etc.)
#' @return List of DE result tables (one per contrast)
#' @export
run_deseq2_de <- function(counts, meta, contrasts_df, de_cfg) {
    # Validate inputs
    if (!is.matrix(counts) && !is.data.frame(counts)) {
        stop("counts must be a matrix or data frame")
    }
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

    # DESeq2 dataset
    dds0 <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts,
        colData   = as.data.frame(meta),
        design    = stats::as.formula(paste0("~ ", factor_col))
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
        contrast <- c(contrasts_df[i, "Factor"], contrasts_df[i, "Numerator"], contrasts_df[i, "Denominator"])
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

    # Extract thresholds from config
    # CRITICAL FIX: Use p_cutoff (standard key) instead of padj_cutoff
    padj_cutoff <- de_cfg$p_cutoff %||% de_cfg$padj_cutoff %||% 0.05
    linear_fc_cutoff <- de_cfg$linear_fc_cutoff %||% 1.5

    # Debug logging
    message(sprintf(
        "[build_rnaseq_summary_df] Using padj_cutoff=%.2f, linear_fc_cutoff=%.2f",
        padj_cutoff, linear_fc_cutoff
    ))

    # Get all genes
    all_genes <- unique(unlist(lapply(de_tables, function(x) x$FeatureID)))

    # Initialize summary data frame
    summary_df <- data.frame(
        FeatureID = all_genes,
        stringsAsFactors = FALSE
    )

    # Add columns for each contrast (using proteomics-compatible naming)
    for (cn in names(de_tables)) {
        tab <- de_tables[[cn]]

        # Match genes
        idx <- match(summary_df$FeatureID, tab$FeatureID)

        # Add linearFC (convert log2FC to linear FC for compatibility)
        fc_col <- paste0("linearFC.imputs.", cn)
        summary_df[[fc_col]] <- 2^tab$log2FoldChange[idx]

        # Add pvalue
        pval_col <- paste0("pvalue.imputs.", cn)
        summary_df[[pval_col]] <- tab$pvalue[idx]

        # Add padj
        padj_col <- paste0("padj.imputs.", cn)
        summary_df[[padj_col]] <- tab$padj[idx]

        # Add pass flag
        pass_col <- paste0(cn, "_pass")
        summary_df[[pass_col]] <- (
            !is.na(tab$padj[idx]) &
                tab$padj[idx] <= padj_cutoff &
                abs(tab$log2FoldChange[idx]) >= log2(linear_fc_cutoff)
        )
    }

    # Add pass_any_contrast column
    pass_cols <- grep("_pass$", colnames(summary_df), value = TRUE)
    if (length(pass_cols) > 0) {
        summary_df$pass_any_contrast <- apply(summary_df[, pass_cols, drop = FALSE], 1, any, na.rm = TRUE)
    } else {
        summary_df$pass_any_contrast <- FALSE
    }

    return(summary_df)
}
