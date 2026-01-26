#' Run DESeq2 differential expression
#'
#' @param expr_filt Filtered counts matrix (integers)
#' @param meta Metadata dataframe
#' @param inputs rna inputs list (must contain contrasts)
#' @param config full config
run_deseq2_de <- function(expr_filt, meta, inputs, config) {
    cfg <- config$modes$rna
    contrasts_df <- inputs$contrasts

    if (is.null(contrasts_df) || nrow(contrasts_df) == 0) {
        stop("No contrasts provided for RNA DE.")
    }

    # Prepare counts (integer)
    counts <- expr_filt
    storage.mode(counts) <- "integer"

    # Factor column
    factor_col <- unique(contrasts_df$Factor)
    if (length(factor_col) != 1) stop("RNA contrasts must use a single Factor.")
    factor_col <- factor_col[[1]]

    if (!factor_col %in% colnames(meta)) stop("meta missing Factor column: ", factor_col)
    meta[[factor_col]] <- as.factor(meta[[factor_col]])

    # DESeq2 dataset
    dds0 <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts,
        colData   = as.data.frame(meta),
        design    = stats::as.formula(paste0("~ ", factor_col))
    )

    # Run DESeq
    dds <- DESeq2::DESeq(dds0)

    # Extract results per contrast
    tables <- list()
    for (i in seq_len(nrow(contrasts_df))) {
        cn <- contrasts_df$Contrast_name[i]
        num <- contrasts_df$Numerator[i]
        den <- contrasts_df$Denominator[i]

        res <- DESeq2::results(
            dds,
            contrast = c(factor_col, num, den),
            independentFiltering = FALSE,
            cooksCutoff = FALSE
        )

        tab <- as.data.frame(res)
        tab$FeatureID <- rownames(tab)
        tab$Contrast <- cn

        # Keep relevant columns
        cols <- intersect(c("FeatureID", "Contrast", "log2FoldChange", "pvalue", "padj"), colnames(tab))
        tab <- tab[, cols, drop = FALSE]
        tables[[cn]] <- tab
    }

    list(
        method = "DESeq2",
        tables = tables,
        contrasts = contrasts_df,
        dds = dds, # Store for legacy export
        info = list(
            design    = paste0("~", factor_col),
            n_genes   = nrow(dds),
            n_samples = ncol(dds)
        )
    )
}

build_rnaseq_summary_df <- function(de_res, config) {
    stopifnot(is.list(de_res), !is.null(de_res$tables), !is.null(de_res$contrasts))
    contrasts_df <- de_res$contrasts

    de_cfg <- config$modes$rna$de %||% list()
    padj_cutoff <- as.numeric(de_cfg$padj_cutoff %||% 0.05)
    linear_fc_cutoff <- as.numeric(de_cfg$linear_fc_cutoff %||% 1.5)

    cn0 <- contrasts_df$Contrast_name[[1]]
    tab0 <- de_res$tables[[cn0]]
    if (is.null(tab0)) stop("RNA DE tables missing first contrast: ", cn0)

    base_ids <- tab0$FeatureID
    if (is.null(base_ids)) stop("FeatureID missing in first contrast table")

    out <- data.frame(FeatureID = base_ids, stringsAsFactors = FALSE, check.names = FALSE)
    pass_counts <- integer(length(base_ids))

    for (i in seq_len(nrow(contrasts_df))) {
        cn <- contrasts_df$Contrast_name[i]
        tab <- de_res$tables[[cn]]

        tab <- tab[match(base_ids, tab$FeatureID), , drop = FALSE]

        log2fc <- tab$log2FoldChange
        pval <- tab$pvalue
        padj <- tab$padj

        linearFC <- ifelse(is.na(log2fc), NA_real_, ifelse(log2fc > 0, 2^log2fc, -1 / (2^log2fc)))
        linearRatio <- ifelse(is.na(log2fc), NA_real_, 2^log2fc)

        pass_dir <- ifelse(!is.na(padj) & (padj <= padj_cutoff) & (abs(linearFC) >= linear_fc_cutoff),
            ifelse(linearFC > 0, "up", "down"), NA_character_
        )

        out[[paste0("sum.pass.", cn)]] <- ifelse(is.na(pass_dir), 0, 1)
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
