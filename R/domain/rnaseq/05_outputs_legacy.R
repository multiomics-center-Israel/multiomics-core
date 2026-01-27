#' Writes intermediate RNA matrices (legacy-like)
write_rnaseq_datasets_legacy <- function(pre, config, out_dir) {
    dirs <- create_legacy_output_dirs(out_dir)
    files <- character(0)

    files <- c(files, save_tsv(as.data.frame(pre$expr_filt, check.names = FALSE), dirs$datasets, "rna_counts_filtered.tsv"))
    files <- c(files, save_tsv(as.data.frame(pre$expr_work, check.names = FALSE), dirs$datasets, sprintf("rna_expr_work_%s.tsv", attr(pre$expr_work, "method") %||% "norm")))

    unique(files)
}

write_rnaseq_outputs_legacy <- function(pre, de_res, inputs, config, out_dir) {
    files <- character(0)
    dirs <- create_legacy_output_dirs(out_dir)

    # 1) datasets
    files <- c(files, write_rnaseq_datasets_legacy(pre, config, out_dir))

    # 2) summary_df
    summary_df <- build_rnaseq_summary_df(de_res, config)
    files <- c(files, save_tsv(summary_df, dirs$datasets, sprintf("deseq2_summary_p%s.tsv", p_tag_generic(config, "rna"))))

    # 3) final results TSV
    if (!is.null(inputs$contrasts)) {
        final_results <- build_final_results_rnaseq(pre = pre, summary_df = summary_df, contrasts_df = inputs$contrasts, row_data = pre$row_data)
        files <- c(files, save_tsv(final_results, dirs$datasets, "final_results.tsv"))

        # 4) Excel outputs
        if (exists("write_final_results_excels_legacy_rna")) {
            files <- c(files, write_final_results_excels_legacy_rna(final_results, pre, config, out_dir))
        } else if (exists("write_final_results_excels_legacy_generic")) {
            # Use generic if wrapper not available
            files <- c(files, write_final_results_excels_legacy_generic(
                final_results = final_results, config = config, out_dir = out_dir, mode = "rna", id_col = "Gene", expr_for_de = pre$expr_work, with_cutoffs = TRUE
            ))
        }
    }
    unique(files)
}


build_final_results_rnaseq <- function(pre, summary_df, contrasts_df, row_data = NULL) {
    expr_df <- as.data.frame(pre$expr_work, check.names = FALSE)
    gene_ids <- summary_df$FeatureID

    base <- data.frame(Gene = gene_ids, stringsAsFactors = FALSE, check.names = FALSE)

    if (!is.null(row_data)) {
        for (col in intersect(c("gene_name", "symbol", "description", "GeneSymbol", "Description"), colnames(row_data))) {
            if (!is.null(rownames(row_data)) && all(gene_ids %in% rownames(row_data))) {
                base[[col]] <- row_data[gene_ids, col]
            }
        }
    }

    if (!is.null(rownames(pre$expr_work))) {
        expr_df2 <- expr_df[match(gene_ids, rownames(pre$expr_work)), , drop = FALSE]
    } else {
        stop("pre$expr_work must have rownames for final results.")
    }
    base <- cbind(base, expr_df2)

    contrast_names <- contrasts_df$Contrast_name
    for (cn in contrast_names) {
        cols <- get_contrast_cols(cn)
        needed <- c(cols$fc, cols$p, cols$padj)
        if (!all(needed %in% colnames(summary_df))) stop("Missing columns in summary_df")
    }

    m <- match(base$Gene, summary_df$FeatureID)
    for (cn in contrast_names) {
        cols <- get_contrast_cols(cn)
        fc_vals <- summary_df[[cols$fc]][m]
        pass_vals <- if (cols$pass %in% colnames(summary_df)) summary_df[[cols$pass]][m] else rep(NA, length(m))

        base[[cols$fc]] <- fc_vals
        base[[cols$p]] <- summary_df[[cols$p]][m]
        base[[cols$padj]] <- summary_df[[cols$padj]][m]
        base[[cols$updown]] <- ifelse(!is.na(pass_vals), ifelse(as.numeric(fc_vals) >= 0, "up", "down"), "")
        base[[cols$manual]] <- NA
    }

    pass_cols <- paste0("pass.imputs.", contrast_names)
    existing_pass_cols <- intersect(pass_cols, colnames(summary_df))

    if (length(existing_pass_cols) > 0) {
        pass_mat <- summary_df[m, existing_pass_cols, drop = FALSE]
        base$pass_any_contrast <- ifelse(rowSums(!is.na(pass_mat)) > 0, 1, NA)
    } else {
        base$pass_any_contrast <- NA
    }
    base
}
