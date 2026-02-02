#' Writes intermediate RNA matrices (legacy-like)
write_rnaseq_datasets_legacy <- function(pre, config, out_dir) {
    dirs <- create_legacy_output_dirs(out_dir)
    files <- character(0)

    # FIX 3: Add 'gene' column as first column (convert rownames to explicit column)
    # FIX 4: Rename output from rna_expr_work_* to rna_norm_*
    counts_filt_df <- as.data.frame(pre$expr_filt, check.names = FALSE)
    counts_filt_df <- cbind(gene = rownames(counts_filt_df), counts_filt_df, stringsAsFactors = FALSE)
    rownames(counts_filt_df) <- NULL

    norm_df <- as.data.frame(pre$expr_work, check.names = FALSE)
    norm_df <- cbind(gene = rownames(norm_df), norm_df, stringsAsFactors = FALSE)
    rownames(norm_df) <- NULL

    files <- c(files, save_tsv(counts_filt_df, dirs$datasets, "rna_counts_filtered.tsv"))
    files <- c(files, save_tsv(norm_df, dirs$datasets, sprintf("rna_norm_%s.tsv", attr(pre$expr_work, "method") %||% "norm")))

    unique(files)
}

write_rnaseq_outputs_legacy <- function(pre, de_res, inputs, config, out_dir, clustering_res = NULL) {
    files <- character(0)
    dirs <- create_legacy_output_dirs(out_dir)

    # 1) datasets
    files <- c(files, write_rnaseq_datasets_legacy(pre, config, out_dir))

    # 2) summary_df
    summary_df <- build_rnaseq_summary_df(de_res$tables, config$modes$rna$de)
    # Ensure ID column is named "Gene" for compatibility with generic export (which expects 'Gene')
    if ("FeatureID" %in% names(summary_df) && !"Gene" %in% names(summary_df)) {
        names(summary_df)[names(summary_df) == "FeatureID"] <- "Gene"
    }

    # ---- Validation / Logging ----
    p_cutoff <- config$modes$rna$de$p_cutoff %||% 0.05
    n_total <- nrow(summary_df)
    n_sig <- if ("pass_any_contrast" %in% names(summary_df)) {
        sum(summary_df$pass_any_contrast == 1, na.rm = TRUE)
    } else {
        NA
    }

    message(sprintf("[RNA-seq Export] Total Genes: %d", n_total))
    message(sprintf("[RNA-seq Export] Significant DE Genes (p < %.2f): %s", p_cutoff, if (is.na(n_sig)) "N/A" else n_sig))

    if (n_total == 0) {
        warning("[RNA-seq Export] CRITICAL: Summary DataFrame is empty!")
    } else if (!is.na(n_sig) && n_sig == 0) {
        warning("[RNA-seq Export] No significant DE genes found across any contrast.")
    }

    files <- c(files, save_tsv(summary_df, dirs$datasets, sprintf("deseq2_summary_p%s.tsv", p_tag_generic(config, "rna"))))

    # 2b) DE summary counts table (up/down/any per contrast)
    if (!is.null(inputs$contrasts)) {
        summary_counts <- build_de_summary_counts(summary_df, inputs$contrasts)
        files <- c(files, save_tsv(summary_counts, dirs$datasets, "de_summary_counts.tsv"))
    }

    # 3) final results TSV
    if (!is.null(inputs$contrasts)) {
        final_results <- build_final_results_rnaseq(pre = pre, summary_df = summary_df, contrasts_df = inputs$contrasts, row_data = pre$row_data)
        files <- c(files, save_tsv(final_results, dirs$datasets, "final_results.tsv"))

        # 4) Excel outputs
        # Extract excel_order from clustering_res if available
        excel_order <- if (!is.null(clustering_res)) clustering_res$excel_order else NULL

        if (exists("write_final_results_excels_legacy_rna")) {
            files <- c(files, write_final_results_excels_legacy_rna(final_results, pre, config, out_dir, excel_order = excel_order))
        } else if (exists("write_final_results_excels_legacy_generic")) {
            # Use generic if wrapper not available
            # Pass RAW counts (expr_filt) and clustering results
            files <- c(files, write_final_results_excels_legacy_generic(
                final_results = final_results,
                config = config,
                out_dir = out_dir,
                mode = "rna",
                id_col = "Gene",
                expr_for_de = pre$expr_filt,  # RAW counts, not normalized
                with_cutoffs = TRUE,
                clustering_res = list(excel_order = excel_order)
            ))
        }
    }
    unique(files)
}


build_final_results_rnaseq <- function(pre, summary_df, contrasts_df, row_data = NULL) {
    # FIX 1: Remove annotation columns (gene_name, symbol, description) from final_results
    # FIX 2: Pass mode="rna" for correct column naming (no ".imputs.")
    # FIX 5: Use RAW counts (expr_filt) not normalized (expr_work)
    build_final_results_generic(
        summary_df = summary_df,
        expr_df = pre$expr_filt,  # FIX 5: RAW counts for final results tables
        contrasts_df = contrasts_df,
        feature_id_col = "Gene",
        annot_cols = NULL,  # FIX 1: No annotation columns in RNA final_results
        row_data = row_data %||% pre$row_data,
        fc_is_signed = TRUE,  # log2FoldChange is signed
        mode = "rna"  # FIX 2: Use RNA column naming
    )
}
