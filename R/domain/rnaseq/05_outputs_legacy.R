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
    build_final_results_generic(
        summary_df = summary_df,
        expr_df = pre$expr_work,
        contrasts_df = contrasts_df,
        feature_id_col = "Gene",
        annot_cols = c(
            "gene_name" = "gene_name",
            "symbol" = "symbol",
            "description" = "description"
        ),
        row_data = row_data %||% pre$row_data,
        fc_is_signed = TRUE # log2FoldChange is signed
    )
}
