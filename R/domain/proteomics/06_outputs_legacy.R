#' Orchestrator: write all proteomics multi-imputation outputs (legacy-compatible)
write_proteomics_multimpute_outputs <- function(pre, de_res, inputs, config, out_dir, excel_order = NULL, write_runs = FALSE) {
    files <- character(0)
    dirs <- create_legacy_output_dirs(out_dir)

    # 1) datasets
    runs_for_datasets <- if (isTRUE(write_runs)) de_res$runs else NULL
    files <- c(files, write_proteomics_datasets_legacy(pre, runs_for_datasets, config, out_dir))

    # 2) summary
    if (!is.null(de_res$summary_df)) {
        files <- c(files, write_limma_multimp_summary_legacy(de_res$summary_df, config, out_dir))
    }

    # 3) wide limma per contrast
    if (!is.null(de_res$runs_de_tables) && length(de_res$runs_de_tables) > 0) {
        contrast_names <- names(de_res$runs_de_tables[[1]])
        for (cn in contrast_names) {
            files <- c(files, write_limma_results_multimp_legacy(de_res = de_res, contrast_name = cn, config = config, out_dir = out_dir))
        }
    }

    # 4) final results TSV
    if (!is.null(inputs$contrasts) && !is.null(de_res$summary_df)) {
        final_results <- build_final_results_proteomics(
            pre = pre,
            summary_df = de_res$summary_df,
            contrasts_df = inputs$contrasts,
            row_data = pre$row_data,
            feature_id_col = config$modes$proteomics$de_table$id_col %||% "FeatureID"
        )
        files <- c(files, save_tsv(final_results, dirs$datasets, "final_results.tsv"))

        # 5) Excel outputs (delegated to dedicated file or assume available)
        if (exists("write_final_results_excels_proteomics")) {
            files <- c(files, write_final_results_excels_proteomics(
                final_results = final_results,
                pre = pre,
                config = config,
                out_dir = out_dir,
                excel_order = excel_order
            ))
        }
    }
    unique(files)
}

write_proteomics_datasets_legacy <- function(pre, runs = NULL, config, out_dir) {
    dirs <- create_legacy_output_dirs(out_dir)
    cfg <- config$modes$proteomics
    files <- character(0)

    files <- c(files, save_tsv(as.data.frame(pre$expr_filt, check.names = FALSE), dirs$datasets, "protein_log2_filtered_unimputed.tsv"))

    fname_imp <- sprintf("protein_log2_filtered_imputed_once_width_%s_shift_%s.tsv", cfg$imputation$width, cfg$imputation$downshift)
    files <- c(files, save_tsv(as.data.frame(pre$expr_imp_single, check.names = FALSE), dirs$datasets, fname_imp))

    if (!is.null(runs) && length(runs) > 0) {
        rep_dir <- file.path(dirs$datasets, "imputed_repetitions")
        ensure_dir(rep_dir)
        for (i in seq_along(runs)) {
            expr_i <- runs[[i]]$expr_imp
            if (is.null(expr_i)) next
            files <- c(files, save_tsv(as.data.frame(expr_i, check.names = FALSE), rep_dir, sprintf("protein_log2_filtered_imputed_%02d.tsv", i)))
        }
    }
    unique(files)
}

write_limma_multimp_summary_legacy <- function(summary_df, config, out_dir) {
    dirs <- create_legacy_output_dirs(out_dir)
    save_tsv(summary_df, dirs$datasets, sprintf("limma_multimp_summary_p%s.tsv", p_tag_generic(config, "proteomics")))
}

write_limma_results_multimp_legacy <- function(de_res, contrast_name, config, out_dir) {
    dirs <- create_legacy_output_dirs(out_dir)
    feature_id_col <- config$modes$proteomics$de_table$id_col %||% "FeatureID"
    wide_df <- build_limma_results_multimp_wide(
        runs_de_tables = de_res$runs_de_tables,
        contrast_name = contrast_name,
        feature_id_col = feature_id_col
    )
    fname <- sprintf("limma_results_multimp_p%s.tsv", p_tag_generic(config, "proteomics"))
    save_tsv(wide_df, dirs$datasets, fname)
}

build_limma_results_multimp_wide <- function(runs_de_tables, contrast_name, stats_cols = c("logFC", "P.Value", "adj.P.Val"), feature_id_col = "FeatureID") {
    stopifnot(length(runs_de_tables) >= 1)
    base <- runs_de_tables[[1]][[contrast_name]]
    if (is.null(base)) stop(sprintf("Contrast '%s' missing in imputation run 1", contrast_name))

    if (!feature_id_col %in% colnames(base)) {
        stop(sprintf("Feature ID column '%s' missing in base table", feature_id_col))
    }

    id_cols <- intersect(c(feature_id_col, "Protein.Names", "Genes", "First.Protein.Description", "Contrast"), colnames(base))
    out <- base[, id_cols, drop = FALSE]

    for (i in seq_along(runs_de_tables)) {
        tab <- runs_de_tables[[i]][[contrast_name]]
        if (is.null(tab)) stop(sprintf("Critical: Contrast '%s' is missing in imputation run %d", contrast_name, i))

        tab <- align_de_table_by_feature_id(tab = tab, ref_ids = out[[feature_id_col]], run_i = i, contrast_name = contrast_name, id_col = feature_id_col)
        stat_block <- tab[, intersect(stats_cols, colnames(tab)), drop = FALSE]
        colnames(stat_block) <- paste0(colnames(stat_block), ".", i)
        out <- cbind(out, stat_block)
    }
    out
}

build_final_results_proteomics <- function(pre, summary_df, contrasts_df, row_data = NULL, feature_id_col = "FeatureID") {
    build_final_results_generic(
        summary_df = summary_df,
        expr_df = pre$expr_filt,
        contrasts_df = contrasts_df,
        feature_id_col = feature_id_col,
        annot_cols = c(
            "Protein.Names" = "Protein.Names",
            "Genes" = "Genes",
            "First.Protein.Description" = "First.Protein.Description"
        ),
        row_data = row_data %||% pre$row_data,
        fc_is_signed = TRUE # linearFC is signed
    )
}
