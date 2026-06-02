#' Writes intermediate RNA matrices (legacy-like)
write_rnaseq_datasets_legacy <- function(pre, config, dirs ) {
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
    files <- c(files, write_rnaseq_datasets_legacy(pre, config, dirs))

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
        final_results <- build_final_results_rnaseq(pre = pre, summary_df = summary_df, contrasts_df = inputs$contrasts, row_data = pre$row_data, config = config)
        files <- c(files, save_tsv(final_results, dirs$datasets, "final_results.tsv"))

        # 4) Excel outputs
        # Extract excel_order from clustering_res if available
        excel_order <- if (!is.null(clustering_res)) clustering_res$excel_order else NULL

        if (exists("write_final_results_excels_legacy_rna")) {
            files <- c(files, write_final_results_excels_legacy_rna(final_results, pre, config, out_dir, excel_order = excel_order))
        } else if (exists("write_final_results_excels_legacy_generic")) {
            # Use generic if wrapper not available
            # Pass RAW counts (expr_filt) and clustering results
            # ID column will be "Gene" (from legacy rename in write_rnaseq_outputs_legacy)
            id_col <- if ("Gene" %in% names(final_results)) "Gene" else "FeatureID"

            # Extract Excel config for enriched layout (annotation rows, sample labels)
            excel_cfg <- config$modes$rna$excel %||% list()
            rna_sample_id_col <- config$modes$rna$effects$samples %||% "SampleID"

            files <- c(files, write_final_results_excels_legacy_generic(
                final_results = final_results,
                config = config,
                out_dir = out_dir,
                mode = "rna",
                id_col = id_col,
                expr_for_de = pre$expr_filt,  # RAW counts, not normalized
                with_cutoffs = TRUE,
                clustering_res = list(excel_order = excel_order),
                sample_meta = pre$meta,
                sample_id_col = rna_sample_id_col,
                annotation_rows = excel_cfg$annotation_rows,
                sample_label_cols = excel_cfg$sample_label_cols
            ))
        }
    }
    unique(files)
}


build_final_results_rnaseq <- function(pre, summary_df, contrasts_df, row_data = NULL,
                                        config = NULL) {
    # FIX 1: Remove annotation columns (gene_name, symbol, description) from final_results
    # FIX 2: Pass mode="rna" for correct column naming (no ".imputs.")
    # FIX 5: Use RAW counts (expr_filt) not normalized (expr_work)
    # FIX 6: Auto-detect ID column (may be "Gene" after legacy rename, or "FeatureID" from build_rnaseq_summary_df)
    id_col <- if ("Gene" %in% names(summary_df)) "Gene" else "FeatureID"

    # Per-group CV on linear CPM. We deliberately use library-size-corrected
    # CPM (2^expr_work would be wrong for VST/preprocessed runs; raw-count CV
    # is confounded by library size), computed from the raw counts directly.
    cv_cols <- build_group_cv_rnaseq(pre, contrasts_df, config)

    build_final_results_generic(
        summary_df = summary_df,
        expr_df = pre$expr_filt,  # FIX 5: RAW counts for final results tables
        contrasts_df = contrasts_df,
        feature_id_col = id_col,  # FIX 6: Auto-detected ID column
        annot_cols = NULL,  # FIX 1: No annotation columns in RNA final_results
        row_data = row_data %||% pre$row_data,
        fc_is_signed = TRUE,  # log2FoldChange is signed
        mode = "rna",  # FIX 2: Use RNA column naming
        cv_cols = cv_cols
    )
}

#' Build per-group CV columns for RNA-seq final results
#'
#' CV is computed on linear CPM (\code{compute_cpm} of the filtered raw counts),
#' which is library-size-corrected and method-agnostic. Returns \code{NULL} when
#' the feature is disabled (\code{modes$rna$excel$group_cv = false}), inputs are
#' unavailable, or the input is not genuine raw counts (see the source-type
#' whitelist below) — since \code{compute_cpm()} divides by column sums assuming
#' counts, running it on already-processed values would yield meaningless CPM.
#'
#' @param pre RNA-seq preprocessing results (uses \code{expr_filt}, \code{meta},
#'   \code{info$source_type}).
#' @param contrasts_df Contrasts table (Factor, Numerator, Denominator).
#' @param config Full pipeline config (for the feature flag + sample-ID column).
#' @return Feature-indexed data.frame of \code{CV.<group>} columns, or NULL.
build_group_cv_rnaseq <- function(pre, contrasts_df, config = NULL) {
    if (is.null(config)) return(NULL)
    enabled <- config$modes$rna$excel$group_cv %||% TRUE
    if (!isTRUE(enabled)) return(NULL)
    if (is.null(pre$expr_filt) || is.null(pre$meta)) return(NULL)

    # WHITELIST: compute CV only when expr_filt is genuine raw counts, so
    # compute_cpm() (counts / colSums) is valid. Anything else -- "preprocessed"
    # (expr_filt holds already-processed values, not counts) or any future
    # non-count source_type -- skips CV + warns. Fail-safe by construction: a new
    # source_type can't silently produce garbage CPM until it is added here.
    raw_count_sources <- c("matrix", "tximport")
    source_type <- pre$info$source_type %||% attr(pre, "source_type") %||% NA_character_
    if (!isTRUE(source_type %in% raw_count_sources)) {
        warning(sprintf(
            "RNA-seq group CV: source_type='%s' is not a known raw-count type (%s); skipping CV (compute_cpm() is only valid on raw counts).",
            source_type, paste(raw_count_sources, collapse = "/")
        ))
        return(NULL)
    }

    sample_id_col <- config$modes$rna$effects$samples %||% "SampleID"
    cpm_linear <- compute_cpm(pre$expr_filt)
    compute_group_cv_columns(
        expr_linear   = cpm_linear,
        sample_meta   = pre$meta,
        sample_id_col = sample_id_col,
        contrasts_df  = contrasts_df
    )
}
