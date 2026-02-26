# 03_mod_rnaseq_de
# mod_rnaseq_de
mod_rnaseq_de <- function(pre, inputs, config, verbose = FALSE) {
    assert_pre_contract(pre, stage = "rna")

    # Check for pre-computed DE tables (skip limma/DESeq2 if provided)
    de_table_files <- config$modes$rna$files$de_table
    has_precomputed <- !is.null(de_table_files) && length(de_table_files) > 0 &&
                       all(nzchar(unlist(de_table_files)))

    if (has_precomputed) {
        message("RNA DE: loading pre-computed DE tables (skipping limma/DESeq2)")
        return(load_precomputed_rna_de(config, contrasts_df = inputs$contrasts))
    }

    # Extract DE configuration
    de_cfg <- config$modes$rna$de

    # Get contrasts (already loaded by load_rna_inputs), auto-generate if missing
    contrasts_df <- inputs$contrasts
    if (is.null(contrasts_df)) {
        group_col <- config$modes$rna$effects$color %||%
                     config$modes$rna$filtering$group_col %||% "Condition"
        contrasts_df <- auto_generate_contrasts(pre$meta, group_col)
    }

    # Dispatch DE method based on source type
    source_type <- attr(pre, "source_type") %||%
                   attr(pre$de_input, "source_type") %||%
                   "matrix"

    if (source_type == "preprocessed") {
        run_limma_rna_de(
            expr = pre$expr_work,
            meta = pre$meta,
            contrasts_df = contrasts_df,
            de_cfg = de_cfg
        )
    } else {
        run_deseq2_de(
            counts = pre$de_input,
            meta = pre$meta,
            contrasts_df = contrasts_df,
            de_cfg = de_cfg
        )
    }
}

#' Auto-generate all pairwise contrasts from metadata
#' @param meta    data.frame with sample metadata
#' @param group_col character, column name defining biological groups
#' @return data.frame with Contrast_name, Factor, Numerator, Denominator
auto_generate_contrasts <- function(meta, group_col) {
    lvls <- sort(unique(as.character(meta[[group_col]])))
    lvls <- lvls[!is.na(lvls) & nzchar(lvls)]
    if (length(lvls) < 2) stop("Cannot auto-generate contrasts: fewer than 2 groups in '", group_col, "'.")
    pairs <- combn(lvls, 2)
    df <- data.frame(
        Contrast_name = apply(pairs, 2, function(p) paste0(p[1], "_vs_", p[2])),
        Factor        = group_col,
        Numerator     = pairs[1, ],
        Denominator   = pairs[2, ],
        stringsAsFactors = FALSE
    )
    message(sprintf("Auto-generated %d pairwise contrast(s) from '%s': %s",
                    nrow(df), group_col, paste(df$Contrast_name, collapse = ", ")))
    df
}
