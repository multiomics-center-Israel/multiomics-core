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

    # Get contrasts (already loaded by load_rna_inputs)
    contrasts_df <- inputs$contrasts
    if (is.null(contrasts_df)) {
        stop("Contrasts not found in inputs. Check config$modes$rna$files$contrasts")
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
