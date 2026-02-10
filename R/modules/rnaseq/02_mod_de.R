# 03_mod_rnaseq_de
# mod_rnaseq_de
mod_rnaseq_de <- function(pre, inputs, config, verbose = FALSE) {
    assert_pre_contract(pre, stage = "rna")

    # Extract DE configuration
    de_cfg <- config$modes$rna$de

    # Get contrasts (already loaded by load_rna_inputs)
    contrasts_df <- inputs$contrasts
    if (is.null(contrasts_df)) {
        stop("Contrasts not found in inputs. Check config$modes$rna$files$contrasts")
    }

    # Delegate actual DE run to domain function
    run_deseq2_de(
        counts = pre$de_input,
        meta = pre$meta,
        contrasts_df = contrasts_df,
        de_cfg = de_cfg
    )
}
