# 03_mod_rnaseq_de
# mod_rnaseq_de
mod_rnaseq_de <- function(pre, inputs, config, verbose = FALSE) {
    assert_pre_contract(pre, stage = "rna")

    # Delegate actual DE run to domain function
    run_deseq2_de(
        expr_filt = pre$expr_filt,
        meta = pre$meta,
        inputs = inputs,
        config = config
    )
}
