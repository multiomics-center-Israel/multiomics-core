pipe_rnaseq <- function() {
    list(
        tar_target(rna_inputs, load_rna_inputs(config)),
        tar_target(rna_pre, preprocess_rna(rna_inputs, config, gene_lengths = NULL, verbose = TRUE)),
        tar_target(rna_out_dir, get_mode_out_dir(run_dir, "rna")),
        tar_target(
            rna_qc_pre_obj,
            mod_rnaseq_qc_pre(
                pre     = rna_pre,
                config  = config,
                out_dir = rna_out_dir
            )
        ),
        tar_target(rna_de_res, mod_rnaseq_de(rna_pre, rna_inputs, config, verbose = TRUE)),
        # Clustering (run before outputs to provide excel_order)
        tar_target(
            rna_clustering_obj,
            mod_rnaseq_clustering(
                pre = rna_pre,
                de_res = rna_de_res,
                config = config,
                out_dir = rna_out_dir
            )
        ),
        # Legacy outputs (TSV files) - now receives clustering_res
        tar_target(
            rna_outputs_legacy,
            write_rnaseq_outputs_legacy(
                pre = rna_pre,
                de_res = rna_de_res,
                inputs = rna_inputs,
                config = config,
                out_dir = file.path(run_dir, "rna"),
                clustering_res = rna_clustering_obj
            ),
            format = "file"
        ),
        tar_target(
            rna_qc_post_obj,
            mod_rnaseq_qc_post(
                pre     = rna_pre,
                de_res  = rna_de_res,
                config  = config,
                out_dir = rna_out_dir
            )
        ),
        tar_target(
            rna_shiny_legacy,
            save_data_to_shiny_legacy_rna(
                pre = rna_pre,
                de_res = rna_de_res,
                inputs = rna_inputs,
                config = config,
                pca_res = rna_qc_pre_obj,
                clustering_res = rna_clustering_obj,
                out_file = file.path(rna_out_dir, "data_to_shiny_legacy_rnaseq.rds")
            ),
            format = "file"
        )
    )
}
