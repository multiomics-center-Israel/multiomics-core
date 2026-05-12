#' Module: Multi-omics data harmonization
#'
#' Loads preprocessed data from individual omics pipelines,
#' builds gene-protein mapping, and constructs a MultiAssayExperiment.
#'
#' @param config Full config object
#' @param rna_pre Preprocessed RNA-seq data (can be NULL)
#' @param prot_pre Preprocessed proteomics data (can be NULL)
#' @param metab_pre Preprocessed metabolomics data (can be NULL)
#' @param out_dir Output directory
#' @return List with: mae, gene_protein_mapping, inputs, validation
mod_multiomics_harmonization <- function(config, rna_pre = NULL, prot_pre = NULL,
                                          metab_pre = NULL, out_dir) {

    message("\n=== Multi-Omics Harmonization ===\n")

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Load multi-omics inputs
    inputs <- load_multiomics_inputs(
        config = config,
        rna_pre = rna_pre,
        prot_pre = prot_pre,
        metab_pre = metab_pre
    )

    # Validate inputs
    validate_multiomics_inputs(inputs, config)

    # Build gene-protein mapping (if RNA + protein are both present)
    gene_protein_mapping <- NULL
    if ("transcriptomics" %in% inputs$omics_available &&
        "proteomics" %in% inputs$omics_available) {

        gene_protein_mapping <- build_gene_protein_mapping(
            rna_data = inputs$transcriptomics,
            prot_data = inputs$proteomics,
            config = config
        )

        # Optionally filter to 1:1 mapping
        if (isTRUE(config$modes$multiomics$require_one_to_one_mapping)) {
            gene_protein_mapping <- filter_to_one_to_one_mapping(gene_protein_mapping)
        }

        # Save mapping
        if (nrow(gene_protein_mapping) > 0) {
            write.csv(gene_protein_mapping,
                      file.path(out_dir, "gene_protein_mapping.csv"),
                      row.names = FALSE)
        }
    }

    # Build MultiAssayExperiment
    mae <- build_mae(
        inputs = inputs,
        config = config,
        gene_protein_mapping = gene_protein_mapping
    )

    # Optionally harmonize features via mapping
    if (!is.null(gene_protein_mapping) && nrow(gene_protein_mapping) > 0) {
        harmonization_res <- harmonize_features_via_mapping(mae, gene_protein_mapping)
        mae <- harmonization_res$mae_harmonized

        message(sprintf(
            "Feature harmonization: %d gene-protein pairs aligned",
            harmonization_res$feature_overlap
        ))
    }

    message(sprintf(
        "MAE constructed: %d samples, %s",
        ncol(mae),
        paste(
            sprintf("%s=%d", names(mae@ExperimentList),
                    sapply(mae@ExperimentList, nrow)),
            collapse = ", "
        )
    ))

    list(
        mae = mae,
        gene_protein_mapping = gene_protein_mapping,
        inputs = inputs,
        omics_available = inputs$omics_available
    )
}
