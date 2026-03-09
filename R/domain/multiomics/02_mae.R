#' Build MultiAssayExperiment from preprocessed omics data
#'
#' Creates a unified MAE object with:
#' - Aligned expression matrices (features in rows, samples in columns)
#' - Sample metadata (colData)
#' - Feature metadata per assay (rowRanges/rowData)
#'
#' @param inputs Multi-omics inputs from load_multiomics_inputs()
#' @param config Full config object
#' @param gene_protein_mapping Optional gene-protein mapping for harmonization
#' @return MultiAssayExperiment object
build_mae <- function(inputs, config, gene_protein_mapping = NULL) {

    if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
        stop("Package 'MultiAssayExperiment' is required for multi-omics integration")
    }
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("Package 'SummarizedExperiment' is required")
    }

    omics <- inputs$omics_available

    # Extract sample metadata from each omics
    sample_maps <- list()
    expr_lists <- list()
    row_data_lists <- list()

    for (om in omics) {
        om_data <- inputs[[om]]

        # Map internal omics names to config keys
        om_config_key <- switch(om,
                                "transcriptomics" = "rna",
                                "proteomics" = "proteomics",
                                "metabolomics" = "metabolomics",
                                om)  # fallback
        om_cfg <- config$modes[[om_config_key]]

        sample_col <- om_cfg$effects$samples %||%
                      om_cfg$id_columns$sample_col %||%
                      "SampleID"

        # Expression matrix (features x samples)
        expr <- om_data$expr_work
        if (is.null(expr)) {
            stop(sprintf("Missing expr_work for %s", om))
        }

        # Sample IDs
        sample_ids <- colnames(expr)

        # Sample metadata
        meta <- om_data$meta
        if (!sample_col %in% colnames(meta)) {
            stop(sprintf("Sample column '%s' not found in %s metadata", sample_col, om))
        }

        # Align metadata to expression columns
        meta <- meta[match(sample_ids, meta[[sample_col]]), , drop = FALSE]
        rownames(meta) <- sample_ids

        sample_maps[[om]] <- meta
        expr_lists[[om]] <- expr
        row_data_lists[[om]] <- om_data$row_data
    }

    # Find common samples across all omics
    all_samples <- lapply(sample_maps, rownames)
    common_samples <- Reduce(intersect, all_samples)

    if (length(common_samples) < 3) {
        stop(
            "Insufficient common samples for MAE construction: ",
            length(common_samples), " (need ≥3)"
        )
    }

    message(sprintf(
        "Building MAE with %d common samples across %d omics layers",
        length(common_samples), length(omics)
    ))

    # Subset expression matrices to common samples
    for (om in omics) {
        expr_lists[[om]] <- expr_lists[[om]][, common_samples, drop = FALSE]
        sample_maps[[om]] <- sample_maps[[om]][common_samples, , drop = FALSE]
    }

    # Merge sample metadata (use first omics as base, add omics-specific columns with prefixes)
    coldata <- sample_maps[[omics[1]]]

    for (i in seq_along(omics)[-1]) {
        om <- omics[i]
        om_meta <- sample_maps[[om]]

        # Add omics-specific columns with prefix
        om_specific_cols <- setdiff(colnames(om_meta), colnames(coldata))
        if (length(om_specific_cols) > 0) {
            om_meta_subset <- om_meta[, om_specific_cols, drop = FALSE]
            colnames(om_meta_subset) <- paste0(om, "_", colnames(om_meta_subset))
            coldata <- cbind(coldata, om_meta_subset)
        }
    }

    # Build SummarizedExperiment objects for each omics
    se_list <- list()
    for (om in omics) {
        expr <- expr_lists[[om]]
        rd <- row_data_lists[[om]]

        # Ensure row_data has feature IDs
        if (is.null(rd)) {
            rd <- data.frame(
                feature_id = rownames(expr),
                row.names = rownames(expr),
                stringsAsFactors = FALSE
            )
        } else {
            if (!identical(rownames(rd), rownames(expr))) {
                # Align row_data to expression
                rd <- rd[match(rownames(expr), rownames(rd)), , drop = FALSE]
                rownames(rd) <- rownames(expr)
            }
        }

        se <- SummarizedExperiment::SummarizedExperiment(
            assays = list(expr = expr),
            rowData = rd,
            colData = coldata[colnames(expr), , drop = FALSE]
        )

        se_list[[om]] <- se
    }

    # Build MAE
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = se_list,
        colData = coldata
    )

    message(sprintf(
        "MAE construction complete:\n  Assays: %s\n  Samples: %d\n  Features: %s",
        paste(names(se_list), collapse = ", "),
        nrow(coldata),
        paste(
            sprintf("%s=%d", names(se_list), sapply(se_list, nrow)),
            collapse = ", "
        )
    ))

    mae
}


#' Extract common features across omics using gene-protein mapping
#'
#' Creates a harmonized feature set where RNA and protein data are matched
#' by gene symbols or mapping table.
#'
#' @param mae MultiAssayExperiment object
#' @param gene_protein_mapping Mapping data frame (gene_id, protein_id)
#' @return List with: mae_harmonized, mapping, feature_overlap
harmonize_features_via_mapping <- function(mae, gene_protein_mapping) {

    if (is.null(gene_protein_mapping) || nrow(gene_protein_mapping) == 0) {
        message("No gene-protein mapping available; returning unharmonized MAE")
        return(list(mae_harmonized = mae, mapping = NULL, feature_overlap = 0))
    }

    # Extract RNA and protein assays
    rna_assay <- mae[["transcriptomics"]]
    prot_assay <- mae[["proteomics"]]

    if (is.null(rna_assay) || is.null(prot_assay)) {
        message("Cannot harmonize: missing RNA or proteomics assay")
        return(list(mae_harmonized = mae, mapping = NULL, feature_overlap = 0))
    }

    rna_ids <- rownames(rna_assay)
    prot_ids <- rownames(prot_assay)

    # Filter mapping to available IDs
    map <- gene_protein_mapping[
        gene_protein_mapping$gene_id %in% rna_ids &
        gene_protein_mapping$protein_id %in% prot_ids,
    ]

    if (nrow(map) == 0) {
        message("No overlapping features found via mapping")
        return(list(mae_harmonized = mae, mapping = map, feature_overlap = 0))
    }

    message(sprintf(
        "Feature harmonization: %d gene-protein pairs",
        nrow(map)
    ))

    # Subset assays to mapped features
    rna_subset <- rna_assay[map$gene_id, ]
    prot_subset <- prot_assay[map$protein_id, ]

    # Create harmonized feature IDs (gene symbols if available)
    harmonized_ids <- paste0("GENE_", seq_len(nrow(map)))
    if ("symbol" %in% colnames(SummarizedExperiment::rowData(rna_subset))) {
        symbols <- SummarizedExperiment::rowData(rna_subset)$symbol
        harmonized_ids <- ifelse(!is.na(symbols) & symbols != "",
                                  symbols, harmonized_ids)
    }

    # Update rownames
    rownames(rna_subset) <- harmonized_ids
    rownames(prot_subset) <- harmonized_ids

    # Update MAE with harmonized assays
    mae_harmonized <- mae
    mae_harmonized[["transcriptomics"]] <- rna_subset
    mae_harmonized[["proteomics"]] <- prot_subset

    list(
        mae_harmonized = mae_harmonized,
        mapping = map,
        feature_overlap = nrow(map)
    )
}
