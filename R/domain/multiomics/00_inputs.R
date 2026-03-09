#' Load preprocessed data from individual omics pipelines
#'
#' Expects that rna_pre, prot_pre, and/or metab_pre targets have already run.
#' When input_mode == "outputs", loads from shiny payload RDS files instead.
#' Returns a named list with available omics data.
#'
#' @param config Full config object
#' @param rna_pre Preprocessed RNA-seq data (from pipe_rnaseq)
#' @param prot_pre Preprocessed proteomics data (from pipe_proteomics)
#' @param metab_pre Preprocessed metabolomics data (from pipe_metabolomics)
#' @return Named list with: omics_available, rna, proteomics, metabolomics
load_multiomics_inputs <- function(config, rna_pre = NULL, prot_pre = NULL, metab_pre = NULL) {

    multi_cfg <- config$modes$multiomics

    # If input_mode == "outputs", load from shiny payload RDS files
    if (identical(multi_cfg$input_mode, "outputs")) {
        return(load_multiomics_from_payloads(config))
    }

    omics_available <- character(0)
    result <- list()

    # Collect available omics
    if (!is.null(rna_pre)) {
        omics_available <- c(omics_available, "transcriptomics")
        result$transcriptomics <- rna_pre
        message(sprintf("RNA-seq data loaded: %d features x %d samples",
                        nrow(rna_pre$expr_work), ncol(rna_pre$expr_work)))
    }

    if (!is.null(prot_pre)) {
        omics_available <- c(omics_available, "proteomics")
        result$proteomics <- prot_pre
        message(sprintf("Proteomics data loaded: %d features x %d samples",
                        nrow(prot_pre$expr_work), ncol(prot_pre$expr_work)))
    }

    if (!is.null(metab_pre)) {
        omics_available <- c(omics_available, "metabolomics")
        result$metabolomics <- metab_pre
        message(sprintf("Metabolomics data loaded: %d features x %d samples",
                        nrow(metab_pre$expr_work), ncol(metab_pre$expr_work)))
    }

    if (length(omics_available) < 2) {
        stop(
            "Multi-omics integration requires at least 2 omics layers. ",
            "Found: ", paste(omics_available, collapse = ", "),
            call. = FALSE
        )
    }

    result$omics_available <- omics_available
    result$n_omics <- length(omics_available)

    message(sprintf("Multi-omics integration ready: %s",
                    paste(omics_available, collapse = " + ")))

    result
}


#' Load multi-omics data from pre-computed shiny payload RDS files
#'
#' Reads shiny payloads, applies sample mapping to unify sample IDs,
#' and converts to the _pre-like format expected by harmonization.
#'
#' @param config Full config object with modes$multiomics$payload_paths
#' @return Named list compatible with load_multiomics_inputs() output
load_multiomics_from_payloads <- function(config) {

    multi_cfg <- config$modes$multiomics
    payload_paths <- multi_cfg$payload_paths
    condition_col <- multi_cfg$condition_column %||%
                     config$design$condition_column %||% "condition"

    if (is.null(payload_paths) || length(payload_paths) == 0) {
        stop("input_mode is 'outputs' but no payload_paths specified in config",
             call. = FALSE)
    }

    # Load sample mapping if provided
    sample_map <- NULL
    mapping_path <- multi_cfg$sample_mapping
    if (!is.null(mapping_path) && file.exists(mapping_path)) {
        sample_map <- read.csv(mapping_path, stringsAsFactors = FALSE)
        message(sprintf("Sample mapping loaded: %d rows from %s",
                        nrow(sample_map), basename(mapping_path)))
    }

    # Map config keys to internal omics names
    key_to_omics <- c(rna = "transcriptomics", proteomics = "proteomics",
                      metabolomics = "metabolomics")

    omics_available <- character(0)
    result <- list()

    for (key in names(payload_paths)) {
        path <- payload_paths[[key]]
        if (is.null(path) || !file.exists(path)) {
            message(sprintf("Payload not found for '%s': %s — skipping", key, path %||% "NULL"))
            next
        }

        payload <- readRDS(path)
        omics_name <- key_to_omics[[key]] %||% key

        # Convert payload to _pre-like format
        pre <- payload_to_pre(payload, key, sample_map, condition_col)
        if (is.null(pre)) next

        omics_available <- c(omics_available, omics_name)
        result[[omics_name]] <- pre
        message(sprintf("%s loaded from payload: %d features x %d samples",
                        omics_name, nrow(pre$expr_work), ncol(pre$expr_work)))
    }

    if (length(omics_available) < 2) {
        stop(
            "Multi-omics integration requires at least 2 omics layers. ",
            "Found: ", paste(omics_available, collapse = ", "),
            call. = FALSE
        )
    }

    result$omics_available <- omics_available
    result$n_omics <- length(omics_available)

    message(sprintf("Multi-omics integration ready (from payloads): %s",
                    paste(omics_available, collapse = " + ")))

    result
}


#' Convert a shiny payload to _pre-like format
#'
#' @param payload Shiny payload list (from readRDS)
#' @param omics_key Config key (rna, proteomics, metabolomics)
#' @param sample_map Optional sample mapping data.frame
#' @param condition_col Condition column name
#' @return List with expr_work, meta, row_data, de_stats
payload_to_pre <- function(payload, omics_key, sample_map = NULL,
                           condition_col = "condition") {

    expr <- payload$expr_norm
    if (is.null(expr)) {
        warning("No expr_norm in payload for ", omics_key)
        return(NULL)
    }

    # Build metadata
    meta <- payload$sample_meta
    if (is.null(meta)) {
        meta <- data.frame(SampleID = colnames(expr), stringsAsFactors = FALSE)
    }

    # Apply sample mapping: rename columns and metadata to unified IDs
    if (!is.null(sample_map) && omics_key %in% colnames(sample_map)) {
        orig_ids <- sample_map[[omics_key]]
        unified_ids <- sample_map$unified_id

        # Match expression columns to mapping
        col_match <- match(colnames(expr), orig_ids)
        matched <- !is.na(col_match)
        if (sum(matched) > 0) {
            expr <- expr[, matched, drop = FALSE]
            colnames(expr) <- unified_ids[col_match[matched]]

            # Rebuild metadata with unified IDs and condition from mapping
            meta <- data.frame(
                SampleID = unified_ids[col_match[matched]],
                stringsAsFactors = FALSE
            )
            if ("condition" %in% colnames(sample_map)) {
                meta[[condition_col]] <- sample_map$condition[col_match[matched]]
            }
            rownames(meta) <- meta$SampleID
            message(sprintf("  %s: %d/%d samples matched via mapping",
                            omics_key, sum(matched), length(orig_ids)))
        }
    }

    # Build row_data from feature_annot if available
    row_data <- payload$feature_annot
    if (is.null(row_data)) {
        row_data <- data.frame(feature_id = rownames(expr),
                               stringsAsFactors = FALSE)
        rownames(row_data) <- rownames(expr)
    }

    # Extract DE stats if available
    de_stats <- payload$de_stats

    list(
        expr_work = expr,
        meta = meta,
        row_data = row_data,
        de_stats = de_stats
    )
}


#' Validate multi-omics inputs before harmonization
#'
#' Checks that all omics have compatible sample metadata and sufficient overlap.
#'
#' @param inputs Output from load_multiomics_inputs()
#' @param config Full config object
#' @return Invisibly returns TRUE if valid, stops otherwise
validate_multiomics_inputs <- function(inputs, config) {

    omics <- inputs$omics_available

    # Extract sample IDs from each omics
    sample_lists <- list()
    for (om in omics) {
        meta <- inputs[[om]]$meta

        # Map internal omics names to config keys
        om_config_key <- switch(om,
                                "transcriptomics" = "rna",
                                "proteomics" = "proteomics",
                                "metabolomics" = "metabolomics",
                                om)  # fallback

        sample_col <- config$modes[[om_config_key]]$effects$samples %||%
                      config$modes[[om_config_key]]$id_columns$sample_col %||%
                      "SampleID"

        if (!sample_col %in% colnames(meta)) {
            stop(sprintf("Sample column '%s' not found in %s metadata",
                         sample_col, om))
        }

        sample_lists[[om]] <- as.character(meta[[sample_col]])
    }

    # Check sample overlap
    common_samples <- Reduce(intersect, sample_lists)

    if (length(common_samples) < 3) {
        stop(
            "Insufficient sample overlap across omics layers. ",
            "Common samples: ", length(common_samples), " (need ≥3)\n",
            paste(
                sprintf("  %s: %d samples", names(sample_lists),
                        sapply(sample_lists, length)),
                collapse = "\n"
            ),
            call. = FALSE
        )
    }

    message(sprintf(
        "Sample validation passed: %d common samples across %d omics layers",
        length(common_samples), length(omics)
    ))

    invisible(TRUE)
}
