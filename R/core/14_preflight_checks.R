#' Pre-flight Data Validation Checks
#'
#' Validates data files before running pipelines.
#' Catches data format issues, sample mismatches, and potential problems early.
#'
#' Uses read_table_auto() from 01_io.R and detect_id_type() /
#' detect_organism_from_ids() from 10_organism_detection.R.

# ==============================================================================
# MAIN PREFLIGHT CHECK FUNCTION
# ==============================================================================

#' Run all pre-flight checks for a pipeline
#'
#' @param config Configuration list (full config from load_config)
#' @param pipeline_type One of "rnaseq", "proteomics"
#' @param verbose Print progress
#' @return List with check results
#' @export
run_preflight_checks <- function(config, pipeline_type, verbose = TRUE) {
    result <- list(
        passed = TRUE,
        errors = character(0),
        warnings = character(0),
        info = character(0),
        stats = list()
    )

    if (verbose) message("\n=== Running Pre-flight Checks ===\n")

    result <- switch(pipeline_type,
        "rnaseq" = preflight_rnaseq(config, result, verbose),
        result
    )

    result$passed <- length(result$errors) == 0

    if (verbose) {
        message("\n=== Pre-flight Check Summary ===")

        if (length(result$errors) > 0) {
            message("\nERRORS (must fix):")
            for (e in result$errors) message(sprintf("  [x] %s", e))
        }

        if (length(result$warnings) > 0) {
            message("\nWARNINGS:")
            for (w in result$warnings) message(sprintf("  [!] %s", w))
        }

        if (length(result$info) > 0) {
            message("\nINFO:")
            for (i in result$info) message(sprintf("  [i] %s", i))
        }

        if (result$passed) {
            message("\n  All pre-flight checks passed")
        } else {
            message("\n  Pre-flight checks FAILED - please fix errors before running")
        }
    }

    result
}

# ==============================================================================
# RNASEQ PREFLIGHT
# ==============================================================================

#' @noRd
preflight_rnaseq <- function(config, result, verbose) {
    # Navigate multiomics-core nested config
    rna_cfg <- config$modes$rna
    base_path <- file.path(config$project$dir, config$paths$raw)

    # --- Counts file ---
    counts <- NULL
    counts_rel <- rna_cfg$files$counts
    if (!is.null(counts_rel) && nzchar(counts_rel)) {
        counts_file <- file.path(base_path, counts_rel)
        if (!file.exists(counts_file)) {
            result$errors <- c(result$errors,
                sprintf("Counts file not found: %s", counts_file))
        } else {
            if (verbose) message("Checking counts file...")
            counts <- tryCatch(read_table_auto(counts_file), error = function(e) {
                result$errors <<- c(result$errors, sprintf("Failed to read counts: %s", e$message))
                NULL
            })

            if (!is.null(counts)) {
                result <- check_expression_matrix(counts, "counts", result, verbose)
                result$stats$n_genes <- nrow(counts)
                result$stats$n_samples_counts <- ncol(counts) - 1L
            }
        }
    } else {
        # Preprocessed data may be provided instead of raw counts
        has_preprocessed <- !is.null(rna_cfg$files$preprocessed_counts) &&
                            nzchar(rna_cfg$files$preprocessed_counts %||% "")
        if (has_preprocessed) {
            result$info <- c(result$info, "No raw counts file; using preprocessed_counts instead")
        } else {
            result$errors <- c(result$errors, "No counts file specified in modes.rna.files.counts")
        }
    }

    # --- Metadata file ---
    meta_rel <- rna_cfg$files$metadata
    metadata <- NULL
    if (!is.null(meta_rel) && nzchar(meta_rel)) {
        meta_file <- file.path(base_path, meta_rel)
        if (!file.exists(meta_file)) {
            result$errors <- c(result$errors,
                sprintf("Metadata file not found: %s", meta_file))
        } else {
            if (verbose) message("Checking metadata file...")
            metadata <- tryCatch(read_table_auto(meta_file), error = function(e) {
                result$errors <<- c(result$errors, sprintf("Failed to read metadata: %s", e$message))
                NULL
            })

            if (!is.null(metadata)) {
                result <- check_metadata(metadata, rna_cfg, result, verbose)
                result$stats$n_samples_metadata <- nrow(metadata)

                if (!is.null(counts)) {
                    result <- check_sample_matching(counts, metadata, rna_cfg, result, verbose)
                }

                result <- check_group_sizes(metadata, rna_cfg, result, verbose)
            }
        }
    } else {
        result$errors <- c(result$errors, "No metadata file specified in modes.rna.files.metadata")
    }

    # --- Gene IDs ---
    if (!is.null(counts)) {
        gene_ids <- .pf_get_feature_ids(counts)
        result <- check_gene_ids(gene_ids, rna_cfg, result, verbose)
    }

    result
}

# ==============================================================================
# INDIVIDUAL CHECK FUNCTIONS
# ==============================================================================

#' Check expression/intensity matrix format
#' @noRd
check_expression_matrix <- function(data, name, result, verbose) {
    if (verbose) message(sprintf("  Matrix dimensions: %d rows x %d columns", nrow(data), ncol(data)))

    numeric_cols <- sapply(data, is.numeric)
    n_numeric <- sum(numeric_cols)

    if (n_numeric == 0) {
        result$errors <- c(result$errors, sprintf("%s matrix has no numeric columns", name))
    } else if (n_numeric < ncol(data) - 1) {
        result$info <- c(result$info,
            sprintf("%s matrix: %d numeric columns (samples), %d non-numeric (ID columns)",
                    name, n_numeric, ncol(data) - n_numeric))
    }

    if (any(numeric_cols)) {
        numeric_data <- data[, numeric_cols, drop = FALSE]
        all_na_rows <- apply(numeric_data, 1, function(x) all(is.na(x)))
        if (sum(all_na_rows) > 0) {
            result$warnings <- c(result$warnings,
                sprintf("%d rows with all missing values in %s matrix", sum(all_na_rows), name))
        }
        all_na_cols <- apply(numeric_data, 2, function(x) all(is.na(x)))
        if (sum(all_na_cols) > 0) {
            result$errors <- c(result$errors,
                sprintf("%d columns with all missing values in %s matrix", sum(all_na_cols), name))
        }
    }

    result
}

#' Check metadata format
#' @noRd
check_metadata <- function(metadata, rna_cfg, result, verbose) {
    if (verbose) message(sprintf("  Metadata: %d samples, %d columns", nrow(metadata), ncol(metadata)))

    sample_col <- rna_cfg$id_columns$sample_col %||%
                  rna_cfg$effects$samples %||%
                  "SampleID"

    if (!(sample_col %in% names(metadata))) {
        potential <- c("sample_id", "SampleID", "Sample_ID", "sample", "ID")
        found <- intersect(potential, names(metadata))
        if (length(found) > 0) {
            result$info <- c(result$info,
                sprintf("Sample column '%s' not found, but found: %s", sample_col, paste(found, collapse = ", ")))
        } else {
            result$warnings <- c(result$warnings,
                sprintf("Sample column '%s' not found in metadata", sample_col))
        }
    }

    group_col <- rna_cfg$filtering$group_col %||% "condition"

    if (!(group_col %in% names(metadata))) {
        potential <- c("condition", "Condition", "group", "Group", "treatment", "Treatment")
        found <- intersect(potential, names(metadata))
        if (length(found) > 0) {
            result$info <- c(result$info,
                sprintf("Group column '%s' not found, but found: %s", group_col, paste(found, collapse = ", ")))
        } else {
            result$warnings <- c(result$warnings,
                sprintf("Group column '%s' not found in metadata", group_col))
        }
    } else {
        groups <- unique(metadata[[group_col]])
        result$info <- c(result$info,
            sprintf("Groups found: %s (%d levels)", paste(groups, collapse = ", "), length(groups)))
    }

    result
}

#' Check sample ID matching between data and metadata
#' @noRd
check_sample_matching <- function(data, metadata, rna_cfg, result, verbose) {
    data_samples <- names(data)
    exclude_patterns <- c("^gene", "^feature", "^protein", "^id$", "^name$", "^symbol$", "^description$", "^entrez")
    is_sample <- !grepl(paste(exclude_patterns, collapse = "|"), data_samples, ignore.case = TRUE)
    data_samples <- data_samples[is_sample]

    sample_col <- rna_cfg$id_columns$sample_col %||%
                  rna_cfg$effects$samples %||%
                  "SampleID"

    if (!(sample_col %in% names(metadata))) {
        sample_col <- names(metadata)[1]
    }

    metadata_samples <- as.character(metadata[[sample_col]])

    common <- intersect(data_samples, metadata_samples)
    only_data <- setdiff(data_samples, metadata_samples)
    only_meta <- setdiff(metadata_samples, data_samples)

    if (verbose) {
        message(sprintf("  Sample matching: %d in both, %d only in data, %d only in metadata",
                        length(common), length(only_data), length(only_meta)))
    }

    if (length(common) == 0) {
        result$errors <- c(result$errors,
            "No matching sample IDs between data and metadata - check ID formats")
        if (length(data_samples) > 0 && length(metadata_samples) > 0) {
            result$info <- c(result$info,
                sprintf("Data sample IDs look like: %s", paste(head(data_samples, 3), collapse = ", ")))
            result$info <- c(result$info,
                sprintf("Metadata sample IDs look like: %s", paste(head(metadata_samples, 3), collapse = ", ")))
        }
    } else {
        if (length(only_data) > 0) {
            result$warnings <- c(result$warnings,
                sprintf("%d samples in data but not in metadata (will be excluded)", length(only_data)))
        }
        if (length(only_meta) > 0) {
            result$info <- c(result$info,
                sprintf("%d samples in metadata but not in data", length(only_meta)))
        }
    }

    result$stats$n_matched_samples <- length(common)
    result
}

#' Check group sizes for statistical validity
#' @noRd
check_group_sizes <- function(metadata, rna_cfg, result, verbose) {
    group_col <- rna_cfg$filtering$group_col %||% "condition"

    if (!(group_col %in% names(metadata))) return(result)

    group_sizes <- table(metadata[[group_col]])
    min_size <- min(group_sizes)

    if (verbose) {
        message("  Group sizes:")
        for (g in names(group_sizes)) {
            message(sprintf("    %s: %d", g, group_sizes[g]))
        }
    }

    critical <- names(group_sizes)[group_sizes < 2]
    if (length(critical) > 0) {
        result$errors <- c(result$errors,
            sprintf("Groups with n < 2: %s. Cannot perform statistical tests!",
                    paste(critical, collapse = ", ")))
    }

    severe <- names(group_sizes)[group_sizes >= 2 & group_sizes < 3]
    if (length(severe) > 0) {
        result$warnings <- c(result$warnings,
            sprintf("Groups with n = 2: %s. Results will be unreliable.",
                    paste(severe, collapse = ", ")))
    }

    low_power <- names(group_sizes)[group_sizes >= 3 & group_sizes < 5]
    if (length(low_power) > 0) {
        result$warnings <- c(result$warnings,
            sprintf("Groups with n < 5: %s. Statistical power is limited.",
                    paste(low_power, collapse = ", ")))
    }

    if (length(group_sizes) >= 2) {
        ratio <- max(group_sizes) / min(group_sizes)
        if (ratio > 5) {
            result$warnings <- c(result$warnings,
                sprintf("Highly unbalanced groups (ratio %.1f:1).", ratio))
        }
    }

    result$stats$group_sizes <- as.list(group_sizes)
    result
}

#' Check gene ID format and type
#' @noRd
check_gene_ids <- function(gene_ids, rna_cfg, result, verbose) {
    # Use detect_id_type / detect_organism_from_ids from 10_organism_detection.R
    tryCatch({
        id_type <- detect_id_type(gene_ids)
        if (verbose) message(sprintf("  Detected ID type: %s", id_type))
        result$stats$detected_id_type <- id_type

        if (has_ensembl_version(gene_ids)) {
            result$info <- c(result$info,
                "Ensembl IDs have version suffixes - will be stripped during annotation")
        }

        org_result <- detect_organism_from_ids(gene_ids)
        if (org_result$organism != "unknown") {
            if (verbose) message(sprintf("  Detected organism: %s (confidence: %s)",
                                          org_result$organism, org_result$confidence))
            result$stats$detected_organism <- org_result$organism

            config_organism <- rna_cfg$annotation$organism
            if (!is.null(config_organism) && config_organism != "auto" &&
                tolower(org_result$organism) != tolower(config_organism) &&
                org_result$confidence == "high") {
                result$warnings <- c(result$warnings,
                    sprintf("Detected organism (%s) differs from config (%s)",
                            org_result$organism, config_organism))
            }
        }
    }, error = function(e) {
        # Silently skip if organism detection unavailable
    })

    n_dup <- sum(duplicated(gene_ids))
    if (n_dup > 0) {
        result$info <- c(result$info, sprintf("%d duplicate feature IDs", n_dup))
    }

    result
}

# ==============================================================================
# INTERNAL HELPERS
# ==============================================================================

#' Get feature IDs from a data frame
#' @noRd
.pf_get_feature_ids <- function(data) {
    first_col <- data[[1]]
    if (is.character(first_col) || is.factor(first_col)) {
        return(as.character(first_col))
    }
    if (!is.null(rownames(data)) && !all(rownames(data) == as.character(seq_len(nrow(data))))) {
        return(rownames(data))
    }
    as.character(first_col)
}
