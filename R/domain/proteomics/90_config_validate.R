#' Validate Proteomics Configuration
validate_proteomics_config <- function(cfg) {
    assert_named_list(cfg, "modes$proteomics")

    # 1. Engine & Scale
    assert_scalar_chr(cfg$engine, "proteomics$engine", allow_null = TRUE)
    assert_one_of(cfg$scale_in, "proteomics$scale_in", c("linear", "log2"), allow_null = TRUE)

    if (!is.null(cfg$transform)) {
        assert_one_of(cfg$transform$target_scale, "proteomics$transform$target_scale", c("log2"), allow_null = TRUE)
    }

    # 2. ID Columns
    assert_named_list(cfg$id_columns, "proteomics$id_columns")
    assert_scalar_chr(cfg$id_columns$protein_id, "id_columns$protein_id")
    assert_scalar_chr(cfg$id_columns$sample_col, "id_columns$sample_col")

    # Ensure consistency between sample_col and map_to (if map_to is used)
    if (!is.null(cfg$id_columns$map_to)) {
        assert_scalar_chr(cfg$id_columns$map_to, "id_columns$map_to")
        if (!identical(cfg$id_columns$sample_col, cfg$id_columns$map_to)) {
            stop("id_columns$sample_col and id_columns$map_to must be identical if map_to is provided")
        }
    }

    # 3. Imputation Settings
    if (!is.null(cfg$imputation) && !is.null(cfg$imputation$method)) {
        assert_one_of(cfg$imputation$method, "imputation$method", c("none", "perseus", "dep2", "qrilc", "minval"))

        # Only validate repetition/consensus params when multi_imputation is not explicitly FALSE
        if (!identical(cfg$imputation$multi_imputation, FALSE)) {
            assert_scalar_num(cfg$imputation$no_repetitions, "imputation$no_repetitions", min_val = 1)
            assert_scalar_num(cfg$imputation$min_no_passed, "imputation$min_no_passed", min_val = 1)

            if (cfg$imputation$min_no_passed > cfg$imputation$no_repetitions) {
                stop("imputation$min_no_passed cannot be greater than imputation$no_repetitions")
            }
        }

        # Method-specific validations
        if (identical(cfg$imputation$method, "perseus")) {
            assert_scalar_num(cfg$imputation$width, "imputation$width", min_val = 0.0001)
            assert_scalar_num(cfg$imputation$downshift, "imputation$downshift", min_val = 0.0001)
        }

        if (identical(cfg$imputation$method, "dep2")) {
            assert_scalar_chr(cfg$imputation$dep2_method, "imputation$dep2_method")
            assert_scalar_num(cfg$imputation$dep2_random_seed, "imputation$dep2_random_seed")
        }

        if (identical(cfg$imputation$method, "qrilc")) {
            assert_scalar_num(cfg$imputation$qrilc_random_seed, "imputation$qrilc_random_seed", allow_null = TRUE)
        }
    }

    # 3b. Filtering Settings
    if (!is.null(cfg$filtering)) {
        if (!is.null(cfg$filtering$remove_contaminants)) {
            assert_scalar_bool(cfg$filtering$remove_contaminants, "filtering$remove_contaminants", allow_null = TRUE)
        }
        if (!is.null(cfg$filtering$contaminant_prefix)) {
            assert_scalar_chr(cfg$filtering$contaminant_prefix, "filtering$contaminant_prefix", allow_null = TRUE)
        }
        if (!is.null(cfg$filtering$min_groups)) {
            assert_scalar_num(cfg$filtering$min_groups, "filtering$min_groups", allow_null = TRUE, min_val = 1)
        }
    }

    # 4. de Settings
    if (!is.null(cfg$de)) {
        if (!is.null(cfg$de$method)) {
            assert_one_of(cfg$de$method, "de$method", c("limma", "ttest", "welch", "anova"), allow_null = TRUE)
        }
        assert_scalar_num(cfg$de$p_cutoff, "de$p_cutoff",
            allow_null = TRUE, min_val = .Machine$double.eps, max_val = 1
        )
        assert_scalar_num(cfg$de$linear_fc_cutoff, "de$linear_fc_cutoff", allow_null = TRUE, min_val = 1)
        if (!is.null(cfg$de$p_adjust_method)) {
            assert_one_of(cfg$de$p_adjust_method, "de$p_adjust_method",
                          c("BH", "bonferroni", "holm", "BY", "fdr", "none"), allow_null = TRUE)
        }
        if (!is.null(cfg$de$paired)) {
            assert_scalar_bool(cfg$de$paired, "de$paired", allow_null = TRUE)
        }
        if (!is.null(cfg$de$pairing_col)) {
            assert_scalar_chr(cfg$de$pairing_col, "de$pairing_col", allow_null = TRUE)
        }
    }

    # 4b. Pathway Settings
    if (!is.null(cfg$pathway)) {
        if (!is.null(cfg$pathway$gsea_ranking)) {
            assert_one_of(cfg$pathway$gsea_ranking, "pathway$gsea_ranking",
                          c("stat", "abs_lfc", "lfc"), allow_null = TRUE)
        }
        if (!is.null(cfg$pathway$pathway_volcano)) {
            assert_scalar_bool(cfg$pathway$pathway_volcano, "pathway$pathway_volcano", allow_null = TRUE)
        }
    }

    # 4c. Pathway databases validation
    if (!is.null(cfg$pathway$databases)) {
        valid_dbs <- c("all", "GO", "GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome")
        for (db in cfg$pathway$databases) {
            if (!db %in% valid_dbs) {
                warning("Unknown pathway database '", db, "'. Valid options: ",
                        paste(valid_dbs, collapse = ", "))
            }
        }
    }

    # 4d. Domain analysis validation
    if (!is.null(cfg$domain_analysis)) {
        assert_scalar_bool(cfg$domain_analysis$enabled, "domain_analysis$enabled", allow_null = TRUE)
    }

    # 5. PPI validation
    ppi <- cfg$ppi
    if (!is.null(ppi)) {
        assert_scalar_bool(ppi$enabled, "ppi$enabled", allow_null = TRUE)
        assert_scalar_num(ppi$significance_threshold, "ppi$significance_threshold",
                          allow_null = TRUE, min_val = 0, max_val = 1)
        assert_scalar_num(ppi$string_score_threshold, "ppi$string_score_threshold",
                          allow_null = TRUE, min_val = 0, max_val = 1000)
        assert_scalar_num(ppi$lfc_threshold, "ppi$lfc_threshold", allow_null = TRUE, min_val = 0)
        assert_scalar_bool(ppi$active_subnetwork, "ppi$active_subnetwork", allow_null = TRUE)
        assert_scalar_bool(ppi$complex_analysis, "ppi$complex_analysis", allow_null = TRUE)
    }

    # 6. Advanced stats validation
    adv <- cfg$advanced_stats
    if (!is.null(adv)) {
        assert_scalar_bool(adv$enabled, "advanced_stats$enabled", allow_null = TRUE)
        assert_scalar_num(adv$bootstrap_n, "advanced_stats$bootstrap_n",
                          allow_null = TRUE, min_val = 100)
        assert_scalar_num(adv$ci_level, "advanced_stats$ci_level",
                          allow_null = TRUE, min_val = 0.5, max_val = 0.999)
        assert_scalar_bool(adv$run_robust_regression, "advanced_stats$run_robust_regression", allow_null = TRUE)
        assert_scalar_bool(adv$run_power_analysis, "advanced_stats$run_power_analysis", allow_null = TRUE)
    }

    # 7. Clustering Settings (optional)
    if (!is.null(cfg$clustering)) {
        validate_clustering_config(cfg$clustering)
    }

    # 6. PPI Settings (optional)
    if (!is.null(cfg$ppi)) {
        assert_scalar_bool(cfg$ppi$enabled, "ppi$enabled", allow_null = TRUE)
        assert_scalar_num(cfg$ppi$string_score_threshold, "ppi$string_score_threshold",
                          allow_null = TRUE, min_val = 0, max_val = 1000)
        assert_scalar_num(cfg$ppi$significance_threshold, "ppi$significance_threshold",
                          allow_null = TRUE, min_val = 0, max_val = 1)
        assert_scalar_num(cfg$ppi$lfc_threshold, "ppi$lfc_threshold",
                          allow_null = TRUE, min_val = 0)
    }

    # 7. Advanced Stats Settings (optional)
    if (!is.null(cfg$advanced_stats)) {
        assert_scalar_bool(cfg$advanced_stats$enabled, "advanced_stats$enabled", allow_null = TRUE)
        assert_scalar_num(cfg$advanced_stats$bootstrap_n, "advanced_stats$bootstrap_n",
                          allow_null = TRUE, min_val = 100)
        assert_scalar_bool(cfg$advanced_stats$compute_effect_size_ci,
                           "advanced_stats$compute_effect_size_ci", allow_null = TRUE)
        assert_scalar_bool(cfg$advanced_stats$run_robust_regression,
                           "advanced_stats$run_robust_regression", allow_null = TRUE)
    }

    # 8. QC Settings (optional)
    if (!is.null(cfg$qc)) {
        assert_scalar_bool(cfg$qc$run_umap, "qc$run_umap", allow_null = TRUE)
        assert_scalar_num(cfg$qc$umap_n_neighbors, "qc$umap_n_neighbors",
                          allow_null = TRUE, min_val = 2)
        assert_scalar_num(cfg$qc$outlier_sd_threshold, "qc$outlier_sd_threshold",
                          allow_null = TRUE, min_val = 0)
    }

    invisible(TRUE)
}

validate_clustering_config <- function(clust_cfg) {
    if (is.null(clust_cfg) || isFALSE(clust_cfg$enabled)) {
        return(invisible(TRUE))
    }

    assert_named_list(clust_cfg, "modes$proteomics$clustering")
    assert_scalar_bool(clust_cfg$enabled, "modes$proteomics$clustering$enabled")

    assert_scalar_num(clust_cfg$min_groups, "modes$proteomics$clustering$min_groups",
        allow_null = TRUE, min_val = 2
    )

    if (!is.null(clust_cfg$de_source)) {
        assert_one_of(clust_cfg$de_source, "modes$proteomics$clustering$de_source", c("any_contrast"), allow_null = TRUE)
    }

    if (!is.null(clust_cfg$steps)) {
        assert_named_list(clust_cfg$steps, "modes$proteomics$clustering$steps")

        if (!is.null(clust_cfg$steps$hierarchical)) {
            h <- clust_cfg$steps$hierarchical
            assert_named_list(h, "clustering$steps$hierarchical")
            assert_scalar_bool(h$enabled, "clustering$steps$hierarchical$enabled", allow_null = TRUE)
            assert_scalar_chr(h$distance, "clustering$steps$hierarchical$distance", allow_null = TRUE)
            assert_scalar_chr(h$linkage, "clustering$steps$hierarchical$linkage", allow_null = TRUE)
            assert_scalar_num(h$k, "clustering$steps$hierarchical$k", allow_null = TRUE, min_val = 2)
        }

        if (!is.null(clust_cfg$steps$partition)) {
            p <- clust_cfg$steps$partition
            assert_named_list(p, "clustering$steps$partition")
            assert_scalar_bool(p$enabled, "clustering$steps$partition$enabled", allow_null = TRUE)
            assert_one_of(p$algorithm, "clustering$steps$partition$algorithm", c("pam", "kmeans", "hclust"), allow_null = TRUE)

            if (!is.null(p$k)) assert_scalar_num(p$k, "clustering$steps$partition$k", min_val = 2)
            if (!is.null(p$k_max)) assert_scalar_num(p$k_max, "clustering$steps$partition$k_max", min_val = 2)
            if (is.null(p$k) && is.null(p$k_max)) stop("clustering$steps$partition must define either 'k' or 'k_max'.")

            assert_scalar_chr(p$distance, "clustering$steps$partition$distance", allow_null = TRUE)
            assert_scalar_num(p$nstart, "clustering$steps$partition$nstart", allow_null = TRUE, min_val = 1)
        }

        if (!is.null(clust_cfg$steps$binary_patterns)) {
            b <- clust_cfg$steps$binary_patterns
            assert_named_list(b, "clustering$steps$binary_patterns")
            assert_scalar_bool(b$enabled, "clustering$steps$binary_patterns$enabled", allow_null = TRUE)
            assert_scalar_num(b$corr_cutoff, "clustering$steps$binary_patterns$corr_cutoff", allow_null = TRUE, min_val = -1, max_val = 1)
            assert_scalar_num(b$counts_cutoff, "clustering$steps$binary_patterns$counts_cutoff", allow_null = TRUE)
            assert_scalar_num(b$counts_cutoff_high, "clustering$steps$binary_patterns$counts_cutoff_high", allow_null = TRUE)
            assert_scalar_num(b$counts_cutoff_low, "clustering$steps$binary_patterns$counts_cutoff_low", allow_null = TRUE)
        }

        if (!is.null(clust_cfg$heatmap)) {
            assert_named_list(clust_cfg$heatmap, "clustering$heatmap")
            assert_scalar_num(clust_cfg$heatmap$max_rows, "clustering$heatmap$max_rows", allow_null = TRUE, min_val = 1)
            assert_scalar_bool(clust_cfg$heatmap$cluster_cols, "clustering$heatmap$cluster_cols", allow_null = TRUE)
        }
        return(invisible(TRUE))
    }

    assert_scalar_chr(clust_cfg$method, "modes$proteomics$clustering$method")
    assert_one_of(clust_cfg$method, "modes$proteomics$clustering$method", c("hierarchical", "kmeans", "pam"))
    invisible(TRUE)
}
