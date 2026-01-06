#' Validate global pipeline configuration
#'
#' Runs fail-fast checks on the YAML-derived config object.
#' Ensures critical parameters exist and have valid types/values before execution.
#'
#' @param config Full config list (from config.yaml).
#' @return Invisibly TRUE if all checks pass; otherwise stops with an error.
validate_config <- function(config) {
  assert_named_list(config, "config")
  
  # Validate basic required fields structure
  assert_named_list(config$paths,  "config$paths")
  assert_named_list(config$modes,  "config$modes")
  assert_named_list(config$params, "config$params")
  
  # Validate specific scalar values
  assert_scalar_chr(config$paths$out, "config$paths$out")
  assert_scalar_num(config$params$seed, "config$params$seed")
  
  # Mode-specific validations (only if the mode is enabled/present)
  if (!is.null(config$modes$proteomics))   validate_proteomics_config(config$modes$proteomics)
  if (!is.null(config$modes$rna))          validate_rna_config(config$modes$rna)
  if (!is.null(config$modes$metabolomics)) validate_metabolomics_config(config$modes$metabolomics)
  if (!is.null(config$modes$lipidomics))   validate_lipidomics_config(config$modes$lipidomics)
  
  invisible(TRUE)
}

# ---- Validation Helpers ----

assert_named_list <- function(x, name) {
  if (is.null(x) || !is.list(x)) stop(sprintf("'%s' must be a list", name))
  invisible(TRUE)
}

assert_scalar_bool <- function(x, name, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(invisible(TRUE))
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop(sprintf("'%s' must be TRUE/FALSE", name))
  }
  invisible(TRUE)
}

assert_scalar_chr <- function(x, name, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(invisible(TRUE))
  if (!is.character(x) || length(x) != 1 || is.na(x) || x == "") {
    stop(sprintf("'%s' must be a non-empty string", name))
  }
  invisible(TRUE)
}

assert_scalar_num <- function(x, name, allow_null = FALSE, min_val = -Inf, max_val = Inf) {
  if (allow_null && is.null(x)) return(invisible(TRUE))
  if (!is.numeric(x) || length(x) != 1 || is.na(x)) {
    stop(sprintf("'%s' must be a number", name))
  }
  if (x < min_val || x > max_val) {
    stop(sprintf("'%s' must be between %s and %s", name, min_val, max_val))
  }
  invisible(TRUE)
}

assert_one_of <- function(x, name, choices, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(invisible(TRUE))
  assert_scalar_chr(x, name)
  if (!(tolower(x) %in% tolower(choices))) {
    stop(sprintf("'%s' must be one of: %s", name, paste(choices, collapse = ", ")))
  }
  invisible(TRUE)
}

# ---- Clustering Validation ----

validate_clustering_config <- function(clust_cfg) {
  # Missing or disabled -> OK
  if (is.null(clust_cfg) || isFALSE(clust_cfg$enabled)) return(invisible(TRUE))
  
  assert_named_list(clust_cfg, "modes$proteomics$clustering")
  assert_scalar_bool(clust_cfg$enabled, "modes$proteomics$clustering$enabled")
  
  # Optional
  assert_scalar_num(clust_cfg$min_groups, "modes$proteomics$clustering$min_groups",
                    allow_null = TRUE, min_val = 2)
  
  # Optional DE source (keep your current logic)
  if (!is.null(clust_cfg$de_source)) {
    assert_one_of(clust_cfg$de_source,
                  "modes$proteomics$clustering$de_source",
                  c("any_contrast"),
                  allow_null = TRUE)
  }
  
  # ---- NEW STYLE: steps ----
  if (!is.null(clust_cfg$steps)) {
    assert_named_list(clust_cfg$steps, "modes$proteomics$clustering$steps")
    
    # hierarchical
    if (!is.null(clust_cfg$steps$hierarchical)) {
      h <- clust_cfg$steps$hierarchical
      assert_named_list(h, "clustering$steps$hierarchical")
      assert_scalar_bool(h$enabled, "clustering$steps$hierarchical$enabled", allow_null = TRUE)
      assert_scalar_chr(h$distance, "clustering$steps$hierarchical$distance", allow_null = TRUE)
      assert_scalar_chr(h$linkage,  "clustering$steps$hierarchical$linkage",  allow_null = TRUE)
      assert_scalar_num(h$k, "clustering$steps$hierarchical$k", allow_null = TRUE, min_val = 2)
    }
    
    # partition
    if (!is.null(clust_cfg$steps$partition)) {
      p <- clust_cfg$steps$partition
      assert_named_list(p, "clustering$steps$partition")
      assert_scalar_bool(p$enabled, "clustering$steps$partition$enabled", allow_null = TRUE)
      
      # algorithm pam/kmeans
      assert_one_of(p$algorithm, "clustering$steps$partition$algorithm",
                    c("pam", "kmeans"), allow_null = TRUE)
      
      # either k fixed OR k_max (for auto selection)
      if (!is.null(p$k)) {
        assert_scalar_num(p$k, "clustering$steps$partition$k", min_val = 2)
      }
      if (!is.null(p$k_max)) {
        assert_scalar_num(p$k_max, "clustering$steps$partition$k_max", min_val = 2)
      }
      if (is.null(p$k) && is.null(p$k_max)) {
        stop("clustering$steps$partition must define either 'k' (fixed) or 'k_max' (auto).")
      }
      
      assert_scalar_chr(p$distance, "clustering$steps$partition$distance", allow_null = TRUE)
      assert_scalar_num(p$nstart, "clustering$steps$partition$nstart", allow_null = TRUE, min_val = 1)
    }
    
    # binary patterns
    if (!is.null(clust_cfg$steps$binary_patterns)) {
      b <- clust_cfg$steps$binary_patterns
      assert_named_list(b, "clustering$steps$binary_patterns")
      assert_scalar_bool(b$enabled, "clustering$steps$binary_patterns$enabled", allow_null = TRUE)
      assert_scalar_num(b$corr_cutoff, "clustering$steps$binary_patterns$corr_cutoff",
                        allow_null = TRUE, min_val = -1, max_val = 1)
      assert_scalar_num(b$counts_cutoff, "clustering$steps$binary_patterns$counts_cutoff",
                        allow_null = TRUE)
    }
    
    # heatmap common settings (optional)
    if (!is.null(clust_cfg$heatmap)) {
      assert_named_list(clust_cfg$heatmap, "clustering$heatmap")
      assert_scalar_num(clust_cfg$heatmap$max_rows, "clustering$heatmap$max_rows",
                        allow_null = TRUE, min_val = 1)
      assert_scalar_bool(clust_cfg$heatmap$cluster_cols, "clustering$heatmap$cluster_cols",
                         allow_null = TRUE)
    }
    
    return(invisible(TRUE))
  }
  
  # ---- OLD STYLE (backwards compatible): method ----
  # If no steps block, fallback to old 'method' validation
  assert_scalar_chr(clust_cfg$method, "modes$proteomics$clustering$method")
  assert_one_of(clust_cfg$method, "modes$proteomics$clustering$method",
                c("hierarchical", "kmeans", "pam"))
  
  invisible(TRUE)
}


# ---- Proteomics Validation ----

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
    assert_one_of(cfg$imputation$method, "imputation$method", c("none", "perseus", "dep2"))
    
    assert_scalar_num(cfg$imputation$no_repetitions, "imputation$no_repetitions", min_val = 1)
    assert_scalar_num(cfg$imputation$min_no_passed, "imputation$min_no_passed", min_val = 1)
    
    if (cfg$imputation$min_no_passed > cfg$imputation$no_repetitions) {
      stop("imputation$min_no_passed cannot be greater than imputation$no_repetitions")
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
  }
  
  # 4. de Settings
  if (!is.null(cfg$de)) {
    assert_scalar_num(cfg$de$p_cutoff, "de$p_cutoff",
                      allow_null = TRUE, min_val = .Machine$double.eps, max_val = 1)
    assert_scalar_num(cfg$de$linear_fc_cutoff, "de$linear_fc_cutoff", allow_null = TRUE, min_val = 1)
  }
  
  invisible(TRUE)
  
  # 5. Clustering Settings (optional)
  if (!is.null(cfg$clustering)) {
    validate_clustering_config(cfg$clustering)
  }
  
}

# --- Future Mode Stubs ---
validate_rna_config <- function(cfg) invisible(TRUE)
validate_metabolomics_config <- function(cfg) invisible(TRUE)
validate_lipidomics_config <- function(cfg) invisible(TRUE)