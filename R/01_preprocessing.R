#' Normalize RNA-seq counts
#'
#' Options:
#' - "TMMlogCPM" – edgeR TMM normalization + logCPM
#' - "VST"       – DESeq2 variance stabilizing transform
#'
#' @param counts matrix/data.frame of raw counts (genes x samples)
#' @param meta   data.frame of sample metadata (required for VST);
#'               rownames(meta) must match colnames(counts)
#' @param method one of c("TMMlogCPM","VST")
#' @param prior.count numeric added before log in logCPM (default 1)
#' @return numeric matrix; attr(., "method") indicates method used
normalize_counts <- function(counts,
                             meta = NULL,
                             method = c("TMMlogCPM", "VST"),
                             prior.count = 1) {
  stopifnot(is.matrix(counts) || is.data.frame(counts))
  counts <- as.matrix(counts)
  
  # match and normalize legacy alias
  method <- match.arg(method)
  
  if (method == "TMMlogCPM") {
    if (!requireNamespace("edgeR", quietly = TRUE))
      stop("Package 'edgeR' is required. Install via BiocManager::install('edgeR').")
    
    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    mat <- edgeR::cpm(dge, log = TRUE, prior.count = prior.count)
    attr(mat, "method") <- "TMMlogCPM"
    return(mat)
  }
  
  # VST branch
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("Package 'DESeq2' is required. Install via BiocManager::install('DESeq2').")
  
  if (is.null(meta))
    stop("For method = 'VST', 'meta' must be provided.")
  
  if (is.null(rownames(meta)) || !identical(colnames(counts), rownames(meta))) {
    stop("rownames(meta) must exactly match colnames(counts) in order and identity for VST.")
  }
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData   = as.data.frame(meta),
    design    = ~ 1
  )
  keep_nonzero <- rowSums(DESeq2::counts(dds)) > 0
  dds <- dds[keep_nonzero, , drop = FALSE]
  if (nrow(dds) == 0)
    stop("No rows with nonzero counts available for VST.")
  
  mat <- tryCatch({
    # Small matrices: use varianceStabilizingTransformation directly
    if (nrow(dds) < 50) {
      vt <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE, fitType = "mean")
      SummarizedExperiment::assay(vt)
    } else {
      # Larger: use vst with sub-sampling guard
      nsub <- min(1000L, nrow(dds))
      vt <- DESeq2::vst(dds, blind = TRUE, nsub = nsub, fitType = "mean")
      SummarizedExperiment::assay(vt)
    }
  }, error = function(e1) {
    message("[VST] Fallback to TMMlogCPM due to error: ", conditionMessage(e1))
    if (!requireNamespace("edgeR", quietly = TRUE))
      stop("Fallback requires 'edgeR'. Install via BiocManager::install('edgeR').")
    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    fb  <- edgeR::cpm(dge, log = TRUE, prior.count = prior.count)
    attr(fb, "method") <- "TMMlogCPM_fallback_from_VST"
    fb
  })
  
  # If we actually got VST values (not fallback), tag accordingly
  if (is.null(attr(mat, "method"))) attr(mat, "method") <- "VST"
  return(mat)
}

# compute CPM 
compute_cpm <- function(counts) {
  counts <- as.matrix(counts)
  lib <- colSums(counts, na.rm = TRUE)
  if (any(lib == 0 | !is.finite(lib))) stop("Zero/invalid library size for CPM.")
  sweep(counts, 2, lib, "/") * 1e6
}

# compute TPM if gene_length exist 
compute_tpm <- function(counts, gene_lengths_bp) {
  stopifnot(!is.null(rownames(counts)))
  stopifnot(all(c("gene", "length") %in% colnames(gene_lengths_bp)))
  
  common <- intersect(rownames(counts), gene_lengths_bp$gene)
  if (length(common) == 0) stop("No overlapping genes for TPM.")
  
  counts <- as.matrix(counts[common, , drop = FALSE])
  gl_kb <- gene_lengths_bp$length[match(common, gene_lengths_bp$gene)] / 1000
  if (any(!is.finite(gl_kb) | gl_kb <= 0)) stop("Invalid gene lengths.")
  
  rpk <- counts / gl_kb
  sf  <- colSums(rpk, na.rm = TRUE)
  if (any(sf == 0 | !is.finite(sf))) stop("Zero/invalid TPM scaling factor.")
  
  tpm <- sweep(rpk, 2, sf, "/") * 1e6
  rownames(tpm) <- common
  tpm
}

#' Determine which features pass the filter based on min_count per condition
#'
#' Returns a logical vector where TRUE means the feature passed the filter.
#'
#' @param expr_mat Numeric matrix/data.frame (features x samples)
#' @param group Condition/factor for each sample
#' @param min_per_group Named numeric vector of min_count thresholds
#'
#' @return Logical vector (TRUE/FALSE per feature)
pass_filter <- function(expr_mat, group, min_per_group) {
  expr_mat <- as.matrix(expr_mat)
  
  if (length(group) != ncol(expr_mat)) {
    stop("Length of `group` must equal number of columns in `expr_mat`.")
  }
  
  group  <- as.character(group)
  groups <- unique(group)
  
  # Allow single threshold without names
  if (length(min_per_group) == 1 && is.null(names(min_per_group))) {
    min_per_group <- setNames(rep(min_per_group, length(groups)), groups)
  }
  
  if (!all(groups %in% names(min_per_group))) {
    stop(
      "min_per_group must have names for all groups. Missing: ",
      paste(setdiff(groups, names(min_per_group)), collapse = ", ")
    )
  }
  
  passes_per_group <- sapply(groups, function(g) {
    cols <- which(group == g)
    if (length(cols) == 0) return(rep(FALSE, nrow(expr_mat)))
    
    sums <- rowSums(!is.na(expr_mat[, cols, drop = FALSE]))
    sums >= min_per_group[[g]]
  })
  
  apply(passes_per_group, 1, any)
}

#' Extract min_count per condition from config
#'
#' Supports all of these config shapes:
#' - min_count: 3
#' - min_count: { default: 3 }
#' - min_count: { default: 3, C: 2 }
#' - min_count: { C: 3, S: 3, SH: 2 }
#'
#' @param min_cfg config$modes$proteomics$filtering$min_count
#' @param groups  character vector of condition values (one per sample)
#'
#' @return Named numeric vector of length unique(groups)
extract_min_count <- function(min_cfg, groups) {
  groups <- unique(as.character(groups))
  
  # Case 0: numeric scalar (e.g. min_count: 3)
  if (is.numeric(min_cfg) && length(min_cfg) == 1 && is.null(names(min_cfg))) {
    return(setNames(rep(min_cfg, length(groups)), groups))
  }
  
  # Case 1: list-like with $default and optional overrides
  if (!is.null(min_cfg$default)) {
    out <- setNames(rep(min_cfg$default, length(groups)), groups)
    
    overrides <- min_cfg[setdiff(names(min_cfg), "default")]
    if (length(overrides) > 0) {
      for (nm in names(overrides)) {
        if (nm %in% names(out)) {
          out[nm] <- overrides[[nm]]
        }
      }
    }
    return(out)
  }
  
  # Case 2: fully named thresholds (no default)
  out <- unlist(min_cfg)
  names(out) <- names(min_cfg)
  out
}


filter_features_dynamic <- function(norm_mat, meta, sample_col, group_col, threshold = 3) {
  stopifnot(sample_col %in% colnames(meta))
  stopifnot(group_col %in% colnames(meta))
  
  # strict alignment (no silent dropping)
  if (!all(colnames(norm_mat) %in% meta[[sample_col]])) {
    miss <- setdiff(colnames(norm_mat), meta[[sample_col]])
    stop("Samples missing in metadata: ", paste(miss, collapse = ", "))
  }
  meta2 <- meta[match(colnames(norm_mat), meta[[sample_col]]), , drop = FALSE]
  
  group_list <- split(meta2[[sample_col]], meta2[[group_col]])
  
  keep_feature <- function(expr_values) {
    for (grp in names(group_list)) {
      samp <- group_list[[grp]]
      sub_vals <- expr_values[samp]
      n_reps <- sum(!is.na(sub_vals))
      min_pass <- ceiling(n_reps / 2)
      if (sum(sub_vals >= threshold, na.rm = TRUE) >= min_pass) return(TRUE)
    }
    FALSE
  }
  
  keep_vec <- apply(norm_mat, 1, keep_feature)
  
  list(
    filtered = norm_mat[keep_vec, , drop = FALSE],
    keep_vec = keep_vec
  )
}
#' Perseus-style imputation (downshifted & narrowed normal distribution)
#'
#' For each sample (column), missing values are imputed from a normal
#' distribution with:
#'   mean = mean(sample, na.rm = TRUE) - downshift * sd(sample, na.rm = TRUE)
#'   sd   = width * sd(sample, na.rm = TRUE)
#'
#' @param expr_mat  Numeric matrix/data.frame (features x samples), usually log2.
#' @param width  Fraction of the sample sd used as sd of the imputation dist.
#' @param downshift How many sd to subtract from the mean for the imputation mean.
#'
#' @return A list with:
#'   - imputed: numeric matrix with imputed values
#'   - imputed_flag: logical matrix, TRUE where a value was imputed
perseus_impute_with_flags <- function(expr_mat, width = 0.3, downshift = 1.8) {
  # Coerce to matrix for safety
  expr_mat <- as.matrix(expr_mat)
  sample_cols <- colnames(expr_mat)
  
  # Logical matrix of missing values (these will be imputed)
  imputed_flag <- is.na(expr_mat)
  
  # Copy for imputation
  imputed <- expr_mat
  
  for (j in seq_len(ncol(imputed))) {
    x <- imputed[, j]
    
    # stop if all values are NA
    if (all(is.na(x))) {
      stop("Imputation failed: sample '", colnames(imputed)[j], "' is all-NA after filtering.")
    }
    
    
    # Compute sd and mean of observed values
    obs <- x[!is.na(x)]
    s <- stats::sd(obs)
    if (!is.finite(s) || s == 0) s <- 1e-8
    m <- mean(obs)
    
    imp_sd   <- width * s
    imp_mean <- m - downshift * s
    
    n_missing <- sum(is.na(x))
    if (n_missing > 0) {
      x[is.na(x)] <- stats::rnorm(n_missing, mean = imp_mean, sd = imp_sd)
    }
    
    imputed[, j] <- x
  }
  
  colnames(imputed)      <- sample_cols
  rownames(imputed_flag) <- rownames(imputed)
  colnames(imputed_flag) <- sample_cols
  
  list(
    imputed      = imputed,
    imputed_flag = imputed_flag
  )
}

make_imputations_proteomics <- function(expr_mat, cfg, verbose = FALSE) {
  
  imp_cfg <- cfg$modes$proteomics$imputation
  n_imputations <- as.integer(imp_cfg$no_repetitions)
  
  seed_base <- cfg$params$seed
  if (is.null(seed_base) || is.na(seed_base)) {
    stop("Missing cfg$params$seed. Set it in the YAML to ensure reproducibility.")
  }
  
  stopifnot(is.matrix(expr_mat))
  stopifnot(!is.na(n_imputations), n_imputations >= 1)
  
  imps <- vector("list", n_imputations)
  
  for (i in seq_len(n_imputations)) {
    if (isTRUE(verbose)) {
      message(sprintf("Imputation: %d / %d", i, n_imputations))
    }
    
    set.seed(as.integer(seed_base) + i)
    
    expr_imp_i <- impute_proteomics_perseus(
      expr_mat,
      cfg = cfg$modes$proteomics,
      return_flags = FALSE
    )
    
    stopifnot(is.matrix(expr_imp_i))
    stopifnot(identical(dim(expr_imp_i), dim(expr_mat)))
    stopifnot(identical(colnames(expr_imp_i), colnames(expr_mat)))
    stopifnot(identical(rownames(expr_imp_i), rownames(expr_mat)))
    
    imps[[i]] <- expr_imp_i
  }
  
  imps
}


as_legacy_pre <- function(pre) {
  pre$expr_filt <- pre$expr_filt_mat %||% pre$expr_filt
  pre$expr_imp  <- pre$expr_imp_single %||% pre$expr_imp
  
  if (is.null(pre$imputation)) pre$imputation <- list()
  if (is.null(pre$imputation$imputed_flag) && !is.null(pre$imputation_qc$imputed_flag)) {
    pre$imputation$imputed_flag <- pre$imputation_qc$imputed_flag
  }
  
  pre
}

#' Build a standardized proteomics expression object (engine-aware)
#'
#' For now this supports DIANN only, but keeps a stable contract so that
#' adding MaxQuant later does not require changing downstream code.
#'
#' Contract:
#'   - assay_log2   : matrix (features × samples), log2 scale (downstream contract)
#'   - assay_linear : optional linear-scale matrix (NULL unless available/needed)
#'   - row_data     : feature annotations aligned to assay rows
#'   - col_data     : sample metadata aligned to assay columns
#'   - info         : engine/scale metadata
#'
#' @param inputs List from load_proteomics_inputs().
#' @param config Full config list.
#'
#' @return List with fields: assay_log2, assay_linear, row_data, col_data, info.
get_proteomics_expression_matrix <- function(inputs, config) {
  
  cfg <- config$modes$proteomics
  
  # Build matrix depending on engine
  if (cfg$engine == "DIANN") {
    
    # DIA-NN → expression matrix (your current builder already returns log2)
    assay_log2 <- get_measurements_per_sample_diann(
      protein    = inputs$protein,
      sample_map = inputs$sample_map,
      meta       = inputs$metadata,
      cfg        = cfg
    )
    
    assay_linear <- NULL
    
    # Feature annotations (row_data)
    row_data <- inputs$protein[ , c(cfg$id_columns$protein_id, unlist(cfg$id_columns$protein_annot)), drop = FALSE]
    
  } else {
    stop(sprintf("Unsupported proteomics engine: %s", cfg$engine))
  }
  
  # Align col_data to assay columns 
  col_data <- inputs$metadata
  # sample_col <- cfg$effects$samples %||% cfg$id_columns$sample_col
  col_data <- align_meta_to_expr(assay_log2, inputs$metadata, cfg)
  
  list(
    assay_log2   = assay_log2,
    assay_linear = assay_linear,
    row_data     = row_data,
    col_data     = col_data,
    info = list(
      mode         = "proteomics",
      engine       = cfg$engine,
      scale_in     = cfg$scale_in %||% "log2",
      target_scale = cfg$transform$target_scale %||% "log2"
    )
  )
}
#' Load omics inputs (generic)
#'
#' Generic loader for any omics mode (proteomics / rna / etc.).
#' Expects:
#'   config$modes[[mode]]$files – a named list of paths
#'   config$modes[[mode]]$engine (optional)
#'
#' @param config list as returned by load_config()
#' @param mode character, e.g. "proteomics", "rna"
#' @return named list of tibbles + optional engine element
load_omics_inputs <- function(config, mode = c("proteomics", "rna")) {
  mode <- match.arg(mode)
  
  # ---- get mode config ----
  cfg <- config$modes[[mode]]
  if (is.null(cfg)) {
    stop("No configuration found for mode '", mode, "' in config$modes$")
  }
  if (is.null(cfg$files) || !length(cfg$files)) {
    stop("config$modes$", mode, "$files is missing or empty.")
  }
  
  files <- cfg$files
  inputs <- list()
  
  # resolve relative -> absolute under project.dir + paths.raw
  for (nm in names(files)) {
    rel <- files[[nm]]
    if (is.null(rel) || !nzchar(rel)) {
      stop("In mode '", mode, "': path for '", nm, "' is missing/empty in cfg$files.")
    }
    
    abs <- resolve_raw_path(config, rel)
    files[[nm]] <- abs
    
    if (!file.exists(abs)) {
      stop("In mode '", mode, "': file not found for '", nm, "': ", abs)
    }
  }
  
  inputs <- lapply(files, function(path) readr::read_csv(path, show_col_types = FALSE))
  names(inputs) <- names(files)
  
  # ---- attach engine if present ----
  if (!is.null(cfg$engine)) {
    inputs$engine <- cfg$engine
  }
  
  inputs
}


read_table_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else {
    readr::read_csv(path, show_col_types = FALSE)
  }
}


#' Load proteomics inputs 
#'
#' @param config list as returned by load_config()
#' @return list (protein, sample_map, meta, contrasts, engine, ...)
load_proteomics_inputs <- function(config) {
  load_omics_inputs(config, mode = "proteomics")
}

#' Load rna inputs
#'
#' @param config list as returned by load_config()
#' @return list (genes, sample_map, meta, contrasts, engine, ...)
load_rna_inputs <- function(config) {
  load_omics_inputs(config, mode = "rna")
}


#' Load pipeline configuration from YAML.
#'
#' @param path Path to a YAML config file (default: "config/config.yml").
#' @return A named list with the parsed configuration.
load_config <- function(path = "config/config.yml") {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Install it with install.packages('yaml').")
  }
  
  cfg <- yaml::read_yaml(path)
  
  cfg$.config_path  <- normalizePath(path, mustWork = FALSE)
  cfg$.config_mtime <- file.info(path)$mtime
  
  return(cfg)
}

# Get effects (color / shape / samples) for a given modality (e.g. "rna", "proteomics")
get_effects <- function(cfg, mode) {
  mode_cfg <- cfg$modes[[mode]]
  if (is.null(mode_cfg)) {
    stop(sprintf("Mode '%s' not found under config$modes.", mode))
  }
  
  eff <- mode_cfg$effects
  if (is.null(eff)) eff <- list()
  
  list(
    color   = eff$color   %||% NULL,
    shape   = eff$shape   %||% NULL,
    samples = eff$samples %||% "SampleID"
  )
}

#' Validate a list of imputed matrices + metadata alignment
#'
#' @param imputations list of numeric matrices (features x samples)
#' @param meta data.frame with one row per sample
#' @param sample_col column name in meta holding sample IDs
#' @param warn_extra_meta warn if meta contains samples not in matrix
#' @param allow_na allow NA values in matrices
validate_imputations <- function(imputations,
                                 meta,
                                 sample_col,
                                 warn_extra_meta = TRUE,
                                 allow_na = FALSE) {
  if (!is.list(imputations) || length(imputations) < 1) {
    stop("imputations must be a non-empty list of matrices.")
  }
  if (!is.data.frame(meta)) stop("meta must be a data.frame.")
  if (is.null(sample_col) || !nzchar(sample_col)) stop("sample_col must be provided.")
  if (!(sample_col %in% colnames(meta))) stop("meta is missing sample_col: ", sample_col)
  
  ref <- imputations[[1]]
  assert_numeric_matrix(ref, "imputations[[1]]")  # rownames+colnames+unique
  
  ref_dim <- dim(ref)
  ref_rn  <- rownames(ref)
  ref_cn  <- colnames(ref)
  
  # ---- validate each imputation matrix ----
  for (i in seq_along(imputations)) {
    x <- imputations[[i]]
    if (!is.matrix(x)) stop(sprintf("Imputation %d is not a matrix.", i))
    if (!is.numeric(x)) stop(sprintf("Imputation %d is not numeric (storage: %s).", i, typeof(x)))
    
    if (!identical(dim(x), ref_dim)) {
      stop(sprintf(
        "Imputation %d has different dimensions. Expected %s, got %s.",
        i, paste(ref_dim, collapse = "x"), paste(dim(x), collapse = "x")
      ))
    }
    if (!identical(rownames(x), ref_rn)) {
      j <- which(rownames(x) != ref_rn)[1]
      stop(sprintf(
        "Imputation %d has different rownames/order at position %d: expected '%s' got '%s'.",
        i, j, ref_rn[j], rownames(x)[j]
      ))
    }
    if (!identical(colnames(x), ref_cn)) {
      j <- which(colnames(x) != ref_cn)[1]
      stop(sprintf(
        "Imputation %d has different colnames/order at position %d: expected '%s' got '%s'.",
        i, j, ref_cn[j], colnames(x)[j]
      ))
    }
    
    if (!isTRUE(allow_na) && anyNA(x)) {
      stop(sprintf("Imputation %d still contains NA values.", i))
    }
    if (!all(is.finite(x[is.na(x) == FALSE]))) {
      stop(sprintf("Imputation %d contains non-finite values (Inf/-Inf/NaN).", i))
    }
  }
  
  # ---- meta checks vs matrix ----
  meta_ids <- as.character(meta[[sample_col]])
  if (anyDuplicated(meta_ids) > 0) {
    dups <- unique(meta_ids[duplicated(meta_ids)])
    stop("meta[[", sample_col, "]] has duplicated sample IDs (e.g. ", paste(head(dups, 3), collapse = ", "), ").")
  }
  
  missing_in_meta <- setdiff(ref_cn, meta_ids)
  if (length(missing_in_meta) > 0) {
    stop(sprintf(
      "meta is missing %d samples present in expression matrix (sample_col='%s'), e.g. %s",
      length(missing_in_meta), sample_col, paste(head(missing_in_meta, 3), collapse = ", ")
    ))
  }
  
  extra_in_meta <- setdiff(meta_ids, ref_cn)
  if (length(extra_in_meta) > 0 && isTRUE(warn_extra_meta)) {
    warning(sprintf(
      "meta contains %d samples not present in expression matrix (sample_col='%s'), e.g. %s",
      length(extra_in_meta), sample_col, paste(head(extra_in_meta, 3), collapse = ", ")
    ))
  }
  
  invisible(TRUE)
}

#' Validate consistency of proteomics imputations and metadata sample IDs
#'
#' Fail-fast validation for runtime objects (not file I/O, not config schema).
#' This is designed to catch silent misalignment bugs before DE/QC steps.
#'
#' @param imputations list of imputed expression matrices (features x samples)
#' @param meta data.frame with one row per sample
#' @param cfg full config list (expects modes$proteomics$id_columns$sample_col)
#' @param warn_extra_meta logical; warn if meta contains samples not present in matrix
#' @param allow_na logical; if TRUE, allow NA values in matrices (useful pre-imputation / NA-QC)
#'
#' @return invisibly TRUE; otherwise stops with informative error

validate_proteomics_imputations <- function(imputations,
                                            meta,
                                            cfg,
                                            warn_extra_meta = TRUE,
                                            allow_na = FALSE,
                                            sample_col = NULL) {
  if (is.null(sample_col) || !nzchar(sample_col)) {
    p <- cfg$modes$proteomics
    sample_col <- p$effects$samples %||% p$id_columns$sample_col %||% "SampleID"
  }
  
  validate_imputations(
    imputations     = imputations,
    meta            = meta,
    sample_col      = sample_col,
    warn_extra_meta = warn_extra_meta,
    allow_na        = allow_na
  )
}

#' Check that required columns exist in a data frame
#'
#' @param df A data.frame / tibble
#' @param required Character vector of required column names
#' @param df_name Optional name of the data frame (for error messages)
check_has_cols <- function(df, required, df_name = deparse(substitute(df))) {
  required <- unlist(required)
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(
      "In table '", df_name, "': missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

#' Check that all values in x are present in y
#'
#' @param x Vector of values that must be contained in y
#' @param y Reference vector
#' @param label_x Name/label for x (used in error message)
#' @param label_y Name/label for y (used in error message)
check_all_in <- function(x, y, label_x = "x", label_y = "y") {
  missing <- setdiff(x, y)
  if (length(missing) > 0) {
    stop(
      "Values in ", label_x, " not found in ", label_y, ": ",
      paste(head(missing, 10), collapse = ", "),
      if (length(missing) > 10) " ... (truncated)" else "",
      call. = FALSE
    )
  }
}

#' Generic omics input validation
#'
#' This function loads inputs for a given omics mode (e.g. proteomics, rna)
#' using \code{load_omics_inputs()}, and then dispatches to a mode-specific
#' validation function.
#'
#' @param config List as returned by \code{load_config()}
#' @param mode   Character, e.g. "proteomics", "rna"
#'
#' @return (invisibly) the loaded inputs (a named list)
validate_omics_inputs <- function(config,
                                  mode = c("proteomics", "rna",
                                           "metabolomics", "lipidomics")) {
  mode <- match.arg(mode)
  
  # Load inputs using the generic loader
  inputs <- load_omics_inputs(config, mode = mode)
  cfg    <- config$modes[[mode]]
  
  # Dispatch to mode-specific validator
  switch(
    mode,
    proteomics  = validate_proteomics_inputs(inputs,  cfg),
    rna         = validate_rna_inputs(inputs,         cfg),
    metabolomics = validate_metabolomics_inputs(inputs, cfg),
    lipidomics   = validate_lipidomics_inputs(inputs,   cfg)
  )
  
  invisible(inputs)
}

#' Validate proteomics inputs (wide matrix: samples in columns)
#'
#' Assumes:
#' - protein: one row per protein, one column per sample (wide format)
#'            with non-sample columns defined in cfg$id_columns
#' - sample_map: maps from raw sample names (map_from) to unified sample IDs (map_to)
#' - metadata: one row per sample, with sample ID = effects$samples
#'
#' @param inputs List returned by load_proteomics_inputs()
#' @param cfg    config$modes$proteomics
#'
#' @return TRUE (invisibly) if validation passes; otherwise stops with an error
validate_proteomics_inputs <- function(inputs, cfg) {
  protein    <- inputs$protein
  sample_map <- inputs$sample_map
  meta       <- inputs$metadata
  contrasts  <- inputs$contrasts
  
  id_cols  <- cfg$id_columns
  eff_cols <- cfg$effects
  
  # 1) Column existence checks ----------------------------------
  
  # protein: must have the protein ID column
  check_has_cols(
    protein,
    id_cols$protein_id,
    df_name = "protein"
  )
  
  # sample_map: must have mapping columns
  check_has_cols(
    sample_map,
    c(id_cols$map_from, id_cols$map_to),
    df_name = "sample_map"
  )
  
  # metadata: must have sample label column + effect columns (color/shape/label)
  meta_required <- unique(c(
    eff_cols$samples,
    eff_cols$color,
    eff_cols$shape
  ))
  meta_required <- meta_required[!is.na(meta_required)]
  
  check_has_cols(
    meta,
    meta_required,
    df_name = "metadata"
  )
  
  # 2) Determine sample IDs from protein matrix -----------------
  
  # Non-sample columns: protein_id + optional annotation columns
  annot_cols <- id_cols$protein_annot
  if (is.null(annot_cols)) {
    annot_cols <- character(0)
  } else {
    annot_cols <- unlist(annot_cols)
  }
  
  non_sample_cols <- unique(c(id_cols$protein_id, annot_cols))
  
  # Wide format: all remaining columns are assumed to be samples
  protein_sample_cols <- setdiff(colnames(protein), non_sample_cols)
  
  # 3) Cross-table consistency checks ---------------------------
  
  # Here we assume:
  # - protein column names correspond to the *raw* sample names
  # - those raw names are stored in sample_map[[map_from]]
  #
  # So: all protein sample columns must exist in sample_map map_from
  check_all_in(
    x = protein_sample_cols,
    y = sample_map[[id_cols$map_from]],
    label_x = "protein sample columns",
    label_y = "sample_map map_from"
  )
  
  # And all unified sample IDs (map_to) must exist in metadata samples column
  meta_sample_ids <- meta[[eff_cols$samples]]
  
  check_all_in(
    x = sample_map[[id_cols$map_to]],
    y = meta_sample_ids,
    label_x = "sample_map map_to",
    label_y = "metadata label"
  )
  
  invisible(TRUE)
}


#' Validate RNA inputs
#'
#' Placeholder / skeleton. To be filled once RNA files and id_columns
#' are fully defined in the config.
#'
#' @param inputs List returned by \code{load_rna_inputs()}
#' @param cfg    The \code{config$modes$rna} sub-list
validate_rna_inputs <- function(inputs, cfg) {
  # Example structure; adjust once your RNA inputs format is finalized
  # counts     <- inputs$counts
  # sample_map <- inputs$sample_map
  # meta       <- inputs$metadata
  # contrasts  <- inputs$contrasts
  #
  # id_cols  <- cfg$id_columns
  # eff_cols <- cfg$effects
  #
  # check_has_cols(counts, id_cols$gene_id, df_name = "counts")
  #
  # (Add consistency checks similar to proteomics)
  
  invisible(TRUE)
}

#' Validate metabolomics inputs
#'
#' Placeholder function. Implement when metabolomics file structure is defined.
validate_metabolomics_inputs <- function(inputs, cfg) {
  # Implement when metabolomics inputs are added to the config
  invisible(TRUE)
}

#' Validate lipidomics inputs
#'
#' Placeholder function. Implement when lipidomics file structure is defined.
validate_lipidomics_inputs <- function(inputs, cfg) {
  # Implement when lipidomics inputs are added to the config
  invisible(TRUE)
}

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
