# R/domain/lipidomics/00_inputs.R
#
# Loading and parsing lipidomics inputs.
# Supports:
#   - "lipidsearch"     : LipidSearch export (Excel or CSV with annotation columns)
#   - "translated_csv"  : Translated/clean CSV (Metabolite col + sample cols, optional Label row)
#   - "processed_wide"  : Already-processed wide table (clean sample columns)
#
# Reuses: resolve_raw_path, read_table_auto, check_has_cols,
#         coerce_df_to_numeric_matrix, assert_numeric_matrix, assert_meta_contract


# ---- public entry point -----------------------------------------------------

#' Load lipidomics inputs from config
#'
#' @param config Full pipeline config (list).
#' @return list with: data (raw data.frame), metadata (data.frame or NULL),
#'         format
load_lipidomics_inputs <- function(config) {
    cfg <- config$modes$lipidomics
    if (is.null(cfg)) stop("No config for mode lipidomics")

    files <- cfg$files
    if (is.null(files$data)) stop("lipidomics config$files$data is required")

    abs_data <- resolve_raw_path(config, files$data)
    if (!file.exists(abs_data)) stop("Lipidomics data file not found: ", abs_data)

    data_df <- read_lipid_file(abs_data, sheet = cfg$input$sheet)

    # Optional metadata
    meta <- NULL
    if (!is.null(files$metadata) && nzchar(files$metadata)) {
        abs_meta <- resolve_raw_path(config, files$metadata)
        if (!file.exists(abs_meta)) stop("Metadata file not found: ", abs_meta)
        meta <- read_table_auto(abs_meta)
    }

    # Handle group-row format: row 1 of data contains condition assignments
    if (isTRUE(cfg$input$has_group_row) && nrow(data_df) > 0) {
        group_row <- data_df[1, , drop = FALSE]
        data_df   <- data_df[-1, , drop = FALSE]
        rownames(data_df) <- NULL

        if (is.null(meta)) {
            id_cfg     <- cfg$id_columns %||% list()
            annot_cols <- id_cfg$annotation_cols %||% character(0)
            # Also exclude the feature name column from sample columns
            name_col <- id_cfg$feature_id_col %||% id_cfg$name_col %||% "Metabolite"
            annot_cols <- unique(c(annot_cols, name_col))
            sample_cols <- setdiff(colnames(data_df), annot_cols)

            cond_col <- cfg$de$condition_column %||%
                        cfg$effects$color %||% "condition"

            meta <- data.frame(
                sample_id = sample_cols,
                stringsAsFactors = FALSE
            )
            meta[[cond_col]] <- as.character(unlist(group_row[1, sample_cols]))
            message(sprintf(
                "lipidomics: built metadata from group row (%d samples, condition='%s')",
                length(sample_cols), cond_col
            ))
        }
    }

    list(
        data     = data_df,
        metadata = meta,
        format   = cfg$input$format %||% "lipidsearch"
    )
}


#' Validate lipidomics config (called by validate_config dispatch)
validate_lipidomics_config <- function(cfg) {
    assert_one_of(cfg$input$format, "input$format",
                  c("lipidsearch", "translated_csv", "processed_wide"))

    if (is.null(cfg$files$data) || !nzchar(cfg$files$data)) {
        stop("lipidomics: files$data is required")
    }

    norm <- cfg$normalization
    if (!is.null(norm)) {
        assert_one_of(norm$sample_norm, "normalization$sample_norm",
                      c("none", "sum", "median", "pqn", "is"),
                      allow_null = TRUE)
        assert_one_of(norm$transform, "normalization$transform",
                      c("none", "log2", "log10", "glog10"),
                      allow_null = TRUE)
        assert_one_of(norm$scaling, "normalization$scaling",
                      c("none", "center", "auto", "pareto", "range"),
                      allow_null = TRUE)
        assert_one_of(norm$na_policy, "normalization$na_policy",
                      c("keep", "zero", "min_half", "lod"),
                      allow_null = TRUE)
    }

    ff <- cfg$feature_filter
    if (!is.null(ff) && !is.null(ff$max_group_missing_frac)) {
        v <- ff$max_group_missing_frac
        if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0 || v > 1) {
            stop("lipidomics: feature_filter$max_group_missing_frac must be a single numeric in [0, 1]")
        }
    }

    invisible(TRUE)
}


# ---- file reader ------------------------------------------------------------

#' Read lipidomics data file (xlsx or csv/tsv)
#' @param path Absolute file path.
#' @param sheet Sheet name or index for xlsx (optional).
#' @return data.frame
read_lipid_file <- function(path, sheet = NULL) {
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("xlsx", "xls")) {
        if (!requireNamespace("readxl", quietly = TRUE)) {
            stop("Package 'readxl' is required to read Excel files.")
        }
        df <- if (!is.null(sheet)) {
            readxl::read_excel(path, sheet = sheet)
        } else {
            readxl::read_excel(path)
        }
        as.data.frame(df)
    } else {
        read_table_auto(path)
    }
}


# ---- format-specific parsers ------------------------------------------------

#' Parse LipidSearch export -> expr_raw + row_data
#'
#' LipidSearch exports have annotation columns (Molecule List Name, Molecule,
#' RT, ClassKey, SubClassKey, FAKey, FAGroupKey, TotalGrade, Precursor Mz)
#' followed by sample intensity columns. Pool/blank columns are included.
#'
#' @param data_df  Raw LipidSearch data.frame.
#' @param cfg      lipidomics mode config.
#' @param meta     Metadata data.frame (optional).
#' @return list(expr_raw, row_data, sample_ids)
parse_lipidsearch <- function(data_df, cfg, meta = NULL) {
    id_cfg <- cfg$id_columns %||% list()

    # LipidSearch annotation columns (configurable)
    default_annot <- c("Molecule List Name", "Molecule",
                       "Average Measured Retention Time",
                       "LipidMolecGroup", "ClassKey", "SubClassKey",
                       "FAKey", "FAGroupKey", "TotalGrade", "Precursor Mz")
    annot_cols <- id_cfg$annotation_cols %||% default_annot

    all_cols <- colnames(data_df)
    # Trim whitespace from column names
    colnames(data_df) <- trimws(colnames(data_df))
    all_cols <- colnames(data_df)

    # Identify annotation vs sample columns
    annot_present <- intersect(annot_cols, all_cols)
    sample_cols <- setdiff(all_cols, annot_present)

    if (length(sample_cols) == 0) {
        stop("No sample columns found after removing annotation columns.")
    }

    # Snapshot pool columns BEFORE exclude_patterns strips them. Captured for
    # downstream pool-CV QC; matches `^Pool` (case-insensitive) on column
    # names regardless of whether sample_filter$exclude_patterns also drops
    # them.
    pool_col_names <- sample_cols[grepl("^Pool", sample_cols,
                                          ignore.case = TRUE)]

    # Optionally filter to specific samples via metadata or config
    exclude_patterns <- cfg$sample_filter$exclude_patterns %||% NULL
    if (!is.null(exclude_patterns)) {
        for (pat in exclude_patterns) {
            sample_cols <- sample_cols[!grepl(pat, sample_cols, ignore.case = TRUE)]
        }
    }

    # If metadata is provided, restrict to matching sample columns
    if (!is.null(meta)) {
        sample_col_name <- cfg$effects$samples %||% "sample_id"
        if (sample_col_name %in% colnames(meta)) {
            meta_ids <- as.character(meta[[sample_col_name]])
            matched <- intersect(sample_cols, meta_ids)
            if (length(matched) > 0) {
                sample_cols <- matched
            }
        }
    }

    # Build feature IDs
    feat_ids <- build_lipid_feature_ids(data_df, id_cfg)

    # Build expression matrix
    expr_df <- data_df[, sample_cols, drop = FALSE]
    expr_raw <- coerce_df_to_numeric_matrix(expr_df, rownames_vec = feat_ids,
                                             name = "lipidsearch_expr")

    # Build pool matrix with the same row identity as expr_raw, suppressing
    # warnings from coerce_df_to_numeric_matrix when fewer than 2 pool cols.
    pool_matrix <- NULL
    if (length(pool_col_names) > 0) {
        pool_df <- data_df[, pool_col_names, drop = FALSE]
        pool_matrix <- tryCatch(
            coerce_df_to_numeric_matrix(pool_df, rownames_vec = feat_ids,
                                         name = "lipidsearch_pool"),
            error = function(e) {
                warning("parse_lipidsearch: pool matrix coercion failed: ",
                        e$message)
                NULL
            }
        )
    }

    # Row data (annotations) — include lipid class info
    row_data <- data_df[, annot_present, drop = FALSE]
    row_data$feature_id <- feat_ids

    # Standardize lipid class column name
    if ("ClassKey" %in% colnames(row_data) && !"lipid_class" %in% colnames(row_data)) {
        row_data$lipid_class <- as.character(row_data$ClassKey)
    }
    if ("SubClassKey" %in% colnames(row_data) && !"lipid_subclass" %in% colnames(row_data)) {
        row_data$lipid_subclass <- as.character(row_data$SubClassKey)
    }

    list(
        expr_raw    = expr_raw,
        row_data    = row_data,
        sample_ids  = sample_cols,
        pool_matrix = pool_matrix
    )
}


#' Parse translated CSV -> expr_raw + row_data
#'
#' Translated CSV has: Metabolite column, then sample columns.
#' Optionally has a Label row (row 1 = group assignments).
#' Lipid names follow shorthand notation (e.g., "Cer 18:1;O2/16:0").
#'
#' @param data_df  Raw data.frame.
#' @param cfg      lipidomics mode config.
#' @param meta     Metadata data.frame (optional).
#' @return list(expr_raw, row_data, sample_ids)
parse_translated_csv <- function(data_df, cfg, meta = NULL) {
    id_cfg <- cfg$id_columns %||% list()
    name_col <- id_cfg$feature_id_col %||% id_cfg$name_col %||% "Metabolite"

    if (!name_col %in% colnames(data_df)) {
        stop("Translated CSV: expected column '", name_col,
             "' not found. Available: ", paste(head(colnames(data_df), 10), collapse = ", "))
    }

    all_cols <- colnames(data_df)
    sample_cols <- setdiff(all_cols, name_col)

    # Remove empty rows
    data_df <- data_df[!is.na(data_df[[name_col]]) & nzchar(data_df[[name_col]]), ]

    feat_ids <- make.unique(as.character(data_df[[name_col]]), sep = "_")

    # Build expression matrix
    expr_df <- data_df[, sample_cols, drop = FALSE]
    expr_raw <- coerce_df_to_numeric_matrix(expr_df, rownames_vec = feat_ids,
                                             name = "translated_csv_expr")

    # Row data: parse lipid class from name
    row_data <- data.frame(
        feature_id = feat_ids,
        Name       = as.character(data_df[[name_col]]),
        stringsAsFactors = FALSE
    )
    row_data$lipid_class <- parse_lipid_class(row_data$Name)

    list(
        expr_raw   = expr_raw,
        row_data   = row_data,
        sample_ids = sample_cols
    )
}


#' Parse processed wide table (reuse metabolomics pattern)
#'
#' @param data_df Processed data.frame (wide).
#' @param cfg     lipidomics mode config.
#' @param meta    Metadata data.frame (REQUIRED).
#' @return list(expr_raw, row_data, sample_ids)
parse_lipid_processed_wide <- function(data_df, cfg, meta) {
    if (is.null(meta)) {
        stop("processed_wide: metadata is required (meta = NULL).")
    }

    sample_col <- cfg$effects$samples %||% "sample_id"
    if (!sample_col %in% colnames(meta)) {
        stop("processed_wide: metadata is missing column '", sample_col, "'.")
    }

    sample_ids <- as.character(meta[[sample_col]])
    sample_ids <- sample_ids[nzchar(sample_ids)]

    df_cols <- colnames(data_df)
    sample_cols <- intersect(df_cols, sample_ids)

    if (length(sample_cols) == 0) {
        stop("processed_wide: none of the metadata sample IDs were found in data columns.")
    }

    missing <- setdiff(sample_ids, df_cols)
    if (length(missing) > 0) {
        stop("processed_wide: samples in metadata missing from data: ",
             paste(missing, collapse = ", "))
    }

    id_cfg <- cfg$id_columns %||% list()
    feat_ids <- build_lipid_feature_ids(data_df, id_cfg)

    expr_df <- data_df[, sample_cols, drop = FALSE]
    expr_raw <- coerce_df_to_numeric_matrix(
        expr_df, rownames_vec = feat_ids, name = "processed_wide_expr"
    )

    row_data <- data_df[, setdiff(df_cols, sample_cols), drop = FALSE]
    row_data$feature_id <- feat_ids

    # Standardize lipid class column
    if (!"lipid_class" %in% colnames(row_data)) {
        if ("ClassKey" %in% colnames(row_data)) {
            row_data$lipid_class <- as.character(row_data$ClassKey)
        } else if ("Name" %in% colnames(row_data)) {
            row_data$lipid_class <- parse_lipid_class(row_data$Name)
        } else {
            row_data$lipid_class <- parse_lipid_class(feat_ids)
        }
    }

    list(
        expr_raw   = expr_raw,
        row_data   = row_data,
        sample_ids = sample_cols
    )
}


# ---- helpers ----------------------------------------------------------------

#' Build feature IDs for lipidomics data
#'
#' Priority: feature_id_col > Molecule > name_col > construct from name+mz+RT
build_lipid_feature_ids <- function(data_df, id_cfg) {
    fid_col <- id_cfg$feature_id_col

    if (!is.null(fid_col) && fid_col %in% colnames(data_df)) {
        ids <- as.character(data_df[[fid_col]])
        if (anyDuplicated(ids) > 0) {
            warning("Duplicate feature IDs in column '", fid_col, "'; making unique.")
            ids <- make.unique(ids, sep = "_")
        }
        return(ids)
    }

    # Try Molecule column (LipidSearch standard)
    if ("Molecule" %in% colnames(data_df)) {
        ids <- as.character(data_df[["Molecule"]])
        ids[is.na(ids) | ids == ""] <- paste0("unknown_", which(is.na(ids) | ids == ""))
        return(make.unique(ids, sep = "_"))
    }

    # Try name_col
    name_col <- id_cfg$name_col %||% "Metabolite"
    if (name_col %in% colnames(data_df)) {
        ids <- as.character(data_df[[name_col]])
        ids[is.na(ids) | ids == ""] <- paste0("unknown_", which(is.na(ids) | ids == ""))
        return(make.unique(ids, sep = "_"))
    }

    # Fallback: construct from Name + mz + RT
    mz_col <- id_cfg$mz_col %||% "Precursor Mz"
    rt_col <- id_cfg$rt_col %||% "Average Measured Retention Time"

    parts <- paste0("lipid_", seq_len(nrow(data_df)))
    if (mz_col %in% colnames(data_df)) {
        mz <- round(as.numeric(data_df[[mz_col]]), 4)
        parts <- paste0(parts, "_mz", mz)
    }
    if (rt_col %in% colnames(data_df)) {
        rt <- round(as.numeric(data_df[[rt_col]]), 2)
        parts <- paste0(parts, "_rt", rt)
    }

    make.unique(parts, sep = "_")
}


#' Detect lipid bond type from lipid name
#'
#' Classifies the bond linkage type based on class-name suffixes
#' (e.g. "PC-O", "PE-P") and chain-level prefixes (e.g. "(O-", "(P-").
#'
#' @param lipid_names Character vector of lipid names.
#' @return Character vector: "Acyl", "Alkyl", or "Alkenyl".
detect_lipid_bond_type <- function(lipid_names) {
    bond_type <- rep("Acyl", length(lipid_names))
    # Alkyl: class names ending in -O or chains with O- prefix
    bond_type[grepl("-O\\b|\\(O-", lipid_names)] <- "Alkyl"
    # Alkenyl: class names ending in -P or chains with P- prefix
    bond_type[grepl("-P\\b|\\(P-", lipid_names)] <- "Alkenyl"
    bond_type
}


#' Parse lipid class from lipid name string
#'
#' Extracts the headgroup/class from shorthand notation:
#'   "Cer 18:1;O2/16:0" -> "Cer"
#'   "PC 34:1"           -> "PC"
#'   "TG 52:3"           -> "TG"
#'   "SM d18:1/16:0"     -> "SM"
#'
#' @param names Character vector of lipid names.
#' @return Character vector of lipid classes.
parse_lipid_class <- function(names) {
    # Extract everything before the first space or parenthesis
    classes <- sub("^([A-Za-z0-9_-]+)\\s.*", "\\1", names)
    classes <- sub("^([A-Za-z0-9_-]+)\\(.*", "\\1", classes)
    classes[is.na(classes) | classes == ""] <- "Unknown"
    classes
}


#' Parse fatty acyl chain information from lipid name
#'
#' Extracts total carbon count and double bond count.
#'
#' @param names Character vector of lipid names.
#' @return data.frame with columns: total_carbons, total_double_bonds
parse_lipid_chains <- function(names) {
    # Match patterns like "34:1", "16:0/18:1", "d18:1/16:0"
    # Sum up all chain:db pairs found
    result <- data.frame(
        total_carbons      = NA_integer_,
        total_double_bonds = NA_integer_,
        stringsAsFactors   = FALSE
    )
    result <- result[rep(1, length(names)), ]
    rownames(result) <- NULL

    for (i in seq_along(names)) {
        nm <- names[i]
        # Find all C:DB patterns (digits:digits)
        matches <- gregexpr("(\\d+):(\\d+)", nm)[[1]]
        if (matches[1] == -1) next

        chain_strs <- regmatches(nm, list(matches))[[1]]
        carbons <- 0L
        dbs <- 0L
        for (cs in chain_strs) {
            parts <- strsplit(cs, ":")[[1]]
            carbons <- carbons + as.integer(parts[1])
            dbs <- dbs + as.integer(parts[2])
        }
        result$total_carbons[i] <- carbons
        result$total_double_bonds[i] <- dbs
    }

    result
}
