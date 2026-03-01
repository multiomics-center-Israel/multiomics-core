# R/domain/metabolomics/00_inputs.R
#
# Loading and parsing metabolomics inputs.
# Supports:
#   - "cd_raw"          : Compound Discoverer wide export (Area: columns)
#   - "processed_wide"  : Already-processed wide table (clean sample columns)
#
# Reuses: resolve_raw_path, read_table_auto, check_has_cols, check_all_in,
#         coerce_df_to_numeric_matrix, assert_numeric_matrix, assert_meta_contract


# ---- public entry point -----------------------------------------------------

#' Load metabolomics inputs from config
#'
#' @param config Full pipeline config (list).
#' @return list with: data (raw data.frame), metadata (data.frame or NULL),
#'         format, feature_id_col, sample_col
load_metabolomics_inputs <- function(config) {
    cfg <- config$modes$metabolomics
    if (is.null(cfg)) stop("No config for mode metabolomics")

    files <- cfg$files
    if (is.null(files$data)) stop("metabolomics config$files$data is required")

    abs_data <- resolve_raw_path(config, files$data)
    if (!file.exists(abs_data)) stop("Metabolomics data file not found: ", abs_data)

    data_df <- read_metab_file(abs_data, sheet = cfg$input$sheet)

    # Optional metadata
    meta <- NULL
    if (!is.null(files$metadata) && nzchar(files$metadata)) {
        abs_meta <- resolve_raw_path(config, files$metadata)
        if (!file.exists(abs_meta)) stop("Metadata file not found: ", abs_meta)
        meta <- read_table_auto(abs_meta)
    }

    # Handle group-row format: row 1 of data contains condition assignments
    # for each sample column (annotation columns are empty/NA in that row).
    if (isTRUE(cfg$input$has_group_row) && nrow(data_df) > 0) {
        group_row <- data_df[1, , drop = FALSE]
        data_df   <- data_df[-1, , drop = FALSE]
        rownames(data_df) <- NULL

        # Build metadata from group row if no external metadata provided
        if (is.null(meta)) {
            id_cfg     <- cfg$id_columns %||% list()
            annot_cols <- id_cfg$annotation_cols %||% character(0)
            sample_cols <- setdiff(colnames(data_df), annot_cols)

            cond_col <- cfg$de$condition_column %||%
                        cfg$effects$color %||% "condition"

            meta <- data.frame(
                sample_id = sample_cols,
                stringsAsFactors = FALSE
            )
            meta[[cond_col]] <- as.character(unlist(group_row[1, sample_cols]))
            message(sprintf(
                "metabolomics: built metadata from group row (%d samples, condition='%s')",
                length(sample_cols), cond_col
            ))
        }
    }

    list(
        data     = data_df,
        metadata = meta,
        format   = cfg$input$format %||% "cd_raw"
    )
}


#' Validate metabolomics config (called by validate_config dispatch)
validate_metabolomics_config <- function(cfg) {
    assert_one_of(cfg$input$format, "input$format",
                  c("cd_raw", "processed_wide", "long"))

    if (is.null(cfg$files$data) || !nzchar(cfg$files$data)) {
        stop("metabolomics: files$data is required")
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
    invisible(TRUE)
}


# ---- file reader ------------------------------------------------------------

#' Read metabolomics data file (xlsx or csv/tsv)
#' @param path Absolute file path.
#' @param sheet Sheet name or index for xlsx (optional).
#' @return data.frame
read_metab_file <- function(path, sheet = NULL) {
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

#' Parse Compound Discoverer raw export → expr_raw + row_data + meta
#'
#' @param data_df Raw CD data.frame.
#' @param cfg     metabolomics mode config.
#' @return list(expr_raw, row_data, sample_map, sample_ids)
parse_cd_raw <- function(data_df, cfg) {
    id_cfg   <- cfg$id_columns
    parse_cfg <- cfg$parsing %||% list()

    area_prefix <- parse_cfg$cd_area_prefix %||% "Area:"

    # 1) Identify Area: columns
    all_cols <- colnames(data_df)
    area_mask <- grepl(paste0("^", escapeRegex(area_prefix)), all_cols)
    area_cols <- all_cols[area_mask]

    if (length(area_cols) == 0) {
        stop(sprintf(
            "No columns matching prefix '%s' found. Available columns: %s",
            area_prefix, paste(head(all_cols, 20), collapse = ", ")
        ))
    }

    # 2) Parse sample IDs from Area: column names
    #    Default regex: "Area:\\s*(.+?)\\.raw\\s*\\(F\\d+\\)"
    #    Captures the sample name between "Area: " and ".raw (Fn)"
    sample_regex <- parse_cfg$cd_sample_regex %||%
        paste0(escapeRegex(area_prefix), "\\s*(.+?)\\.raw\\s*\\(F\\d+\\)")

    sample_ids <- vapply(area_cols, function(col) {
        m <- regmatches(col, regexec(sample_regex, col))[[1]]
        if (length(m) >= 2) m[2] else NA_character_
    }, character(1), USE.NAMES = FALSE)

    failed <- is.na(sample_ids)
    if (all(failed)) {
        stop("Could not parse any sample IDs from Area: columns using regex: ", sample_regex,
             "\nExample columns: ", paste(head(area_cols, 5), collapse = ", "))
    }
    if (any(failed)) {
        warning("Could not parse sample ID from: ",
                paste(area_cols[failed], collapse = ", "), " — dropping these columns.")
        area_cols  <- area_cols[!failed]
        sample_ids <- sample_ids[!failed]
    }

    # Check for duplicate sample IDs
    if (anyDuplicated(sample_ids) > 0) {
        dups <- unique(sample_ids[duplicated(sample_ids)])
        stop("Duplicate sample IDs parsed from Area: columns: ",
             paste(dups, collapse = ", "))
    }

    # 3) Build sample_map (original CD name → clean sample_id)
    sample_map <- data.frame(
        cd_column = area_cols,
        sample_id = sample_ids,
        stringsAsFactors = FALSE
    )

    # 4) Build expression matrix
    expr_df <- data_df[, area_cols, drop = FALSE]
    colnames(expr_df) <- sample_ids

    # 5) Build feature_id
    feat_ids <- build_feature_ids(data_df, id_cfg)

    expr_raw <- coerce_df_to_numeric_matrix(expr_df, rownames_vec = feat_ids,
                                             name = "cd_raw_expr")

    # 6) Row data (annotations)
    annot_cols <- setdiff(all_cols[!area_mask], character(0))
    row_data <- data_df[, annot_cols, drop = FALSE]
    row_data$feature_id <- feat_ids

    list(
        expr_raw    = expr_raw,
        row_data    = row_data,
        sample_map  = sample_map,
        sample_ids  = sample_ids
    )
}


#' Parse processed wide table → expr_raw + row_data (META-ONLY)
#'
#' This version relies ONLY on metadata to identify sample columns.
#' It does NOT use config heuristics / annotation_cols / numeric inference.
#'
#' Required:
#'   - meta must be provided
#'   - meta must contain the sample id column (cfg$effects$samples or "sample_id")
#'   - processed data must contain columns matching those sample IDs
#'
#' @param data_df Processed data.frame (wide).
#' @param cfg     metabolomics mode config.
#' @param meta    Metadata data.frame (REQUIRED).
#' @return list(expr_raw, row_data, sample_ids)
#' Parse processed wide table → expr_raw + row_data (META-ONLY, minimal)
parse_processed_wide <- function(data_df, cfg, meta) {
  
  if (is.null(meta)) {
    stop("processed_wide: metadata is required (meta = NULL).")
  }
  
  sample_col <- cfg$effects$samples %||% "sample_id"
  if (!sample_col %in% colnames(meta)) {
    stop("processed_wide: metadata is missing column '", sample_col, "'.")
  }
  
  sample_ids <- as.character(meta[[sample_col]])
  sample_ids <- sample_ids[nzchar(sample_ids)]
  
  if (anyDuplicated(sample_ids)) {
    stop(
      "processed_wide: duplicated sample IDs in metadata: ",
      paste(unique(sample_ids[duplicated(sample_ids)]), collapse = ", ")
    )
  }
  
  df_cols <- colnames(data_df)
  sample_cols <- intersect(df_cols, sample_ids)
  
  missing <- setdiff(sample_ids, df_cols)
  extra   <- setdiff(sample_cols, sample_ids)
  
  if (length(sample_cols) == 0) {
    stop(
      "processed_wide: none of the metadata sample IDs were found in data columns.\n",
      "Example metadata IDs: ", paste(head(sample_ids, 5), collapse = ", ")
    )
  }
  
  if (length(missing) > 0) {
    stop(
      "processed_wide: samples in metadata missing from data: ",
      paste(missing, collapse = ", ")
    )
  }
  
  # Build feature IDs
  feat_ids <- build_feature_ids(data_df, cfg$id_columns)
  
  expr_df <- data_df[, sample_cols, drop = FALSE]
  expr_raw <- coerce_df_to_numeric_matrix(
    expr_df,
    rownames_vec = feat_ids,
    name = "processed_wide_expr"
  )
  
  # Everything else = annotation
  row_data <- data_df[, setdiff(df_cols, sample_cols), drop = FALSE]
  row_data$feature_id <- feat_ids
  
  list(
    expr_raw   = expr_raw,
    row_data   = row_data,
    sample_ids = sample_cols
  )
}




# ---- helpers ----------------------------------------------------------------

#' Build feature IDs from config rules
#'
#' If feature_id_col is specified and exists, use it directly.
#' Otherwise construct from Name + mz + RT.
build_feature_ids <- function(data_df, id_cfg) {
    fid_col <- id_cfg$feature_id_col

    if (!is.null(fid_col) && fid_col %in% colnames(data_df)) {
        ids <- as.character(data_df[[fid_col]])
        if (anyDuplicated(ids) > 0) {
            warning("Duplicate feature IDs in column '", fid_col, "'; appending row index.")
            ids <- make.unique(ids, sep = "_")
        }
        return(ids)
    }

    # Construct from Name / mz / RT
    name_col <- id_cfg$name_col %||% "Name"
    mz_col   <- id_cfg$mz_col   %||% "m/z"
    rt_col   <- id_cfg$rt_col   %||% "RT [min]"

    parts <- character(nrow(data_df))

    if (name_col %in% colnames(data_df)) {
        nm <- as.character(data_df[[name_col]])
        nm[is.na(nm) | nm == ""] <- "Unknown"
        parts <- nm
    }

    if (mz_col %in% colnames(data_df)) {
        mz <- round(as.numeric(data_df[[mz_col]]), 4)
        parts <- paste0(parts, "_mz", mz)
    }

    if (rt_col %in% colnames(data_df)) {
        rt <- round(as.numeric(data_df[[rt_col]]), 2)
        parts <- paste0(parts, "_rt", rt)
    }

    if (all(parts == "" | parts == "Unknown")) {
        parts <- paste0("feature_", seq_len(nrow(data_df)))
    }

    make.unique(parts, sep = "_")
}




#' Build minimal metadata when no metadata file is provided
#'
#' Infers is_QC and is_blank flags from sample names.
build_minimal_meta <- function(sample_ids) {
    meta <- data.frame(
        sample_id = sample_ids,
        is_QC     = grepl("^QC", sample_ids, ignore.case = TRUE),
        is_blank  = grepl("^Blank", sample_ids, ignore.case = TRUE),
        stringsAsFactors = FALSE
    )

    # Create a simple group variable for QC plots
    meta$sample_type <- ifelse(meta$is_blank, "Blank",
                        ifelse(meta$is_QC, "QC", "Sample"))
    meta
}


#' Escape special regex characters in a string
escapeRegex <- function(string) {
    gsub("([\\[\\]\\{\\}\\(\\)\\^\\$\\.\\*\\+\\?\\|\\\\])", "\\\\\\1", string)
}
