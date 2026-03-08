# R/domain/metabolomics/00_inputs.R
#
# Loading and parsing metabolomics inputs.
# Supports:
#   - "cd_raw"          : Compound Discoverer wide export (Area: columns)
#   - "processed_wide"  : Already-processed wide table (clean sample columns)
#   - "multi_level"     : Batch-load multiple Level_*.xlsx from a directory;
#                         each file is parsed with a shared per-file format
#                         (cfg$input[["level_format"]]: "cd_raw" | "processed_wide")
#                         and merged into a single expr_raw / row_data contract.
#
# Reuses: resolve_raw_path, read_table_auto, check_has_cols, check_all_in,
#         coerce_df_to_numeric_matrix, assert_numeric_matrix, assert_meta_contract
# NOTE on strict indexing:
#   All config key accesses that could be subject to R's partial-matching rules
#   (e.g. files[["data"]] vs files[["data_dir"]]) use [[ ]] with an explicit
#   is.null() guard.  Never use $ for these keys.


# ---- public entry point -----------------------------------------------------

#' Load metabolomics inputs from config
#'
#' For single-file formats ("cd_raw", "processed_wide", "long"),
#' \code{inp$data} is a single data.frame — identical to the previous behaviour.
#'
#' For the "multi_level" format, \code{inp$data} is a named list of data.frames,
#' one per level file (e.g., \code{list(Level_1 = df, Level_2 = df, ...)}).
#' \code{mod_met_raw()} remains the single dispatch point that handles the
#' structural difference.
#'
#' @param config Full pipeline config (list).
#' @return list with: data (data.frame or named list), metadata (data.frame or NULL),
#'         format
load_metabolomics_inputs <- function(config) {
    cfg <- config$modes$metabolomics
    if (is.null(cfg)) stop("No config for mode metabolomics")

    files <- cfg$files
    fmt   <- cfg$input[["format"]] %||% "cd_raw"

    # ── multi_level: read directory of level files ───────────────────────────
    if (fmt == "multi_level") {
        data_dir <- files[["data_dir"]]
        if (is.null(data_dir) || !nzchar(data_dir))
            stop("metabolomics: files$data_dir is required when format = 'multi_level'")

        abs_dir <- resolve_raw_path(config, data_dir)
        if (!dir.exists(abs_dir))
            stop("Metabolomics data directory not found: ", abs_dir)

        pattern   <- cfg$input[["level_pattern"]] %||% "\\.xlsx$"
        data_list <- read_multi_level_dir(abs_dir, pattern = pattern,
                                          sheet = cfg$input[["sheet"]])

        meta       <- .load_optional_metadata(config, files)
        sample_map <- .load_optional_sample_map(config, files)

        return(list(
            data       = data_list,   # named list: level_name -> data.frame
            metadata   = meta,
            sample_map = sample_map,
            format     = fmt
        ))
    }

    # ── single-file path (cd_raw / processed_wide / long) ───────────────────
    data_path <- files[["data"]]
    if (is.null(data_path) || !nzchar(data_path))
        stop("metabolomics config$files$data is required")

    abs_data <- resolve_raw_path(config, data_path)
    if (!file.exists(abs_data)) stop("Metabolomics data file not found: ", abs_data)

    data_df <- read_metab_file(abs_data, sheet = cfg$input[["sheet"]])

    meta       <- .load_optional_metadata(config, files)
    sample_map <- .load_optional_sample_map(config, files)

    # Handle group-row format: row 1 of data contains condition assignments
    # for each sample column (annotation columns are empty/NA in that row).
    if (isTRUE(cfg$input$has_group_row) && nrow(data_df) > 0) {
        group_row <- data_df[1, , drop = FALSE]
        data_df   <- data_df[-1, , drop = FALSE]
        rownames(data_df) <- NULL

        # Build metadata from group row if no external metadata provided
        if (is.null(meta)) {
            id_cfg      <- cfg$id_columns %||% list()
            annot_cols  <- id_cfg$annotation_cols %||% character(0)
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
        data       = data_df,
        metadata   = meta,
        sample_map = sample_map,
        format     = fmt
    )
}

# Internal: read optional metadata file (shared by single-file and multi_level paths)
.load_optional_metadata <- function(config, files) {
    meta_path <- files[["metadata"]]
    if (is.null(meta_path) || !nzchar(meta_path)) return(NULL)
    abs_meta <- resolve_raw_path(config, meta_path)
    if (!file.exists(abs_meta)) stop("Metadata file not found: ", abs_meta)
    read_table_auto(abs_meta)
}

# Internal: read optional sample_map file (shared by single-file and multi_level paths)
.load_optional_sample_map <- function(config, files) {
    sm_path <- files[["sample_map"]]
    if (is.null(sm_path) || !nzchar(sm_path)) return(NULL)
    abs_sm <- resolve_raw_path(config, sm_path)
    if (!file.exists(abs_sm)) stop("Sample map file not found: ", abs_sm)
    read_table_auto(abs_sm)
}


#' Validate metabolomics config (called by validate_config dispatch)
#' Uses [[ ]] strict indexing throughout to prevent R's partial-matching from
#' resolving files[["data"]] to files[["data_dir"]] (or vice versa) when only
#' one of the two keys is present.
validate_metabolomics_config <- function(cfg) {
    assert_one_of(cfg$input[["format"]], "input$format",
                  c("cd_raw", "processed_wide", "long", "multi_level"))

    fmt          <- cfg$input[["format"]] %||% "cd_raw"
    data_val     <- cfg$files[["data"]]
    data_dir_val <- cfg$files[["data_dir"]]

    has_data     <- !is.null(data_val)     && nzchar(data_val)
    has_data_dir <- !is.null(data_dir_val) && nzchar(data_dir_val)
    # Mutual exclusivity: files$data and files$data_dir cannot both be set
    if (has_data && has_data_dir) {
        stop("metabolomics: files$data and files$data_dir are mutually exclusive; ",
             "set only one depending on format ('multi_level' uses data_dir, ",
             "all others use data)")
    }

    if (fmt == "multi_level") {
        if (!has_data_dir)
            stop("metabolomics: files$data_dir is required when format = 'multi_level'")
        assert_one_of(cfg$input[["level_format"]], "input$level_format",
                      c("cd_raw", "processed_wide"))
    } else {
        if (!has_data)
            stop("metabolomics: files$data is required")
    }

    # If sample_map is provided, map_from and map_to must also be set
    if (!is.null(cfg$files[["sample_map"]]) && nzchar(cfg$files[["sample_map"]] %||% "")) {
        if (is.null(cfg$id_columns$map_from) || is.null(cfg$id_columns$map_to))
            stop("metabolomics: id_columns$map_from and id_columns$map_to are required ",
                 "when files$sample_map is set")
    }

    norm <- cfg$normalization
    if (!is.null(norm)) {
        assert_one_of(norm$sample_norm, "normalization$sample_norm",
                      c("none", "sum", "median", "pqn", "is"),
                      allow_null = TRUE)
        assert_one_of(norm$transform, "normalization$transform",
                      c("none", "log2", "log10"),
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

    df_cols     <- colnames(data_df)
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
    orig_id  <- attr(feat_ids, "original_id")

    expr_df  <- data_df[, sample_cols, drop = FALSE]
    expr_raw <- coerce_df_to_numeric_matrix(
        expr_df,
        rownames_vec = feat_ids,
        name = "processed_wide_expr"
    )

    # Everything else = annotation
    row_data <- data_df[, setdiff(df_cols, sample_cols), drop = FALSE]
    row_data$feature_id  <- feat_ids
    row_data$original_id <- orig_id

    list(
        expr_raw   = expr_raw,
        row_data   = row_data,
        sample_ids = sample_cols
    )
}


# ---- multi_level helpers ----------------------------------------------------

#' Enumerate and read all level files from a directory
#'
#' Files are sorted alphanumerically so Level_1 < Level_2 < Level_10.
#' Each file is read via \code{read_metab_file()}.
#'
#' @param dir_path Absolute path to the directory containing level files.
#' @param pattern  Regex to filter filenames (default: \code{"\\.xlsx$"}).
#' @param sheet    Excel sheet name/index (passed to \code{read_metab_file}).
#' @return Named list of \code{list(data_df, level_name, file_path)}, one per file.
#'         Names equal the level names (filename without extension).
read_multi_level_dir <- function(dir_path, pattern = "\\.xlsx$", sheet = NULL) {
    paths <- sort(list.files(dir_path, pattern = pattern,
                             full.names = TRUE, recursive = FALSE))
    if (length(paths) == 0)
        stop("No files matching pattern '", pattern, "' found in: ", dir_path)

    level_names <- tools::file_path_sans_ext(basename(paths))

    result <- lapply(seq_along(paths), function(i) {
        list(
            data_df    = read_metab_file(paths[i], sheet = sheet),
            level_name = level_names[i],
            file_path  = paths[i]
        )
    })
    names(result) <- level_names
    result
}


#' Parse a multi-level directory input into the canonical contract
#'
#' Dispatches each per-level data frame to \code{parse_cd_raw()} or
#' \code{parse_processed_wide()} (controlled by \code{cfg$input[["level_format"]]}),
#' validates that all levels share the same sample set, then calls
#' \code{merge_level_parsed()} to produce a single merged output whose structure
#' is identical to a single-file parse result.
#'
#' @param level_data_list Named list returned by \code{read_multi_level_dir()}.
#'        Names are used as level labels (e.g., \code{"Level_1"}).
#' @param cfg  metabolomics mode config.
#' @param meta Optional metadata data.frame (required for \code{"processed_wide"}).
#' @return \code{list(expr_raw, row_data, sample_ids, sample_map)} — identical
#'         contract to \code{parse_cd_raw()}.
parse_multi_level <- function(level_data_list, cfg, meta) {

    level_format <- cfg$input[["level_format"]]
    if (is.null(level_format) || !nzchar(level_format))
        stop("parse_multi_level: cfg$input$level_format is required")

    level_names <- names(level_data_list)

    # Parse each level file with the shared per-file format
    parsed_levels <- lapply(level_data_list, function(item) {
        switch(level_format,
            cd_raw         = parse_cd_raw(item$data_df, cfg),
            processed_wide = parse_processed_wide(item$data_df, cfg, meta),
            stop("parse_multi_level: unsupported level_format '", level_format, "'")
        )
    })
    names(parsed_levels) <- level_names

    # Validate that all levels expose the same sample set
    ref_ids <- sort(parsed_levels[[1]]$sample_ids)
    for (i in seq_along(parsed_levels)[-1]) {
        lv_ids <- sort(parsed_levels[[i]]$sample_ids)
        if (!identical(ref_ids, lv_ids)) {
            stop(sprintf(
                paste0("parse_multi_level: sample mismatch between '%s' and '%s'.\n",
                       "  '%s' samples: %s\n",
                       "  '%s' samples: %s"),
                level_names[1], level_names[i],
                level_names[1], paste(ref_ids, collapse = ", "),
                level_names[i], paste(lv_ids,  collapse = ", ")
            ))
        }
    }

    merge_level_parsed(parsed_levels, level_names)
}


#' Merge per-level parse results into a single contract-conformant output
#'
#' Data-contract guarantees:
#' \enumerate{
#'   \item \code{expr_raw} rows are prefixed with \code{"<Level>__"} to ensure
#'         cross-level uniqueness.  IDs are in \code{RT[rt]_MZ[mz]} format using
#'         the raw (unrounded) coordinates from \code{build_feature_ids()};
#'         no rounding is applied here.
#'   \item \code{row_data} column order is always:
#'         \code{c("feature_id", "Source_File", "feature_id_orig", <remaining annotation cols>)}.
#'         Remaining columns are the deterministic union of annotation columns
#'         across all levels (sorted by first appearance via \code{Reduce(union, ...)}).
#'   \item \code{feature_id_orig} is mandatory: it preserves the pre-prefix
#'         identifier so downstream code can always recover the raw source ID.
#'   \item Annotation columns absent in a given level are filled with \code{NA}.
#'   \item \code{sample_ids} is the reference column order from level 1 (all
#'         levels are reordered to match before merging).
#' }
#'
#' @param parsed_levels Named list of per-level parse results (each a list with
#'        keys \code{expr_raw}, \code{row_data}, \code{sample_ids}, \code{sample_map}).
#' @param level_names   Character vector of level labels, same length and order
#'        as \code{parsed_levels}.
#' @return \code{list(expr_raw, row_data, sample_ids, sample_map)}.
merge_level_parsed <- function(parsed_levels, level_names) {
    ref_ids <- parsed_levels[[1]]$sample_ids   # canonical column order (level 1)

    # Collect all annotation column names across all levels.
    # "feature_id" is included (it exists in every p$row_data); the fixed_cols
    # block below re-positions it at the front alongside the two new columns.
    all_annot_cols <- Reduce(union, lapply(parsed_levels, function(p) colnames(p$row_data)))

    # Deterministic column ordering:
    #   1. feature_id  (RT[rt]_MZ[mz], no prefix)
    #   2. level_id    (lowercase level label, e.g. "level_1")
    #   3. identification_level (integer parsed from file/level name) / original_id
    #   4. remaining annotation columns in union order
    fixed_cols      <- c("feature_id", "level_id", "identification_level", "original_id")
    remaining_cols  <- setdiff(all_annot_cols, fixed_cols)
    final_col_order <- c(fixed_cols, remaining_cols)

    expr_list <- vector("list", length(parsed_levels))
    rd_list   <- vector("list", length(parsed_levels))
    sm_list   <- vector("list", length(parsed_levels))

    for (i in seq_along(parsed_levels)) {
        p  <- parsed_levels[[i]]
        lv <- level_names[i]

        # ── Expression matrix ────────────────────────────────────────────────
        # Reorder columns to canonical ref_ids order; prefix rownames with the
        # level label to guarantee cross-level uniqueness after rbind.
        expr_mat                <- p$expr_raw[, ref_ids, drop = FALSE]
        rownames(expr_mat)      <- paste0(tolower(lv), "__", rownames(expr_mat))
        expr_list[[i]]          <- expr_mat

        # ── row_data ─────────────────────────────────────────────────────────
        rd                      <- p$row_data
        rd$level_id             <- tolower(lv)   # e.g. "level_1"
        rd$identification_level <- suppressWarnings(
            as.integer(sub("^Level_(\\d+)$", "\\1", lv, perl = TRUE))
        )

        # Fill annotation columns absent in this level
        for (col in setdiff(all_annot_cols, colnames(rd))) {
            rd[[col]] <- NA
        }

        # Apply deterministic column ordering
        rd_list[[i]] <- rd[, final_col_order, drop = FALSE]

        # ── sample_map (cd_raw only; NULL for processed_wide) ────────────────
        sm_list[[i]] <- p$sample_map
    }

    expr_raw <- do.call(rbind, expr_list)
    row_data <- do.call(rbind, rd_list)
    rownames(row_data) <- paste0(row_data$level_id, "__", row_data$feature_id)

    # ── Deduplication: one representative per feature_id ─────────────────────
    # Sort: lowest identification_level first (NA treated as Inf = worst),
    # then by original row position for deterministic tie-breaking.
    # After sorting, the first occurrence of each feature_id is "best".
    row_data$.sort_level <- ifelse(
        is.na(row_data$identification_level),
        Inf,
        as.numeric(row_data$identification_level)
    )
    row_data$.row_idx <- seq_len(nrow(row_data))

    o          <- order(row_data$feature_id,
                        row_data$.sort_level,
                        row_data$.row_idx)
    rd_sorted  <- row_data[o, , drop = FALSE]
    is_dup     <- duplicated(rd_sorted$feature_id)
    kept_rows  <- rd_sorted[!is_dup, , drop = FALSE]
    dropped    <- rd_sorted[ is_dup, , drop = FALSE]

    # ── Audit log with provenance columns ────────────────────────────────────
    dup_log <- NULL
    if (nrow(dropped) > 0) {
        kept_lookup <- data.frame(
            feature_id                = kept_rows$feature_id,
            kept_level_id             = kept_rows$level_id,
            kept_identification_level = kept_rows$identification_level,
            stringsAsFactors          = FALSE
        )
        dup_log <- merge(dropped, kept_lookup, by = "feature_id", all.x = TRUE)

        dup_log$drop_reason <- ifelse(
            is.na(dup_log$identification_level) &
                !is.na(dup_log$kept_identification_level),
            "missing_level_lost_to_numeric",
            ifelse(
                !is.na(dup_log$identification_level) &
                    !is.na(dup_log$kept_identification_level) &
                    dup_log$identification_level == dup_log$kept_identification_level,
                "tie_same_level_kept_first",
                "higher_identification_level"
            )
        )

        # Lead columns first, then remaining annotation columns, drop helpers
        log_lead <- c("feature_id", "identification_level", "original_id", "Name",
                      "drop_reason", "kept_level_id", "kept_identification_level")
        log_lead    <- intersect(log_lead, colnames(dup_log))
        log_rest    <- setdiff(colnames(dup_log),
                               c(log_lead, ".sort_level", ".row_idx"))
        dup_log     <- dup_log[, c(log_lead, log_rest), drop = FALSE]
    }

    # Strip helper columns; rebuild row_data and filter expr_raw
    kept_rows$.sort_level <- NULL
    kept_rows$.row_idx    <- NULL
    row_data              <- kept_rows
    rownames(row_data)    <- paste0(row_data$level_id, "__", row_data$feature_id)
    expr_raw              <- expr_raw[rownames(row_data), , drop = FALSE]

    n_kept    <- nrow(row_data)
    n_dropped <- if (is.null(dup_log)) 0L else nrow(dup_log)
    message("Deduplication complete: Kept ", n_kept, " unique features, removed ",
            n_dropped, " duplicates (see duplicate_log attribute for details).")

    # Deduplicated sample_map (identical across levels for cd_raw)
    non_null_maps <- Filter(Negate(is.null), sm_list)
    sample_map    <- if (length(non_null_maps) > 0) unique(do.call(rbind, non_null_maps)) else NULL

    out <- list(
        expr_raw        = expr_raw,
        row_data        = row_data,
        sample_ids      = ref_ids,
        sample_map      = sample_map,
        duplicate_log   = dup_log
    )
    attr(out, "duplicate_log") <- dup_log
    out
}


# ---- helpers ----------------------------------------------------------------

#' Build feature IDs from config rules
#'
#' Constructs IDs in the form \code{RT[rt]_MZ[mz]} using the raw (unrounded)
#' numeric values from \code{rt_col} and \code{mz_col}.  Per-row fallback to
#' \code{name_col} or a generic index when either coordinate is missing.
#'
#' When \code{feature_id_col} is present the unmodified source string is
#' attached as an attribute on the returned vector:
#' \describe{
#'   \item{\code{original_id}}{Character vector of the raw source strings.}
#' }
build_feature_ids <- function(data_df, id_cfg) {
    name_col <- id_cfg$name_col %||% "Name"
    mz_col   <- id_cfg$mz_col   %||% "m/z"
    rt_col   <- id_cfg$rt_col   %||% "RT [min]"
    fid_col  <- id_cfg$feature_id_col

    nr     <- nrow(data_df)
    has_mz <- mz_col %in% colnames(data_df)
    has_rt <- rt_col %in% colnames(data_df)
    has_nm <- name_col %in% colnames(data_df)

    # Vectorised RT[rt]_MZ[mz] builder; per-row fallback when a coordinate is NA
    make_rt_mz_ids <- function() {
        mz_vals <- if (has_mz) as.numeric(data_df[[mz_col]]) else rep(NA_real_, nr)
        rt_vals <- if (has_rt) as.numeric(data_df[[rt_col]]) else rep(NA_real_, nr)
        both_ok <- !is.na(mz_vals) & !is.na(rt_vals)

        fallback <- if (has_nm) {
            nm <- as.character(data_df[[name_col]])
            ifelse(is.na(nm) | nm == "", paste0("feature_", seq_len(nr)), nm)
        } else {
            paste0("feature_", seq_len(nr))
        }

        ifelse(both_ok,
               paste0("RT", as.character(rt_vals), "_MZ", as.character(mz_vals)),
               fallback)
    }

    if (!is.null(fid_col) && fid_col %in% colnames(data_df)) {
        original_id <- as.character(data_df[[fid_col]])
        ids <- make.unique(make_rt_mz_ids(), sep = "_dup")
        attr(ids, "original_id") <- original_id
        return(ids)
    }

    # Constructed path
    original_id <- if (has_nm) {
        nm <- as.character(data_df[[name_col]])
        ifelse(is.na(nm) | nm == "", paste0("feature_", seq_len(nr)), nm)
    } else {
        paste0("feature_", seq_len(nr))
    }

    ids <- make.unique(make_rt_mz_ids(), sep = "_dup")
    attr(ids, "original_id") <- original_id
    ids
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
