# R/domain/metabolomics/06_enrichment.R
#
# Pathway enrichment analysis for metabolomics:
#   1. QEA  (Quantitative Enrichment Analysis) via globaltest       [self-contained]
#   2. ssGSEA via GSVA + Wilcoxon rank-sum test                     [competitive]
#   3. ORA  (Over-Representation Analysis) via Fisher's exact test  [competitive]
#   4. GSEA (Gene-Set Enrichment Analysis) via fgsea, ranked by log2FC [competitive]
#   5. Self-contained set test via limma::fry/mroast (rotation)     [self-contained]
#
# Operates on the standard pre-processing contract (expr_work, meta, row_data).
# ORA, GSEA and the self-contained test additionally require DE results (de_res).
#
# On competitive vs self-contained: the competitive methods (ORA/GSEA/ssGSEA)
# ask "is this set more changed than the rest of the features?" and lose power
# when a large fraction of features move (there is nothing left to be a contrast
# against). The self-contained methods (QEA, and the limma rotation test added
# here) ask "did this set change at all?", so they still report signal in that
# regime. The rotation test reuses the SAME limma model as DE (same design,
# same variance moderation), which QEA — run per-contrast without moderation —
# does not, so it stays powered at small n.
# Reuses: %||%


# ==== GMT PARSING =============================================================

#' Read a GMT file into a named list of gene/compound sets
#'
#' @param gmt_file       Path to GMT file.
#' @param include_descriptions Logical. If TRUE, also return pathway descriptions.
#' @return If include_descriptions=TRUE: list(sets, descriptions).
#'         Otherwise: named list of character vectors.
read_gmt_list <- function(gmt_file, include_descriptions = FALSE) {
    # Return empty containers when the file does not exist
    if (!file.exists(gmt_file)) {
        if (include_descriptions) return(list(sets = list(), descriptions = character(0)))
        return(list())
    }

    lines    <- readLines(gmt_file, warn = FALSE)
    gmt_list <- list()
    desc_map <- character(0)

    # Each GMT line is: pathway_id <TAB> description <TAB> member1 <TAB> ...
    for (line in lines) {
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) < 3) next

        pathway_id   <- parts[1]
        pathway_desc <- parts[2]
        members      <- parts[3:length(parts)]
        # Drop empty strings that can appear from trailing tabs
        members      <- members[nzchar(members)]

        if (length(members) > 0) {
            gmt_list[[pathway_id]] <- members
            if (nzchar(pathway_desc)) {
                desc_map[[pathway_id]] <- pathway_desc
            }
        }
    }

    if (include_descriptions) {
        return(list(sets = gmt_list, descriptions = desc_map))
    }
    gmt_list
}


#' Create "ID - Description" labels from pathway IDs and description map
#'
#' Concatenates each pathway ID with its human-readable description from
#' desc_map. Falls back to the bare ID when no description exists.
#'
#' @param ids      Character vector of pathway IDs.
#' @param desc_map Named character vector mapping IDs to descriptions.
#' @return Character vector of formatted labels, same length as ids.
make_pathway_labels <- function(ids, desc_map) {
    vapply(ids, function(id) {
        desc <- desc_map[[id]]
        if (!is.null(desc) && nzchar(desc)) paste0(id, " - ", desc) else id
    }, character(1), USE.NAMES = FALSE)
}


# ==== COMPOUND ID MAPPING ====================================================

#' Read an HMDB → KEGG mapping table from a TSV
#'
#' Returns a named character vector keyed by HMDB. NULL when the path is
#' missing/invalid or the required \code{HMDB}/\code{KEGG} columns are absent,
#' so callers can skip the mapping step without branching on file existence.
#'
#' @param path Path to a tab-separated mapping file.
#' @return Named character vector (KEGG values, HMDB names), or NULL.
read_hmdb_kegg_map <- function(path) {
    if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
    df <- tryCatch(
        readr::read_delim(path, delim = "\t",
                          col_types = readr::cols(), show_col_types = FALSE),
        error = function(e) NULL
    )
    if (is.null(df) || !all(c("HMDB", "KEGG") %in% colnames(df))) return(NULL)
    stats::setNames(as.character(df$KEGG), as.character(df$HMDB))
}

#' Map feature IDs to compound identifiers for enrichment
#'
#' Extracts compound names from row_data and optionally maps HMDB to KEGG.
#'
#' @param row_data     data.frame with feature annotations.
#' @param expr_mat     Expression matrix (features x samples).
#' @param mapping_file Optional path to HMDB-to-KEGG mapping TSV.
#' @param annotated_only Logical (default TRUE). When the data carries KEGG
#'   compound IDs, restrict the universe to the KEGG-annotated features so the
#'   background matches a KEGG-based GMT's namespace (MetaboAnalyst-style). It is
#'   a no-op when no KEGG ids are present, so name-only datasets are unaffected.
#'   Not exposed to config on purpose — this is the fixed, correct default.
#' @return list(expr_mapped, compound_names, feature_map) with remapped rownames.
map_compounds_for_enrichment <- function(row_data, expr_mat, mapping_file = NULL,
                                         annotated_only = TRUE) {
    # Align the expression matrix to the annotation rows. row_data can be a
    # filtered subset of expr_mat (e.g. missingness filtering drops features, so
    # expr_raw ends up with more rows than row_data); matching by feature id
    # keeps compound_names (taken from row_data) and the matrix rows in lockstep.
    # Without this, the logical filtering below indexes a mismatched matrix and
    # errors, which the module's tryCatch turns into a silent NULL for every
    # enrichment method.
    if (!is.null(row_data) && "feature_id" %in% colnames(row_data) &&
        !is.null(rownames(expr_mat))) {
        common <- intersect(as.character(row_data$feature_id), rownames(expr_mat))
        row_data <- row_data[match(common, as.character(row_data$feature_id)), ,
                             drop = FALSE]
        expr_mat <- expr_mat[common, , drop = FALSE]
    }

    # Extract compound names from row_data
    if ("Name" %in% colnames(row_data)) {
        compound_names <- as.character(row_data$Name)
    } else if ("Molecule" %in% colnames(row_data)) {
        compound_names <- as.character(row_data$Molecule)
    } else {
        compound_names <- rownames(expr_mat)
    }

    # HMDB-to-KEGG mapping: replace compound names with KEGG IDs when a
    # mapping file is provided, so that IDs match KEGG-based GMT pathways.
    map_vec <- read_hmdb_kegg_map(mapping_file)
    if (!is.null(map_vec)) {
        # Locate HMDB IDs — try common column-name variants, then parse feature_id
        hmdb_ids <- NULL
        hmdb_col <- intersect(c("HMDB", "HMDB_ID", "HMDB ID", "hmdb_id"),
                              colnames(row_data))[1]
        if (!is.na(hmdb_col)) {
            hmdb_ids <- as.character(row_data[[hmdb_col]])
        } else if ("feature_id" %in% colnames(row_data)) {
            # Try extracting HMDB from feature_id (HMDB|Name format)
            hmdb_ids <- sub("\\|.*$", "", as.character(row_data$feature_id))
        }

        if (!is.null(hmdb_ids)) {
            mapped_kegg <- map_vec[hmdb_ids]
            has_map <- !is.na(mapped_kegg) & mapped_kegg != ""
            compound_names[has_map] <- mapped_kegg[has_map]
            message("enrichment: mapped ", sum(has_map), " features to KEGG IDs")
        }
    }

    # Prefer explicit KEGG compound IDs when the annotation carries them: a
    # KEGG column is unambiguous, whereas compound names are spelled
    # inconsistently across software. Applied AFTER the HMDB->KEGG step so a
    # direct KEGG ID wins over a lookup. Cells may hold "C00031", "cpd:C00031"
    # or "C00031;C00267" — take the first KEGG compound ID present.
    kegg_col <- intersect(c("KEGG", "KEGG_ID", "KEGG ID", "kegg", "kegg_id"),
                          colnames(row_data))[1]
    if (!is.na(kegg_col)) {
        kegg_raw <- as.character(row_data[[kegg_col]])
        has_kegg <- grepl("C[0-9]{5}", kegg_raw)
        kegg_id  <- sub(".*?(C[0-9]{5}).*", "\\1", kegg_raw, perl = TRUE)
        compound_names[has_kegg] <- kegg_id[has_kegg]
        message("enrichment: using KEGG IDs from '", kegg_col, "' for ",
                sum(has_kegg), " features")
    }

    # Optionally restrict the universe to features carrying a KEGG compound ID.
    # With a KEGG-based GMT, name-only features can never match a pathway, so
    # keeping them only inflates the background; MetaboAnalyst-style ORA uses the
    # KEGG-mappable set as the reference. Detect KEGG IDs by their Cxxxxx shape
    # (set above from the KEGG column or the HMDB->KEGG mapping).
    if (isTRUE(annotated_only)) {
        is_kegg <- grepl("^C[0-9]{5}$", trimws(compound_names))
        # No-op when the data carries no KEGG ids at all (e.g. a name-only
        # dataset with a name-based GMT), so the background is never emptied.
        if (any(is_kegg)) {
            n_drop <- sum(!is_kegg & !is.na(compound_names))
            compound_names[!is_kegg] <- NA_character_
            message("enrichment: restricting background to ", sum(is_kegg),
                    " KEGG-annotated features (dropped ", n_drop, " name-only).")
        }
    }

    # Per-feature id -> compound map, keyed by the feature ids callers hold
    # (ORA/GSEA look up by feature_id). compound_names is row-aligned to
    # row_data, so key on row_data$feature_id when available (falling back to
    # the matrix rownames). ORA/GSEA use this to translate their feature ids to
    # compounds; deriving that positionally after the filtering/dedup below
    # would misalign, so build it here beforehand.
    feat_key <- if ("feature_id" %in% colnames(row_data) &&
                    length(row_data$feature_id) == length(compound_names)) {
        as.character(row_data$feature_id)
    } else {
        rownames(expr_mat)
    }
    feature_map <- stats::setNames(compound_names, feat_key)

    # Filter valid compound names
    valid <- !is.na(compound_names) & nzchar(compound_names) & compound_names != "NA"
    mat_use   <- expr_mat[valid, , drop = FALSE]
    names_use <- compound_names[valid]

    # Deduplicate (keep first occurrence)
    dup <- duplicated(names_use)
    mat_use <- mat_use[!dup, , drop = FALSE]
    rownames(mat_use) <- names_use[!dup]

    list(
        expr_mapped    = mat_use,
        compound_names = names_use[!dup],
        feature_map    = feature_map
    )
}


#' Translate HMDB IDs in a GMT list to KEGG using a mapping file
#'
#' Replaces HMDB compound identifiers within each pathway set with their
#' corresponding KEGG IDs. Compounds without a mapping are kept as-is.
#' Returns the original list unchanged if the mapping file is missing or
#' does not contain the required HMDB/KEGG columns.
#'
#' @param gmt_list     Named list of character vectors (pathway -> compound IDs).
#' @param mapping_file Path to a TSV with HMDB and KEGG columns (or NULL).
#' @return Named list with the same structure, HMDB IDs replaced where possible.
translate_gmt_hmdb_to_kegg <- function(gmt_list, mapping_file) {
    # No-op when mapping is absent or file does not exist
    if (is.null(mapping_file) || !file.exists(mapping_file)) return(gmt_list)

    # Read the HMDB-to-KEGG lookup table
    mapping_df <- readr::read_delim(mapping_file, delim = "\t",
                                    col_types = readr::cols(),
                                    show_col_types = FALSE)

    # Skip if the required columns are missing
    if (!("HMDB" %in% colnames(mapping_df) && "KEGG" %in% colnames(mapping_df))) {
        return(gmt_list)
    }

    # Build named vector for vectorised replacement in each pathway set
    hmdb2kegg <- stats::setNames(mapping_df$KEGG, mapping_df$HMDB)
    lapply(gmt_list, function(cpds) {
        mapped <- hmdb2kegg[cpds]
        # Keep original ID when no KEGG mapping is found
        ifelse(!is.na(mapped) & mapped != "", mapped, cpds)
    })
}


# ==== QEA VIA GLOBALTEST =====================================================

#' Run Quantitative Enrichment Analysis (QEA) using globaltest
#'
#' Orchestrates QEA across one or more pathway libraries provided as GMT files.
#' Each GMT file basename becomes a library name. Libraries listed in config
#' without a matching GMT file are warned about and skipped. Results from all
#' libraries are combined into a single table with FDR correction across all
#' pathways.
#'
#' @param pre     Preprocessing results (must include expr_work/expr_log, meta,
#'   row_data, info). QEA runs on the same analysis matrix as DE — see below.
#' @param config  Full pipeline config (reads modes$metabolomics$enrichment).
#' @return list(table, per_library_tables, method) or NULL if disabled/unavailable.
run_metabolomics_qea <- function(pre, config) {
    cfg      <- config$modes$metabolomics
    enr_cfg  <- cfg$enrichment %||% list()

    # Resolve gmt_file / mapping_file like the data files: absolute paths are
    # kept as-is, relative paths are located under paths$raw.
    enr_cfg$gmt_file     <- resolve_input_path(config, enr_cfg$gmt_file)
    enr_cfg$mapping_file <- resolve_input_path(config, enr_cfg$mapping_file)

    # Guard: skip when enrichment is disabled in config
    if (!isTRUE(enr_cfg$run_enrichment)) {
        message("metabolomics QEA: disabled in config — skipping")
        return(NULL)
    }

    # Guard: skip when globaltest is not installed
    if (!requireNamespace("globaltest", quietly = TRUE)) {
        message("metabolomics QEA: globaltest not available — skipping")
        return(NULL)
    }

    # Resolve which metadata columns define condition and sample identity,
    # cascading through enrichment > DE > effects > defaults.
    condition_col <- enr_cfg$condition_column %||%
                     cfg$de$condition_column %||%
                     cfg$effects$color %||% "sample_type"
    sample_col <- cfg$effects$samples %||% "sample_id"

    # ---- Determine libraries to run ----
    # Merge explicit library names from config with labels derived from GMT
    # file basenames. The gmt_lookup maps each label to its file path.
    libraries <- enr_cfg$libraries
    if (is.null(libraries)) {
        lib_single <- enr_cfg$library
        libraries <- if (!is.null(lib_single)) list(lib_single) else list()
    }

    gmt_files <- enr_cfg$gmt_file
    gmt_lookup <- character(0)
    if (!is.null(gmt_files)) {
        gmt_files <- unlist(gmt_files)
        for (gf in gmt_files) {
            if (file.exists(gf)) {
                gmt_label <- tools::file_path_sans_ext(basename(gf))
                libraries <- c(libraries, gmt_label)
                gmt_lookup[[gmt_label]] <- gf
            } else {
                warning("GMT file not found: ", gf)
            }
        }
    }

    mapping_file <- enr_cfg$mapping_file
    # standardize (globaltest::gt) rescales predictors so each metabolite gets
    # equal baseline weight. It changes pathway statistics, so treat it as a
    # sensitivity setting, not an automatic correction. Validate loudly — a
    # mistyped value must not silently coerce to FALSE.
    standardize <- (enr_cfg$qea %||% list())$standardize %||% FALSE
    if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize)) {
        stop("enrichment.qea.standardize must be TRUE or FALSE.")
    }

    # ---- Prepare QEA data ----
    # Map features to compound IDs and transpose into the samples-by-compounds
    # layout expected by globaltest. The result is written to a temp TSV.
    #
    # Matrix source: QEA now starts from the SAME upstream analysis matrix, scale,
    # and biological-sample set as DE — via .metab_de_matrix_condition() (the
    # helper the self-contained test uses) — instead of pre$expr_raw. This makes
    # the QEA-vs-DE comparison fair (globaltest is scale/parametrization
    # sensitive). Note: map_compounds_for_enrichment() then applies enrichment-
    # specific compound mapping, KEGG restriction and dedup, so the FINAL feature
    # set is not identical to DE's. Behavior change vs the old expr_raw path:
    # QC/blank samples excluded, values on the log-normalized scale.
    de_in  <- .metab_de_matrix_condition(pre, config)
    mapped <- map_compounds_for_enrichment(pre$row_data, de_in$mat, mapping_file)
    if (nrow(mapped$expr_mapped) < 3) {
        message("metabolomics QEA: too few compounds — skipping")
        return(NULL)
    }

    # Build globaltest input format: rows=samples, cols=compounds.
    # Prepend Sample and Group columns required by run_qea_gmt_internal().
    # Align conditions against the ALREADY-FILTERED metadata (de_in$meta), matched
    # by sample id (not position), and guard: a stale sample id or a bad
    # condition_column must fail loudly here, not surface as NAs deep inside
    # globaltest.
    mat_t      <- as.data.frame(t(mapped$expr_mapped), check.names = FALSE)
    sample_ids <- rownames(mat_t)
    idx        <- match(sample_ids, de_in$meta[[sample_col]])
    if (anyNA(idx)) {
        stop("metabolomics QEA: analysis-matrix samples could not be matched to ",
             "metadata via '", sample_col, "'.")
    }
    conditions <- droplevels(factor(de_in$meta[[condition_col]][idx]))
    if (anyNA(conditions)) {
        stop("metabolomics QEA: condition column '", condition_col,
             "' has missing values after sample alignment.")
    }
    if (nlevels(conditions) < 2L) {
        stop("metabolomics QEA: at least two condition levels are required ",
             "(column '", condition_col, "').")
    }
    df_t <- cbind(
        data.frame(Sample = sample_ids, Group = as.character(conditions),
                   stringsAsFactors = FALSE),
        mat_t
    )

    # Write to a temp directory so each library reads from the same file
    tmp_dir  <- tempfile("metabo_qea_")
    dir.create(tmp_dir, recursive = TRUE)
    data_file <- file.path(tmp_dir, "combined_data.txt")
    utils::write.table(df_t, data_file, sep = "\t", row.names = FALSE,
                        quote = FALSE)

    message("QEA data: ", nrow(df_t), " samples x ", ncol(df_t) - 2, " compounds")

    # ---- Warn about library names with no GMT file ----
    unresolved <- setdiff(libraries, names(gmt_lookup))
    if (length(unresolved) > 0) {
        warning("QEA: no GMT file for library name(s): ",
                paste(unresolved, collapse = ", "), " — skipping these")
        libraries <- intersect(libraries, names(gmt_lookup))
    }

    if (length(libraries) == 0) {
        message("metabolomics QEA: no libraries with GMT files — skipping")
        unlink(tmp_dir, recursive = TRUE)
        return(NULL)
    }

    # ---- Run QEA for each library ----
    # Dispatch globaltest independently per library; failures are caught so
    # one bad GMT does not abort the entire enrichment step.
    all_results <- list()

    for (lib in libraries) {
        result <- tryCatch(
            run_qea_gmt_internal(data_file, gmt_lookup[[lib]], mapping_file,
                                 standardize = standardize),
            error = function(e) {
                warning("QEA failed for ", lib, ": ", e$message)
                NULL
            }
        )
        if (!is.null(result)) all_results[[lib]] <- result
    }

    # Clean up temp data file
    unlink(tmp_dir, recursive = TRUE)

    if (length(all_results) == 0) {
        message("metabolomics QEA: no results from any library")
        return(NULL)
    }

    # Combine per-library tables and apply global FDR correction
    combined_df <- do.call(rbind, lapply(names(all_results), function(lib) {
        df <- all_results[[lib]]
        df$library <- lib
        df
    }))
    combined_df$FDR <- stats::p.adjust(combined_df$raw_p, method = "fdr")
    combined_df <- combined_df[order(combined_df$FDR), ]

    n_sig <- sum(combined_df$FDR < 0.05, na.rm = TRUE)
    message("QEA complete: ", nrow(combined_df), " pathways, ",
            n_sig, " with FDR < 0.05")

    list(
        table              = combined_df,
        per_library_tables = all_results,
        method             = "globaltest_qea"
    )
}


# ---- QEA internals ----------------------------------------------------------

#' Thin seam around globaltest::gt (isolates the standardize wiring)
#'
#' Exists so tests can assert that \code{standardize} is threaded through to
#' \code{globaltest::gt()} without depending on gt's numeric output.
#'
#' @param response   Response factor of sample groups.
#' @param X          Samples x compounds matrix.
#' @param subsets    Named list of compound sets (pathways).
#' @param standardize Logical, forwarded to \code{globaltest::gt()}.
#' @return A \code{gt.object}.
.run_globaltest <- function(response, X, subsets, standardize = FALSE) {
    globaltest::gt(response, X, subsets = subsets, standardize = standardize)
}

#' Run QEA via globaltest on a single GMT file
#'
#' Parses the GMT, optionally translates HMDB IDs to KEGG, filters to pathways
#' with at least 2 matching compounds in the data, and runs globaltest::gt().
#' Returns a data.frame of pathway-level p-values and hit counts.
#'
#' @param data_file    Path to the tab-delimited QEA data (Sample, Group, compounds).
#' @param gmt_file     Path to GMT file defining pathway sets.
#' @param mapping_file Path to HMDB-to-KEGG mapping TSV (or NULL).
#' @param standardize  Logical, passed to \code{globaltest::gt()}. When TRUE,
#'   each compound is standardized to equal baseline weight (recommended when the
#'   features' relative scales are arbitrary). Default FALSE preserves gt()'s
#'   default.
#' @return data.frame(pathway, raw_p, hits) or NULL if no testable pathways.
run_qea_gmt_internal <- function(data_file, gmt_file, mapping_file,
                                 standardize = FALSE) {
    if (!requireNamespace("globaltest", quietly = TRUE)) {
        stop("Package 'globaltest' required for GMT enrichment.")
    }

    # Parse pathway sets and their descriptions from the GMT file
    gmt_parsed <- read_gmt_list(gmt_file, include_descriptions = TRUE)
    gmt_list   <- gmt_parsed$sets
    desc_map   <- gmt_parsed$descriptions

    if (length(gmt_list) == 0) return(NULL)

    # Harmonise compound IDs: translate HMDB to KEGG so GMT members match data
    gmt_list <- translate_gmt_hmdb_to_kegg(gmt_list, mapping_file)

    # Attach human-readable descriptions to pathway names for output tables
    names(gmt_list) <- make_pathway_labels(names(gmt_list), desc_map)

    # Read the samples-by-compounds data written by run_metabolomics_qea()
    df <- utils::read.table(data_file, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, check.names = FALSE)
    response <- factor(df$Group)
    X <- as.matrix(df[, -c(1, 2), drop = FALSE])

    # Keep only pathways where at least 2 members are present in the data.
    # Match case-insensitively so "GLUCOSE" in data matches "Glucose" in GMT.
    avail_lower <- tolower(colnames(X))
    avail_map   <- stats::setNames(colnames(X), avail_lower)
    subsets <- lapply(gmt_list, function(cpds) {
        matched <- avail_map[tolower(cpds)]
        unname(matched[!is.na(matched)])
    })
    keep <- vapply(subsets, length, integer(1)) >= 2L
    subsets <- subsets[keep]

    if (length(subsets) == 0) {
        message("QEA GMT: no pathways with >= 2 matching compounds")
        return(NULL)
    }

    message("QEA GMT: ", length(subsets), " / ", length(gmt_list),
            " pathways with >= 2 matching compounds")

    # Run the global test: tests whether compound profiles in each pathway
    # are collectively associated with the condition factor. Wrapped in
    # .run_globaltest() so the standardize wiring is a single, testable seam.
    res_gt <- tryCatch(
        .run_globaltest(response, X, subsets = subsets, standardize = standardize),
        error = function(e) { warning("globaltest error: ", e$message); NULL }
    )
    if (is.null(res_gt)) return(NULL)

    # Extract the result table and return as a tidy data.frame
    res_tbl <- globaltest::result(res_gt)

    data.frame(
        pathway = rownames(res_tbl),
        raw_p   = res_tbl[, "p-value"],
        hits    = res_tbl[, "#Cov"],
        stringsAsFactors = FALSE
    )
}


# ==== ssGSEA VIA GSVA ========================================================

#' Run ssGSEA pathway enrichment via GSVA + Wilcoxon
#'
#' Computes single-sample GSEA scores for each pathway and sample, then
#' applies a Wilcoxon rank-sum test per pathway to compare two conditions.
#' Requires exactly two condition levels for statistical testing; if more
#' are present, scores are returned without p-values.
#'
#' @param pre     Preprocessing results (must include expr_work/expr_log, meta,
#'   row_data, info). Scores are computed on the same analysis matrix as DE.
#' @param config  Full pipeline config (reads modes$metabolomics$enrichment).
#' @return list(table, scores, method) or NULL if disabled/unavailable.
run_metabolomics_ssgsea <- function(pre, config) {
    cfg     <- config$modes$metabolomics
    enr_cfg <- cfg$enrichment %||% list()

    # Resolve gmt_file / mapping_file like the data files: absolute paths are
    # kept as-is, relative paths are located under paths$raw.
    enr_cfg$gmt_file     <- resolve_input_path(config, enr_cfg$gmt_file)
    enr_cfg$mapping_file <- resolve_input_path(config, enr_cfg$mapping_file)

    # Guard: skip when enrichment is disabled in config
    if (!isTRUE(enr_cfg$run_enrichment)) {
        message("metabolomics ssGSEA: disabled — skipping")
        return(NULL)
    }

    # Guard: skip when GSVA is not installed
    if (!requireNamespace("GSVA", quietly = TRUE)) {
        message("metabolomics ssGSEA: GSVA not available — skipping")
        return(NULL)
    }

    # Resolve condition and sample columns (same cascade as QEA)
    condition_col <- enr_cfg$condition_column %||%
                     cfg$de$condition_column %||%
                     cfg$effects$color %||% "sample_type"
    sample_col <- cfg$effects$samples %||% "sample_id"
    mapping_file <- enr_cfg$mapping_file

    # ---- Build expression matrix with compound IDs ----
    # Map features to enrichment-ready IDs (Name/HMDB/KEGG). Source matrix: the
    # DE analysis matrix via .metab_de_matrix_condition(), not pre$expr_raw.
    # Unlike QEA, ssGSEA is RANK-BASED within each sample, so a monotonic
    # transform (log) alone would barely move the scores — the alignment here is
    # for consistency with DE on sample inclusion (QC/blank excluded), the
    # normalization actually applied, and the feature universe, not for scale per
    # se. (GSVA's cross-sample score normalization also makes the sample set
    # matter.)
    de_in    <- .metab_de_matrix_condition(pre, config)
    mapped   <- map_compounds_for_enrichment(pre$row_data, de_in$mat, mapping_file)
    expr_mat <- as.matrix(mapped$expr_mapped)
    if (nrow(expr_mat) < 2) {
        message("metabolomics ssGSEA: too few features — skipping")
        return(NULL)
    }

    # ---- Load gene sets from GMT files ----
    # Parse each GMT, translate IDs, label pathways, and accumulate
    gene_sets <- list()
    gmt_files <- unlist(enr_cfg$gmt_file)
    if (!is.null(gmt_files)) {
        for (gf in gmt_files) {
            if (!file.exists(gf)) next
            gmt_parsed <- read_gmt_list(gf, include_descriptions = TRUE)
            gmt <- gmt_parsed$sets
            desc_map <- gmt_parsed$descriptions
            gmt <- translate_gmt_hmdb_to_kegg(gmt, mapping_file)
            names(gmt) <- make_pathway_labels(names(gmt), desc_map)
            gene_sets <- c(gene_sets, gmt)
            message("ssGSEA: loaded ", length(gmt), " sets from ", basename(gf))
        }
    }

    if (length(gene_sets) == 0) {
        message("metabolomics ssGSEA: no gene sets available — skipping")
        return(NULL)
    }

    # Discard pathway sets with fewer than 2 members present in the data
    # (case-insensitive match)
    avail_map <- stats::setNames(rownames(expr_mat), tolower(rownames(expr_mat)))
    gene_sets <- lapply(gene_sets, function(cpds) {
        m <- avail_map[tolower(cpds)]
        unname(m[!is.na(m)])
    })
    gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) >= 2L]

    if (length(gene_sets) == 0) {
        message("ssGSEA: no gene sets with >= 2 matching compounds")
        return(NULL)
    }
    message("ssGSEA: ", length(gene_sets), " gene sets with >= 2 members")

    # ---- Run ssGSEA ----
    # Try the new ssgseaParam API first (GSVA >= 1.50); fall back to the
    # legacy gsva() interface for older versions.
    scores <- tryCatch({
        param <- GSVA::ssgseaParam(expr_mat, gene_sets)
        GSVA::gsva(param)
    }, error = function(e1) {
        tryCatch(
            GSVA::gsva(expr_mat, gene_sets, method = "ssgsea", verbose = FALSE),
            error = function(e2) { warning("ssGSEA failed: ", e2$message); NULL }
        )
    })

    if (is.null(scores) || nrow(scores) == 0) {
        message("ssGSEA: no results")
        return(NULL)
    }
    message("ssGSEA: computed scores for ", nrow(scores), " pathways x ",
            ncol(scores), " samples")

    # ---- Wilcoxon test per pathway ----
    # Align condition labels to score matrix columns via sample IDs
    meta <- pre$meta
    conditions <- factor(meta[[condition_col]][
        match(colnames(scores), meta[[sample_col]])
    ])
    cond_levels <- levels(conditions)

    # Wilcoxon rank-sum requires exactly two groups; return raw scores otherwise
    if (length(cond_levels) != 2) {
        message("ssGSEA: Wilcoxon test requires exactly 2 conditions, found ",
                length(cond_levels))
        return(list(table = NULL, scores = scores, method = "ssgsea"))
    }

    # Split sample indices by condition for the two-group comparison
    grp1_idx <- which(conditions == cond_levels[1])
    grp2_idx <- which(conditions == cond_levels[2])

    # Test each pathway independently; NA on failure (e.g. tied ranks)
    pvalues <- vapply(seq_len(nrow(scores)), function(i) {
        tryCatch(
            stats::wilcox.test(scores[i, grp1_idx], scores[i, grp2_idx])$p.value,
            error = function(e) NA_real_
        )
    }, numeric(1))

    # FDR correction across all pathways and per-group mean scores
    fdr   <- stats::p.adjust(pvalues, method = "fdr")
    mean1 <- rowMeans(scores[, grp1_idx, drop = FALSE], na.rm = TRUE)
    mean2 <- rowMeans(scores[, grp2_idx, drop = FALSE], na.rm = TRUE)

    # Assemble result table with significance flag and per-condition means
    results_df <- data.frame(
        pathway    = rownames(scores),
        p_value    = pvalues,
        FDR        = fdr,
        score_diff = mean2 - mean1,
        significant = !is.na(fdr) & fdr < 0.05,
        stringsAsFactors = FALSE
    )
    colnames(results_df)[colnames(results_df) == "score_diff"] <- "score_diff"
    results_df[[paste0("mean_", cond_levels[1])]] <- mean1
    results_df[[paste0("mean_", cond_levels[2])]] <- mean2
    results_df <- results_df[order(results_df$p_value), ]

    n_sig <- sum(results_df$significant, na.rm = TRUE)
    message("ssGSEA: ", nrow(results_df), " pathways, ", n_sig,
            " significant (FDR < 0.05)")

    list(
        table  = results_df,
        scores = scores,
        method = "ssgsea_wilcoxon"
    )
}


# ==== CONTRAST SELECTION (shared by ORA / GSEA) ==============================

#' Resolve which DE contrast an enrichment method should use
#'
#' Returns \code{contrast} when it names an existing \code{de_res$de_tables}
#' entry; otherwise falls back to the first contrast (warning when a non-NULL
#' name was not found). Callers guarantee \code{de_res$de_tables} is non-empty.
#'
#' @param contrast Requested contrast name, or NULL for the first.
#' @param de_res   DE results carrying a named \code{de_tables} list.
#' @param method   Short method label used only in the warning message.
#' @return A single contrast name (a key of \code{de_res$de_tables}).
.resolve_enrichment_contrast <- function(contrast, de_res, method) {
    available <- names(de_res$de_tables)
    if (is.null(contrast)) return(available[1])
    if (contrast %in% available) return(contrast)
    warning("metabolomics ", method, ": contrast '", contrast,
            "' not found; using '", available[1], "'.")
    available[1]
}


# ==== ORA VIA FISHER'S EXACT TEST =============================================

#' Run Over-Representation Analysis (ORA) using Fisher's exact test
#'
#' Tests whether significant features from DE are over-represented in each
#' pathway defined by GMT files.  The foreground is the set of significant
#' features (adj.P.Val < p_cutoff AND |logFC| >= log2fc_cutoff); the
#' background is the full set of measured features.
#'
#' @param pre     Preprocessing results (expr_raw, meta, row_data).
#' @param de_res  DE results from run_metabolomics_de() (must contain de_tables).
#' @param config  Full pipeline config.
#' @param contrast Optional contrast name (a key of \code{de_res$de_tables}) to
#'   test. Defaults to the first contrast; an unknown name falls back to the
#'   first with a warning.
#' @return list(table, contrast, method) or NULL.
run_metabolomics_ora <- function(pre, de_res, config, contrast = NULL) {
    cfg     <- config$modes$metabolomics
    enr_cfg <- cfg$enrichment %||% list()

    # Resolve gmt_file / mapping_file like the data files: absolute paths are
    # kept as-is, relative paths are located under paths$raw.
    enr_cfg$gmt_file     <- resolve_input_path(config, enr_cfg$gmt_file)
    enr_cfg$mapping_file <- resolve_input_path(config, enr_cfg$mapping_file)

    if (!isTRUE(enr_cfg$run_enrichment)) {
        message("metabolomics ORA: disabled — skipping")
        return(NULL)
    }

    # ---- Determine significance thresholds ----
    de_cfg       <- cfg$de %||% list()
    p_cutoff     <- de_cfg$p_cutoff %||% 0.05
    if (!is.null(de_cfg$logfc_cutoff)) {
        lfc_cutoff <- de_cfg$logfc_cutoff
    } else {
        lfc_cutoff <- log2(de_cfg$linear_fc_cutoff %||% 1.5)
    }
    mapping_file <- enr_cfg$mapping_file

    # ---- Build background (all measured features mapped to IDs) ----
    mapped <- map_compounds_for_enrichment(pre$row_data, pre$expr_raw, mapping_file)
    bg_ids <- mapped$compound_names
    if (length(bg_ids) < 5) {
        message("metabolomics ORA: too few background compounds — skipping")
        return(NULL)
    }

    # ---- Build foreground (significant features) ----
    if (is.null(de_res) || is.null(de_res$de_tables) ||
        length(de_res$de_tables) == 0) {
        message("metabolomics ORA: no DE results — skipping")
        return(NULL)
    }

    contrast_name <- .resolve_enrichment_contrast(contrast, de_res, "ORA")
    de_tbl <- de_res$de_tables[[contrast_name]]

    # Significance filtering: choose p-value column via config.
    # Defaults to "P.Value" (raw) — MetaboAnalyst-style ORA convention.
    # Set enr_cfg$pvalue_column = "adj.P.Val" for FDR-controlled filtering.
    pval_col <- enr_cfg$pvalue_column %||% "P.Value"
    if (!pval_col %in% colnames(de_tbl)) {
      stop("metabolomics ORA: pvalue_column '", pval_col,
           "' not found in DE table. Available columns: ",
           paste(colnames(de_tbl), collapse = ", "))
    }
    
    sig_mask <- !is.na(de_tbl[[pval_col]]) & de_tbl[[pval_col]] < p_cutoff &
      !is.na(de_tbl$logFC) & abs(de_tbl$logFC) >= lfc_cutoff
    sig_features <- de_tbl$feature_id[sig_mask]
    
    if (length(sig_features) == 0) {
      message("metabolomics ORA: no significant features at ", pval_col, "<", p_cutoff,
              " |logFC|>=", round(lfc_cutoff, 3), " — skipping")
      return(NULL)
    }

    # Map significant feature IDs to the compound namespace via the aligned
    # per-feature map (keyed by feature id), then keep those present in bg.
    fg_ids <- unique(na.omit(mapped$feature_map[sig_features]))
    fg_ids <- fg_ids[fg_ids %in% bg_ids]

    message("ORA: ", length(fg_ids), " foreground / ", length(bg_ids),
            " background compounds (contrast: ", contrast_name, ")")

    if (length(fg_ids) < 2) {
        message("metabolomics ORA: fewer than 2 mapped foreground compounds — skipping")
        return(NULL)
    }

    # ---- Load GMT sets ----
    gene_sets <- list()
    gmt_files <- unlist(enr_cfg$gmt_file)
    if (is.null(gmt_files)) {
        message("metabolomics ORA: no GMT files — skipping")
        return(NULL)
    }
    for (gf in gmt_files) {
        if (!file.exists(gf)) next
        gmt_parsed <- read_gmt_list(gf, include_descriptions = TRUE)
        gmt <- gmt_parsed$sets
        desc_map <- gmt_parsed$descriptions
        gmt <- translate_gmt_hmdb_to_kegg(gmt, mapping_file)
        names(gmt) <- make_pathway_labels(names(gmt), desc_map)
        lib_label <- tools::file_path_sans_ext(basename(gf))
        for (nm in names(gmt)) {
            gene_sets[[nm]] <- gmt[[nm]]
            attr(gene_sets[[nm]], "library") <- lib_label
        }
        message("ORA: loaded ", length(gmt), " sets from ", basename(gf))
    }

    if (length(gene_sets) == 0) {
        message("metabolomics ORA: no gene sets — skipping")
        return(NULL)
    }

    # Filter to sets with >= 2 members in background (case-insensitive)
    bg_map <- stats::setNames(bg_ids, tolower(bg_ids))
    gene_sets_filt <- lapply(gene_sets, function(cpds) {
        m <- bg_map[tolower(cpds)]
        unname(m[!is.na(m)])
    })
    keep <- vapply(gene_sets_filt, length, integer(1)) >= 2L
    gene_sets_filt <- gene_sets_filt[keep]
    gene_sets      <- gene_sets[keep]

    if (length(gene_sets_filt) == 0) {
        message("ORA: no pathways with >= 2 matching background compounds")
        return(NULL)
    }
    message("ORA: testing ", length(gene_sets_filt), " pathways")

    # ---- Fisher's exact test per pathway ----
    N <- length(bg_ids)   # total background
    K <- length(fg_ids)   # total foreground (significant)

    ora_rows <- lapply(names(gene_sets_filt), function(pw) {
        pw_members <- gene_sets_filt[[pw]]
        n  <- length(pw_members)                          # pathway size in bg
        k  <- sum(fg_ids %in% pw_members)                 # overlap fg ∩ pathway

        # 2x2 contingency: fg-in-pw, fg-not-in-pw, bg-in-pw-not-fg, bg-not-in-pw-not-fg
        mat <- matrix(c(k, K - k, n - k, N - K - (n - k)), nrow = 2)
        mat[mat < 0] <- 0  # safety

        pval <- tryCatch(
            stats::fisher.test(mat, alternative = "greater")$p.value,
            error = function(e) NA_real_
        )

        lib <- attr(gene_sets[[pw]], "library") %||% NA_character_

        data.frame(
            pathway     = pw,
            overlap     = k,
            pathway_size = n,
            fg_size     = K,
            bg_size     = N,
            fold_enrichment = (k / K) / (n / N),
            raw_p       = pval,
            library     = lib,
            stringsAsFactors = FALSE
        )
    })

    result_df <- do.call(rbind, ora_rows)
    result_df$FDR <- stats::p.adjust(result_df$raw_p, method = "fdr")
    result_df$overlap_genes <- vapply(names(gene_sets_filt), function(pw) {
        paste(fg_ids[fg_ids %in% gene_sets_filt[[pw]]], collapse = ";")
    }, character(1))
    result_df <- result_df[order(result_df$FDR), ]

    n_sig <- sum(result_df$FDR < 0.05, na.rm = TRUE)
    message("ORA complete: ", nrow(result_df), " pathways, ",
            n_sig, " with FDR < 0.05")

    list(
        table    = result_df,
        contrast = contrast_name,
        method   = "fisher_ora"
    )
}


# ==== GSEA VIA FGSEA ==========================================================

#' Run Gene-Set Enrichment Analysis (GSEA) using fgsea
#'
#' Ranks all features by log2 fold-change from DE results and tests pathway
#' enrichment using the fgsea algorithm.
#'
#' @param pre     Preprocessing results.
#' @param de_res  DE results from run_metabolomics_de().
#' @param config  Full pipeline config.
#' @param contrast Optional contrast name (a key of \code{de_res$de_tables}) to
#'   rank. Defaults to the first contrast; an unknown name falls back to the
#'   first with a warning.
#' @return list(table, ranks, contrast, method) or NULL.
run_metabolomics_gsea <- function(pre, de_res, config, contrast = NULL) {
    cfg     <- config$modes$metabolomics
    enr_cfg <- cfg$enrichment %||% list()

    # Resolve gmt_file / mapping_file like the data files: absolute paths are
    # kept as-is, relative paths are located under paths$raw.
    enr_cfg$gmt_file     <- resolve_input_path(config, enr_cfg$gmt_file)
    enr_cfg$mapping_file <- resolve_input_path(config, enr_cfg$mapping_file)

    if (!isTRUE(enr_cfg$run_enrichment)) {
        message("metabolomics GSEA: disabled — skipping")
        return(NULL)
    }

    if (!requireNamespace("fgsea", quietly = TRUE)) {
        message("metabolomics GSEA: fgsea not available — skipping")
        return(NULL)
    }

    if (is.null(de_res) || is.null(de_res$de_tables) ||
        length(de_res$de_tables) == 0) {
        message("metabolomics GSEA: no DE results — skipping")
        return(NULL)
    }

    mapping_file <- enr_cfg$mapping_file

    # ---- Build ranked list from DE log2FC ----
    contrast_name <- .resolve_enrichment_contrast(contrast, de_res, "GSEA")
    de_tbl <- de_res$de_tables[[contrast_name]]

    # Map feature IDs to the compound namespace via the aligned per-feature map.
    mapped <- map_compounds_for_enrichment(pre$row_data, pre$expr_raw, mapping_file)
    de_tbl$compound_id <- mapped$feature_map[de_tbl$feature_id]
    de_tbl <- de_tbl[!is.na(de_tbl$compound_id) & nzchar(de_tbl$compound_id), ]
    de_tbl <- de_tbl[!is.na(de_tbl$logFC), ]

    if (nrow(de_tbl) < 5) {
        message("metabolomics GSEA: too few mapped compounds — skipping")
        return(NULL)
    }

    # Ranking metric: log2 fold-change
    ranks <- stats::setNames(de_tbl$logFC, de_tbl$compound_id)

    # Deduplicate: keep entry with largest absolute logFC
    if (any(duplicated(names(ranks)))) {
        ord <- order(abs(ranks), decreasing = TRUE)
        ranks <- ranks[ord]
        ranks <- ranks[!duplicated(names(ranks))]
    }

    ranks <- sort(ranks, decreasing = TRUE)
    message("GSEA: ranked list of ", length(ranks),
            " compounds by log2FC (contrast: ", contrast_name, ")")

    # ---- Load GMT sets ----
    gene_sets <- list()
    gmt_files <- unlist(enr_cfg$gmt_file)
    if (is.null(gmt_files)) {
        message("metabolomics GSEA: no GMT files — skipping")
        return(NULL)
    }
    for (gf in gmt_files) {
        if (!file.exists(gf)) next
        gmt_parsed <- read_gmt_list(gf, include_descriptions = TRUE)
        gmt <- gmt_parsed$sets
        desc_map <- gmt_parsed$descriptions
        gmt <- translate_gmt_hmdb_to_kegg(gmt, mapping_file)
        names(gmt) <- make_pathway_labels(names(gmt), desc_map)
        lib_label <- tools::file_path_sans_ext(basename(gf))
        for (nm in names(gmt)) {
            gene_sets[[nm]] <- gmt[[nm]]
            attr(gene_sets[[nm]], "library") <- lib_label
        }
        message("GSEA: loaded ", length(gmt), " sets from ", basename(gf))
    }

    if (length(gene_sets) == 0) {
        message("metabolomics GSEA: no gene sets — skipping")
        return(NULL)
    }

    # Filter to sets with >= 2 members present in ranked list
    # (case-insensitive match — map GMT compound names to actual ranks keys)
    avail_map <- stats::setNames(names(ranks), tolower(names(ranks)))
    gene_sets_filt <- lapply(gene_sets, function(cpds) {
        m <- avail_map[tolower(cpds)]
        unname(m[!is.na(m)])
    })
    keep <- vapply(gene_sets_filt, length, integer(1)) >= 2L
    gene_sets_filt <- gene_sets_filt[keep]
    gene_sets      <- gene_sets[keep]

    if (length(gene_sets_filt) == 0) {
        message("GSEA: no pathways with >= 2 matching compounds in ranked list")
        return(NULL)
    }
    message("GSEA: testing ", length(gene_sets_filt), " pathways")

    # ---- Run fgsea ----
    fgsea_res <- tryCatch(
        fgsea::fgsea(pathways = gene_sets_filt, stats = ranks,
                     minSize = 2, maxSize = 500, nPermSimple = 10000),
        error = function(e) {
            warning("fgsea failed: ", e$message)
            NULL
        }
    )

    if (is.null(fgsea_res) || nrow(fgsea_res) == 0) {
        message("GSEA: no results from fgsea")
        return(NULL)
    }

    result_df <- data.frame(
        pathway        = fgsea_res$pathway,
        pval           = fgsea_res$pval,
        FDR            = fgsea_res$padj,
        NES            = fgsea_res$NES,
        ES             = fgsea_res$ES,
        pathway_size   = fgsea_res$size,
        stringsAsFactors = FALSE
    )

    result_df$library <- vapply(result_df$pathway, function(pw) {
        attr(gene_sets[[pw]], "library") %||% NA_character_
    }, character(1))

    result_df$leading_edge <- vapply(seq_len(nrow(fgsea_res)), function(i) {
        le <- fgsea_res$leadingEdge[[i]]
        if (is.null(le)) return("")
        paste(le, collapse = ";")
    }, character(1))

    result_df <- result_df[order(result_df$FDR), ]

    n_sig <- sum(result_df$FDR < 0.05, na.rm = TRUE)
    message("GSEA complete: ", nrow(result_df), " pathways, ",
            n_sig, " with FDR < 0.05")

    list(
        table    = result_df,
        ranks    = ranks,
        contrast = contrast_name,
        method   = "fgsea"
    )
}


# ==== SELF-CONTAINED SET TEST VIA LIMMA ROTATION ==============================

#' Rebuild the limma DE analysis matrix + condition factor
#'
#' Enrichment methods that must run on the SAME data as DE (the self-contained
#' rotation test, QEA, and ssGSEA) use this to obtain the same expression matrix
#' (pre-scaling \code{expr_log} when a variance scaling was applied, otherwise
#' \code{expr_work}), the same biological-sample filtering, and the same
#' condition factor. It mirrors the matrix/condition setup in
#' \code{run_metabolomics_de()} (\code{03_differential.R}) without re-running DE;
#' that function remains the source of truth — if its setup changes, update this
#' to match.
#'
#' PARITY CONTRACT: this is a deliberate duplication of DE's matrix/scale/sample
#' selection, not an independent policy — \code{run_metabolomics_de()} does not
#' call this helper, so the two can drift. Any change to DE's matrix selection
#' (the expr_log-vs-expr_work branch or the biological-sample filtering) MUST be
#' mirrored here. A single canonical helper used by both is the eventual fix, but
#' is out of scope for this change.
#'
#' @param pre    Preprocessing results (expr_work/expr_log, meta, info).
#' @param config Full pipeline config.
#' @return list(mat, condition, meta): the features x biological-samples matrix,
#'   the condition factor (aligned to matrix columns), and the aligned metadata.
.metab_de_matrix_condition <- function(pre, config) {
    cfg    <- config$modes$metabolomics
    de_cfg <- cfg$de %||% list()

    condition_col <- de_cfg$condition_column %||% cfg$effects$color %||% "sample_type"
    sample_col    <- cfg$effects$samples %||% "sample_id"

    # Variance-scaling (auto/pareto/range) distorts within-group variance, so DE
    # tests on the pre-scaling matrix (expr_log). Keep this branch identical to
    # run_metabolomics_de().
    scaling_used <- pre$info$normalization$scaling %||% "none"
    mat <- if (scaling_used %in% c("auto", "pareto", "range")) {
        pre$expr_log %||% pre$expr_work
    } else {
        pre$expr_work
    }

    meta <- pre$meta[match(colnames(mat), pre$meta[[sample_col]]), , drop = FALSE]
    bio  <- filter_to_biological(mat, meta, condition_col, sample_col,
                                 label = "metabolomics self-contained",
                                 qc_flag_column = cfg$qc$qc_flag_column)
    list(mat = bio$mat, condition = bio$condition, meta = bio$meta)
}


#' Run a self-contained pathway set test via limma rotation (fry / mroast)
#'
#' Fifth enrichment method. Unlike ORA/GSEA/ssGSEA (competitive: "is this set
#' more changed than the rest?"), this is a SELF-CONTAINED test ("did this set
#' change at all?"), so it retains power when a large fraction of features move.
#' It reuses the SAME limma design and contrast as DE (via
#' \code{.metab_de_matrix_condition()} + \code{limma::makeContrasts}), the same
#' compound->KEGG mapping as the other methods
#' (\code{map_compounds_for_enrichment()}), and the same GMT pathway definitions.
#'
#' The GMT is converted to row indices with \code{limma::ids2indices()} and the
#' test is run with \code{limma::fry()} (default, deterministic — no seed needed)
#' or \code{limma::mroast()} (rotation-based, stochastic — seeded via
#' \code{withr::with_seed()}). Both directional and mixed p-values/FDR are
#' returned; the directional result ("coordinated up/down") is the interpretable
#' one. This test is more permissive than the competitive methods, so its hits
#' are hypothesis-generating, not a correction of the other methods.
#'
#' @param pre      Preprocessing results (expr_work/expr_log, meta, row_data, info).
#' @param config   Full pipeline config (reads modes$metabolomics$enrichment).
#' @param contrast_str  Contrast in "Numerator - Denominator" form (levels of the
#'   DE condition column), used to build the limma contrast — same as DE.
#' @param contrast_label Human-readable contrast label for the output table.
#' @return list(table, contrast, method) or NULL if disabled/unavailable/empty.
run_metabolomics_selfcontained <- function(pre, config, contrast_str,
                                           contrast_label = contrast_str) {
    cfg     <- config$modes$metabolomics
    enr_cfg <- cfg$enrichment %||% list()

    # Resolve gmt_file / mapping_file like the data files (see the other methods).
    enr_cfg$gmt_file     <- resolve_input_path(config, enr_cfg$gmt_file)
    enr_cfg$mapping_file <- resolve_input_path(config, enr_cfg$mapping_file)

    if (!isTRUE(enr_cfg$run_enrichment)) {
        message("metabolomics self-contained: enrichment disabled — skipping")
        return(NULL)
    }

    # Opt-in: this method runs only when explicitly enabled, since it is more
    # permissive than the competitive methods and is hypothesis-generating.
    sc_cfg <- enr_cfg$selfcontained %||% list()
    if (!isTRUE(sc_cfg$enabled)) {
        message("metabolomics self-contained: not enabled ",
                "(enrichment.selfcontained.enabled) — skipping")
        return(NULL)
    }

    if (!requireNamespace("limma", quietly = TRUE)) {
        message("metabolomics self-contained: limma not available — skipping")
        return(NULL)
    }

    # fry (deterministic) by default; mroast (rotation) is opt-in and seeded.
    sc_method <- tolower(sc_cfg$method %||% "fry")
    if (!sc_method %in% c("fry", "mroast")) {
        warning("metabolomics self-contained: unknown method '", sc_method,
                "'; falling back to 'fry'.")
        sc_method <- "fry"
    }

    mapping_file <- enr_cfg$mapping_file

    # ---- Same limma model as DE: matrix, condition, design, contrast ----
    de_in <- .metab_de_matrix_condition(pre, config)
    if (nlevels(de_in$condition) < 2) {
        message("metabolomics self-contained: fewer than 2 condition levels — skipping")
        return(NULL)
    }
    condition <- de_in$condition
    design <- stats::model.matrix(~ 0 + condition)
    colnames(design) <- levels(condition)

    contrast_matrix <- tryCatch(
        limma::makeContrasts(contrasts = contrast_str, levels = design),
        error = function(e) {
            warning("metabolomics self-contained: cannot build contrast '",
                    contrast_str, "' from condition levels (",
                    paste(levels(condition), collapse = ", "), "): ", e$message)
            NULL
        }
    )
    if (is.null(contrast_matrix)) return(NULL)
    # fry/mroast expect a numeric contrast vector of length ncol(design); a
    # single "Num - Den" string gives a one-column makeContrasts matrix.
    contrast_vec <- contrast_matrix[, 1]

    # ---- Map features to compound IDs (same universe/GMT as the other methods) ----
    mapped   <- map_compounds_for_enrichment(pre$row_data, de_in$mat, mapping_file)
    expr_mat <- as.matrix(mapped$expr_mapped)
    if (nrow(expr_mat) < 3) {
        message("metabolomics self-contained: too few mapped compounds — skipping")
        return(NULL)
    }
    # map_compounds_for_enrichment subsets rows only, so columns still align to
    # the design rows (biological samples in condition order).

    # ---- Load + translate GMT sets, tracking their source library ----
    gene_sets <- list()
    gmt_files <- unlist(enr_cfg$gmt_file)
    if (is.null(gmt_files)) {
        message("metabolomics self-contained: no GMT files — skipping")
        return(NULL)
    }
    for (gf in gmt_files) {
        if (!file.exists(gf)) next
        gmt_parsed <- read_gmt_list(gf, include_descriptions = TRUE)
        gmt <- translate_gmt_hmdb_to_kegg(gmt_parsed$sets, mapping_file)
        names(gmt) <- make_pathway_labels(names(gmt), gmt_parsed$descriptions)
        lib_label <- tools::file_path_sans_ext(basename(gf))
        for (nm in names(gmt)) {
            gene_sets[[nm]] <- gmt[[nm]]
            attr(gene_sets[[nm]], "library") <- lib_label
        }
        message("self-contained: loaded ", length(gmt), " sets from ", basename(gf))
    }
    if (length(gene_sets) == 0) {
        message("metabolomics self-contained: no gene sets — skipping")
        return(NULL)
    }

    # ---- Convert GMT to row indices and keep sets with >= 2 members present ----
    index <- limma::ids2indices(gene_sets, identifiers = rownames(expr_mat))
    index <- index[vapply(index, length, integer(1)) >= 2L]
    if (length(index) == 0) {
        message("self-contained: no pathways with >= 2 matching compounds")
        return(NULL)
    }
    message("self-contained: testing ", length(index), " pathways (", sc_method, ")")

    # ---- Run the rotation test ----
    res_raw <- tryCatch({
        if (sc_method == "mroast") {
            # Rotation-based → stochastic; seed for reproducibility (repo convention).
            n_rot <- sc_cfg$n_rotations %||% 9999
            seed  <- sc_cfg$seed %||% 42
            withr::with_seed(seed, limma::mroast(
                y = expr_mat, index = index, design = design,
                contrast = contrast_vec, nrot = n_rot
            ))
        } else {
            # fry approximates an infinite-rotation mroast → deterministic.
            limma::fry(y = expr_mat, index = index, design = design,
                       contrast = contrast_vec)
        }
    }, error = function(e) {
        warning("metabolomics self-contained (", sc_method, ") failed: ", e$message)
        NULL
    })
    if (is.null(res_raw) || nrow(res_raw) == 0) {
        message("self-contained: no results")
        return(NULL)
    }

    # ---- Tidy to the shared enrichment schema ----
    # fry/mroast both return: NGenes, Direction, PValue, FDR, PValue.Mixed,
    # FDR.Mixed (rownames = pathway). The directional PValue/FDR is the
    # interpretable "coordinated up/down" result; the mixed one is retained too.
    lib_lookup <- vapply(names(gene_sets), function(nm) {
        attr(gene_sets[[nm]], "library") %||% NA_character_
    }, character(1))

    result_df <- data.frame(
        pathway      = rownames(res_raw),
        n_hits       = res_raw$NGenes,
        direction    = as.character(res_raw$Direction),
        PValue       = res_raw$PValue,
        FDR          = res_raw$FDR,
        PValue_mixed = res_raw[["PValue.Mixed"]],
        FDR_mixed    = res_raw[["FDR.Mixed"]],
        library      = unname(lib_lookup[rownames(res_raw)]),
        stringsAsFactors = FALSE
    )
    result_df <- result_df[order(result_df$FDR), ]

    n_sig <- sum(result_df$FDR < 0.05, na.rm = TRUE)
    message("self-contained complete: ", nrow(result_df), " pathways, ",
            n_sig, " with directional FDR < 0.05 (contrast: ", contrast_label, ")")

    list(
        table    = result_df,
        contrast = contrast_label,
        method   = paste0("limma_", sc_method)
    )
}


# ==== ENRICHMENT PLOTS ========================================================

#' Enrichment barplot (for QEA or ssGSEA results)
#'
#' Plots the top pathways as horizontal bars of -log10(FDR). When a "library"
#' column is present, bars are colored by database source; otherwise a
#' continuous gradient is used. A dashed line marks FDR = 0.05.
#'
#' @param enrich_df data.frame with pathway, FDR columns (and optional library).
#' @param top_n     Number of top pathways to show.
#' @param title     Plot title.
#' @return ggplot object or NULL if input is empty.
plot_enrichment_barplot <- function(enrich_df, top_n = 20,
                                    title = "Pathway Enrichment") {
    if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)

    # Select the top pathways and compute -log10(FDR) for the y-axis
    top_df <- utils::head(enrich_df, top_n)
    top_df$neg_log10_fdr <- -log10(pmax(top_df$FDR, 1e-20))

    # Truncate long pathway names and order for a horizontal bar layout
    top_df$pathway_short <- ifelse(
        nchar(top_df$pathway) > 50,
        paste0(substr(top_df$pathway, 1, 47), "..."),
        top_df$pathway
    )
    top_df$pathway_short <- factor(top_df$pathway_short,
                                    levels = rev(top_df$pathway_short))

    # Color by database source when a library column exists, else by FDR gradient
    has_library <- "library" %in% colnames(top_df)

    if (has_library) {
        top_df$lib_label <- toupper(gsub("_pathway$", "", top_df$library))
        p <- ggplot2::ggplot(top_df, ggplot2::aes(x = pathway_short,
                                                    y = neg_log10_fdr,
                                                    fill = lib_label)) +
            ggplot2::geom_col() +
            ggplot2::labs(fill = "Database")
    } else {
        p <- ggplot2::ggplot(top_df, ggplot2::aes(x = pathway_short,
                                                    y = neg_log10_fdr,
                                                    fill = neg_log10_fdr)) +
            ggplot2::geom_col() +
            ggplot2::scale_fill_gradient(low = "steelblue", high = "firebrick",
                                          name = "-log10(FDR)")
    }

    # Add FDR = 0.05 threshold line and style the plot
    p +
        ggplot2::coord_flip() +
        ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                             color = "grey40") +
        ggplot2::labs(title = title,
                       subtitle = paste0("Top ", nrow(top_df), " pathways"),
                       x = NULL, y = "-log10(FDR)") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 9),
            plot.title  = ggplot2::element_text(face = "bold")
        )
}


#' ssGSEA pathway boxplots for significant pathways
#'
#' Creates faceted boxplots comparing ssGSEA enrichment scores between
#' conditions for the top significant pathways (FDR < 0.05). Individual
#' sample points are overlaid with jitter.
#'
#' @param scores       Pathway scores matrix (pathways x samples).
#' @param conditions   Factor of sample conditions (same order as score columns).
#' @param results_df   ssGSEA results data.frame with a "significant" logical flag.
#' @param max_pathways Maximum number of significant pathways to display.
#' @return ggplot object or NULL if no pathways are significant.
plot_ssgsea_boxplots <- function(scores, conditions, results_df,
                                  max_pathways = 12) {
    # Select the top significant pathways up to max_pathways
    sig_pathways <- results_df$pathway[results_df$significant == TRUE]
    if (length(sig_pathways) == 0) return(NULL)

    sig_pathways <- utils::head(sig_pathways, max_pathways)
    sig_scores   <- scores[sig_pathways, , drop = FALSE]

    # Reshape score matrix into long format for ggplot2
    plot_data <- do.call(rbind, lapply(sig_pathways, function(pw) {
        data.frame(
            pathway   = pw,
            score     = sig_scores[pw, ],
            condition = as.character(conditions),
            stringsAsFactors = FALSE
        )
    }))

    # Truncate long pathway names for facet labels
    plot_data$pathway_short <- ifelse(
        nchar(plot_data$pathway) > 40,
        paste0(substr(plot_data$pathway, 1, 37), "..."),
        plot_data$pathway
    )

    # Draw boxplots with jittered sample points, one facet per pathway
    ggplot2::ggplot(plot_data, ggplot2::aes(x = condition, y = score,
                                             fill = condition)) +
        ggplot2::geom_boxplot(outlier.shape = 16, outlier.size = 1.5) +
        ggplot2::geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
        ggplot2::facet_wrap(~pathway_short, scales = "free_y") +
        ggplot2::labs(
            title = "ssGSEA Pathway Scores by Condition",
            subtitle = paste0("Top ", length(sig_pathways),
                              " significant pathways (FDR < 0.05)"),
            x = NULL, y = "ssGSEA Enrichment Score", fill = "Condition"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            strip.text = ggplot2::element_text(size = 8),
            plot.title = ggplot2::element_text(face = "bold")
        )
}


#' GSEA NES barplot
#'
#' Horizontal barplot of Normalized Enrichment Scores (NES) for top pathways,
#' colored by direction (up/down) with FDR < 0.05 threshold marking.
#'
#' @param gsea_df  data.frame from run_metabolomics_gsea() with pathway, NES, FDR.
#' @param top_n    Number of top pathways to display (by FDR).
#' @param title    Plot title.
#' @return ggplot object or NULL.
plot_gsea_nes_barplot <- function(gsea_df, top_n = 20,
                                   title = "GSEA — Normalized Enrichment Scores") {
    if (is.null(gsea_df) || nrow(gsea_df) == 0) return(NULL)

    top_df <- utils::head(gsea_df[order(gsea_df$FDR), ], top_n)
    top_df$pathway_short <- ifelse(
        nchar(top_df$pathway) > 50,
        paste0(substr(top_df$pathway, 1, 47), "..."),
        top_df$pathway
    )
    top_df$pathway_short <- factor(top_df$pathway_short,
                                    levels = rev(top_df$pathway_short))
    top_df$direction <- ifelse(top_df$NES > 0, "Up", "Down")

    ggplot2::ggplot(top_df, ggplot2::aes(x = pathway_short, y = NES,
                                          fill = direction)) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = c(Up = "firebrick", Down = "steelblue"),
                                    name = "Direction") +
        ggplot2::geom_hline(yintercept = 0, color = "grey30") +
        ggplot2::labs(title = title,
                       subtitle = paste0("Top ", nrow(top_df),
                                         " pathways by FDR"),
                       x = NULL, y = "NES") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 9),
            plot.title  = ggplot2::element_text(face = "bold")
        )
}


#' ORA dot plot
#'
#' Dot plot showing fold-enrichment vs pathway, sized by overlap count and
#' colored by −log10(FDR).
#'
#' @param ora_df  data.frame from run_metabolomics_ora().
#' @param top_n   Number of top pathways to display (by FDR).
#' @param title   Plot title.
#' @return ggplot object or NULL.
plot_ora_dotplot <- function(ora_df, top_n = 20,
                              title = "ORA — Pathway Over-Representation") {
    if (is.null(ora_df) || nrow(ora_df) == 0) return(NULL)

    top_df <- utils::head(ora_df[order(ora_df$FDR), ], top_n)
    top_df$neg_log10_fdr <- -log10(pmax(top_df$FDR, 1e-20))
    top_df$pathway_short <- ifelse(
        nchar(top_df$pathway) > 50,
        paste0(substr(top_df$pathway, 1, 47), "..."),
        top_df$pathway
    )
    top_df$pathway_short <- factor(top_df$pathway_short,
                                    levels = rev(top_df$pathway_short))

    ggplot2::ggplot(top_df, ggplot2::aes(x = fold_enrichment,
                                          y = pathway_short,
                                          size = overlap,
                                          color = neg_log10_fdr)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_gradient(low = "steelblue", high = "firebrick",
                                       name = "-log10(FDR)") +
        ggplot2::scale_size_continuous(name = "Overlap", range = c(2, 8)) +
        ggplot2::geom_vline(xintercept = 1, linetype = "dashed",
                             color = "grey40") +
        ggplot2::labs(title = title,
                       subtitle = paste0("Top ", nrow(top_df), " pathways by FDR"),
                       x = "Fold Enrichment", y = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 9),
            plot.title  = ggplot2::element_text(face = "bold")
        )
}


#' ORA lollipop plot
#'
#' Lollipop plot of −log10(FDR) per pathway, colored by fold-enrichment.
#'
#' @param ora_df  data.frame from run_metabolomics_ora().
#' @param top_n   Number of top pathways to display.
#' @param title   Plot title.
#' @return ggplot object or NULL.
plot_ora_lollipop <- function(ora_df, top_n = 20,
                               title = "ORA — Pathway Enrichment") {
    if (is.null(ora_df) || nrow(ora_df) == 0) return(NULL)

    top_df <- utils::head(ora_df[order(ora_df$FDR), ], top_n)
    top_df$neg_log10_fdr <- -log10(pmax(top_df$FDR, 1e-20))
    top_df$pathway_short <- ifelse(
        nchar(top_df$pathway) > 50,
        paste0(substr(top_df$pathway, 1, 47), "..."),
        top_df$pathway
    )
    top_df$pathway_short <- factor(top_df$pathway_short,
                                    levels = rev(top_df$pathway_short))

    ggplot2::ggplot(top_df, ggplot2::aes(x = neg_log10_fdr,
                                          y = pathway_short)) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = neg_log10_fdr,
                                            yend = pathway_short),
                               color = "grey60") +
        ggplot2::geom_point(ggplot2::aes(size = overlap,
                                          color = fold_enrichment)) +
        ggplot2::scale_color_gradient(low = "steelblue", high = "firebrick",
                                       name = "Fold Enrichment") +
        ggplot2::scale_size_continuous(name = "Overlap", range = c(2, 7)) +
        ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed",
                             color = "grey40") +
        ggplot2::labs(title = title,
                       subtitle = paste0("Top ", nrow(top_df), " pathways by FDR"),
                       x = "-log10(FDR)", y = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 9),
            plot.title  = ggplot2::element_text(face = "bold")
        )
}


#' GSEA lollipop plot
#'
#' Lollipop plot of NES per pathway, colored by direction and sized by
#' −log10(FDR).
#'
#' @param gsea_df  data.frame from run_metabolomics_gsea().
#' @param top_n    Number of top pathways to display.
#' @param title    Plot title.
#' @return ggplot object or NULL.
plot_gsea_lollipop <- function(gsea_df, top_n = 20,
                                title = "GSEA — Pathway Enrichment") {
    if (is.null(gsea_df) || nrow(gsea_df) == 0) return(NULL)

    top_df <- utils::head(gsea_df[order(gsea_df$FDR), ], top_n)
    top_df$neg_log10_fdr <- -log10(pmax(top_df$FDR, 1e-20))
    top_df$pathway_short <- ifelse(
        nchar(top_df$pathway) > 50,
        paste0(substr(top_df$pathway, 1, 47), "..."),
        top_df$pathway
    )
    top_df$pathway_short <- factor(top_df$pathway_short,
                                    levels = rev(top_df$pathway_short))
    top_df$direction <- ifelse(top_df$NES > 0, "Up", "Down")

    ggplot2::ggplot(top_df, ggplot2::aes(x = NES, y = pathway_short)) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES,
                                            yend = pathway_short),
                               color = "grey60") +
        ggplot2::geom_point(ggplot2::aes(size = neg_log10_fdr,
                                          color = direction)) +
        ggplot2::scale_color_manual(values = c(Up = "firebrick",
                                                Down = "steelblue"),
                                     name = "Direction") +
        ggplot2::scale_size_continuous(name = "-log10(FDR)", range = c(2, 7)) +
        ggplot2::geom_vline(xintercept = 0, color = "grey30") +
        ggplot2::labs(title = title,
                       subtitle = paste0("Top ", nrow(top_df), " pathways by FDR"),
                       x = "NES", y = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 9),
            plot.title  = ggplot2::element_text(face = "bold")
        )
}


#' QEA lollipop plot
#'
#' Lollipop plot of −log10(FDR) per pathway, sized by number of hits.
#'
#' @param qea_df  data.frame from run_metabolomics_qea() with pathway, FDR, hits.
#' @param top_n   Number of top pathways to display.
#' @param title   Plot title.
#' @return ggplot object or NULL.
plot_qea_lollipop <- function(qea_df, top_n = 20,
                               title = "QEA — Pathway Enrichment") {
    if (is.null(qea_df) || nrow(qea_df) == 0) return(NULL)

    top_df <- utils::head(qea_df[order(qea_df$FDR), ], top_n)
    top_df$neg_log10_fdr <- -log10(pmax(top_df$FDR, 1e-20))
    top_df$pathway_short <- ifelse(
        nchar(top_df$pathway) > 50,
        paste0(substr(top_df$pathway, 1, 47), "..."),
        top_df$pathway
    )
    top_df$pathway_short <- factor(top_df$pathway_short,
                                    levels = rev(top_df$pathway_short))

    has_library <- "library" %in% colnames(top_df)

    p <- ggplot2::ggplot(top_df, ggplot2::aes(x = neg_log10_fdr,
                                                y = pathway_short)) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = neg_log10_fdr,
                                            yend = pathway_short),
                               color = "grey60")

    if (has_library) {
        top_df$lib_label <- toupper(gsub("_pathway$", "", top_df$library))
        p <- p +
            ggplot2::geom_point(data = top_df,
                                 ggplot2::aes(size = hits, color = lib_label)) +
            ggplot2::labs(color = "Database")
    } else {
        p <- p +
            ggplot2::geom_point(ggplot2::aes(size = hits),
                                 color = "firebrick")
    }

    p +
        ggplot2::scale_size_continuous(name = "Hits", range = c(2, 7)) +
        ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed",
                             color = "grey40") +
        ggplot2::labs(title = title,
                       subtitle = paste0("Top ", nrow(top_df), " pathways by FDR"),
                       x = "-log10(FDR)", y = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 9),
            plot.title  = ggplot2::element_text(face = "bold")
        )
}


#' ssGSEA lollipop plot
#'
#' Lollipop plot of −log10(FDR) per pathway, colored by enrichment direction
#' (based on score_diff) and sized by absolute score difference.
#'
#' @param ssgsea_df data.frame from run_metabolomics_ssgsea().
#' @param top_n     Number of top pathways to display.
#' @param title     Plot title.
#' @return ggplot object or NULL.
plot_ssgsea_lollipop <- function(ssgsea_df, top_n = 20,
                                  title = "ssGSEA — Pathway Enrichment") {
    if (is.null(ssgsea_df) || nrow(ssgsea_df) == 0) return(NULL)
    if (!"FDR" %in% colnames(ssgsea_df)) return(NULL)

    top_df <- utils::head(ssgsea_df[order(ssgsea_df$FDR), ], top_n)
    top_df$neg_log10_fdr <- -log10(pmax(top_df$FDR, 1e-20))
    top_df$pathway_short <- ifelse(
        nchar(top_df$pathway) > 50,
        paste0(substr(top_df$pathway, 1, 47), "..."),
        top_df$pathway
    )
    top_df$pathway_short <- factor(top_df$pathway_short,
                                    levels = rev(top_df$pathway_short))

    if ("score_diff" %in% colnames(top_df)) {
        top_df$direction <- ifelse(top_df$score_diff > 0, "Up", "Down")
        top_df$abs_diff  <- abs(top_df$score_diff)
    } else {
        top_df$direction <- "N/A"
        top_df$abs_diff  <- 1
    }

    ggplot2::ggplot(top_df, ggplot2::aes(x = neg_log10_fdr,
                                          y = pathway_short)) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = neg_log10_fdr,
                                            yend = pathway_short),
                               color = "grey60") +
        ggplot2::geom_point(ggplot2::aes(size = abs_diff,
                                          color = direction)) +
        ggplot2::scale_color_manual(values = c(Up = "firebrick",
                                                Down = "steelblue",
                                                "N/A" = "grey50"),
                                     name = "Direction") +
        ggplot2::scale_size_continuous(name = "|Score Diff|", range = c(2, 7)) +
        ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed",
                             color = "grey40") +
        ggplot2::labs(title = title,
                       subtitle = paste0("Top ", nrow(top_df), " pathways by FDR"),
                       x = "-log10(FDR)", y = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 9),
            plot.title  = ggplot2::element_text(face = "bold")
        )
}


#' Self-contained set-test lollipop plot
#'
#' Lollipop plot of −log10(directional FDR) per pathway from the limma rotation
#' test, colored by coordinated direction (Up/Down) and sized by the number of
#' matched compounds (n_hits).
#'
#' @param sc_df  data.frame from run_metabolomics_selfcontained() with pathway,
#'   FDR, direction, n_hits.
#' @param top_n  Number of top pathways to display (by FDR).
#' @param title  Plot title.
#' @return ggplot object or NULL.
plot_selfcontained_lollipop <- function(sc_df, top_n = 20,
                                        title = "Self-contained set test") {
    if (is.null(sc_df) || nrow(sc_df) == 0) return(NULL)
    if (!"FDR" %in% colnames(sc_df)) return(NULL)

    top_df <- utils::head(sc_df[order(sc_df$FDR), ], top_n)
    top_df$neg_log10_fdr <- -log10(pmax(top_df$FDR, 1e-20))
    top_df$pathway_short <- ifelse(
        nchar(top_df$pathway) > 50,
        paste0(substr(top_df$pathway, 1, 47), "..."),
        top_df$pathway
    )
    top_df$pathway_short <- factor(top_df$pathway_short,
                                    levels = rev(top_df$pathway_short))

    direction <- if ("direction" %in% colnames(top_df)) {
        as.character(top_df$direction)
    } else {
        rep("N/A", nrow(top_df))
    }
    top_df$direction <- ifelse(direction %in% c("Up", "Down"), direction, "N/A")
    top_df$n_hits <- if ("n_hits" %in% colnames(top_df)) top_df$n_hits else 1

    ggplot2::ggplot(top_df, ggplot2::aes(x = neg_log10_fdr, y = pathway_short)) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = neg_log10_fdr,
                                            yend = pathway_short),
                               color = "grey60") +
        ggplot2::geom_point(ggplot2::aes(size = n_hits, color = direction)) +
        ggplot2::scale_color_manual(values = c(Up = "firebrick",
                                                Down = "steelblue",
                                                "N/A" = "grey50"),
                                     name = "Direction") +
        ggplot2::scale_size_continuous(name = "Compounds", range = c(2, 7)) +
        ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed",
                             color = "grey40") +
        ggplot2::labs(title = title,
                       subtitle = paste0("Top ", nrow(top_df),
                                         " pathways by directional FDR"),
                       x = "-log10(directional FDR)", y = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 9),
            plot.title  = ggplot2::element_text(face = "bold")
        )
}
