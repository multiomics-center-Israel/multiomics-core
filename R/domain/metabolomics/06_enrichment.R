# R/domain/metabolomics/06_enrichment.R
#
# Pathway enrichment analysis for metabolomics:
#   1. QEA (Quantitative Enrichment Analysis) via globaltest
#   2. ssGSEA via GSVA + Wilcoxon rank-sum test
#
# Operates on the standard pre-processing contract (expr_work, meta, row_data).
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

#' Map feature IDs to compound identifiers for enrichment
#'
#' Extracts compound names from row_data and optionally maps HMDB to KEGG.
#'
#' @param row_data     data.frame with feature annotations.
#' @param expr_mat     Expression matrix (features x samples).
#' @param mapping_file Optional path to HMDB-to-KEGG mapping TSV.
#' @return list(expr_mapped, compound_names) with remapped rownames.
map_compounds_for_enrichment <- function(row_data, expr_mat, mapping_file = NULL) {
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
    if (!is.null(mapping_file) && file.exists(mapping_file)) {
        mapping_df <- readr::read_delim(mapping_file, delim = "\t",
                                        col_types = readr::cols(),
                                        show_col_types = FALSE)

        # Locate HMDB IDs — try dedicated column first, then parse feature_id
        hmdb_ids <- NULL
        if ("HMDB" %in% colnames(row_data)) {
            hmdb_ids <- as.character(row_data$HMDB)
        } else if ("feature_id" %in% colnames(row_data)) {
            # Try extracting HMDB from feature_id (HMDB|Name format)
            hmdb_ids <- sub("\\|.*$", "", as.character(row_data$feature_id))
        }

        # Vectorised lookup: build named vector HMDB->KEGG then index
        if (!is.null(hmdb_ids) &&
            "HMDB" %in% colnames(mapping_df) &&
            "KEGG" %in% colnames(mapping_df)) {
            map_vec <- stats::setNames(mapping_df$KEGG, mapping_df$HMDB)
            mapped_kegg <- map_vec[hmdb_ids]
            has_map <- !is.na(mapped_kegg) & mapped_kegg != ""
            compound_names[has_map] <- mapped_kegg[has_map]
            message("enrichment: mapped ", sum(has_map), " features to KEGG IDs")
        }
    }

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
        compound_names = names_use[!dup]
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
#' @param pre     Preprocessing results (must include expr_raw, meta, row_data).
#' @param config  Full pipeline config (reads modes$metabolomics$enrichment).
#' @return list(table, per_library_tables, method) or NULL if disabled/unavailable.
run_metabolomics_qea <- function(pre, config) {
    cfg      <- config$modes$metabolomics
    enr_cfg  <- cfg$enrichment %||% list()

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

    # ---- Prepare QEA data ----
    # Map features to compound IDs and transpose into the samples-by-compounds
    # layout expected by globaltest. The result is written to a temp TSV.
    mapped <- map_compounds_for_enrichment(pre$row_data, pre$expr_raw,
                                            mapping_file)
    if (nrow(mapped$expr_mapped) < 3) {
        message("metabolomics QEA: too few compounds — skipping")
        return(NULL)
    }

    # Build globaltest input format: rows=samples, cols=compounds.
    # Prepend Sample and Group columns required by run_qea_gmt_internal().
    meta <- pre$meta
    mat_t <- as.data.frame(t(mapped$expr_mapped), check.names = FALSE)
    conditions <- meta[[condition_col]][match(rownames(mat_t), meta[[sample_col]])]
    df_t <- cbind(
        data.frame(Sample = rownames(mat_t), Group = as.character(conditions),
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
            run_qea_gmt_internal(data_file, gmt_lookup[[lib]], mapping_file),
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

#' Run QEA via globaltest on a single GMT file
#'
#' Parses the GMT, optionally translates HMDB IDs to KEGG, filters to pathways
#' with at least 2 matching compounds in the data, and runs globaltest::gt().
#' Returns a data.frame of pathway-level p-values and hit counts.
#'
#' @param data_file    Path to the tab-delimited QEA data (Sample, Group, compounds).
#' @param gmt_file     Path to GMT file defining pathway sets.
#' @param mapping_file Path to HMDB-to-KEGG mapping TSV (or NULL).
#' @return data.frame(pathway, raw_p, hits) or NULL if no testable pathways.
run_qea_gmt_internal <- function(data_file, gmt_file, mapping_file) {
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

    # Keep only pathways where at least 2 members are present in the data
    available <- colnames(X)
    subsets <- lapply(gmt_list, function(cpds) cpds[cpds %in% available])
    keep <- vapply(subsets, length, integer(1)) >= 2L
    subsets <- subsets[keep]

    if (length(subsets) == 0) {
        message("QEA GMT: no pathways with >= 2 matching compounds")
        return(NULL)
    }

    message("QEA GMT: ", length(subsets), " / ", length(gmt_list),
            " pathways with >= 2 matching compounds")

    # Run the global test: tests whether compound profiles in each pathway
    # are collectively associated with the condition factor.
    res_gt <- tryCatch(
        globaltest::gt(response, X, subsets = subsets),
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
#' @param pre     Preprocessing results (must include expr_raw, meta, row_data).
#' @param config  Full pipeline config (reads modes$metabolomics$enrichment).
#' @return list(table, scores, method) or NULL if disabled/unavailable.
run_metabolomics_ssgsea <- function(pre, config) {
    cfg     <- config$modes$metabolomics
    enr_cfg <- cfg$enrichment %||% list()

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
    # Map features to enrichment-ready IDs (Name/HMDB/KEGG)
    mapped <- map_compounds_for_enrichment(pre$row_data, pre$expr_raw,
                                            mapping_file)
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
    available <- rownames(expr_mat)
    gene_sets <- lapply(gene_sets, function(cpds) cpds[cpds %in% available])
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
