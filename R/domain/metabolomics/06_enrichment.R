# R/domain/metabolomics/06_enrichment.R
#
# Pathway enrichment analysis for metabolomics:
#   1. QEA (Quantitative Enrichment Analysis) via MetaboAnalystR / globaltest
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
    if (!file.exists(gmt_file)) {
        if (include_descriptions) return(list(sets = list(), descriptions = character(0)))
        return(list())
    }

    lines    <- readLines(gmt_file, warn = FALSE)
    gmt_list <- list()
    desc_map <- character(0)

    for (line in lines) {
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) < 3) next

        pathway_id   <- parts[1]
        pathway_desc <- parts[2]
        members      <- parts[3:length(parts)]
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

    # HMDB-to-KEGG mapping (if provided)
    if (!is.null(mapping_file) && file.exists(mapping_file)) {
        mapping_df <- readr::read_delim(mapping_file, delim = "\t",
                                        col_types = readr::cols(),
                                        show_col_types = FALSE)

        hmdb_ids <- NULL
        if ("HMDB" %in% colnames(row_data)) {
            hmdb_ids <- as.character(row_data$HMDB)
        } else if ("feature_id" %in% colnames(row_data)) {
            # Try extracting HMDB from feature_id (HMDB|Name format)
            hmdb_ids <- sub("\\|.*$", "", as.character(row_data$feature_id))
        }

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
translate_gmt_hmdb_to_kegg <- function(gmt_list, mapping_file) {
    if (is.null(mapping_file) || !file.exists(mapping_file)) return(gmt_list)

    mapping_df <- readr::read_delim(mapping_file, delim = "\t",
                                    col_types = readr::cols(),
                                    show_col_types = FALSE)

    if (!("HMDB" %in% colnames(mapping_df) && "KEGG" %in% colnames(mapping_df))) {
        return(gmt_list)
    }

    hmdb2kegg <- stats::setNames(mapping_df$KEGG, mapping_df$HMDB)
    lapply(gmt_list, function(cpds) {
        mapped <- hmdb2kegg[cpds]
        ifelse(!is.na(mapped) & mapped != "", mapped, cpds)
    })
}


# ==== QEA VIA METABOANALYSTR =================================================

#' Run QEA enrichment using MetaboAnalystR + globalTest
#'
#' @param pre     Preprocessing results.
#' @param config  Full pipeline config.
#' @return list(table, per_library_tables, method) or NULL if disabled/unavailable.
run_metabolomics_qea <- function(pre, config) {
    cfg      <- config$modes$metabolomics
    enr_cfg  <- cfg$enrichment %||% list()

    if (!isTRUE(enr_cfg$run_enrichment)) {
        message("metabolomics QEA: disabled in config — skipping")
        return(NULL)
    }

    if (!requireNamespace("MetaboAnalystR", quietly = TRUE)) {
        message("metabolomics QEA: MetaboAnalystR not available — skipping")
        return(NULL)
    }

    condition_col <- enr_cfg$condition_column %||%
                     cfg$de$condition_column %||%
                     cfg$effects$color %||% "sample_type"
    sample_col <- cfg$effects$samples %||% "sample_id"

    # ---- Determine libraries to run ----
    libraries <- enr_cfg$libraries
    if (is.null(libraries)) {
        lib_single <- enr_cfg$library
        libraries <- if (!is.null(lib_single)) list(lib_single) else list("smpdb_pathway")
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
    mapped <- map_compounds_for_enrichment(pre$row_data, pre$expr_raw,
                                            mapping_file)
    if (nrow(mapped$expr_mapped) < 3) {
        message("metabolomics QEA: too few compounds — skipping")
        return(NULL)
    }

    # Build MetaboAnalystR input format: rows=samples, cols=compounds
    meta <- pre$meta
    mat_t <- as.data.frame(t(mapped$expr_mapped), check.names = FALSE)
    conditions <- meta[[condition_col]][match(rownames(mat_t), meta[[sample_col]])]
    df_t <- cbind(
        data.frame(Sample = rownames(mat_t), Group = as.character(conditions),
                   stringsAsFactors = FALSE),
        mat_t
    )

    tmp_dir  <- tempfile("metabo_qea_")
    dir.create(tmp_dir, recursive = TRUE)
    data_file <- file.path(tmp_dir, "combined_data.txt")
    utils::write.table(df_t, data_file, sep = "\t", row.names = FALSE,
                        quote = FALSE)

    message("QEA data: ", nrow(df_t), " samples x ", ncol(df_t) - 2, " compounds")

    # ---- Run QEA for each library ----
    all_results <- list()

    for (lib in libraries) {
        if (lib %in% names(gmt_lookup)) {
            result <- tryCatch(
                run_qea_gmt_internal(data_file, gmt_lookup[[lib]], mapping_file),
                error = function(e) {
                    warning("QEA GMT failed for ", lib, ": ", e$message)
                    NULL
                }
            )
        } else {
            result <- tryCatch(
                run_qea_metaboanalyst_internal(data_file, lib, enr_cfg),
                error = function(e) {
                    warning("QEA failed for ", lib, ": ", e$message)
                    NULL
                }
            )
        }
        if (!is.null(result)) all_results[[lib]] <- result
    }

    unlink(tmp_dir, recursive = TRUE)

    if (length(all_results) == 0) {
        message("metabolomics QEA: no results from any library")
        return(NULL)
    }

    # Combine results
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

#' Run QEA via MetaboAnalystR for a single built-in library
run_qea_metaboanalyst_internal <- function(data_file, library_name, enr_cfg) {
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)

    row_norm   <- enr_cfg$metaboanalyst_row_norm   %||% "NULL"
    trans_norm <- enr_cfg$metaboanalyst_trans_norm  %||% "LogNorm"
    scale_norm <- enr_cfg$metaboanalyst_scale_norm  %||% "MeanCenter"

    mSet <- MetaboAnalystR::InitDataObjects("conc", "msetqea", FALSE)
    mSet <- MetaboAnalystR::Read.TextData(mSet, data_file, "rowu", "disc")
    mSet <- MetaboAnalystR::SanityCheckData(mSet)
    mSet <- MetaboAnalystR::ReplaceMin(mSet)
    mSet <- MetaboAnalystR::PreparePrenormData(mSet)
    mSet <- MetaboAnalystR::Normalization(mSet, row_norm, trans_norm, scale_norm,
                                           "S10T0", ratio = FALSE, ratioNum = 20)

    mSet <- MetaboAnalystR::CrossReferencing(mSet, "name")
    mSet <- MetaboAnalystR::CreateMappingResultTable(mSet)
    mSet <- MetaboAnalystR::SetMetabolomeFilter(mSet, FALSE)
    mSet <- MetaboAnalystR::SetCurrentMsetLib(mSet, library_name, 2)
    mSet <- MetaboAnalystR::CalculateGlobalTestScore(mSet)

    qea_mat <- mSet$analSet$qea.mat
    if (is.null(qea_mat) || nrow(qea_mat) == 0) return(NULL)

    data.frame(
        pathway = rownames(qea_mat),
        raw_p   = as.numeric(qea_mat[, "Raw p"]),
        hits    = as.numeric(qea_mat[, "Hits"]),
        stringsAsFactors = FALSE
    )
}


#' Run QEA via globaltest on a custom GMT file
run_qea_gmt_internal <- function(data_file, gmt_file, mapping_file) {
    if (!requireNamespace("globaltest", quietly = TRUE)) {
        stop("Package 'globaltest' required for GMT enrichment.")
    }

    gmt_parsed <- read_gmt_list(gmt_file, include_descriptions = TRUE)
    gmt_list   <- gmt_parsed$sets
    desc_map   <- gmt_parsed$descriptions

    if (length(gmt_list) == 0) return(NULL)

    # Translate HMDB to KEGG in GMT if mapping available
    gmt_list <- translate_gmt_hmdb_to_kegg(gmt_list, mapping_file)

    # Label pathways with descriptions
    names(gmt_list) <- make_pathway_labels(names(gmt_list), desc_map)

    # Read data
    df <- utils::read.table(data_file, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, check.names = FALSE)
    response <- factor(df$Group)
    X <- as.matrix(df[, -c(1, 2), drop = FALSE])

    # Filter to pathways with >= 2 matching compounds
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

    res_gt <- tryCatch(
        globaltest::gt(response, X, subsets = subsets),
        error = function(e) { warning("globaltest error: ", e$message); NULL }
    )
    if (is.null(res_gt)) return(NULL)

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
#' @param pre     Preprocessing results.
#' @param config  Full pipeline config.
#' @return list(table, scores, method) or NULL if disabled/unavailable.
run_metabolomics_ssgsea <- function(pre, config) {
    cfg     <- config$modes$metabolomics
    enr_cfg <- cfg$enrichment %||% list()

    if (!isTRUE(enr_cfg$run_enrichment)) {
        message("metabolomics ssGSEA: disabled — skipping")
        return(NULL)
    }

    if (!requireNamespace("GSVA", quietly = TRUE)) {
        message("metabolomics ssGSEA: GSVA not available — skipping")
        return(NULL)
    }

    condition_col <- enr_cfg$condition_column %||%
                     cfg$de$condition_column %||%
                     cfg$effects$color %||% "sample_type"
    sample_col <- cfg$effects$samples %||% "sample_id"
    mapping_file <- enr_cfg$mapping_file

    # ---- Build expression matrix with compound IDs ----
    mapped <- map_compounds_for_enrichment(pre$row_data, pre$expr_raw,
                                            mapping_file)
    expr_mat <- as.matrix(mapped$expr_mapped)
    if (nrow(expr_mat) < 2) {
        message("metabolomics ssGSEA: too few features — skipping")
        return(NULL)
    }

    # ---- Load gene sets from GMT files ----
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

    # Filter to sets with >= 2 members in data
    available <- rownames(expr_mat)
    gene_sets <- lapply(gene_sets, function(cpds) cpds[cpds %in% available])
    gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) >= 2L]

    if (length(gene_sets) == 0) {
        message("ssGSEA: no gene sets with >= 2 matching compounds")
        return(NULL)
    }
    message("ssGSEA: ", length(gene_sets), " gene sets with >= 2 members")

    # ---- Run ssGSEA ----
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
    meta <- pre$meta
    conditions <- factor(meta[[condition_col]][
        match(colnames(scores), meta[[sample_col]])
    ])
    cond_levels <- levels(conditions)

    if (length(cond_levels) != 2) {
        message("ssGSEA: Wilcoxon test requires exactly 2 conditions, found ",
                length(cond_levels))
        return(list(table = NULL, scores = scores, method = "ssgsea"))
    }

    grp1_idx <- which(conditions == cond_levels[1])
    grp2_idx <- which(conditions == cond_levels[2])

    pvalues <- vapply(seq_len(nrow(scores)), function(i) {
        tryCatch(
            stats::wilcox.test(scores[i, grp1_idx], scores[i, grp2_idx])$p.value,
            error = function(e) NA_real_
        )
    }, numeric(1))

    fdr   <- stats::p.adjust(pvalues, method = "fdr")
    mean1 <- rowMeans(scores[, grp1_idx, drop = FALSE], na.rm = TRUE)
    mean2 <- rowMeans(scores[, grp2_idx, drop = FALSE], na.rm = TRUE)

    results_df <- data.frame(
        pathway    = rownames(scores),
        p_value    = pvalues,
        FDR        = fdr,
        score_diff = mean2 - mean1,
        significant = !is.na(fdr) & fdr < 0.05,
        stringsAsFactors = FALSE
    )
    colnames(results_df)[colnames(results_df) == "score_diff"] <- "score_diff"
    # Add per-condition means with descriptive names
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
#' @param enrich_df data.frame with pathway, FDR columns.
#' @param top_n     Number of top pathways to show.
#' @param title     Plot title.
#' @return ggplot object.
plot_enrichment_barplot <- function(enrich_df, top_n = 20,
                                    title = "Pathway Enrichment") {
    if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)

    top_df <- utils::head(enrich_df, top_n)
    top_df$neg_log10_fdr <- -log10(pmax(top_df$FDR, 1e-20))
    top_df$pathway_short <- ifelse(
        nchar(top_df$pathway) > 50,
        paste0(substr(top_df$pathway, 1, 47), "..."),
        top_df$pathway
    )
    top_df$pathway_short <- factor(top_df$pathway_short,
                                    levels = rev(top_df$pathway_short))

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
#' @param scores       Pathway scores matrix (pathways x samples).
#' @param conditions   Factor of sample conditions.
#' @param results_df   ssGSEA results data.frame with significant flag.
#' @param max_pathways Maximum pathways to display.
#' @return ggplot object or NULL.
plot_ssgsea_boxplots <- function(scores, conditions, results_df,
                                  max_pathways = 12) {
    sig_pathways <- results_df$pathway[results_df$significant == TRUE]
    if (length(sig_pathways) == 0) return(NULL)

    sig_pathways <- utils::head(sig_pathways, max_pathways)
    sig_scores   <- scores[sig_pathways, , drop = FALSE]

    plot_data <- do.call(rbind, lapply(sig_pathways, function(pw) {
        data.frame(
            pathway   = pw,
            score     = sig_scores[pw, ],
            condition = as.character(conditions),
            stringsAsFactors = FALSE
        )
    }))

    plot_data$pathway_short <- ifelse(
        nchar(plot_data$pathway) > 40,
        paste0(substr(plot_data$pathway, 1, 37), "..."),
        plot_data$pathway
    )

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
