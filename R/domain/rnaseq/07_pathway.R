#' RNA-seq Pathway / Enrichment Analysis — Domain Helpers
#'
#' Domain-specific pathway utilities for RNA-seq.
#' Core enrichment functions (run_pathway_analysis, save_pathway_results,
#' generate_pathway_plots, lookup_go_term_names, add_pathway_names) live
#' in R/core/09_enrichment.R and are shared across omics types.

# ==============================================================================
# PATHWAY VOLCANO DATA
# ==============================================================================

#' Build pathway-colored volcano data
#'
#' Creates a data frame tagging each gene with its enriched pathway memberships,
#' where genes are colored by pathway membership on a volcano plot.
#'
#' @param de_table DE table with FeatureID, log2FoldChange, pvalue, padj columns
#' @param pathway_results Pathway results from run_pathway_analysis()
#' @return data.frame: FeatureID, log2FC, neg_log10_pvalue, padj, pathways
build_pathway_volcano_data <- function(de_table, pathway_results) {
  if (is.null(de_table) || nrow(de_table) == 0) return(NULL)
  if (is.null(pathway_results)) return(NULL)
  
  # Cap -log10(p) at the 99.5th percentile to remove extreme outliers
  # that would compress the bulk of the data to the bottom of the plot
  neg_log10_p <- -log10(de_table$pvalue + 1e-300)
  cap <- quantile(neg_log10_p[is.finite(neg_log10_p)], 0.995, na.rm = TRUE)
  cap <- max(cap, 10)
  neg_log10_p <- pmin(neg_log10_p, cap)
  
  volcano_df <- data.frame(
    FeatureID = de_table$FeatureID,
    log2FC = de_table$log2FoldChange,
    neg_log10_pvalue = neg_log10_p,
    padj = de_table$padj,
    stringsAsFactors = FALSE
  )
  
  pathway_memberships <- character(nrow(volcano_df))
  
  # Flatten: collect all data.frames from the nested structure, keeping the name
  # of the result they came from. Since gene-set collections are scored one per
  # GMT file, two of them can carry the same pathway id and label; without the
  # source, the later one's leading edge would replace the earlier one's under a
  # single key and that membership would vanish from the export.
  all_result_dfs <- list()
  result_sources <- character(0)
  for (top_name in names(pathway_results)) {
    top_elem <- pathway_results[[top_name]]
    if (is.data.frame(top_elem)) {
      all_result_dfs <- c(all_result_dfs, list(top_elem))
      result_sources <- c(result_sources, top_name)
    } else if (is.list(top_elem)) {
      for (sub_name in names(top_elem)) {
        sub_elem <- top_elem[[sub_name]]
        if (is.data.frame(sub_elem)) {
          all_result_dfs <- c(all_result_dfs, list(sub_elem))
          result_sources <- c(result_sources, sub_name)
        }
      }
    }
  }

  enriched_sets <- list()
  for (i in seq_along(all_result_dfs)) {
    res <- all_result_dfs[[i]]
    src <- result_sources[i]
    # Only a label already taken by different genes is qualified by its source,
    # so labels stay as they were except where one would otherwise be lost.
    set_key <- function(label, genes) {
      if (!is.null(enriched_sets[[label]]) && !identical(enriched_sets[[label]], genes)) {
        paste0(label, " (", src, ")")
      } else {
        label
      }
    }
    if ("leadingEdge" %in% colnames(res)) {
      sig <- res[!is.na(res$padj) & res$padj < 0.05, , drop = FALSE]
      for (r in seq_len(nrow(sig))) {
        pw_label <- if ("pathway_name" %in% colnames(sig) && nzchar(sig$pathway_name[r] %||% "")) {
          sig$pathway_name[r]
        } else {
          sig$pathway[r]
        }
        genes <- trimws(strsplit(as.character(sig$leadingEdge[[r]]), ",")[[1]])
        if (length(genes) > 0) enriched_sets[[set_key(pw_label, genes)]] <- genes
      }
    } else if ("geneID" %in% colnames(res)) {
      sig <- res[!is.na(res$p.adjust) & res$p.adjust < 0.05, , drop = FALSE]
      for (r in seq_len(nrow(sig))) {
        pw_label <- sig$Description[r] %||% sig$ID[r]
        genes <- unlist(strsplit(as.character(sig$geneID[r]), "/"))
        if (length(genes) > 0) enriched_sets[[set_key(pw_label, genes)]] <- genes
      }
    }
  }
  
  if (length(enriched_sets) == 0) return(volcano_df)
  
  # Tag each gene with its pathway memberships (top 5 per gene)
  for (i in seq_len(nrow(volcano_df))) {
    fid <- volcano_df$FeatureID[i]
    member_of <- names(enriched_sets)[vapply(enriched_sets, function(gs) fid %in% gs, logical(1))]
    if (length(member_of) > 5) member_of <- member_of[1:5]
    pathway_memberships[i] <- paste(member_of, collapse = "; ")
  }
  
  volcano_df$pathways <- pathway_memberships
  volcano_df
}

# ==============================================================================
# GSEA SETTINGS
# ==============================================================================

#' GSEA engines the offline RNA path can run
#' @noRd
.RNA_GSEA_ENGINES <- c("fgsea", "clusterProfiler")

#' Resolve the RNA GSEA settings
#'
#' The one reader of the offline RNA GSEA engine and gene-set size limits, and
#' of the enrichment seed, which the online path also takes from here. The
#' online path keeps its own size limits (`modes.rna.pathway.min_size` /
#' `max_size`, default 10 / 500), read where it calls run_pathway_analysis().
#'
#' - `engine`: `modes.rna.enrichment.gsea_engine`, default "fgsea". "fgsea"
#'   runs the shared fgsea call (run_fgsea_core()) and keeps every tested set;
#'   "clusterProfiler" reproduces earlier offline runs (clusterProfiler::GSEA,
#'   non-significant sets dropped).
#' - `min_size` / `max_size`: `modes.rna.enrichment.gsea.min_size` /
#'   `max_size`, defaulting to the limits of the earlier offline runs: 4, and
#'   NULL, which run_gsea_all() reads as the number of unique genes in each
#'   database's TERM2GENE table (no upper limit in practice). The
#'   clusterProfiler engine always uses those defaults.
#' - `seed`: `params.seed`, default 1.
#'
#' @param config Full pipeline config.
#' @return list(engine, min_size, max_size, seed); `max_size` may be NULL.
resolve_rna_gsea_settings <- function(config) {
  rna  <- config$modes$rna %||% list()
  enr  <- rna$enrichment %||% list()
  # [[ ]] not $: enr$gsea would partially match enr$gsea_engine when the
  # gsea block is absent and hand back the engine string.
  gsea <- enr[["gsea", exact = TRUE]] %||% list()

  engine <- enr$gsea_engine %||% "fgsea"
  # Case-insensitive, as validate_rna_config() accepts it.
  hit <- if (is.character(engine) && length(engine) == 1L) {
    match(tolower(engine), tolower(.RNA_GSEA_ENGINES))
  } else NA_integer_
  if (is.na(hit)) {
    stop("modes.rna.enrichment.gsea_engine must be one of: ",
         paste(.RNA_GSEA_ENGINES, collapse = ", "), " (got: ",
         paste(format(engine), collapse = ", "), ").", call. = FALSE)
  }
  engine <- .RNA_GSEA_ENGINES[hit]

  list(
    engine   = engine,
    min_size = gsea$min_size %||% 4,
    max_size = gsea$max_size,
    seed     = config$params$seed %||% 1L
  )
}
