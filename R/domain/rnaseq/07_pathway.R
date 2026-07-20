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
#' @return data.frame: FeatureID, log2FC, simple_log2FC, neg_log10_pvalue, padj, pathways
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
    # Naive group-mean log-ratio companion (present for the DESeq2 path); carried
    # so downstream consumers (e.g. the multiomics RNA-protein plots) can use it.
    simple_log2FC = if ("simple_log2FC" %in% names(de_table)) de_table$simple_log2FC else NA_real_,
    neg_log10_pvalue = neg_log10_p,
    padj = de_table$padj,
    stringsAsFactors = FALSE
  )
  
  pathway_memberships <- character(nrow(volcano_df))
  
  # Flatten: collect all data.frames from the nested structure
  all_result_dfs <- list()
  for (top_name in names(pathway_results)) {
    top_elem <- pathway_results[[top_name]]
    if (is.data.frame(top_elem)) {
      all_result_dfs <- c(all_result_dfs, list(top_elem))
    } else if (is.list(top_elem)) {
      for (sub_name in names(top_elem)) {
        sub_elem <- top_elem[[sub_name]]
        if (is.data.frame(sub_elem)) {
          all_result_dfs <- c(all_result_dfs, list(sub_elem))
        }
      }
    }
  }
  
  enriched_sets <- list()
  for (res in all_result_dfs) {
    if ("leadingEdge" %in% colnames(res)) {
      sig <- res[!is.na(res$padj) & res$padj < 0.05, , drop = FALSE]
      for (r in seq_len(nrow(sig))) {
        pw_label <- if ("pathway_name" %in% colnames(sig) && nzchar(sig$pathway_name[r] %||% "")) {
          sig$pathway_name[r]
        } else {
          sig$pathway[r]
        }
        genes <- trimws(strsplit(as.character(sig$leadingEdge[[r]]), ",")[[1]])
        if (length(genes) > 0) enriched_sets[[pw_label]] <- genes
      }
    } else if ("geneID" %in% colnames(res)) {
      sig <- res[!is.na(res$p.adjust) & res$p.adjust < 0.05, , drop = FALSE]
      for (r in seq_len(nrow(sig))) {
        pw_label <- sig$Description[r] %||% sig$ID[r]
        genes <- unlist(strsplit(as.character(sig$geneID[r]), "/"))
        if (length(genes) > 0) enriched_sets[[pw_label]] <- genes
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