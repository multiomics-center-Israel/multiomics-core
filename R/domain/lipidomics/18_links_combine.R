# R/domain/lipidomics/18_links_combine.R
#
# Combine pathway- and literature-based evidence into final per-layer
# ranked tables.
#
# Inputs:
#   pathway_links_list  list of per-layer tables from score_pathway_links()
#   literature_res      list returned by run_pubmed_literature_links() with
#                       $papers and $summary; may be NULL when literature
#                       step is disabled or yielded no hits.
#
# Output per layer (long table, ranked by combined_score):
#   layer, contrast, lipid_class, lipid_id, lipid_name, lipid_logFC,
#   lipid_padj, feature_id, feature_name, feature_logFC, feature_padj,
#   sign_concordant, n_pathways, pathway_ids (semicolon-joined),
#   n_papers, top_pmid, combined_score


combine_lipid_links <- function(pathway_links_list, literature_res = NULL,
                                  pathway_weight = 1.0,
                                  literature_weight = 0.5) {
    if (is.null(pathway_links_list) || length(pathway_links_list) == 0) {
        return(NULL)
    }
    lit_summary <- literature_res$summary
    lit_papers  <- literature_res$papers

    out_list <- list()
    for (lay in names(pathway_links_list)) {
        pw <- pathway_links_list[[lay]]
        if (is.null(pw) || nrow(pw) == 0) next

        # Aggregate to one row per (lipid, feature, contrast)
        pair_key <- paste(pw$layer, pw$contrast, pw$lipid_id,
                           pw$feature_id, sep = "\x1f")

        # Identity columns (first-row representatives per pair)
        first_idx <- !duplicated(pair_key)
        keys <- pw[first_idx,
                    c("layer", "contrast", "lipid_id", "lipid_name",
                      "lipid_class", "lipid_logFC", "lipid_padj",
                      "feature_id", "feature_name", "feature_logFC",
                      "feature_padj", "sign_concordant"),
                    drop = FALSE]

        n_pw_v  <- tapply(pw$kegg_pathway_id, pair_key,
                           function(x) length(unique(x)))
        sum_v   <- tapply(pw$pathway_link_score, pair_key,
                           function(x) sum(x, na.rm = TRUE))
        pids_v  <- tapply(pw$kegg_pathway_id, pair_key,
                           function(x) paste(sort(unique(x)),
                                              collapse = ";"))

        agg <- keys
        keys_pair <- paste(keys$layer, keys$contrast, keys$lipid_id,
                            keys$feature_id, sep = "\x1f")
        agg$n_pathways        <- as.integer(n_pw_v[keys_pair])
        agg$pathway_link_score <- as.numeric(sum_v[keys_pair])
        agg$pathway_ids        <- as.character(pids_v[keys_pair])

        # Literature support
        agg$n_papers <- 0L
        agg$top_pmid <- NA_character_
        if (!is.null(lit_summary) && nrow(lit_summary) > 0) {
            lit_lay <- lit_summary[lit_summary$layer == lay, , drop = FALSE]
            if (nrow(lit_lay) > 0) {
                # Key-based lookup preserves agg row order; merge(all.x=TRUE)
                # previously reordered its result, mis-assigning n_papers.
                key_str <- paste(agg$contrast, agg$lipid_class,
                                  agg$feature_name, sep = "||")
                pap_str <- paste(lit_lay$contrast, lit_lay$lipid_class,
                                  lit_lay$feature_name, sep = "||")
                n_papers <- vapply(key_str, function(k) {
                    hits <- lit_lay$n_papers[pap_str == k]
                    if (length(hits) == 0) 0L else as.integer(hits[1])
                }, integer(1))
                agg$n_papers <- n_papers
            }
        }
        if (!is.null(lit_papers) && nrow(lit_papers) > 0) {
            lit_lay <- lit_papers[lit_papers$layer == lay, , drop = FALSE]
            if (nrow(lit_lay) > 0) {
                key_str <- paste(agg$contrast, agg$lipid_class,
                                  agg$feature_name, sep = "||")
                pap_str <- paste(lit_lay$contrast, lit_lay$lipid_class,
                                  lit_lay$feature_name, sep = "||")
                first_pmid <- vapply(key_str, function(k) {
                    hits <- lit_lay$pmid[pap_str == k]
                    if (length(hits) == 0) NA_character_ else hits[1]
                }, character(1))
                agg$top_pmid <- first_pmid
            }
        }

        agg$combined_score <- (pathway_weight * agg$pathway_link_score) +
            (literature_weight * log1p(agg$n_papers))
        agg <- agg[order(-agg$combined_score), , drop = FALSE]
        out_list[[lay]] <- agg
    }
    out_list
}
