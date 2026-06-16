# R/domain/lipidomics/16_links_scoring.R
#
# Score pathway-based links between changed lipids and DE features in the
# external omics layers.
#
# Inputs:
#   lipid_pw_long    long table from lipids_to_pathways() — one row per
#                    (sig_lipid, contrast, kegg_pathway_id)
#   ext_long_list    list of harmonized DE tables (rnaseq / proteomics /
#                    metabolomics) from load_external_omics_for_links()
#   pw_members       long table from fetch_kegg_pathway_members()
#
# Outputs:
#   list of per-layer ranked tables. Each row is a (lipid, feature, contrast,
#   pathway_id) triple plus a `pathway_link_score`. Multiple pathway hits
#   for the same lipid-feature pair are kept as separate rows; aggregation
#   happens in 18_links_combine.R.


score_pathway_links <- function(lipid_pw_long, ext_long_list, pw_members,
                                  fdr_cutoff_layer = 0.25,
                                  abs_logfc_min = 0.1) {
    if (is.null(lipid_pw_long) || nrow(lipid_pw_long) == 0) return(NULL)

    # Index pathway members by id
    pw <- pw_members
    pw_genes_by_id <- pw[pw$member_type == "gene", , drop = FALSE]
    pw_cpd_by_id   <- pw[pw$member_type == "compound", , drop = FALSE]

    out_list <- list()

    # --- RNAseq + Proteomics (gene members): match on Entrez ID --------------
    for (lay in c("rnaseq", "proteomics")) {
        ext <- ext_long_list[[lay]]
        if (is.null(ext) || nrow(ext) == 0) next

        # Filter to DE-passing features for this layer
        kept <- ext[!is.na(ext$p_adj) & ext$p_adj < fdr_cutoff_layer &
                     !is.na(ext$logFC) & abs(ext$logFC) >= abs_logfc_min, ,
                     drop = FALSE]
        if (nrow(kept) == 0) {
            message("[links] ", lay, ": no DE features pass filter")
            next
        }

        # Determine match column on the layer side: Entrez for both layers
        # (RNAseq stores Entrez in entrez_id; proteomics has Entrez via
        # wormbase_id-to-entrez? — we already have wormbase_id; rely on
        # pathway-side member_id == NCBI gene id and proteomics has
        # wormbase_id only. We approximate by joining proteomics on
        # wormbase_id matched to Entrez via the rnaseq mapping when possible.
        # Simpler approach: KEGG cel pathway member_id IS an Entrez Gene ID;
        # proteomics has uniprot+wormbase but not entrez. So we map via the
        # RNAseq table's gene_id (Wormbase) -> entrez_id, build a small
        # wb->entrez lookup, then join.

        if (lay == "rnaseq") {
            kept$.entrez_match <- as.character(kept$entrez_id)
        } else {
            # proteomics: need wormbase_id -> entrez_id via rnaseq table
            rna <- ext_long_list$rnaseq
            wb2entrez <- NULL
            if (!is.null(rna) && nrow(rna) > 0) {
                wb2entrez <- unique(rna[, c("feature_id", "entrez_id"),
                                          drop = FALSE])
                wb2entrez <- wb2entrez[!is.na(wb2entrez$entrez_id) &
                                        nzchar(wb2entrez$entrez_id), ,
                                        drop = FALSE]
            }
            if (is.null(wb2entrez) || nrow(wb2entrez) == 0) {
                # No way to map proteomics to KEGG genes — skip silently
                next
            }
            kept <- merge(kept, wb2entrez, by.x = "wormbase_id",
                          by.y = "feature_id", all.x = TRUE,
                          suffixes = c("", ".rna"))
            kept$.entrez_match <- ifelse(is.na(kept$entrez_id.rna),
                                          NA_character_,
                                          as.character(kept$entrez_id.rna))
            kept$entrez_id.rna <- NULL
        }

        kept <- kept[!is.na(kept$.entrez_match) & nzchar(kept$.entrez_match), ,
                     drop = FALSE]
        if (nrow(kept) == 0) next

        # Join lipid-pathway long with pathway-gene members, then with kept.
        # Use the entrez_id column on pw_members (added by fetcher via
        # keggConv) so the match against the layer-side Entrez ID is exact.
        join_cols <- c("kegg_pathway_id", "member_id", "member_name")
        if ("entrez_id" %in% colnames(pw_genes_by_id)) {
            join_cols <- c(join_cols, "entrez_id")
        }
        pw_join <- merge(
            lipid_pw_long,
            pw_genes_by_id[, join_cols, drop = FALSE],
            by = "kegg_pathway_id", all.x = FALSE
        )
        if (nrow(pw_join) == 0) next

        # Drop rows where the pathway gene has no Entrez mapping
        if ("entrez_id" %in% colnames(pw_join)) {
            pw_join <- pw_join[!is.na(pw_join$entrez_id) &
                                 nzchar(pw_join$entrez_id), , drop = FALSE]
            if (nrow(pw_join) == 0) next
            join_left <- c("entrez_id", "contrast")
        } else {
            join_left <- c("member_id", "contrast")
        }

        # Join layer features by Entrez ID; require matching contrast.
        joined <- merge(
            pw_join,
            kept,
            by.x = join_left,
            by.y = c(".entrez_match", "contrast"),
            all = FALSE,
            suffixes = c(".lipid", ".feat")
        )
        if (nrow(joined) == 0) next

        out_list[[lay]] <- .pkg_link_score_rows(joined, layer = lay)
    }

    # --- Metabolomics (compound members): match on KEGG cpd id --------------
    ext_m <- ext_long_list$metabolomics
    if (!is.null(ext_m) && nrow(ext_m) > 0 && nrow(pw_cpd_by_id) > 0) {
        kept <- ext_m[!is.na(ext_m$p_adj) & ext_m$p_adj < fdr_cutoff_layer &
                       !is.na(ext_m$logFC) &
                       abs(ext_m$logFC) >= abs_logfc_min &
                       !is.na(ext_m$kegg_cpd) & nzchar(ext_m$kegg_cpd), ,
                       drop = FALSE]
        if (nrow(kept) > 0) {
            pw_join <- merge(
                lipid_pw_long,
                pw_cpd_by_id[, c("kegg_pathway_id", "member_id",
                                  "member_name")],
                by = "kegg_pathway_id", all.x = FALSE
            )
            if (nrow(pw_join) > 0) {
                joined <- merge(
                    pw_join,
                    kept,
                    by.x = c("member_id", "contrast"),
                    by.y = c("kegg_cpd", "contrast"),
                    all = FALSE,
                    suffixes = c(".lipid", ".feat")
                )
                if (nrow(joined) > 0) {
                    out_list$metabolomics <- .pkg_link_score_rows(
                        joined, layer = "metabolomics"
                    )
                }
            }
        }
    }

    out_list
}


# Compute the per-row link score and tidy column names.
.pkg_link_score_rows <- function(joined, layer) {
    # Sign concordance bonus when lipid and feature move the same direction
    sign_conc <- sign(joined$lipid_logFC) == sign(joined$logFC)
    sign_conc[is.na(sign_conc)] <- FALSE
    score <- abs(joined$lipid_logFC) * abs(joined$logFC) *
             ifelse(sign_conc, 1.0, 0.5)

    data.frame(
        layer                = layer,
        contrast             = joined$contrast,
        lipid_id             = joined$feature_id.lipid %||%
                                  joined$feature_id,
        lipid_name           = joined$feature_name.lipid %||%
                                  joined$feature_name,
        lipid_class          = joined$ClassKey,
        lipid_logFC          = joined$lipid_logFC,
        lipid_padj           = joined$lipid_padj,
        kegg_pathway_id      = joined$kegg_pathway_id,
        pathway_name         = joined$pathway_name,
        feature_id           = joined$feature_id.feat %||% joined$feature_id,
        feature_name         = joined$feature_name.feat %||%
                                  joined$feature_name,
        feature_kegg_member  = joined$member_name,
        feature_logFC        = joined$logFC,
        feature_padj         = joined$p_adj,
        sign_concordant      = sign_conc,
        pathway_link_score   = score,
        stringsAsFactors = FALSE
    )
}
