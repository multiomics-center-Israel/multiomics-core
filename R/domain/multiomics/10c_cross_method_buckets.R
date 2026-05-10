# =============================================================================
# Cross-Method Bucket Enrichment (MOFA Factor x DIABLO Component story)
# =============================================================================
# For a chosen MOFA factor + DIABLO component, partition each omics layer's
# top features into three buckets:
#   - MOFA-only   (top in MOFA factor, not in DIABLO component)
#   - DIABLO-only (top in DIABLO component, not in MOFA factor)
#   - Shared      (top in both)
#
# For each bucket, run hypergeometric ORA against the configured KEGG GMT
# (and optional GO GMT), then identify pathways that are enriched across
# multiple omics layers within the same bucket -- the cross-omics story.
#
# Outputs (under <out_dir>/buckets/):
#   buckets_<view>.csv                 -- per-feature bucket assignment
#   bucket_enrichment_<view>_<bucket>_<set>.csv -- ORA results
#   cross_omics_pathway_convergence_<bucket>_<set>.csv -- pathways x omics
#   plots/buckets_overlap_<view>.png   -- 2-set Venn (MOFA vs DIABLO)
#   plots/cross_omics_convergence_<set>.png -- pathway x omics x bucket heatmap
# =============================================================================

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

# Ordered fallback for picking an external ID column (highest priority first)
.bucket_external_id <- function(df) {
    pick <- intersect(c("original_name", "name", "ext_id", "feature"),
                      colnames(df))[1]
    df[[pick]]
}

# Build a per-view GENE_xxx -> external id lookup from MOFA + DIABLO files.
# Returns a named character vector. Names = MOFA/DIABLO feature ids (e.g.
# "GENE_3262"); values = external IDs (e.g. "GL50803_006289" or "KAE…").
.bucket_view_id_map <- function(view, mo_root) {
    map <- character(0)
    diablo_top  <- file.path(mo_root, "integration", "diablo",
                             sprintf("diablo_top_features_%s.csv", view))
    mofa_w      <- file.path(mo_root, "integration", "mofa",
                             sprintf("mofa_weights_%s.csv", view))
    mofa_top    <- file.path(mo_root, "integration", "mofa",
                             "mofa_top_features.csv")
    add_pairs <- function(feats, names) {
        ok <- !is.na(names) & nzchar(names) & !is.na(feats) & nzchar(feats)
        stats::setNames(names[ok], feats[ok])
    }
    if (file.exists(diablo_top)) {
        d <- utils::read.csv(diablo_top, stringsAsFactors = FALSE)
        if (all(c("feature", "original_name") %in% colnames(d))) {
            map <- c(map, add_pairs(d$feature, d$original_name))
        }
    }
    if (file.exists(mofa_w)) {
        m <- utils::read.csv(mofa_w, stringsAsFactors = FALSE)
        if (all(c("feature_id", "original_name") %in% colnames(m))) {
            fid <- sub("_(transcriptomics|proteomics|metabolomics)$", "",
                       m$feature_id)
            map <- c(map, add_pairs(fid, m$original_name))
        }
    }
    if (file.exists(mofa_top)) {
        mt <- utils::read.csv(mofa_top, stringsAsFactors = FALSE)
        mt <- mt[mt$view == view, , drop = FALSE]
        if (all(c("feature_id", "original_name") %in% colnames(mt))) {
            fid <- sub("_(transcriptomics|proteomics|metabolomics)$", "",
                       mt$feature_id)
            map <- c(map, add_pairs(fid, mt$original_name))
        }
    }
    map <- map[!duplicated(names(map))]
    # Normalize to gene-id space when working with transcriptomics: some
    # upstream files emit protein KAE… IDs even on the transcriptomics view.
    # Use gene_protein_mapping.csv to translate KAE… → GL50803_ so all
    # transcriptomics ext_ids live in the same namespace. Also recover the
    # GL50803_ identity for genes whose `original_name` was dropped during
    # integration writing — by convention GENE_n corresponds to row n of
    # gene_protein_mapping.csv (verified against the resolved rows).
    if (view == "transcriptomics") {
        gpm_f <- file.path(mo_root, "gene_protein_mapping.csv")
        if (file.exists(gpm_f)) {
            gpm <- utils::read.csv(gpm_f, stringsAsFactors = FALSE)
            if (all(c("gene_id", "protein_id") %in% colnames(gpm))) {
                # KAE… → GL50803_
                if (length(map) > 0) {
                    lut <- stats::setNames(gpm$gene_id, gpm$protein_id)
                    kae <- grepl("^KAE", map)
                    hit <- map[kae] %in% names(lut)
                    if (any(hit)) {
                        new_vals <- map[kae]
                        new_vals[hit] <- unname(lut[new_vals[hit]])
                        map[kae] <- new_vals
                    }
                }
                # Row-index recovery: GENE_n -> gpm$gene_id[n]. Add for any
                # GENE_n in [1..nrow(gpm)] not already in the map, and overwrite
                # entries that are still GENE_xxx-shaped (didn't resolve above).
                n <- nrow(gpm)
                ord <- function(x) {
                    suppressWarnings(as.integer(sub("^GENE_(\\d+).*", "\\1", x)))
                }
                # Any names whose value is still "GENE_..."-like → recover
                still_gene <- grepl("^GENE_", map)
                if (any(still_gene)) {
                    keys <- names(map)[still_gene]
                    idx <- ord(keys)
                    valid <- !is.na(idx) & idx >= 1 & idx <= n
                    if (any(valid)) map[still_gene][valid] <- gpm$gene_id[idx[valid]]
                }
                # Add map entries for GENE_1..GENE_n that aren't in the map yet
                missing_keys <- setdiff(sprintf("GENE_%d", seq_len(n)), names(map))
                if (length(missing_keys) > 0) {
                    add <- stats::setNames(gpm$gene_id[ord(missing_keys)], missing_keys)
                    map <- c(map, add)
                }
            }
        }
    }
    map
}

# Translate any remaining GENE_xxx ids using the view id-map. ids that already
# look like external IDs are left unchanged.
.bucket_translate_via_map <- function(ids, id_map) {
    out <- ids
    needs <- grepl("^GENE_", out)
    keys <- sub("_(transcriptomics|proteomics|metabolomics)$", "", out)
    hit <- needs & keys %in% names(id_map)
    out[hit] <- unname(id_map[keys[hit]])
    out
}

# Build a display-name table for the report. Joins the bucket ext_ids against
# the concordance table (which carries gene symbols + protein descriptions).
.bucket_display_names <- function(ext_ids, view, mo_root) {
    res <- data.frame(ext_id = ext_ids,
                      display_name = ext_ids,
                      symbol = NA_character_,
                      description = NA_character_,
                      stringsAsFactors = FALSE)
    conc_f <- file.path(mo_root, "concordance",
                        "concordance_transcriptomics_vs_proteomics.csv")
    if (!file.exists(conc_f)) return(res)
    conc <- utils::read.csv(conc_f, stringsAsFactors = FALSE)
    pick_col <- function(df, names) {
        intersect(names, colnames(df))[1]
    }
    sym_col  <- pick_col(conc, c("gene_symbol", "Genes", "Gene_symbol"))
    name_col <- pick_col(conc, c("Protein.Names", "Protein_Names"))
    desc_col <- pick_col(conc, c("First.Protein.Description",
                                 "Protein_Description", "description"))
    if (view == "transcriptomics" && "gene_id" %in% colnames(conc)) {
        m <- match(res$ext_id, conc$gene_id)
        if (!is.null(sym_col))  res$symbol      <- conc[[sym_col]][m]
        if (!is.null(desc_col)) res$description <- conc[[desc_col]][m]
    } else if (view == "proteomics" && "protein_id" %in% colnames(conc)) {
        m <- match(res$ext_id, conc$protein_id)
        if (!is.null(sym_col))  res$symbol      <- conc[[sym_col]][m]
        if (!is.null(desc_col)) res$description <- conc[[desc_col]][m]
    }
    # Compose human-readable display. Concordance often packs gene_id and
    # description into a single "Genes" field already (e.g. "GL50803_006289 |
    # IscU"); avoid duplicating it.
    has_sym  <- !is.na(res$symbol) & nzchar(res$symbol) & res$symbol != "|"
    has_desc <- !is.na(res$description) & nzchar(res$description)
    sym_has_desc <- has_sym & has_desc &
        mapply(function(s, d) grepl(d, s, fixed = TRUE),
               res$symbol, res$description, USE.NAMES = FALSE)
    res$display_name <- ifelse(sym_has_desc, res$symbol,
        ifelse(has_sym & has_desc,
               sprintf("%s | %s", res$symbol, res$description),
        ifelse(has_sym,  res$symbol,
        ifelse(has_desc, sprintf("%s | %s", res$ext_id, res$description),
               res$ext_id))))
    res
}

# Build a per-feature DE stats lookup for one view (log2FC, pvalue, padj, sig).
# Sources:
#   transcriptomics → concordance_transcriptomics_vs_proteomics.csv
#                     (logFC_rna / padj_rna / pvalue.<contrast>_rna)
#   proteomics      → concordance_transcriptomics_vs_proteomics.csv
#                     (logFC_prot / padj_prot / pvalue.imputs.<contrast>)
#   metabolomics    → upstream payload rds (de_final_table); fall back to
#                     a metabolomics DE CSV if configured.
.bucket_de_stats <- function(ext_ids, view, mo_root, config) {
    res <- data.frame(ext_id = ext_ids,
                      log2FC  = NA_real_,
                      pvalue  = NA_real_,
                      padj    = NA_real_,
                      sig     = NA,
                      stringsAsFactors = FALSE)
    pick_col <- function(df, pat) {
        cols <- grep(pat, colnames(df), value = TRUE, ignore.case = TRUE)
        if (length(cols) == 0) NULL else cols[1]
    }
    if (view %in% c("transcriptomics", "proteomics")) {
        conc_f <- file.path(mo_root, "concordance",
                            "concordance_transcriptomics_vs_proteomics.csv")
        if (!file.exists(conc_f)) return(res)
        conc <- utils::read.csv(conc_f, stringsAsFactors = FALSE,
                                check.names = FALSE)
        key_col <- if (view == "transcriptomics") "gene_id" else "protein_id"
        if (!key_col %in% colnames(conc)) return(res)
        m <- match(res$ext_id, conc[[key_col]])
        if (view == "transcriptomics") {
            fc_col   <- pick_col(conc, "^logFC_rna$")
            padj_col <- pick_col(conc, "^padj_rna$")
            pval_col <- pick_col(conc, "^pvalue\\.[^_]+_rna|^pvalue\\..*_vs_.*$")
        } else {
            fc_col   <- pick_col(conc, "^logFC_prot$|^linearFC\\.imputs")
            padj_col <- pick_col(conc, "^padj_prot$|^padj\\.imputs")
            pval_col <- pick_col(conc, "^pvalue\\.imputs|^pvalue_prot")
        }
        if (!is.null(fc_col))   res$log2FC <- suppressWarnings(as.numeric(conc[[fc_col]][m]))
        if (!is.null(padj_col)) res$padj   <- suppressWarnings(as.numeric(conc[[padj_col]][m]))
        if (!is.null(pval_col)) res$pvalue <- suppressWarnings(as.numeric(conc[[pval_col]][m]))
    } else if (view == "metabolomics") {
        # Try upstream payload first
        metab_payload <- config$inputs$payloads$metabolomics %||%
                         config$global$inputs$payloads$metabolomics
        de_df <- NULL
        if (!is.null(metab_payload) && file.exists(metab_payload)) {
            tryCatch({
                p <- readRDS(metab_payload)
                # Prefer the full de_stats (covers all compounds), fall back
                # to de_final_table (filtered to significant only).
                if (!is.null(p$de_stats) && is.data.frame(p$de_stats)) {
                    de_df <- p$de_stats
                } else if (!is.null(p$de_final_table)) {
                    de_df <- p$de_final_table
                }
            }, error = function(e) NULL)
        }
        # Fallback: configured CSV path
        if (is.null(de_df)) {
            de_csv <- config$modes$multiomics$consensus$metab_de_csv %||%
                      config$modes$multiomics$enrichment$metab_de_csv
            if (!is.null(de_csv) && file.exists(de_csv)) {
                de_df <- utils::read.csv(de_csv, stringsAsFactors = FALSE,
                                          check.names = FALSE)
            }
        }
        if (!is.null(de_df)) {
            id_col <- intersect(c("feature_id", "FeatureID", "HMDB_ID", "id"),
                                colnames(de_df))[1]
            name_col <- intersect(c("Name", "name", "metabolite", "compound"),
                                  colnames(de_df))[1]
            if (!is.null(id_col) || !is.null(name_col)) {
                # Try id then name
                m <- if (!is.null(id_col)) match(res$ext_id, de_df[[id_col]]) else rep(NA_integer_, nrow(res))
                if (!is.null(name_col)) {
                    nm_match <- match(res$ext_id, de_df[[name_col]])
                    m[is.na(m)] <- nm_match[is.na(m)]
                }
                fc_col   <- pick_col(de_df, "linearFC|log2FoldChange|logFC")
                padj_col <- pick_col(de_df, "adj\\.P\\.Val|padj")
                pval_col <- pick_col(de_df, "^P\\.Value|pvalue|^pvalue\\.")
                if (!is.null(fc_col)) {
                    fc <- suppressWarnings(as.numeric(de_df[[fc_col]][m]))
                    # linearFC is fold-change, not log2; convert if needed
                    if (grepl("linearFC", fc_col, ignore.case = TRUE)) {
                        fc <- ifelse(fc > 0, log2(fc),
                                     ifelse(fc < 0, -log2(-fc), 0))
                    }
                    res$log2FC <- fc
                }
                if (!is.null(padj_col)) res$padj   <- suppressWarnings(as.numeric(de_df[[padj_col]][m]))
                if (!is.null(pval_col)) res$pvalue <- suppressWarnings(as.numeric(de_df[[pval_col]][m]))
            }
        }
    }
    res$sig <- !is.na(res$padj) & res$padj < 0.05
    res
}

# Build a feature -> external id lookup using available files. Returns a
# named character vector (names = MOFA/DIABLO feature ids, e.g. GENE_3182).
.bucket_feature_id_map <- function(view, mo_root) {
    map <- character(0)
    diablo_f <- file.path(mo_root, "integration", "diablo",
                          sprintf("diablo_top_features_%s.csv", view))
    mofa_w_f <- file.path(mo_root, "integration", "mofa",
                          sprintf("mofa_weights_%s.csv", view))
    if (file.exists(diablo_f)) {
        d <- utils::read.csv(diablo_f, stringsAsFactors = FALSE)
        if (nrow(d) > 0 && all(c("feature", "original_name") %in% colnames(d))) {
            ok <- !is.na(d$original_name) & nzchar(d$original_name)
            map <- c(map, stats::setNames(d$original_name[ok], d$feature[ok]))
        }
    }
    if (file.exists(mofa_w_f)) {
        m <- utils::read.csv(mofa_w_f, stringsAsFactors = FALSE)
        # MOFA feature_ids may carry a "_<view>" suffix
        if (nrow(m) > 0 && "feature_id" %in% colnames(m) && "original_name" %in% colnames(m)) {
            fid <- sub("_(transcriptomics|proteomics|metabolomics)$", "", m$feature_id)
            ok <- !is.na(m$original_name) & nzchar(m$original_name)
            map <- c(map, stats::setNames(m$original_name[ok], fid[ok]))
        }
    }
    map[!duplicated(names(map))]
}

# Translate KAE… protein IDs to GL50803_ gene IDs via gene_protein_mapping.csv
.bucket_protein_to_gene <- function(ids, gpm_df) {
    if (is.null(gpm_df) || nrow(gpm_df) == 0) return(ids)
    if (!all(c("gene_id", "protein_id") %in% colnames(gpm_df))) return(ids)
    lut <- stats::setNames(gpm_df$gene_id, gpm_df$protein_id)
    out <- ids
    hit <- ids %in% names(lut)
    out[hit] <- lut[ids[hit]]
    out
}

# Reuse load_gmt_file from 07_enrichment.R if available, else inline
.bucket_load_gmt <- function(gmt_file) {
    if (exists("load_gmt_file", mode = "function")) return(load_gmt_file(gmt_file))
    if (!file.exists(gmt_file)) return(NULL)
    lines <- readLines(gmt_file, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]
    out <- list(); nm <- list()
    for (l in lines) {
        f <- strsplit(l, "\t")[[1]]
        if (length(f) < 3) next
        out[[f[1]]] <- f[3:length(f)]
        nm[[f[1]]] <- f[2]
    }
    attr(out, "pathway_names") <- nm
    out
}

# Hypergeometric ORA: returns data.frame
.bucket_ora <- function(hits, universe, gmt, min_size = 3, max_size = 500) {
    hits     <- unique(hits[!is.na(hits) & nzchar(hits)])
    universe <- unique(universe[!is.na(universe) & nzchar(universe)])
    if (length(hits) < 2 || length(universe) < 10 || length(gmt) == 0) return(NULL)
    pw_names <- attr(gmt, "pathway_names")
    N <- length(universe); k <- length(hits)
    rows <- list()
    for (pw in names(gmt)) {
        in_pw <- intersect(gmt[[pw]], universe)
        m <- length(in_pw)
        if (m < min_size || m > max_size) next
        x <- length(intersect(hits, in_pw))
        if (x < 1) next
        # P(X >= x) under hypergeometric (white = m, black = N-m, draws = k)
        p <- stats::phyper(x - 1, m, N - m, k, lower.tail = FALSE)
        rows[[length(rows) + 1L]] <- data.frame(
            pathway_id   = pw,
            pathway_name = if (!is.null(pw_names)) pw_names[[pw]] %||% pw else pw,
            set_size     = m,
            hits_in_set  = x,
            n_hits       = k,
            n_universe   = N,
            pvalue       = p,
            gene_ratio   = sprintf("%d/%d", x, k),
            bg_ratio     = sprintf("%d/%d", m, N),
            leading_ids  = paste(intersect(hits, in_pw), collapse = ";"),
            stringsAsFactors = FALSE
        )
    }
    if (length(rows) == 0) return(NULL)
    out <- do.call(rbind, rows)
    out$padj <- stats::p.adjust(out$pvalue, method = "BH")
    out[order(out$pvalue), ]
}

# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------

#' Compute MOFA/DIABLO buckets and enrich per omics, then find cross-omics
#' pathway convergence within each bucket.
#'
#' @param mo_root multiomics output root (contains integration/, consensus/)
#' @param config full pipeline config
#' @param out_dir directory where buckets/ subdir will be written. Defaults to
#'   <mo_root>/consensus/consensus/buckets if NULL.
#' @return list with buckets, enrichments, convergence (also written to disk)
run_cross_method_buckets <- function(mo_root, config, out_dir = NULL) {
    message("=== Cross-Method Bucket Enrichment ===")

    cfg_b <- config$modes$multiomics$consensus %||% list()
    factor_pick    <- cfg_b$bucket_factor    %||% "Factor1"
    component_pick <- cfg_b$bucket_component %||% "comp1"
    top_n          <- cfg_b$bucket_top_n     %||% 50L
    do_enrich      <- cfg_b$run_bucket_enrichment %||% TRUE
    kegg_gmt       <- config$modes$multiomics$enrichment$kegg_gmt_file
    proteins_gmt   <- config$modes$multiomics$enrichment$proteins_gmt_file
    go_gmt         <- config$modes$multiomics$enrichment$go_gmt_file

    if (is.null(out_dir)) out_dir <- file.path(mo_root, "consensus", "consensus", "buckets")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    plots_dir <- file.path(out_dir, "plots")
    dir.create(plots_dir, showWarnings = FALSE)

    diablo_dir <- file.path(mo_root, "integration", "diablo")
    mofa_dir   <- file.path(mo_root, "integration", "mofa")
    if (!dir.exists(diablo_dir) || !dir.exists(mofa_dir)) {
        message("  Need both DIABLO and MOFA outputs; skipping bucket analysis.")
        return(NULL)
    }

    # --- Read MOFA top features (already saved by MOFA module) ---
    mofa_top_f <- file.path(mofa_dir, "mofa_top_features.csv")
    if (!file.exists(mofa_top_f)) {
        message("  mofa_top_features.csv missing; skipping.")
        return(NULL)
    }
    mofa_all <- utils::read.csv(mofa_top_f, stringsAsFactors = FALSE)
    if (!factor_pick %in% mofa_all$factor) {
        message("  Requested factor '", factor_pick, "' not found in MOFA results; skipping.")
        return(NULL)
    }

    # --- Gene-protein mapping (proteomics → gene_id) ---
    gpm_f <- file.path(mo_root, "gene_protein_mapping.csv")
    gpm <- if (file.exists(gpm_f)) utils::read.csv(gpm_f, stringsAsFactors = FALSE) else NULL

    views <- intersect(c("transcriptomics", "proteomics", "metabolomics"),
                       unique(mofa_all$view))
    if (length(views) == 0) views <- c("transcriptomics", "proteomics", "metabolomics")

    bucket_tables <- list()
    enrich_tables <- list()  # nested: [[set]][[view]][[bucket]] -> df

    for (view in views) {
        diablo_f <- file.path(diablo_dir,
                              sprintf("diablo_top_features_%s.csv", view))
        if (!file.exists(diablo_f)) {
            message("  No DIABLO file for ", view, "; skipping view.")
            next
        }
        diablo_v <- utils::read.csv(diablo_f, stringsAsFactors = FALSE)
        if (!component_pick %in% diablo_v$component) {
            message("  Component '", component_pick, "' not in DIABLO ", view, "; skipping view.")
            next
        }

        # Restrict to chosen factor / component
        m <- mofa_all[mofa_all$factor    == factor_pick   & mofa_all$view == view, ,
                      drop = FALSE]
        d <- diablo_v[diablo_v$component == component_pick, , drop = FALSE]

        # Top-N by absolute weight / loading (already sorted, but enforce)
        if ("abs_weight"  %in% colnames(m)) m <- m[order(-m$abs_weight), ]
        if ("abs_loading" %in% colnames(d)) d <- d[order(-d$abs_loading), ]
        m <- utils::head(m, top_n)
        d <- utils::head(d, top_n)

        # External IDs (e.g. GL50803_… or HMDB0000…). MOFA's transcriptomics
        # original_name is often NA, so back-fill via a per-view id map built
        # from MOFA weights + DIABLO loadings.
        id_map <- .bucket_view_id_map(view, mo_root)
        m_ids <- .bucket_external_id(m)
        d_ids <- .bucket_external_id(d)
        # Where original_name was NA, fall back to feature_id and translate
        if ("feature_id" %in% colnames(m)) {
            blank <- is.na(m_ids) | !nzchar(m_ids)
            m_ids[blank] <- sub("_(transcriptomics|proteomics|metabolomics)$",
                                "", m$feature_id[blank])
        }
        m_ids <- .bucket_translate_via_map(m_ids, id_map)
        d_ids <- .bucket_translate_via_map(d_ids, id_map)

        shared      <- intersect(m_ids, d_ids)
        mofa_only   <- setdiff(m_ids, d_ids)
        diablo_only <- setdiff(d_ids, m_ids)

        bt <- data.frame(
            view        = view,
            ext_id      = c(shared, mofa_only, diablo_only),
            bucket      = c(rep("Shared", length(shared)),
                            rep("MOFA_only", length(mofa_only)),
                            rep("DIABLO_only", length(diablo_only))),
            mofa_weight = NA_real_,
            diablo_loading = NA_real_,
            stringsAsFactors = FALSE
        )
        if (nrow(m) > 0) {
            mw <- stats::setNames(m$weight  %||% m$abs_weight,  m_ids)
            bt$mofa_weight[bt$ext_id %in% names(mw)] <- mw[bt$ext_id[bt$ext_id %in% names(mw)]]
        }
        if (nrow(d) > 0) {
            dl <- stats::setNames(d$loading %||% d$abs_loading, d_ids)
            bt$diablo_loading[bt$ext_id %in% names(dl)] <- dl[bt$ext_id[bt$ext_id %in% names(dl)]]
        }
        # Attach human-readable symbol + description from concordance table
        dn <- .bucket_display_names(bt$ext_id, view, mo_root)
        bt$symbol       <- dn$symbol
        bt$description  <- dn$description
        bt$display_name <- dn$display_name
        # Attach per-layer DE stats (log2FC, pvalue, padj, sig)
        de <- .bucket_de_stats(bt$ext_id, view, mo_root, config)
        bt$log2FC <- de$log2FC
        bt$pvalue <- de$pvalue
        bt$padj   <- de$padj
        bt$sig    <- de$sig
        # Reorder for readability
        col_order <- c("view", "bucket", "ext_id", "symbol", "description",
                       "display_name", "mofa_weight", "diablo_loading",
                       "log2FC", "pvalue", "padj", "sig")
        bt <- bt[, intersect(col_order, colnames(bt)), drop = FALSE]
        utils::write.csv(bt, file.path(out_dir, sprintf("buckets_%s.csv", view)),
                         row.names = FALSE)
        bucket_tables[[view]] <- bt

        message(sprintf("  %s: Shared=%d  MOFA-only=%d  DIABLO-only=%d",
                        view, length(shared), length(mofa_only), length(diablo_only)))

        # Simple 2-set Venn-style bar (no extra package required)
        png(file.path(plots_dir, sprintf("buckets_overlap_%s.png", view)),
            width = 700, height = 420, res = 110)
        op <- graphics::par(mar = c(4, 8, 3, 1))
        on.exit(graphics::par(op), add = TRUE)
        graphics::barplot(
            c(MOFA_only = length(mofa_only),
              Shared    = length(shared),
              DIABLO_only = length(diablo_only)),
            horiz = TRUE, las = 1,
            col = c("#4C9F70", "#7B5BA1", "#E07B39"),
            main = sprintf("%s — MOFA %s vs DIABLO %s (top %d each)",
                           view, factor_pick, component_pick, top_n),
            xlab = "n features"
        )
        dev.off()

        if (!do_enrich) next

        # Build universe (all measured features in this view, in external ID space)
        mofa_w_f <- file.path(mofa_dir, sprintf("mofa_weights_%s.csv", view))
        univ_ids <- character(0)
        if (file.exists(mofa_w_f)) {
            mw_all <- utils::read.csv(mofa_w_f, stringsAsFactors = FALSE)
            if ("original_name" %in% colnames(mw_all)) {
                univ_ids <- unique(stats::na.omit(mw_all$original_name))
            }
            # If MOFA weights lack original_name (e.g. transcriptomics), use the
            # feature_id minus suffix and translate via gpm later.
            if (length(univ_ids) < 10 && "feature_id" %in% colnames(mw_all)) {
                univ_ids <- unique(sub("_(transcriptomics|proteomics|metabolomics)$",
                                       "", mw_all$feature_id))
            }
        }

        # Proteomics: if a protein-ID GMT (e.g. KAE…) is configured, use raw KAE
        # space directly; otherwise translate KAE→GL50803_ to match a gene GMT.
        use_protein_gmt <- view == "proteomics" &&
            !is.null(proteins_gmt) && file.exists(proteins_gmt)
        if (view == "proteomics" && !use_protein_gmt && !is.null(gpm)) {
            univ_ids <- unique(.bucket_protein_to_gene(univ_ids, gpm))
        }
        # Transcriptomics: if MOFA universe came in as GENE_xxx (no original_name),
        # translate via id-map to GL50803_
        if (view == "transcriptomics" && length(univ_ids) > 0 &&
            any(grepl("^GENE_", univ_ids))) {
            id_map <- .bucket_feature_id_map("transcriptomics", mo_root)
            tr <- univ_ids
            hit <- tr %in% names(id_map)
            tr[hit] <- id_map[tr[hit]]
            univ_ids <- tr
        }

        # Translate bucket IDs to the same external space as universe / GMT
        translate <- function(ids) {
            if (view == "proteomics" && !use_protein_gmt && !is.null(gpm)) {
                ids <- .bucket_protein_to_gene(ids, gpm)
            }
            ids
        }

        # Pick gene-set GMT(s) per view
        gmt_paths <- list()
        if (view == "transcriptomics") {
            if (!is.null(kegg_gmt) && file.exists(kegg_gmt)) gmt_paths$kegg <- kegg_gmt
            if (!is.null(go_gmt)   && file.exists(go_gmt))   gmt_paths$go   <- go_gmt
        } else if (view == "proteomics") {
            if (use_protein_gmt) {
                gmt_paths$kegg <- proteins_gmt
            } else if (!is.null(kegg_gmt) && file.exists(kegg_gmt)) {
                gmt_paths$kegg <- kegg_gmt
            }
            if (!is.null(go_gmt) && file.exists(go_gmt)) gmt_paths$go <- go_gmt
        } else if (view == "metabolomics") {
            metab_gmt <- config$modes$multiomics$enrichment$metab_gmt_file %||%
                         config$modes$multiomics$enrichment$metabolomics_gmt_file
            if (!is.null(metab_gmt) && file.exists(metab_gmt)) gmt_paths$kegg <- metab_gmt
        }

        if (length(gmt_paths) == 0) {
            message("    No GMT configured for ", view, "; skipping enrichment.")
            next
        }

        for (set_name in names(gmt_paths)) {
            gmt <- .bucket_load_gmt(gmt_paths[[set_name]])
            if (is.null(gmt) || length(gmt) == 0) next
            for (bk in c("Shared", "MOFA_only", "DIABLO_only")) {
                ids <- bt$ext_id[bt$bucket == bk]
                if (length(ids) < 2) next
                ids <- translate(ids)
                ora <- .bucket_ora(ids, univ_ids, gmt)
                if (is.null(ora)) next
                ora$view   <- view
                ora$bucket <- bk
                ora$set    <- set_name
                f <- file.path(out_dir,
                               sprintf("bucket_enrichment_%s_%s_%s.csv",
                                       view, bk, set_name))
                utils::write.csv(ora, f, row.names = FALSE)
                enrich_tables[[set_name]][[view]][[bk]] <- ora
            }
        }
    }

    # --- Cross-omics pathway convergence within each bucket ---
    # Always include every view we processed as a column in the convergence
    # tables/heatmaps, even if it had zero ORA hits. NA = no hits in that view.
    all_views <- views
    convergence <- list()
    for (set_name in names(enrich_tables)) {
        for (bk in c("Shared", "MOFA_only", "DIABLO_only")) {
            rows <- list()
            for (view in all_views) {
                df <- enrich_tables[[set_name]][[view]][[bk]]
                if (is.null(df)) {
                    # Insert a sentinel row so the wide pivot creates a column
                    # for this view even if it has no ORA hits.
                    rows[[length(rows) + 1L]] <- data.frame(
                        pathway_id = "__none__",
                        pathway_name = "__none__",
                        set_size = NA_integer_, hits_in_set = NA_integer_,
                        n_hits = NA_integer_, n_universe = NA_integer_,
                        pvalue = NA_real_, gene_ratio = NA_character_,
                        bg_ratio = NA_character_, leading_ids = NA_character_,
                        padj = NA_real_, view = view, bucket = bk,
                        set = set_name, stringsAsFactors = FALSE)
                } else {
                    rows[[length(rows) + 1L]] <- df
                }
            }
            if (length(rows) == 0) next
            all_df <- do.call(rbind, rows)
            # Drop the sentinel rows now that the column structure is preserved
            all_df <- all_df[all_df$pathway_id != "__none__" |
                             is.na(all_df$pathway_id), , drop = FALSE]
            if (nrow(all_df) == 0) next
            # Merge cross-omics by pathway_name (pathway IDs differ between
            # species KEGG codes, e.g. gla* vs hsa*, but names are stable).
            all_df$pathway_name_norm <- tolower(trimws(all_df$pathway_name))
            wide <- stats::reshape(
                all_df[, c("pathway_name_norm", "pathway_name", "view",
                           "pvalue", "padj")],
                idvar = c("pathway_name_norm", "pathway_name"),
                timevar = "view", direction = "wide"
            )
            # Collapse duplicate pathway_name rows (different IDs, same name)
            stat_cols <- grep("^(pvalue|padj)\\.", colnames(wide), value = TRUE)
            agg <- stats::aggregate(
                wide[, stat_cols, drop = FALSE],
                by = list(pathway_name_norm = wide$pathway_name_norm),
                FUN = function(x) suppressWarnings(min(x, na.rm = TRUE))
            )
            for (cc in stat_cols) agg[[cc]][!is.finite(agg[[cc]])] <- NA_real_
            name_lookup <- tapply(wide$pathway_name, wide$pathway_name_norm,
                                  function(x) x[1])
            agg$pathway_name <- name_lookup[agg$pathway_name_norm]
            wide <- agg
            pval_cols <- grep("^pvalue\\.", colnames(wide), value = TRUE)
            padj_cols <- grep("^padj\\.",   colnames(wide), value = TRUE)
            wide$n_omics_present  <- rowSums(!is.na(wide[, pval_cols, drop = FALSE]))
            wide$n_omics_pval05   <- rowSums(wide[, pval_cols, drop = FALSE] < 0.05,
                                             na.rm = TRUE)
            wide$n_omics_padj10   <- rowSums(wide[, padj_cols, drop = FALSE] < 0.10,
                                             na.rm = TRUE)
            wide$min_pvalue <- apply(wide[, pval_cols, drop = FALSE], 1,
                                     function(x) suppressWarnings(min(x, na.rm = TRUE)))
            wide$min_pvalue[!is.finite(wide$min_pvalue)] <- NA_real_
            wide$pathway_name_norm <- NULL
            wide <- wide[, c("pathway_name", pval_cols, padj_cols,
                             "n_omics_present", "n_omics_pval05",
                             "n_omics_padj10", "min_pvalue")]
            wide <- wide[order(-wide$n_omics_padj10, -wide$n_omics_pval05,
                               wide$min_pvalue), ]
            f <- file.path(out_dir,
                           sprintf("cross_omics_pathway_convergence_%s_%s.csv",
                                   bk, set_name))
            utils::write.csv(wide, f, row.names = FALSE)
            convergence[[paste(set_name, bk, sep = "_")]] <- wide

            # Heatmap of top convergent pathways: show all where ≥2 omics
            # present (visible), and ≥1 hits p<0.05 (real signal). Falls back
            # to top 15 by min p-value if nothing strict survives.
            shared_pw <- wide[wide$n_omics_present >= 2 & wide$n_omics_pval05 >= 1,
                              , drop = FALSE]
            if (nrow(shared_pw) < 5) shared_pw <- utils::head(wide, 15)
            # Force every configured view as a column, in a stable order, even
            # when the view contributed no enrichments — the column appears
            # with NAs so the heatmap shows all 3 layers consistently.
            for (v in all_views) {
                col <- paste0("pvalue.", v)
                if (!col %in% colnames(shared_pw)) shared_pw[[col]] <- NA_real_
            }
            view_order_cols <- paste0("pvalue.", all_views)
            if (nrow(shared_pw) >= 1) {
                top <- utils::head(shared_pw, 25)
                mat <- as.matrix(-log10(pmax(top[, view_order_cols, drop = FALSE], 1e-300)))
                # NA → 0 for plotting (no hit ⇒ contributes 0 to -log10p heatmap)
                mat[is.na(mat)] <- 0
                colnames(mat) <- tools::toTitleCase(all_views)
                rownames(mat) <- substr(top$pathway_name, 1, 60)
                # ggplot2 + scale_fill_gradient: clean labels, no clipping
                df_long <- data.frame(
                    pathway = factor(rep(rownames(mat), times = ncol(mat)),
                                     levels = rev(rownames(mat))),
                    omics   = factor(rep(colnames(mat), each = nrow(mat)),
                                     levels = colnames(mat)),
                    score   = as.numeric(mat),
                    stringsAsFactors = FALSE)
                p <- ggplot2::ggplot(df_long,
                        ggplot2::aes(x = omics, y = pathway, fill = score)) +
                    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
                    ggplot2::scale_fill_gradientn(
                        colours = c("#FFFFFF", "#FED976", "#E31A1C"),
                        limits = c(0, max(mat, na.rm = TRUE) + 0.001),
                        name = "-log10(p)") +
                    ggplot2::labs(
                        title = sprintf("Cross-omics convergence — %s (%s)",
                                        gsub("_", "-", bk), toupper(set_name)),
                        x = NULL, y = NULL) +
                    ggplot2::theme_minimal(base_size = 11) +
                    ggplot2::theme(
                        axis.text.x = ggplot2::element_text(
                            angle = 30, hjust = 1, size = 11),
                        axis.text.y = ggplot2::element_text(size = 9),
                        plot.title  = ggplot2::element_text(face = "bold"),
                        panel.grid  = ggplot2::element_blank(),
                        plot.margin = ggplot2::margin(8, 8, 14, 8))
                ggplot2::ggsave(
                    filename = file.path(plots_dir,
                                         sprintf("cross_omics_convergence_%s_%s.png",
                                                 bk, set_name)),
                    plot = p,
                    width = 8, height = max(4, 0.32 * nrow(top) + 1.5),
                    dpi = 150)
            }
        }
    }

    message("Cross-method bucket analysis complete; outputs in: ", out_dir)
    invisible(list(
        buckets       = bucket_tables,
        enrichments   = enrich_tables,
        convergence   = convergence,
        factor_used   = factor_pick,
        component_used = component_pick,
        top_n         = top_n
    ))
}
