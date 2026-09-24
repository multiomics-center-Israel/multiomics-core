# =============================================================================
# Enzyme-metabolite pairs: assembling the table
# =============================================================================
#
# One row per enzyme-metabolite pair: a protein that changed in the proteomics
# layer, a metabolite the metabolomics layer measured, and the KEGG annotation
# that says the enzyme class can act on that compound. It is a lookup, not a
# test -- nothing here is scored, and the report says so.
#
# KEGG access and equation parsing live in 07f_enzyme_metabolite.R.


#' Which features of one DE table count as hits
#'
#' The standardized DE tables carry no pass flag, so the rule is the one the
#' ORA path already applies: the adjusted p-value, falling back to the raw one
#' for a contrast where nothing clears adjustment. The rule that was used is
#' returned beside the flags rather than left to be guessed from the numbers.
#'
#' @param tab Standardized DE table with \code{pvalue} and \code{padj}.
#' @param cutoff Significance cutoff.
#' @return List with \code{hit} (logical, never NA) and \code{rule} ("padj" or
#'   "raw p").
#' @examples
#' de_hit_flags(data.frame(pvalue = c(0.001, 0.4), padj = c(0.2, 0.9)))$rule
de_hit_flags <- function(tab, cutoff = 0.05) {
    padj <- suppressWarnings(as.numeric(tab$padj))
    pval <- suppressWarnings(as.numeric(tab$pvalue))
    hit <- !is.na(padj) & padj < cutoff
    if (any(hit)) return(list(hit = hit, rule = "padj"))
    list(hit = !is.na(pval) & pval < cutoff, rule = "raw p")
}


#' Gene symbol for each proteomics feature id
#'
#' The DE tables key on whatever \code{extract_proteomics_de_tables()} resolved
#' the protein to, which is the expression matrix's row name;
#' \code{extract_protein_symbols()} reads symbols off \code{row_data} in the
#' same order. Pairing them by position is only valid because the two are one
#' per protein, so the lengths are checked rather than assumed.
#'
#' @param harmonization_res Harmonization result.
#' @param config Full config object.
#' @return Named character vector, feature id to symbol, or NULL.
#' @keywords internal
.proteomics_symbol_lookup <- function(harmonization_res, config) {
    pre <- harmonization_res$inputs$proteomics
    if (is.null(pre) || is.null(pre$row_data) || is.null(pre$expr_work)) return(NULL)
    symbols <- extract_protein_symbols(pre$row_data, config$modes$proteomics)
    ids <- rownames(pre$expr_work)
    if (is.null(symbols) || is.null(ids) || length(symbols) != length(ids)) return(NULL)
    stats::setNames(as.character(symbols), ids)
}


#' Pair the two layers' contrasts by their canonical key
#'
#' Proteomics and metabolomics spell one comparison differently often enough
#' that matching on the literal name loses every pair.
#'
#' @param prot_names,metab_names Contrast names of each layer.
#' @return Data frame with \code{prot} and \code{metab}, one row per matched
#'   contrast; zero rows when none match.
#' @keywords internal
.match_contrasts <- function(prot_names, metab_names) {
    if (length(prot_names) == 0 || length(metab_names) == 0) {
        return(data.frame(prot = character(0), metab = character(0),
                          stringsAsFactors = FALSE))
    }
    metab_keys <- normalize_contrast_key(metab_names)
    rows <- list()
    for (p in prot_names) {
        hit <- metab_names[metab_keys == normalize_contrast_key(p)]
        if (length(hit) == 0) next
        rows[[length(rows) + 1]] <- data.frame(prot = p, metab = hit[1],
                                               stringsAsFactors = FALSE)
    }
    if (length(rows) == 0) {
        return(data.frame(prot = character(0), metab = character(0),
                          stringsAsFactors = FALSE))
    }
    do.call(rbind, rows)
}


#' Collapse a character vector into one deterministic cell
#'
#' @param x Character vector.
#' @return Sorted unique values joined by ";", or NA when there are none.
#' @keywords internal
.join_unique <- function(x) {
    x <- sort(unique(x[!is.na(x) & nzchar(x)]))
    if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}


#' Enzyme-metabolite pairs from the two layers and KEGG
#'
#' Builds, for each contrast the two layers share, the pairs of a proteomics
#' feature and a measured metabolite that KEGG links through an enzyme class:
#' protein to Entrez to KEGG gene to EC number, EC to reaction to compound.
#' A pair survives when the metabolite was measured, is not a currency
#' compound, and shares a KEGG pathway with the enzyme's gene -- each of the
#' three switchable in the config.
#'
#' One row per (contrast, protein, EC, compound). A pair that KEGG links
#' through several reactions or pathways keeps them in one cell, joined by ";",
#' so a row stays one enzyme-metabolite pair and the table can be counted.
#'
#' Returns NULL, with a message, whenever the annotation cannot be built:
#' no KEGG code for the organism, no OrgDb, a layer missing, no shared
#' contrast, or KEGG unreachable. None of those is an error -- the pairs are a
#' lookup beside the results, and the run does not depend on them.
#'
#' @param de_results Named list of DE results per omics layer.
#' @param harmonization_res Harmonization result.
#' @param config Full config object.
#' @param out_dir Directory whose \code{metabolomics/} subdirectory holds the
#'   compound-pathway cache; KEGG link tables are cached beside it.
#' @param cutoff Significance cutoff for the hit flags.
#' @return Data frame of pairs, or NULL.
build_enzyme_metabolite_pairs <- function(de_results, harmonization_res, config,
                                          out_dir, cutoff = 0.05) {
    cfg <- enzyme_metabolite_config(config)
    if (!cfg$enabled) {
        message("  Enzyme-metabolite pairs disabled in config")
        return(NULL)
    }
    if (is.null(de_results$proteomics) || is.null(de_results$metabolomics)) {
        message("  Enzyme-metabolite pairs need both a proteomics and a ",
                "metabolomics layer")
        return(NULL)
    }

    organism <- config$global$organism
    kegg_org <- resolve_kegg_org_code(organism)
    org_db <- get_organism_db(organism)
    if (is.null(kegg_org) || is.null(org_db)) {
        message("  Enzyme-metabolite pairs need a KEGG organism code and an ",
                "annotation database; skipping for ", organism %||% "no organism")
        return(NULL)
    }

    prot_de <- extract_de_tables(de_results$proteomics, "proteomics", harmonization_res)
    metab_de <- extract_de_tables(de_results$metabolomics, "metabolomics", harmonization_res)
    contrasts <- .match_contrasts(names(prot_de), names(metab_de))
    if (nrow(contrasts) == 0) {
        message("  No contrast is shared by the proteomics and metabolomics layers")
        return(NULL)
    }

    cache_dir <- file.path(out_dir, "kegg_cache")
    ann <- .enzyme_metabolite_annotation(prot_de, metab_de, harmonization_res,
                                         kegg_org, org_db, out_dir, cache_dir)
    if (is.null(ann)) return(NULL)

    symbols <- .proteomics_symbol_lookup(harmonization_res, config)
    rows <- list()
    for (i in seq_len(nrow(contrasts))) {
        rows[[i]] <- .pairs_for_contrast(
            prot_tab = prot_de[[contrasts$prot[i]]],
            metab_tab = metab_de[[contrasts$metab[i]]],
            contrast = contrasts$prot[i], ann = ann, cfg = cfg,
            symbols = symbols, cutoff = cutoff)
    }
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0) {
        message("  No enzyme-metabolite pair survived the filters")
        return(NULL)
    }
    pairs <- do.call(rbind, rows)

    roles <- fetch_reaction_roles(unlist(strsplit(pairs$reaction, ";")), cache_dir)
    pairs$role <- .roles_for_pairs(pairs, roles)

    front <- c("contrast", "protein", "gene_symbol", "ec", "enzyme_log2fc",
               "enzyme_p", "enzyme_padj", "enzyme_hit", "enzyme_hit_rule",
               "compound", "metabolite", "metabolite_log2fc", "metabolite_p",
               "metabolite_padj", "metabolite_hit", "metabolite_hit_rule",
               "role", "reaction", "pathway_id", "pathway_name", "same_direction",
               "note")
    pairs <- pairs[, intersect(front, names(pairs)), drop = FALSE]
    # Pairs whose metabolite cleared FDR lead, then the rest of the pairs, then
    # the enzymes nothing could be paired with; enzyme FDR orders within each.
    tier <- ifelse(is.na(pairs$compound), 3L,
                   ifelse(!is.na(pairs$metabolite_padj) & pairs$metabolite_padj < 0.05,
                          1L, 2L))
    pairs <- pairs[order(pairs$contrast, tier, pairs$enzyme_padj,
                         pairs$metabolite_padj, pairs$gene_symbol, pairs$protein,
                         pairs$ec, pairs$compound, na.last = TRUE), , drop = FALSE]
    rownames(pairs) <- NULL
    n_pairs <- sum(!is.na(pairs$compound))
    message("  Enzyme-metabolite table: ", n_pairs, " pair(s) and ",
            nrow(pairs) - n_pairs, " enzyme(s) with nothing measured to pair, over ",
            nrow(contrasts), " contrast(s)")
    pairs
}


#' The annotation every contrast reuses
#'
#' Fetched once: it does not depend on which contrast is being assembled, and
#' each piece is a bulk KEGG call or a cache read.
#'
#' @param prot_de,metab_de Standardized DE tables per layer.
#' @param harmonization_res Harmonization result.
#' @param kegg_org KEGG organism code.
#' @param org_db OrgDb name for the organism.
#' @param out_dir Cross-enrichment output directory.
#' @param cache_dir Directory for the KEGG link caches.
#' @return List with \code{protein_ec} (protein to EC), \code{ec_compound}
#'   (EC to compound with reaction), \code{gene_pathways} (EC to normalized
#'   pathway id), \code{compound_pathways} and \code{pathway_names}; NULL when
#'   a required piece is missing.
#' @keywords internal
.enzyme_metabolite_annotation <- function(prot_de, metab_de, harmonization_res,
                                          kegg_org, org_db, out_dir, cache_dir) {
    id_map <- tryCatch(
        map_feature_ids_to_entrez(prot_de, "proteomics", harmonization_res, org_db),
        error = function(e) NULL)
    if (is.null(id_map) || nrow(id_map) == 0) {
        message("  No proteomics feature mapped to a gene id")
        return(NULL)
    }
    kegg_genes <- tryCatch(
        convert_entrez_to_kegg(unique(id_map$ENTREZID), kegg_org),
        error = function(e) NULL)
    if (is.null(kegg_genes) || length(kegg_genes) == 0) {
        message("  No gene id converted to a KEGG gene")
        return(NULL)
    }
    id_map$kegg_gene <- unname(kegg_genes[id_map$ENTREZID])
    id_map <- id_map[!is.na(id_map$kegg_gene), , drop = FALSE]

    gene_ec <- kegg_link_table("enzyme", kegg_org, cache_dir)
    ec_rn <- kegg_link_table("reaction", "enzyme", cache_dir)
    rn_cpd <- kegg_link_table("compound", "reaction", cache_dir)
    gene_path <- kegg_link_table("pathway", kegg_org, cache_dir)
    if (is.null(gene_ec) || is.null(ec_rn) || is.null(rn_cpd)) {
        message("  KEGG enzyme annotation unavailable; skipping ",
                "enzyme-metabolite pairs")
        return(NULL)
    }

    # Joined on the bare gene id: the two sides prefix it differently, see
    # kegg_gene_key().
    id_map$gene_key <- kegg_gene_key(id_map$kegg_gene, kegg_org)
    gene_ec$gene_key <- kegg_gene_key(gene_ec$from, kegg_org)
    protein_ec <- merge(id_map[, c("feature_id", "kegg_gene", "gene_key")],
                        gene_ec[, c("gene_key", "to")], by = "gene_key")
    names(protein_ec)[names(protein_ec) == "to"] <- "ec"
    if (nrow(protein_ec) == 0) {
        message("  No measured protein carries an EC number")
        return(NULL)
    }

    # Renamed before merging: both link tables call their columns from/to, and
    # merging them on a mismatched pair leaves two columns named "to".
    ec_rn <- data.frame(ec = ec_rn$from, reaction = ec_rn$to,
                        stringsAsFactors = FALSE)
    ec_rn <- ec_rn[ec_rn$ec %in% protein_ec$ec, , drop = FALSE]
    rn_cpd <- data.frame(reaction = rn_cpd$from, compound = rn_cpd$to,
                         stringsAsFactors = FALSE)
    ec_compound <- merge(ec_rn, rn_cpd, by = "reaction")

    metab_map <- tryCatch(map_metabolite_ids_to_kegg(metab_de, harmonization_res),
                          error = function(e) NULL)
    if (is.null(metab_map) || nrow(metab_map) == 0) {
        message("  No metabolite mapped to a KEGG compound")
        return(NULL)
    }

    cpd_path <- get_kegg_compound_pathways(file.path(out_dir, "metabolomics"))
    gene_pathways <- NULL
    if (!is.null(gene_path) && nrow(gene_path) > 0) {
        gene_path$gene_key <- kegg_gene_key(gene_path$from, kegg_org)
        gp <- merge(protein_ec[, c("gene_key", "ec")],
                    gene_path[, c("gene_key", "to")], by = "gene_key")
        gene_pathways <- unique(data.frame(
            ec = gp$ec, pathway = normalize_pathway_join_key(gp$to, kegg_org),
            stringsAsFactors = FALSE))
    }

    list(
        protein_ec = unique(protein_ec[, c("feature_id", "ec")]),
        ec_compound = unique(ec_compound),
        metab_map = metab_map,
        gene_pathways = gene_pathways,
        compound_pathways = if (is.null(cpd_path)) NULL else unique(data.frame(
            compound = cpd_path$compound,
            pathway = normalize_pathway_join_key(cpd_path$pathway, kegg_org),
            stringsAsFactors = FALSE)),
        pathway_names = if (is.null(cpd_path)) NULL else stats::setNames(
            cpd_path$name, normalize_pathway_join_key(cpd_path$pathway, kegg_org))
    )
}


#' Pairs for one matched contrast
#'
#' @param prot_tab,metab_tab Standardized DE tables of the two layers.
#' @param contrast Contrast name to record, as the proteomics layer spells it.
#' @param ann Output of \code{.enzyme_metabolite_annotation()}.
#' @param cfg Output of \code{enzyme_metabolite_config()}.
#' @param symbols Named feature-id-to-symbol vector, or NULL.
#' @param cutoff Significance cutoff.
#' @return Data frame of pairs without the \code{role} column, or NULL.
#' @keywords internal
.pairs_for_contrast <- function(prot_tab, metab_tab, contrast, ann, cfg,
                                symbols, cutoff) {
    prot_hits <- de_hit_flags(prot_tab, cutoff)
    prot_tab$hit <- prot_hits$hit
    metab_hits <- de_hit_flags(metab_tab, cutoff)
    metab_tab$hit <- metab_hits$hit

    enzymes <- if (cfg$enzyme_hits_only) {
        prot_tab[prot_tab$hit, , drop = FALSE]
    } else prot_tab
    if (nrow(enzymes) == 0) return(NULL)

    with_ec <- merge(enzymes, ann$protein_ec, by = "feature_id")
    pairs <- if (nrow(with_ec) == 0) with_ec else merge(with_ec, ann$ec_compound, by = "ec")

    if (cfg$drop_currency_metabolites) {
        pairs <- pairs[!pairs$compound %in% ENZYME_METABOLITE_CURRENCY_COMPOUNDS, ,
                       drop = FALSE]
    }

    measured <- merge(ann$metab_map, metab_tab, by = "feature_id")
    if (nrow(pairs) > 0) {
        pairs <- merge(pairs, measured, by.x = "compound", by.y = "KEGG_CPD",
                       suffixes = c("_enzyme", "_metabolite"))
    }
    if (nrow(pairs) > 0) {
        pairs$pathway <- .shared_pathways(pairs$ec, pairs$compound, ann)
        if (cfg$require_shared_pathway) {
            pairs <- pairs[!is.na(pairs$pathway), , drop = FALSE]
        }
    }

    key <- paste(pairs$feature_id_enzyme, pairs$ec, pairs$compound, sep = "\r")
    collapsed <- lapply(split(seq_len(nrow(pairs)), key), function(idx) {
        r <- pairs[idx[1], , drop = FALSE]
        data.frame(
            contrast = contrast,
            protein = r$feature_id_enzyme,
            gene_symbol = if (is.null(symbols)) NA_character_
                          else unname(symbols[r$feature_id_enzyme]),
            ec = r$ec,
            enzyme_log2fc = r$log2fc_enzyme, enzyme_p = r$pvalue_enzyme,
            enzyme_padj = r$padj_enzyme, enzyme_hit = r$hit_enzyme,
            enzyme_hit_rule = prot_hits$rule,
            compound = r$compound,
            metabolite = r$feature_id_metabolite,
            metabolite_log2fc = r$log2fc_metabolite,
            metabolite_p = r$pvalue_metabolite,
            metabolite_padj = r$padj_metabolite,
            metabolite_hit = r$hit_metabolite,
            metabolite_hit_rule = metab_hits$rule,
            reaction = .join_unique(pairs$reaction[idx]),
            pathway_id = r$pathway,
            pathway_name = .pathway_labels(r$pathway, ann$pathway_names),
            same_direction = !is.na(r$log2fc_enzyme) & !is.na(r$log2fc_metabolite) &
                sign(r$log2fc_enzyme) == sign(r$log2fc_metabolite),
            stringsAsFactors = FALSE)
    })
    paired <- if (nrow(pairs) == 0) NULL else do.call(rbind, collapsed)
    if (!is.null(paired)) paired$note <- NA_character_
    if (!isTRUE(cfg$list_unpaired_enzymes)) return(paired)

    # Every enzyme that changed is listed, even when nothing it acts on was
    # measured: a reader asking "what happened to this enzyme's metabolites?"
    # needs to see that the question has no answer here, not an absent row.
    left <- enzymes[!enzymes$feature_id %in% paired$protein, , drop = FALSE]
    if (nrow(left) == 0) return(paired)
    ec_of <- split(with_ec$ec, with_ec$feature_id)
    unpaired <- lapply(seq_len(nrow(left)), function(i) {
        id <- left$feature_id[i]
        ecs <- .join_unique(ec_of[[id]])
        data.frame(
            contrast = contrast, protein = id,
            gene_symbol = if (is.null(symbols)) NA_character_ else unname(symbols[id]),
            ec = ecs,
            enzyme_log2fc = left$log2fc[i], enzyme_p = left$pvalue[i],
            enzyme_padj = left$padj[i], enzyme_hit = left$hit[i],
            enzyme_hit_rule = prot_hits$rule,
            compound = NA_character_, metabolite = NA_character_,
            metabolite_log2fc = NA_real_, metabolite_p = NA_real_,
            metabolite_padj = NA_real_, metabolite_hit = NA,
            metabolite_hit_rule = metab_hits$rule,
            reaction = NA_character_, pathway_id = NA_character_,
            pathway_name = NA_character_, same_direction = NA,
            note = if (is.na(ecs)) "no EC number in KEGG"
                   else if (cfg$require_shared_pathway)
                       "no measured metabolite in a shared pathway"
                   else "no measured metabolite",
            stringsAsFactors = FALSE)
    })
    rbind(paired, do.call(rbind, unpaired))
}


#' Pathways an EC and a compound both sit in
#'
#' @param ec,compound Vectors of equal length.
#' @param ann Output of \code{.enzyme_metabolite_annotation()}.
#' @return Character vector of ";"-joined pathway ids, NA where the two share
#'   none or a membership table is missing.
#' @keywords internal
.shared_pathways <- function(ec, compound, ann) {
    if (is.null(ann$gene_pathways) || is.null(ann$compound_pathways)) {
        return(rep(NA_character_, length(ec)))
    }
    by_ec <- split(ann$gene_pathways$pathway, ann$gene_pathways$ec)
    by_cpd <- split(ann$compound_pathways$pathway, ann$compound_pathways$compound)
    vapply(seq_along(ec), function(i) {
        .join_unique(intersect(by_ec[[ec[i]]] %||% character(0),
                               by_cpd[[compound[i]]] %||% character(0)))
    }, character(1))
}


#' Readable names for ";"-joined pathway ids
#'
#' @param ids Character vector of joined pathway ids, possibly NA.
#' @param names_map Named vector of pathway names, or NULL.
#' @return Character vector of joined names, NA where nothing is known.
#' @keywords internal
.pathway_labels <- function(ids, names_map) {
    if (is.null(names_map)) return(rep(NA_character_, length(ids)))
    vapply(ids, function(id) {
        if (is.na(id)) return(NA_character_)
        .join_unique(unname(names_map[strsplit(id, ";", fixed = TRUE)[[1]]]))
    }, character(1), USE.NAMES = FALSE)
}


#' The role each pair's compound plays in the reactions behind it
#'
#' @param pairs Pair table with \code{reaction} and \code{compound}.
#' @param roles Output of \code{fetch_reaction_roles()}, or NULL.
#' @return Character vector: "substrate", "product", "substrate;product" where
#'   the reactions disagree, or "unknown" where no equation was read.
#' @keywords internal
.roles_for_pairs <- function(pairs, roles) {
    if (is.null(roles) || nrow(roles) == 0) return(rep("unknown", nrow(pairs)))
    key <- paste(roles$reaction, roles$compound, sep = "\r")
    by_key <- split(roles$role, key)
    vapply(seq_len(nrow(pairs)), function(i) {
        rns <- strsplit(pairs$reaction[i], ";", fixed = TRUE)[[1]]
        found <- unlist(by_key[paste(rns, pairs$compound[i], sep = "\r")],
                        use.names = FALSE)
        out <- .join_unique(found)
        if (is.na(out)) "unknown" else out
    }, character(1))
}


#' Write the enzyme-metabolite pairs, and clear a previous run's file
#'
#' The report includes the table on \code{file.exists()} alone, so a rerun that
#' produces no pairs must leave no file behind.
#'
#' @param pairs Output of \code{build_enzyme_metabolite_pairs()}, or NULL.
#' @param out_dir Directory to write into.
#' @return The path written, or NULL.
write_enzyme_metabolite_pairs <- function(pairs, out_dir) {
    path <- file.path(out_dir, "enzyme_metabolite_pairs.tsv")
    if (is.null(pairs) || !is.data.frame(pairs) || nrow(pairs) == 0) {
        if (file.exists(path)) unlink(path)
        return(NULL)
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.table(pairs, path, sep = "\t", quote = FALSE, row.names = FALSE,
                       na = "")
    path
}
