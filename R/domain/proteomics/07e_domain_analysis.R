# =============================================================================
# Protein Domain & Functional Feature Analysis
# =============================================================================
# Queries InterPro/Pfam domain annotations and UniProt functional features
# (binding sites, active sites, PTMs) for significant proteins.

# -----------------------------------------------------------------------------
# Main Orchestrator
# -----------------------------------------------------------------------------

#' Run protein domain and functional feature analysis
#'
#' @param de_res   DE results list (with summary_df)
#' @param config   Full config list
#' @param out_dir  Output root directory
#' @return List with domain enrichment and feature annotation results
#' @export
run_domain_analysis <- function(de_res, config, out_dir) {
    message("=== Protein Domain & Functional Feature Analysis ===")

    cfg <- config$modes$proteomics
    domain_cfg <- cfg$domain_analysis %||% list()

    summary_df <- de_res$summary_df
    if (is.null(summary_df)) {
        message("No summary_df available. Skipping domain analysis.")
        return(NULL)
    }

    # Extract significant proteins across all contrasts
    padj_cols <- grep("^padj\\.imputs\\.", colnames(summary_df), value = TRUE)
    fc_cols   <- grep("^linearFC\\.imputs\\.", colnames(summary_df), value = TRUE)

    de_cutoff <- cfg$de$p_cutoff %||% 0.05
    fc_cutoff <- cfg$de$linear_fc_cutoff %||% 1.5

    sig_mask <- rep(FALSE, nrow(summary_df))
    for (i in seq_along(padj_cols)) {
        padj_vals <- as.numeric(summary_df[[padj_cols[i]]])
        fc_vals   <- as.numeric(summary_df[[fc_cols[min(i, length(fc_cols))]]])
        sig_mask <- sig_mask | (!is.na(padj_vals) & padj_vals < de_cutoff &
                                !is.na(fc_vals) & abs(fc_vals) >= fc_cutoff)
    }

    # summary_df is keyed by the DE FeatureID column, not the raw protein_id
    # column (cfg$id_columns$protein_id, e.g. "Protein.Group"). Using the latter
    # silently looked up a non-existent column -> 0 significant -> always skipped.
    id_col <- cfg$de_table$id_col %||% "FeatureID"
    if (!id_col %in% colnames(summary_df)) {
        id_col <- intersect(c("FeatureID", cfg$id_columns$protein_id), colnames(summary_df))[1]
    }
    if (is.na(id_col) || !nzchar(id_col) || !id_col %in% colnames(summary_df)) {
        stop("Domain analysis: could not determine ID column in summary_df. Available: ",
             paste(colnames(summary_df), collapse = ", "))
    }
    sig_proteins <- summary_df[[id_col]][sig_mask]

    # Genes may be "SYMBOL | description" (custom-mapped) and/or ";"-joined;
    # extract the bare leading gene symbol for UniProt/InterPro queries.
    .first_gene_symbol <- function(g) {
        if (is.na(g) || !nzchar(g)) return(NA_character_)
        trimws(strsplit(strsplit(g, ";", fixed = TRUE)[[1]][1], "|", fixed = TRUE)[[1]][1])
    }

    # Map to gene symbols if available
    if ("Genes" %in% colnames(summary_df)) {
        gene_map <- setNames(
            vapply(summary_df[["Genes"]], .first_gene_symbol, character(1)),
            summary_df[[id_col]]
        )
        sig_genes <- gene_map[sig_proteins]
        sig_genes <- sig_genes[!is.na(sig_genes) & nzchar(sig_genes)]
    } else {
        sig_genes <- sig_proteins
    }

    all_genes <- if ("Genes" %in% colnames(summary_df)) {
        g <- vapply(summary_df[["Genes"]], .first_gene_symbol, character(1))
        g[!is.na(g) & nzchar(g)]
    } else {
        summary_df[[id_col]]
    }

    message("  Significant proteins: ", length(sig_genes),
            " / Background: ", length(all_genes))

    if (length(sig_genes) < 5) {
        message("Too few significant proteins for domain analysis. Skipping.")
        return(NULL)
    }

    domain_dir <- file.path(out_dir, "domain_analysis")
    dir.create(domain_dir, recursive = TRUE, showWarnings = FALSE)

    organism <- cfg$annotation$organism %||% "Homo sapiens"

    # 1. InterPro/Pfam domain enrichment
    domain_enrichment <- run_domain_enrichment(
        sig_genes, all_genes, organism, domain_dir
    )

    # 2. UniProt functional features (binding sites, active sites, PTMs)
    uniprot_features <- query_uniprot_features(
        sig_genes, organism, domain_dir
    )

    # 3. Generate domain plots
    figures <- generate_domain_plots(domain_enrichment, uniprot_features, domain_dir)

    results <- list(
        domain_enrichment = domain_enrichment,
        uniprot_features  = uniprot_features,
        figures           = figures,
        n_sig_proteins    = length(sig_genes)
    )

    saveRDS(results, file.path(domain_dir, "domain_analysis_results.rds"))
    message("Domain analysis complete!")
    results
}

# -----------------------------------------------------------------------------
# InterPro / Pfam Domain Enrichment
# -----------------------------------------------------------------------------

#' Test for enrichment of protein domains among significant proteins
#'
#' Uses biomaRt to fetch InterPro/Pfam annotations, then runs Fisher's exact
#' test for over-representation of each domain in significant vs background.
#'
#' @param sig_genes Character vector of significant gene symbols
#' @param all_genes Character vector of all (background) gene symbols
#' @param organism Organism name
#' @param output_dir Output directory
#' @return Data frame of enriched domains
run_domain_enrichment <- function(sig_genes, all_genes, organism, output_dir) {
    message("Running InterPro/Pfam domain enrichment...")

    domain_annotations <- NULL

    # Strategy 1: biomaRt
    if (requireNamespace("biomaRt", quietly = TRUE)) {
        domain_annotations <- tryCatch({
            dataset <- .domain_ensembl_dataset(organism)
            if (is.null(dataset)) stop("Unknown Ensembl dataset for: ", organism)

            mart <- biomaRt::useMart("ensembl", dataset = dataset)

            annotations <- biomaRt::getBM(
                attributes = c("external_gene_name", "interpro", "interpro_description",
                               "pfam", "pfam_start", "pfam_end"),
                filters = "external_gene_name",
                values = unique(all_genes),
                mart = mart
            )

            if (nrow(annotations) == 0) stop("No domain annotations returned")
            annotations
        }, error = function(e) {
            message("  biomaRt domain lookup failed: ", e$message)
            NULL
        })
    }

    # Strategy 2: UniProt REST API fallback for InterPro domains
    if (is.null(domain_annotations) || nrow(domain_annotations) == 0) {
        domain_annotations <- tryCatch({
            fetch_interpro_from_uniprot(all_genes, organism)
        }, error = function(e) {
            message("  UniProt InterPro fallback failed: ", e$message)
            NULL
        })
    }

    if (is.null(domain_annotations) || nrow(domain_annotations) == 0) {
        message("  No domain annotations available. Skipping domain enrichment.")
        return(NULL)
    }

    # Run enrichment for InterPro domains
    interpro_results <- NULL
    if ("interpro" %in% colnames(domain_annotations)) {
        ipro <- domain_annotations[!is.na(domain_annotations$interpro) &
                                    nzchar(domain_annotations$interpro), ]
        if (nrow(ipro) > 0) {
            interpro_results <- .domain_ora(
                sig_genes, all_genes,
                gene_col = "external_gene_name",
                domain_col = "interpro",
                desc_col = "interpro_description",
                annotations = ipro,
                label = "InterPro"
            )
        }
    }

    # Run enrichment for Pfam domains
    pfam_results <- NULL
    if ("pfam" %in% colnames(domain_annotations)) {
        pfam <- domain_annotations[!is.na(domain_annotations$pfam) &
                                    nzchar(domain_annotations$pfam), ]
        if (nrow(pfam) > 0) {
            pfam_results <- .domain_ora(
                sig_genes, all_genes,
                gene_col = "external_gene_name",
                domain_col = "pfam",
                desc_col = NULL,
                annotations = pfam,
                label = "Pfam"
            )
        }
    }

    # Combine and save
    combined <- rbind(interpro_results, pfam_results)
    if (!is.null(combined) && nrow(combined) > 0) {
        combined <- combined[order(combined$pvalue), ]
        write.csv(combined, file.path(output_dir, "domain_enrichment.csv"), row.names = FALSE)
        n_sig <- sum(combined$padj < 0.05, na.rm = TRUE)
        message("  Domain enrichment: ", nrow(combined), " domains tested, ",
                n_sig, " significant (padj < 0.05)")
    }

    list(annotations = domain_annotations, enrichment = combined)
}

#' Over-representation analysis for domains
#' @keywords internal
.domain_ora <- function(sig_genes, all_genes, gene_col, domain_col,
                        desc_col, annotations, label) {
    domain_gene_map <- split(annotations[[gene_col]], annotations[[domain_col]])
    domain_gene_map <- domain_gene_map[lengths(domain_gene_map) >= 3]

    if (length(domain_gene_map) == 0) return(NULL)

    results <- lapply(names(domain_gene_map), function(dom_id) {
        dom_genes <- unique(intersect(domain_gene_map[[dom_id]], all_genes))
        sig_in_dom <- length(intersect(sig_genes, dom_genes))
        if (sig_in_dom < 2) return(NULL)

        sig_not_dom <- length(sig_genes) - sig_in_dom
        dom_not_sig <- length(dom_genes) - sig_in_dom
        neither <- length(all_genes) - length(sig_genes) - dom_not_sig

        mat <- matrix(c(sig_in_dom, dom_not_sig, sig_not_dom, max(0, neither)), nrow = 2)
        test <- fisher.test(mat, alternative = "greater")

        desc <- if (!is.null(desc_col)) {
            idx <- which(annotations[[domain_col]] == dom_id)[1]
            annotations[[desc_col]][idx]
        } else dom_id

        data.frame(
            domain_id = dom_id,
            domain_description = desc %||% dom_id,
            database = label,
            domain_size = length(dom_genes),
            overlap = sig_in_dom,
            pvalue = test$p.value,
            odds_ratio = as.numeric(test$estimate),
            genes = paste(intersect(sig_genes, dom_genes), collapse = ";"),
            stringsAsFactors = FALSE
        )
    })

    results_df <- do.call(rbind, results[!sapply(results, is.null)])
    if (!is.null(results_df) && nrow(results_df) > 0) {
        results_df$padj <- p.adjust(results_df$pvalue, method = "BH")
    }
    results_df
}

#' Map organism name to Ensembl dataset
#' @keywords internal
.domain_ensembl_dataset <- function(organism) {
    datasets <- c(
        "Homo sapiens"              = "hsapiens_gene_ensembl",
        "Mus musculus"              = "mmusculus_gene_ensembl",
        "Rattus norvegicus"         = "rnorvegicus_gene_ensembl",
        "Danio rerio"              = "drerio_gene_ensembl",
        "Drosophila melanogaster"  = "dmelanogaster_gene_ensembl",
        "Caenorhabditis elegans"   = "celegans_gene_ensembl",
        "Saccharomyces cerevisiae" = "scerevisiae_gene_ensembl",
        "Sus scrofa"               = "sscrofa_gene_ensembl",
        "Bos taurus"               = "btaurus_gene_ensembl",
        "Gallus gallus"            = "ggallus_gene_ensembl"
    )
    organism_norm <- normalize_organism_name(organism)
    unname(datasets[organism_norm])
}

# -----------------------------------------------------------------------------
# UniProt Functional Features (Binding Sites, Active Sites, PTMs)
# -----------------------------------------------------------------------------

#' Query UniProt REST API for functional features of significant proteins
#'
#' @param genes Character vector of gene symbols
#' @param organism Organism name
#' @param output_dir Output directory
#' @return Data frame of protein features
query_uniprot_features <- function(genes, organism, output_dir) {
    message("Querying UniProt for functional features...")

    if (length(genes) > 200) {
        message("  Limiting to top 200 proteins for UniProt query")
        genes <- genes[1:200]
    }

    taxid <- .uniprot_taxid(organism)
    if (is.null(taxid)) {
        message("  Unknown organism for UniProt. Skipping.")
        return(NULL)
    }

    # Query UniProt REST API for features
    all_features <- list()

    # Batch genes in groups of 25 to avoid URL length limits
    gene_batches <- split(genes, ceiling(seq_along(genes) / 25))

    for (batch in gene_batches) {
        gene_query <- paste(batch, collapse = "+OR+gene_exact:")
        url <- paste0(
            "https://rest.uniprot.org/uniprotkb/search?query=",
            "(gene_exact:", gene_query, ")+AND+(organism_id:", taxid, ")",
            "&fields=accession,gene_names,ft_binding,ft_act_site,ft_site,",
            "ft_domain,ft_region,ft_motif,ft_mod_res",
            "&format=tsv&size=500"
        )

        batch_result <- tryCatch({
            resp <- utils::read.delim(url(url), stringsAsFactors = FALSE,
                                       quote = "", fill = TRUE)
            resp
        }, error = function(e) {
            # Silently skip failed batches
            NULL
        })

        if (!is.null(batch_result) && nrow(batch_result) > 0) {
            all_features <- c(all_features, list(batch_result))
        }

        Sys.sleep(0.5)  # Rate limit
    }

    if (length(all_features) == 0) {
        message("  No UniProt features retrieved.")
        return(NULL)
    }

    features_df <- do.call(rbind, all_features)
    features_df <- unique(features_df)

    # Parse feature columns into a long-format summary
    feature_summary <- .parse_uniprot_features(features_df)

    if (!is.null(feature_summary) && nrow(feature_summary) > 0) {
        write.csv(feature_summary, file.path(output_dir, "uniprot_features.csv"),
                  row.names = FALSE)
        write.csv(features_df, file.path(output_dir, "uniprot_features_raw.csv"),
                  row.names = FALSE)
        message("  Retrieved features for ", nrow(features_df), " proteins")

        # Summarize feature types
        if ("feature_type" %in% colnames(feature_summary)) {
            type_counts <- as.data.frame(table(feature_summary$feature_type),
                                          stringsAsFactors = FALSE)
            colnames(type_counts) <- c("feature_type", "count")
            type_counts <- type_counts[order(-type_counts$count), ]
            write.csv(type_counts, file.path(output_dir, "feature_type_summary.csv"),
                      row.names = FALSE)
            message("  Feature types found: ",
                    paste(head(type_counts$feature_type, 5), collapse = ", "))
        }
    }

    list(raw = features_df, summary = feature_summary)
}

#' Parse UniProt TSV feature columns into long-format summary
#' @keywords internal
.parse_uniprot_features <- function(features_df) {
    if (is.null(features_df) || nrow(features_df) == 0) return(NULL)

    feature_cols <- c(
        "Binding.site"        = "Binding Site",
        "Active.site"         = "Active Site",
        "Site"                = "Functional Site",
        "Domain..FT."         = "Domain",
        "Region"              = "Region",
        "Motif"               = "Motif",
        "Modified.residue"    = "PTM"
    )

    # Identify which columns exist (TSV header names vary)
    available <- intersect(names(feature_cols), colnames(features_df))
    if (length(available) == 0) {
        # Try alternate naming
        alt_map <- c(
            "Binding site" = "Binding Site", "Active site" = "Active Site",
            "Site" = "Functional Site", "Domain [FT]" = "Domain",
            "Region" = "Region", "Motif" = "Motif",
            "Modified residue" = "PTM"
        )
        available <- intersect(names(alt_map), colnames(features_df))
        if (length(available) == 0) return(NULL)
        feature_cols <- alt_map
    }

    gene_col <- if ("Gene.Names" %in% colnames(features_df)) "Gene.Names"
                else if ("Gene Names" %in% colnames(features_df)) "Gene Names"
                else if ("Entry" %in% colnames(features_df)) "Entry"
                else colnames(features_df)[1]

    rows <- list()
    for (fc in available) {
        label <- feature_cols[fc]
        col_data <- features_df[[fc]]
        for (i in seq_along(col_data)) {
            val <- col_data[i]
            if (!is.na(val) && nzchar(trimws(val))) {
                gene <- trimws(strsplit(as.character(features_df[[gene_col]][i]), " ")[[1]][1])
                # Each annotation entry may have multiple features separated by "; "
                entries <- strsplit(val, ";\\s*")[[1]]
                for (entry in entries) {
                    rows <- c(rows, list(data.frame(
                        gene = gene,
                        feature_type = label,
                        feature_detail = trimws(entry),
                        stringsAsFactors = FALSE
                    )))
                }
            }
        }
    }

    if (length(rows) == 0) return(NULL)
    do.call(rbind, rows)
}

#' Map organism to NCBI taxonomy ID for UniProt queries
#' @keywords internal
.uniprot_taxid <- function(organism) {
    taxids <- c(
        "Homo sapiens" = "9606", "Mus musculus" = "10090",
        "Rattus norvegicus" = "10116", "Danio rerio" = "7955",
        "Drosophila melanogaster" = "7227", "Caenorhabditis elegans" = "6239",
        "Saccharomyces cerevisiae" = "559292", "Sus scrofa" = "9823",
        "Bos taurus" = "9913", "Gallus gallus" = "9031",
        "Giardia lamblia" = "5741"
    )
    organism_norm <- normalize_organism_name(organism)
    unname(taxids[organism_norm])
}

#' Fetch InterPro annotations from UniProt as fallback
#' @keywords internal
fetch_interpro_from_uniprot <- function(genes, organism) {
    taxid <- .uniprot_taxid(organism)
    if (is.null(taxid)) return(NULL)

    if (length(genes) > 100) genes <- genes[1:100]

    gene_query <- paste(genes, collapse = "+OR+gene_exact:")
    url <- paste0(
        "https://rest.uniprot.org/uniprotkb/search?query=",
        "(gene_exact:", gene_query, ")+AND+(organism_id:", taxid, ")",
        "&fields=accession,gene_names,xref_interpro,xref_pfam",
        "&format=tsv&size=500"
    )

    resp <- utils::read.delim(url(url), stringsAsFactors = FALSE, quote = "", fill = TRUE)
    if (nrow(resp) == 0) return(NULL)

    # Parse into the expected format
    gene_col <- if ("Gene.Names" %in% colnames(resp)) "Gene.Names"
                else if ("Gene Names" %in% colnames(resp)) "Gene Names"
                else colnames(resp)[1]

    rows <- list()
    for (i in seq_len(nrow(resp))) {
        gene <- trimws(strsplit(as.character(resp[[gene_col]][i]), " ")[[1]][1])

        # InterPro
        ipro_col <- grep("interpro", colnames(resp), ignore.case = TRUE, value = TRUE)[1]
        if (!is.null(ipro_col) && !is.na(resp[[ipro_col]][i]) && nzchar(resp[[ipro_col]][i])) {
            for (ipro in strsplit(resp[[ipro_col]][i], ";\\s*")[[1]]) {
                rows <- c(rows, list(data.frame(
                    external_gene_name = gene, interpro = trimws(ipro),
                    interpro_description = "", pfam = NA_character_,
                    pfam_start = NA, pfam_end = NA, stringsAsFactors = FALSE
                )))
            }
        }

        # Pfam
        pfam_col <- grep("pfam", colnames(resp), ignore.case = TRUE, value = TRUE)[1]
        if (!is.null(pfam_col) && !is.na(resp[[pfam_col]][i]) && nzchar(resp[[pfam_col]][i])) {
            for (pf in strsplit(resp[[pfam_col]][i], ";\\s*")[[1]]) {
                rows <- c(rows, list(data.frame(
                    external_gene_name = gene, interpro = NA_character_,
                    interpro_description = NA_character_, pfam = trimws(pf),
                    pfam_start = NA, pfam_end = NA, stringsAsFactors = FALSE
                )))
            }
        }
    }

    if (length(rows) == 0) return(NULL)
    do.call(rbind, rows)
}

# -----------------------------------------------------------------------------
# Visualizations
# -----------------------------------------------------------------------------

generate_domain_plots <- function(domain_enrichment, uniprot_features, output_dir) {
    plots_dir <- file.path(output_dir, "plots")
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
    figures <- list()

    # Domain enrichment dotplot
    if (!is.null(domain_enrichment$enrichment)) {
        enr <- domain_enrichment$enrichment
        sig_enr <- enr[!is.na(enr$padj) & enr$padj < 0.05, ]
        if (nrow(sig_enr) > 0) {
            sig_enr <- head(sig_enr[order(sig_enr$padj), ], 25)
            sig_enr$neg_log10_padj <- -log10(sig_enr$padj + 1e-300)
            sig_enr$label <- ifelse(nchar(sig_enr$domain_description) > 50,
                                     paste0(substr(sig_enr$domain_description, 1, 47), "..."),
                                     sig_enr$domain_description)

            p <- ggplot2::ggplot(sig_enr,
                     ggplot2::aes(x = neg_log10_padj,
                                  y = reorder(label, neg_log10_padj),
                                  size = overlap, color = database)) +
                ggplot2::geom_point() +
                ggplot2::scale_size_continuous(range = c(2, 8), name = "Overlap") +
                ggplot2::labs(title = "Enriched Protein Domains (DE proteins)",
                              x = "-log10(adjusted p-value)", y = NULL) +
                ggplot2::theme_minimal() +
                ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

            fig_path <- file.path(plots_dir, "domain_enrichment_dotplot.png")
            ggplot2::ggsave(fig_path, p, width = 10, height = max(4, nrow(sig_enr) * 0.3 + 2))
            figures$domain_dotplot <- fig_path
        }
    }

    # Feature type bar chart
    if (!is.null(uniprot_features$summary) && nrow(uniprot_features$summary) > 0) {
        type_counts <- as.data.frame(table(uniprot_features$summary$feature_type),
                                      stringsAsFactors = FALSE)
        colnames(type_counts) <- c("feature_type", "count")
        type_counts <- type_counts[order(-type_counts$count), ]

        p <- ggplot2::ggplot(type_counts,
                 ggplot2::aes(x = reorder(feature_type, count), y = count)) +
            ggplot2::geom_bar(stat = "identity", fill = "#4a90d9") +
            ggplot2::coord_flip() +
            ggplot2::labs(title = "UniProt Functional Features in DE Proteins",
                          x = NULL, y = "Number of annotations") +
            ggplot2::theme_minimal()

        fig_path <- file.path(plots_dir, "feature_type_barchart.png")
        ggplot2::ggsave(fig_path, p, width = 8, height = 5)
        figures$feature_types <- fig_path
    }

    figures
}
