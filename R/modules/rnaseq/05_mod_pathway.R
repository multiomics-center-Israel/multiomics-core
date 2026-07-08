#' Build a pure-compute ORA worker with a minimal captured environment
#'
#' Returns a `function(job)` that computes cluster ORA for one job. Defining it
#' here (not inside mod_rnaseq_pathway) bounds the closure's environment to just
#' the arguments below, so future.apply serializes only these to parallel
#' workers — not unrelated large objects (e.g. the expression matrix in `pre`)
#' that would otherwise be captured from the module frame. The worker does pure
#' computation only: no file I/O, no plotting, no messages.
#'
#' @param gene_lists Output of build_gene_lists().
#' @param local_tables Output of load_local_pathway_tables().
#' @param pval_cutoff,padj_method ORA thresholds passed to run_cluster_ora_compute().
#' @param go_simplify,orgdb GO-simplify controls (gated; default off).
#' @return A function(job) -> list(job, db_type, result), where result is the
#'   4-element list from run_cluster_ora_compute() (or list() if not significant).
#' @noRd
.make_ora_worker <- function(gene_lists, local_tables, pval_cutoff, padj_method,
                             go_simplify, orgdb) {
    force(gene_lists); force(local_tables); force(pval_cutoff)
    force(padj_method); force(go_simplify); force(orgdb)
    function(job) {
        clusters <- gene_lists[[job$clust_method]][[job$clust_round]]
        tbl      <- local_tables[[job$db_name]]
        db_type  <- if (grepl("^GO", job$db_name)) "GO" else "KEGG"
        db_ont   <- if (db_type == "GO") sub("^GO_", "", job$db_name) else NULL

        res <- run_cluster_ora_compute(
            clusters      = clusters,
            TERM2GENE     = tbl$TERM2GENE,
            TERM2NAME     = tbl$TERM2NAME,
            type          = db_type,
            pvalueCutoff  = pval_cutoff,
            pAdjustMethod = padj_method,
            go_simplify   = go_simplify,
            orgdb         = orgdb,
            ont           = db_ont
        )
        list(job = job, db_type = db_type, result = res)
    }
}

#' RNA-seq Pathway Analysis Module
#'
#' Orchestrates organism detection, gene annotation, gene set loading,
#' and pathway/enrichment analysis (fGSEA + ORA).
#'
#' @param de_res  Result from mod_rnaseq_de() — list(dds, tables)
#' @param pre     Preprocessed data list (with expr_filt, meta)
#' @param config  Full pipeline config
#' @param out_dir Output directory for the RNA mode (e.g. .../rna)
#' @param clustering_res Result from mod_rnaseq_clustering(), or NULL. Used only
#'   by the local (offline) enrichment path: when provided (single-omics mode)
#'   cluster-based ORA runs; when NULL (multiomics mode) ORA is skipped with a
#'   warning and GSEA still runs. Ignored by the online fallback.
#' @return List with annotation, pathway_results, and plot_files
#' @export
mod_rnaseq_pathway <- function(de_res, pre, config, out_dir, clustering_res = NULL) {

    rna_cfg <- config$modes$rna
    ann_cfg <- rna_cfg$annotation %||% list()
    pw_cfg  <- rna_cfg$pathway
    enr_cfg <- rna_cfg$enrichment %||% list()

    # Skip entirely if pathway analysis is disabled
    if (is.null(pw_cfg) || isFALSE(pw_cfg$enabled)) {
        message("Pathway analysis disabled in config (pathway.enabled: false)")
        return(list(annotation = NULL, pathway_results = list(), plot_files = list()))
    }

    de_tables <- de_res$tables
    if (is.null(de_tables) || length(de_tables) == 0) {
        warning("No DE tables available for pathway analysis")
        return(list(annotation = NULL, pathway_results = list(), plot_files = list()))
    }

    # ------------------------------------------------------------------
    # 1. Collect all gene IDs from DE tables
    # ------------------------------------------------------------------
    all_gene_ids <- unique(unlist(lapply(de_tables, function(x) x$FeatureID)))

    # Strip "Gene:" prefix if present (e.g. "Gene:WBGene00000001" -> "WBGene00000001")
    has_prefix <- grepl("^Gene:", all_gene_ids)
    if (mean(has_prefix) > 0.5) {
        message("Stripping 'Gene:' prefix from feature IDs for annotation/pathway")
        # Also strip in each DE table for consistency
        for (cn in names(de_tables)) {
            de_tables[[cn]]$FeatureID <- sub("^Gene:", "", de_tables[[cn]]$FeatureID)
        }
        all_gene_ids <- sub("^Gene:", "", all_gene_ids)
    }

    # ------------------------------------------------------------------
    # Route: local (offline, table-driven) enrichment vs. online fallback.
    # The local path activates ONLY when enrichment.annotation_dir is set;
    # if it is unset/empty the existing online behavior below runs unchanged.
    # ------------------------------------------------------------------
    annotation_dir <- enr_cfg$annotation_dir
    if (!is.null(annotation_dir) && nzchar(annotation_dir)) {
        return(.run_local_enrichment(
            de_tables      = de_tables,
            feature_ids    = all_gene_ids,
            enr_cfg        = enr_cfg,
            config         = config,
            out_dir        = out_dir,
            clustering_res = clustering_res
        ))
    }

    # ------------------------------------------------------------------
    # 2. Auto-detect organism
    # ------------------------------------------------------------------
    organism <- ann_cfg$organism %||% "auto"
    if (tolower(organism) == "auto") {
        org_detect <- detect_organism_from_ids(all_gene_ids)
        if (org_detect$confidence %in% c("high", "medium")) {
            organism <- org_detect$organism
            message("Auto-detected organism: ", organism,
                    " (", org_detect$confidence, " confidence)")
        } else {
            warning("Could not auto-detect organism. Pathway analysis may be limited.")
            organism <- "unknown"
        }
    } else {
        organism <- normalize_organism_name(organism)
    }

    # ------------------------------------------------------------------
    # 3. Annotate genes (custom -> biomaRt -> OrgDb fallback)
    # ------------------------------------------------------------------
    # Ensembl gene IDs from many quantifiers carry a version suffix
    # (e.g. "ENSG00000290825.2"). biomaRt / org.Hs.eg.db key on the
    # unversioned form, so look up against stripped IDs and then map
    # the result back to the original versioned IDs that downstream
    # tables (DE results, report) actually carry.
    has_ensembl_version <- grepl("^ENS[A-Z]*G[0-9]+\\.[0-9]+$", all_gene_ids)
    if (any(has_ensembl_version)) {
        lookup_ids <- sub("\\.[0-9]+$", "", all_gene_ids)
        versioned_for <- setNames(all_gene_ids, lookup_ids)
    } else {
        lookup_ids <- all_gene_ids
        versioned_for <- NULL
    }

    annotation_result <- NULL
    if (!isTRUE(ann_cfg$skip_annotation)) {
        anno_config <- list(
            organism = organism,
            annotation = list(
                skip_annotation = FALSE,
                custom_mapping_file = ann_cfg$custom_mapping_file,
                id_type = ann_cfg$id_type %||% "auto",
                fallback_chain = c("custom", "biomart", "orgdb", "keggrest")
            )
        )
        annotation_result <- annotate_genes_v2(lookup_ids, anno_config, verbose = TRUE)
    }

    annotation_df <- if (!is.null(annotation_result)) annotation_result$annotation else NULL

    # Restore versioned IDs so the saved gene_annotation.csv matches the
    # FeatureID column in DE tables (the report joins on this key).
    if (!is.null(annotation_df) && !is.null(versioned_for)) {
        annotation_df$gene_id <- unname(versioned_for[annotation_df$gene_id])
    }

    # ------------------------------------------------------------------
    # 4. Load gene sets (GO, KEGG, and/or custom GMT)
    # ------------------------------------------------------------------
    databases <- pw_cfg$databases %||% c("GO", "KEGG")
    gmt_file  <- pw_cfg$gmt_file

    gene_sets <- load_gene_sets(
        organism         = organism,
        pathway_database = databases,
        gmt_file         = gmt_file,
        annotation       = annotation_df
    )

    if (length(gene_sets) == 0) {
        message("No gene sets loaded. Skipping pathway analysis.")
        return(list(annotation = annotation_result,
                    pathway_results = list(),
                    plot_files = list()))
    }

    # ------------------------------------------------------------------
    # 5. Run fGSEA / ORA per contrast
    # ------------------------------------------------------------------
    pw_method  <- pw_cfg$method   %||% "fgsea"
    pw_min     <- pw_cfg$min_size %||% 10
    pw_max     <- pw_cfg$max_size %||% 500

    # TODO(simplify-go): GO term simplification was wired here via simplify_go_results
    # (commit 4564b09, dropped by merge 29ffe3e). Restore via cluster_enrichment_terms()
    # in R/core/09_enrichment.R, which has correct score/sim_matrix alignment.
    pathway_results <- run_pathway_analysis(
        de_tables          = de_tables,
        gene_sets          = gene_sets,
        annotation         = annotation_df,
        method             = pw_method,
        min_size           = pw_min,
        max_size           = pw_max
    )

    # ------------------------------------------------------------------
    # 6. Save CSV results + dotplot PNGs
    # ------------------------------------------------------------------
    enrich_dir <- file.path(out_dir, "Enrichment")
    plot_dir   <- file.path(enrich_dir, "plots")

    saved_files <- save_pathway_results(pathway_results, enrich_dir)
    plot_files  <- generate_pathway_plots(pathway_results, plot_dir)

    # Save annotation table
    if (!is.null(annotation_df)) {
        anno_file <- file.path(enrich_dir, "gene_annotation.csv")
        write.csv(annotation_df, anno_file, row.names = FALSE)
        message("Saved gene annotation to: ", anno_file)
    }

    # Build pathway-colored volcano data
    pw_volcano_enabled <- isTRUE(pw_cfg$pathway_volcano) || is.null(pw_cfg$pathway_volcano)
    if (pw_volcano_enabled && length(pathway_results) > 0) {
        message("Building pathway-colored volcano data for RNA-seq...")
        for (cn in names(de_tables)) {
            volcano_data <- build_pathway_volcano_data(de_tables[[cn]], pathway_results)
            if (!is.null(volcano_data)) {
                volcano_file <- file.path(enrich_dir, sprintf("pathway_volcano_data_%s.csv", cn))
                write.csv(volcano_data, volcano_file, row.names = FALSE)
                message("  Saved pathway volcano data: ", volcano_file)
            }
        }
    }

    list(
        annotation      = annotation_result,
        pathway_results = pathway_results,
        plot_files      = plot_files
    )
}


# ==============================================================================
# LOCAL ENRICHMENT PATH (offline, table-driven) — enrichment migration v2
# ==============================================================================

#' Internal: run enrichment using local precomputed KEGG/GO tables
#'
#' Offline path activated by enrichment.annotation_dir. Runs cluster-based ORA
#' (when clustering_res is available) and multi-method GSEA. Returns the same
#' shape as mod_rnaseq_pathway: list(annotation, pathway_results, plot_files),
#' with pathway_results a (possibly empty) named list of data.frames carrying a
#' padj column for downstream compatibility.
#'
#' @param de_tables Named list of per-contrast DE tables.
#' @param feature_ids All unique feature IDs (for the overlap guard).
#' @param enr_cfg The config$modes$rna$enrichment section.
#' @param config Full pipeline config.
#' @param out_dir Output directory for the RNA mode.
#' @param clustering_res Result from mod_rnaseq_clustering(), or NULL.
#' @return list(annotation, pathway_results, plot_files)
#' @noRd
.run_local_enrichment <- function(de_tables, feature_ids, enr_cfg, config,
                                  out_dir, clustering_res = NULL) {

    message("\n=== Local Enrichment (offline, table-driven) ===\n")

    annotation_dir <- enr_cfg$annotation_dir

    # Resolve a relative annotation_dir against the project dir if needed.
    if (!dir.exists(annotation_dir)) {
        project_dir <- config$project$dir
        if (!is.null(project_dir) && nzchar(project_dir)) {
            candidate <- file.path(project_dir, annotation_dir)
            if (dir.exists(candidate)) annotation_dir <- candidate
        }
    }

    # ------------------------------------------------------------------
    # 1. Load local pathway tables (missing DBs are skipped + warned inside)
    # ------------------------------------------------------------------
    databases <- enr_cfg$databases %||% c("KEGG", "GO_BP", "GO_MF", "GO_CC")

    message("Loading local pathway tables from: ", annotation_dir)
    local_tables <- load_local_pathway_tables(
        annotation_dir = annotation_dir,
        databases      = databases,
        feature_ids    = feature_ids
    )

    if (length(local_tables) == 0) {
        warning("No local pathway tables loaded. Returning empty enrichment result.")
        return(list(annotation = NULL, pathway_results = list(), plot_files = list()))
    }

    enrich_dir  <- file.path(out_dir, "Enrichment")
    pval_cutoff <- enr_cfg$pvalue_cutoff %||% 0.05
    padj_method <- enr_cfg$padj_method   %||% "fdr"
    rna_de_cfg  <- config$modes$rna$de %||% list()
    go_simplify <- isTRUE(enr_cfg$go_simplify)
    orgdb       <- enr_cfg$orgdb
    # Single control for enrichment parallelism (ORA + GSEA). <=1 == sequential.
    workers     <- enr_cfg$workers %||% 1
    # Reproducibility: the project's single seed (params$seed) is the source of
    # truth for enrichment RNG. Threaded into run_enrichment_jobs() as an explicit
    # future.seed, it makes GSEA identical across worker counts AND independent
    # rebuilds. Falls back to 1L if params$seed is unset.
    enr_seed    <- config$params$seed %||% 1L

    pathway_results <- list()
    plot_files      <- list()

    # ------------------------------------------------------------------
    # 2. Cluster-based ORA across all gene-list methods
    # ------------------------------------------------------------------
    gene_lists <- build_gene_lists(
        de_tables      = de_tables,
        clustering_res = clustering_res,
        p_cutoff       = rna_de_cfg$p_cutoff %||% 0.05,
        lfc_cutoff     = log2(rna_de_cfg$linear_fc_cutoff %||% 1.5)
    )

    if (length(gene_lists) > 0) {
        message("\n--- Cluster-based ORA ---")

        # Build a flat, method-agnostic ORA job list: one job per
        # (collection x round x database). Each job is an independent enrichment
        # unit; jobs run through the shared run_enrichment_jobs() orchestration
        # layer (parallel compute when workers > 1). Build order (collection ->
        # round -> database) is preserved by run_enrichment_jobs() and by the
        # serial assembly loop below, so ora_results key order is identical to
        # the previous sequential implementation.
        ora_jobs <- list()
        for (clust_method in names(gene_lists)) {
            for (clust_round in names(gene_lists[[clust_method]])) {
                clusters <- gene_lists[[clust_method]][[clust_round]]
                if (length(clusters) == 0) next
                for (db_name in names(local_tables)) {
                    ora_jobs[[length(ora_jobs) + 1]] <- list(
                        clust_method = clust_method,
                        clust_round  = clust_round,
                        db_name      = db_name
                    )
                }
            }
        }

        if (length(ora_jobs) > 0) {
            message("  ", length(ora_jobs), " ORA jobs to run",
                    if (workers > 1) paste0(" (", workers, " workers)") else " (sequential)")

            # Build the pure-compute worker with a MINIMAL captured environment
            # (only the gene sets, tables, and scalar params it uses) — see
            # .make_ora_worker(). This keeps future from serializing unrelated
            # large objects (e.g. the expression matrix in `pre`) to workers.
            run_one_ora_job <- .make_ora_worker(
                gene_lists   = gene_lists,
                local_tables = local_tables,
                pval_cutoff  = pval_cutoff,
                padj_method  = padj_method,
                go_simplify  = go_simplify,
                orgdb        = orgdb
            )

            ora_job_results <- run_enrichment_jobs(ora_jobs, run_one_ora_job, workers, seed = enr_seed)

            # Serial assembly + file writing (deterministic; never in a worker).
            ora_results <- list()
            for (jr in ora_job_results) {
                job         <- jr$job
                res         <- jr$result
                result_base <- paste0(job$db_name, "_", job$clust_method, "_", job$clust_round)

                if (length(res) == 0) {
                    message("    ORA: ", result_base, " — no significant enrichment")
                    next
                }

                write_cluster_ora_outputs(
                    res,
                    outDir      = file.path(enrich_dir, "ORA", job$db_name),
                    file_name   = result_base,
                    type        = jr$db_type,
                    maxCategory = enr_cfg$max_terms_in_dotplot %||% 20
                )

                ora_results <- .store_ora_result(
                    ora_results, res[[3]], paste0(result_base, "_ora"))
                ora_results <- .store_ora_result(
                    ora_results, res[[4]], paste0(result_base, "_ora_simplify"))

                n_terms <- if (!is.null(res[[3]])) nrow(res[[3]]) else 0
                message("    ORA: ", result_base, " — ", n_terms, " enriched terms")
            }

            if (length(ora_results) > 0) pathway_results[["cluster_ora"]] <- ora_results
        }
    } else {
        message("Cluster-based ORA requires clustering results. ",
                "Enable clustering in config to run ORA. GSEA will proceed.")
    }

    # ------------------------------------------------------------------
    # 3. GSEA (multiple ranking methods; parallel when workers > 1)
    # ------------------------------------------------------------------
    message("\n--- GSEA ---")
    ranked_genes <- build_ranked_gene_lists(de_tables)

    gsea_dir  <- file.path(enrich_dir, "GSEA")
    gsea_pval <- enr_cfg$gsea_pvalue_cutoff %||% pval_cutoff
    gsea_padj <- enr_cfg$gsea_padj_method   %||% padj_method
    # `workers` is defined once above and controls both ORA and GSEA.
    # OFF by default (legacy parity); only emitted on demand when enabled.
    per_pathway_artifacts <- isTRUE(enr_cfg$gsea_per_pathway_artifacts)

    gsea_out <- run_gsea_all(
        ranked_genes          = ranked_genes,
        local_tables          = local_tables,
        pvalueCutoff          = gsea_pval,
        pAdjustMethod         = gsea_padj,
        output_dir            = gsea_dir,
        workers               = workers,
        per_pathway_artifacts = per_pathway_artifacts,
        max_terms_in_dotplot  = enr_cfg$max_terms_in_dotplot %||% 20,
        seed                  = enr_seed
    )

    for (contrast in names(gsea_out$results)) {
        if (is.null(pathway_results[[contrast]])) {
            pathway_results[[contrast]] <- gsea_out$results[[contrast]]
        } else {
            pathway_results[[contrast]] <- c(
                pathway_results[[contrast]], gsea_out$results[[contrast]])
        }
    }
    plot_files <- c(plot_files, gsea_out$plot_files)

    message("\n=== Local enrichment complete ===")
    if (length(gene_lists) == 0) message("  ORA: skipped (no gene lists available)")

    list(
        annotation      = NULL,
        pathway_results = pathway_results,
        plot_files      = plot_files
    )
}


#' Store an ORA result table in the accumulator with downstream-compatible columns
#'
#' Adds `padj` (from p.adjust) and `pathway` (from Description) columns when
#' missing, so collect_pipeline_stats()/extract_enrichment_df() can consume it.
#'
#' @param ora_results Current accumulator list.
#' @param df Data.frame from run_cluster_ora(), or NULL.
#' @param key Storage key.
#' @return Updated ora_results list.
#' @noRd
.store_ora_result <- function(ora_results, df, key) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(ora_results)

    if ("p.adjust" %in% colnames(df) && !"padj" %in% colnames(df)) {
        df$padj <- df$p.adjust
    }
    if ("Description" %in% colnames(df) && !"pathway" %in% colnames(df)) {
        df$pathway <- df$Description
    }
    ora_results[[key]] <- df
    ora_results
}
