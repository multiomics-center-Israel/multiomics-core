#' Build a pure-compute ORA worker with a minimal captured environment
#'
#' Returns a `function(job)` that computes cluster ORA for one job. Defining it
#' here (not inside mod_rnaseq_pathway) bounds the closure's environment to just
#' the three scalar thresholds below — it captures NO large objects. Everything
#' the job needs travels IN the job: its `clusters` (gene-list for this
#' collection/round) and its OWN single database's `term2gene`/`term2name`. The
#' worker never captures the whole multi-database `local_tables` (nor the
#' expression matrix in `pre`), so future serializes at most one database's
#' tables per job. Pure computation only: no file I/O, no plotting, no messages.
#'
#' The returned value carries a SLIM `job` (only the coordinate fields the serial
#' assembly stage needs) — NOT the heavy `clusters`/`term2gene`/`term2name` — so
#' the large tables are not serialized back from parallel workers to the master.
#'
#' @param pval_cutoff,padj_method ORA thresholds passed to run_cluster_ora_compute().
#' @param orgdb Reserved for future IC-based GO simplify measures (or NULL). GO
#'   simplify itself is mandatory for GO results (Wang/GO.db) — see
#'   run_cluster_ora_compute().
#' @return A function(job) -> list(job, db_type, result), where `job` is the slim
#'   coordinate record and result is the 4-element list from
#'   run_cluster_ora_compute() (or list() if not significant).
#' @noRd
.make_ora_worker <- function(pval_cutoff, padj_method, orgdb) {
    force(pval_cutoff); force(padj_method); force(orgdb)
    function(job) {
        db_type <- if (grepl("^GO", job$db_name)) "GO" else "KEGG"
        db_ont  <- if (db_type == "GO") sub("^GO_", "", job$db_name) else NULL

        res <- run_cluster_ora_compute(
            clusters      = job$clusters,
            TERM2GENE     = job$term2gene,
            TERM2NAME     = job$term2name,
            type          = db_type,
            pvalueCutoff  = pval_cutoff,
            pAdjustMethod = padj_method,
            orgdb         = orgdb,
            ont           = db_ont
        )
        # Return only the coordinate fields (not the heavy per-job payload) so
        # parallel workers don't serialize the tables back to the master.
        list(job = list(clust_method = job$clust_method,
                        clust_round  = job$clust_round,
                        db_name      = job$db_name),
             db_type = db_type, result = res)
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
#' @param clustering_res Result from mod_rnaseq_clustering(), or NULL. Its
#'   NULL-ness is the single-omics/multiomics discriminator, NOT a
#'   clustering-enabled flag: the single-omics pipeline always passes the
#'   clustering object (a non-NULL list even when clustering is disabled),
#'   while the multiomics pipeline passes NULL. Used only by the local (offline)
#'   enrichment path. When non-NULL (single-omics): contrast-based ORA always
#'   runs from the DE tables, and cluster-based ORA additionally runs when the
#'   object actually carries clusters. When NULL (multiomics): all ORA is
#'   skipped with a warning and only GSEA runs. Ignored by the online fallback.
#' @return List with annotation, pathway_results, and plot_files
#' @export
mod_rnaseq_pathway <- function(de_res, pre, config, out_dir, clustering_res = NULL) {

    rna_cfg <- config$modes$rna
    ann_cfg <- rna_cfg$annotation %||% list()
    pw_cfg  <- rna_cfg$pathway
    enr_cfg <- rna_cfg$enrichment %||% list()

    # Decide whether to run enrichment. Skip when either switch is explicitly
    # false (`pathway.enabled: false` or `enrichment.enabled: false`), OR when the
    # config omits BOTH blocks entirely — a legacy/minimal RNA config never ran
    # enrichment, and proceeding would introduce surprise network annotation / I/O
    # (pinned by tests/testthat/test-config-defaults.R). A standalone `enrichment:`
    # block still opts in even without a legacy `pathway:` block (e.g. the shipped
    # template), and pre-existing configs that carry a `pathway:` block keep running.
    # Returns the standard empty-safe shape so the rest of the RNA pipeline (and
    # all downstream consumers of pathway_results) continue without failure.
    pathway_off    <- isFALSE(pw_cfg$enabled)
    enrichment_off <- isFALSE(enr_cfg$enabled)
    both_absent    <- is.null(pw_cfg) && is.null(rna_cfg$enrichment)
    if (both_absent || pathway_off || enrichment_off) {
        reason <- if (both_absent) "no pathway/enrichment config"
                  else if (pathway_off) "pathway.enabled: false"
                  else "enrichment.enabled: false"
        message("RNA enrichment disabled (", reason,
                "); skipping ORA/GSEA and continuing the pipeline.")
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
            clustering_res = clustering_res,
            pre            = pre
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
#' @param pre Preprocessed data list (expr_work normalized matrix + meta). Used
#'   only in the serial output stage to build the rich per-pathway core-gene
#'   tables and (opt-in) per-pathway expression heatmaps; never sent to workers.
#' @return list(annotation, pathway_results, plot_files)
#' @noRd
.run_local_enrichment <- function(de_tables, feature_ids, enr_cfg, config,
                                  out_dir, clustering_res = NULL, pre = NULL) {

    message("\n=== Local Enrichment (offline, table-driven) ===\n")

    # Resolve annotation_dir. Absolute paths (Unix "/...", Windows "C:\..."/UNC)
    # are used as given; RELATIVE paths are ALWAYS resolved against the configured
    # project directory (config$project$dir) — never a same-named directory that
    # happens to sit under the process working directory. The resolved path is the
    # one validated and reported below.
    # Expand a leading "~" before anything else: it must happen before the
    # absolute-path test (an unexpanded "~/..." looks relative and would wrongly
    # be joined onto project_dir) and before the file.path() join for the
    # relative case.
    annotation_dir_raw <- path.expand(enr_cfg$annotation_dir)
    is_abs      <- grepl("^(/|\\\\|[A-Za-z]:)", annotation_dir_raw)
    project_dir <- config$project$dir
    if (!is.null(project_dir) && nzchar(project_dir)) project_dir <- path.expand(project_dir)
    annotation_dir <- if (is_abs || is.null(project_dir) || !nzchar(project_dir)) {
        annotation_dir_raw
    } else {
        file.path(project_dir, annotation_dir_raw)
    }

    if (!dir.exists(annotation_dir)) {
        stop("Local/offline enrichment was requested ",
             "(modes.rna.enrichment.annotation_dir = \"", annotation_dir_raw, "\"), ",
             "but the resolved annotation directory does not exist:\n  ",
             normalizePath(annotation_dir, winslash = "/", mustWork = FALSE), "\n",
             "Provide the local annotation tables at that path, or set ",
             "modes.rna.enrichment.annotation_dir: \"\" to use the online ",
             "enrichment workflow.", call. = FALSE)
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
        return(list(annotation = NULL, pathway_results = list(), plot_files = list(),
                    enrichment_manifest = .empty_enrichment_manifest(),
                    enrichment_index    = .empty_enrichment_index(),
                    gsea_rankings       = list(),
                    pathway_membership  = list()))
    }

    enrich_dir  <- file.path(out_dir, "Enrichment")
    pval_cutoff <- enr_cfg$pvalue_cutoff %||% 0.05
    padj_method <- enr_cfg$padj_method   %||% "fdr"
    rna_de_cfg  <- config$modes$rna$de %||% list()
    # GO simplify is MANDATORY for GO results and NOT configurable: the compute
    # layer always collapses redundant GO terms with the Wang measure over GO.db
    # (no OrgDb, no network). Any legacy `enrichment.go_simplify` field is ignored.
    # If GO.db/GOSemSim are missing it warns and keeps the unsimplified GO ORA —
    # the pipeline never fails. `orgdb` stays reserved for future IC-based measures.
    orgdb       <- enr_cfg$orgdb
    # Single control for enrichment parallelism (ORA + GSEA). <=1 == sequential.
    workers     <- enr_cfg$workers %||% 1
    # Reproducibility: the project's single seed (params$seed) is the source of
    # truth for enrichment RNG. Threaded into run_enrichment_jobs() as an explicit
    # future.seed, it makes GSEA identical across worker counts AND independent
    # rebuilds. Falls back to 1L if params$seed is unset.
    enr_seed    <- config$params$seed %||% 1L
    # Enrichment output toggles. Rich-by-default: all default TRUE except the
    # high-volume per-pathway heatmaps (opt-in). Set any to false to disable.
    #   modes.rna.enrichment.plots.<name> and .gsea.per_pathway_artifacts
    plots_cfg           <- enr_cfg$plots %||% list()
    gsea_cfg            <- enr_cfg$gsea  %||% list()
    do_dotplot          <- isTRUE(plots_cfg$dotplot             %||% TRUE)
    do_ridgeplot        <- isTRUE(plots_cfg$ridgeplot           %||% TRUE)
    do_ridgeplot_all    <- isTRUE(plots_cfg$ridgeplot_all_genes %||% TRUE)
    do_pathway_heatmaps <- isTRUE(plots_cfg$pathway_heatmaps    %||% FALSE)
    do_shared_genes     <- isTRUE(plots_cfg$shared_genes        %||% TRUE)
    # per-pathway gseaplot2 + rich core-gene tables. New key gsea.per_pathway_artifacts
    # (default TRUE); the legacy flat key enrichment.gsea_per_pathway_artifacts is
    # still honored for backward compatibility. `%||%` chains on NULL, so an
    # explicit `false` at either key is respected.
    per_pathway_artifacts <- isTRUE(
        gsea_cfg$per_pathway_artifacts %||%
        enr_cfg$gsea_per_pathway_artifacts %||% TRUE)

    pathway_results <- list()
    plot_files      <- list()
    # Stage 3A: availability-manifest + storage-index accumulators. ORA rows are
    # built here (before empty units are dropped); GSEA rows arrive from
    # run_gsea_all() below. Combined into rna_pathway_res$enrichment_manifest /
    # $enrichment_index at return.
    ora_manifest_rows <- list()
    ora_index_rows    <- list()

    # ------------------------------------------------------------------
    # 2. Cluster-based ORA across all gene-list methods
    # ------------------------------------------------------------------
    # Contract (see the clustering_res @param): this guard is the multiomics ORA
    # skip. In multiomics mode (clustering_res = NULL) all ORA is skipped and only
    # GSEA runs; without it build_gene_lists() would still build the contrast-derived
    # collections from the DE tables alone (they need no clustering) and ORA would
    # run contrary to the documented behavior. In single-omics mode clustering_res is
    # always a non-NULL list, so build_gene_lists() runs: contrast-based ORA always,
    # cluster-based ORA additionally when the object carries clusters.
    gene_lists <- if (is.null(clustering_res)) {
        list()
    } else {
        build_gene_lists(
            de_tables      = de_tables,
            clustering_res = clustering_res,
            p_cutoff       = rna_de_cfg$p_cutoff %||% 0.05,
            lfc_cutoff     = log2(rna_de_cfg$linear_fc_cutoff %||% 1.5)
        )
    }

    if (length(gene_lists) > 0) {
        message("\n--- Cluster-based ORA ---")

        # Build a flat, method-agnostic ORA job list: one job per
        # (collection x round x database). Each job is an independent enrichment
        # unit; jobs run through the shared run_enrichment_jobs() orchestration
        # layer (parallel compute when workers > 1). Build order (collection ->
        # round -> database) is preserved by run_enrichment_jobs() and by the
        # serial assembly loop below, so ora_results key order is identical to
        # the previous sequential implementation.
        # Each job carries EVERYTHING it needs: its `clusters` (tiny) and its OWN
        # single database's TERM2GENE/TERM2NAME. `local_tables[[db]]$TERM2GENE` is
        # assigned by REFERENCE (copy-on-write) — jobs sharing a database point at
        # the same object in the master, so this does NOT duplicate the annotation
        # tables in memory; only the single database in play is serialized when a
        # job is dispatched to a worker. The worker (see .make_ora_worker) thus
        # captures no multi-database `local_tables`.
        ora_jobs <- list()
        for (clust_method in names(gene_lists)) {
            for (clust_round in names(gene_lists[[clust_method]])) {
                clusters <- gene_lists[[clust_method]][[clust_round]]
                if (length(clusters) == 0) next
                for (db_name in names(local_tables)) {
                    ora_jobs[[length(ora_jobs) + 1]] <- list(
                        clust_method = clust_method,
                        clust_round  = clust_round,
                        db_name      = db_name,
                        clusters     = clusters,
                        term2gene    = local_tables[[db_name]]$TERM2GENE,
                        term2name    = local_tables[[db_name]]$TERM2NAME
                    )
                }
            }
        }

        if (length(ora_jobs) > 0) {
            message("  ", length(ora_jobs), " ORA jobs to run",
                    if (workers > 1) paste0(" (", workers, " workers)") else " (sequential)")

            # Build the pure-compute worker. It captures ONLY the scalar
            # thresholds — everything else (clusters + this job's single-database
            # tables) rides in each job (see .make_ora_worker()), so future never
            # serializes the multi-database `local_tables` or the expression
            # matrix in `pre`.
            run_one_ora_job <- .make_ora_worker(
                pval_cutoff  = pval_cutoff,
                padj_method  = padj_method,
                orgdb        = orgdb
            )

            # ORA is deterministic (enricher + Wang GO-simplify use no RNG), so at
            # workers <= 1 we take the true base::lapply path
            # (prefer_lapply_when_sequential = TRUE): no future global discovery,
            # no future.globals.maxSize guard. workers > 1 still uses future.
            ora_job_results <- run_enrichment_jobs(
                ora_jobs, run_one_ora_job, workers, seed = enr_seed,
                prefer_lapply_when_sequential = TRUE)

            # Serial assembly + file writing (deterministic; never in a worker).
            ora_results <- list()
            for (jr in ora_job_results) {
                job         <- jr$job
                res         <- jr$result
                result_base <- paste0(job$db_name, "_", job$clust_method, "_", job$clust_round)
                storage_key <- paste0(result_base, "_ora")
                simplify_key <- paste0(result_base, "_ora_simplify")
                has_simpl   <- identical(jr$db_type, "GO")
                main_df     <- if (length(res) >= 3) res[[3]] else NULL
                clusters    <- gene_lists[[job$clust_method]][[job$clust_round]]

                # Availability-manifest rows for EVERY evaluated unit (captured
                # before the empty-unit `next` below, so evaluated-but-empty
                # clusters/patterns/contrasts get n_significant = 0 rows).
                ora_manifest_rows[[length(ora_manifest_rows) + 1]] <-
                    .build_ora_manifest_rows(
                        collection = job$clust_method, round = job$clust_round,
                        db_name = job$db_name, main_df = main_df,
                        has_simplify = has_simpl, clusters = clusters,
                        storage_key = storage_key)

                if (length(res) == 0) {
                    message("    ORA: ", result_base, " — no significant enrichment")
                    next
                }

                unit_dir <- ora_unit_dir(
                    file.path(enrich_dir, "ORA"),
                    job$db_name, job$clust_method, job$clust_round)
                write_cluster_ora_outputs(
                    res,
                    unit_dir     = unit_dir,
                    type         = jr$db_type,
                    maxCategory  = enr_cfg$max_terms_in_dotplot %||% 20,
                    dotplot      = do_dotplot,
                    shared_genes = do_shared_genes
                )

                ora_results <- .store_ora_result(
                    ora_results, res[[3]], storage_key)
                ora_results <- .store_ora_result(
                    ora_results, res[[4]], simplify_key)

                # Storage-index row (table-level) only when a table was stored.
                if (storage_key %in% names(ora_results)) {
                    ora_index_rows[[length(ora_index_rows) + 1]] <- data.frame(
                        analysis = "ORA", database = job$db_name,
                        group = job$clust_method, item = job$clust_round,
                        container = "cluster_ora", storage_key = storage_key,
                        has_simplify = has_simpl,
                        simplify_key = if (simplify_key %in% names(ora_results))
                            simplify_key else NA_character_,
                        stringsAsFactors = FALSE)
                }

                n_terms <- if (!is.null(res[[3]])) nrow(res[[3]]) else 0
                message("    ORA: ", result_base, " — ", n_terms, " enriched terms")
            }

            if (length(ora_results) > 0) pathway_results[["cluster_ora"]] <- ora_results
        }
    } else if (is.null(clustering_res)) {
        # Multiomics mode: ORA is intentionally skipped (see the clustering_res
        # @param contract); GSEA still runs below.
        message("Cluster-based ORA skipped (multiomics mode: no clustering ",
                "result); running GSEA only.")
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
    # `per_pathway_artifacts` and the plot toggles are resolved once at the top.

    # Rich per-pathway core-gene tables + expression heatmaps need a per-gene
    # context (expression + DE stats + z-scores) and the normalized expression
    # matrix + sample annotation. Built ONCE here in the serial stage from
    # de_tables + `pre`; passed only into run_gsea_all()'s serial writer (never
    # into a worker). NULL when `pre` is unavailable -> falls back to the minimal
    # gene/rank_value table and skips heatmaps (fail-soft).
    gene_context   <- NULL
    heatmap_expr   <- NULL
    annotation_col <- NULL
    if ((per_pathway_artifacts || do_pathway_heatmaps) && !is.null(pre)) {
        gene_context <- tryCatch(
            .build_gene_context(de_tables, pre, rna_de_cfg),
            error = function(e) {
                warning("Could not build per-pathway gene context: ",
                        conditionMessage(e)); NULL
            })
        if (do_pathway_heatmaps && !is.null(pre$expr_work)) {
            heatmap_expr   <- pre$expr_work
            annotation_col <- tryCatch(
                build_heatmap_annotation_col(pre$meta, config$modes$rna),
                error = function(e) NULL)
        }
    }

    gsea_out <- run_gsea_all(
        ranked_genes          = ranked_genes,
        local_tables          = local_tables,
        pvalueCutoff          = gsea_pval,
        pAdjustMethod         = gsea_padj,
        output_dir            = gsea_dir,
        workers               = workers,
        per_pathway_artifacts = per_pathway_artifacts,
        max_terms_in_dotplot  = enr_cfg$max_terms_in_dotplot %||% 20,
        dotplot               = do_dotplot,
        ridgeplot             = do_ridgeplot,
        ridgeplot_all_genes   = do_ridgeplot_all,
        pathway_heatmaps      = do_pathway_heatmaps,
        gene_context          = gene_context,
        expr_mat              = heatmap_expr,
        annotation_col        = annotation_col,
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

    # ------------------------------------------------------------------
    # 4. Stage 3A: preserve enrichment metadata + biological primitives that
    #    were previously discarded. These are ADDITIVE top-level siblings; the
    #    existing keys (annotation/pathway_results/plot_files) are unchanged, and
    #    no biological result, threshold, or output file is affected. Integer
    #    encoding against a shared gene_index is a Stage 3B (export) concern.
    # ------------------------------------------------------------------
    ora_manifest <- if (length(ora_manifest_rows) > 0)
        do.call(rbind, ora_manifest_rows) else .empty_enrichment_manifest()
    ora_index    <- if (length(ora_index_rows) > 0)
        do.call(rbind, ora_index_rows) else .empty_enrichment_index()

    enrichment_manifest <- rbind(
        ora_manifest, gsea_out$manifest %||% .empty_enrichment_manifest())
    enrichment_index <- rbind(
        ora_index, gsea_out$index %||% .empty_enrichment_index())

    # Exact ranked vectors used by GSEA (shared across databases — verified in
    # run_gsea_all: job$ranked = ranked_genes[[method]][[contrast]] for every db,
    # no db-specific filtering). Returned verbatim; never recomputed downstream.
    gsea_rankings <- ranked_genes

    # GSEA pathway membership: significant-pathway union per db from the exact
    # local TERM2GENE tables GSEA used.
    pathway_membership <- .build_gsea_pathway_membership(gsea_out$results, local_tables)

    message("\n=== Local enrichment complete ===")
    if (length(gene_lists) == 0) message("  ORA: skipped (no gene lists available)")

    list(
        annotation          = NULL,
        pathway_results     = pathway_results,
        plot_files          = plot_files,
        enrichment_manifest = enrichment_manifest,
        enrichment_index    = enrichment_index,
        gsea_rankings       = gsea_rankings,
        pathway_membership  = pathway_membership
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


#' Build availability-manifest rows for one ORA unit (collection x round x db)
#'
#' Produces one manifest row per USER-FACING item, capturing evaluated-but-empty
#' units (n_significant = 0) that are otherwise dropped from pathway_results:
#'   - partition / binary_patterns: one row per evaluated cluster/pattern label
#'     (from the evaluated gene-list, NOT from the result), n_significant counted
#'     from the result's `Cluster` column (0 when the label is absent);
#'   - contrasts / contrasts_wo_direction / all_DE: one row for the round
#'     (contrast or `any_contrast`), n_significant = total rows in the unit table.
#'
#' @param collection Collection name (clust_method).
#' @param round Round within the collection (contrast / any_contrast / k / best).
#' @param db_name Database.
#' @param main_df The stored ORA result table for this unit (res[[3]]), or NULL.
#' @param has_simplify Logical (TRUE for GO databases).
#' @param clusters The evaluated gene-list for this unit (named vector
#'   feature -> label); its unique labels are the evaluated cluster/pattern set.
#' @param storage_key The `<db>_<collection>_<round>_ora` key.
#' @return A data.frame of manifest rows (schema of .empty_enrichment_manifest()).
#' @noRd
.build_ora_manifest_rows <- function(collection, round, db_name, main_df,
                                     has_simplify, clusters, storage_key) {
    # ORA compute is not wrapped in a per-unit tryCatch (unlike GSEA), so an ORA
    # error propagates and aborts the run loudly rather than being silently
    # recorded — there is no "failed" ORA manifest state to misrepresent. Every
    # ORA row is therefore a SUCCESSFUL evaluation: "significant" (n_sig > 0) or
    # "empty" (n_sig == 0).
    mk <- function(item, n_sig) data.frame(
        analysis = "ORA", database = db_name, group = collection,
        item = as.character(item), evaluated = TRUE,
        status = if (n_sig > 0) "significant" else "empty",
        n_significant = as.integer(n_sig), has_simplify = has_simplify,
        storage_key = storage_key, stringsAsFactors = FALSE)

    if (collection %in% c("partition", "binary_patterns")) {
        labels <- sort(unique(as.character(clusters)))
        if (length(labels) == 0) return(.empty_enrichment_manifest())
        rows <- lapply(labels, function(lab) {
            n <- if (!is.null(main_df) && "Cluster" %in% names(main_df)) {
                sum(as.character(main_df$Cluster) == lab)
            } else 0L
            mk(lab, n)
        })
        do.call(rbind, rows)
    } else {
        n <- if (!is.null(main_df)) nrow(main_df) else 0L
        mk(round, n)
    }
}


#' Preserve GSEA pathway membership (union of significant GSEA pathways per db)
#'
#' From the assembled GSEA result tables, collects the union of pathway IDs that
#' are significant in at least one unit PER database, then subsets the exact
#' local TERM2GENE table used by GSEA to those pathways. This is the natural
#' (character) representation; integer encoding against a shared gene_index is a
#' Stage 3B (export) concern. No ORA memberships are added (ORA keeps `geneID`),
#' and membership is stored once per database (never per ranking/contrast).
#'
#' @param gsea_results The nested `contrast -> key -> data.frame` GSEA results.
#' @param local_tables Output of load_local_pathway_tables() (per-db TERM2GENE).
#' @return list(database -> list(pathway_id -> character gene IDs)); empty list
#'   when no GSEA pathways are significant.
#' @noRd
.build_gsea_pathway_membership <- function(gsea_results, local_tables) {
    ids_by_db <- list()
    for (contrast in names(gsea_results)) {
        for (key in names(gsea_results[[contrast]])) {
            df <- gsea_results[[contrast]][[key]]
            if (is.null(df) || !all(c("ID", "database") %in% names(df)) ||
                nrow(df) == 0) next
            db <- as.character(df$database[1])
            ids_by_db[[db]] <- unique(c(ids_by_db[[db]], as.character(df$ID)))
        }
    }

    membership <- list()
    for (db in names(ids_by_db)) {
        t2g <- local_tables[[db]]$TERM2GENE
        if (is.null(t2g) || !all(c("term", "gene") %in% names(t2g))) next
        keep <- as.character(t2g$term) %in% ids_by_db[[db]]
        if (!any(keep)) next
        membership[[db]] <- split(as.character(t2g$gene[keep]),
                                  as.character(t2g$term[keep]))
    }
    membership
}


#' Build a per-gene context table for rich per-pathway core-gene outputs
#'
#' Assembles one row per feature carrying all-contrast DE statistics (reusing
#' \code{build_rnaseq_summary_df()} so no DE-column logic is duplicated) plus the
#' normalized per-sample expression and per-gene z-scores. Keyed by the feature
#' ID in the first column so the core-gene writer can left-join by core gene.
#' This is the offline-path equivalent of the columns in \code{final_results.tsv}
#' (annotation columns are omitted — RNA final_results carries none either).
#'
#' @param de_tables Named list of per-contrast DE tables.
#' @param pre Preprocessed data list (uses \code{pre$expr_work}).
#' @param de_cfg The \code{modes.rna.de} config section.
#' @return A data.frame (first column = FeatureID) or NULL on failure/empty.
#' @noRd
.build_gene_context <- function(de_tables, pre, de_cfg = list()) {
    summ <- tryCatch(build_rnaseq_summary_df(de_tables, de_cfg),
                     error = function(e) NULL)
    if (is.null(summ) || !is.data.frame(summ) || nrow(summ) == 0) return(NULL)

    id_col <- if ("FeatureID" %in% names(summ)) "FeatureID" else names(summ)[1]
    summ   <- summ[, c(id_col, setdiff(names(summ), id_col)), drop = FALSE]

    expr <- pre$expr_work
    if (!is.null(expr)) {
        expr <- as.matrix(expr)
        # Match on the `Gene:`-normalized key on BOTH sides: DE FeatureIDs may have
        # had the prefix stripped upstream while expr_work rownames keep it (or vice
        # versa). Normalizing both aligns them regardless. `summ`'s own (original)
        # IDs are preserved in the output; only the match key is normalized.
        # Ambiguity guard: if two expr rows normalize to the SAME key (e.g. both
        # "Gene:X" and "X" present) we cannot know which is X's expression, so we
        # EXCLUDE those keys from matching (their expression stays NA) rather than
        # silently pick one. Unique matches are unaffected.
        expr_keys <- sub("^Gene:", "", rownames(expr))
        dup <- duplicated(expr_keys) | duplicated(expr_keys, fromLast = TRUE)
        if (any(dup)) {
            ambig <- unique(expr_keys[dup])
            warning(".build_gene_context(): ", length(ambig), " gene ID(s) are ",
                    "ambiguous after Gene: normalization (e.g. ",
                    paste(utils::head(ambig, 5L), collapse = ", "),
                    "); their expression is left NA to avoid a wrong match.")
            expr_keys[dup] <- NA_character_   # ambiguous keys never match a gene
        }
        idx  <- match(sub("^Gene:", "", summ[[id_col]]), expr_keys)
        sub  <- expr[idx, , drop = FALSE]           # NA rows for unmatched genes
        # Per-gene z-scores across samples (scale() works column-wise -> transpose).
        z <- tryCatch(t(scale(t(sub))), error = function(e) NULL)
        expr_df <- as.data.frame(sub, stringsAsFactors = FALSE, check.names = FALSE)
        summ <- cbind(summ, expr_df)
        if (!is.null(z)) {
            colnames(z) <- paste0(colnames(expr), ".zscore")
            summ <- cbind(summ, as.data.frame(z, stringsAsFactors = FALSE,
                                              check.names = FALSE))
        }
    }
    rownames(summ) <- NULL
    summ
}
