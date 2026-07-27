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
#' @param orgdb Reserved for future IC-based GO simplify measures (or NULL). GO
#'   simplify itself is mandatory for GO results (Wang/GO.db) — see
#'   run_cluster_ora_compute().
#' @return A function(job) -> list(job, db_type, result), where result is the
#'   4-element list from run_cluster_ora_compute() (or list() if not significant).
#' @noRd
.make_ora_worker <- function(gene_lists, local_tables, pval_cutoff, padj_method,
                             orgdb) {
    force(gene_lists); force(local_tables); force(pval_cutoff)
    force(padj_method); force(orgdb)
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
    annotation_dir_raw <- enr_cfg$annotation_dir
    is_abs      <- grepl("^(/|\\\\|[A-Za-z]:)", annotation_dir_raw)
    project_dir <- config$project$dir
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
        return(list(annotation = NULL, pathway_results = list(), plot_files = list()))
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
