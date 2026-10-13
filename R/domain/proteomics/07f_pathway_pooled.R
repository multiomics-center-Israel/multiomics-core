#' Proteomics enrichment on a Rubin-pooled ranking
#'
#' A second enrichment, next to the legacy one in 07b_pathway.R: proteins are
#' ranked by the z-score of the limma statistic combined across imputation draws
#' with Rubin's rules (pool_limma_imputations()), collapsed to one protein per
#' gene, and tested with fgsea and limma::cameraPR on the same gene sets. Output
#' goes to Enrichment/pooled_ranking/ only, so the legacy tables, plots and the
#' contrast-keyed prot_pathway_res consumed by multiomics are untouched.

#' Output directory of the pooled-ranking enrichment
#'
#' One resolver for the producer, the report and the Methods text.
#'
#' @param out_dir Proteomics output directory.
#' @return Path to \code{Enrichment/pooled_ranking} under \code{out_dir}.
proteomics_pooled_enrichment_dir <- function(out_dir) {
    file.path(out_dir, "Enrichment", "pooled_ranking")
}

#' Build per-contrast gene-level rankings from pooled limma statistics
#'
#' @param de_res Output of mod_proteomics_de() with limma \code{runs}.
#' @param pre Preprocessed proteomics object; \code{expr_filt} supplies the
#'   pre-imputation missing-value mask.
#' @param config Full config list.
#' @return Named list (by contrast) of data.frames, one row per gene (or per
#'   protein when annotation is skipped), with \code{FeatureID},
#'   \code{ProteinID}, \code{stat} (= \code{z_pooled}), \code{log2FoldChange}
#'   (= \code{pooled_logFC}), \code{legacy_stat} (the legacy fgsea ranking value
#'   of the same protein, NA when unavailable) and every column of
#'   pool_limma_imputations().
build_pooled_ranking_tables <- function(de_res, pre, config) {
    cfg <- config$modes$proteomics
    runs <- de_res$runs
    if (length(runs) == 0 || is.null(runs[[1]]$fit2)) {
        stop("build_pooled_ranking_tables(): the DE result has no limma fits to pool ",
             "(de.method must be 'limma').", call. = FALSE)
    }
    skip_symbol_remap <- isTRUE((cfg$annotation %||% list())$skip_annotation)
    observed <- if (!is.null(pre$expr_filt)) !is.na(pre$expr_filt) else NULL
    contrasts <- colnames(runs[[1]]$fit2$coefficients)

    tables <- lapply(contrasts, function(cn) {
        pooled <- pool_limma_imputations(runs, cn, observed = observed)
        pooled$stat <- pooled$z_pooled
        pooled$log2FoldChange <- pooled$pooled_logFC

        legacy <- tryCatch(
            extract_de_table_for_pathway(de_res$summary_df, cn, config),
            error = function(e) NULL)
        pooled$legacy_stat <- if (is.null(legacy)) NA_real_ else
            legacy$stat[match(pooled$FeatureID, legacy$FeatureID)]

        ann_tbl <- runs[[1]]$de_tables[[cn]]
        if (!skip_symbol_remap && !is.null(ann_tbl)) {
            # map_proteins_to_gene_symbols() reads Genes keyed by the configured
            # ID column, which the per-run DE table carries.
            pooled <- map_proteins_to_gene_symbols(pooled, ann_tbl, config, priority = "observed")
        }
        # Unmapped (no Genes column, or annotation skipped): one row per protein.
        if (is.null(pooled$ProteinID)) pooled$ProteinID <- pooled$FeatureID
        pooled
    })
    names(tables) <- normalize_contrast_name(contrasts)
    tables
}

#' Diagnostic plots for one pooled ranking
#'
#' Two panels written as separate PNGs: the fraction of missing information
#' against the fraction of imputed values, and the rank of the pooled z-score
#' against the rank of the legacy ranking statistic.
#'
#' @param ranking One element of build_pooled_ranking_tables().
#' @param contrast Contrast name, used in titles and file names.
#' @param plots_dir Directory for the PNG files.
#' @return Character vector of written file paths.
plot_pooled_ranking_diagnostics <- function(ranking, contrast, plots_dir) {
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
    clean <- gsub("[^a-zA-Z0-9_-]", "_", contrast)
    files <- character(0)

    d1 <- ranking[is.finite(ranking$fmi) & is.finite(ranking$frac_imputed), , drop = FALSE]
    if (nrow(d1) > 0) {
        p1 <- ggplot2::ggplot(d1, ggplot2::aes(x = frac_imputed, y = fmi)) +
            ggplot2::geom_jitter(width = 0.01, height = 0, alpha = 0.4, size = 1) +
            ggplot2::labs(title = paste("Fraction of missing information:", contrast),
                          x = "Fraction of imputed values (contrast samples)",
                          y = "Fraction of missing information") +
            ggplot2::theme_minimal()
        f1 <- file.path(plots_dir, paste0("diagnostic_fmi_", clean, ".png"))
        ggplot2::ggsave(f1, p1, width = 7, height = 5)
        files <- c(files, f1)
    }

    d2 <- ranking[is.finite(ranking$z_pooled) & is.finite(ranking$legacy_stat), , drop = FALSE]
    if (nrow(d2) > 2) {
        rho <- stats::cor(d2$z_pooled, d2$legacy_stat, method = "spearman")
        d2$rank_pooled <- rank(-d2$z_pooled)
        d2$rank_legacy <- rank(-d2$legacy_stat)
        p2 <- ggplot2::ggplot(d2, ggplot2::aes(x = rank_legacy, y = rank_pooled)) +
            ggplot2::geom_point(alpha = 0.4, size = 1) +
            ggplot2::labs(title = paste("Pooled z vs legacy ranking:", contrast),
                          subtitle = sprintf("Spearman correlation %.3f", rho),
                          x = "Rank by legacy statistic",
                          y = "Rank by pooled z") +
            ggplot2::theme_minimal()
        f2 <- file.path(plots_dir, paste0("diagnostic_rank_", clean, ".png"))
        ggplot2::ggsave(f2, p2, width = 7, height = 5)
        files <- c(files, f2)
    }
    files
}

#' Run fgsea and cameraPR on the Rubin-pooled ranking
#'
#' @param de_res Output of mod_proteomics_de() with limma \code{runs}.
#' @param pre Preprocessed proteomics object.
#' @param config Full config list.
#' @param out_dir Proteomics output directory.
#' @return list(ranking, results, settings): \code{ranking} is
#'   build_pooled_ranking_tables(); \code{results} is a run_pathway_analysis()
#'   list with \code{<db>_fgsea} and \code{<db>_camerapr} tables (empty when no
#'   gene set loaded); \code{settings} records \code{m},
#'   \code{inter_gene_cor}, \code{use_ranks}, \code{collapse_rule} and
#'   \code{df_com_rule}.
run_proteomics_pooled_enrichment <- function(de_res, pre, config, out_dir) {
    cfg <- config$modes$proteomics
    pw_cfg <- cfg$pathway %||% list()
    cam_cfg <- (pw_cfg$pooled_ranking %||% list())$camerapr %||% list()
    inter_gene_cor <- cam_cfg$inter_gene_cor %||% 0.01
    use_ranks <- isTRUE(cam_cfg$use_ranks)

    message("=== Proteomics enrichment on the pooled ranking ===")
    ranking <- build_pooled_ranking_tables(de_res, pre, config)

    pooled_dir <- proteomics_pooled_enrichment_dir(out_dir)
    plots_dir <- file.path(pooled_dir, "plots")
    dir.create(pooled_dir, recursive = TRUE, showWarnings = FALSE)
    for (cn in names(ranking)) {
        tbl <- ranking[[cn]]
        tbl <- tbl[order(-tbl$z_pooled, na.last = TRUE), , drop = FALSE]
        clean <- gsub("[^a-zA-Z0-9_-]", "_", cn)
        utils::write.csv(tbl, file.path(pooled_dir, paste0("pooled_ranking_", clean, ".csv")),
                         row.names = FALSE)
        plot_pooled_ranking_diagnostics(tbl, cn, plots_dir)
    }

    settings <- list(
        m              = length(de_res$runs),
        inter_gene_cor = inter_gene_cor,
        use_ranks      = use_ranks,
        collapse_rule  = "most_observed_then_AveExpr",
        df_com_rule    = pooled_df_com_rule()
    )

    results <- list()
    gs_info <- prepare_proteomics_gene_sets(ranking, config)
    if (!is.null(gs_info)) {
        results <- run_pathway_analysis(
            de_tables               = ranking,
            gene_sets               = gs_info$gene_sets,
            annotation              = gs_info$annotation,
            method                  = c("fgsea", "camerapr"),
            min_size                = pw_cfg$min_size %||% 10,
            max_size                = pw_cfg$max_size %||% 500,
            seed                    = config$params$seed %||% 1L,
            camerapr_inter_gene_cor = inter_gene_cor,
            camerapr_use_ranks      = use_ranks
        )
        save_pathway_results(results, pooled_dir)
        generate_pathway_plots(results, plots_dir)
    }

    message("Pooled-ranking enrichment complete.")
    list(ranking = ranking, results = results, settings = settings)
}
