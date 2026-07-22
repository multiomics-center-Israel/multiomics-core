#' Metabolomics Commentary Module
#'
#' Generates AI-powered or data-driven commentary for metabolomics figures.
#' Adapted from the lipidomics commentary module.
#'
#' @param de_res      Result from metabolomics DE
#' @param qc_pre_obj  Result from metabolomics QC pre
#' @param config      Full pipeline config
#' @param out_dir     Metabolomics output directory
#' @return Path to commentary_all.json (character, format = "file")
mod_metabolomics_commentary <- function(de_res, qc_pre_obj = NULL, config, out_dir) {

    metab_cfg    <- config$modes$metabolomics
    comment_cfg  <- metab_cfg$commentary %||% list()

    # Skip if commentary is not enabled
    if (!isTRUE(comment_cfg$enabled) || identical(comment_cfg$backend, "none")) {
        message("Metabolomics commentary disabled in config")
        commentary_dir <- file.path(out_dir, "commentary")
        dir.create(commentary_dir, recursive = TRUE, showWarnings = FALSE)
        out_file <- file.path(commentary_dir, "commentary_all.json")
        jsonlite::write_json(list(), out_file, auto_unbox = TRUE)
        return(out_file)
    }

    # 1. Discover all figures
    figures_tbl <- discover_metabolomics_figures(out_dir, config)

    if (nrow(figures_tbl) == 0) {
        message("No metabolomics figures found for commentary.")
        commentary_dir <- file.path(out_dir, "commentary")
        dir.create(commentary_dir, recursive = TRUE, showWarnings = FALSE)
        out_file <- file.path(commentary_dir, "commentary_all.json")
        jsonlite::write_json(list(), out_file, auto_unbox = TRUE)
        return(out_file)
    }

    # 2. Build enriched data for richer commentary
    extra_data <- build_metabolomics_extra_data(de_res, qc_pre_obj, config)

    # 3. Override config so generate_all_commentary reads commentary settings
    config_for_commentary <- config
    if (is.null(config_for_commentary$modes$rna))
        config_for_commentary$modes$rna <- list()
    config_for_commentary$modes$rna$commentary <- comment_cfg

    # 4. Generate commentary
    commentary_dir <- file.path(out_dir, "commentary")
    commentary_list <- generate_all_commentary(
        figures_tbl = figures_tbl,
        config      = config_for_commentary,
        de_summary  = NULL,
        extra_data  = extra_data,
        output_dir  = commentary_dir
    )

    out_file <- file.path(commentary_dir, "commentary_all.json")
    message("Metabolomics commentary complete: ", length(commentary_list), " figures processed")
    out_file
}


#' Discover metabolomics figures for commentary
#'
#' Scans the output directory for PNG plots and returns a table of figures.
#'
#' @param out_dir  Metabolomics output directory
#' @param config   Full pipeline config (unused, kept for interface consistency)
#' @return data.frame with columns: figure_id, filepath, plot_type, section, etc.
discover_metabolomics_figures <- function(out_dir, config = NULL) {

    diag_dir <- file.path(out_dir, "Diagnostic_plots")
    figures <- list()

    add_if_exists <- function(id, filename, dir, plot_type, section, title, desc,
                               x_axis = NA, y_axis = NA) {
        fp <- file.path(dir, filename)
        if (file.exists(fp)) {
            figures[[id]] <<- data.frame(
                figure_id   = id,
                filepath    = fp,
                plot_type   = plot_type,
                section     = section,
                title       = title,
                description = desc,
                x_axis      = x_axis,
                y_axis      = y_axis,
                contrast    = NA_character_,
                stringsAsFactors = FALSE
            )
        }
    }

    # ---- QC: PCA ----
    if (dir.exists(diag_dir)) {
        pca_files <- list.files(diag_dir, pattern = "^PCA_.*\\.png$", full.names = TRUE)
        for (pf in pca_files) {
            bn <- gsub("\\.png$", "", basename(pf))
            fid <- tolower(gsub("[^a-zA-Z0-9]", "_", bn))
            figures[[fid]] <- data.frame(
                figure_id = fid, filepath = pf, plot_type = "scatter",
                section = "QC", title = gsub("\\.", " ", bn),
                description = "PCA of metabolomics samples",
                x_axis = "PC1", y_axis = "PC2", contrast = NA_character_,
                stringsAsFactors = FALSE
            )
        }

        # ---- QC: Density ----
        dens_files <- list.files(diag_dir, pattern = "^intensity_density.*\\.png$", full.names = TRUE)
        for (df in dens_files) {
            bn <- gsub("\\.png$", "", basename(df))
            fid <- gsub("intensity_", "", bn)
            is_raw <- grepl("raw", bn)
            figures[[fid]] <- data.frame(
                figure_id = fid, filepath = df, plot_type = "density",
                section = "QC",
                title = paste(if (is_raw) "Raw" else "Normalized", "Intensity Density"),
                description = paste("Per-sample density of", if (is_raw) "raw" else "normalized", "metabolite intensities"),
                x_axis = "Intensity", y_axis = "Density", contrast = NA_character_,
                stringsAsFactors = FALSE
            )
        }

        # ---- QC: Boxplots ----
        box_files <- list.files(diag_dir, pattern = "^intensity_boxplot.*\\.png$", full.names = TRUE)
        for (bf in box_files) {
            bn <- gsub("\\.png$", "", basename(bf))
            fid <- gsub("intensity_", "", bn)
            is_raw <- grepl("raw", bn)
            figures[[fid]] <- data.frame(
                figure_id = fid, filepath = bf, plot_type = "boxplot",
                section = "QC",
                title = paste(if (is_raw) "Raw" else "Normalized", "Intensity Boxplot"),
                description = paste("Per-sample boxplot of", if (is_raw) "raw" else "normalized", "metabolite intensities"),
                x_axis = "Sample", y_axis = "Intensity", contrast = NA_character_,
                stringsAsFactors = FALSE
            )
        }

        # ---- QC: Heatmaps ----
        hm_files <- list.files(diag_dir, pattern = "^sample_(distance|correlation)_heatmap.*\\.png$", full.names = TRUE)
        for (hf in hm_files) {
            bn <- gsub("\\.png$", "", basename(hf))
            fid <- bn
            figures[[fid]] <- data.frame(
                figure_id = fid, filepath = hf, plot_type = "heatmap",
                section = "QC", title = gsub("_", " ", tools::toTitleCase(bn)),
                description = "Sample-level distance or correlation heatmap",
                x_axis = "Samples", y_axis = "Samples", contrast = NA_character_,
                stringsAsFactors = FALSE
            )
        }

        # ---- DE: Volcano, MA, P-value hist ----
        volcano_files <- list.files(diag_dir, pattern = "^volcano_.*\\.png$", full.names = TRUE)
        for (vf in volcano_files) {
            cn <- gsub("^volcano_|\\.png$", "", basename(vf))
            fid <- paste0("volcano_", cn)
            figures[[fid]] <- data.frame(
                figure_id = fid, filepath = vf, plot_type = "volcano",
                section = "Differential Expression",
                title = paste("Volcano Plot:", gsub("_", " ", cn)),
                description = "Statistical significance vs fold change for all metabolites",
                x_axis = "log2 Fold Change", y_axis = "-log10(p-value)",
                contrast = cn, stringsAsFactors = FALSE
            )
        }

        ma_files <- list.files(diag_dir, pattern = "^MA_.*\\.png$", full.names = TRUE)
        for (mf in ma_files) {
            cn <- gsub("^MA_|\\.png$", "", basename(mf))
            fid <- paste0("ma_", cn)
            figures[[fid]] <- data.frame(
                figure_id = fid, filepath = mf, plot_type = "ma_plot",
                section = "Differential Expression",
                title = paste("MA Plot:", gsub("_", " ", cn)),
                description = "Mean-difference plot for all metabolites",
                x_axis = "Average Expression", y_axis = "log2 Fold Change",
                contrast = cn, stringsAsFactors = FALSE
            )
        }

        # ---- Feature Selection: RF, PLS-DA ----
        add_if_exists("rf_importance", "rf_importance.png", diag_dir,
                      "bar", "Feature Selection", "Random Forest Importance",
                      "Top features by Random Forest variable importance",
                      "Importance", "Feature")

        add_if_exists("plsda_scores", "plsda_scores.png", diag_dir,
                      "scatter", "Feature Selection", "PLS-DA Scores",
                      "PLS-DA scores plot showing sample separation",
                      "Component 1", "Component 2")

        add_if_exists("plsda_vip", "plsda_vip.png", diag_dir,
                      "bar", "Feature Selection", "PLS-DA VIP Scores",
                      "Variable Importance in Projection scores from PLS-DA",
                      "VIP Score", "Feature")
    }

    if (length(figures) == 0) return(data.frame())
    do.call(rbind, figures)
}


#' Build enriched data for metabolomics commentary context
#'
#' @param de_res      DE results object
#' @param qc_pre_obj  QC pre results
#' @param config      Full pipeline config
#' @return Named list with domain info, top metabolites, PCA variance, etc.
build_metabolomics_extra_data <- function(de_res, qc_pre_obj = NULL, config) {
    metab_cfg <- config$modes$metabolomics
    extra <- list()

    extra$domain <- "metabolomics"

    # --- Experimental groups ---
    tryCatch({
        group_col <- metab_cfg$de$condition_column %||%
                     metab_cfg$effects$color %||% "condition"
        extra$group_col <- group_col
        if (!is.null(de_res$summary_df)) {
            contrasts <- metab_cfg$de$contrasts
            if (is.list(contrasts)) contrasts <- unlist(contrasts)
            extra$contrasts <- contrasts
        }
    }, error = function(e) message("  Could not extract group info: ", e$message))

    # --- PCA variance ---
    tryCatch({
        if (!is.null(qc_pre_obj)) {
            var_expl <- qc_pre_obj$objects$var_expl
            if (!is.null(var_expl) && length(var_expl) >= 2) {
                extra$pca_variance <- as.numeric(var_expl[1:min(3, length(var_expl))])
            }
        }
    }, error = function(e) message("  Could not extract PCA variance: ", e$message))

    # --- Methods ---
    norm_info <- metab_cfg$preprocessing %||% list()
    norm_parts <- character(0)
    if (!is.null(norm_info$sample_norm) && tolower(norm_info$sample_norm) != "none")
        norm_parts <- c(norm_parts, norm_info$sample_norm)
    if (!is.null(norm_info$transform) && tolower(norm_info$transform) != "none")
        norm_parts <- c(norm_parts, norm_info$transform)
    if (!is.null(norm_info$scaling) && tolower(norm_info$scaling) != "none")
        norm_parts <- c(norm_parts, norm_info$scaling)
    extra$norm_method <- if (length(norm_parts) > 0) paste(norm_parts, collapse = " + ") else "none"

    extra$de_method <- metab_cfg$de$method %||% "limma"
    extra$de_cutoffs <- list(
        p_cutoff  = metab_cfg$de$p_cutoff %||% 0.05,
        fc_cutoff = metab_cfg$de$linear_fc_cutoff %||% 1.5
    )

    # --- Top DE metabolites per contrast ---
    extra$top_genes <- list()
    tryCatch({
        de_tables <- de_res$de_tables %||% de_res$tables
        if (!is.null(de_tables)) {
            for (cn in names(de_tables)) {
                tbl <- de_tables[[cn]]
                if (!is.data.frame(tbl)) next

                pval_col <- intersect(c("P.Value", "pvalue"), names(tbl))[1]
                fc_col <- intersect(c("logFC", "log2FC", "log2FoldChange"), names(tbl))[1]
                if (is.na(pval_col) || is.na(fc_col)) next

                p_cut <- extra$de_cutoffs$p_cutoff
                fc_cut <- log2(extra$de_cutoffs$fc_cutoff)
                sig <- !is.na(tbl[[pval_col]]) & tbl[[pval_col]] <= p_cut &
                       abs(tbl[[fc_col]]) >= fc_cut

                if (sum(sig) == 0) next

                sig_df <- data.frame(
                    gene   = tbl$feature_id[sig],
                    log2FC = tbl[[fc_col]][sig],
                    padj   = tbl[[pval_col]][sig],
                    stringsAsFactors = FALSE
                )

                up_df <- sig_df[sig_df$log2FC > 0, ]
                dn_df <- sig_df[sig_df$log2FC < 0, ]
                up_df <- up_df[order(-abs(up_df$log2FC)), ]
                dn_df <- dn_df[order(-abs(dn_df$log2FC)), ]
                top_combined <- rbind(
                    if (nrow(up_df) > 0) cbind(head(up_df, 5), direction = "up") else NULL,
                    if (nrow(dn_df) > 0) cbind(head(dn_df, 5), direction = "down") else NULL
                )

                if (!is.null(top_combined) && nrow(top_combined) > 0) {
                    extra$top_genes[[cn]] <- top_combined
                    message("  Top DE metabolites for ", cn, ": ",
                            min(nrow(up_df), 5), " up + ", min(nrow(dn_df), 5), " down")
                }
            }
        }
    }, error = function(e) message("  Could not extract top DE metabolites: ", e$message))

    extra
}
