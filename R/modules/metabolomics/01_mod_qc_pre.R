# R/modules/metabolomics/01_mod_qc_pre.R
#
# Metabolomics pre-normalization and post-normalization QC module.
# Produces diagnostics to support normalization method choice.
#
# Reuses core plotting: qc_pca_scatter, qc_pca_3d, norm_boxplot,
#   qc_omic_density (works for any intensity matrix),
#   qc_sample_distance_heatmap,
#   plot_density_overlay, save_heatmap_to_file, save_tsv


#' Metabolomics pre-DE QC module
#'
#' @param pre     List from preprocess_metabolomics().
#' @param config  Full config.
#' @param out_dir Output directory for this mode.
#' @return list(files, plots, objects, missingness, normalization_eval)
mod_metabolomics_qc_pre <- function(pre, config, out_dir) {
    stage <- "metabolomics"
    stopifnot(is.character(out_dir), length(out_dir) == 1)
    assert_pre_contract(pre, stage = stage)

    dirs <- create_legacy_output_dirs(out_dir)
    out_qc <- dirs$diagnostic_plots
    out_ds <- dirs$datasets
    cfg <- config$modes$metabolomics

    plots <- list()
    files <- character(0)

    # Determine which matrix to use for QC (expr_work = normalized)
    expr_qc <- pre$expr_work

    # ---- Missingness summary ----
    miss_info <- pre$info$missingness
    if (!is.null(miss_info)) {
        miss_df <- data.frame(
            sample_id = names(miss_info$per_sample),
            pct_missing = round(100 * miss_info$per_sample, 2),
            stringsAsFactors = FALSE
        )
        f_miss_sample <- save_tsv(miss_df, out_ds, "missingness_per_sample.tsv")
        files <- c(files, f_miss_sample)

        feat_miss_df <- data.frame(
            feature_id = names(miss_info$per_feature),
            pct_missing = round(100 * miss_info$per_feature, 2),
            stringsAsFactors = FALSE
        )
        f_miss_feat <- save_tsv(feat_miss_df, out_ds, "missingness_per_feature.tsv")
        files <- c(files, f_miss_feat)
    }

    # ---- Normalization evaluation table ----
    if (!is.null(pre$normalization_eval)) {
        f_eval <- save_tsv(pre$normalization_eval, out_ds,
                           "normalization_method_comparison.tsv")
        files <- c(files, f_eval)
    }

    # ---- Normalization info summary ----
    norm_info <- pre$info$normalization
    if (!is.null(norm_info)) {
        info_df <- data.frame(
            parameter = names(norm_info),
            value = as.character(unlist(norm_info)),
            stringsAsFactors = FALSE
        )
        f_norm_info <- save_tsv(info_df, out_ds, "normalization_applied.tsv")
        files <- c(files, f_norm_info)
    }

    # ---- Build normalization label for titles ----
    norm_label <- build_norm_label(pre$info$normalization)

    # ---- Build sample subsets: "all" and "noQC" ----
    sample_col <- cfg$effects$samples %||% "sample_id"
    subsets <- build_qc_subsets(expr_qc, pre$expr_filt, pre$meta, sample_col)

    # ---- PCA plots ----
    color_config <- cfg$effects$color
    if (is.null(color_config)) color_config <- list("sample_type")
    if (!is.list(color_config)) color_config <- as.list(color_config)
    n_colors <- length(color_config)
    is_multi <- n_colors > 1

    for (s in subsets[1]) {
        tag   <- s$tag
        label <- s$label

        for (ci in 1) {
            color_var <- as.character(color_config[[ci]])

            if (!color_var %in% colnames(s$meta)) {
                warning("PCA: color column '", color_var, "' not in metadata; skipping.")
                next
            }

            cfg_temp <- cfg
            cfg_temp$effects$color <- color_var

            color_suffix <- ""

            f_pca12 <- file.path(out_qc, paste0("PCA_PC1.vs.PC2", color_suffix, tag, ".png"))
            p12 <- qc_pca_scatter(s$expr_work, s$meta, cfg_temp, pcs = c(1, 2),
                                   out_file = f_pca12)
            files <- c(files, f_pca12)
            plots[[paste0("pca_1_2", color_suffix, tag)]] <- p12

            f_pca13 <- file.path(out_qc, paste0("PCA_PC1.vs.PC3", color_suffix, tag, ".png"))
            p13 <- qc_pca_scatter(s$expr_work, s$meta, cfg_temp, pcs = c(1, 3),
                                   out_file = f_pca13)
            files <- c(files, f_pca13)
            plots[[paste0("pca_1_3", color_suffix, tag)]] <- p13

            # 3D PCA only for the first color variable, all-samples subset
            if (ci == 1 && tag == "" && ncol(s$expr_work) >= 10) {
                f_pca3d <- file.path(out_qc, "PCA_3D.html")
                p_3d <- qc_pca_3d(s$expr_work, s$meta, cfg_temp, out_file = f_pca3d)
                files <- c(files, f_pca3d)
                if (!is.null(p_3d)) plots$pca_3d <- p_3d
            }
        }
    }

    # ---- Density / distribution plots ----
    cfg_primary <- cfg
    cfg_primary$effects$color <- as.character(color_config[[1]])

    for (s in subsets[1]) {
        tag   <- s$tag       # "" or "_noQC"
        label <- s$label     # "" or " [excl. QC]"

        # -- Density: raw (log2-transformed for visualization) --
        f_dens_raw <- file.path(out_qc, paste0("intensity_density_raw", tag, ".png"))
        expr_filt_log2 <- log2(s$expr_filt + 1)
        p_dens_raw <- qc_omic_density(
            expr_filt_log2, s$meta, cfg_primary,
            out_file = f_dens_raw,
            title = paste0("Density: log2(raw intensities) (before normalization)", label)
        )
        files <- c(files, f_dens_raw)
        plots[[paste0("density_raw", tag)]] <- p_dens_raw

        # -- Density: log2(raw) -- (retained for plots list; not saved to disk)
        expr_log2_raw <- log2(pmax(s$expr_filt, 1))
        p_dens_log2 <- qc_omic_density(
            expr_log2_raw, s$meta, cfg_primary,
            out_file = NULL,
            title = paste0("Density: log2(raw intensities)", label)
        )
        plots[[paste0("density_log2raw", tag)]] <- p_dens_log2

        # -- Density: normalized --
        f_dens_norm <- file.path(out_qc, paste0("intensity_density_normalized", tag, ".png"))
        p_dens_norm <- qc_omic_density(
            s$expr_work, s$meta, cfg_primary,
            out_file = f_dens_norm,
            title = paste0("Density: ", norm_label, label)
        )
        files <- c(files, f_dens_norm)
        plots[[paste0("density_norm", tag)]] <- p_dens_norm

        # -- Boxplot: raw (log2-transformed for visualization) --
        f_box_raw <- file.path(out_qc, paste0("intensity_boxplot_raw", tag, ".png"))
        expr_filt_log2_box <- log2(s$expr_filt + 1)
        p_box_raw <- norm_boxplot(
            expr_filt_log2_box, s$meta, cfg_primary,
            out_file = f_box_raw,
            title = paste0("Boxplot: log2(raw intensities)", label)
        )
        files <- c(files, f_box_raw)
        plots[[paste0("boxplot_raw", tag)]] <- p_box_raw

        # -- Boxplot: normalized --
        f_box <- file.path(out_qc, paste0("intensity_boxplot_normalized", tag, ".png"))
        p_box <- norm_boxplot(
            s$expr_work, s$meta, cfg_primary,
            out_file = f_box,
            title = paste0("Boxplot: ", norm_label, label)
        )
        files <- c(files, f_box)
        plots[[paste0("boxplot", tag)]] <- p_box
    }

    # ---- Heatmaps (sample-level) ----
    n_samples <- ncol(expr_qc)
    max_hm <- cfg$qc$thresholds$max_samples_for_heatmaps %||% 120

    if (n_samples <= max_hm) {
        annot_cols <- build_heatmap_annotation(pre$meta, cfg)

        # Distance heatmap: all samples (with QC)
        s_all <- subsets[[1]]
        f_dist_wqc <- file.path(out_qc, "sample_distance_metab_heatmap_w_qc.png")
        ph_dist_wqc <- qc_sample_distance_heatmap(s_all$expr_work, s_all$meta, cfg_primary,
                                    out_file = f_dist_wqc, annot_cols = annot_cols)
        files <- c(files, f_dist_wqc)
        plots[["dist_heatmap"]] <- ph_dist_wqc

        # Distance heatmap: without QC (only when QC samples exist)
        if (length(subsets) >= 2) {
            s_noqc <- subsets[[2]]
            f_dist <- file.path(out_qc, "sample_distance_metab_heatmap.png")
            ph_dist <- qc_sample_distance_heatmap(s_noqc$expr_work, s_noqc$meta, cfg_primary,
                                        out_file = f_dist, annot_cols = annot_cols)
            files <- c(files, f_dist)
            plots[["dist_heatmap_noQC"]] <- ph_dist
        }

        # Expression heatmap: all samples
        f_hm <- file.path(out_qc, "samples_metab_heatmap.png")
        hm_clusters <- wrap_qc_heatmap(s_all$expr_work, s_all$meta, cfg_primary,
                                       stage = stage, out_file = f_hm)
        files <- c(files, f_hm)
        plots[["heatmap_clusters"]] <- hm_clusters

        # Expression heatmap: without column clustering
        f_hm_nocol <- file.path(out_qc, "samples_metab_heatmap_wo_col.png")
        hm_nocol <- wrap_qc_heatmap(s_all$expr_work, s_all$meta, cfg_primary,
                                    stage = stage, out_file = f_hm_nocol,
                                    cluster_cols = FALSE)
        files <- c(files, f_hm_nocol)
        plots[["heatmap_nocol"]] <- hm_nocol

        # Correlation heatmap: without QC
        if (length(subsets) >= 2) {
            s_noqc <- subsets[[2]]
            f_cor <- file.path(out_qc, "sample_correlation_metab_heatmap.png")
            p_cor <- qc_sample_correlation_heatmap(s_noqc$expr_work, s_noqc$meta,
                                                   cfg_primary, out_file = f_cor,
                                                   annot_cols = annot_cols,
                                                   adjust_scale = TRUE)
            files <- c(files, f_cor)
            plots[["correlation_noQC"]] <- p_cor
        }
    }

    # ---- Extract PCA objects for downstream use ----
    pca_obj <- NULL
    scores <- NULL
    var_expl <- NULL

    # Prefer the all-samples PCA (no _noQC tag) for downstream objects
    pca_keys <- grep("^pca_1_2", names(plots), value = TRUE)
    first_pca_key <- pca_keys[!grepl("_noQC$", pca_keys)][1]
    if (is.na(first_pca_key)) first_pca_key <- pca_keys[1]
    if (!is.na(first_pca_key)) {
        p_main <- plots[[first_pca_key]]
        pca_obj  <- attr(p_main, "pca_result")
        scores   <- attr(p_main, "scores")
        var_expl <- attr(p_main, "var_expl")
    }

    list(
        plots  = plots,
        files  = unique(files),
        objects = list(
            norm_log_counts_pca = pca_obj,
            pca_scores = assert_pca_scores(scores, context = "metabolomics QC"),
            var_expl  = var_expl,
            color     = cfg$effects$color,
            shape     = cfg$effects$shape %||% NULL
        ),
        missingness       = miss_info,
        normalization_eval = pre$normalization_eval
    )
}


# ---- helpers ----------------------------------------------------------------

#' Build heatmap annotation column names from config
#'
#' Uses cfg$effects$heatmap_annotations if set; otherwise falls back to
#' the primary color variable.
#' @return Character vector of column names, or NULL (which means use default).
build_heatmap_annotation <- function(meta, cfg) {
    annot_cols <- cfg$effects$heatmap_annotations
    if (!is.null(annot_cols)) {
        annot_cols <- unlist(annot_cols)
        # Filter to columns that actually exist
        annot_cols <- intersect(annot_cols, colnames(meta))
        if (length(annot_cols) > 0) return(annot_cols)
    }
    NULL
}


#' Build a human-readable normalization label from the applied info
#'
#' @param norm_applied List with sample_norm, transform, scaling.
#' @return Character string, e.g. "PQN + log2 + pareto".
build_norm_label <- function(norm_applied) {
    if (is.null(norm_applied)) return("normalized intensities")
    parts <- character(0)
    sn <- norm_applied$sample_norm
    if (!is.null(sn) && tolower(sn) != "none") parts <- c(parts, toupper(sn))
    tr <- norm_applied$transform
    if (!is.null(tr) && tolower(tr) != "none") parts <- c(parts, tr)
    sc <- norm_applied$scaling
    if (!is.null(sc) && tolower(sc) != "none") parts <- c(parts, sc)
    if (length(parts) == 0) return("intensities (no normalization)")
    paste(parts, collapse = " + ")
}


#' Build sample subsets for QC: all samples and (if QCs exist) without QC
#'
#' QC samples are identified by treatment == "QC" (case-insensitive) in metadata.
#' Returns a list of subset descriptors, each with: tag, label, expr_work,
#' expr_filt, meta.
#'
#' @param expr_work  Normalized expression matrix.
#' @param expr_filt  Filtered (raw) expression matrix.
#' @param meta       Metadata data.frame.
#' @param sample_col Column name for sample IDs.
#' @return list of subset descriptors.
build_qc_subsets <- function(expr_work, expr_filt, meta, sample_col) {
    all_subset <- list(
        tag       = "",
        label     = "",
        expr_work = expr_work,
        expr_filt = expr_filt,
        meta      = meta
    )

    subsets <- list(all_subset)

    # Identify QC samples via the treatment column
    if ("treatment" %in% colnames(meta)) {
        is_qc <- tolower(trimws(as.character(meta$treatment))) == "qc"
        if (any(is_qc) && !all(is_qc)) {
            keep_ids <- as.character(meta[[sample_col]][!is_qc])
            subsets[[2]] <- list(
                tag       = "_noQC",
                label     = " [excl. QC]",
                expr_work = expr_work[, keep_ids, drop = FALSE],
                expr_filt = expr_filt[, keep_ids, drop = FALSE],
                meta      = meta[!is_qc, , drop = FALSE]
            )
        }
    }

    subsets
}
