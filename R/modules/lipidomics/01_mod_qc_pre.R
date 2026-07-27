# R/modules/lipidomics/01_mod_qc_pre.R
#
# Lipidomics pre/post-normalization QC module.
# Mirrors the metabolomics QC module, reusing core plotting functions.


#' Lipidomics pre-DE QC module
#'
#' @param pre     List from preprocess_lipidomics().
#' @param config  Full config.
#' @param out_dir Output directory for this mode.
#' @return list(files, plots, objects, missingness, normalization_eval)
mod_lipidomics_qc_pre <- function(pre, config, out_dir) {
    stage <- "lipidomics"
    stopifnot(is.character(out_dir), length(out_dir) == 1)
    assert_pre_contract(pre, stage = stage)

    dirs <- create_legacy_output_dirs(out_dir)
    out_qc <- dirs$diagnostic_plots
    out_ds <- dirs$datasets
    cfg <- config$modes$lipidomics

    plots <- list()
    files <- character(0)

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

    norm_label <- build_norm_label(pre$info$normalization)

    # ---- Build sample subsets ----
    sample_col <- cfg$effects$samples %||% "sample_id"
    subsets <- build_qc_subsets(expr_qc, pre$expr_filt, pre$meta, sample_col,
                                 expr_log = pre$expr_log)

    # ---- PCA plots ----
    color_config <- cfg$effects$color
    if (is.null(color_config)) color_config <- list("sample_type")
    if (!is.list(color_config)) color_config <- as.list(color_config)
    n_colors <- length(color_config)
    is_multi <- n_colors > 1

    for (s in subsets) {
        tag   <- s$tag
        for (ci in seq_along(color_config)) {
            color_var <- as.character(color_config[[ci]])
            if (!color_var %in% colnames(s$meta)) next

            cfg_temp <- cfg
            cfg_temp$effects$color <- color_var
            color_suffix <- if (is_multi) paste0("_by_", color_var) else ""

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
        }
    }

    # ---- Density / distribution plots ----
    cfg_primary <- cfg
    cfg_primary$effects$color <- as.character(color_config[[1]])

    for (s in subsets) {
        tag   <- s$tag
        label <- s$label

        # Density: raw (log2)
        f_dens_raw <- file.path(out_qc, paste0("intensity_density_raw", tag, ".png"))
        expr_filt_log2 <- log2(s$expr_filt + 1)
        p_dens_raw <- qc_omic_density(
            expr_filt_log2, s$meta, cfg_primary,
            out_file = f_dens_raw,
            title = paste0("Density: log2(raw intensities)", label)
        )
        files <- c(files, f_dens_raw)
        plots[[paste0("density_raw", tag)]] <- p_dens_raw

        # Density: transformed (log2 + normalization, no scaling)
        if (!is.null(s$expr_log)) {
            f_dens_trans <- file.path(out_qc, paste0("intensity_density_transformed", tag, ".png"))
            p_dens_trans <- qc_omic_density(
                s$expr_log, s$meta, cfg_primary,
                out_file = f_dens_trans,
                title = paste0("Density: log2 + normalized", label)
            )
            files <- c(files, f_dens_trans)
            plots[[paste0("density_transformed", tag)]] <- p_dens_trans
        }

        # Density: fully normalized (log2 + normalization + scaling)
        f_dens_norm <- file.path(out_qc, paste0("intensity_density_normalized", tag, ".png"))
        p_dens_norm <- qc_omic_density(
            s$expr_work, s$meta, cfg_primary,
            out_file = f_dens_norm,
            title = paste0("Density: ", norm_label, label)
        )
        files <- c(files, f_dens_norm)
        plots[[paste0("density_norm", tag)]] <- p_dens_norm

        # Boxplot: raw (log2)
        f_box_raw <- file.path(out_qc, paste0("intensity_boxplot_raw", tag, ".png"))
        p_box_raw <- norm_boxplot(
            expr_filt_log2, s$meta, cfg_primary,
            out_file = f_box_raw,
            title = paste0("Boxplot: log2(raw intensities)", label),
            y_label = "log2(raw intensity)"
        )
        files <- c(files, f_box_raw)
        plots[[paste0("boxplot_raw", tag)]] <- p_box_raw

        # Boxplot: transformed (log2 + normalization, no scaling)
        if (!is.null(s$expr_log)) {
            f_box_trans <- file.path(out_qc, paste0("intensity_boxplot_transformed", tag, ".png"))
            p_box_trans <- norm_boxplot(
                s$expr_log, s$meta, cfg_primary,
                out_file = f_box_trans,
                title = paste0("Boxplot: log2 + normalized", label),
                y_label = "log2(normalized intensity)"
            )
            files <- c(files, f_box_trans)
            plots[[paste0("boxplot_transformed", tag)]] <- p_box_trans
        }

        # Boxplot: fully normalized (log2 + normalization + scaling)
        f_box <- file.path(out_qc, paste0("intensity_boxplot_normalized", tag, ".png"))
        p_box <- norm_boxplot(
            s$expr_work, s$meta, cfg_primary,
            out_file = f_box,
            title = paste0("Boxplot: ", norm_label, label),
            y_label = "Normalized intensity"
        )
        files <- c(files, f_box)
        plots[[paste0("boxplot", tag)]] <- p_box
    }

    # ---- Heatmaps ----
    n_samples <- ncol(expr_qc)
    max_hm <- cfg$qc$thresholds$max_samples_for_heatmaps %||% 120

    if (n_samples <= max_hm) {
        annot_cols <- build_heatmap_annotation(pre$meta, cfg)

        for (s in subsets) {
            tag <- s$tag

            f_dist <- file.path(out_qc, paste0("sample_distance_heatmap", tag, ".png"))
            qc_sample_distance_heatmap(s$expr_work, s$meta, cfg_primary,
                                        out_file = f_dist, annot_cols = annot_cols)
            files <- c(files, f_dist)

            f_cor <- file.path(out_qc, paste0("sample_correlation_heatmap", tag, ".png"))
            qc_sample_correlation_heatmap(s$expr_work, s$meta, cfg_primary,
                                           out_file = f_cor, annot_cols = annot_cols)
            files <- c(files, f_cor)
        }
    }

    # ---- Extract PCA objects ----
    pca_obj <- NULL
    scores <- NULL
    var_expl <- NULL

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
            pca_scores = assert_pca_scores(scores, context = "lipidomics QC"),
            var_expl  = var_expl,
            color     = cfg$effects$color,
            shape     = cfg$effects$shape %||% NULL
        ),
        missingness       = miss_info,
        normalization_eval = pre$normalization_eval
    )
}
