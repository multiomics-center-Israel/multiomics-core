#' RNA-seq pre-DE QC module
#'
#' QC plots that depend only on preprocessing outputs.
#' Mirrors Proteomics logic using core plotting functions.
#'
#' @param pre List returned by preprocess_rna()
#' @param config Full config list
#' @param out_dir Output root directory
#' @return list(plots, files, objects, ...)
mod_rnaseq_qc_pre <- function(pre, config, out_dir) {
    stage <- "rna"
    assert_pre_contract(pre, stage = stage)

    dirs <- create_legacy_output_dirs(out_dir, create = TRUE)
    out_qc <- dirs$diagnostic_plots
    cfg <- config$modes$rna

    plots <- list()
    files <- character(0)

    # Use working expression matrix (normalized log-counts) for QC
    mat <- pre$expr_work
    meta <- pre$meta

    # --- Validate color/shape configuration ---
    validate_effects_config(cfg$effects)

    # --- Determine if adaptive plotting is enabled ---
    adaptive_enabled <- if (!is.null(cfg$qc) && !is.null(cfg$qc$adaptive_plots)) {
        isTRUE(cfg$qc$adaptive_plots)
    } else {
        TRUE  # Default: adaptive plotting ON
    }

    # --- Get thresholds (use defaults if not configured) ---
    thresh <- if (!is.null(cfg$qc)) cfg$qc$thresholds else NULL
    thresh <- thresh %||% list()
    min_samples_pca3d <- thresh$min_samples_for_pca3d %||% 10
    max_samples_heatmaps <- thresh$max_samples_for_heatmaps %||% 120
    max_samples_expr_heatmap <- thresh$max_samples_for_expr_heatmap %||% 60

    n_samples <- nrow(meta)

    # --- Parse color configuration (string OR array) ---
    color_config <- cfg$effects$color
    if (is.null(color_config)) {
        stop("Config error: effects$color is required but missing")
    }

    # Check if color is multi-value (more defensive)
    n_colors <- length(color_config)
    is_multi <- n_colors > 1

    color_values <- if (isTRUE(is_multi)) {
        as.character(color_config)
    } else {
        as.character(color_config[1])
    }
    primary_color <- color_values[1]

    # ---------- PCA plots (one per color value) ----------
    # PCA 1v2 (always generated for each color)
    for (i in seq_along(color_values)) {
        color_var <- color_values[i]

        # Create temporary config with this color
        cfg_temp <- cfg
        cfg_temp$effects$color <- color_var

        # Determine filename
        if (i == 1) {
            f_pca12 <- file.path(out_qc, "PCA_PC1.vs.PC2.png")
        } else {
            safe_name <- gsub("[^A-Za-z0-9_]", "_", color_var)
            f_pca12 <- file.path(out_qc, sprintf("PCA_PC1.vs.PC2_by_%s.png", safe_name))
        }

        p12 <- qc_pca_scatter(mat, meta, cfg_temp, pcs = c(1, 2), out_file = f_pca12)
        files <- c(files, f_pca12)

        # Store first plot for object extraction later
        if (i == 1) {
            plots$pca_1_2 <- p12
            primary_pca <- p12
        } else {
            plots[[sprintf("pca_1_2_by_%s", color_var)]] <- p12
        }
    }

    # PCA 1v3 (always generated, all colors)
    for (i in seq_along(color_values)) {
        color_var <- color_values[i]
        cfg_temp <- cfg
        cfg_temp$effects$color <- color_var

        if (i == 1) {
            f_pca13 <- file.path(out_qc, "PCA_PC1.vs.PC3.png")
        } else {
            safe_name <- gsub("[^A-Za-z0-9_]", "_", color_var)
            f_pca13 <- file.path(out_qc, sprintf("PCA_PC1.vs.PC3_by_%s.png", safe_name))
        }

        p13 <- qc_pca_scatter(mat, meta, cfg_temp, pcs = c(1, 3), out_file = f_pca13)
        files <- c(files, f_pca13)

        if (i == 1) {
            plots$pca_1_3 <- p13
        } else {
            plots[[sprintf("pca_1_3_by_%s", color_var)]] <- p13
        }
    }

    # PCA 3D (conditional: requires minimum sample count, primary color only)
    should_generate_pca3d <- !isTRUE(adaptive_enabled) || (n_samples >= min_samples_pca3d)
    if (should_generate_pca3d) {
        f_pca3d <- file.path(out_qc, "PCA_3D.html")
        cfg_temp <- cfg
        cfg_temp$effects$color <- primary_color
        p_3d <- qc_pca_3d(mat, meta, cfg_temp, out_file = f_pca3d)
        files <- c(files, f_pca3d)
        if (!is.null(p_3d)) plots$pca_3d <- p_3d
    } else if (isTRUE(adaptive_enabled)) {
        message(sprintf("Skipping 3D PCA (%d samples < %d threshold)",
                        n_samples, min_samples_pca3d))
    }

    # ---------- Density (primary color only) ----------
    cfg_temp <- cfg
    cfg_temp$effects$color <- primary_color
    f_hist <- file.path(out_qc, "rna_histograms_summary.png")
    p_dens <- qc_omic_density(
        mat, meta, cfg_temp,
        out_file = f_hist,
        title = "Density plot of normalized counts"
    )
    files <- c(files, f_hist)
    plots$density <- p_dens

    # ---------- Boxplots (primary color only) ----------
    f_box <- file.path(out_qc, "norm_boxplot.png")
    p_bp <- norm_boxplot(mat, meta, cfg_temp, out_file = f_box)
    files <- c(files, f_box)
    plots$boxplot <- p_bp

    # --- Get heatmap annotation columns ---
    annot_cols <- cfg$effects$heatmap_annotations %||% NULL
    # Fallback to primary color if no annotations specified
    if (is.null(annot_cols)) {
        annot_cols <- NULL  # Will use default single-column in heatmap functions
    }

    # --- Get heatmap visualization settings ---
    hm_viz <- cfg$qc$heatmap_viz %||% list()
    show_labels <- hm_viz$show_sample_labels %||% TRUE  # Default TRUE to show sample names
    cluster_samples <- hm_viz$cluster_samples %||% TRUE
    adjust_scale <- hm_viz$adjust_correlation_scale %||% TRUE

    # ---------- Sample distance heatmap (conditional: skip for large datasets) ----------
    should_generate_dist_heatmap <- !isTRUE(adaptive_enabled) || (n_samples <= max_samples_heatmaps)
    if (should_generate_dist_heatmap) {
        f_dist <- file.path(out_qc, "sample_distance_heatmap.png")
        # Use primary color for heatmap
        cfg_temp <- cfg
        cfg_temp$effects$color <- primary_color
        ph_dist <- qc_sample_distance_heatmap(
            mat, meta, cfg_temp,
            annot_cols = annot_cols,
            out_file = f_dist,
            show_labels = show_labels,
            cluster_samples = cluster_samples
        )
        files <- c(files, f_dist)
        plots$dist_heatmap <- ph_dist
    } else if (isTRUE(adaptive_enabled)) {
        message(sprintf("Skipping distance heatmap (%d samples > %d threshold)",
                        n_samples, max_samples_heatmaps))
    }

    # ---------- Correlation heatmap (conditional: skip for large datasets) ----------
    should_generate_corr_heatmap <- !isTRUE(adaptive_enabled) || (n_samples <= max_samples_heatmaps)
    if (should_generate_corr_heatmap) {
        f_cor <- file.path(out_qc, "sample_correlation_heatmap.png")
        # Use primary color for heatmap
        cfg_temp <- cfg
        cfg_temp$effects$color <- primary_color
        p_cor <- qc_sample_correlation_heatmap(
            mat, meta, cfg_temp,
            annot_cols = annot_cols,
            out_file = f_cor,
            show_labels = show_labels,
            cluster_samples = cluster_samples,
            adjust_scale = adjust_scale
        )
        files <- c(files, f_cor)
        plots$correlation <- p_cor
    } else if (isTRUE(adaptive_enabled)) {
        message(sprintf("Skipping correlation heatmap (%d samples > %d threshold)",
                        n_samples, max_samples_heatmaps))
    }

    # ---------- Expression heatmaps (conditional: skip for large datasets) ----------
    should_generate_expr_heatmap <- !isTRUE(adaptive_enabled) || (n_samples <= max_samples_expr_heatmap)
    if (should_generate_expr_heatmap) {
        # Prepare annotations
        if (!is.null(annot_cols)) {
            # Use configured annotation columns
            annot <- meta[, annot_cols, drop = FALSE]
            rownames(annot) <- meta[[cfg$effects$samples]]
        } else {
            # Fallback to primary color
            annot <- data.frame(
                Condition = meta[[primary_color]],
                row.names = meta[[cfg$effects$samples]]
            )
        }

        # Heatmap with column clustering
        f_hm <- file.path(out_qc, "samples_rna_heatmap.png")
        # Use primary color for heatmap
        cfg_temp <- cfg
        cfg_temp$effects$color <- primary_color
        hm_clusters <- wrap_qc_heatmap(mat, meta, cfg_temp, stage = stage, out_file = f_hm)
        files <- c(files, f_hm)
        plots$heatmap_clusters <- hm_clusters

        # Heatmap without column clustering
        f_hm_nocol <- file.path(out_qc, "samples_rna_heatmap_wo_col.png")
        p_hm_nocol <- plot_heatmap_core(
            expr_mat = mat,
            annotation_col = annot,
            title = "QC: Sample RNA Expression",
            max_rows = 2000,
            cluster_rows = TRUE,
            cluster_cols = FALSE
        )
        save_heatmap_to_file(p_hm_nocol, f_hm_nocol)
        files <- c(files, f_hm_nocol)
        plots$heatmap_nocol <- p_hm_nocol
    } else if (isTRUE(adaptive_enabled)) {
        message(sprintf("Skipping expression heatmap (%d samples > %d threshold)",
                        n_samples, max_samples_expr_heatmap))
    }

    # ---------- Extract PCA objects from plot attributes (primary PCA plot) ----------
    pca_obj <- attr(primary_pca, "pca_result")
    scores <- attr(primary_pca, "scores")
    var_expl <- attr(primary_pca, "var_expl")

    eff_color <- cfg$effects$color %||% NULL
    eff_shape <- cfg$effects$shape %||% NULL

    # objects list (Proteomics style)
    objs <- list(
        norm_log_counts_pca = pca_obj,
        pca_scores = assert_pca_scores(scores, context = "rnaseq QC"),
        var_expl = var_expl,
        color = eff_color,
        shape = eff_shape
    )

    list(
        files = unique(files),
        plots = plots,
        objects = objs,
        # Keep these at top level for compatibility with current build_data_to_shiny_legacy_rna
        # (until that function is updated to look inside 'objects')
        norm_log_counts_pca = pca_obj,
        pca_scores = assert_pca_scores(scores, context = "rnaseq QC")
    )
}

#' Validate effects configuration for color/shape conflicts
#'
#' Ensures that:
#' 1. Both color and shape cannot be arrays/lists
#' 2. Color and shape cannot overlap (same column in both)
#'
#' @param effects The effects configuration list
validate_effects_config <- function(effects) {
    if (is.null(effects)) return(invisible(NULL))

    color <- effects$color
    shape <- effects$shape

    # Check if both are arrays/lists (more defensive)
    n_color <- if (!is.null(color)) length(color) else 0
    n_shape <- if (!is.null(shape)) length(shape) else 0

    color_is_multi <- n_color > 1
    shape_is_multi <- n_shape > 1

    # Use explicit TRUE checks to ensure scalar logicals
    if (isTRUE(color_is_multi) && isTRUE(shape_is_multi)) {
        stop(
            "Invalid effects configuration: Both 'color' and 'shape' cannot be arrays/lists.\n",
            "Please specify only one as an array for multiple plots.\n",
            "Current config:\n",
            "  color: ", paste(color, collapse = ", "), "\n",
            "  shape: ", paste(shape, collapse = ", ")
        )
    }

    # Check for overlap
    if (!is.null(color) && !is.null(shape)) {
        color_values <- if (is.list(color)) unlist(color) else as.character(color)
        shape_values <- if (is.list(shape)) unlist(shape) else as.character(shape)

        overlap <- intersect(color_values, shape_values)
        n_overlap <- length(overlap)
        if (n_overlap > 0) {
            stop(
                "Invalid effects configuration: 'color' and 'shape' cannot overlap.\n",
                "The following column(s) appear in both: ", paste(overlap, collapse = ", "), "\n",
                "Current config:\n",
                "  color: ", paste(color_values, collapse = ", "), "\n",
                "  shape: ", paste(shape_values, collapse = ", ")
            )
        }
    }

    invisible(NULL)
}
