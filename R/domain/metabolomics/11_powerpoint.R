# R/domain/metabolomics/11_powerpoint.R
#
# Generate a PowerPoint presentation summarizing metabolomics pipeline results.
# Each slide includes an AI-generated bottom line via Claude Code CLI,
# following the same pattern as 10_pipeline_summary.R and 12_commentary.R.
#
# Requires: officer, flextable (optional, for tables)
# Depends:  R/core/16_pptx_helpers.R (generate_slide_bottom_line, pptx_enabled)


# ==============================================================================
# TEXT FORMATTING HELPERS
# ==============================================================================

# Slide dimensions (inches) — standard widescreen 10x7.5
SLIDE_W <- 10
SLIDE_H <- 7.5

# Colors
COL_TITLE   <- "#1B2A4A"
COL_BODY    <- "#333333"
COL_ACCENT  <- "#2563EB"
COL_SUBTLE  <- "#6B7280"
COL_UP      <- "#16A34A"
COL_DOWN    <- "#DC2626"
COL_BL_BG   <- "#F0F4FF"

# Font properties
fp_slide_title <- function() {
    officer::fp_text(font.size = 20, bold = TRUE, color = COL_TITLE,
                      font.family = "Calibri")
}

fp_body <- function(sz = 11) {
    officer::fp_text(font.size = sz, color = COL_BODY, font.family = "Calibri")
}

fp_body_bold <- function(sz = 11) {
    officer::fp_text(font.size = sz, bold = TRUE, color = COL_BODY,
                      font.family = "Calibri")
}

fp_label <- function() {
    officer::fp_text(font.size = 9, color = COL_SUBTLE, font.family = "Calibri")
}

fp_stat_value <- function(color = COL_ACCENT) {
    officer::fp_text(font.size = 22, bold = TRUE, color = color,
                      font.family = "Calibri")
}

fp_stat_label <- function() {
    officer::fp_text(font.size = 9, color = COL_SUBTLE, font.family = "Calibri")
}

fp_bottom_line <- function() {
    officer::fp_text(font.size = 9, italic = TRUE, color = COL_TITLE,
                      font.family = "Calibri")
}

fp_bl_prefix <- function() {
    officer::fp_text(font.size = 9, bold = TRUE, italic = TRUE, color = COL_ACCENT,
                      font.family = "Calibri")
}

fp_footer <- function() {
    officer::fp_text(font.size = 7, color = COL_SUBTLE, font.family = "Calibri")
}


#' Add a figure slide with title, image, optional subtitle, and AI bottom line
#'
#' @param pptx    officer pptx object
#' @param title   Slide title
#' @param img     Path to PNG
#' @param subtitle Optional text below title (stats, etc.)
#' @param bl      AI bottom line text
#' @param footer  Footer text
#' @return Updated pptx object
add_figure_slide <- function(pptx, title, img, subtitle = NULL, bl = "",
                              footer = "") {
    pptx <- officer::add_slide(pptx, layout = "Blank", master = "Office Theme")

    # Title bar
    pptx <- officer::ph_with(pptx,
        value = officer::fpar(officer::ftext(title, fp_slide_title())),
        location = officer::ph_location(left = 0.5, top = 0.3, width = 9, height = 0.5))

    # Thin accent line under title
    pptx <- officer::ph_with(pptx,
        value = officer::block_list(officer::fpar(officer::ftext("", fp_label()))),
        location = officer::ph_location(left = 0.5, top = 0.75, width = 9, height = 0.02))

    # Subtitle (stats line)
    top_img <- 0.95
    if (!is.null(subtitle) && nzchar(subtitle)) {
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(officer::ftext(subtitle, fp_label())),
            location = officer::ph_location(left = 0.5, top = 0.85, width = 9, height = 0.3))
        top_img <- 1.15
    }

    # Image — centered, with good aspect ratio
    img_w <- 7.5
    img_h <- 4.3
    img_left <- (SLIDE_W - img_w) / 2
    if (file.exists(img)) {
        pptx <- officer::ph_with(pptx,
            value = officer::external_img(img, width = img_w, height = img_h),
            location = officer::ph_location(left = img_left, top = top_img,
                                             width = img_w, height = img_h))
    }

    # AI bottom line
    if (nzchar(bl)) {
        bl_parts <- officer::fpar(
            officer::ftext("AI Summary: ", fp_bl_prefix()),
            officer::ftext(bl, fp_bottom_line())
        )
        pptx <- officer::ph_with(pptx, value = bl_parts,
            location = officer::ph_location(left = 0.5, top = 5.7, width = 9, height = 1.1))
    }

    # Footer
    if (nzchar(footer)) {
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(officer::ftext(footer, fp_footer())),
            location = officer::ph_location(left = 0.3, top = 7.1, width = 9.4, height = 0.3))
    }

    pptx
}


# ==============================================================================
# POWERPOINT GENERATION
# ==============================================================================

#' Generate a metabolomics summary PowerPoint
#'
#' @param pre             Preprocessed data list.
#' @param qc_res          QC results list (with $plots and $files).
#' @param de_res          DE results list (with $summary_df, $de_tables).
#' @param feature_sel_res Feature selection results (rf, plsda) or NULL.
#' @param enrichment_res  Enrichment results list or NULL.
#' @param config          Full pipeline config.
#' @param out_dir         Output directory for this mode.
#' @return Character path to the generated .pptx file.
generate_metabolomics_pptx <- function(pre, qc_res, de_res, feature_sel_res,
                                        enrichment_res, config, out_dir) {

    # Check both metabolomics and lipidomics mode names
    mode_key <- if (!is.null(config$modes$metabolomics)) "metabolomics" else "lipidomics"
    if (!pptx_enabled(config, mode_key)) {
        message("PowerPoint generation disabled — skipping")
        return(character(0))
    }

    if (!requireNamespace("officer", quietly = TRUE)) {
        warning("officer package not available — skipping PowerPoint generation")
        return(character(0))
    }

    # Detect actual mode: lipidomics reuses this function
    mode_name <- basename(out_dir)
    is_lipidomics <- grepl("lipid", mode_name, ignore.case = TRUE)
    domain_label <- if (is_lipidomics) "Lipidomics" else "Metabolomics"
    feature_label <- if (is_lipidomics) "Lipids" else "Metabolites"
    feature_label_lc <- tolower(feature_label)

    message("Generating ", domain_label, " PowerPoint presentation...")

    metab_cfg  <- config$modes$metabolomics %||% config$modes$lipidomics
    de_cfg     <- metab_cfg$de %||% list()
    norm_cfg   <- metab_cfg$normalization %||% list()
    diag_dir   <- file.path(out_dir, "Diagnostic_plots")

    # Project info
    project_name <- config$project$name %||% paste(domain_label, "Analysis")
    analyst      <- gsub("_", " ", config$project$analyst %||% "")
    date_str     <- format(Sys.Date(), "%B %d, %Y")
    cond_col     <- de_cfg$condition_column %||% metab_cfg$effects$color %||% "sample_type"

    # Group info
    groups <- if (!is.null(pre$meta) && cond_col %in% colnames(pre$meta)) {
        table(pre$meta[[cond_col]])
    } else NULL

    # DE stats
    n_features <- if (!is.null(pre$expr_filt)) nrow(pre$expr_filt) else NA
    n_samples  <- if (!is.null(pre$expr_filt)) ncol(pre$expr_filt) else NA
    n_de       <- if (!is.null(de_res$summary_df) && "pass_any_contrast" %in% colnames(de_res$summary_df)) {
        sum(de_res$summary_df$pass_any_contrast == 1, na.rm = TRUE)
    } else 0
    de_method  <- de_res$method %||% "unknown"

    # Thresholds
    p_cut <- de_cfg$p_cutoff %||% 0.05
    lfc_cut <- if (!is.null(de_cfg$logfc_cutoff)) de_cfg$logfc_cutoff else log2(de_cfg$linear_fc_cutoff %||% 1.5)

    # Contrast info
    contrasts <- names(de_res$de_tables)

    # Name resolver
    name_lut <- NULL
    if (!is.null(pre$row_data) && "Name" %in% colnames(pre$row_data)) {
        name_lut <- stats::setNames(as.character(pre$row_data$Name),
                                     as.character(pre$row_data$feature_id))
    }
    resolve <- function(ids) {
        if (!is.null(name_lut)) ifelse(is.na(name_lut[ids]), ids, name_lut[ids])
        else ids
    }

    # Footer
    footer_txt <- paste0("multiomics-core | ", analyst, " | ", date_str,
                          "   |   Generated by AI-assisted pipeline")

    # ---- Create presentation ----
    pptx <- officer::read_pptx()

    # ================================================================
    # SLIDE 1: Title
    # ================================================================
    title_context <- sprintf(
        "Title slide for %s analysis. Project: %s. Analyst: %s. "
        , domain_label, project_name, analyst)
    if (!is.null(groups)) {
        title_context <- paste0(title_context,
            sprintf("Groups: %s. ", paste(sprintf("%s (n=%d)", names(groups), as.integer(groups)), collapse = ", ")))
    }
    title_context <- paste0(title_context,
        sprintf("%d %s, %d samples, %d significant DE.", n_features, feature_label_lc, n_samples, n_de))
    title_bl <- generate_slide_bottom_line(title_context, config = config, mode = "metabolomics")

    pptx <- officer::add_slide(pptx, layout = "Blank", master = "Office Theme")

    # Main title
    pptx <- officer::ph_with(pptx,
        value = officer::fpar(
            officer::ftext(paste(domain_label, "Analysis"), officer::fp_text(
                font.size = 28, bold = TRUE, color = COL_TITLE, font.family = "Calibri"))
        ),
        location = officer::ph_location(left = 1, top = 1.5, width = 8, height = 0.8))

    # Project name
    pptx <- officer::ph_with(pptx,
        value = officer::fpar(
            officer::ftext(project_name, officer::fp_text(
                font.size = 18, color = COL_ACCENT, font.family = "Calibri"))
        ),
        location = officer::ph_location(left = 1, top = 2.4, width = 8, height = 0.5))

    # Analyst / date line
    pptx <- officer::ph_with(pptx,
        value = officer::fpar(
            officer::ftext(paste0(analyst, "  |  ", date_str), officer::fp_text(
                font.size = 12, color = COL_SUBTLE, font.family = "Calibri"))
        ),
        location = officer::ph_location(left = 1, top = 3.1, width = 8, height = 0.4))

    # Key stats row
    groups_str <- if (!is.null(groups)) {
        paste(sprintf("%s (n=%d)", names(groups), as.integer(groups)), collapse = "  vs  ")
    } else "N/A"
    stats_text <- sprintf("%d %s  |  %d samples  |  %s  |  %d DE %s",
                           n_features, feature_label_lc, n_samples, groups_str, n_de, feature_label_lc)
    pptx <- officer::ph_with(pptx,
        value = officer::fpar(officer::ftext(stats_text, fp_body(11))),
        location = officer::ph_location(left = 1, top = 4.0, width = 8, height = 0.4))

    # AI bottom line
    if (nzchar(title_bl)) {
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(
                officer::ftext("AI Summary: ", fp_bl_prefix()),
                officer::ftext(title_bl, fp_bottom_line())
            ),
            location = officer::ph_location(left = 1, top = 5.0, width = 8, height = 1.0))
    }

    pptx <- officer::ph_with(pptx,
        value = officer::fpar(officer::ftext(footer_txt, fp_footer())),
        location = officer::ph_location(left = 0.3, top = 7.1, width = 9.4, height = 0.3))

    # ================================================================
    # SLIDE 2: Study Design
    # ================================================================
    design_context <- sprintf(
        "Study design. %d samples, %d %s. Groups: %s. DE method: %s. Normalization: %s + %s.",
        n_samples, n_features, feature_label_lc, groups_str, de_method,
        norm_cfg$transform %||% "glog10", norm_cfg$scaling %||% "auto")
    design_bl <- generate_slide_bottom_line(design_context, config = config, mode = "metabolomics")

    pptx <- officer::add_slide(pptx, layout = "Blank", master = "Office Theme")
    pptx <- officer::ph_with(pptx,
        value = officer::fpar(officer::ftext("Study Design", fp_slide_title())),
        location = officer::ph_location(left = 0.5, top = 0.3, width = 9, height = 0.5))

    # Design info as clean formatted text
    design_items <- list(
        list(label = "Samples", value = as.character(n_samples)),
        list(label = feature_label, value = as.character(n_features)),
        list(label = "Groups", value = groups_str),
        list(label = "DE Method", value = de_method),
        list(label = "Transform", value = norm_cfg$transform %||% "glog10"),
        list(label = "Scaling", value = norm_cfg$scaling %||% "auto"),
        list(label = "P-value cutoff", value = as.character(p_cut)),
        list(label = "log2(FC) cutoff", value = as.character(lfc_cut)),
        list(label = paste("DE", feature_label_lc), value = as.character(n_de))
    )

    # Build formatted block list
    design_fpars <- lapply(design_items, function(item) {
        officer::fpar(
            officer::ftext(paste0(item$label, ":  "), fp_body_bold(12)),
            officer::ftext(item$value, fp_body(12))
        )
    })
    pptx <- officer::ph_with(pptx,
        value = do.call(officer::block_list, design_fpars),
        location = officer::ph_location(left = 1.5, top = 1.2, width = 7, height = 4.2))

    if (nzchar(design_bl)) {
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(
                officer::ftext("AI Summary: ", fp_bl_prefix()),
                officer::ftext(design_bl, fp_bottom_line())
            ),
            location = officer::ph_location(left = 0.5, top = 5.7, width = 9, height = 1.1))
    }
    pptx <- officer::ph_with(pptx,
        value = officer::fpar(officer::ftext(footer_txt, fp_footer())),
        location = officer::ph_location(left = 0.3, top = 7.1, width = 9.4, height = 0.3))

    # ================================================================
    # SLIDE 3: PCA
    # ================================================================
    pca_png <- file.path(diag_dir, "PCA_PC1.vs.PC2.png")
    if (file.exists(pca_png)) {
        pca_context <- sprintf(
            "PCA plot (PC1 vs PC2). %d samples, groups: %s.", n_samples, groups_str)
        pca_bl <- generate_slide_bottom_line(pca_context, image_path = pca_png, config = config, mode = "metabolomics")
        pptx <- add_figure_slide(pptx, "PCA: PC1 vs PC2", pca_png,
                                  subtitle = sprintf("%d samples | %s", n_samples, groups_str),
                                  bl = pca_bl, footer = footer_txt)
    }

    # ================================================================
    # SLIDE 4+: Volcano (per contrast)
    # ================================================================
    for (ctr in contrasts) {
        volcano_png <- file.path(diag_dir, paste0("volcano_", ctr, ".png"))
        if (!file.exists(volcano_png)) next

        ctr_tbl <- de_res$de_tables[[ctr]]
        n_up <- sum(!is.na(ctr_tbl$P.Value) & ctr_tbl$P.Value < p_cut &
                     !is.na(ctr_tbl$logFC) & ctr_tbl$logFC >= lfc_cut, na.rm = TRUE)
        n_down <- sum(!is.na(ctr_tbl$P.Value) & ctr_tbl$P.Value < p_cut &
                       !is.na(ctr_tbl$logFC) & ctr_tbl$logFC <= -lfc_cut, na.rm = TRUE)

        # Top hits
        sig_idx <- which(!is.na(ctr_tbl$P.Value) & ctr_tbl$P.Value < p_cut &
                          abs(ctr_tbl$logFC) >= lfc_cut)
        top_names <- ""
        if (length(sig_idx) > 0) {
            ordered_idx <- sig_idx[order(abs(ctr_tbl$logFC[sig_idx]), decreasing = TRUE)]
            top_ids <- ctr_tbl$feature_id[utils::head(ordered_idx, 8)]
            top_names <- paste(resolve(top_ids), collapse = ", ")
        }

        vol_stats <- sprintf(
            '{"contrast":"%s","up":%d,"down":%d,"total_tested":%d,"top":"%s"}',
            gsub("_vs_", " vs ", ctr), n_up, n_down, nrow(ctr_tbl), top_names)
        vol_context <- sprintf(
            "Volcano plot: %s. %d up, %d down (p<%.2g, |log2FC|>=%.1f). Top: %s.",
            gsub("_vs_", " vs ", ctr), n_up, n_down, p_cut, lfc_cut, top_names)
        vol_bl <- generate_slide_bottom_line(vol_context, stats_json = vol_stats,
                                              image_path = volcano_png, config = config, mode = "metabolomics")

        subtitle <- sprintf("Up: %d  |  Down: %d  |  Total DE: %d  (p < %s, |log2FC| >= %s)",
                             n_up, n_down, n_up + n_down, p_cut, lfc_cut)
        pptx <- add_figure_slide(pptx,
            paste("Volcano Plot:", gsub("_vs_", " vs ", ctr)),
            volcano_png, subtitle = subtitle, bl = vol_bl, footer = footer_txt)
    }

    # ================================================================
    # SLIDE: MA plot (per contrast)
    # ================================================================
    for (ctr in contrasts) {
        ma_png <- file.path(diag_dir, paste0("MA_", ctr, ".png"))
        if (!file.exists(ma_png)) next

        ma_context <- sprintf("MA plot: %s. log2(FC) vs average expression.",
                               gsub("_vs_", " vs ", ctr))
        ma_bl <- generate_slide_bottom_line(ma_context, image_path = ma_png, config = config, mode = "metabolomics")
        pptx <- add_figure_slide(pptx,
            paste("MA Plot:", gsub("_vs_", " vs ", ctr)),
            ma_png, bl = ma_bl, footer = footer_txt)
    }

    # ================================================================
    # SLIDE: Top DE metabolites table (per contrast)
    # ================================================================
    for (ctr in contrasts) {
        ctr_tbl <- de_res$de_tables[[ctr]]
        if (is.null(ctr_tbl) || nrow(ctr_tbl) == 0) next

        sig <- ctr_tbl[!is.na(ctr_tbl$P.Value) & ctr_tbl$P.Value < p_cut &
                        abs(ctr_tbl$logFC) >= lfc_cut, , drop = FALSE]
        sig <- sig[order(sig$P.Value), ]
        top_n <- min(15, nrow(sig))
        if (top_n == 0) next

        top_tbl <- sig[seq_len(top_n), , drop = FALSE]
        display_df <- data.frame(
            Metabolite  = resolve(as.character(top_tbl$feature_id)),
            `log2(FC)`  = round(top_tbl$logFC, 2),
            `P.Value`   = signif(top_tbl$P.Value, 3),
            `adj.P.Val` = signif(top_tbl$adj.P.Val, 3),
            check.names = FALSE, stringsAsFactors = FALSE
        )

        tbl_context <- sprintf(
            "Top %d DE %s for %s. Top: %s (log2FC=%.2f, p=%.2e). Total DE: %d.",
            top_n, feature_label_lc, gsub("_vs_", " vs ", ctr),
            display_df$Metabolite[1], display_df$`log2(FC)`[1],
            display_df$P.Value[1], nrow(sig))
        tbl_bl <- generate_slide_bottom_line(tbl_context, config = config, mode = "metabolomics")

        pptx <- officer::add_slide(pptx, layout = "Blank", master = "Office Theme")
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(officer::ftext(
                paste("Top DE", paste0(feature_label, ":"), gsub("_vs_", " vs ", ctr)),
                fp_slide_title())),
            location = officer::ph_location(left = 0.5, top = 0.3, width = 9, height = 0.5))
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(officer::ftext(
                sprintf("%d significant %s (p < %s, |log2FC| >= %s)", nrow(sig), feature_label_lc, p_cut, lfc_cut),
                fp_label())),
            location = officer::ph_location(left = 0.5, top = 0.85, width = 9, height = 0.3))

        # Table on left, AI sidebar on right to avoid overlap
        tbl_w <- if (nzchar(tbl_bl)) 6.0 else 9.0
        if (requireNamespace("flextable", quietly = TRUE)) {
            ft <- flextable::flextable(display_df)
            ft <- flextable::fontsize(ft, size = 8, part = "body")
            ft <- flextable::fontsize(ft, size = 9, part = "header")
            ft <- flextable::bold(ft, part = "header")
            ft <- flextable::color(ft, color = "#FFFFFF", part = "header")
            ft <- flextable::bg(ft, bg = COL_TITLE, part = "header")
            ft <- flextable::align(ft, j = 2:4, align = "right", part = "all")
            ft <- flextable::autofit(ft)
            ft <- flextable::padding(ft, padding = 3, part = "all")
            pptx <- officer::ph_with(pptx, value = ft,
                location = officer::ph_location(left = 0.5, top = 1.2, width = tbl_w, height = 5.5))
        } else {
            tbl_text <- paste(
                apply(display_df, 1, function(r) paste(r, collapse = "    ")),
                collapse = "\n")
            pptx <- officer::ph_with(pptx,
                value = officer::fpar(officer::ftext(tbl_text, fp_body(9))),
                location = officer::ph_location(left = 0.5, top = 1.2, width = tbl_w, height = 5.5))
        }

        if (nzchar(tbl_bl)) {
            pptx <- officer::ph_with(pptx,
                value = officer::fpar(
                    officer::ftext("AI Summary: ", fp_bl_prefix()),
                    officer::ftext(tbl_bl, fp_bottom_line())),
                location = officer::ph_location(left = 6.8, top = 1.2, width = 2.9, height = 5.5))
        }
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(officer::ftext(footer_txt, fp_footer())),
            location = officer::ph_location(left = 0.3, top = 7.1, width = 9.4, height = 0.3))
    }

    # ================================================================
    # SLIDE: DE Heatmap
    # ================================================================
    heatmap_pngs <- list.files(diag_dir, pattern = "^heatmap_top(10|25)_DE\\.png$", full.names = TRUE)
    if (length(heatmap_pngs) > 0) {
        hm_png <- heatmap_pngs[grepl("top25", heatmap_pngs)]
        if (length(hm_png) == 0) hm_png <- heatmap_pngs[1]
        hm_png <- hm_png[1]

        hm_context <- sprintf("Heatmap of top DE %s with hierarchical clustering.", feature_label_lc)
        hm_bl <- generate_slide_bottom_line(hm_context, image_path = hm_png, config = config, mode = "metabolomics")
        pptx <- add_figure_slide(pptx, paste("DE", feature_label, "Heatmap"), hm_png,
                                  subtitle = "Z-scored glog10 intensities | hierarchical clustering",
                                  bl = hm_bl, footer = footer_txt)
    }

    # ================================================================
    # SLIDE: P-value histogram (per contrast)
    # ================================================================
    for (ctr in contrasts) {
        phist_png <- file.path(diag_dir, paste0("pvalue_hist_", ctr, ".png"))
        if (!file.exists(phist_png)) next

        phist_context <- sprintf("P-value distribution for %s.", gsub("_vs_", " vs ", ctr))
        phist_bl <- generate_slide_bottom_line(phist_context, image_path = phist_png, config = config, mode = "metabolomics")
        pptx <- add_figure_slide(pptx,
            paste("P-value Distribution:", gsub("_vs_", " vs ", ctr)),
            phist_png, bl = phist_bl, footer = footer_txt)
    }

    # ================================================================
    # SLIDE: Random Forest importance
    # ================================================================
    rf_png <- file.path(diag_dir, "rf_importance.png")
    if (file.exists(rf_png) && !is.null(feature_sel_res$rf)) {
        rf_imp <- feature_sel_res$rf$importance
        rf_top_names <- character(0)
        if (!is.null(rf_imp) && is.data.frame(rf_imp) && nrow(rf_imp) > 0) {
            imp_col <- if ("MeanDecreaseAccuracy" %in% colnames(rf_imp)) {
                "MeanDecreaseAccuracy"
            } else if ("importance" %in% colnames(rf_imp)) {
                "importance"
            } else setdiff(colnames(rf_imp), "feature_id")[1]
            sorted_idx <- order(rf_imp[[imp_col]], decreasing = TRUE)
            top_rf <- rf_imp[utils::head(sorted_idx, 5), , drop = FALSE]
            top_ids <- if ("feature_id" %in% colnames(top_rf)) as.character(top_rf$feature_id) else rownames(top_rf)
            rf_top_names <- resolve(top_ids)
        }

        rf_context <- sprintf("Random Forest importance. Top features: %s.",
                               paste(rf_top_names, collapse = ", "))
        rf_bl <- generate_slide_bottom_line(rf_context, image_path = rf_png, config = config, mode = "metabolomics")
        pptx <- add_figure_slide(pptx, "Random Forest: Variable Importance", rf_png,
                                  subtitle = paste("Top:", paste(rf_top_names, collapse = ", ")),
                                  bl = rf_bl, footer = footer_txt)
    }

    # ================================================================
    # SLIDE: PLS-DA scores
    # ================================================================
    plsda_png <- file.path(diag_dir, "plsda_scores.png")
    if (file.exists(plsda_png) && !is.null(feature_sel_res$plsda)) {
        plsda_context <- sprintf("PLS-DA scores. Supervised separation of %s.", groups_str)
        plsda_bl <- generate_slide_bottom_line(plsda_context, image_path = plsda_png, config = config, mode = "metabolomics")
        pptx <- add_figure_slide(pptx, "PLS-DA: Scores Plot", plsda_png,
                                  subtitle = sprintf("Supervised classification | %s", groups_str),
                                  bl = plsda_bl, footer = footer_txt)
    }

    # ================================================================
    # SLIDE: PLS-DA VIP
    # ================================================================
    vip_png <- file.path(diag_dir, "plsda_vip_scores.png")
    if (file.exists(vip_png) && !is.null(feature_sel_res$plsda)) {
        vip_context <- "PLS-DA VIP scores. Features with VIP > 1 are important discriminators."
        vip_bl <- generate_slide_bottom_line(vip_context, image_path = vip_png, config = config, mode = "metabolomics")
        pptx <- add_figure_slide(pptx, "PLS-DA: VIP Scores", vip_png,
                                  subtitle = "Variable Importance in Projection | VIP > 1 threshold",
                                  bl = vip_bl, footer = footer_txt)
    }

    # ================================================================
    # SLIDE: Sample correlation
    # ================================================================
    cor_png <- file.path(diag_dir, "sample_correlation_heatmap.png")
    if (file.exists(cor_png)) {
        cor_context <- "Sample-to-sample Pearson correlation after normalization."
        cor_bl <- generate_slide_bottom_line(cor_context, image_path = cor_png, config = config, mode = "metabolomics")
        pptx <- add_figure_slide(pptx, "Sample Correlation Heatmap", cor_png,
                                  bl = cor_bl, footer = footer_txt)
    }

    # ================================================================
    # SLIDE: Enrichment (only if results exist)
    # ================================================================
    enrich_dir <- file.path(out_dir, "Enrichment")
    if (dir.exists(enrich_dir)) {
        enrich_pngs <- list.files(enrich_dir, pattern = "\\.png$",
                                   full.names = TRUE, recursive = TRUE)
        for (ep in enrich_pngs) {
            ep_label <- gsub("_", " ", tools::file_path_sans_ext(basename(ep)))
            ep_label <- tools::toTitleCase(ep_label)
            ep_context <- sprintf("Pathway enrichment: %s.", ep_label)
            ep_bl <- generate_slide_bottom_line(ep_context, image_path = ep, config = config, mode = "metabolomics")
            pptx <- add_figure_slide(pptx, paste("Enrichment:", ep_label), ep,
                                      bl = ep_bl, footer = footer_txt)
        }
    }

    # ---- Save to Results root (parent of mode dir) ----
    results_root <- dirname(out_dir)
    out_file <- file.path(results_root, paste0(tolower(domain_label), "_summary.pptx"))
    print(pptx, target = out_file)
    message("PowerPoint saved: ", out_file)
    out_file
}
