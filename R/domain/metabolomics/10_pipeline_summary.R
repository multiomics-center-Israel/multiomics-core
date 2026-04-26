#' Metabolomics Pipeline Summary HTML Generation
#'
#' Generates a dark-themed HTML summary of the metabolomics pipeline workflow.
#' Data is extracted from pipeline objects and config. Descriptions polished
#' by Claude Code CLI with R fallback.
#'
#' Follows the pattern of R/domain/rnaseq/10_pipeline_summary.R.

# ==============================================================================
# DATA COLLECTION
# ==============================================================================

#' Collect all metabolomics pipeline statistics for the summary
#'
#' @param config          Full pipeline config
#' @param pre             Preprocessed data list (expr_filt, meta, etc.)
#' @param de_res          DE results list
#' @param feature_sel_res Feature selection results (rf, plsda)
#' @param enrichment_res  Enrichment results list
#' @return Named list with all stats organized by section
collect_metab_pipeline_stats <- function(config, pre, de_res,
                                         feature_sel_res, enrichment_res) {

    metab_cfg <- config$modes$metabolomics
    norm_cfg  <- metab_cfg$normalization %||% list()

    # --- Project info ---
    analyst <- gsub("_", " ", config$project$analyst %||% "")

    project <- list(
        name     = config$project$name %||% "Metabolomics Analysis",
        analyst  = analyst,
        date     = format(Sys.Date(), "%B %d, %Y"),
        date_iso = format(Sys.Date(), "%Y-%m-%d")
    )

    # --- Feature & sample overview ---
    n_samples   <- if (!is.null(pre$expr_filt)) ncol(pre$expr_filt) else NA
    n_features_raw <- if (!is.null(pre$expr_raw)) nrow(pre$expr_raw) else
                      if (!is.null(pre$expr_filt)) nrow(pre$expr_filt) else NA
    n_features  <- if (!is.null(pre$expr_filt)) nrow(pre$expr_filt) else NA

    group_col <- metab_cfg$de$condition_column %||%
                 metab_cfg$effects$color %||% "sample_type"
    n_groups <- NA
    groups <- character()
    if (!is.null(pre$meta) && !is.null(group_col) && group_col %in% names(pre$meta)) {
        groups <- unique(as.character(pre$meta[[group_col]]))
        n_groups <- length(groups)
    }

    # --- DE statistics ---
    # Single source of truth: use pipeline pass columns (pass.{cn}) from summary_df.
    # Fallback: recompute from p-value + FC thresholds.
    n_de_total <- 0; n_de_up <- 0; n_de_down <- 0
    de_contrasts <- list()

    p_cut  <- metab_cfg$de$p_cutoff %||% 0.05
    fc_lin <- metab_cfg$de$linear_fc_cutoff %||% 1.0
    log2_fc <- if (fc_lin > 1) log2(fc_lin) else 0

    if (!is.null(de_res$summary_df)) {
        sdf <- de_res$summary_df
        pass_cols <- grep("^pass\\.", names(sdf), value = TRUE)
        pass_cols <- setdiff(pass_cols, "pass_any_contrast")

        if (length(pass_cols) > 0) {
            for (pcol in pass_cols) {
                cn <- sub("^pass\\.", "", pcol)
                is_sig <- !is.na(sdf[[pcol]]) & sdf[[pcol]] == 1
                fc_col <- paste0("linearFC.", cn)
                if (fc_col %in% names(sdf)) {
                    fc_vals <- as.numeric(sdf[[fc_col]])
                    up <- sum(is_sig & fc_vals > 0, na.rm = TRUE)
                    dn <- sum(is_sig & fc_vals < 0, na.rm = TRUE)
                } else {
                    up <- 0; dn <- 0
                }
                de_contrasts[[cn]] <- list(
                    name = cn, total = sum(is_sig),
                    up = up, down = dn, tested = nrow(sdf)
                )
                n_de_total <- n_de_total + sum(is_sig)
                n_de_up <- n_de_up + up; n_de_down <- n_de_down + dn
            }
        }
    }

    # Fallback: use per-contrast DE tables
    if (length(de_contrasts) == 0) {
        de_tables <- de_res$de_tables %||% de_res$tables
        if (!is.null(de_tables) && length(de_tables) > 0) {
            for (cn in names(de_tables)) {
                tbl <- de_tables[[cn]]
                if (!is.data.frame(tbl)) next
                padj_col <- intersect(c("P.Value", "pvalue", "padj", "adj.P.Val"), names(tbl))[1]
                lfc_col  <- intersect(c("log2FoldChange", "logFC", "log2FC"), names(tbl))[1]
                if (is.na(padj_col)) next
                sig <- !is.na(tbl[[padj_col]]) & tbl[[padj_col]] <= p_cut
                if (log2_fc > 0 && !is.na(lfc_col)) sig <- sig & (abs(tbl[[lfc_col]]) >= log2_fc)
                up <- if (!is.na(lfc_col)) sum(sig & tbl[[lfc_col]] > 0, na.rm = TRUE) else 0
                dn <- if (!is.na(lfc_col)) sum(sig & tbl[[lfc_col]] < 0, na.rm = TRUE) else 0
                de_contrasts[[cn]] <- list(
                    name = cn, total = sum(sig, na.rm = TRUE),
                    up = up, down = dn, tested = sum(!is.na(tbl[[padj_col]]))
                )
                n_de_total <- n_de_total + sum(sig, na.rm = TRUE)
                n_de_up <- n_de_up + up; n_de_down <- n_de_down + dn
            }
        }
    }

    # --- Feature selection statistics ---
    rf_top_n <- 0
    plsda_top_n <- 0

    if (!is.null(feature_sel_res$rf)) {
        rf_res <- feature_sel_res$rf
        rf_imp <- rf_res$importance_df %||% rf_res$importance
        if (!is.null(rf_imp) && is.data.frame(rf_imp)) {
            rf_top_n <- metab_cfg$rf$top_n %||% min(20, nrow(rf_imp))
        }
    }
    if (!is.null(feature_sel_res$plsda)) {
        plsda_res <- feature_sel_res$plsda
        plsda_vip <- plsda_res$vip_df %||% plsda_res$vip
        if (!is.null(plsda_vip) && is.data.frame(plsda_vip)) {
            plsda_top_n <- metab_cfg$plsda$vip_top_n %||% min(15, nrow(plsda_vip))
        }
    }

    # --- Enrichment statistics ---
    n_enriched_total <- 0
    if (!is.null(enrichment_res) && is.list(enrichment_res)) {
        count_sig <- function(x) {
            if (is.data.frame(x) && "padj" %in% names(x)) {
                return(sum(!is.na(x$padj) & x$padj < 0.05, na.rm = TRUE))
            }
            if (is.list(x)) return(sum(vapply(x, count_sig, integer(1))))
            0L
        }
        n_enriched_total <- count_sig(enrichment_res)
    }

    # --- Subtitle ---
    subtitle <- project$name
    if (n_groups >= 2) {
        subtitle <- sprintf("%s &mdash; %s",
                            project$name,
                            paste(groups, collapse = " vs "))
    }

    list(
        project  = project,
        subtitle = subtitle,
        overview = list(
            n_samples = n_samples %||% 0, n_groups = n_groups %||% 0,
            groups = groups, n_features = n_features %||% 0,
            n_features_fmt = format(n_features %||% 0, big.mark = ","),
            n_de_features_fmt = format(n_de_total, big.mark = ",")
        ),
        de = list(total = n_de_total, up = n_de_up, down = n_de_down,
                  contrasts = de_contrasts),
        feature_selection = list(
            rf_enabled = isTRUE(metab_cfg$rf$run_rf),
            rf_top_n = rf_top_n,
            plsda_enabled = isTRUE(metab_cfg$plsda$run_plsda),
            plsda_top_n = plsda_top_n
        ),
        enrichment = list(total = n_enriched_total),
        methods = list(
            input_format = metab_cfg$input$format %||% "cd_raw",
            sample_norm  = metab_cfg$preprocessing$chosen_norm %||% norm_cfg$sample_norm %||% "pqn",
            transform    = norm_cfg$transform %||% "log2",
            scaling      = norm_cfg$scaling %||% "none",
            de_method    = metab_cfg$de$method %||% "limma",
            p_cutoff     = metab_cfg$de$p_cutoff %||% 0.05,
            fc_cutoff    = metab_cfg$de$linear_fc_cutoff %||% 1.0
        ),
        filtering = list(
            features_before = n_features_raw %||% 0,
            features_after = n_features %||% 0,
            before_fmt = format(n_features_raw %||% 0, big.mark = ","),
            after_fmt  = format(n_features %||% 0, big.mark = ","),
            pct_retained = if (!is.na(n_features_raw) && n_features_raw > 0 && !is.na(n_features))
                round(100 * n_features / n_features_raw, 1) else NA
        ),
        feature_flags = list(
            rf_enabled          = isTRUE(metab_cfg$rf$run_rf),
            plsda_enabled       = isTRUE(metab_cfg$plsda$run_plsda),
            enrichment_enabled  = isTRUE(metab_cfg$enrichment$run_enrichment)
        )
    )
}


# ==============================================================================
# CLAUDE CODE CLI — generate the complete HTML via `claude --print`
# ==============================================================================

#' Generate metabolomics pipeline summary HTML using Claude Code CLI
#'
#' @param stats       Pipeline stats from collect_metab_pipeline_stats()
#' @param comment_cfg Commentary config section (for model selection)
#' @return Complete HTML document string, or NULL on failure
generate_metab_summary_with_claude <- function(stats, comment_cfg) {

    if (Sys.which("claude") == "") {
        message("  'claude' CLI not found in PATH.")
        return(NULL)
    }

    model <- comment_cfg$claude_code_model %||% "sonnet"
    stats_json <- jsonlite::toJSON(stats, auto_unbox = TRUE, pretty = TRUE)

    tmp_stats <- tempfile(fileext = ".json")
    on.exit(unlink(tmp_stats), add = TRUE)
    writeLines(stats_json, tmp_stats)

    prompt <- paste0(
        "Read the pipeline statistics JSON file at '", tmp_stats, "' and generate a COMPLETE, self-contained HTML document ",
        "for a metabolomics pipeline workflow summary page.\n\n",
        "DESIGN REQUIREMENTS:\n",
        "- Dark theme: background #0a0e17, cards #111827, border #1e2d45\n",
        "- Google Fonts: JetBrains Mono (monospace elements) + DM Sans (body text)\n",
        "- Gradient h1 text (white to emerald-green #34d399), subtle background grid, two blurred glow circles\n",
        "- Vertical timeline with gradient connector line on the left\n",
        "- Each step: node circle (emoji + step number) + card (title, phase badge, description, tags)\n",
        "- Cards have hover effects (translateX + border glow + shadow)\n",
        "- Steps animate in with fadeUp (staggered delays)\n",
        "- Responsive: mobile breakpoint at 640px\n\n",
        "STRUCTURE:\n",
        "Header: badge 'multiomics-core pipeline', h1 'Metabolomics Analysis Workflow', subtitle from stats, ",
        "analyst + date, meta row (Samples, Groups, Features, DE Features) with actual numbers from stats.\n\n",
        "Pipeline sections:\n",
        "1. 'Data Processing' section:\n",
        "   - Step 01 Data Import (phase-input, cyan, emoji inbox) - input format from stats.methods.input_format\n",
        "   - Step 02 Feature Filtering (phase-filter, orange, emoji magnifier) - before/after feature counts\n",
        "   - Step 03 Normalization (phase-norm, rose, emoji scales) - sample_norm, transform, scaling methods\n",
        "2. 'Statistical Analysis' section:\n",
        "   - Step 04 Quality Control (phase-qc, green, emoji checkmark) - PCA, heatmaps, sample QC\n",
        "   - Step 05 Differential Expression (phase-de, violet, emoji chart) - DE method, sig feature counts\n",
        "   - Step 06 Feature Selection (phase-fs, amber, emoji star) - RF importance, PLS-DA VIP (if enabled)\n",
        "   - Step 07 Pathway Enrichment (phase-enrich, blue, emoji compass) - enrichment results (if enabled)\n",
        "   - Step 08 Report & Outputs (phase-out, green, emoji folder) - HTML report, Excel, CSV\n\n",
        "Footer: analyst + date.\n\n",
        "Write SCIENTIFICALLY PRECISE descriptions. Use HTML entities for special chars.\n",
        "Output ONLY the complete HTML document - no markdown fences, no explanations."
    )

    cmd <- sprintf(
        "unset CLAUDECODE; claude --print --model %s --no-session-persistence %s",
        model, shQuote(prompt)
    )

    result <- tryCatch({
        raw <- system(cmd, intern = TRUE, timeout = 180)
        html <- paste(raw, collapse = "\n")

        html <- gsub("^```html\\s*\n?", "", html)
        html <- gsub("\n?```\\s*$", "", html)

        if (grepl("<!DOCTYPE html>", html, ignore.case = TRUE) &&
            grepl("phase-de", html) &&
            grepl("pipeline", html)) {
            message("  Claude Code generated metabolomics pipeline summary HTML")
            html
        } else {
            message("  Claude output did not match expected HTML structure. Using R fallback.")
            NULL
        }
    }, error = function(e) {
        message("  Claude Code call failed: ", e$message, ". Using R fallback.")
        NULL
    })

    result
}


# ==============================================================================
# R FALLBACK — template-based HTML generation
# ==============================================================================

#' Build HTML for a single pipeline step card
#' @noRd
build_metab_step_html <- function(num, emoji, title, phase_class, phase_label,
                                   description, tags, stats_html = "") {
    tags_html <- paste(
        sprintf('            <span class="tag">%s</span>', tags),
        collapse = "\n"
    )
    sprintf(
'      <div class="step %s">
        <div class="node">
          <div class="node-circle">%s</div>
          <div class="node-num">%02d</div>
        </div>
        <div class="card">
          <div class="card-title">
            %s
            <span class="card-phase">%s</span>
          </div>
          <div class="card-desc">%s</div>
%s          <div class="card-tags">
%s
          </div>
        </div>
      </div>',
        phase_class, emoji, num,
        title, phase_label, description,
        if (nzchar(stats_html)) paste0(stats_html, "\n") else "",
        tags_html
    )
}

#' Build stats row HTML
#' @noRd
build_metab_stats_row <- function(...) {
    items <- list(...)
    html_items <- vapply(items, function(item) {
        if (isTRUE(item$arrow)) {
            '<div class="stat" style="color:var(--text-dim)">&rarr;</div>'
        } else {
            extra <- if (!is.null(item$margin)) sprintf(' style="margin-left:%s"', item$margin) else ""
            sprintf('<div class="stat"%s><span class="stat-val" style="color:%s">%s</span><span class="stat-label">%s</span></div>',
                    extra, item$color %||% "var(--text)", item$value, item$label)
        }
    }, character(1))
    sprintf('          <div class="stats-row">\n%s\n          </div>',
            paste(paste0("            ", html_items), collapse = "\n"))
}

#' Build a single flowchart node HTML
#' @noRd
metab_fc_node <- function(icon, title, detail, variant) {
    sprintf(
'<div class="fc-node %s">
  <div class="fc-node-icon">%s</div>
  <div class="fc-node-title">%s</div>
  <div class="fc-node-detail">%s</div>
</div>', variant, icon, title, detail)
}

#' Build vertical connector arrow
#' @noRd
metab_fc_arrow <- function() {
    '<div class="fc-connector"></div>'
}

#' Build T-split (stem down + horizontal bar + drop lines with arrows)
#' @noRd
metab_fc_t_split <- function(n_cols) {
    cls <- sprintf("fc-cols-%d", n_cols)
    drops <- paste(rep(
        '<div class="fc-drop-wrap"><div class="fc-drop-line"></div></div>',
        n_cols), collapse = "\n    ")
    sprintf(
'<div class="fc-t-split">
  <div class="fc-t-stem"></div>
  <div class="fc-branch-drops %s">
    <div class="fc-t-bar"></div>
    %s
  </div>
</div>', cls, drops)
}

#' Build merge connector (rise lines + horizontal bar + stem down with arrow)
#' @noRd
metab_fc_merge <- function(n_cols) {
    cls <- sprintf("fc-cols-%d", n_cols)
    rises <- paste(rep(
        '<div class="fc-rise-wrap"><div class="fc-rise-line"></div></div>',
        n_cols), collapse = "\n    ")
    sprintf(
'<div class="fc-merge-rises %s">
  <div class="fc-t-bar"></div>
  %s
</div>
<div class="fc-merge-stem"></div>', cls, rises)
}

#' Generate the metabolomics analysis overview flowchart HTML
#'
#' Builds a visual data-flow diagram from pipeline stats and feature flags.
#' Nodes are conditionally included based on which pipeline steps were enabled.
#'
#' @param stats Pipeline stats from collect_metab_pipeline_stats()
#' @return HTML string for the flowchart section
#' @noRd
generate_metab_flowchart_html <- function(stats) {

    ff  <- stats$feature_flags
    esc <- function(x) if (requireNamespace("htmltools", quietly = TRUE)) htmltools::htmlEscape(x) else x

    parts <- list()

    # Section header
    parts <- c(parts, list(
        '<div class="fc-section">',
        '  <div class="fc-section-title">Analysis Overview</div>',
        '  <div class="fc-section-sub">Data flow through the pipeline</div>',
        '  <div class="fc-flow">'
    ))

    # --- Input Data ---
    input_detail <- sprintf("%s samples &middot; %s features",
                            stats$overview$n_samples,
                            format(stats$filtering$features_before, big.mark = ","))
    parts <- c(parts, list(
        metab_fc_node("\U0001F4E5", "Input Data", input_detail, "fc-node-input"),
        metab_fc_arrow()
    ))

    # --- Feature Filtering ---
    filt_detail <- sprintf("%s &rarr; %s features (%.1f%% retained)",
                           stats$filtering$before_fmt, stats$filtering$after_fmt,
                           stats$filtering$pct_retained %||% 0)
    parts <- c(parts, list(
        metab_fc_node("\U0001F50D", "Feature Filtering", filt_detail, "fc-node-filter"),
        metab_fc_arrow()
    ))

    # --- Normalization ---
    norm_parts <- c()
    if (stats$methods$sample_norm != "none") norm_parts <- c(norm_parts, toupper(stats$methods$sample_norm))
    if (stats$methods$transform != "none")   norm_parts <- c(norm_parts, stats$methods$transform)
    if (stats$methods$scaling != "none")     norm_parts <- c(norm_parts, stats$methods$scaling)
    norm_detail <- if (length(norm_parts) > 0) paste(norm_parts, collapse = " &middot; ") else "none"
    parts <- c(parts, list(
        metab_fc_node("\u2696\uFE0F", "Normalization", esc(norm_detail), "fc-node-norm"),
        metab_fc_arrow()
    ))

    # --- QC Diagnostics (always present) ---
    parts <- c(parts, list(
        metab_fc_node("\u2705", "QC Diagnostics", "PCA &middot; Heatmaps &middot; Sample QC", "fc-node-qc")
    ))

    # --- Downstream branches: DE (always), Feature Selection (conditional), Enrichment (conditional) ---
    downstream_nodes <- list()

    # DE is always present
    de_method <- toupper(stats$methods$de_method)
    de_detail <- sprintf("%s &middot; %s DE features", de_method,
                         format(stats$de$total, big.mark = ","))
    downstream_nodes <- c(downstream_nodes, list(
        metab_fc_node("\U0001F4CA", "Differential Expression", de_detail, "fc-node-de")
    ))

    # Feature Selection (conditional)
    if (isTRUE(ff$rf_enabled) || isTRUE(ff$plsda_enabled)) {
        fs_parts <- c()
        if (isTRUE(ff$rf_enabled))    fs_parts <- c(fs_parts, "RF")
        if (isTRUE(ff$plsda_enabled)) fs_parts <- c(fs_parts, "PLS-DA")
        fs_detail <- paste(fs_parts, collapse = " &middot; ")
        downstream_nodes <- c(downstream_nodes, list(
            metab_fc_node("\u2B50", "Feature Selection", fs_detail, "fc-node-fs")
        ))
    }

    # Enrichment (conditional)
    if (isTRUE(ff$enrichment_enabled)) {
        enrich_detail <- sprintf("%s pathways", stats$enrichment$total)
        downstream_nodes <- c(downstream_nodes, list(
            metab_fc_node("\U0001F9ED", "Enrichment", enrich_detail, "fc-node-enrich")
        ))
    }

    n_downstream <- length(downstream_nodes)

    if (n_downstream == 1) {
        # Single node (just DE) -- no split needed
        parts <- c(parts, list(
            metab_fc_arrow(),
            downstream_nodes[[1]]
        ))
    } else {
        # Multi-branch T-split / merge
        cols_cls <- sprintf("fc-cols-%d", n_downstream)
        row_html <- sprintf('<div class="fc-branch-row %s">%s</div>',
                            cols_cls, paste(downstream_nodes, collapse = ""))
        parts <- c(parts, list(
            metab_fc_t_split(n_downstream),
            row_html,
            metab_fc_merge(n_downstream)
        ))
    }

    # --- Report & Export ---
    parts <- c(parts, list(
        metab_fc_arrow(),
        metab_fc_node("\U0001F4C1", "Report &amp; Export", "HTML &middot; Excel &middot; CSV", "fc-node-report")
    ))

    # Close section
    parts <- c(parts, list(
        '  </div>',
        '</div>'
    ))

    paste(parts, collapse = "\n")
}


#' Generate metabolomics pipeline body HTML from stats (R fallback)
#' @noRd
generate_metab_summary_body_r <- function(stats) {

    esc <- function(x) if (requireNamespace("htmltools", quietly = TRUE)) htmltools::htmlEscape(x) else x

    # -- Data Processing --
    s_proc <- '      <div class="section-label">Data Processing <span class="team-badge team-multiomics">Multi-omics Core</span></div>\n'

    # Step 1: Data Import
    input_desc <- sprintf("Metabolomics data loaded from %s format input file.",
                           toupper(stats$methods$input_format))
    step1 <- build_metab_step_html(1, "\U0001F4E5", "Data Import", "phase-input", "input",
        input_desc, c(stats$methods$input_format))

    # Step 2: Feature Filtering — render dynamically based on whether
    # filtering actually changed the feature count.
    n_before <- stats$filtering$features_before %||% 0
    n_after  <- stats$filtering$features_after  %||% 0
    filtering_ran <- n_before > 0 && n_after > 0 && n_before != n_after
    filt_stats <- ""
    if (filtering_ran) {
        filt_desc <- "Low-abundance and unreliable features removed. Features retained based on detection frequency across samples."
        filt_stats <- build_metab_stats_row(
            list(value = stats$filtering$before_fmt, label = "before", color = "var(--accent-orange)"),
            list(arrow = TRUE),
            list(value = stats$filtering$after_fmt, label = "after", color = "var(--accent-green)"),
            list(value = "", label = sprintf("%.1f%% retained", stats$filtering$pct_retained %||% 0),
                 margin = "8px", color = "var(--text)"))
        filt_tags <- c("detection filter")
    } else if (n_before > 0) {
        filt_desc <- sprintf("No filtering applied; all %s features retained for downstream analysis.",
                             format(n_before, big.mark = ","))
        filt_tags <- c("no filtering")
    } else {
        filt_desc <- "Filtering step was skipped."
        filt_tags <- c("no filtering")
    }
    step2 <- build_metab_step_html(2, "\U0001F50D", "Feature Filtering", "phase-filter", "processing",
        filt_desc, filt_tags, filt_stats)

    # Step 3: Normalization
    # Human-readable normalization method names
    norm_display_names <- c(
        pqn = "PQN (Probabilistic Quotient Normalization)",
        tss = "TSS (Total Sum Scaling)",
        median = "Median normalization",
        eigenms = "EigenMS (SVD-based bias removal, ProteoMM)",
        eigenms_forced = "EigenMS Forced (SVD-based, fixed eigentrend removal)",
        sum = "Sum normalization",
        none = "None"
    )
    norm_method <- stats$methods$sample_norm
    norm_label <- norm_display_names[norm_method] %||% toupper(norm_method)

    norm_parts <- c()
    if (norm_method != "none")
        norm_parts <- c(norm_parts, sprintf("Sample normalization: %s", norm_label))
    if (stats$methods$transform != "none")
        norm_parts <- c(norm_parts, sprintf("%s transformation", stats$methods$transform))
    if (stats$methods$scaling != "none")
        norm_parts <- c(norm_parts, sprintf("%s scaling", stats$methods$scaling))
    if (length(norm_parts) == 0) norm_parts <- "No normalization applied"
    norm_desc <- paste(norm_parts, collapse = ". ") |> paste0(".")
    norm_tags <- c(toupper(norm_method), stats$methods$transform, stats$methods$scaling)
    norm_tags <- norm_tags[norm_tags != "none" & norm_tags != "NONE"]
    if (length(norm_tags) == 0) norm_tags <- "raw"
    step3 <- build_metab_step_html(3, "\u2696", "Normalization", "phase-norm", "processing",
        norm_desc, norm_tags)

    # -- Statistical Analysis --
    s_stat <- '      <div class="section-label">Statistical Analysis <span class="team-badge team-multiomics">Multi-omics Core</span></div>\n'

    # Step 4: QC
    step4 <- build_metab_step_html(4, "\u2713", "Quality Control", "phase-qc", "analysis",
        "Sample-level QC including PCA, hierarchical clustering, correlation analysis, and sample distribution assessment.",
        c("PCA", "Heatmaps", "Correlation"))

    # Step 5: DE
    de_label <- toupper(stats$methods$de_method)
    de_desc <- sprintf("%s differential expression test. Significance: p-value &le; %s and |linear FC| &ge; %s.",
        de_label, stats$methods$p_cutoff, stats$methods$fc_cutoff)
    de_stats <- build_metab_stats_row(
        list(value = format(stats$de$total, big.mark = ","), label = "DE features", color = "var(--accent-violet)"),
        list(value = format(stats$de$up, big.mark = ","), label = "up", color = "var(--accent-green)"),
        list(value = format(stats$de$down, big.mark = ","), label = "down", color = "var(--accent-rose)"))
    de_tags <- c(de_label, sprintf("p \u2264 %s", stats$methods$p_cutoff),
                 sprintf("|FC| \u2265 %s", stats$methods$fc_cutoff))
    step5 <- build_metab_step_html(5, "\U0001F4CA", "Differential Expression", "phase-de", "analysis",
        de_desc, de_tags, de_stats)

    # Step 6: Feature Selection
    fs_parts <- c()
    if (isTRUE(stats$feature_flags$rf_enabled))
        fs_parts <- c(fs_parts, sprintf("Random Forest (top %d features)", stats$feature_selection$rf_top_n))
    if (isTRUE(stats$feature_flags$plsda_enabled))
        fs_parts <- c(fs_parts, sprintf("PLS-DA VIP (top %d features)", stats$feature_selection$plsda_top_n))
    fs_desc <- if (length(fs_parts) > 0) paste(fs_parts, collapse = ". ") |> paste0(".")
               else "Feature selection not enabled."
    fs_tags <- c()
    if (isTRUE(stats$feature_flags$rf_enabled)) fs_tags <- c(fs_tags, "Random Forest")
    if (isTRUE(stats$feature_flags$plsda_enabled)) fs_tags <- c(fs_tags, "PLS-DA")
    if (length(fs_tags) == 0) fs_tags <- "skipped"
    step6 <- build_metab_step_html(6, "\u2B50", "Feature Selection", "phase-fs", "analysis",
        fs_desc, fs_tags)

    # Step 7: Pathway Enrichment
    enrich_desc <- if (isTRUE(stats$feature_flags$enrichment_enabled))
        sprintf("Pathway enrichment analysis performed. %d significant pathways detected (adj. p &lt; 0.05).",
                stats$enrichment$total)
    else "Pathway enrichment analysis not enabled for this run."
    enrich_tags <- if (isTRUE(stats$feature_flags$enrichment_enabled))
        c("QEA", "ssGSEA", "ORA", "GSEA", "padj < 0.05") else c("skipped")
    step7 <- build_metab_step_html(7, "\U0001F9ED", "Pathway Enrichment", "phase-enrich", "analysis",
        enrich_desc, enrich_tags)

    # Step 8: Report
    report_tag <- '<a href="report_metabolomics.html" class="tag tag-link">HTML report &rarr;</a>'
    step8 <- build_metab_step_html(8, "\U0001F4C1", "Report &amp; Outputs", "phase-out", "output",
        "Interactive HTML report with QC, PCA, volcano plots, feature importance, and enrichment results. Excel workbooks with DE results and normalized data.",
        c("report_metabolomics.html", "DE tables", "Normalized data"))
    step8 <- sub("(</div>\\s*</div>\\s*</div>)$",
                 paste0("\n            ", report_tag, "\\1"), step8)

    paste0(s_proc, "\n", step1, "\n\n", step2, "\n\n", step3, "\n\n",
           s_stat, "\n", step4, "\n\n", step5, "\n\n", step6, "\n\n", step7, "\n\n", step8)
}


# ==============================================================================
# HTML DOCUMENT ASSEMBLY
# ==============================================================================

#' @noRd
get_metab_pipeline_summary_css <- function() {
'  :root {
    --bg: #0a0e17; --surface: #111827; --surface-2: #1a2234; --border: #1e2d45;
    --text: #e2e8f0; --text-dim: #7a8ba8;
    --accent-wet: #38bdf8; --accent-green: #34d399; --accent-amber: #fbbf24;
    --accent-rose: #fb7185; --accent-violet: #a78bfa; --accent-cyan: #22d3ee;
    --accent-orange: #fb923c; --accent-blue: #60a5fa;
  }
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body { background: var(--bg); color: var(--text); font-family: "DM Sans", sans-serif; min-height: 100vh; overflow-x: hidden; }
  .bg-grid { position: fixed; inset: 0; background-image: linear-gradient(rgba(52,211,153,0.03) 1px, transparent 1px), linear-gradient(90deg, rgba(52,211,153,0.03) 1px, transparent 1px); background-size: 48px 48px; pointer-events: none; z-index: 0; }
  .bg-glow { position: fixed; width: 600px; height: 600px; border-radius: 50%; filter: blur(120px); opacity: 0.15; pointer-events: none; z-index: 0; }
  .bg-glow-1 { top: -200px; left: -100px; background: var(--accent-green); }
  .bg-glow-2 { bottom: -200px; right: -100px; background: var(--accent-violet); }
  .container { position: relative; z-index: 1; max-width: 1100px; margin: 0 auto; padding: 48px 24px 80px; }
  .header { text-align: center; margin-bottom: 56px; }
  .header-badge { display: inline-flex; align-items: center; gap: 8px; background: rgba(52,211,153,0.08); border: 1px solid rgba(52,211,153,0.2); border-radius: 100px; padding: 6px 18px; font-family: "JetBrains Mono", monospace; font-size: 11px; letter-spacing: 1.5px; text-transform: uppercase; color: var(--accent-green); margin-bottom: 20px; }
  .header-badge .dot { width: 6px; height: 6px; background: var(--accent-green); border-radius: 50%; animation: pulse 2s ease-in-out infinite; }
  @keyframes pulse { 0%,100% { opacity:1; box-shadow: 0 0 0 0 rgba(52,211,153,0.5); } 50% { opacity:0.7; box-shadow: 0 0 0 6px rgba(52,211,153,0); } }
  h1 { font-size: 36px; font-weight: 700; letter-spacing: -0.5px; background: linear-gradient(135deg, #e2e8f0, #34d399); -webkit-background-clip: text; -webkit-text-fill-color: transparent; margin-bottom: 10px; }
  .header-sub { color: var(--text-dim); font-size: 15px; }
  .header-author { color: var(--text-dim); font-size: 13px; margin-top: 6px; font-family: "JetBrains Mono", monospace; }
  .header-meta { display: flex; justify-content: center; gap: 28px; margin-top: 18px; font-family: "JetBrains Mono", monospace; font-size: 12px; color: var(--text-dim); }
  .header-meta span { display: flex; align-items: center; gap: 6px; }
  .header-meta .val { color: var(--accent-green); font-weight: 500; }
  .pipeline { position: relative; display: flex; flex-direction: column; gap: 0; }
  .pipeline::before { content: ""; position: absolute; left: 40px; top: 40px; bottom: 40px; width: 2px; background: linear-gradient(to bottom, var(--accent-green), var(--accent-amber), var(--accent-violet), var(--accent-green)); opacity: 0.3; }
  .step { display: flex; gap: 24px; align-items: flex-start; padding: 16px 0; position: relative; opacity: 0; transform: translateY(20px); animation: fadeUp 0.5s ease forwards; }
  .step:nth-child(1){animation-delay:.1s} .step:nth-child(2){animation-delay:.2s} .step:nth-child(3){animation-delay:.3s}
  .step:nth-child(4){animation-delay:.4s} .step:nth-child(5){animation-delay:.5s} .step:nth-child(6){animation-delay:.6s}
  .step:nth-child(7){animation-delay:.7s} .step:nth-child(8){animation-delay:.8s}
  @keyframes fadeUp { to { opacity: 1; transform: translateY(0); } }
  .node { flex-shrink: 0; width: 80px; display: flex; flex-direction: column; align-items: center; gap: 6px; }
  .node-circle { width: 44px; height: 44px; border-radius: 14px; display: flex; align-items: center; justify-content: center; font-size: 18px; position: relative; z-index: 2; border: 2px solid; transition: transform 0.2s, box-shadow 0.2s; }
  .step:hover .node-circle { transform: scale(1.1); }
  .node-num { font-family: "JetBrains Mono", monospace; font-size: 10px; font-weight: 600; color: var(--text-dim); }
  .card { flex: 1; background: var(--surface); border: 1px solid var(--border); border-radius: 16px; padding: 22px 26px; transition: border-color 0.2s, box-shadow 0.2s, transform 0.2s; }
  .step:hover .card { transform: translateX(4px); border-color: rgba(52,211,153,0.25); box-shadow: 0 4px 30px rgba(0,0,0,0.3); }
  .card-title { font-size: 16px; font-weight: 600; margin-bottom: 6px; display: flex; align-items: center; gap: 10px; }
  .card-phase { font-family: "JetBrains Mono", monospace; font-size: 9px; text-transform: uppercase; letter-spacing: 1.2px; padding: 2px 8px; border-radius: 4px; font-weight: 500; }
  .card-desc { color: var(--text-dim); font-size: 13.5px; line-height: 1.6; margin-bottom: 12px; }
  .card-tags { display: flex; flex-wrap: wrap; gap: 6px; }
  .tag { font-family: "JetBrains Mono", monospace; font-size: 11px; padding: 3px 10px; border-radius: 6px; background: rgba(255,255,255,0.04); border: 1px solid rgba(255,255,255,0.06); color: var(--text-dim); white-space: nowrap; }
  .phase-input .node-circle { background: rgba(34,211,238,0.1); border-color: rgba(34,211,238,0.35); color: var(--accent-cyan); }
  .phase-input .card-phase { background: rgba(34,211,238,0.1); color: var(--accent-cyan); }
  .phase-filter .node-circle { background: rgba(251,146,60,0.1); border-color: rgba(251,146,60,0.35); color: var(--accent-orange); }
  .phase-filter .card-phase { background: rgba(251,146,60,0.1); color: var(--accent-orange); }
  .phase-norm .node-circle { background: rgba(251,113,133,0.1); border-color: rgba(251,113,133,0.35); color: var(--accent-rose); }
  .phase-norm .card-phase { background: rgba(251,113,133,0.1); color: var(--accent-rose); }
  .phase-qc .node-circle { background: rgba(52,211,153,0.1); border-color: rgba(52,211,153,0.35); color: var(--accent-green); }
  .phase-qc .card-phase { background: rgba(52,211,153,0.1); color: var(--accent-green); }
  .phase-de .node-circle { background: rgba(167,139,250,0.1); border-color: rgba(167,139,250,0.35); color: var(--accent-violet); }
  .phase-de .card-phase { background: rgba(167,139,250,0.1); color: var(--accent-violet); }
  .phase-fs .node-circle { background: rgba(251,191,36,0.1); border-color: rgba(251,191,36,0.35); color: var(--accent-amber); }
  .phase-fs .card-phase { background: rgba(251,191,36,0.1); color: var(--accent-amber); }
  .phase-enrich .node-circle { background: rgba(96,165,250,0.1); border-color: rgba(96,165,250,0.35); color: var(--accent-blue); }
  .phase-enrich .card-phase { background: rgba(96,165,250,0.1); color: var(--accent-blue); }
  .phase-out .node-circle { background: rgba(52,211,153,0.1); border-color: rgba(52,211,153,0.35); color: var(--accent-green); }
  .phase-out .card-phase { background: rgba(52,211,153,0.1); color: var(--accent-green); }
  .stats-row { display: flex; gap: 10px; margin-top: 10px; flex-wrap: wrap; }
  .stat { display: flex; align-items: baseline; gap: 5px; font-family: "JetBrains Mono", monospace; font-size: 12px; }
  .stat-val { font-weight: 600; font-size: 16px; }
  .stat-label { color: var(--text-dim); font-size: 11px; }
  .section-label { display: flex; align-items: center; gap: 12px; padding: 20px 0 8px 96px; font-family: "JetBrains Mono", monospace; font-size: 10px; text-transform: uppercase; letter-spacing: 2px; color: var(--text-dim); opacity: 0; animation: fadeUp 0.5s ease forwards; }
  .section-label .team-badge { font-size: 9px; letter-spacing: 1px; padding: 2px 10px; border-radius: 100px; border: 1px solid; text-transform: none; font-weight: 500; }
  .team-multiomics { color: var(--accent-violet); border-color: rgba(167,139,250,0.3) !important; background: rgba(167,139,250,0.08); }
  .section-label::after { content: ""; flex: 1; height: 1px; background: linear-gradient(to right, var(--border), transparent); }
  .section-label:nth-of-type(1){animation-delay:.05s} .section-label:nth-of-type(2){animation-delay:.35s}
  .header-cta { margin-top: 22px; }
  .header-cta a { display: inline-flex; align-items: center; gap: 8px; padding: 10px 24px; border-radius: 100px; background: rgba(52,211,153,0.1); border: 1px solid rgba(52,211,153,0.3); color: var(--accent-green); font-family: "JetBrains Mono", monospace; font-size: 12px; font-weight: 500; letter-spacing: 0.5px; text-decoration: none; transition: background 0.2s, border-color 0.2s, box-shadow 0.2s; }
  .header-cta a:hover { background: rgba(52,211,153,0.18); border-color: rgba(52,211,153,0.5); box-shadow: 0 0 20px rgba(52,211,153,0.15); }
  .header-cta a .arrow { transition: transform 0.2s; }
  .header-cta a:hover .arrow { transform: translateX(3px); }
  .tag-link { text-decoration: none; color: var(--accent-green); border-color: rgba(52,211,153,0.2) !important; transition: background 0.2s, border-color 0.2s; }
  .tag-link:hover { background: rgba(52,211,153,0.1) !important; border-color: rgba(52,211,153,0.35) !important; }
  .footer { text-align: center; margin-top: 48px; padding-top: 24px; border-top: 1px solid var(--border); font-size: 12px; color: var(--text-dim); font-family: "JetBrains Mono", monospace; }
  @media (max-width: 640px) { .container{padding:24px 12px 48px} h1{font-size:24px} .pipeline::before{left:24px} .node{width:48px} .node-circle{width:36px;height:36px;font-size:14px;border-radius:10px} .card{padding:16px} .card-title{font-size:14px} .header-meta{flex-wrap:wrap;justify-content:center;gap:12px} .section-label{padding-left:64px} }

  /* === Analysis Overview Flowchart === */
  .fc-section {
    border-top: 1px solid var(--border); margin-top: 56px; padding-top: 40px;
    opacity: 0; animation: fadeUp 0.6s ease forwards; animation-delay: 1s;
  }
  .fc-section-title {
    text-align: center; font-size: 20px; font-weight: 600; margin-bottom: 8px;
    background: linear-gradient(135deg, #e2e8f0, #34d399);
    -webkit-background-clip: text; -webkit-text-fill-color: transparent;
  }
  .fc-section-sub {
    text-align: center; font-size: 12px; color: var(--text-dim);
    font-family: "JetBrains Mono", monospace; letter-spacing: 1px;
    text-transform: uppercase; margin-bottom: 32px;
  }
  .fc-flow { display: flex; flex-direction: column; align-items: center; gap: 0; }

  .fc-node {
    width: 240px; background: var(--surface); border: 1px solid var(--border);
    border-radius: 14px; padding: 16px 18px; text-align: center;
    transition: border-color 0.2s, box-shadow 0.2s, transform 0.2s;
    position: relative;
  }
  .fc-node:hover { transform: translateY(-3px); box-shadow: 0 8px 32px rgba(0,0,0,0.4); }
  .fc-node-icon {
    width: 38px; height: 38px; border-radius: 50%; display: flex;
    align-items: center; justify-content: center; font-size: 16px;
    margin: 0 auto 8px; border: 2px solid;
  }
  .fc-node-title {
    font-size: 13px; font-weight: 600; margin-bottom: 4px; color: var(--text);
  }
  .fc-node-detail {
    font-family: "JetBrains Mono", monospace; font-size: 10px;
    color: var(--text-dim); line-height: 1.5;
  }

  /* Connector arrow (vertical) */
  .fc-connector {
    width: 2px; height: 28px; position: relative;
    background: linear-gradient(to bottom, rgba(52,211,153,0.3), rgba(167,139,250,0.3));
  }
  .fc-connector::after {
    content: ""; position: absolute; bottom: -5px; left: 50%;
    transform: translateX(-50%);
    border-left: 5px solid transparent; border-right: 5px solid transparent;
    border-top: 6px solid rgba(167,139,250,0.4);
  }

  /* T-split connector (one line down, then horizontal bar, then drop lines) */
  .fc-t-split { display: flex; flex-direction: column; align-items: center; }
  .fc-t-stem {
    width: 2px; height: 16px;
    background: linear-gradient(to bottom, rgba(52,211,153,0.3), rgba(167,139,250,0.3));
  }
  .fc-branch-drops { display: flex; justify-content: center; position: relative; }
  .fc-branch-drops .fc-t-bar {
    position: absolute; top: 0; height: 2px;
    background: linear-gradient(to right, rgba(52,211,153,0.3), rgba(167,139,250,0.3));
  }
  .fc-branch-drops .fc-t-bar {
    left: 50%; transform: translateX(-50%);
  }
  .fc-cols-2 .fc-t-bar { width: calc(180px + 24px); }
  .fc-cols-3 .fc-t-bar { width: calc(360px + 48px); }
  .fc-drop-line {
    width: 2px; height: 16px;
    background: linear-gradient(to bottom, rgba(52,211,153,0.3), rgba(167,139,250,0.3));
    position: relative;
  }
  .fc-drop-line::after {
    content: ""; position: absolute; bottom: -5px; left: 50%;
    transform: translateX(-50%);
    border-left: 4px solid transparent; border-right: 4px solid transparent;
    border-top: 5px solid rgba(167,139,250,0.4);
  }

  /* Branch row (grid of nodes side by side) */
  .fc-branch-row { display: grid; gap: 24px; justify-content: center; }
  .fc-branch-row .fc-node { width: 180px; }
  .fc-branch-row.fc-cols-2 { grid-template-columns: 180px 180px; }
  .fc-branch-row.fc-cols-3 { grid-template-columns: 180px 180px 180px; }
  .fc-drop-wrap { display: flex; flex-direction: column; align-items: center; }

  /* Merge connector (reverse T: rise lines + horizontal bar + stem down) */
  .fc-merge-rises { display: flex; justify-content: center; position: relative; }
  .fc-merge-rises .fc-t-bar {
    position: absolute; bottom: 0; height: 2px;
    background: linear-gradient(to right, rgba(52,211,153,0.3), rgba(167,139,250,0.3));
    left: 50%; transform: translateX(-50%);
  }
  .fc-merge-rises.fc-cols-2 .fc-t-bar { width: calc(180px + 24px); }
  .fc-merge-rises.fc-cols-3 .fc-t-bar { width: calc(360px + 48px); }
  .fc-rise-line {
    width: 2px; height: 16px;
    background: linear-gradient(to top, rgba(52,211,153,0.3), rgba(167,139,250,0.3));
  }
  .fc-rise-wrap { display: flex; flex-direction: column; align-items: center; }
  .fc-merge-stem {
    width: 2px; height: 16px;
    background: linear-gradient(to bottom, rgba(52,211,153,0.3), rgba(167,139,250,0.3));
    position: relative;
  }
  .fc-merge-stem::after {
    content: ""; position: absolute; bottom: -5px; left: 50%;
    transform: translateX(-50%);
    border-left: 5px solid transparent; border-right: 5px solid transparent;
    border-top: 6px solid rgba(167,139,250,0.4);
  }

  /* Color variants for flowchart nodes */
  .fc-node-input { border-color: rgba(34,211,238,0.35); }
  .fc-node-input:hover { border-color: rgba(34,211,238,0.6); box-shadow: 0 8px 32px rgba(34,211,238,0.15); }
  .fc-node-input .fc-node-icon { background: rgba(34,211,238,0.1); border-color: rgba(34,211,238,0.4); color: var(--accent-cyan); }
  .fc-node-filter { border-color: rgba(251,146,60,0.35); }
  .fc-node-filter:hover { border-color: rgba(251,146,60,0.6); box-shadow: 0 8px 32px rgba(251,146,60,0.15); }
  .fc-node-filter .fc-node-icon { background: rgba(251,146,60,0.1); border-color: rgba(251,146,60,0.4); color: var(--accent-orange); }
  .fc-node-norm { border-color: rgba(251,113,133,0.35); }
  .fc-node-norm:hover { border-color: rgba(251,113,133,0.6); box-shadow: 0 8px 32px rgba(251,113,133,0.15); }
  .fc-node-norm .fc-node-icon { background: rgba(251,113,133,0.1); border-color: rgba(251,113,133,0.4); color: var(--accent-rose); }
  .fc-node-qc { border-color: rgba(52,211,153,0.35); }
  .fc-node-qc:hover { border-color: rgba(52,211,153,0.6); box-shadow: 0 8px 32px rgba(52,211,153,0.15); }
  .fc-node-qc .fc-node-icon { background: rgba(52,211,153,0.1); border-color: rgba(52,211,153,0.4); color: var(--accent-green); }
  .fc-node-de { border-color: rgba(167,139,250,0.35); }
  .fc-node-de:hover { border-color: rgba(167,139,250,0.6); box-shadow: 0 8px 32px rgba(167,139,250,0.15); }
  .fc-node-de .fc-node-icon { background: rgba(167,139,250,0.1); border-color: rgba(167,139,250,0.4); color: var(--accent-violet); }
  .fc-node-fs { border-color: rgba(251,191,36,0.35); }
  .fc-node-fs:hover { border-color: rgba(251,191,36,0.6); box-shadow: 0 8px 32px rgba(251,191,36,0.15); }
  .fc-node-fs .fc-node-icon { background: rgba(251,191,36,0.1); border-color: rgba(251,191,36,0.4); color: var(--accent-amber); }
  .fc-node-enrich { border-color: rgba(96,165,250,0.35); }
  .fc-node-enrich:hover { border-color: rgba(96,165,250,0.6); box-shadow: 0 8px 32px rgba(96,165,250,0.15); }
  .fc-node-enrich .fc-node-icon { background: rgba(96,165,250,0.1); border-color: rgba(96,165,250,0.4); color: var(--accent-blue); }
  .fc-node-report { border-color: rgba(52,211,153,0.35); }
  .fc-node-report:hover { border-color: rgba(52,211,153,0.6); box-shadow: 0 8px 32px rgba(52,211,153,0.15); }
  .fc-node-report .fc-node-icon { background: rgba(52,211,153,0.1); border-color: rgba(52,211,153,0.4); color: var(--accent-green); }

  @media (max-width: 640px) {
    .fc-branch-row { grid-template-columns: 1fr !important; gap: 12px; }
    .fc-branch-row .fc-node { width: 240px; margin: 0 auto; }
    .fc-branch-drops .fc-t-bar, .fc-merge-rises .fc-t-bar { display: none; }
    .fc-branch-drops, .fc-merge-rises { flex-direction: column; align-items: center; gap: 0; }
    .fc-drop-wrap, .fc-rise-wrap { width: auto; }
    .fc-node { width: 260px; }
  }'
}

#' Wrap body in full HTML document
#' @noRd
wrap_metab_html_document <- function(body_html, stats) {
    css <- get_metab_pipeline_summary_css()
    esc <- function(x) if (requireNamespace("htmltools", quietly = TRUE)) htmltools::htmlEscape(x) else x

    sprintf(
'<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Metabolomics Pipeline Workflow &mdash; %s</title>
<link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@300;400;500;600&family=DM+Sans:wght@300;400;500;600;700&display=swap" rel="stylesheet">
<style>
%s
</style>
</head>
<body>
  <div class="bg-grid"></div>
  <div class="bg-glow bg-glow-1"></div>
  <div class="bg-glow bg-glow-2"></div>
  <div class="container">
    <div class="header">
      <div class="header-badge"><span class="dot"></span> multiomics-core pipeline</div>
      <h1>Metabolomics Analysis Workflow</h1>
      <div class="header-sub">%s</div>
      <div class="header-author">Bioinformatician: %s &nbsp;&middot;&nbsp; %s</div>
      <div class="header-meta">
        <span>Samples <span class="val">%s</span></span>
        <span>Groups <span class="val">%s</span></span>
        <span>Features <span class="val">%s</span></span>
        <span>DE Features <span class="val">%s</span></span>
      </div>
      <div class="header-cta">
        <a href="report_metabolomics.html">View Full Interactive Report <span class="arrow">&rarr;</span></a>
      </div>
    </div>
    <div class="pipeline">
%s
    </div>
%s
    <div class="footer">%s &middot; %s</div>
  </div>
</body>
</html>',
        esc(stats$project$name), css, stats$subtitle,
        esc(stats$project$analyst), stats$project$date,
        stats$overview$n_samples, stats$overview$n_groups,
        stats$overview$n_features_fmt, stats$overview$n_de_features_fmt,
        body_html,
        generate_metab_flowchart_html(stats),
        esc(stats$project$analyst), stats$project$date_iso
    )
}


# ==============================================================================
# MAIN ENTRY POINT
# ==============================================================================

#' Generate the metabolomics pipeline summary HTML report
#'
#' @param config          Full pipeline config
#' @param pre             Preprocessed data list
#' @param de_res          DE results list
#' @param feature_sel_res Feature selection results
#' @param enrichment_res  Enrichment results list
#' @param run_dir         Run output directory
#' @return Path to metabolomics_pipeline_summary.html
#' @export
generate_metab_pipeline_summary <- function(config, pre, de_res,
                                             feature_sel_res, enrichment_res,
                                             run_dir) {

    message("=== Generating Metabolomics Pipeline Summary ===")

    # 1. Collect stats
    stats <- collect_metab_pipeline_stats(config, pre, de_res,
                                          feature_sel_res, enrichment_res)

    # 2. Try Claude Code CLI, fall back to R template
    metab_cfg   <- config$modes$metabolomics
    comment_cfg <- metab_cfg$commentary %||% list()
    body_html   <- NULL

    if (isTRUE(comment_cfg$enabled) && (comment_cfg$backend %||% "none") == "claude-code") {
        body_html <- generate_metab_summary_with_claude(stats, comment_cfg)
    }

    if (is.null(body_html)) {
        message("  Using R-generated metabolomics pipeline summary")
        body_html <- generate_metab_summary_body_r(stats)
        full_html <- wrap_metab_html_document(body_html, stats)
    } else {
        full_html <- body_html
    }

    # 3. Save directly to Results root (avoids duplicate inside metabolomics/)
    results_root <- dirname(run_dir)
    out_file <- file.path(results_root, "metabolomics_pipeline_summary.html")
    writeLines(full_html, out_file, useBytes = TRUE)
    message("Metabolomics pipeline summary saved to: ", out_file)

    out_file
}
