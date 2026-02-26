#' Pipeline Summary HTML Generation
#'
#' Generates a dark-themed HTML summary of the complete RNA-seq pipeline
#' workflow, from wet lab through downstream analysis. Data is extracted
#' from pipeline objects and config. Descriptions polished by Claude Code CLI.
#'
#' Design template by Dor Shalev (Azrieli Technion Genomics Center).

# ==============================================================================
# DATA COLLECTION
# ==============================================================================

#' Collect all pipeline statistics for the summary
#'
#' @param config      Full pipeline config
#' @param pre         Preprocessed data list (expr_raw, expr_filt, meta)
#' @param de_res      DE results list (summary_df, tables)
#' @param pathway_res Pathway results list
#' @return Named list with all stats organized by section
collect_pipeline_stats <- function(config, pre, de_res, pathway_res) {

    rna_cfg  <- config$modes$rna
    tech_rep <- rna_cfg$technical_report %||% list()

    # --- Project info ---
    analyst <- gsub("_", " ", config$project$analyst %||% "")

    project <- list(
        name     = config$project$name %||% "RNA-seq Analysis",
        analyst  = analyst,
        date     = format(Sys.Date(), "%B %d, %Y"),
        date_iso = format(Sys.Date(), "%Y-%m-%d"),
        facility = tech_rep$facility %||% ""
    )

    # --- Sample & gene overview ---
    n_samples   <- if (!is.null(pre$expr_filt)) ncol(pre$expr_filt) else NA
    n_genes_raw <- if (!is.null(pre$expr_raw)) nrow(pre$expr_raw) else NA
    n_genes     <- if (!is.null(pre$expr_filt)) nrow(pre$expr_filt) else NA

    group_col <- rna_cfg$filtering$group_col %||% rna_cfg$effects$color
    n_groups <- NA
    groups <- character()
    if (!is.null(pre$meta) && !is.null(group_col) && group_col %in% names(pre$meta)) {
        groups <- unique(as.character(pre$meta[[group_col]]))
        n_groups <- length(groups)
    }

    # --- DE statistics ---
    # Extract from de_res$tables (per-contrast data frames with padj, log2FoldChange)
    # or from de_res$summary_df (wide format with padj.ContrastName columns)
    n_de_total <- 0; n_de_up <- 0; n_de_down <- 0
    de_contrasts <- list()

    p_cut  <- rna_cfg$de$p_cutoff %||% 0.05
    fc_lin <- rna_cfg$de$linear_fc_cutoff %||% 1.0
    log2_fc <- if (fc_lin > 1) log2(fc_lin) else 0

    if (!is.null(de_res$tables) && length(de_res$tables) > 0) {
        # Primary path: per-contrast tables (each has padj, log2FoldChange)
        for (cn in names(de_res$tables)) {
            tbl <- de_res$tables[[cn]]
            if (!is.data.frame(tbl) || !("padj" %in% names(tbl))) next

            sig <- !is.na(tbl$padj) & tbl$padj <= p_cut
            if (log2_fc > 0 && "log2FoldChange" %in% names(tbl)) {
                sig <- sig & (abs(tbl$log2FoldChange) >= log2_fc)
            }

            up <- if ("log2FoldChange" %in% names(tbl)) sum(sig & tbl$log2FoldChange > 0, na.rm = TRUE) else NA
            dn <- if ("log2FoldChange" %in% names(tbl)) sum(sig & tbl$log2FoldChange < 0, na.rm = TRUE) else NA

            de_contrasts[[cn]] <- list(
                name = cn, total = sum(sig, na.rm = TRUE),
                up = up %||% 0, down = dn %||% 0, tested = sum(!is.na(tbl$padj))
            )
            n_de_total <- n_de_total + sum(sig, na.rm = TRUE)
            n_de_up <- n_de_up + (up %||% 0); n_de_down <- n_de_down + (dn %||% 0)
        }
    } else if (!is.null(de_res$summary_df)) {
        # Fallback: wide-format summary_df (padj.ContrastName columns)
        sdf <- de_res$summary_df
        padj_cols <- grep("^padj\\.", names(sdf), value = TRUE)

        for (pc in padj_cols) {
            cn <- sub("^padj\\.", "", pc)
            lcol <- paste0("log2FoldChange.", cn)
            if (!(lcol %in% names(sdf))) next

            sig <- !is.na(sdf[[pc]]) & sdf[[pc]] <= p_cut
            if (log2_fc > 0) sig <- sig & (abs(sdf[[lcol]]) >= log2_fc)

            up <- sum(sig & sdf[[lcol]] > 0, na.rm = TRUE)
            dn <- sum(sig & sdf[[lcol]] < 0, na.rm = TRUE)

            de_contrasts[[cn]] <- list(
                name = cn, total = sum(sig, na.rm = TRUE),
                up = up, down = dn, tested = sum(!is.na(sdf[[pc]]))
            )
            n_de_total <- n_de_total + sum(sig, na.rm = TRUE)
            n_de_up <- n_de_up + up; n_de_down <- n_de_down + dn
        }
    }

    # --- Pathway statistics ---
    # Structure: pathway_results -> ContrastName -> Database_method -> data.frame
    n_pathways_total <- 0; n_pathways_up <- 0; n_pathways_down <- 0
    if (!is.null(pathway_res$pathway_results)) {
        count_from_df <- function(df) {
            if (!is.data.frame(df) || !"padj" %in% names(df)) return(NULL)
            sig_pw <- df[!is.na(df$padj) & df$padj < 0.05, , drop = FALSE]
            n_up <- if ("NES" %in% names(sig_pw)) sum(sig_pw$NES > 0, na.rm = TRUE) else 0
            n_dn <- if ("NES" %in% names(sig_pw)) sum(sig_pw$NES < 0, na.rm = TRUE) else 0
            list(total = nrow(sig_pw), up = n_up, down = n_dn)
        }

        for (pr_name in names(pathway_res$pathway_results)) {
            pr <- pathway_res$pathway_results[[pr_name]]
            if (is.data.frame(pr)) {
                # Direct data frame (flat structure)
                res <- count_from_df(pr)
            } else if (is.list(pr)) {
                # Nested: list of database results per contrast
                res <- list(total = 0, up = 0, down = 0)
                for (sub_name in names(pr)) {
                    sub_res <- count_from_df(pr[[sub_name]])
                    if (!is.null(sub_res)) {
                        res$total <- res$total + sub_res$total
                        res$up    <- res$up + sub_res$up
                        res$down  <- res$down + sub_res$down
                    }
                }
            } else {
                res <- NULL
            }
            if (!is.null(res)) {
                n_pathways_total <- n_pathways_total + res$total
                n_pathways_up    <- n_pathways_up + res$up
                n_pathways_down  <- n_pathways_down + res$down
            }
        }
    }

    # --- Subtitle from organism + groups ---
    organism <- rna_cfg$annotation$organism %||% ""
    subtitle <- project$name
    if (nzchar(organism) && n_groups >= 2) {
        genus <- sub("\\s.*", "", organism)
        subtitle <- sprintf("%s &mdash; %s %s",
                            project$name, genus,
                            paste(groups, collapse = " vs "))
    }

    list(
        project  = project,
        subtitle = subtitle,
        overview = list(
            n_samples = n_samples %||% 0, n_groups = n_groups %||% 0,
            groups = groups, n_genes = n_genes %||% 0,
            n_genes_fmt = format(n_genes %||% 0, big.mark = ","),
            n_de_genes_fmt = format(n_de_total, big.mark = ",")
        ),
        de = list(total = n_de_total, up = n_de_up, down = n_de_down,
                  contrasts = de_contrasts),
        pathway = list(total = n_pathways_total, up = n_pathways_up,
                       down = n_pathways_down),
        methods = list(
            normalization = rna_cfg$normalization$method %||% "TMM",
            prior_count   = rna_cfg$normalization$prior.count %||% 1,
            filtering_mode = rna_cfg$filtering$mode %||% "adaptive",
            de_method = rna_cfg$de$method %||% "deseq2",
            de_mode   = rna_cfg$de$deseq_mode %||% "default",
            p_cutoff  = rna_cfg$de$p_cutoff %||% 0.05,
            fc_cutoff = rna_cfg$de$linear_fc_cutoff %||% 1.0,
            pathway_method = rna_cfg$pathway$method %||% "fgsea",
            pathway_dbs = rna_cfg$pathway$databases %||% c("GO", "KEGG"),
            organism = organism
        ),
        technical = list(
            library_prep = tech_rep$library_prep %||% "",
            sequencing   = tech_rep$sequencing %||% "",
            upstream     = tech_rep$upstream_bioinformatics %||% "",
            de_note      = tech_rep$upstream_de_note %||% ""
        ),
        filtering = list(
            genes_before = n_genes_raw %||% 0, genes_after = n_genes %||% 0,
            before_fmt = format(n_genes_raw %||% 0, big.mark = ","),
            after_fmt  = format(n_genes %||% 0, big.mark = ","),
            pct_retained = if (!is.na(n_genes_raw) && n_genes_raw > 0 && !is.na(n_genes))
                round(100 * n_genes / n_genes_raw, 1) else NA
        ),
        feature_flags = list(
            batch_enabled      = isTRUE(rna_cfg$batch_correction$enabled),
            batch_method       = rna_cfg$batch_correction$method %||% "combat_seq",
            pathway_enabled    = isTRUE(rna_cfg$pathway$enabled),
            clustering_enabled = isTRUE(rna_cfg$clustering$enabled)
        )
    )
}


# ==============================================================================
# CLAUDE CODE CLI — generate the complete HTML via `claude --print`
# ==============================================================================

#' Generate pipeline summary HTML using Claude Code CLI
#'
#' Sends the collected stats to Claude via `claude --print` and asks it to
#' generate the complete pipeline summary HTML document.
#'
#' @param stats       Pipeline stats from collect_pipeline_stats()
#' @param comment_cfg Commentary config section (for model selection)
#' @return Complete HTML document string, or NULL on failure
generate_summary_with_claude <- function(stats, comment_cfg) {

    if (Sys.which("claude") == "") {
        message("  'claude' CLI not found in PATH.")
        return(NULL)
    }

    model <- comment_cfg$claude_code_model %||% "sonnet"
    stats_json <- jsonlite::toJSON(stats, auto_unbox = TRUE, pretty = TRUE)

    # Write stats to a temp file so we don't hit shell argument limits
    tmp_stats <- tempfile(fileext = ".json")
    on.exit(unlink(tmp_stats), add = TRUE)
    writeLines(stats_json, tmp_stats)

    # Build the prompt — Claude generates the full HTML document
    prompt <- paste0(
        "Read the pipeline statistics JSON file at '", tmp_stats, "' and generate a COMPLETE, self-contained HTML document ",
        "for an RNA-seq pipeline workflow summary page.\n\n",
        "DESIGN REQUIREMENTS:\n",
        "- Dark theme: background #0a0e17, cards #111827, border #1e2d45\n",
        "- Google Fonts: JetBrains Mono (monospace elements) + DM Sans (body text)\n",
        "- Gradient h1 text (white to sky-blue), subtle background grid, two blurred glow circles\n",
        "- Vertical timeline with gradient connector line on the left\n",
        "- Each step: node circle (emoji + step number) + card (title, phase badge, description, tags)\n",
        "- Cards have hover effects (translateX + border glow + shadow)\n",
        "- Steps animate in with fadeUp (staggered delays)\n",
        "- Responsive: mobile breakpoint at 640px\n\n",
        "STRUCTURE:\n",
        "Header: badge 'multiomics-core pipeline', h1 'RNA-seq Analysis Workflow', subtitle from stats, ",
        "analyst + date, meta row (Samples, Groups, Genes, DE Genes) with actual numbers from stats.\n\n",
        "Pipeline sections:\n",
        "1. 'Wet Lab' section (team badge 'ATGC', cyan):\n",
        "   - Step 01 Library Preparation (phase-wet, blue, emoji DNA) — use stats.technical.library_prep for description\n",
        "   - Step 02 Sequencing (phase-seq, cyan, emoji lightning) — use stats.technical.sequencing\n",
        "2. 'Upstream Bioinformatics' section (team badge 'ATGC', cyan):\n",
        "   - Step 03 QC & Trimming (phase-qc, green, emoji checkmark) — extract QC/trimming info from stats.technical.upstream\n",
        "   - Step 04 Alignment & Quantification (phase-align, amber, emoji target) — extract alignment info from stats.technical.upstream\n",
        "3. 'Downstream Analysis' section (team badge 'Multi-omics Core', violet):\n",
        "   - Step 05 Gene Filtering (phase-filter, orange, emoji magnifier) — include stats-row with before/after gene counts\n",
        "   - Step 06 Normalization (phase-norm, rose, emoji scales) — describe the normalization method from stats.methods\n",
        "   - Step 07 Differential Expression (phase-de, violet, emoji chart) — include stats-row with DE/up/down counts\n",
        "   - Step 08 Pathway Enrichment (phase-enrich, blue, emoji compass) — include stats-row with pathway counts\n",
        "   - Step 09 Report & Outputs (phase-out, green, emoji folder)\n\n",
        "Footer: facility name + date.\n\n",
        "COLOR SCHEME per phase class:\n",
        "phase-wet: rgba(56,189,248,*), phase-seq: rgba(34,211,238,*), phase-qc: rgba(52,211,153,*), ",
        "phase-align: rgba(251,191,36,*), phase-filter: rgba(251,146,60,*), phase-norm: rgba(251,113,133,*), ",
        "phase-de: rgba(167,139,250,*), phase-enrich: rgba(96,165,250,*), phase-out: rgba(52,211,153,*)\n\n",
        "Tags should be tool names/versions extracted from the technical report text, or method parameters.\n",
        "Stats rows use class='stat' with 'stat-val' (large colored number) and 'stat-label' (dim text).\n\n",
        "Write SCIENTIFICALLY PRECISE descriptions. Use HTML entities for special chars.\n",
        "Output ONLY the complete HTML document — no markdown fences, no explanations."
    )

    cmd <- sprintf(
        "unset CLAUDECODE; claude --print --model %s --no-session-persistence %s",
        model, shQuote(prompt)
    )

    result <- tryCatch({
        raw <- system(cmd, intern = TRUE, timeout = 180)
        html <- paste(raw, collapse = "\n")

        # Strip markdown fences if Claude wrapped the output
        html <- gsub("^```html\\s*\n?", "", html)
        html <- gsub("\n?```\\s*$", "", html)

        # Validate: must be valid HTML with expected elements
        if (grepl("<!DOCTYPE html>", html, ignore.case = TRUE) &&
            grepl("phase-de", html) &&
            grepl("pipeline", html)) {
            message("  Claude Code generated pipeline summary HTML")
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
build_step_html <- function(num, emoji, title, phase_class, phase_label,
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
build_stats_row <- function(...) {
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
fc_node <- function(icon, title, detail, variant) {
    sprintf(
'<div class="fc-node %s">
  <div class="fc-node-icon">%s</div>
  <div class="fc-node-title">%s</div>
  <div class="fc-node-detail">%s</div>
</div>', variant, icon, title, detail)
}

#' Build vertical connector arrow
#' @noRd
fc_arrow <- function() {
    '<div class="fc-connector"></div>'
}

#' Build T-split (stem down + horizontal bar + drop lines with arrows)
#' @noRd
fc_t_split <- function(n_cols) {
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
fc_merge <- function(n_cols) {
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

#' Generate the analysis overview flowchart HTML
#'
#' Builds a visual data-flow diagram from pipeline stats and feature flags.
#' Nodes are conditionally included based on which pipeline steps were enabled.
#'
#' @param stats Pipeline stats from collect_pipeline_stats()
#' @return HTML string for the flowchart section
#' @noRd
generate_flowchart_html <- function(stats) {

    ff <- stats$feature_flags
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
    input_detail <- sprintf("%s samples &middot; %s genes",
                            stats$overview$n_samples, format(stats$filtering$genes_before, big.mark = ","))
    parts <- c(parts, list(
        fc_node("\U0001F4E5", "Input Data", input_detail, "fc-node-input"),
        fc_arrow()
    ))

    # --- Gene Filtering ---
    filt_detail <- sprintf("%s &rarr; %s genes (%.1f%% retained)",
                           stats$filtering$before_fmt, stats$filtering$after_fmt,
                           stats$filtering$pct_retained %||% 0)
    parts <- c(parts, list(
        fc_node("\U0001F50D", "Gene Filtering", filt_detail, "fc-node-filter"),
        fc_arrow()
    ))

    # --- Normalization ---
    norm_detail <- esc(stats$methods$normalization)
    parts <- c(parts, list(
        fc_node("\u2696\uFE0F", "Normalization", norm_detail, "fc-node-norm")
    ))

    # --- Batch Correction (optional) ---
    if (isTRUE(ff$batch_enabled)) {
        batch_detail <- esc(ff$batch_method)
        parts <- c(parts, list(
            fc_arrow(),
            fc_node("\U0001F504", "Batch Correction", batch_detail, "fc-node-batch")
        ))
    }

    # --- QC / DE split (always 2 branches) ---
    de_method <- toupper(stats$methods$de_method)
    de_detail <- sprintf("%s &middot; %s DE genes", de_method,
                         format(stats$de$total, big.mark = ","))

    parts <- c(parts, list(
        fc_t_split(2),
        sprintf('<div class="fc-branch-row fc-cols-2">%s%s</div>',
                fc_node("\u2705", "Quality Control", "PCA &middot; Heatmaps", "fc-node-qc"),
                fc_node("\U0001F4CA", "Differential Expression", de_detail, "fc-node-de")),
        fc_merge(2)
    ))

    # --- Downstream branches: Volcano (always), Pathway (if enabled), Clustering (if enabled) ---
    downstream_nodes <- list()
    downstream_nodes <- c(downstream_nodes, list(
        fc_node("\U0001F30B", "Volcano &amp; MA Plots", "Visualization", "fc-node-viz")
    ))
    if (isTRUE(ff$pathway_enabled)) {
        pw_detail <- sprintf("%s &middot; %s pathways",
                             toupper(stats$methods$pathway_method),
                             stats$pathway$total)
        downstream_nodes <- c(downstream_nodes, list(
            fc_node("\U0001F9ED", "Pathway Enrichment", pw_detail, "fc-node-pathway")
        ))
    }
    if (isTRUE(ff$clustering_enabled)) {
        downstream_nodes <- c(downstream_nodes, list(
            fc_node("\U0001F4A0", "Clustering", "Gene clusters", "fc-node-cluster")
        ))
    }

    n_downstream <- length(downstream_nodes)

    if (n_downstream == 1) {
        # Single node — no split needed, show inline
        parts <- c(parts, list(
            fc_arrow(),
            downstream_nodes[[1]]
        ))
    } else {
        # Multi-branch split
        cols_cls <- sprintf("fc-cols-%d", n_downstream)
        row_html <- sprintf('<div class="fc-branch-row %s">%s</div>',
                            cols_cls, paste(downstream_nodes, collapse = ""))
        parts <- c(parts, list(
            fc_t_split(n_downstream),
            row_html,
            fc_merge(n_downstream)
        ))
    }

    # --- Report & Export ---
    parts <- c(parts, list(
        fc_arrow(),
        fc_node("\U0001F4C1", "Report &amp; Export", "HTML &middot; Excel &middot; CSV", "fc-node-report")
    ))

    # Close section
    parts <- c(parts, list(
        '  </div>',
        '</div>'
    ))

    paste(parts, collapse = "\n")
}

#' Extract tool tags from text
#' @noRd
extract_tags_from_text <- function(text, pattern = NULL, fallback = character()) {
    if (is.null(text) || !nzchar(text)) return(fallback)
    if (is.null(pattern))
        pattern <- "([A-Za-z][A-Za-z0-9_-]+)\\s*[v(]+\\s*v?([0-9]+\\.[0-9.]+)"
    matches <- gregexpr(pattern, text, perl = TRUE)
    extracted <- regmatches(text, matches)[[1]]
    if (length(extracted) == 0) return(fallback)
    tags <- gsub("\\(v?", "v", extracted)
    tags <- gsub("[)]", "", tags)
    tags <- gsub("\\s+v", " v", tags)
    if (length(tags) == 0) fallback else tags
}

#' Generate pipeline body HTML from stats (R fallback)
#' @noRd
generate_summary_body_r <- function(stats) {

    esc <- function(x) if (requireNamespace("htmltools", quietly = TRUE)) htmltools::htmlEscape(x) else x

    # -- Wet Lab --
    s_wet <- '      <div class="section-label">Wet Lab <span class="team-badge team-atgc">ATGC</span></div>\n'

    lib_desc <- if (nzchar(stats$technical$library_prep)) esc(stats$technical$library_prep)
                else sprintf("%d RNA-seq libraries constructed. QC performed prior to pooling.", stats$overview$n_samples)
    lib_tags <- extract_tags_from_text(stats$technical$library_prep,
                                        fallback = c(sprintf("%d libraries", stats$overview$n_samples)))
    step1 <- build_step_html(1, "\U0001F9EC", "Library Preparation", "phase-wet", "wet lab", lib_desc, lib_tags)

    seq_desc <- if (nzchar(stats$technical$sequencing)) esc(stats$technical$sequencing) else "Illumina sequencing performed."
    seq_tags <- extract_tags_from_text(stats$technical$sequencing, fallback = c("Illumina"))
    step2 <- build_step_html(2, "\u26A1", "Sequencing", "phase-seq", "sequencing", seq_desc, seq_tags)

    # -- Upstream --
    s_up <- '      <div class="section-label">Upstream Bioinformatics <span class="team-badge team-atgc">ATGC</span></div>\n'

    up_text <- stats$technical$upstream
    qc_desc <- ""; align_desc <- ""
    if (nzchar(up_text)) {
        parts <- strsplit(up_text, "(?=Reference:|Alignment:)", perl = TRUE)[[1]]
        if (length(parts) >= 2) { qc_desc <- esc(trimws(parts[1])); align_desc <- esc(trimws(paste(parts[-1], collapse = " "))) }
        else { qc_desc <- esc(up_text); align_desc <- "Reads aligned to reference genome." }
    } else { qc_desc <- "Read quality assessed and adapters trimmed."; align_desc <- "Reads aligned and quantified." }

    qc_tags <- extract_tags_from_text(up_text,
        pattern = "(FastQC|MultiQC|Cutadapt|CUTADAPT|Trimmomatic)\\s*[v(]+\\s*v?([0-9]+\\.[0-9.]+)",
        fallback = c("QC", "Trimming"))
    step3 <- build_step_html(3, "\u2713", "QC &amp; Trimming", "phase-qc", "upstream", qc_desc, qc_tags)

    align_tags <- extract_tags_from_text(up_text,
        pattern = "(STAR|HISAT2|Salmon|kallisto|gffread|featureCounts)\\s*[v(]+\\s*v?([0-9]+\\.[0-9.]+)",
        fallback = c("Alignment"))
    if (nzchar(stats$methods$organism)) align_tags <- c(align_tags, stats$methods$organism)
    step4 <- build_step_html(4, "\U0001F3AF", "Alignment &amp; Quantification", "phase-align", "upstream", align_desc, align_tags)

    # -- Downstream --
    s_down <- '      <div class="section-label">Downstream Analysis <span class="team-badge team-multiomics">Multi-omics Core</span></div>\n'

    # Filtering
    filt_desc <- switch(stats$methods$filtering_mode,
        "adaptive" = "Adaptive density-based threshold (KDE on log2-CPM) selects optimal CPM cutoff. Genes passing in at least one biological group are retained.",
        "deseq2_only" = "DESeq2 independent filtering applied. Low-count genes removed by DESeq2&rsquo;s automatic threshold.",
        "fixed" = "Fixed CPM threshold applied. Genes above threshold in at least one group retained.",
        "Gene filtering applied to remove lowly expressed genes.")
    filt_stats <- ""
    if (stats$filtering$genes_before > 0 && stats$filtering$genes_after > 0) {
        filt_stats <- build_stats_row(
            list(value = stats$filtering$before_fmt, label = "before", color = "var(--accent-orange)"),
            list(arrow = TRUE),
            list(value = stats$filtering$after_fmt, label = "after", color = "var(--accent-green)"),
            list(value = "", label = sprintf("%.1f%% retained", stats$filtering$pct_retained %||% 0), margin = "8px", color = "var(--text)"))
    }
    step5 <- build_step_html(5, "\U0001F50D", "Gene Filtering", "phase-filter", "downstream", filt_desc, c(stats$methods$filtering_mode), filt_stats)

    # Normalization
    nm <- stats$methods$normalization
    norm_desc <- switch(nm,
        "TMMlogCPM" = sprintf("TMM scaling factors (edgeR) correct library composition, then log2 CPM (prior count = %s) for stabilized values.", stats$methods$prior_count),
        "VST" = "Variance Stabilizing Transformation (DESeq2) stabilizes variance across expression range.",
        sprintf("Normalization via %s.", nm))
    norm_tags <- if (nm == "TMMlogCPM") c("TMM scaling", "log2 CPM", sprintf("prior.count = %s", stats$methods$prior_count))
                 else if (nm == "VST") c("VST", "DESeq2") else nm
    step6 <- build_step_html(6, "\u2696", "Normalization", "phase-norm", "downstream", norm_desc, norm_tags)

    # DE
    de_label <- toupper(stats$methods$de_method)
    de_desc <- sprintf("%s %s test (%s mode). Significance: adj. p-value &le; %s and |linear FC| &ge; %s. FDR via Benjamini-Hochberg.",
        de_label, if (stats$methods$de_mode == "default") "Wald" else stats$methods$de_mode,
        stats$methods$de_mode, stats$methods$p_cutoff, stats$methods$fc_cutoff)
    de_stats <- build_stats_row(
        list(value = format(stats$de$total, big.mark = ","), label = "DE genes", color = "var(--accent-violet)"),
        list(value = format(stats$de$up, big.mark = ","), label = "up", color = "var(--accent-green)"),
        list(value = format(stats$de$down, big.mark = ","), label = "down", color = "var(--accent-rose)"))
    de_tags <- c(de_label, sprintf("padj \u2264 %s", stats$methods$p_cutoff), sprintf("|FC| \u2265 %s", stats$methods$fc_cutoff))
    step7 <- build_step_html(7, "\U0001F4CA", "Differential Expression", "phase-de", "downstream", de_desc, de_tags, de_stats)

    # Pathway
    pw_dbs <- paste(stats$methods$pathway_dbs, collapse = " and ")
    pw_desc <- sprintf("%s ranks all genes by DE statistic and tests enrichment against %s. Pathways with adj. p &lt; 0.05 reported.",
        toupper(stats$methods$pathway_method), pw_dbs)
    pw_stats <- ""
    if (stats$pathway$total > 0) {
        pw_items <- list(list(value = as.character(stats$pathway$total), label = "pathways", color = "var(--accent-blue)"))
        if (stats$pathway$up > 0) pw_items <- c(pw_items, list(list(value = as.character(stats$pathway$up), label = "up", color = "var(--accent-green)")))
        if (stats$pathway$down > 0) pw_items <- c(pw_items, list(list(value = as.character(stats$pathway$down), label = "down", color = "var(--accent-rose)")))
        pw_stats <- do.call(build_stats_row, pw_items)
    }
    pw_tags <- c(toupper(stats$methods$pathway_method), stats$methods$pathway_dbs, "padj < 0.05")
    step8 <- build_step_html(8, "\U0001F9ED", "Pathway Enrichment", "phase-enrich", "downstream", pw_desc, pw_tags, pw_stats)

    # Report — "HTML report" tag is a link to the full interactive report
    report_tag <- '<a href="report_rnaseq.html" class="tag tag-link">HTML report &rarr;</a>'
    step9 <- build_step_html(9, "\U0001F4C1", "Report &amp; Outputs", "phase-out", "output",
        "Interactive HTML report with QC, PCA, heatmaps, volcano &amp; MA plots. Excel workbooks with DE results, normalized counts, and enrichment tables.",
        c("Final_results_DE.xlsx", "Normalized counts", "Enrichment CSVs"),
        stats_html = "")
    # Inject the link tag at the end of the tags div
    step9 <- sub("(</div>\\s*</div>\\s*</div>)$",
                 paste0("\n            ", report_tag, "\\1"), step9)

    paste0(s_wet, "\n", step1, "\n\n", step2, "\n\n",
           s_up, "\n", step3, "\n\n", step4, "\n\n",
           s_down, "\n", step5, "\n\n", step6, "\n\n", step7, "\n\n", step8, "\n\n", step9)
}


# ==============================================================================
# HTML DOCUMENT ASSEMBLY
# ==============================================================================

#' @noRd
get_pipeline_summary_css <- function() {
'  :root {
    --bg: #0a0e17; --surface: #111827; --surface-2: #1a2234; --border: #1e2d45;
    --text: #e2e8f0; --text-dim: #7a8ba8;
    --accent-wet: #38bdf8; --accent-green: #34d399; --accent-amber: #fbbf24;
    --accent-rose: #fb7185; --accent-violet: #a78bfa; --accent-cyan: #22d3ee;
    --accent-orange: #fb923c; --accent-blue: #60a5fa;
  }
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body { background: var(--bg); color: var(--text); font-family: "DM Sans", sans-serif; min-height: 100vh; overflow-x: hidden; }
  .bg-grid { position: fixed; inset: 0; background-image: linear-gradient(rgba(56,189,248,0.03) 1px, transparent 1px), linear-gradient(90deg, rgba(56,189,248,0.03) 1px, transparent 1px); background-size: 48px 48px; pointer-events: none; z-index: 0; }
  .bg-glow { position: fixed; width: 600px; height: 600px; border-radius: 50%; filter: blur(120px); opacity: 0.15; pointer-events: none; z-index: 0; }
  .bg-glow-1 { top: -200px; left: -100px; background: var(--accent-wet); }
  .bg-glow-2 { bottom: -200px; right: -100px; background: var(--accent-violet); }
  .container { position: relative; z-index: 1; max-width: 1100px; margin: 0 auto; padding: 48px 24px 80px; }
  .header { text-align: center; margin-bottom: 56px; }
  .header-badge { display: inline-flex; align-items: center; gap: 8px; background: rgba(56,189,248,0.08); border: 1px solid rgba(56,189,248,0.2); border-radius: 100px; padding: 6px 18px; font-family: "JetBrains Mono", monospace; font-size: 11px; letter-spacing: 1.5px; text-transform: uppercase; color: var(--accent-wet); margin-bottom: 20px; }
  .header-badge .dot { width: 6px; height: 6px; background: var(--accent-green); border-radius: 50%; animation: pulse 2s ease-in-out infinite; }
  @keyframes pulse { 0%,100% { opacity:1; box-shadow: 0 0 0 0 rgba(52,211,153,0.5); } 50% { opacity:0.7; box-shadow: 0 0 0 6px rgba(52,211,153,0); } }
  h1 { font-size: 36px; font-weight: 700; letter-spacing: -0.5px; background: linear-gradient(135deg, #e2e8f0, #38bdf8); -webkit-background-clip: text; -webkit-text-fill-color: transparent; margin-bottom: 10px; }
  .header-sub { color: var(--text-dim); font-size: 15px; }
  .header-author { color: var(--text-dim); font-size: 13px; margin-top: 6px; font-family: "JetBrains Mono", monospace; }
  .header-meta { display: flex; justify-content: center; gap: 28px; margin-top: 18px; font-family: "JetBrains Mono", monospace; font-size: 12px; color: var(--text-dim); }
  .header-meta span { display: flex; align-items: center; gap: 6px; }
  .header-meta .val { color: var(--accent-wet); font-weight: 500; }
  .pipeline { position: relative; display: flex; flex-direction: column; gap: 0; }
  .pipeline::before { content: ""; position: absolute; left: 40px; top: 40px; bottom: 40px; width: 2px; background: linear-gradient(to bottom, var(--accent-wet), var(--accent-green), var(--accent-amber), var(--accent-rose), var(--accent-violet)); opacity: 0.3; }
  .step { display: flex; gap: 24px; align-items: flex-start; padding: 16px 0; position: relative; opacity: 0; transform: translateY(20px); animation: fadeUp 0.5s ease forwards; }
  .step:nth-child(1){animation-delay:.1s} .step:nth-child(2){animation-delay:.2s} .step:nth-child(3){animation-delay:.3s}
  .step:nth-child(4){animation-delay:.4s} .step:nth-child(5){animation-delay:.5s} .step:nth-child(6){animation-delay:.6s}
  .step:nth-child(7){animation-delay:.7s} .step:nth-child(8){animation-delay:.8s} .step:nth-child(9){animation-delay:.9s}
  @keyframes fadeUp { to { opacity: 1; transform: translateY(0); } }
  .node { flex-shrink: 0; width: 80px; display: flex; flex-direction: column; align-items: center; gap: 6px; }
  .node-circle { width: 44px; height: 44px; border-radius: 14px; display: flex; align-items: center; justify-content: center; font-size: 18px; position: relative; z-index: 2; border: 2px solid; transition: transform 0.2s, box-shadow 0.2s; }
  .step:hover .node-circle { transform: scale(1.1); }
  .node-num { font-family: "JetBrains Mono", monospace; font-size: 10px; font-weight: 600; color: var(--text-dim); }
  .card { flex: 1; background: var(--surface); border: 1px solid var(--border); border-radius: 16px; padding: 22px 26px; transition: border-color 0.2s, box-shadow 0.2s, transform 0.2s; }
  .step:hover .card { transform: translateX(4px); border-color: rgba(56,189,248,0.25); box-shadow: 0 4px 30px rgba(0,0,0,0.3); }
  .card-title { font-size: 16px; font-weight: 600; margin-bottom: 6px; display: flex; align-items: center; gap: 10px; }
  .card-phase { font-family: "JetBrains Mono", monospace; font-size: 9px; text-transform: uppercase; letter-spacing: 1.2px; padding: 2px 8px; border-radius: 4px; font-weight: 500; }
  .card-desc { color: var(--text-dim); font-size: 13.5px; line-height: 1.6; margin-bottom: 12px; }
  .card-tags { display: flex; flex-wrap: wrap; gap: 6px; }
  .tag { font-family: "JetBrains Mono", monospace; font-size: 11px; padding: 3px 10px; border-radius: 6px; background: rgba(255,255,255,0.04); border: 1px solid rgba(255,255,255,0.06); color: var(--text-dim); white-space: nowrap; }
  .phase-wet .node-circle { background: rgba(56,189,248,0.1); border-color: rgba(56,189,248,0.35); color: var(--accent-wet); }
  .phase-wet .card-phase { background: rgba(56,189,248,0.1); color: var(--accent-wet); }
  .phase-seq .node-circle { background: rgba(34,211,238,0.1); border-color: rgba(34,211,238,0.35); color: var(--accent-cyan); }
  .phase-seq .card-phase { background: rgba(34,211,238,0.1); color: var(--accent-cyan); }
  .phase-qc .node-circle { background: rgba(52,211,153,0.1); border-color: rgba(52,211,153,0.35); color: var(--accent-green); }
  .phase-qc .card-phase { background: rgba(52,211,153,0.1); color: var(--accent-green); }
  .phase-align .node-circle { background: rgba(251,191,36,0.1); border-color: rgba(251,191,36,0.35); color: var(--accent-amber); }
  .phase-align .card-phase { background: rgba(251,191,36,0.1); color: var(--accent-amber); }
  .phase-filter .node-circle { background: rgba(251,146,60,0.1); border-color: rgba(251,146,60,0.35); color: var(--accent-orange); }
  .phase-filter .card-phase { background: rgba(251,146,60,0.1); color: var(--accent-orange); }
  .phase-norm .node-circle { background: rgba(251,113,133,0.1); border-color: rgba(251,113,133,0.35); color: var(--accent-rose); }
  .phase-norm .card-phase { background: rgba(251,113,133,0.1); color: var(--accent-rose); }
  .phase-de .node-circle { background: rgba(167,139,250,0.1); border-color: rgba(167,139,250,0.35); color: var(--accent-violet); }
  .phase-de .card-phase { background: rgba(167,139,250,0.1); color: var(--accent-violet); }
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
  .team-atgc { color: var(--accent-cyan); border-color: rgba(34,211,238,0.3) !important; background: rgba(34,211,238,0.08); }
  .team-multiomics { color: var(--accent-violet); border-color: rgba(167,139,250,0.3) !important; background: rgba(167,139,250,0.08); }
  .section-label::after { content: ""; flex: 1; height: 1px; background: linear-gradient(to right, var(--border), transparent); }
  .section-label:nth-of-type(1){animation-delay:.05s} .section-label:nth-of-type(2){animation-delay:.35s} .section-label:nth-of-type(3){animation-delay:.55s}
  .header-cta { margin-top: 22px; }
  .header-cta a { display: inline-flex; align-items: center; gap: 8px; padding: 10px 24px; border-radius: 100px; background: rgba(56,189,248,0.1); border: 1px solid rgba(56,189,248,0.3); color: var(--accent-wet); font-family: "JetBrains Mono", monospace; font-size: 12px; font-weight: 500; letter-spacing: 0.5px; text-decoration: none; transition: background 0.2s, border-color 0.2s, box-shadow 0.2s; }
  .header-cta a:hover { background: rgba(56,189,248,0.18); border-color: rgba(56,189,248,0.5); box-shadow: 0 0 20px rgba(56,189,248,0.15); }
  .header-cta a .arrow { transition: transform 0.2s; }
  .header-cta a:hover .arrow { transform: translateX(3px); }
  .tag-link { text-decoration: none; color: var(--accent-wet); border-color: rgba(56,189,248,0.2) !important; transition: background 0.2s, border-color 0.2s; }
  .tag-link:hover { background: rgba(56,189,248,0.1) !important; border-color: rgba(56,189,248,0.35) !important; }
  .footer { text-align: center; margin-top: 48px; padding-top: 24px; border-top: 1px solid var(--border); font-size: 12px; color: var(--text-dim); font-family: "JetBrains Mono", monospace; }
  @media (max-width: 640px) { .container{padding:24px 12px 48px} h1{font-size:24px} .pipeline::before{left:24px} .node{width:48px} .node-circle{width:36px;height:36px;font-size:14px;border-radius:10px} .card{padding:16px} .card-title{font-size:14px} .header-meta{flex-wrap:wrap;justify-content:center;gap:12px} .section-label{padding-left:64px} }

  /* === Analysis Overview Flowchart === */
  .fc-section {
    border-top: 1px solid var(--border); margin-top: 56px; padding-top: 40px;
    opacity: 0; animation: fadeUp 0.6s ease forwards; animation-delay: 1s;
  }
  .fc-section-title {
    text-align: center; font-size: 20px; font-weight: 600; margin-bottom: 8px;
    background: linear-gradient(135deg, #e2e8f0, #38bdf8);
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
    background: linear-gradient(to bottom, rgba(56,189,248,0.3), rgba(167,139,250,0.3));
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
    background: linear-gradient(to bottom, rgba(56,189,248,0.3), rgba(167,139,250,0.3));
  }
  .fc-branch-drops { display: flex; justify-content: center; position: relative; }
  .fc-branch-drops .fc-t-bar {
    position: absolute; top: 0; height: 2px;
    background: linear-gradient(to right, rgba(56,189,248,0.3), rgba(167,139,250,0.3));
  }
  .fc-cols-2 .fc-t-bar { left: 90px; right: 90px; width: calc(180px + 24px); margin: 0 auto; }
  .fc-cols-3 .fc-t-bar { left: 90px; right: 90px; width: calc(360px + 48px); margin: 0 auto; }
  .fc-branch-drops .fc-t-bar {
    left: 50%; transform: translateX(-50%);
  }
  .fc-cols-2 .fc-t-bar { width: calc(180px + 24px); }
  .fc-cols-3 .fc-t-bar { width: calc(360px + 48px); }
  .fc-drop-line {
    width: 2px; height: 16px;
    background: linear-gradient(to bottom, rgba(56,189,248,0.3), rgba(167,139,250,0.3));
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

  /* Merge connector (reverse T: rise lines + horizontal bar + stem up) */
  .fc-merge-rises { display: flex; justify-content: center; position: relative; }
  .fc-merge-rises .fc-t-bar {
    position: absolute; bottom: 0; height: 2px;
    background: linear-gradient(to right, rgba(56,189,248,0.3), rgba(167,139,250,0.3));
    left: 50%; transform: translateX(-50%);
  }
  .fc-cols-2 .fc-t-bar, .fc-merge-rises.fc-cols-2 .fc-t-bar { width: calc(180px + 24px); }
  .fc-cols-3 .fc-t-bar, .fc-merge-rises.fc-cols-3 .fc-t-bar { width: calc(360px + 48px); }
  .fc-rise-line {
    width: 2px; height: 16px;
    background: linear-gradient(to top, rgba(56,189,248,0.3), rgba(167,139,250,0.3));
  }
  .fc-rise-wrap { display: flex; flex-direction: column; align-items: center; }
  .fc-merge-stem {
    width: 2px; height: 16px;
    background: linear-gradient(to bottom, rgba(56,189,248,0.3), rgba(167,139,250,0.3));
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
  .fc-node-batch { border-color: rgba(251,191,36,0.35); }
  .fc-node-batch:hover { border-color: rgba(251,191,36,0.6); box-shadow: 0 8px 32px rgba(251,191,36,0.15); }
  .fc-node-batch .fc-node-icon { background: rgba(251,191,36,0.1); border-color: rgba(251,191,36,0.4); color: var(--accent-amber); }
  .fc-node-qc { border-color: rgba(52,211,153,0.35); }
  .fc-node-qc:hover { border-color: rgba(52,211,153,0.6); box-shadow: 0 8px 32px rgba(52,211,153,0.15); }
  .fc-node-qc .fc-node-icon { background: rgba(52,211,153,0.1); border-color: rgba(52,211,153,0.4); color: var(--accent-green); }
  .fc-node-de { border-color: rgba(167,139,250,0.35); }
  .fc-node-de:hover { border-color: rgba(167,139,250,0.6); box-shadow: 0 8px 32px rgba(167,139,250,0.15); }
  .fc-node-de .fc-node-icon { background: rgba(167,139,250,0.1); border-color: rgba(167,139,250,0.4); color: var(--accent-violet); }
  .fc-node-viz { border-color: rgba(56,189,248,0.35); }
  .fc-node-viz:hover { border-color: rgba(56,189,248,0.6); box-shadow: 0 8px 32px rgba(56,189,248,0.15); }
  .fc-node-viz .fc-node-icon { background: rgba(56,189,248,0.1); border-color: rgba(56,189,248,0.4); color: var(--accent-wet); }
  .fc-node-pathway { border-color: rgba(96,165,250,0.35); }
  .fc-node-pathway:hover { border-color: rgba(96,165,250,0.6); box-shadow: 0 8px 32px rgba(96,165,250,0.15); }
  .fc-node-pathway .fc-node-icon { background: rgba(96,165,250,0.1); border-color: rgba(96,165,250,0.4); color: var(--accent-blue); }
  .fc-node-cluster { border-color: rgba(251,191,36,0.35); }
  .fc-node-cluster:hover { border-color: rgba(251,191,36,0.6); box-shadow: 0 8px 32px rgba(251,191,36,0.15); }
  .fc-node-cluster .fc-node-icon { background: rgba(251,191,36,0.1); border-color: rgba(251,191,36,0.4); color: var(--accent-amber); }
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
wrap_html_document <- function(body_html, stats) {
    css <- get_pipeline_summary_css()
    esc <- function(x) if (requireNamespace("htmltools", quietly = TRUE)) htmltools::htmlEscape(x) else x

    sprintf(
'<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>RNA-seq Pipeline Workflow &mdash; %s</title>
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
      <h1>RNA-seq Analysis Workflow</h1>
      <div class="header-sub">%s</div>
      <div class="header-author">Bioinformatician: %s &nbsp;&middot;&nbsp; %s</div>
      <div class="header-meta">
        <span>Samples <span class="val">%s</span></span>
        <span>Groups <span class="val">%s</span></span>
        <span>Genes <span class="val">%s</span></span>
        <span>DE Genes <span class="val">%s</span></span>
      </div>
      <div class="header-cta">
        <a href="report_rnaseq.html">View Full Interactive Report <span class="arrow">&rarr;</span></a>
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
        stats$overview$n_genes_fmt, stats$overview$n_de_genes_fmt,
        body_html,
        generate_flowchart_html(stats),
        esc(stats$project$facility), stats$project$date_iso
    )
}


# ==============================================================================
# MAIN ENTRY POINT
# ==============================================================================

#' Generate the pipeline summary HTML report
#'
#' @param config      Full pipeline config
#' @param pre         Preprocessed data list
#' @param de_res      DE results list
#' @param pathway_res Pathway results list
#' @param run_dir     Run output directory
#' @return Path to pipeline_summary.html
#' @export
generate_pipeline_summary <- function(config, pre, de_res, pathway_res, run_dir) {

    message("=== Generating Pipeline Summary ===")

    # 1. Collect stats
    stats <- collect_pipeline_stats(config, pre, de_res, pathway_res)

    # 2. Try Claude Code CLI, fall back to R template
    rna_cfg     <- config$modes$rna
    comment_cfg <- rna_cfg$commentary %||% list()
    body_html   <- NULL

    if (isTRUE(comment_cfg$enabled) && (comment_cfg$backend %||% "none") == "claude-code") {
        body_html <- generate_summary_with_claude(stats, comment_cfg)
    }

    if (is.null(body_html)) {
        message("  Using R-generated pipeline summary")
        body_html <- generate_summary_body_r(stats)
        full_html <- wrap_html_document(body_html, stats)
    } else {
        # Claude returned the full document
        full_html <- body_html
    }

    # 3. Save
    out_file <- file.path(run_dir, "rna_pipeline_summary.html")
    writeLines(full_html, out_file, useBytes = TRUE)
    message("Pipeline summary saved to: ", out_file)
    out_file
}
