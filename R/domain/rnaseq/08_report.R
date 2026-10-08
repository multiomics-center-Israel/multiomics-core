#' RNA-seq Report Rendering
#'
#' Copies the canonical Rmd template into the run directory and renders
#' a self-contained HTML report.

#' Render the RNA-seq analysis report
#'
#' @param run_dir  The results run directory (e.g. outputs/project/Results_.../rna)
#' @param config   Full pipeline config list
#' @param config_file Path to the original YAML config file. Retained for the
#'   pipeline call sites; the config snapshot is serialized from \code{config}
#'   rather than copied from this path, so that it matches what the
#'   \code{execution_info_files} target writes.
#' @return Path to the rendered HTML file (character, format = "file")
#' @export
render_rnaseq_report <- function(run_dir, config, config_file = NULL) {

    # Locate the canonical Rmd template shipped with multiomics-core.
    template_path <- system.file("report_template.Rmd",
                                  package = "multiomics.core",
                                  mustWork = FALSE)

    if (!nzchar(template_path) || !file.exists(template_path)) {
        template_path <- file.path(getwd(), "R", "domain", "rnaseq",
                                   "report_template.Rmd")
    }

    if (!file.exists(template_path)) {
        warning("Report template not found at: ", template_path,
                ". Skipping report generation.")
        return(NA_character_)
    }

    # The Rmd and its HTML output live in the *parent* results directory
    # (alongside the pptx / pipeline summary), not inside the rna/ sub-folder.
    # The template already references paths as run_dir/rna/... internally.
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    parent_dir <- if (basename(run_dir) %in% c("rna", "rnaseq")) dirname(run_dir) else run_dir
    dest_rmd <- file.path(parent_dir, "report_rnaseq.Rmd")
    file.copy(template_path, dest_rmd, overwrite = TRUE)

    # Write execution_info/config_used.yaml (needed by the template), keeping a
    # copy in both parent and rna/ subdirectories. Rewritten on every render
    # rather than only when absent: the rna/ copy is written by nothing else, so
    # the old guard left the first run's snapshot in place forever and config
    # edits appeared to do nothing.
    #
    # Always serialized from `config`, never copied from the source YAML. The
    # parent path is the same file the execution_info_files target owns, and that
    # target writes it with this same yaml::write_yaml(config, ...). Copying the
    # raw YAML here would put different content at a tracked path -- load_config()
    # adds .config_path and .config_mtime, which the source file does not carry --
    # leaving the target outdated the moment the report finished.
    for (edir in unique(c(file.path(parent_dir, "execution_info"),
                          file.path(run_dir, "execution_info")))) {
        dir.create(edir, recursive = TRUE, showWarnings = FALSE)
        yaml::write_yaml(config, file.path(edir, "config_used.yaml"))
    }

    # Render into the parent results directory
    out_html <- file.path(parent_dir, "report_rnaseq.html")

    # Ensure any stale report from a previous run is gone, as the proteomics and
    # multiomics renderers do. The config snapshot above is now rewritten every
    # time, so a render that fails afterwards would otherwise leave last run's
    # HTML sitting beside this run's snapshot, reading as though that report had
    # been produced from these settings. It also makes the file.exists(out_html)
    # check below meaningful rather than satisfied by the leftover.
    if (file.exists(out_html)) try(file.remove(out_html), silent = TRUE)

    message("Rendering RNA-seq report to: ", out_html)

    tryCatch({
        rmarkdown::render(
            input       = dest_rmd,
            output_file = out_html,
            output_dir  = parent_dir,
            quiet       = TRUE,
            envir       = new.env(parent = globalenv())
        )
    }, error = function(e) {
        warning("RNA-seq report rendering failed: ", e$message,
                "\nTemplate: ", dest_rmd,
                "\nCheck the Rmd for errors.")
    })

    rmarkdown::render(
        input       = dest_rmd,
        output_file = out_html,
        output_dir  = run_dir,
        quiet       = TRUE,
        envir       = new.env(parent = globalenv())
    )

    if (file.exists(out_html)) {
        message("Report rendered successfully: ", out_html)
    } else {
        warning("Report rendering did not produce output file.")
    }

    out_html
}

#' DE gene sets per contrast, for the report's overlap plots
#'
#' Applies the report's DE rule (padj at or below the cutoff and |log2FC| at
#' or above the cutoff) to every contrast in the final results table. Kept out
#' of the template so the All / Up / Down tabs share one rule.
#'
#' @param de_data Final results table (one row per gene, columns
#'   \code{padj.<contrast>} and \code{log2FoldChange.<contrast>} or
#'   \code{linearFC.<contrast>}).
#' @param contrast_names Contrast names, in display order.
#' @param padj_cut Adjusted p-value cutoff.
#' @param lfc_cut Absolute log2 fold-change cutoff.
#' @param direction \code{"all"}, \code{"up"} or \code{"down"}.
#' @return Named list of gene vectors, one per contrast with both columns
#'   present; names have underscores replaced by spaces.
collect_de_gene_sets <- function(de_data, contrast_names, padj_cut, lfc_cut,
                                 direction = c("all", "up", "down")) {
    direction <- match.arg(direction)
    gene_col <- if ("GeneName" %in% names(de_data)) "GeneName" else names(de_data)[1]
    sets <- list()
    for (cn in contrast_names) {
        padj_col <- grep(paste0("^padj\\.", cn, "$"), names(de_data), value = TRUE)[1]
        lfc_col <- grep(paste0("^(log2FoldChange|linearFC)\\.", cn, "$"), names(de_data), value = TRUE)[1]
        if (is.na(padj_col) || is.na(lfc_col)) next
        padj_vals <- as.numeric(de_data[[padj_col]])
        lfc_vals <- as.numeric(de_data[[lfc_col]])
        if (grepl("^linearFC", lfc_col)) lfc_vals <- log2(abs(lfc_vals)) * sign(lfc_vals)
        passes_fc <- switch(direction,
            all  = abs(lfc_vals) >= lfc_cut,
            up   = lfc_vals >= lfc_cut,
            down = lfc_vals <= -lfc_cut)
        is_de <- !is.na(padj_vals) & padj_vals <= padj_cut & !is.na(passes_fc) & passes_fc
        sets[[gsub("_", " ", cn)]] <- de_data[[gene_col]][is_de]
    }
    sets
}

#' Draw the overlap of DE gene sets: a Venn diagram or an UpSet plot
#'
#' A Venn diagram stays readable up to three sets. Beyond that most of its
#' regions are tiny or empty, so an UpSet plot is drawn instead: one bar per
#' combination of contrasts that actually shares genes, largest first.
#'
#' @param gene_sets Named list of gene vectors, from
#'   \code{collect_de_gene_sets()}.
#' @param title Plot title.
#' @param fill_high Fill colour for the Venn regions and UpSet bars.
#' @param max_venn Largest number of sets drawn as a Venn diagram.
#' @param max_combinations Most combinations shown in the UpSet plot.
#' @param use_venn Whether ggVennDiagram is available.
#' @return Invisibly, \code{"venn"}, \code{"upset"} or \code{"none"}: what was drawn.
draw_de_overlap <- function(gene_sets, title, fill_high = "steelblue", max_venn = 3,
                            max_combinations = 25, use_venn = TRUE) {
    total <- length(unique(unlist(gene_sets)))
    if (length(gene_sets) < 2 || total == 0) {
        cat("No DE genes to compare across contrasts.\n")
        return(invisible("none"))
    }

    if (length(gene_sets) <= max_venn && isTRUE(use_venn)) {
        names(gene_sets) <- gsub(" vs ", "\nvs ", names(gene_sets))
        p <- ggVennDiagram::ggVennDiagram(gene_sets, label_alpha = 0, label = "both",
                                          label_percent_digit = 1, set_size = 3.5) +
            ggplot2::scale_fill_gradient(low = "white", high = fill_high) +
            ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.25)) +
            ggplot2::labs(title = title, subtitle = paste0("Total unique genes: ", total)) +
            ggplot2::theme(legend.position = "none",
                           plot.title = ggplot2::element_text(hjust = 0),
                           plot.subtitle = ggplot2::element_text(hjust = 0),
                           plot.margin = ggplot2::margin(10, 40, 10, 40))
        print(p)
        return(invisible("venn"))
    }

    m <- ComplexHeatmap::make_comb_mat(gene_sets)
    m <- m[ComplexHeatmap::comb_size(m) > 0]
    n_comb <- length(ComplexHeatmap::comb_size(m))
    keep <- utils::head(order(ComplexHeatmap::comb_size(m), decreasing = TRUE), max_combinations)
    m <- m[keep]
    ht <- ComplexHeatmap::UpSet(
        m,
        set_order = seq_along(gene_sets),
        comb_order = order(ComplexHeatmap::comb_size(m), decreasing = TRUE),
        comb_col = fill_high,
        top_annotation = ComplexHeatmap::upset_top_annotation(m, add_numbers = TRUE,
                                                              gp = grid::gpar(fill = fill_high)),
        right_annotation = ComplexHeatmap::upset_right_annotation(m, add_numbers = TRUE,
                                                                  gp = grid::gpar(fill = fill_high)),
        row_names_max_width = ComplexHeatmap::max_text_width(names(gene_sets)),
        column_title = sprintf("%s (total unique genes: %d; %s)", title, total,
                               if (n_comb <= max_combinations) "all combinations shown"
                               else sprintf("largest %d combinations shown", max_combinations))
    )
    # Left padding so long contrast names are not clipped at the device edge.
    ComplexHeatmap::draw(ht, padding = grid::unit(c(2, 8, 2, 2), "mm"))
    invisible("upset")
}

#' Write a pheatmap object to a PNG for embedding in the report
#'
#' Printing a pheatmap inside a report chunk draws on whatever device is
#' current. In the pipeline's R session that was not knitr's device, so the
#' Top DE heatmaps went to a stray Rplots.pdf and the report showed nothing.
#' Drawing onto a PNG device opened and closed here does not depend on that
#' state; the chunk then embeds the file with \code{knitr::include_graphics()}.
#'
#' @param ph A pheatmap object (\code{pheatmap(..., silent = TRUE)}).
#' @param file PNG path to write; a temporary file by default.
#' @param width,height Size in inches.
#' @param res Resolution in pixels per inch.
#' @return The PNG path, invisibly.
write_pheatmap_png <- function(ph, file = tempfile(fileext = ".png"),
                               width = 10, height = 8, res = 100) {
    grDevices::png(file, width = width, height = height, units = "in", res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
    invisible(file)
}
