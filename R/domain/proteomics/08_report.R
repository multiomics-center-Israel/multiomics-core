#' Proteomics Report Rendering
#'
#' Copies the canonical Rmd template into the run directory and renders
#' a self-contained HTML report.

#' Render the proteomics analysis report
#'
#' @param run_dir  The results run directory (e.g. outputs/project/Results_...)
#' @param config   Full pipeline config list
#' @param config_file Path to the original YAML config file. Retained for the
#'   pipeline call sites; the config snapshot is serialized from \code{config}
#'   rather than copied from this path, so that it matches what the
#'   \code{execution_info_files} target writes.
#' @param report_type Type of report: "detailed" (default) or "short"
#' @return Path to the rendered HTML file (character, format = "file")
#' @export
render_proteomics_report <- function(run_dir, config, config_file = NULL, report_type = "detailed") {

    # Locate the canonical Rmd template shipped with multiomics-core
    template_path <- system.file("report_template_proteomics.Rmd",
                                  package = "multiomics.core",
                                  mustWork = FALSE)

    # Fallback: template lives alongside this source file
    if (!nzchar(template_path) || !file.exists(template_path)) {
        src_dir <- system.file("R", "domain", "proteomics",
                                package = "multiomics.core",
                                mustWork = FALSE)
        if (!nzchar(src_dir)) {
            # When running via targets::tar_source(), try the working directory first,
            # then fall back to the project root from config
            src_dir <- file.path("R", "domain", "proteomics")
            if (!file.exists(file.path(src_dir, "report_template_proteomics.Rmd"))) {
                proj_dir <- config$project$dir %||% "."
                src_dir <- file.path(proj_dir, "R", "domain", "proteomics")
            }
        }
        template_path <- file.path(src_dir, "report_template_proteomics.Rmd")
    }

    if (!file.exists(template_path)) {
        warning("Proteomics report template not found at: ", template_path,
                ". Skipping report generation.")
        return(NA_character_)
    }


    # Ensure run_dir exists.  The Rmd and its HTML output live in the *parent*
    # results directory (alongside the pptx / pipeline summary), not inside the
    # proteomics/ sub-folder.
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    parent_dir <- if (basename(run_dir) == "proteomics") dirname(run_dir) else run_dir
    dest_rmd <- file.path(parent_dir, "report_proteomics.Rmd")
    file.copy(template_path, dest_rmd, overwrite = TRUE)

    # Write execution_info/config_used.yaml (needed by the template). The template
    # looks for it relative to its own location, so keep a copy in the parent
    # results dir as well as the proteomics/ subdir. Rewritten on every render
    # rather than only when absent: the proteomics/ copy is written by nothing
    # else, so the old guard left the first run's snapshot in place forever and
    # config edits appeared to do nothing.
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
    out_html <- file.path(parent_dir, "report_proteomics.html")

    # Ensure any stale report / stub from a previous run is gone so the
    # post-render existence check is meaningful (not fooled by leftovers).
    if (file.exists(out_html)) try(file.remove(out_html), silent = TRUE)

    # On Windows, `Rscript --vanilla` strips RSTUDIO_PANDOC and the user's PATH
    # may not include pandoc. Probe common install locations and register the
    # first hit with rmarkdown so the render succeeds.
    ensure_pandoc_available <- function() {
        if (isTRUE(rmarkdown::pandoc_available(version = "1.12.3"))) return(TRUE)

        # Repo-bundled location (written by tools/install_pandoc.R) always comes first.
        proj_dir <- config$project$dir %||% getwd()
        candidates <- file.path(proj_dir, "tools", "pandoc")
        if (.Platform$OS.type == "windows") {
            pf   <- Sys.getenv("ProgramFiles",       "C:/Program Files")
            pf86 <- Sys.getenv("ProgramFiles(x86)",  "C:/Program Files (x86)")
            lad  <- Sys.getenv("LOCALAPPDATA",       file.path(Sys.getenv("USERPROFILE"), "AppData/Local"))
            candidates <- c(
                candidates,
                file.path(pf,   "RStudio", "resources", "app", "bin", "quarto", "bin", "tools"),
                file.path(pf,   "RStudio", "resources", "app", "bin", "pandoc"),
                file.path(pf,   "RStudio", "bin", "quarto", "bin", "tools"),
                file.path(pf,   "RStudio", "bin", "pandoc"),
                file.path(pf86, "RStudio", "bin", "pandoc"),
                file.path(pf,   "Pandoc"),
                file.path(pf86, "Pandoc"),
                file.path(lad,  "Pandoc")
            )
        } else {
            candidates <- c(candidates, "/usr/bin", "/usr/local/bin", "/opt/homebrew/bin")
        }
        candidates <- candidates[dir.exists(candidates)]

        for (d in candidates) {
            exe <- if (.Platform$OS.type == "windows") "pandoc.exe" else "pandoc"
            if (file.exists(file.path(d, exe))) {
                tryCatch({
                    rmarkdown::find_pandoc(dir = d)
                    if (isTRUE(rmarkdown::pandoc_available(version = "1.12.3"))) {
                        message("Using pandoc from: ", d)
                        return(TRUE)
                    }
                }, error = function(e) NULL)
            }
        }
        FALSE
    }

    if (!ensure_pandoc_available()) {
        warning("pandoc >= 1.12.3 not found. Install pandoc or RStudio Desktop, ",
                "or set RSTUDIO_PANDOC to a folder containing pandoc.exe.")
    }

    message("Rendering proteomics report (", report_type, ") to: ", out_html)

    render_err      <- NULL
    render_err_full <- NULL
    tryCatch({
        rmarkdown::render(
            input       = dest_rmd,
            output_file = out_html,

            output_dir  = parent_dir,
            quiet       = FALSE,  # surface inner knitr errors to the pipeline log
            envir       = new.env(parent = globalenv())
        )
    }, error = function(e) {
        render_err      <<- e$message
        render_err_full <<- paste(c(e$message,
                                     capture.output(traceback(max.lines = 40))),
                                   collapse = "\n")
        # Print to the live pipeline log so the user can see WHY it failed.
        message("\n==== Proteomics Rmd render ERROR ====")
        message(e$message)
        message("====================================\n")
    })


    # Fallback: if the render did not produce an output file (a setup chunk or
    # the pandoc step crashed), write a minimal HTML stub so downstream targets
    # and the user-facing run flow still see *something* rather than failing.
    if (!file.exists(out_html)) {
        warning("Proteomics report rendering did not produce output file.",
                if (!is.null(render_err)) paste0(" Cause: ", render_err) else "")

        # Drop the full error + traceback next to the stub so the user can send it.
        err_log <- file.path(parent_dir, "report_proteomics_error.log")
        writeLines(render_err_full %||% (render_err %||% "(error not captured)"),
                   err_log)
        message("Full render error + traceback saved to: ", err_log)

        # Diagnose common root causes for a friendlier stub message.
        err_msg <- render_err %||% "(not captured)"
        hint <- if (grepl("pandoc", err_msg, ignore.case = TRUE)) {
            paste0("<p><strong>Likely fix:</strong> pandoc is not installed or not on PATH. ",
                   "Install <a href='https://pandoc.org/installing.html'>Pandoc</a>, ",
                   "or install RStudio Desktop (which bundles pandoc), then re-run the pipeline.</p>")
        } else if (grepl("HTTP|network|connection|Ensembl|biomart", err_msg, ignore.case = TRUE)) {
            paste0("<p><strong>Likely fix:</strong> a network-dependent step (biomaRt / STRING / KEGG) could not reach the internet. ",
                   "Check connectivity or disable pathway/PPI and re-run.</p>")
        } else {
            paste0("<p>Full traceback saved to <code>report_proteomics_error.log</code> in this folder.</p>")
        }

        stub <- paste0(
            "<!doctype html><html><head><meta charset='utf-8'>",
            "<title>Proteomics Report (partial)</title>",
            "<style>body{font-family:system-ui,sans-serif;max-width:780px;margin:40px auto;padding:0 20px;color:#222}",
            "h1{color:#b91c1c}code{background:#f3f4f6;padding:2px 6px;border-radius:3px;font-size:0.95em}",
            "pre{background:#0d1220;color:#e7eefc;padding:14px;border-radius:6px;overflow-x:auto;font-size:13px;white-space:pre-wrap}",
            ".box{background:#fef3c7;border-left:4px solid #d97706;padding:14px 18px;border-radius:4px;margin:18px 0}",
            "a{color:#2563eb}",
            "</style></head><body>",
            "<h1>Proteomics report did not render</h1>",
            "<div class='box'>The pipeline finished, but the HTML report template failed to render. ",
            "All upstream analysis outputs (CSV tables, plots, PPTX, payloads) are still available under ",
            "the results directory.</div>",
            hint,
            "<p><strong>Render error:</strong></p>",
            "<pre>", htmltools::htmlEscape(err_msg), "</pre>",
            "</body></html>"
        )
        writeLines(stub, out_html)
        message("Wrote fallback stub report: ", out_html)
    } else {
        message("Proteomics report rendered successfully: ", out_html)
    }

    out_html
}

#' PCA feature-set panels for the report's protein-subset dropdown
#'
#' Lists the panel images for the report's feature-set selector, in the order
#' the dropdown shows them: all proteins, then the top-variable sets from
#' largest to smallest. Kept out of the template so the order and labels are
#' tested against the code the report actually runs.
#'
#' The complete-case panel (\code{PCA_robust.png}) is deliberately not listed.
#' It answers a different question from "how many proteins, ranked by variance",
#' and the report shows it as its own section beside the main PCA rather than as
#' an entry here — so enabling this selector cannot draw it a second time.
#'
#' @param diag_dir Directory holding the proteomics diagnostic plots.
#' @return Data frame with columns \code{key} (dropdown value), \code{path} and
#'   \code{label}, one row per panel whose image exists; zero rows if none do.
list_pca_feature_panels <- function(diag_dir) {
    panels <- data.frame(key = character(0), path = character(0),
                         label = character(0), stringsAsFactors = FALSE)
    add <- function(panels, key, path, label) {
        if (!file.exists(path)) return(panels)
        rbind(panels, data.frame(key = key, path = path, label = label,
                                 stringsAsFactors = FALSE))
    }

    # All proteins is the main PC1-vs-PC2 plot, not a panel of its own: it is
    # already the PCA of the full matrix, so the QC module writes it once and
    # this maps it to the dropdown's "all" entry.
    panels <- add(panels, "all", file.path(diag_dir, "PCA_PC1.vs.PC2.png"), "All proteins")

    top_files <- list.files(diag_dir, pattern = "^PCA_top[0-9]+\\.png$", full.names = TRUE)
    n_top <- as.numeric(sub("^PCA_top([0-9]+)\\.png$", "\\1", basename(top_files)))
    for (i in order(n_top, decreasing = TRUE)) {
        panels <- add(panels, sprintf("top%d", n_top[i]), top_files[i],
                      sprintf("Top %s variable proteins", format(n_top[i], big.mark = ",")))
    }

    panels
}

#' Sample-subset PCA panels for the report's Subsets section
#'
#' Only files named \code{PCA_subset_<name>.png} count, so the dropdown's
#' feature-set panels (all, top-N, complete-case) and the PC1/PC3 plots are
#' never listed as sample subsets.
#'
#' @param diag_dir Directory holding the proteomics diagnostic plots.
#' @return Data frame with columns \code{name} (the subset name taken from the
#'   filename) and \code{path}, sorted by name; zero rows if there are none.
list_pca_subset_panels <- function(diag_dir) {
    paths <- sort(list.files(diag_dir, pattern = "^PCA_subset_.+\\.png$", full.names = TRUE))
    data.frame(name = sub("^PCA_subset_(.+)\\.png$", "\\1", basename(paths)),
               path = paths, stringsAsFactors = FALSE)
}

# =============================================================================
# Reader-facing Methods text
# =============================================================================
# Kept out of the template for the same reason list_pca_feature_panels() is:
# the wording is then tested against the code the report actually runs, rather
# than living as prose nobody can assert on.
#
# Two rules hold throughout, and the tests pin both:
#
#   * Nothing is stated for a step that did not run. Every sentence below is
#     reached only when the configuration it describes is active.
#   * No default is introduced that disagrees with runtime. imputation$method is
#     read with the dispatcher's own default (impute_proteomics() in
#     03_imputation.R), never with "none" -- a config omitting the key imputes,
#     so a Methods section defaulting to "none" would describe a run that did
#     not happen. The other readers in this repo still default "none"; that
#     divergence is tracked separately and is deliberately not changed here.

#' Describe the min-count filter, which may be per group
#'
#' \code{extract_min_count()} accepts a bare number, a \code{default} plus
#' per-group overrides, or per-group values alone. The phrase has to survive all
#' three without claiming a single scalar threshold.
#'
#' Overrides are described as \emph{configured}, not applied:
#' \code{extract_min_count()} keeps an override only when its name is among the
#' analysed groups (\code{02_filtering.R:185-191}), so a group removed by
#' \code{sample_filter} has a configured threshold that never took effect. This
#' formatter sees the config alone and cannot tell the two apart.
#'
#' @param min_cfg The \code{filtering$min_count} config value, or NULL.
#' @return A character scalar naming the threshold(s).
methods_min_count_phrase <- function(min_cfg) {
    if (is.null(min_cfg)) return("at least 1 sample")
    if (is.numeric(min_cfg) && length(min_cfg) == 1 && is.null(names(min_cfg))) {
        return(sprintf("at least %s sample(s)", format(min_cfg)))
    }
    default_v <- min_cfg$default
    overrides <- min_cfg[setdiff(names(min_cfg), "default")]
    ov_txt <- if (length(overrides) > 0) {
        paste(sprintf("%s: %s", names(overrides), unlist(overrides)), collapse = ", ")
    } else ""
    if (!is.null(default_v) && nzchar(ov_txt)) {
        sprintf("at least %s sample(s), with configured per-group overrides (%s)",
                format(default_v), ov_txt)
    } else if (!is.null(default_v)) {
        sprintf("at least %s sample(s)", format(default_v))
    } else if (nzchar(ov_txt)) {
        sprintf("the configured per-group thresholds %s", ov_txt)
    } else {
        "at least 1 sample"
    }
}

#' The configured imputation method, resolved as the dispatcher resolves it
#'
#' @param imp_cfg The \code{imputation} config block, or NULL.
#' @return Lower-case method name; "perseus" is folded to "perseus_like" and an
#'   absent key resolves to the dispatcher's default, not to "none".
methods_imputation_method <- function(imp_cfg) {
    m <- tolower(as.character(imp_cfg$method %||% "perseus_like"))
    if (identical(m, "perseus")) "perseus_like" else m
}

#' Whether the configured imputation repeats identically
#'
#' make_imputations_proteomics() seeds each run separately, but for a
#' deterministic method every run returns the same matrix. Saying "independent"
#' of those would be false, so the consensus text branches on this.
#'
#' @param imp_cfg The \code{imputation} config block.
#' @return TRUE when repeated runs are identical by construction.
methods_imputation_is_deterministic <- function(imp_cfg) {
    m <- methods_imputation_method(imp_cfg)
    if (m %in% c("none", "minval")) return(TRUE)
    if (identical(m, "dep2")) {
        return(!identical(tolower(as.character(imp_cfg$dep2_method %||% "MinDet")), "minprob"))
    }
    FALSE
}

#' One sentence describing how missing values were filled
#'
#' @param imp_cfg The \code{imputation} config block.
#' @return A character scalar, or NULL when the method is "none".
methods_imputation_phrase <- function(imp_cfg) {
    m <- methods_imputation_method(imp_cfg)
    if (identical(m, "none")) return(NULL)
    if (identical(m, "perseus_like")) {
        return(sprintf(paste0(
            "Missing values were imputed by drawing, independently per sample, from a normal ",
            "distribution shifted down by %s and narrowed to %s of that sample's observed ",
            "standard deviation."),
            format(imp_cfg$downshift %||% 1.8), format(imp_cfg$width %||% 0.3)))
    }
    if (identical(m, "dep2")) {
        if (identical(tolower(as.character(imp_cfg$dep2_method %||% "MinDet")), "minprob")) {
            return(paste0("Missing values were drawn from a narrow normal distribution centred ",
                          "on the 1st percentile of each sample's observed values."))
        }
        return(paste0("Missing values were replaced by the 1st percentile of each sample's ",
                      "observed values."))
    }
    if (identical(m, "qrilc")) {
        return(paste0("Missing values were imputed by quantile regression imputation of ",
                      "left-censored data (QRILC)."))
    }
    if (identical(m, "minval")) {
        return("Missing values were replaced by a per-sample minimum-based value.")
    }
    sprintf("Missing values were imputed using the %s method.", m)
}

#' One sentence describing the fitted differential-abundance model
#'
#' Only for internally fitted modes. Precomputed DE is handled by the caller,
#' which must not describe a model this pipeline did not fit.
#'
#' @param de_cfg The \code{de} config block.
#' @param method The method \code{mod_proteomics_de()} reported.
#' @return A character scalar.
methods_de_phrase <- function(de_cfg, method) {
    m <- tolower(as.character(method))
    if (identical(m, "limma")) {
        s <- paste0("Moderated t-tests were fitted with the *limma* package, using empirical ",
                    "Bayes shrinkage of the per-protein variance estimates.")
        bc <- de_cfg$block_col
        if (!is.null(bc) && nzchar(bc)) {
            # run_limma_proteomics() refits without blocking when the consensus is
            # non-finite (05_de_summary.R:342-346), and that outcome never reaches
            # the report. The sentence therefore states both branches rather than
            # claiming the estimate was incorporated.
            s <- paste0(s, sprintf(paste0(
                " When %s was configured, within-block correlation was estimated with ",
                "limma::duplicateCorrelation(). A finite estimate was incorporated in the ",
                "linear-model fit; if the estimate was non-finite, the model was fitted ",
                "without blocking."), bc))
        }
        return(s)
    }
    if (identical(m, "limma_percontrast")) {
        return(paste0("Each requested two-group comparison was fitted separately with *limma* ",
                      "and empirical Bayes moderation. Within each comparison, a protein was ",
                      "retained only where it was observed above the imputation floor in at ",
                      "least one replicate of either group."))
    }
    if (m %in% c("ttest", "welch")) {
        pcol <- de_cfg$pairing_col
        if (isTRUE(de_cfg$paired) && !is.null(pcol) && nzchar(pcol)) {
            return(sprintf(paste0("Paired t-tests were applied per protein, with samples ",
                                  "matched on %s."), pcol))
        }
        return(sprintf("Two-sample t-tests assuming %s variances were applied per protein.",
                       if (identical(m, "ttest")) "equal" else "unequal"))
    }
    if (identical(m, "anova")) {
        return(paste0("A one-way ANOVA across all groups was fitted per protein, and the ",
                      "requested pairwise comparisons were obtained from Tukey's honest ",
                      "significant difference test."))
    }
    sprintf("Differential abundance was assessed using the %s method.", method)
}

#' Reader-facing Methods blocks for the proteomics report
#'
#' Builds every Methods paragraph from the configuration and from the DE mode
#' that actually ran, so the section describes the run rather than the defaults.
#'
#' The top-level branch is \code{de_method == "precomputed"}:
#' \code{mod_proteomics_de()} returns that before any internal imputation or
#' model fitting, so in that mode no internal test, no fdrtool correction and no
#' multi-imputation consensus may be described.
#'
#' @param config Full pipeline config, as written to \code{config_used.yaml}.
#' @param de_method The method \code{mod_proteomics_de()} reported.
#' @param n_included Proteins carried into the differential analysis, or NA to
#'   omit that sentence.
#' @param run_dir Run directory, read only for \code{execution_info/git_commit.txt}.
#' @param pca_cc_file Path to the complete-case PCA image, or NULL. The
#'   sensitivity-view sentence is emitted only when this file exists, because
#'   the panel is skipped on several legitimate runs.
#' @param hier_file Path to the hierarchical clustering heatmap, or NULL. The
#'   clustering sentence needs this as well as the enable flags, for the same
#'   reason.
#' @return Named list of character scalars: \code{data_processing},
#'   \code{missing_values}, \code{differential}, \code{consensus},
#'   \code{quality_control}, \code{downstream}, \code{software},
#'   \code{results_oneliner} and \code{results_thresholds}. Blocks that do not
#'   apply are NULL.
build_proteomics_methods_text <- function(config, de_method = "limma",
                                          n_included = NA_integer_, run_dir = NULL,
                                          pca_cc_file = NULL, hier_file = NULL) {
    cfg <- config$modes$proteomics %||% list()
    de_cfg <- cfg$de %||% list()
    imp_cfg <- cfg$imputation %||% list()
    filt_cfg <- cfg$filtering %||% list()

    is_precomputed <- identical(tolower(as.character(de_method)), "precomputed")
    imp_method <- methods_imputation_method(imp_cfg)
    p_adjust <- de_cfg$p_adjust_method %||% "BH"
    p_adjust_none <- identical(tolower(as.character(p_adjust)), "none")
    p_cut <- de_cfg$p_cutoff %||% 0.05
    fc_cut <- de_cfg$linear_fc_cutoff %||% 1.5

    # Resolved before the differential block so the use_adj_for_pass1 sentence
    # can be stated once: the consensus paragraph already names which p-value
    # the per-run call used, so repeating it above would say the same thing
    # twice in adjacent paragraphs.
    n_reps <- suppressWarnings(as.integer(imp_cfg$no_repetitions %||% NA))
    multi_on <- isTRUE(imp_cfg$multi_imputation %||% TRUE) &&
                !is.na(n_reps) && n_reps > 1L && !is_precomputed

    # ---- Data processing ----------------------------------------------------
    dp <- sprintf("Proteins were quantified using %s.", cfg$engine %||% "the configured search engine")
    # Resolved as get_proteomics_expression_matrix() resolves it, including the
    # branch: files$is_logtransformed is consulted only on the preprocessed path
    # (01_expression.R:122-134). The DIA-NN branch ignores that legacy flag and
    # defaults an absent scale_in to "linear" before taking log2
    # (01_expression.R:152-155), so the fallback must not be applied there.
    is_preprocessed <- identical(tolower(as.character(cfg$input$format %||% "")),
                                 "preprocessed")
    scale_in <- cfg$scale_in %||% (
        if (is_preprocessed && isTRUE(cfg$files$is_logtransformed)) "log2" else "linear")
    if (identical(tolower(as.character(scale_in)), "linear")) {
        dp <- paste(dp, "Linear-scale intensities were transformed to the log2 scale before",
                    "downstream analysis.")
    }
    # Stated only when the flag is on: filter_contaminants() honours
    # remove_contaminants and skips the step when it is FALSE.
    if (isTRUE(filt_cfg$remove_contaminants %||% TRUE)) {
        dp <- paste(dp, sprintf(
            "Entries whose identifier begins with %s were removed as contaminants.",
            filt_cfg$contaminant_prefix %||% "cRAP-"))
    }
    dp <- paste(dp, sprintf("Proteins were retained when observed in %s in at least %s group(s).",
                            methods_min_count_phrase(filt_cfg$min_count),
                            format(filt_cfg$min_groups %||% 1)))

    norm_m <- tolower(as.character(cfg$normalization$method %||% "none"))
    dp <- paste(dp, if (identical(norm_m, "median")) {
        "Sample intensities were median-centred on the log2 scale."
    } else {
        "No between-sample normalization was applied."
    })

    bc_m <- tolower(as.character(cfg$batch_correction$method %||% "none"))
    # get_proteomics_batch_config() resolves enabled as bc$enabled, defaulting to
    # method != "none". An explicit enabled: false with a method set therefore
    # runs no correction at all, so the method alone cannot gate this sentence.
    bc_default_on <- !identical(bc_m, "none")
    bc_on <- isTRUE(cfg$batch_correction$enabled %||% bc_default_on)
    if (bc_on && !identical(bc_m, "none")) {
        bc_name <- if (identical(bc_m, "combat")) "ComBat" else
                   if (identical(bc_m, "probatch")) "proBatch" else bc_m
        dp <- paste(dp, sprintf(paste0(
            "Batch effects were corrected with %s. The correction was estimated and applied on ",
            "the single-imputed complete matrix."), bc_name))
        # mod_proteomics_de() returns the loaded tables before
        # make_imputations_proteomics() is reached, so in precomputed mode
        # nothing is re-imputed for the model -- and saying otherwise would
        # contradict the missing-values paragraph below.
        dp <- paste(dp, if (is_precomputed) {
            paste("The original missingness pattern was then restored to the filtered matrix.",
                  "Quality control and clustering use the corrected complete matrix.")
        } else {
            paste("The original missingness pattern was then restored to the filtered matrix,",
                  "which is re-imputed for differential analysis, while quality-control and",
                  "clustering use the corrected complete matrix.")
        })
    }

    # ---- Missing values -----------------------------------------------------
    mv <- methods_imputation_phrase(imp_cfg)
    if (!is.null(mv)) {
        if (is_precomputed) {
            mv <- paste(mv, paste0(
                "In precomputed-DE mode, pipeline imputation is not used to fit the ",
                "differential-abundance model; it may still be used for quality-control ",
                "visualisation and, when enabled, batch correction."))
        }
        mv <- paste(mv, "Where measured and model-input values are both presented, they are",
                    "labelled separately.")
    }

    # ---- Differential protein abundance -------------------------------------
    if (is_precomputed) {
        # No internal model, and no re-adjustment: the loader takes the supplied
        # adjusted p-value and only falls back to BH -- specifically BH, not
        # de$p_adjust_method -- when the table carries none.
        de_txt <- paste0(
            "Differential-abundance statistics were loaded from precomputed result tables ",
            "supplied with the project configuration. The upstream statistical model is not ",
            "inferred by this pipeline and is not described here. ",
            "The loader uses the supplied adjusted p-value when available; if only raw ",
            "p-values are provided, adjusted p-values are calculated using the ",
            "Benjamini-Hochberg procedure. ",
            sprintf(paste0("Proteins were reported as differentially abundant where the ",
                           "adjusted p-value fell below %s and |linear fold change| was at ",
                           "least %s."), p_cut, fc_cut))
    } else {
        m <- tolower(as.character(de_method))
        de_txt <- methods_de_phrase(de_cfg, de_method)
        # fdrtool is implemented in run_limma_proteomics() and
        # run_limma_percontrast_proteomics() only; the flag being TRUE does not
        # mean a t-test or ANOVA run applied it.
        if (isTRUE(de_cfg$fdrtool_correction) && m %in% c("limma", "limma_percontrast")) {
            de_txt <- paste(de_txt, paste0(
                "Before adjustment, p-values were re-estimated from the moderated t-statistics ",
                "against an empirical null distribution (*fdrtool*)."))
        }
        # ANOVA states its own adjustment, because the values being adjusted are
        # Tukey's rather than the raw per-protein p-values.
        #
        # de$p_adjust_method accepts "none" (90_config_validate.R:99-101), and
        # p.adjust(method = "none") returns the input untouched. Calling that
        # "adjusted with the none procedure" would describe a correction that
        # did not happen, and the cutoff sentence would then call raw p-values
        # adjusted -- so both branch on it.
        de_txt <- paste(de_txt, if (p_adjust_none) {
            if (identical(m, "anova")) {
                "Those Tukey p-values were not further adjusted across proteins."
            } else {
                "P-values were not adjusted for multiple testing."
            }
        } else if (identical(m, "anova")) {
            sprintf(paste0("Those Tukey p-values were then adjusted across proteins using the ",
                           "%s procedure."), p_adjust)
        } else {
            sprintf("P-values were adjusted with the %s procedure.", p_adjust)
        })
        de_txt <- paste(de_txt, sprintf(paste0(
            "Proteins were reported as differentially abundant at %s p <= %s and ",
            "|linear fold change| >= %s."),
            if (p_adjust_none) "unadjusted" else "adjusted", p_cut, fc_cut))
        # summarize_limma_mult_imputation() reads this with isTRUE(), so an
        # absent key means the per-run call used the raw p-value. Stated here
        # only when no consensus paragraph follows to say it.
        if (!multi_on && !isTRUE(de_cfg$use_adj_for_pass1)) {
            de_txt <- paste(de_txt, "The per-run significance call used the unadjusted p-value.")
        }
    }
    # This is a row count of the summary table, and the wording has to stay
    # neutral about what produced those rows. In precomputed mode the rows come
    # from the upstream tables, not from this pipeline's filtering. Under
    # limma_percontrast each comparison applies its own observed/floor filter
    # and the dropped proteins are re-expanded as NA purely for alignment
    # (05_de_summary.R:595-599), so a summary row is not evidence that the
    # protein was fitted in any comparison.
    if (!is.na(n_included)) {
        de_txt <- paste(de_txt, sprintf(
            if (is_precomputed) {
                "%s proteins were present in the precomputed result tables."
            } else {
                "%s proteins are reported in the differential-abundance results table."
            },
            format(n_included, big.mark = ",")))
    }

    # ---- Multiple-imputation consensus --------------------------------------
    cons <- NULL
    if (multi_on) {
        min_passed <- suppressWarnings(as.integer(imp_cfg$min_no_passed %||% NA))
        cons <- sprintf(paste0(
            "The differential analysis was repeated across %d configured imputation runs. %s ",
            "Within each run, a protein passed when its %s p-value met the cutoff and its ",
            "|fold change| met the threshold. A protein was called differentially abundant ",
            "when it passed in at least %s of the %d runs and its summarised adjusted p-value ",
            "also met the cutoff."),
            n_reps,
            if (methods_imputation_is_deterministic(imp_cfg)) {
                "These runs are identical by construction, since the imputation is deterministic."
            } else {
                "Each run was drawn under its own seed."
            },
            if (isTRUE(de_cfg$use_adj_for_pass1)) "adjusted" else "unadjusted",
            format(min_passed), n_reps)
        cons <- paste(cons, sprintf(paste0(
            "The reported p-value and adjusted p-value are the %s quantile of the per-run ",
            "values. The reported fold change is the mean of the per-run ratios on the linear ",
            "scale; the reported log2 fold change is the logarithm of that mean, not the mean ",
            "of the per-run log2 fold changes."),
            if (!is.na(min_passed) && n_reps > 0) sprintf("%g", min_passed / n_reps) else "configured"))
    }

    # ---- Quality control ----------------------------------------------------
    # Both sentences below describe figures, and a configuration flag does not
    # mean the figure exists. 01_mod_qc_pre.R skips the complete-case panel when
    # there are no missingness flags, when fewer than three proteins are
    # complete cases, or when the plot call fails; mod_proteomics_clustering()
    # returns before the hierarchical step when fewer than two DE features are
    # in the matrix. Both skip with a message rather than erroring, so the
    # artefact is the only honest witness -- the same reason the complete-case
    # chunk in the template tests for its file instead of reasoning about config.
    qc <- paste0("Principal component analysis was computed on the quality-control expression ",
                 "matrix.")
    if (!is.null(pca_cc_file) && file.exists(pca_cc_file)) {
        qc <- paste(qc, "A complete-case panel restricted to proteins observed in every",
                    "included sample is shown as a sensitivity view.")
    }
    hier <- cfg$clustering$steps$hierarchical
    # clustering_run_flags() resolves this as isTRUE(steps$hierarchical$enabled
    # %||% TRUE) (09_clustering.R:626), so a steps block with no enabled key
    # runs the step. Reading it with a bare isTRUE() here would drop the
    # sentence from a run whose heatmap is sitting right there.
    if (isTRUE(cfg$clustering$enabled) && isTRUE(hier$enabled %||% TRUE) &&
        !is.null(hier_file) && file.exists(hier_file)) {
        qc <- paste(qc, sprintf(paste0(
            "When hierarchical clustering was enabled, row-z-scored values were clustered ",
            "using %s distance and %s linkage."),
            hier$distance %||% "euclidean", hier$linkage %||% "complete"))
    }

    # ---- Downstream analyses ------------------------------------------------
    # Carried over from the previous Methods text rather than dropped: both
    # sentences are reader-facing Methods content and both were already gated on
    # their enabled flag. Only the wording is config-driven now.
    down <- character(0)
    if (isTRUE(cfg$pathway$enabled)) {
        pw <- tolower(as.character(cfg$pathway$method %||% "both"))
        pw_desc <- if (identical(pw, "fgsea")) {
            "gene-set enrichment analysis"
        } else if (identical(pw, "ora")) {
            "over-representation analysis"
        } else if (identical(pw, "both")) {
            "gene-set enrichment analysis and over-representation analysis"
        } else {
            sprintf("the %s method", pw)
        }
        down <- c(down, sprintf("Pathway enrichment was performed by %s against %s.",
                                pw_desc,
                                paste(cfg$pathway$databases %||% c("GO", "KEGG", "Reactome"),
                                      collapse = ", ")))
    }
    if (isTRUE(cfg$ppi$enabled)) {
        down <- c(down, paste0("Protein-protein interaction networks were constructed using ",
                               "the STRING database."))
    }
    down <- if (length(down) > 0) paste(down, collapse = " ") else NULL

    # ---- Software -----------------------------------------------------------
    sw <- sprintf("Analyses were performed in R %s using the multiomics-core pipeline",
                  paste(R.version$major, R.version$minor, sep = "."))
    sha <- NULL
    if (!is.null(run_dir)) {
        sha_file <- file.path(run_dir, "execution_info", "git_commit.txt")
        if (file.exists(sha_file)) {
            sha <- tryCatch(substr(trimws(readLines(sha_file, warn = FALSE)[1]), 1, 8),
                            error = function(e) NULL)
        }
    }
    # An empty or whitespace-only git_commit.txt reads back as NA here, via
    # readLines()[1] on a zero-length result. The commit is optional provenance,
    # so an unusable file is treated as absent rather than printed.
    sw <- if (!is.null(sha) && !is.na(sha) && nzchar(sha)) {
        sprintf("%s (commit %s).", sw, sha)
    } else {
        paste0(sw, ".")
    }

    # ---- Results-section pointer --------------------------------------------
    one_liner <- if (is_precomputed) {
        "Differential-abundance statistics were loaded from precomputed tables; see Methods."
    } else {
        sprintf(paste0("Differential protein abundance was assessed with %s; see Methods for ",
                       "the full description."), de_method)
    }

    # The thresholds belong with the one-liner in Results, and they are built
    # here so the comparator cannot drift from the one the differential block
    # states. The two modes genuinely differ: the internal summary calls a
    # protein significant at padj <= cutoff (05_de_summary.R:117, 209), while
    # load_precomputed_proteomics_de() uses a strict < (05_de_summary.R:895).
    thresholds <- if (is_precomputed) {
        sprintf(paste0("Proteins were reported as differentially abundant where the adjusted ",
                       "p-value fell below %s and |linear fold change| was at least %s."),
                p_cut, fc_cut)
    } else {
        sprintf(paste0("Proteins were reported as differentially abundant at %s ",
                       "p-value $\\leq$ %s and |linear fold change| $\\geq$ %s."),
                if (p_adjust_none) "unadjusted" else "adjusted", p_cut, fc_cut)
    }

    list(data_processing = dp, missing_values = mv, differential = de_txt,
         consensus = cons, quality_control = qc, downstream = down, software = sw,
         results_oneliner = one_liner, results_thresholds = thresholds)
}
