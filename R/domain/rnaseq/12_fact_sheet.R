#' Assemble the RNA-seq results fact sheet
#'
#' Collects the headline numbers of a run and pairs each with the artifact it
#' can be checked against. Every value is read or derived from the files already
#' written to \code{out_dir}, so the sheet describes the results as they exist
#' on disk rather than as they existed in memory.
#'
#' Sections, in order: cohort and filtering, differential expression, sample
#' structure, gene-set enrichment, run provenance. Any section whose inputs are
#' missing is skipped with a warning rather than failing the run.
#'
#' @param out_dir The RNA-seq results directory (the one holding Datasets/).
#' @param config Full pipeline config.
#' @param pre Optional preprocessing result; used for the per-group library
#'   counts, which cannot be recovered from the matrices alone.
#' @param inputs Optional loaded inputs; used for the contrast names.
#' @param run_dir Optional run root, for the execution_info provenance rows.
#' @return A data.frame with columns claim, value, source_file.
build_rnaseq_fact_sheet <- function(out_dir, config, pre = NULL, inputs = NULL,
                                    run_dir = NULL) {
    dirs <- create_legacy_output_dirs(out_dir, create = FALSE)
    de_cfg <- config$modes$rna$de %||% list()
    # Exactly the resolution the DE code uses (R/domain/rnaseq/04_de_summary.R:100,148).
    # A legacy config carrying only padj_cutoff would otherwise have its calls
    # made at that value while the sheet claimed 0.05 -- the one kind of error a
    # fact sheet must not make.
    p_cut <- de_cfg$p_cutoff %||% de_cfg$padj_cutoff %||% 0.05
    fc_cut <- de_cfg$linear_fc_cutoff %||% 1.5

    counts <- read_artifact_or_null(file.path(dirs$datasets, "rna_counts_filtered.tsv"))
    final <- read_artifact_or_null(file.path(dirs$datasets, "final_results.tsv"))
    de_counts <- .read_de_summary_counts(out_dir, dirs)
    norm_path <- list.files(dirs$datasets, pattern = "^rna_norm_.*\\.tsv$", full.names = TRUE)[1]
    norm <- read_artifact_or_null(norm_path)
    contrast_names <- if (!is.null(inputs$contrasts)) as.character(inputs$contrasts$Contrast_name) else character(0)

    bind_facts(
        .facts_rnaseq_cohort(counts, final, pre, inputs, config, p_cut, fc_cut),
        .facts_rnaseq_de(final, de_counts, de_cfg),
        .facts_rnaseq_structure(norm, norm_path, pre, config),
        .facts_rnaseq_enrichment(dirs$enrichment, p_cut, contrast_names),
        .facts_run_provenance(run_dir)
    )
}

#' Read the DE summary counts, whichever schema this run wrote
#'
#' Two writers produce a file of this name with different headers: the mode-root
#' copy uses contrast/up/down/total, the Datasets copy uses Name/up/down/any.
#' Rather than guess, take the first that parses and normalise the names.
#'
#' The Datasets copy is tried first on purpose. It is written by
#' \code{rna_outputs_legacy}, which this target depends on, so it is always this
#' run's. The mode-root copy comes from the Shiny export, which runs
#' independently: on a rerun into an existing output directory it can still hold
#' the previous run's counts, and the sheet would then record them as current.
#'
#' @param out_dir The RNA-seq results directory.
#' @param dirs Result of \code{create_legacy_output_dirs()}.
#' @return A data.frame with columns contrast, up, down, total, or NULL.
#' @keywords internal
.read_de_summary_counts <- function(out_dir, dirs) {
    for (p in c(file.path(dirs$datasets, "de_summary_counts.tsv"),
                file.path(out_dir, "de_summary_counts.tsv"))) {
        tab <- read_artifact_or_null(p)
        if (is.null(tab)) next
        names(tab)[names(tab) == "Name"] <- "contrast"
        names(tab)[names(tab) == "any"] <- "total"
        if (all(c("contrast", "up", "down", "total") %in% names(tab))) {
            # both writers append an all-contrasts row under a different label
            return(tab[!tab$contrast %in% c("any", "pass_any"), , drop = FALSE])
        }
    }
    NULL
}

#' Cohort, filtering and analysis settings rows
#' @keywords internal
.facts_rnaseq_cohort <- function(counts, final, pre, inputs, config, p_cut, fc_cut) {
    sample_cols <- if (!is.null(counts)) setdiff(names(counts), names(counts)[1]) else NULL

    group_desc <- NULL
    if (!is.null(pre$meta) && !is.null(inputs$contrasts) && "Factor" %in% names(inputs$contrasts)) {
        fcol <- unique(as.character(inputs$contrasts$Factor))[1]
        if (!is.na(fcol) && fcol %in% names(pre$meta)) {
            tab <- table(as.character(pre$meta[[fcol]]))
            group_desc <- paste(sprintf("%s: %d", names(tab), as.integer(tab)), collapse = "; ")
        }
    }

    list(
        fact("contrasts tested",
             if (!is.null(inputs$contrasts)) paste(inputs$contrasts$Contrast_name, collapse = "; ") else NULL,
             "config contrasts file"),
        fact("libraries analysed", if (!is.null(sample_cols)) length(sample_cols) else NULL,
             "Datasets/rna_counts_filtered.tsv"),
        fact("libraries per group", group_desc, "sample metadata as loaded"),
        fact("genes after expression filtering", if (!is.null(counts)) nrow(counts) else NULL,
             "Datasets/rna_counts_filtered.tsv"),
        fact("genes receiving an adjusted p-value",
             if (!is.null(final)) sum(!is.na(.fact_col(final, "padj"))) else NULL,
             "Datasets/final_results.tsv"),
        fact("significance rule",
             sprintf("adjusted p <= %s and |linear fold change| >= %s", p_cut, fc_cut),
             "../execution_info/config_used.yaml"),
        fact("DESeq2 mode", config$modes$rna$de$deseq_mode %||% "default",
             "../execution_info/config_used.yaml")
    )
}

#' Differential expression rows
#' @keywords internal
.facts_rnaseq_de <- function(final, de_counts, de_cfg) {
    rows <- list()

    if (!is.null(de_counts) && nrow(de_counts) > 0) {
        rows <- c(rows, lapply(seq_len(nrow(de_counts)), function(i) {
            fact(sprintf("differentially expressed genes, %s", de_counts$contrast[i]),
                 sprintf("%s total (%s up, %s down)",
                         de_counts$total[i], de_counts$up[i], de_counts$down[i]),
                 "de_summary_counts.tsv")
        }))
    }

    if (is.null(final)) return(rows)

    mean_cols <- grep("^Mean\\.", names(final), value = TRUE)

    # One set of rows per contrast. final_results holds pvalue.<contrast>,
    # padj.<contrast> and linearFC.<contrast> side by side; reading only the
    # first of each made every claim below silently about contrast 1, and mixed
    # its fold changes with the any-contrast significance flag.
    per_contrast <- unlist(lapply(
        .fact_contrasts(final, c("pvalue", "padj", "linearFC")),
        function(cn) .facts_rnaseq_de_contrast(final, cn, mean_cols)
    ), recursive = FALSE, use.names = FALSE)

    c(rows, per_contrast, list(.fact_onoff(final, mean_cols)))
}

#' Differential expression rows for one contrast
#' @keywords internal
.facts_rnaseq_de_contrast <- function(final, cn, mean_cols) {
    p    <- .fact_col(final, "pvalue", cn)
    padj <- .fact_col(final, "padj", cn)
    lin  <- .fact_col(final, "linearFC", cn)

    # Denominator for the by-chance comparator: genes that HAVE a raw p-value.
    # DESeq2's independent filtering leaves padj NA for plenty of genes that
    # still carry one, and those genes are counted in the observed figure, so
    # using non-NA padj here understated what chance alone would produce.
    n_raw <- sum(!is.na(p))

    list(
        fact(sprintf("genes at raw p < 0.05, %s", cn),
             if (length(p)) sum(p < 0.05, na.rm = TRUE) else NULL,
             "Datasets/final_results.tsv"),
        fact(sprintf("genes expected at raw p < 0.05 by chance, %s", cn),
             if (n_raw > 0) round(0.05 * n_raw) else NULL,
             sprintf("derived: 0.05 x %s genes with a raw p-value", n_raw)),
        fact(sprintf("smallest adjusted p-value, %s", cn),
             if (any(!is.na(padj))) signif(min(padj, na.rm = TRUE), 3) else NULL,
             "Datasets/final_results.tsv"),
        fact(sprintf("|linear fold change| among differentially expressed genes detected in both groups, %s", cn),
             .fact_de_fc_range(final, lin, mean_cols, cn),
             "Datasets/final_results.tsv")
    )
}

#' Count differentially expressed genes with a zero group mean
#' @keywords internal
.fact_onoff <- function(final, mean_cols) {
    if (length(mean_cols) < 2 || !"pass_any_contrast" %in% names(final)) return(NULL)
    de <- final[!is.na(final$pass_any_contrast), mean_cols, drop = FALSE]
    if (nrow(de) == 0) return(NULL)
    zero <- rowSums(as.matrix(de) == 0, na.rm = TRUE) > 0
    fact("differentially expressed genes with a zero group mean",
         sprintf("%d of %d; the fold change is unbounded for these", sum(zero), nrow(de)),
         paste0("Datasets/final_results.tsv (", paste(mean_cols, collapse = ", "), ")"))
}

#' Fold-change range over genes detected in both groups, for one contrast
#'
#' Uses that contrast's own \code{<contrast>_pass} flag
#' (\code{R/domain/rnaseq/04_de_summary.R}). Pairing the fold changes of one
#' contrast with \code{pass_any_contrast} would fold in genes significant only
#' somewhere else.
#' @keywords internal
.fact_de_fc_range <- function(final, lin, mean_cols, cn = NULL) {
    if (length(lin) == 0) return(NULL)
    pass_col <- if (!is.null(cn) && paste0(cn, "_pass") %in% names(final)) {
        paste0(cn, "_pass")
    } else if ("pass_any_contrast" %in% names(final)) {
        "pass_any_contrast"
    } else {
        return(NULL)
    }
    keep <- !is.na(final[[pass_col]])
    if (length(mean_cols) >= 2) {
        keep <- keep & rowSums(as.matrix(final[, mean_cols, drop = FALSE]) == 0, na.rm = TRUE) == 0
    }
    fmt_range(abs(lin[keep]))
}

#' Sample-structure rows: PCA variance and between-library correlation
#' @keywords internal
.facts_rnaseq_structure <- function(norm, norm_path, pre, config) {
    if (is.null(norm) || ncol(norm) < 3) return(NULL)
    m <- as.matrix(norm[, -1, drop = FALSE])
    if (!is.numeric(m) || nrow(m) < 2) return(NULL)
    src <- if (is.na(norm_path)) "normalised matrix" else file.path("Datasets", basename(norm_path))

    # The preprocessed-input path keeps every numeric feature and can carry NA
    # or non-finite values into the normalised TSV, where QC's own PCA removes
    # them before plotting. Without the same step prcomp() aborts -- dropping
    # the PCA row with no explanation -- and cor() returns NA, which reaches the
    # sheet as a literal "NA and NA" claim.
    n_genes_all <- nrow(m)
    m <- m[is.finite(rowSums(m)), , drop = FALSE]
    if (nrow(m) < 2) return(NULL)
    src <- if (nrow(m) == n_genes_all) {
        src
    } else {
        sprintf("%s (recomputed on the %d of %d genes with no missing values)",
                src, nrow(m), n_genes_all)
    }

    pc <- tryCatch(stats::prcomp(t(m), center = TRUE, scale. = FALSE), error = function(e) NULL)
    pc_row <- NULL
    if (!is.null(pc)) {
        v <- pc$sdev^2
        v <- 100 * v / sum(v)
        pc_row <- fact("variance explained by PC1 and PC2",
                       sprintf("%.2f%% and %.2f%%", v[1], v[2]),
                       src)
    }

    cm <- stats::cor(m, method = "pearson")
    off <- cm[upper.tri(cm)]

    grp_row <- NULL
    fcol <- config$modes$rna$effects$color %||% NULL
    fcol <- if (is.list(fcol)) as.character(fcol[[1]]) else as.character(fcol)[1]
    sample_col <- config$modes$rna$effects$samples %||% "SampleID"
    if (!is.null(pre$meta) && !is.na(fcol) && fcol %in% names(pre$meta) &&
        sample_col %in% names(pre$meta)) {
        g <- as.character(pre$meta[[fcol]])[match(colnames(m), as.character(pre$meta[[sample_col]]))]
        if (!anyNA(g)) {
            same <- outer(g, g, "==")
            within <- cm[upper.tri(cm) & same]
            between <- cm[upper.tri(cm) & !same]
            if (length(within) && length(between)) {
                grp_row <- fact("mean correlation within group and between groups",
                                sprintf("%.3f and %.3f", mean(within), mean(between)),
                                src)
            }
        }
    }

    list(
        pc_row,
        fact("pearson correlation between libraries", fmt_range(off, 3),
             src),
        grp_row
    )
}

#' Gene-set enrichment rows
#'
#' Two layouts produce enrichment results and the sheet has to read both. The
#' online path writes flat \code{pathway_<contrast>_<db>_fgsea.csv} files at the
#' Enrichment root; the local path (\code{modes.rna.enrichment.annotation_dir})
#' writes \code{GSEA/<db>/ranking_by_<method>/<contrast>/results.csv} and the
#' matching nest under \code{ORA/}. Looking only for the flat names dropped the
#' whole enrichment section from local runs that had in fact produced results.
#' @keywords internal
.facts_rnaseq_enrichment <- function(enrich_dir, p_cut, contrast_names = character(0)) {
    if (is.null(enrich_dir) || !dir.exists(enrich_dir)) return(NULL)
    rows <- c(
        .facts_enrichment_flat(enrich_dir, p_cut, contrast_names),
        .facts_enrichment_nested(enrich_dir, p_cut)
    )
    if (length(rows) == 0) return(NULL)
    rows
}

#' Adjusted p-values from an enrichment table, whichever column holds them
#' @keywords internal
.fact_enrich_padj <- function(tab) {
    col <- intersect(c("padj", "p.adjust", "adj.P.Val"), names(tab))
    if (length(col) == 0) return(NULL)
    suppressWarnings(as.numeric(tab[[col[1]]]))
}

#' Enrichment rows from the flat online layout
#' @keywords internal
.facts_enrichment_flat <- function(enrich_dir, p_cut, contrast_names) {
    fgsea <- list.files(enrich_dir, pattern = "_fgsea\\.csv$", full.names = TRUE)
    if (length(fgsea) == 0) return(list())

    lapply(fgsea, function(f) {
        tab <- read_artifact_or_null(f, sep = ",")
        padj <- if (is.null(tab)) NULL else .fact_enrich_padj(tab)
        if (is.null(padj)) return(NULL)
        # Strip the known contrast name rather than guess where it ends: a
        # lazy quantifier is not portable here, and contrast names contain
        # underscores, so any positional split gets it wrong. The contrast is
        # kept and named in the claim: one file exists per contrast and per
        # database, so dropping it collapsed several different results into
        # repeated claims reading only "(KEGG)".
        stem <- sub("_fgsea\\.csv$", "", sub("^pathway_", "", basename(f)))
        collection <- stem
        contrast <- NA_character_
        for (cn in contrast_names) {
            if (startsWith(stem, paste0(cn, "_"))) {
                collection <- sub(paste0("^", cn, "_"), "", stem)
                contrast <- cn
                break
            }
        }
        label <- if (is.na(contrast)) collection else paste0(collection, ", ", contrast)

        ora_files <- list.files(
            enrich_dir,
            pattern = paste0("^", sub("_fgsea\\.csv$", "", basename(f)), "_ora_(up|down)\\.csv$"),
            full.names = TRUE
        )
        ora_sig <- if (length(ora_files)) {
            sum(vapply(ora_files, function(o) {
                t2 <- read_artifact_or_null(o, sep = ",")
                p2 <- if (is.null(t2)) NULL else .fact_enrich_padj(t2)
                if (is.null(p2)) return(0L)
                sum(!is.na(p2) & p2 <= p_cut)
            }, integer(1)))
        } else NA_integer_

        fact(sprintf("gene sets tested, ranked and over-representation hits (%s)", label),
             sprintf("%d tested; %d ranked at adjusted p <= %s; %s over-represented",
                     nrow(tab), sum(!is.na(padj) & padj <= p_cut), p_cut,
                     if (is.na(ora_sig)) "no files for" else ora_sig),
             paste0("Enrichment/", basename(f), " and the matching _ora_up/_ora_down files"))
    })
}

#' Enrichment rows from the nested local-annotation layout
#'
#' One row per result unit, labelled by the directories that identify it, so a
#' database/ranking/contrast combination is never merged with another.
#' @keywords internal
.facts_enrichment_nested <- function(enrich_dir, p_cut) {
    unlist(lapply(c("GSEA", "ORA"), function(analysis) {
        root <- file.path(enrich_dir, analysis)
        if (!dir.exists(root)) return(list())
        files <- list.files(root, pattern = "^results\\.csv$",
                            full.names = TRUE, recursive = TRUE)
        lapply(files, function(f) {
            tab <- read_artifact_or_null(f, sep = ",")
            padj <- if (is.null(tab)) NULL else .fact_enrich_padj(tab)
            if (is.null(padj)) return(NULL)
            # The unit directory path IS the identity of the result: for GSEA
            # <db>/ranking_by_<method>/<contrast>, for ORA <db>/<collection>/<round>.
            # Trimmed by length rather than by pattern: root is a real path and
            # would be read as a regex, where a "." matches any character.
            rel <- substring(f, nchar(root) + 2L)
            unit <- gsub("/", ", ", dirname(rel), fixed = TRUE)
            fact(sprintf("gene sets tested and hit at adjusted p <= %s (%s: %s)",
                         p_cut, analysis, unit),
                 sprintf("%d tested; %d significant", nrow(tab),
                         sum(!is.na(padj) & padj <= p_cut)),
                 file.path("Enrichment", analysis, rel))
        })
    }), recursive = FALSE, use.names = FALSE)
}

#' Run provenance rows read from execution_info
#' @keywords internal
.facts_run_provenance <- function(run_dir) {
    if (is.null(run_dir)) return(NULL)
    info <- file.path(run_dir, "execution_info")
    rd <- function(f) {
        p <- file.path(info, f)
        if (!file.exists(p)) return(NULL)
        trimws(paste(readLines(p, warn = FALSE), collapse = " "))
    }
    list(
        fact("run produced at", rd("timestamp.txt"), "../execution_info/timestamp.txt"),
        fact("pipeline commit", substr(rd("git_commit.txt") %||% "", 1, 12),
             "../execution_info/git_commit.txt"),
        fact("config used", rd("config_path.txt"), "../execution_info/config_path.txt")
    )
}

#' Pull a statistic column, for a named contrast or the first one present
#'
#' @param df The results table.
#' @param prefix Statistic prefix, e.g. "padj".
#' @param contrast Contrast name; NULL takes the first matching column.
#' @return Numeric vector, or numeric(0) when the column is absent.
#' @keywords internal
.fact_col <- function(df, prefix, contrast = NULL) {
    col <- if (!is.null(contrast)) {
        paste0(prefix, ".", contrast)
    } else {
        grep(paste0("^", prefix, "\\."), names(df), value = TRUE)[1]
    }
    if (is.na(col) || !col %in% names(df)) return(numeric(0))
    suppressWarnings(as.numeric(df[[col]]))
}

#' Contrast names carrying every one of the given statistic columns
#'
#' Derived from the table rather than the config, so the rows describe what was
#' actually written. Order follows the columns as written.
#'
#' @param df The results table.
#' @param prefixes Statistic prefixes every contrast must have.
#' @return Character vector of contrast names, possibly empty.
#' @keywords internal
.fact_contrasts <- function(df, prefixes) {
    per_prefix <- lapply(prefixes, function(p) {
        sub(paste0("^", p, "\\."), "", grep(paste0("^", p, "\\."), names(df), value = TRUE))
    })
    Reduce(intersect, per_prefix)
}
