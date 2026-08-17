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
    p_cut <- de_cfg$p_cutoff %||% 0.05
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
#' @param out_dir The RNA-seq results directory.
#' @param dirs Result of \code{create_legacy_output_dirs()}.
#' @return A data.frame with columns contrast, up, down, total, or NULL.
#' @keywords internal
.read_de_summary_counts <- function(out_dir, dirs) {
    for (p in c(file.path(out_dir, "de_summary_counts.tsv"),
                file.path(dirs$datasets, "de_summary_counts.tsv"))) {
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
             "execution_info/config_used.yaml"),
        fact("DESeq2 mode", config$modes$rna$de$deseq_mode %||% "default",
             "execution_info/config_used.yaml")
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

    p <- .fact_col(final, "pvalue")
    padj <- .fact_col(final, "padj")
    lin <- .fact_col(final, "linearFC")
    mean_cols <- grep("^Mean\\.", names(final), value = TRUE)

    tested <- sum(!is.na(padj))
    c(rows, list(
        fact("genes at raw p < 0.05", if (length(p)) sum(p < 0.05, na.rm = TRUE) else NULL,
             "Datasets/final_results.tsv"),
        fact("genes expected at raw p < 0.05 by chance",
             if (tested > 0) round(0.05 * tested) else NULL,
             sprintf("derived: 0.05 x %s tested genes", tested)),
        fact("smallest adjusted p-value",
             if (any(!is.na(padj))) signif(min(padj, na.rm = TRUE), 3) else NULL,
             "Datasets/final_results.tsv"),
        .fact_onoff(final, mean_cols),
        fact("|linear fold change| among differentially expressed genes detected in both groups",
             .fact_de_fc_range(final, lin, mean_cols),
             "Datasets/final_results.tsv")
    ))
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

#' Fold-change range over genes detected in both groups
#' @keywords internal
.fact_de_fc_range <- function(final, lin, mean_cols) {
    if (length(lin) == 0 || !"pass_any_contrast" %in% names(final)) return(NULL)
    keep <- !is.na(final$pass_any_contrast)
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

    pc <- tryCatch(stats::prcomp(t(m), center = TRUE, scale. = FALSE), error = function(e) NULL)
    pc_row <- NULL
    if (!is.null(pc)) {
        v <- pc$sdev^2
        v <- 100 * v / sum(v)
        pc_row <- fact("variance explained by PC1 and PC2",
                       sprintf("%.2f%% and %.2f%%", v[1], v[2]),
                       paste0(src, " (recomputed)"))
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
                                paste0(src, " (recomputed)"))
            }
        }
    }

    list(
        pc_row,
        fact("pearson correlation between libraries", fmt_range(off, 3),
             paste0(src, " (recomputed)")),
        grp_row
    )
}

#' Gene-set enrichment rows, one per collection
#' @keywords internal
.facts_rnaseq_enrichment <- function(enrich_dir, p_cut, contrast_names = character(0)) {
    if (is.null(enrich_dir) || !dir.exists(enrich_dir)) return(NULL)
    fgsea <- list.files(enrich_dir, pattern = "_fgsea\\.csv$", full.names = TRUE)
    if (length(fgsea) == 0) return(NULL)

    lapply(fgsea, function(f) {
        tab <- read_artifact_or_null(f, sep = ",")
        if (is.null(tab) || !"padj" %in% names(tab)) return(NULL)
        # Strip the known contrast name rather than guess where it ends: a
        # lazy quantifier is not portable here, and contrast names contain
        # underscores, so any positional split gets it wrong.
        collection <- sub("_fgsea\\.csv$", "", sub("^pathway_", "", basename(f)))
        for (cn in contrast_names) {
            collection <- sub(paste0("^", cn, "_"), "", collection)
        }

        ora_files <- list.files(
            enrich_dir,
            pattern = paste0("^", sub("_fgsea\\.csv$", "", basename(f)), "_ora_(up|down)\\.csv$"),
            full.names = TRUE
        )
        ora_sig <- if (length(ora_files)) {
            sum(vapply(ora_files, function(o) {
                t2 <- read_artifact_or_null(o, sep = ",")
                if (is.null(t2) || !"padj" %in% names(t2)) return(0L)
                sum(!is.na(t2$padj) & t2$padj <= p_cut)
            }, integer(1)))
        } else NA_integer_

        fact(sprintf("gene sets tested, ranked and over-representation hits (%s)", collection),
             sprintf("%d tested; %d ranked at adjusted p <= %s; %s over-represented",
                     nrow(tab), sum(!is.na(tab$padj) & tab$padj <= p_cut), p_cut,
                     if (is.na(ora_sig)) "no files for" else ora_sig),
             paste0("Enrichment/", basename(f), " and the matching _ora_up/_ora_down files"))
    })
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
        fact("run produced at", rd("timestamp.txt"), "execution_info/timestamp.txt"),
        fact("pipeline commit", substr(rd("git_commit.txt") %||% "", 1, 12),
             "execution_info/git_commit.txt"),
        fact("config used", rd("config_path.txt"), "execution_info/config_path.txt")
    )
}

#' Pull the first column whose name starts with a statistic prefix
#' @keywords internal
.fact_col <- function(df, prefix) {
    hit <- grep(paste0("^", prefix, "\\."), names(df), value = TRUE)
    if (length(hit) == 0) return(numeric(0))
    suppressWarnings(as.numeric(df[[hit[1]]]))
}
