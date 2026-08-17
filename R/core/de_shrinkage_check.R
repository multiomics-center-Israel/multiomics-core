#' Detect fold-change collapse by comparing the model estimate to the group means
#'
#' The reported \code{log2FC} is a model coefficient; \code{log2FC_from_means}
#' is the same quantity read straight off the exported group means, with no
#' model behind it. When a fit misbehaves — over-aggressive shrinkage, a
#' mis-specified design — the first collapses towards zero while the second
#' keeps carrying the effect the data holds. In a volcano that shows up as a
#' vertical stripe of points at x = 0, some of them highly significant, which is
#' otherwise easy to read as "nothing changed" rather than "the estimates were
#' flattened".
#'
#' This check turns that picture into numbers the pipeline can act on.
#'
#' Two symptoms are measured per contrast:
#' \describe{
#'   \item{shrinkage}{The median of \code{|log2FC| / |log2FC_from_means|} over
#'     features whose model-free estimate is large enough to divide by. Near 1
#'     means the model estimate tracks the data; well below 1 means it does not.
#'     The floor matters: for a null feature both numbers are ~0 and their ratio
#'     is meaningless, so those features are excluded rather than allowed to
#'     dominate the median.}
#'   \item{the x = 0 stripe}{The fraction of significant features whose
#'     \code{log2FC} is flat while their model-free estimate is not. These are
#'     exactly the points that draw the vertical line in the volcano.}
#' }
#'
#' @param de_stats Data frame of DE statistics (the final results table or the
#'   per-omics summary): must carry the \code{log2FC}, \code{log2FC_from_means}
#'   and adjusted p-value columns for each contrast.
#' @param contrasts Character vector of contrast names to check.
#' @param mode Omics mode, used to resolve column names via
#'   \code{\link{get_contrast_cols}} (e.g. "rna", "proteomics").
#' @param p_cutoff Adjusted p-value cutoff defining "significant".
#' @param naive_floor Minimum \code{|log2FC_from_means|} for a feature to enter
#'   the ratio. Features below it are null-ish, where the ratio is unstable.
#' @param flat_tol \code{|log2FC|} at or below this counts as flattened to zero.
#' @param ratio_warn Median ratio at or below this raises the shrinkage alert.
#'   The default sits between ordinary DESeq2 \code{betaPrior} behaviour
#'   (~0.85 in our fixtures) and a genuine collapse (~0).
#' @param stripe_warn Fraction of significant-but-flat features at or above
#'   which the x = 0 stripe alert is raised.
#' @return A data.frame with one row per checked contrast and columns
#'   \code{contrast}, \code{n_considered}, \code{median_ratio},
#'   \code{frac_flat}, \code{n_significant}, \code{frac_sig_flat},
#'   \code{flag} ("ok"/"shrunk"/"collapsed"). Contrasts whose columns are
#'   missing are skipped. Returns NULL when nothing could be checked.
check_log2fc_shrinkage <- function(de_stats,
                                   contrasts,
                                   mode = "rna",
                                   p_cutoff = 0.05,
                                   naive_floor = 0.5,
                                   flat_tol = 0.01,
                                   ratio_warn = 0.5,
                                   stripe_warn = 0.1) {
    if (is.null(de_stats) || !is.data.frame(de_stats) || nrow(de_stats) == 0) return(NULL)
    if (is.null(contrasts) || length(contrasts) == 0) return(NULL)

    de_cols <- names(de_stats)

    rows <- lapply(contrasts, function(cn) {
        cols <- get_contrast_cols(cn, mode = mode)
        if (is.null(cols$log2fc) || is.null(cols$log2fc_means)) return(NULL)
        if (!all(c(cols$log2fc, cols$log2fc_means) %in% de_cols)) return(NULL)

        model <- as.numeric(de_stats[[cols$log2fc]])
        naive <- as.numeric(de_stats[[cols$log2fc_means]])
        padj  <- if (cols$padj %in% de_cols) as.numeric(de_stats[[cols$padj]]) else rep(NA_real_, length(model))

        usable <- is.finite(model) & is.finite(naive)
        considered <- usable & abs(naive) >= naive_floor

        # Ratio only where the denominator is meaningful (see naive_floor).
        median_ratio <- if (any(considered)) {
            stats::median(abs(model[considered]) / abs(naive[considered]))
        } else NA_real_

        frac_flat <- if (any(considered)) {
            mean(abs(model[considered]) <= flat_tol)
        } else NA_real_

        # The volcano stripe: significant, flat estimate, real effect in the data
        sig <- usable & !is.na(padj) & padj <= p_cutoff
        sig_flat <- sig & abs(model) <= flat_tol & abs(naive) >= naive_floor
        frac_sig_flat <- if (any(sig)) mean(sig_flat[sig]) else NA_real_

        flag <- "ok"
        if (!is.na(frac_flat) && frac_flat >= stripe_warn) {
            flag <- "collapsed"
        } else if (!is.na(frac_sig_flat) && frac_sig_flat >= stripe_warn) {
            flag <- "collapsed"
        } else if (!is.na(median_ratio) && median_ratio <= ratio_warn) {
            flag <- "shrunk"
        }

        data.frame(
            contrast      = cn,
            n_considered  = sum(considered),
            median_ratio  = median_ratio,
            frac_flat     = frac_flat,
            n_significant = sum(sig),
            frac_sig_flat = frac_sig_flat,
            flag          = flag,
            stringsAsFactors = FALSE
        )
    })

    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0) return(NULL)

    do.call(rbind, rows)
}

#' Raise the analyst-facing alert for a shrinkage check
#'
#' Emits one warning per flagged contrast, naming what to look at and where.
#' Split from \code{\link{check_log2fc_shrinkage}} so the measurement stays pure
#' and testable while the side effect lives on its own.
#'
#' @param check Result of \code{\link{check_log2fc_shrinkage}}.
#' @param mode Omics mode, quoted back in the message.
#' @return \code{check}, invisibly.
warn_log2fc_shrinkage <- function(check, mode = "rna") {
    if (is.null(check) || nrow(check) == 0) return(invisible(check))

    flagged <- check[check$flag != "ok", , drop = FALSE]
    if (nrow(flagged) == 0) {
        message(sprintf(
            "[%s] Fold-change check: no shrinkage detected (median |log2FC| / |log2FC_from_means| = %s across %d contrast(s)).",
            mode, paste(signif(check$median_ratio, 3), collapse = ", "), nrow(check)
        ))
        return(invisible(check))
    }

    for (i in seq_len(nrow(flagged))) {
        r <- flagged[i, ]
        headline <- if (identical(r$flag, "collapsed")) {
            sprintf("fold changes appear COLLAPSED towards zero (%.0f%% of the %d features with a real effect have |log2FC| <= 0.01)",
                    100 * r$frac_flat, r$n_considered)
        } else {
            sprintf("fold changes appear SHRUNK (median |log2FC| is only %.0f%% of |log2FC_from_means| across %d features)",
                    100 * r$median_ratio, r$n_considered)
        }

        stripe <- if (!is.na(r$frac_sig_flat) && r$frac_sig_flat > 0) {
            sprintf(" %.0f%% of the %d significant features sit at log2FC ~ 0, which draws a vertical stripe at x = 0 in the volcano.",
                    100 * r$frac_sig_flat, r$n_significant)
        } else ""

        warning(sprintf(
            "[%s] Contrast '%s': %s.%s Compare the log2FC and log2FC_from_means columns in Final_results, and check the DE settings (for RNA-seq, de$deseq_mode; any log2FC shrinkage applied upstream). Full numbers: Datasets/log2fc_shrinkage_check.tsv",
            mode, r$contrast, headline, stripe
        ), call. = FALSE)
    }

    invisible(check)
}

#' Run the shrinkage check, alert the analyst, and write the numbers to disk
#'
#' Convenience wrapper for the export layer: measure, warn, persist. The TSV is
#' written unconditionally so a clean run is on record too, not only a bad one.
#'
#' @param de_stats Final results / DE summary table.
#' @param contrasts_df Contrasts table with a \code{Contrast_name} column.
#' @param mode Omics mode ("rna", "proteomics", ...).
#' @param p_cutoff Adjusted p-value cutoff defining "significant".
#' @param out_dir Directory to write \code{log2fc_shrinkage_check.tsv} into.
#'   NULL skips writing.
#' @return Path to the written TSV, or character(0) when nothing was written.
run_log2fc_shrinkage_check <- function(de_stats, contrasts_df, mode = "rna",
                                       p_cutoff = 0.05, out_dir = NULL) {
    if (is.null(contrasts_df) || !"Contrast_name" %in% names(contrasts_df)) {
        return(character(0))
    }

    check <- check_log2fc_shrinkage(
        de_stats  = de_stats,
        contrasts = as.character(contrasts_df$Contrast_name),
        mode      = mode,
        p_cutoff  = p_cutoff
    )
    if (is.null(check)) return(character(0))

    warn_log2fc_shrinkage(check, mode = mode)

    if (is.null(out_dir)) return(character(0))
    save_tsv(check, out_dir, "log2fc_shrinkage_check.tsv")
}
