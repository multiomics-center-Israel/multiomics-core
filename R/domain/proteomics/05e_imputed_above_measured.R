# Imputed-above-measured rule for proteomics (one resolver).
#
# flag_imputed_above_measured_proteomics() decides, per contrast, whether a
# row's imputed values sit above its measured ones. Two consumers read it:
# the exported imputed_above_measured.<contrast> column, on the exported draw
# (build_final_results_proteomics()), and, when
# de.exclude_imputed_above_measured is on, the DE step, which leaves flagged
# rows untested in that contrast on every draw
# (exclude_imputed_above_measured_proteomics()).

#' Label written in imputed_above_measured.<contrast> for flagged rows
IMPUTED_ABOVE_MEASURED_LABEL <- "FLAG: Imputation value higher than measured"

#' Flag features whose imputed values sit above their measured ones
#'
#' The imputation methods offered for proteomics all treat a missing value as
#' low-abundance and fill it from the low end of each sample's distribution. For
#' a protein that is itself near the bottom of that distribution, the filled-in
#' values can land above the ones that were actually measured -- and then a
#' protein never detected in one group comes out higher there, with a fold
#' change and p-value to match. This flags those rows so a reader can filter
#' them. On its own it changes no statistic; with
#' \code{de.exclude_imputed_above_measured: true} the DE step calls this same
#' function on each imputation run to leave the flagged rows untested
#' (\code{\link{exclude_imputed_above_measured_proteomics}}), so the exported
#' column and the exclusion cannot disagree about a row of the exported draw.
#'
#' Per contrast, the samples of the numerator and denominator groups are pooled
#' and the mean of the imputed cells is compared with the mean of the measured
#' cells. Both are read from \code{expr_model}, the matrix limma was fitted on
#' and the one exported as the \code{.norm} block, so the flag describes the
#' numbers the reader sees. Under \code{multi_imputation: true} that is the
#' first run; other runs drew different values and may flag differently.
#' The DE step passes each run's own matrix here instead.
#'
#' @param pre Proteomics preprocessing results (uses \code{expr_filt}, whose
#'   NAs mark the imputed cells, and \code{meta}).
#' @param expr_model Numeric matrix (features x samples) limma was fitted on.
#' @param contrasts_df Contrasts table (Contrast_name, Factor, Numerator,
#'   Denominator).
#' @param config Full pipeline config (sample-ID column).
#' @return Feature-indexed data.frame with one logical column per
#'   \code{Contrast_name}: TRUE when the imputed cells average above the
#'   measured ones, FALSE when not, NA when the contrast's samples hold no
#'   imputed or no measured value. NULL when an input is missing or the
#'   matrices do not share their features and samples.
flag_imputed_above_measured_proteomics <- function(pre, expr_model, contrasts_df,
                                                   config = NULL) {
    if (is.null(pre$expr_filt) || is.null(pre$meta) || is.null(expr_model)) return(NULL)

    prot_cfg <- config$modes$proteomics %||% list()
    sample_id_col <- prot_cfg$effects$samples %||%
        prot_cfg$id_columns$sample_col %||% "SampleID"

    measured <- as.matrix(pre$expr_filt)
    model <- as.matrix(expr_model)
    if (!all(rownames(measured) %in% rownames(model)) ||
        !all(colnames(measured) %in% colnames(model))) {
        warning("flag_imputed_above_measured_proteomics(): the model matrix does not ",
                "cover every feature and sample of expr_filt, so the ",
                "imputed_above_measured columns were left out.", call. = FALSE)
        return(NULL)
    }
    model <- model[rownames(measured), colnames(measured), drop = FALSE]

    imputed <- is.na(measured)
    imp_vals <- model
    imp_vals[!imputed] <- NA_real_
    obs_vals <- model
    obs_vals[imputed] <- NA_real_

    # Per-group sums and counts through the shared group resolver, so a
    # contrast's samples are exactly the ones behind Mean. and N.observed.
    per_group <- function(m, fn, prefix) {
        compute_group_stat_columns(expr = m, sample_meta = pre$meta,
                                   sample_id_col = sample_id_col,
                                   contrasts_df = contrasts_df,
                                   stat_fn = fn, prefix = prefix)
    }
    row_sum <- function(x) rowSums(x, na.rm = TRUE)
    row_n <- function(x) rowSums(!is.na(x))
    blocks <- list(s_imp = per_group(imp_vals, row_sum, "s_imp."),
                   n_imp = per_group(imp_vals, row_n, "n_imp."),
                   s_obs = per_group(obs_vals, row_sum, "s_obs."),
                   n_obs = per_group(obs_vals, row_n, "n_obs."))
    if (any(vapply(blocks, is.null, logical(1)))) return(NULL)

    flags <- lapply(seq_len(nrow(contrasts_df)), function(i) {
        grp <- c(as.character(contrasts_df$Numerator[i]),
                 as.character(contrasts_df$Denominator[i]))
        pooled <- function(b) {
            rowSums(as.matrix(blocks[[b]][, paste0(b, ".", grp), drop = FALSE]))
        }
        n_imp <- pooled("n_imp")
        n_obs <- pooled("n_obs")
        flag <- pooled("s_imp") / n_imp > pooled("s_obs") / n_obs
        flag[is.na(n_imp) | is.na(n_obs) | n_imp == 0 | n_obs == 0] <- NA
        unname(flag)
    })
    names(flags) <- as.character(contrasts_df$Contrast_name)

    out <- as.data.frame(flags, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(out) <- rownames(measured)
    out
}

#' Column marking, in a per-run DE table, the rows the exclusion left untested
IMPUTED_ABOVE_MEASURED_EXCLUDED_COL <- "imputed_above_measured_excluded"

#' Leave rows flagged by the imputed-above-measured rule untested, per contrast
#'
#' Called by \code{mod_proteomics_de()} between the per-run model fits and
#' \code{\link{summarize_limma_mult_imputation}} when
#' \code{de.exclude_imputed_above_measured} is on. The fits are the ones the
#' DE method made without the option (for pooled limma, every protein and every
#' group, so the variance moderation is unchanged); only the testing is
#' narrowed. For each imputation run and each contrast, the
#' rows that \code{\link{flag_imputed_above_measured_proteomics}} flags on that
#' run's own imputed matrix get \code{P.Value} and \code{adj.P.Val} set to NA,
#' and the adjusted p-value is recomputed with \code{de.p_adjust_method} over
#' the rows of that contrast that still have a p-value. A row excluded in one
#' contrast keeps its p-value in every other contrast. \code{logFC} and the
#' other columns are left as fitted.
#'
#' Because this runs before the summary, a run in which a row was excluded
#' counts as a run in which it did not pass. Each table also gets a logical
#' column \code{imputed_above_measured_excluded}, from which
#' \code{\link{summarize_limma_mult_imputation}} counts an excluded run as
#' \code{P.Value} = \code{adj.P.Val} = 1 in the summarised quantiles, unless
#' the row was excluded in every run: then it has no summarised p-value. With
#' one imputation run this is the plain rule: flagged rows have no p-value and
#' cannot pass.
#'
#' @param runs_de_tables List over imputation runs of named lists of
#'   per-contrast DE tables (each with the feature ID column, \code{P.Value}
#'   and \code{adj.P.Val}), as \code{mod_proteomics_de()} builds them.
#' @param imputations List of imputed matrices (features x samples), one per
#'   run, in the same order as \code{runs_de_tables}.
#' @param pre Proteomics preprocessing results (\code{expr_filt}, whose NAs
#'   mark the imputed cells, and \code{meta}).
#' @param contrasts_df Contrasts table the runs were fitted with.
#' @param config Full pipeline config (feature ID column, p-value adjustment
#'   method, sample-ID column).
#' @return A list with \code{runs_de_tables} (the input with the flagged rows'
#'   p-values removed, the adjusted p-values recomputed and the
#'   \code{imputed_above_measured_excluded} column added) and \code{excluded}:
#'   a list holding \code{n_runs} and \code{n_excluded}, an integer matrix
#'   (features x contrasts, named by feature ID and contrast) counting the runs
#'   in which each row was left untested.
exclude_imputed_above_measured_proteomics <- function(runs_de_tables, imputations,
                                                      pre, contrasts_df, config) {
    if (length(runs_de_tables) != length(imputations)) {
        stop("exclude_imputed_above_measured_proteomics(): got ",
             length(runs_de_tables), " runs of DE tables but ", length(imputations),
             " imputed matrices; each run must come with the matrix it was fitted on.",
             call. = FALSE)
    }
    p_cfg <- config$modes$proteomics %||% list()
    id_col <- (p_cfg$de_table %||% list())$id_col %||% "FeatureID"
    adj_method <- p_cfg$de$p_adjust_method %||% "BH"

    contrasts <- names(runs_de_tables[[1]])
    ref_ids <- as.character(runs_de_tables[[1]][[contrasts[1]]][[id_col]])
    n_excluded <- matrix(0L, nrow = length(ref_ids), ncol = length(contrasts),
                         dimnames = list(ref_ids, contrasts))

    for (i in seq_along(runs_de_tables)) {
        flags <- flag_imputed_above_measured_proteomics(pre, imputations[[i]],
                                                        contrasts_df, config)
        # The option asks for rows to be left untested; carrying on without the
        # flag would report every row as tested while claiming the rule ran.
        if (is.null(flags)) {
            stop("de.exclude_imputed_above_measured is on, but the imputed-above-",
                 "measured flag could not be computed for imputation run ", i,
                 " (see the warning above, or check that the contrasts' Factor ",
                 "column and the sample-ID column are in the metadata).\n",
                 "  Fix the inputs, or set de.exclude_imputed_above_measured: false.",
                 call. = FALSE)
        }
        for (cn in names(runs_de_tables[[i]])) {
            if (!cn %in% colnames(flags)) {
                stop("exclude_imputed_above_measured_proteomics(): contrast '", cn,
                     "' has DE results but no row in contrasts_df, so its ",
                     "imputed-above-measured flag is unknown.", call. = FALSE)
            }
            tab <- runs_de_tables[[i]][[cn]]
            ids <- as.character(tab[[id_col]])
            at <- match(ids, rownames(flags))
            if (anyNA(at)) {
                stop("exclude_imputed_above_measured_proteomics(): ", sum(is.na(at)),
                     " feature(s) in the DE table of contrast '", cn, "' (run ", i,
                     ") are not rows of expr_filt, e.g. '", ids[is.na(at)][1], "'.",
                     call. = FALSE)
            }
            hit <- flags[[cn]][at] %in% TRUE
            if (any(hit)) {
                tab$P.Value[hit] <- NA_real_
                # p.adjust() leaves NAs in place and sizes the correction on the
                # non-missing values, so this is the adjustment over the rows of
                # this contrast that are still tested.
                tab$adj.P.Val <- stats::p.adjust(tab$P.Value, method = adj_method)
            }
            # Read by summarize_limma_mult_imputation(), which has to tell an
            # excluded run apart from a row the method itself left untested.
            tab[[IMPUTED_ABOVE_MEASURED_EXCLUDED_COL]] <- hit
            runs_de_tables[[i]][[cn]] <- tab
            row_at <- match(ids, ref_ids)
            n_excluded[row_at[hit], cn] <- n_excluded[row_at[hit], cn] + 1L
        }
    }

    list(runs_de_tables = runs_de_tables,
         excluded = list(n_runs = length(runs_de_tables), n_excluded = n_excluded))
}

#' Rows the exclusion left untested in a contrast
#'
#' The rows enrichment drops from a contrast when
#' \code{de.exclude_imputed_above_measured} is on. The legacy enrichment
#' (\code{in_every_run = TRUE}) drops the rows excluded in every run, which are
#' exactly the rows left with no summarised p-value; a row excluded in only some
#' runs keeps a summarised p-value (its excluded runs counted as p = 1) and
#' stays. The pooled-ranking enrichment (\code{in_every_run = FALSE}) drops a
#' row excluded in any run, because Rubin pooling averages every draw's fit and
#' has no way to count a draw as untested.
#'
#' @param excluded The \code{imputed_above_measured_excluded} element of the
#'   proteomics DE result (\code{excluded} from
#'   \code{\link{exclude_imputed_above_measured_proteomics}}), or NULL.
#' @param contrast Contrast name as it appears in the summary columns.
#' @param in_every_run TRUE for rows excluded in every run, FALSE for rows
#'   excluded in at least one.
#' @return Character vector of feature IDs; empty when nothing was excluded or
#'   the exclusion did not run.
imputed_above_measured_excluded_ids <- function(excluded, contrast, in_every_run = TRUE) {
    if (is.null(excluded) || is.null(excluded$n_excluded)) return(character(0))
    m <- excluded$n_excluded
    col <- match(normalize_contrast_name(contrast), normalize_contrast_name(colnames(m)))
    if (is.na(col)) return(character(0))
    need <- if (isTRUE(in_every_run)) excluded$n_runs else 1L
    rownames(m)[m[, col] >= need]
}
