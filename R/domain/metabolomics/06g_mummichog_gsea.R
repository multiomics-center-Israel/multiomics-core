# R/domain/metabolomics/06g_mummichog_gsea.R
#
# GSEA for untargeted MS peaks-to-pathways, run in EmpiricalCompound space —
# MetaboAnalyst's complementary analysis to mummichog's over-representation test
# (ORA). It is ADDITIVE: the pinned mummichog 2.7.0 ORA (06c/05b) and everything
# derived from it are left exactly as they are, including the ORA input contract
# (which sends logFC as mummichog's `statistic` column).
#
# The two analyses answer different questions:
#
#   * ORA  — takes only the features/ECs passing the p-value cutoff, asks whether
#            they over-populate a pathway, and reports an empirical permutation
#            p-value with an overlap / detected-pathway-size ratio.
#   * GSEA — takes the COMPLETE ranked EC list with no significance cutoff and
#            asks whether a pathway's ECs sit systematically at one end of that
#            ranking (ES / NES / p / padj), catching coordinated directional
#            shifts that never reach per-feature significance.
#
# ---------------------------------------------------------------------------
# Statistical semantics, verified against MetaboAnalystR (xia-lab/MetaboAnalystR)
# ---------------------------------------------------------------------------
# Read from `R/peaks_to_function.R` and `R/util_fgsea.R` at HEAD:
#
#  1. RANKING STATISTIC — `PreparePeakTable4PSEA()` ranks on the t-score: the
#     limma moderated `t` when the limma path is taken, otherwise the t-test /
#     ANOVA statistic. It is NOT the fold change. Our per-contrast limma tables
#     carry that moderated `t` in their `statistic` column (see `de_limma()`), so
#     that is the preferred metric here; `logFC` is used only when a DE table has
#     no usable statistic (e.g. a pre-computed table loaded from disk). The
#     choice is always recorded in the result and in the report provenance.
#
#  2. EC SPACE — with retention time and v2 (our configuration), MetaboAnalyst
#     performs GSEA in EmpiricalCompound space (`.compute.mummichog.RT.fgsea()`).
#     The EC score is built in two documented steps
#     (`peaks_to_function.R`, the `mumRT && version=="v2"` branch):
#       a. features sharing the same m/z are merged by `rank.metric`
#          (MetaboAnalyst's default is `"mean"`), then
#       b. `ec.exp.vec <- unlist(lapply(ec_exp_dict, max))` — an EC's score is
#          the **signed maximum** of its member features' scores. Not the mean,
#          not `max(abs())`.
#     `mmc_ec_scores()` implements exactly that.
#
#  3. GENE SETS — `current.mset <- mSetObj$pathways$emp_cpds`: pathway ->
#     detected EmpiricalCompounds. `mmc_pathway_ec_sets()` builds the same thing
#     by intersecting each EC's candidate compounds with the pathway's compounds.
#
#  4. ENGINE — MetaboAnalystR's `my.fgsea()` is a vendored copy of fgsea's
#     *simple* (fixed-permutation) scheme: `fgsea::calcGseaStat` for the
#     enrichment score, `calcGseaStatCumulativeBatch` permutations, the two-sided
#     `min((1+leEs)/(1+leZero), (1+geEs)/(1+geZero))` p-value, BH adjustment, and
#     `NES = ES / mean(same-sign permuted ES)`. We call `fgsea::fgseaSimple()`,
#     which IS that algorithm, rather than standing up a second GSEA engine or
#     using the multilevel estimator (`fgsea::fgsea()`, as our compound-space
#     `run_metabolomics_gsea()` does) whose p-values have different semantics.
#
#     Two deliberate departures from the vendored copy, both documented:
#       * We rank on the signed EC scores directly. MetaboAnalystR first collapses
#         tied scores into a single ranked position and passes tie-group indices;
#         with continuous moderated-t scores ties are rare, and collapsing them
#         would silently merge distinct ECs.
#       * We do not reproduce its `stats <- abs(stats); stats <- sort(stats,
#         decreasing = TRUE)` step, which (for signed input) re-orders the score
#         vector after the pathway positions were computed against the
#         signed order.
#     Both keep ES/NES sign meaningful, which the report relies on.
#
# NES sign is carried, never flattened to "up"/"down": positive NES = enrichment
# toward the positive end of the ranking statistic for the contrast as written
# (numerator vs denominator), negative NES = toward the negative end.
#
# Dependencies (R): fgsea (already used by run_metabolomics_gsea), ggplot2,
# withr, optionally ggrepel. All in renv.lock.


# ==== ranking statistic =====================================================

#' Choose and extract the signed GSEA ranking statistic from a DE table
#'
#' Prefers the moderated `t` statistic (our limma tables expose it as
#' `statistic`; a literal `t` column is also accepted) because that is what
#' MetaboAnalyst's peaks-to-pathways GSEA ranks on, and falls back to `logFC`
#' when no usable statistic is present. The ORA input contract is untouched —
#' this function is only ever used for GSEA.
#'
#' @param de_table One contrast's DE table (needs `feature_id` and at least one
#'   of `statistic` / `t` / `logFC`).
#' @return A list with `values` (named numeric vector, keyed by `feature_id`,
#'   finite values only), `metric` (`"moderated_t"` or `"logFC"`) and `label`
#'   (a human-readable provenance string). `NULL` when nothing is usable.
mmc_gsea_rank_metric <- function(de_table) {
  if (is.null(de_table) || !is.data.frame(de_table) || nrow(de_table) == 0 ||
      !"feature_id" %in% names(de_table)) {
    return(NULL)
  }
  stat_col <- .mmc_find_col(de_table, c("statistic", "t"))
  usable <- function(col) {
    if (is.null(col) || !col %in% names(de_table)) return(FALSE)
    v <- suppressWarnings(as.numeric(de_table[[col]]))
    any(is.finite(v))
  }

  if (usable(stat_col)) {
    v      <- suppressWarnings(as.numeric(de_table[[stat_col]]))
    metric <- "moderated_t"
    label  <- paste0("moderated t statistic (DE column '", stat_col, "')")
  } else if (usable("logFC")) {
    v      <- suppressWarnings(as.numeric(de_table$logFC))
    metric <- "logFC"
    label  <- "log2 fold change (no t statistic in the DE table)"
  } else {
    return(NULL)
  }

  keep <- is.finite(v)
  if (!any(keep)) return(NULL)
  list(
    values = stats::setNames(v[keep], as.character(de_table$feature_id)[keep]),
    metric = metric,
    label  = label
  )
}


# ==== EmpiricalCompound scores ==============================================

#' Build the ranked EmpiricalCompound score vector (MetaboAnalyst v2 semantics)
#'
#' Two steps, reproducing MetaboAnalystR's EC-space construction exactly:
#'
#' \enumerate{
#'   \item features that share the same m/z are merged with `rank_metric`
#'     (MetaboAnalyst's default is `"mean"`), then
#'   \item each EmpiricalCompound takes the **signed maximum** of its member
#'     features' merged scores.
#' }
#'
#' The signed max is MetaboAnalyst's rule (`unlist(lapply(ec_exp_dict, max))`),
#' not an aggregation we chose: no mean, no `max(abs())`.
#'
#' @param ec_features  Per-(EC, feature) frame from `read_mummichog_ec_features()`
#'   (needs `EID`, `feature_id`, `mz`).
#' @param feature_scores Named numeric vector of per-feature ranking statistics,
#'   keyed by `feature_id` (from `mmc_gsea_rank_metric()$values`).
#' @param rank_metric  How to merge features sharing one m/z: `"mean"` (default,
#'   MetaboAnalyst's default), `"min"`, `"max"` or `"median"`.
#' @return A named numeric vector of EC scores, sorted decreasing; empty named
#'   numeric when nothing can be scored.
mmc_ec_scores <- function(ec_features, feature_scores,
                          rank_metric = c("mean", "min", "max", "median")) {
  rank_metric <- match.arg(rank_metric)
  out <- stats::setNames(numeric(0), character(0))
  if (is.null(ec_features) || !is.data.frame(ec_features) ||
      nrow(ec_features) == 0 || length(feature_scores) == 0) {
    return(out)
  }
  score <- unname(feature_scores[as.character(ec_features$feature_id)])
  # Keyed on the m/z string mummichog echoed back, so the grouping matches its
  # own notion of "the same m/z" without float-equality surprises. A missing m/z
  # becomes its own group rather than being dropped by the formula interface.
  mz <- as.character(ec_features$mz)
  na_mz <- is.na(mz)
  if (any(na_mz)) mz[na_mz] <- paste0("<na>", which(na_mz))

  df <- data.frame(
    EID   = as.character(ec_features$EID),
    mz    = mz,
    score = score,
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$score) & !is.na(df$EID), , drop = FALSE]
  if (nrow(df) == 0) return(out)

  # -- step 1: merge scores of features sharing one m/z ----------------------
  agg <- switch(rank_metric,
                mean   = mean,
                min    = min,
                max    = max,
                median = stats::median)
  merged <- stats::aggregate(score ~ mz, data = df, FUN = agg)
  df$score <- merged$score[match(df$mz, merged$mz)]

  # -- step 2: EC score = signed max over its member features ---------------
  # Distinct (EC, m/z) pairs only: an m/z matching the same EC through several
  # candidate compounds must not weight the max differently.
  uniq <- unique(df[, c("EID", "mz", "score")])
  ec   <- stats::aggregate(score ~ EID, data = uniq, FUN = max)
  sort(stats::setNames(ec$score, ec$EID), decreasing = TRUE)
}

#' Build pathway gene-set equivalents: pathway -> detected EmpiricalCompounds
#'
#' An EmpiricalCompound belongs to a pathway when ANY of its candidate compounds
#' is a compound of that pathway — the same relation MetaboAnalyst holds as
#' `mSetObj$pathways$emp_cpds`. Only ECs present in `ec_ids` (the ranked
#' universe) are kept, so set sizes describe what was actually detected.
#'
#' @param model_pathways Named list: pathway -> compound ids (from
#'   `read_mummichog_model_pathways()$pathways`).
#' @param ec_candidates  Long EC -> candidate frame from
#'   `read_mummichog_ec_candidates()`.
#' @param ec_ids         Character vector of ECs in the ranked universe.
#' @return Named list: pathway -> character vector of detected EC ids. Pathways
#'   with no detected EC are dropped.
mmc_pathway_ec_sets <- function(model_pathways, ec_candidates, ec_ids) {
  if (length(model_pathways) == 0 || is.null(ec_candidates) ||
      nrow(ec_candidates) == 0 || length(ec_ids) == 0) {
    return(list())
  }
  cand <- ec_candidates[ec_candidates$EID %in% ec_ids, , drop = FALSE]
  if (nrow(cand) == 0) return(list())
  # compound id -> ECs that carry it as a candidate
  cpd_to_ec <- split(cand$EID, cand$compound_id)

  sets <- lapply(model_pathways, function(cpds) {
    hits <- cpd_to_ec[intersect(cpds, names(cpd_to_ec))]
    unique(unlist(hits, use.names = FALSE))
  })
  sets[vapply(sets, length, integer(1)) > 0L]
}


# ==== the GSEA run ==========================================================

#' Describe what a positive NES means for a contrast, in the contrast's own terms
#'
#' `"LL_vs_HL"` (our `make_contrast_label()` shape: numerator_vs_denominator)
#' becomes "positive NES = enrichment toward higher <metric> in LL relative to
#' HL". A label that does not parse gets the neutral statistic-based wording
#' instead of an invented group direction.
#'
#' @param contrast    Contrast label, or NULL.
#' @param metric_label Human-readable ranking-metric label.
#' @return A single character scalar.
mmc_gsea_direction_note <- function(contrast, metric_label) {
  parts <- if (is.null(contrast) || !nzchar(contrast)) character(0) else
    strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) == 2L && all(nzchar(parts))) {
    sprintf(paste0("positive NES = enrichment toward higher %s in %s relative ",
                   "to %s; negative NES = the opposite direction of the same ",
                   "contrast (%s vs %s)"),
            metric_label, parts[[1]], parts[[2]], parts[[1]], parts[[2]])
  } else {
    sprintf(paste0("positive NES = enrichment toward the positive end of the ",
                   "ranked %s; negative NES = toward the negative end"),
            metric_label)
  }
}

#' Run MetaboAnalyst-style GSEA on the mummichog EmpiricalCompound universe
#'
#' Uses the FULL ranked EC list — no feature significance cutoff — so this is a
#' genuinely different question from the ORA the pinned engine already reported.
#' The engine's own results are read-only inputs here; nothing is recomputed or
#' overwritten.
#'
#' @param files    One contrast's mummichog output files.
#' @param de_table That contrast's DE table (for the signed ranking statistic).
#' @param model    Model index from `read_mummichog_model_pathways()` /
#'   `mmc_load_model_pathways()`.
#' @param contrast Contrast label, used for the direction wording.
#' @param min_size,max_size EC-set size bounds passed to fgsea. `min_size = 2`
#'   matches this pipeline's compound-space `run_metabolomics_gsea()`;
#'   MetaboAnalystR leaves fgsea's `minSize = 1`, which admits single-EC
#'   "pathways" whose enrichment score is degenerate. Raise or lower it
#'   deliberately — it filters which pathways are testable, not how they score.
#' @param n_perm   Permutations for the fixed-permutation scheme.
#' @param seed     RNG seed — fgsea permutes, so the run is seeded to keep the
#'   report reproducible across machines and reruns.
#' @param rank_metric How features sharing one m/z are merged (see
#'   `mmc_ec_scores()`).
#' @return `NULL` when GSEA cannot run (fgsea unavailable, no model, no ECs, no
#'   usable statistic, no EC set of at least `min_size`). Otherwise a list:
#'   \describe{
#'     \item{table}{One row per pathway with `Pathway`,
#'       `Pathway size (model compounds)`, `Detected ECs`, `Hits used in GSEA`,
#'       `ES`, `NES`, `P.Value`, `padj`, `Direction`,
#'       `Leading-edge EmpiricalCompounds`, sorted by ascending `P.Value`.}
#'     \item{ranks}{The ranked EC score vector actually used.}
#'     \item{metric, metric_label}{Which ranking statistic was chosen, and why.}
#'     \item{direction_note}{NES orientation for this exact contrast.}
#'     \item{contrast, n_ec, n_pathways, n_perm, seed}{Provenance.}
#'   }
run_mummichog_gsea <- function(files, de_table, model, contrast = NULL,
                               min_size = 2, max_size = 500,
                               n_perm = 1000, seed = 42,
                               rank_metric = "mean") {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    message("mummichog GSEA: fgsea not available — skipping")
    return(NULL)
  }
  if (is.null(model) || length(model$pathways) == 0) {
    message("mummichog GSEA: no metabolic model pathways available — skipping")
    return(NULL)
  }
  rank <- mmc_gsea_rank_metric(de_table)
  if (is.null(rank)) {
    message("mummichog GSEA: no usable ranking statistic in the DE table — skipping")
    return(NULL)
  }

  ec_feat <- read_mummichog_ec_features(files)
  ec_cand <- read_mummichog_ec_candidates(files)
  if (nrow(ec_feat) == 0 || nrow(ec_cand) == 0) {
    message("mummichog GSEA: no EmpiricalCompound tables among the mummichog ",
            "outputs — skipping")
    return(NULL)
  }

  ranks <- mmc_ec_scores(ec_feat, rank$values, rank_metric = rank_metric)
  if (length(ranks) < min_size) {
    message("mummichog GSEA: only ", length(ranks),
            " scored EmpiricalCompounds — skipping")
    return(NULL)
  }
  sets <- mmc_pathway_ec_sets(model$pathways, ec_cand, names(ranks))
  sets <- sets[vapply(sets, length, integer(1)) >= min_size]
  if (length(sets) == 0) {
    message("mummichog GSEA: no pathway has >= ", min_size,
            " detected EmpiricalCompounds — skipping")
    return(NULL)
  }

  # fgseaSimple = the fixed-permutation scheme MetaboAnalystR vendored. Seeded,
  # because permutation p-values would otherwise drift between reruns.
  res <- tryCatch(
    withr::with_seed(seed, withCallingHandlers(
      fgsea::fgseaSimple(pathways = sets, stats = ranks, nperm = n_perm,
                         minSize = min_size, maxSize = max_size,
                         scoreType = "std"),
      # Tied EC scores are expected (features can share a statistic) and
      # MetaboAnalystR warns about them too; mute just that one so it does not
      # bury a real warning, and let everything else through.
      warning = function(w) {
        if (grepl("ties in the preranked stats", conditionMessage(w),
                  fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )),
    error = function(e) {
      warning("mummichog GSEA: fgseaSimple failed: ", conditionMessage(e),
              call. = FALSE)
      NULL
    }
  )
  if (is.null(res) || nrow(res) == 0) {
    message("mummichog GSEA: fgsea returned no results")
    return(NULL)
  }

  # Pathway size in MODEL compounds (what the pathway contains) vs detected ECs
  # (what this dataset actually saw) vs the set size fgsea scored — three
  # different numbers, all reported rather than conflated.
  model_size <- vapply(as.character(res$pathway),
                       function(p) as.numeric(length(model$pathways[[p]])),
                       numeric(1), USE.NAMES = FALSE)
  detected <- vapply(as.character(res$pathway),
                     function(p) as.numeric(length(sets[[p]])),
                     numeric(1), USE.NAMES = FALSE)

  lead <- vapply(seq_len(nrow(res)), function(i) {
    le <- res$leadingEdge[[i]]
    if (is.null(le) || length(le) == 0) NA_character_ else paste(le, collapse = ";")
  }, character(1))

  dir_note <- mmc_gsea_direction_note(contrast, rank$label)
  parts <- if (is.null(contrast) || !nzchar(contrast)) character(0) else
    strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
  direction <- if (length(parts) == 2L && all(nzchar(parts))) {
    ifelse(res$NES >= 0, paste("toward", parts[[1]]), paste("toward", parts[[2]]))
  } else {
    ifelse(res$NES >= 0, "toward positive statistic", "toward negative statistic")
  }

  tbl <- data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    "Pathway"                          = as.character(res$pathway),
    "Pathway size (model compounds)"   = model_size,
    "Detected ECs"                     = detected,
    "Hits used in GSEA"                = as.numeric(res$size),
    "ES"                               = round(as.numeric(res$ES), 4),
    "NES"                              = round(as.numeric(res$NES), 4),
    "P.Value"                          = as.numeric(res$pval),
    "padj"                             = as.numeric(res$padj),
    "Direction"                        = direction,
    "Leading-edge EmpiricalCompounds"  = lead
  )
  tbl <- tbl[order(tbl[["P.Value"]], na.last = TRUE), , drop = FALSE]

  message("mummichog GSEA: ", nrow(tbl), " pathways over ", length(ranks),
          " EmpiricalCompounds, ranked by ", rank$label)

  list(
    table          = tbl,
    ranks          = ranks,
    metric         = rank$metric,
    metric_label   = rank$label,
    direction_note = dir_note,
    contrast       = contrast,
    n_ec           = length(ranks),
    n_pathways     = nrow(tbl),
    n_perm         = n_perm,
    seed           = seed
  )
}


# ==== plot ==================================================================

#' MetaboAnalyst-style GSEA summary scatter for peaks-to-pathways
#'
#' Reproduces the encoding of MetaboAnalystR's `PlotPeaks2Paths()` GSEA branch
#' (`peaks_to_function.R`):
#'
#' \itemize{
#'   \item x = `NES`
#'   \item y = `-log10(GSEA p-value)`
#'   \item colour = `NES` on a diverging scale centred at 0, in MetaboAnalyst's
#'     own colours (`#458B00` -> `#fffee0` -> `#7f0000`)
#'   \item point size = `sqrt(-log10(p))`, which is literally what MetaboAnalyst
#'     computes there (`radi.vec <- sqrt(abs(y))`, with `y` the -log10 p already
#'     on the vertical axis). It is therefore a redundant re-encoding of
#'     significance and carries no additional biological quantity — the legend
#'     says exactly that instead of inventing a meaning for it.
#'   \item the most significant pathways labelled with `ggrepel`
#' }
#'
#' Two additions to the reference, both non-encoding: the points are drawn as
#' outlined circles (`shape = 21`) so the near-white middle of the diverging
#' scale stays visible — matching the outlined bubbles used elsewhere in this
#' report — and a dashed p-cutoff line plus a NES = 0 guide are drawn, as in the
#' other pathway plots in this pipeline.
#'
#' @param gsea_df  The `table` from `run_mummichog_gsea()` (needs `Pathway`,
#'   `NES`, `P.Value`).
#' @param title    Plot title.
#' @param subtitle Plot subtitle (e.g. the contrast). Optional.
#' @param p_cutoff Significance cutoff for the dashed reference line.
#' @param label_top Number of most-significant pathways to label (needs
#'   `ggrepel`; silently unlabelled when it is unavailable).
#' @return A ggplot object, or `NULL` when there is nothing to plot (so the
#'   report section hides gracefully).
plot_mummichog_gsea_scatter <- function(gsea_df, title, subtitle = NULL,
                                        p_cutoff = 0.05, label_top = 8) {
  if (is.null(gsea_df) || !is.data.frame(gsea_df) || nrow(gsea_df) == 0) {
    return(NULL)
  }
  if (!all(c("Pathway", "NES", "P.Value") %in% names(gsea_df))) return(NULL)

  df <- data.frame(
    pathway = as.character(gsea_df$Pathway),
    NES     = suppressWarnings(as.numeric(gsea_df$NES)),
    p       = suppressWarnings(as.numeric(gsea_df[["P.Value"]])),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$NES) & is.finite(df$p), , drop = FALSE]
  if (nrow(df) == 0) return(NULL)

  # -log10 with a floor so an exact-zero p can't blow up the axis.
  df$neg_log10_p <- -log10(pmax(df$p, .Machine$double.eps))
  # MetaboAnalyst's radi.vec, verbatim: sqrt of the -log10 p on the y axis.
  df$significance_size <- sqrt(df$neg_log10_p)

  lim <- max(abs(df$NES), na.rm = TRUE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$NES, y = .data$neg_log10_p)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey80") +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff),
                        linetype = "dashed", colour = "grey50") +
    ggplot2::annotate("text", x = Inf, y = -log10(p_cutoff),
                      label = sprintf("p = %g", p_cutoff),
                      hjust = 1.05, vjust = -0.5, size = 3, colour = "grey40") +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$significance_size, fill = .data$NES),
      shape = 21, colour = "grey25", stroke = 0.4, alpha = 0.9) +
    ggplot2::scale_fill_gradient2(
      low = "#458B00", mid = "#fffee0", high = "#7f0000", midpoint = 0,
      limits = c(-lim, lim), name = "NES") +
    # ASCII legend text on purpose: this label is rendered by the graphics
    # device on PNG/PDF export, where a non-ASCII glyph warns and degrades on a
    # minimal font setup.
    ggplot2::scale_size_continuous(
      range = c(1, 5), name = "Significance\nsqrt(-log10 p)") +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      x        = "NES  (normalized enrichment score)",
      y        = "-log10(GSEA p-value)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5))

  # Label the most significant pathways. ggrepel is optional across the repo, so
  # guard it and simply skip labels when it is not installed.
  if (label_top > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    lab_df <- df[order(df$p), , drop = FALSE]
    lab_df <- lab_df[seq_len(min(label_top, nrow(lab_df))), , drop = FALSE]
    p <- p + ggrepel::geom_text_repel(
      data = lab_df,
      ggplot2::aes(label = .data$pathway),
      size = 3, max.overlaps = Inf, box.padding = 0.5,
      min.segment.length = 0, seed = 42)
  }
  p
}
