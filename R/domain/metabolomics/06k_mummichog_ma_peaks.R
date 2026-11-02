# R/domain/metabolomics/06k_mummichog_ma_peaks.R
#
# MetaboAnalyst-compatible mummichog, component 2: the peak list and m/z
# matching, including genuine mixed polarity.
#
# In MetaboAnalystR's mixed mode the positive- and negative-ion features form
# ONE peak list: one p-value ordering, one de-duplication, one retention-time
# ranking and one reference list. Each feature is then matched against the ion
# forms of ITS OWN polarity, and the two match tables are stacked into one. This
# file reproduces those steps; empirical compounds (component 3) and the pathway
# statistics (component 4) build on its output.
#
# ===========================================================================
# BEHAVIOURAL REFERENCE — pinned upstream version
# ===========================================================================
#     repo:   xia-lab/MetaboAnalystR
#     commit: 398476ae2a0c996e390925fa62adc46b5ecde334
#     file:   R/peaks_to_function.R
#
#   mmc_ma_prepare_peaks()   SanityCheckMummichogData()  L629-L983
#                            + the RT ranking at the top of .search.compoundLib()
#                              (L1629), see "Random tie-breaking" below
#   mmc_ma_input_mzlist()    SetMummichogPval()          L1042-L1098
#   mmc_ma_match_peaks()     .search.compoundLib()       L1605-L1806 (matching,
#                            up to and including .drop.currency.matches())
#
# Scope: the retention-time path with version "v2" — the path that builds
# empirical compounds and applies the primary-ion rule. A peak list without
# retention times takes a different, compound-level path upstream; it is not
# reproduced and is refused.
#
# Reproduced as upstream does it, including what may surprise:
#   * m/z kept only when 50 < m/z < 2001 (upstream's own bound, wider than
#     python mummichog's 2000);
#   * the list is ordered by p-value (stable for ties), and a feature whose m/z
#     EXACTLY equals an earlier one is dropped — across polarities, so of a
#     positive and a negative feature sharing an m/z only the smaller p
#     survives;
#   * m/z values become character strings (15 significant digits), the keys
#     every later upstream lookup uses;
#   * tolerance ppm * 1e-6 * m/z, strict `<`, over the floor-1, floor, floor+1
#     bins of the library tree;
#   * currency compounds (upstream's default list) removed from the match
#     table, unless that would remove every match.
#
# Where we are stricter than upstream (refusals, never different numbers):
#   * mixed mode needs each feature's polarity spelled exactly "positive" or
#     "negative" (mmc_feature_polarity() produces these); upstream silently
#     treats any other value, including a missing column, as positive or
#     negative;
#   * any feature lacking an m/z, retention time, p-value or statistic stops
#     the run; upstream drops such rows with a message. Our DE tables carry
#     every value, so a gap means something upstream of here went wrong.
#
# Random tie-breaking: upstream ranks retention times with
# rank(ties.method = "random"). The RNG state it draws from is the one
# SanityCheckMummichogData() leaves behind: set.seed(123), then any draws of its
# m/z de-duplication loop; nothing between the two calls draws. We perform both
# steps in one withr::with_seed(123, ...) block, which gives the same ranks in a
# session using R's default RNG kinds without touching the caller's RNG state.


# Upstream's default currency compounds (peaks_to_function.R L4405-L4409).
.MMC_MA_CURRENCY <- c("C00001", "C00080", "C00007", "C00006", "C00005",
                      "C00003", "C00004", "C00002", "C00013", "C00008",
                      "C00009", "C00011", "G11113", "", "H2O", "H+",
                      "Oxygen", "NADP+", "NADPH", "NAD+", "NADH", "ATP",
                      "Pyrophosphate", "ADP", "CO2")

#' Prepare a peak list the way MetaboAnalystR's SanityCheckMummichogData() does
#'
#' Orders the features by p-value, keeps 50 < m/z < 2001, drops exact duplicate
#' m/z (keeping the smaller p), turns m/z into the character keys upstream uses,
#' and ranks retention times over all features (see the file header for the
#' seeded tie-breaking). In `"mixed"` mode each feature keeps its own polarity;
#' in `"positive"` / `"negative"` mode every feature is treated as that
#' polarity, as upstream does.
#'
#' @param peaks Data frame with one row per feature and columns `feature_id`,
#'   `mz`, `rt`, `p_value`, `statistic`, and (mixed mode) `polarity`
#'   (`"positive"` / `"negative"`).
#' @param mode `"mixed"`, `"positive"` or `"negative"`.
#' @param rt_frac Fraction used for the retention-time tolerances (upstream
#'   default 0.02).
#' @param rt_tol Retention-time tolerance; `NA` (upstream default) derives it as
#'   `max(rt) * rt_frac`.
#' @param min_features Smallest peak list analysed (upstream: 100).
#' @return A list with `table` (the kept features in upstream order:
#'   `feature_id`, `p_value`, `mz` (character key), `statistic`, `rt`,
#'   `positive`), `ref_mzlist` (character), `pos_inx` (logical), `ret_time` and
#'   `rt_rank` (numeric, named by m/z key), `expr_dic` (statistic named by m/z
#'   key), `rt_tol`, `rt_frac`, `mode`, and `removed` (counts dropped by range
#'   and by duplicate m/z).
#' @examples
#' pk <- data.frame(feature_id = paste0("F", 1:3), mz = c(181.07, 179.06, 181.07),
#'                  rt = c(1, 1.1, 5), p_value = c(0.01, 0.2, 0.5),
#'                  statistic = c(3, -1, 0.2),
#'                  polarity = c("positive", "negative", "negative"))
#' prep <- mmc_ma_prepare_peaks(pk, "mixed", min_features = 1)
#' prep$table
mmc_ma_prepare_peaks <- function(peaks, mode = c("mixed", "positive", "negative"),
                                 rt_frac = 0.02, rt_tol = NA,
                                 min_features = 100) {
  mode <- match.arg(mode)
  need <- c("feature_id", "mz", "rt", "p_value", "statistic")
  if (mode == "mixed") need <- c(need, "polarity")
  miss <- setdiff(need, names(peaks))
  if (length(miss) > 0L) {
    .mmc_stop("the MetaboAnalyst-compatible peak list needs column(s) ",
              paste(miss, collapse = ", "), " (mode '", mode, "').",
              if ("rt" %in% miss) paste0(" Retention time is required: only ",
                "MetaboAnalyst's retention-time (empirical compound) path is ",
                "reproduced.") else "")
  }
  num <- function(x, nm) {
    v <- suppressWarnings(as.numeric(x))
    bad <- is.na(v) | !is.finite(v)
    if (any(bad)) {
      .mmc_stop(sum(bad), " feature(s) have no usable '", nm, "' (e.g. ",
                paste(utils::head(as.character(peaks$feature_id[bad]), 5),
                      collapse = ", "), "). Every feature needs an m/z, ",
                "retention time, p-value and statistic; none is dropped or ",
                "filled silently.")
    }
    v
  }
  mz   <- num(peaks$mz, "mz")
  rt   <- num(peaks$rt, "rt")
  pval <- num(peaks$p_value, "p_value")
  tsc  <- num(peaks$statistic, "statistic")

  positive <- switch(mode,
    positive = rep(TRUE, nrow(peaks)),
    negative = rep(FALSE, nrow(peaks)),
    mixed = {
      pol <- as.character(peaks$polarity)
      bad <- is.na(pol) | !pol %in% c("positive", "negative")
      if (any(bad)) {
        .mmc_stop(sum(bad), " feature(s) have a polarity other than ",
                  "\"positive\" or \"negative\" (e.g. ",
                  paste(utils::head(unique(pol[bad]), 5), collapse = " | "),
                  "). Mixed mode matches each feature against its own ",
                  "polarity's ion forms, so none is assumed.")
      }
      pol == "positive"
    })

  tbl <- data.frame(feature_id = as.character(peaks$feature_id),
                    p_value = pval, mz = mz, statistic = tsc, rt = rt,
                    positive = positive, stringsAsFactors = FALSE)

  # SanityCheckMummichogData(): order by p (order() is stable for ties) ...
  tbl <- tbl[order(tbl$p_value), , drop = FALSE]
  # ... trim to 50 < m/z < 2001 ...
  in_range <- tbl$mz > 50 & tbl$mz < 2001
  n_range  <- sum(!in_range)
  tbl <- tbl[in_range, , drop = FALSE]
  # ... and drop exact duplicate m/z, whatever their polarity.
  dup   <- duplicated(tbl$mz)
  n_dup <- sum(dup)
  tbl <- tbl[!dup, , drop = FALSE]
  rownames(tbl) <- NULL

  if (nrow(tbl) < min_features) {
    .mmc_stop("only ", nrow(tbl), " feature(s) remain after the m/z range and ",
              "duplicate filters; MetaboAnalyst refuses peak lists under ",
              min_features, " features as too few for enrichment.")
  }

  # The character keys, made unique the way upstream does (only distinct
  # doubles printing alike can collide), then the seeded RT ranking.
  keyed <- withr::with_seed(123, {
    mzs <- as.character(tbl$mz)
    mzs_unq <- mzs[duplicated(mzs)]
    while (length(mzs_unq) > 0) {
      mzs[duplicated(mzs)] <- vapply(mzs_unq, function(x) {
        paste0(x, sample(1:9, 1, replace = TRUE))
      }, FUN.VALUE = character(1L))
      mzs_unq <- mzs[duplicated(mzs)]
    }
    ret_time <- stats::setNames(tbl$rt, mzs)
    list(mzs = mzs, ret_time = ret_time,
         rt_rank = rank(ret_time, ties.method = "random"))
  })
  tbl$mz <- keyed$mzs

  if (is.na(rt_tol)) rt_tol <- max(keyed$ret_time) * rt_frac

  list(
    table      = tbl,
    ref_mzlist = keyed$mzs,
    pos_inx    = tbl$positive,
    ret_time   = keyed$ret_time,
    rt_rank    = keyed$rt_rank,
    expr_dic   = stats::setNames(tbl$statistic, keyed$mzs),
    rt_tol     = rt_tol,
    rt_frac    = rt_frac,
    mode       = mode,
    removed    = c(out_of_range = n_range, duplicate_mz = n_dup)
  )
}

#' Select the significant peaks the way MetaboAnalystR's SetMummichogPval() does
#'
#' One cutoff over the whole peak list — both polarities together — with
#' upstream's inclusive `p <= cutoff`. Upstream refuses fewer than 11
#' significant features, or fewer than 1% of the list; so does this.
#'
#' @param prep   Output of `mmc_ma_prepare_peaks()`.
#' @param cutoff p-value cutoff.
#' @return Character vector of the significant m/z keys, in peak-list order.
#' @examples
#' \dontrun{
#' sig <- mmc_ma_input_mzlist(prep, 0.05)
#' }
mmc_ma_input_mzlist <- function(prep, cutoff) {
  if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff)) {
    .mmc_stop("the MetaboAnalyst-compatible cutoff must be a single number; ",
              "got '", paste(format(cutoff), collapse = ", "), "'.")
  }
  sig <- prep$ref_mzlist[prep$table$p_value <= cutoff]
  n <- length(sig)
  part <- round(100 * n / length(prep$ref_mzlist), 2)
  if (n == 0L) {
    .mmc_stop("no feature has p <= ", cutoff, ". Adjust p_cutoff.")
  }
  if (n < 11 || part < 1) {
    .mmc_stop("only ", n, " feature(s) (", part, "%) have p <= ", cutoff,
              "; MetaboAnalyst refuses fewer than 11 significant features or ",
              "under 1% of the list as too few for enrichment. Adjust p_cutoff.")
  }
  sig
}

#' Match a prepared peak list against a MetaboAnalyst library by polarity
#'
#' Port of the matching half of MetaboAnalystR's `.search.compoundLib()` on the
#' retention-time path: positive-ion features against the library's positive ion
#' forms, negative-ion features against its negative forms, both tables stacked
#' (positive first), then currency compounds removed. Row order, column order
#' and the character encoding of every value are upstream's.
#'
#' @param prep Output of `mmc_ma_prepare_peaks()`.
#' @param lib  Library from `mmc_ma_library_from_model()`.
#' @param ppm  Instrument mass accuracy in ppm.
#' @param exclude_currency Drop currency-compound matches (upstream default TRUE).
#' @return Data frame of character columns `Query.Mass`, `Matched.Compound`,
#'   `Matched.Form`, `Mass.Diff`, `Retention.Time`, `RT.Rank`, one row per
#'   (feature, compound, ion form) match.
#' @examples
#' \dontrun{
#' matched <- mmc_ma_match_peaks(prep, lib, ppm = 10)
#' }
mmc_ma_match_peaks <- function(prep, lib, ppm, exclude_currency = TRUE) {
  if (!is.numeric(ppm) || length(ppm) != 1L || is.na(ppm) || ppm <= 0) {
    .mmc_stop("ppm must be a single positive number; got '",
              paste(format(ppm), collapse = ", "), "'.")
  }
  ref_mzlist <- as.numeric(prep$ref_mzlist)
  pos_inx <- prep$pos_inx

  one_polarity <- function(inx, mz_mat, tree) {
    mzs   <- ref_mzlist[inx]
    rts   <- prep$ret_time[inx]
    ranks <- prep$rt_rank[inx]
    tols  <- ppm * 1e-06 * mzs
    self  <- floor(mzs)
    all_mzs <- cbind(self - 1, self, self + 1)
    forms <- colnames(mz_mat)
    out <- list()
    for (i in seq_along(mzs)) {
      mz <- mzs[i]
      rt <- rts[i]
      rt_rank <- ranks[i]
      pos_all <- as.numeric(unique(unlist(tree[all_mzs[i, ]])))
      for (pos in pos_all) {
        id <- lib$cpd.lib$id[pos]
        mw_all <- mz_mat[pos, ]
        diffs <- abs(mw_all - mz)
        hit_inx <- which(diffs < tols[i])
        for (h in hit_inx) {
          out[[length(out) + 1L]] <- c(i, id, mz, rt, rt_rank, mw_all[h],
                                       forms[h], diffs[h])
        }
      }
    }
    out
  }

  res_p <- if (prep$mode != "negative") {
    one_polarity(pos_inx, lib$cpd.lib$adducts$positive, lib$cpd.tree$positive)
  } else list()
  res_n <- if (prep$mode != "positive") {
    one_polarity(!pos_inx, lib$cpd.lib$adducts$negative, lib$cpd.tree$negative)
  } else list()
  if (length(res_p) + length(res_n) == 0L) {
    .mmc_stop("no feature matched any library compound within ", ppm, " ppm.")
  }

  as_frame <- function(x) {
    data.frame(matrix(unlist(x), nrow = length(x), byrow = TRUE),
               stringsAsFactors = FALSE)
  }
  matched <- if (length(res_p) > 0L && length(res_n) > 0L) {
    rbind(as_frame(res_p), as_frame(res_n))
  } else if (length(res_n) > 0L) {
    as_frame(res_n)
  } else {
    as_frame(res_p)
  }
  matched <- matched[, c(3, 2, 7, 8, 4, 5)]
  colnames(matched) <- c("Query.Mass", "Matched.Compound", "Matched.Form",
                         "Mass.Diff", "Retention.Time", "RT.Rank")

  if (isTRUE(exclude_currency)) {
    drop <- matched[, "Matched.Compound"] %in%
      .MMC_MA_CURRENCY[nzchar(.MMC_MA_CURRENCY)]
    # Upstream keeps the table when currency accounts for every match.
    if (any(drop) && !all(drop)) {
      matched <- matched[!drop, , drop = FALSE]
      rownames(matched) <- NULL
    }
  }
  matched
}
