# =============================================================================
# Feature-correlation: "which metabolites behave like this one?"
# =============================================================================
# Answers a single question, on demand: given one feature, which other features
# rise and fall with it across the samples?
#
# This is deliberately a small, pure library rather than a pipeline stage. The
# feature of interest is chosen interactively in the Shiny GUI, so it cannot be
# known when the {targets} plan is built. The Shiny payload already carries
# everything needed (expr_norm, sample_meta, feature_annot, clust_partition), so
# nothing here changes the payload contract.
#
# Both Pearson and Spearman are returned for every partner. Pearson is the
# primary measure -- expr_work is already log2, and build_clustering_distance()
# also uses Pearson, so the ranking and the clustering agree on what "similar"
# means. Spearman rides along as a robustness check: when the two agree the hit
# is convincing, and when they diverge that is itself informative (an outlier,
# or a monotonic but non-linear relationship).
# =============================================================================


# ---- internal: pairwise-complete Pearson ------------------------------------

#' Pairwise-complete Pearson correlation of one vector against every matrix row
#'
#' Computes the correlation sufficient statistics directly over each pair's
#' shared (mutually non-missing) samples. Doing the algebra by hand rather than
#' calling \code{stats::cor(use = "pairwise.complete.obs")} buys two things we
#' need: the per-pair sample count and the per-pair variances fall out of the
#' same pass, so the zero-variance guard is inherently pairwise-complete aware;
#' and because \code{stats::cor()} is never called it cannot emit its "the
#' standard deviation is zero" warning, so we never have to reach for a blanket
#' \code{suppressWarnings()} that would also hide genuine warnings.
#'
#' A pair is undefined -- and returns \code{NA} rather than a number -- when
#' either side is constant *on the samples the pair shares*, even if that
#' feature varies perfectly well elsewhere in the matrix.
#'
#' @param x Numeric vector, the query profile (may contain NA).
#' @param mat Numeric matrix (features x samples), aligned to \code{x} by
#'   column position (may contain NA).
#' @param tol Relative tolerance for calling a shared-subset variance zero.
#'   Compared against the raw sum of squares, so it is scale-free.
#' @return list with \code{r} (correlation, NA where undefined) and \code{n}
#'   (shared sample count), both of length \code{nrow(mat)}.
#' @keywords internal
.pairwise_pearson <- function(x, mat, tol = 1e-9) {
  x <- as.numeric(x)
  mat <- as.matrix(mat)

  ox <- as.numeric(!is.na(x))
  xz <- x
  xz[is.na(xz)] <- 0

  om <- matrix(as.numeric(!is.na(mat)), nrow = nrow(mat))
  mz <- mat
  mz[is.na(mz)] <- 0

  n   <- as.vector(om %*% ox)
  sx  <- as.vector(om %*% (ox * xz))
  sxx <- as.vector(om %*% (ox * xz^2))
  sy  <- as.vector(mz %*% ox)
  syy <- as.vector((mz^2) %*% ox)
  sxy <- as.vector(mz %*% (ox * xz))

  # Guard the division before computing the centred moments.
  ok_n <- n > 2
  safe_n <- ifelse(ok_n, n, NA_real_)

  var_x <- sxx - sx^2 / safe_n
  var_y <- syy - sy^2 / safe_n
  cov_xy <- sxy - sx * sy / safe_n

  # Relative zero test: a constant profile leaves only floating-point dust
  # behind, which an `== 0` test would happily divide by.
  flat_x <- !is.finite(var_x) | var_x <= tol * pmax(abs(sxx), .Machine$double.eps)
  flat_y <- !is.finite(var_y) | var_y <= tol * pmax(abs(syy), .Machine$double.eps)

  r <- cov_xy / sqrt(var_x * var_y)
  r[!ok_n | flat_x | flat_y | !is.finite(r)] <- NA_real_

  # Floating point can nudge a perfect correlation just past 1.
  r <- pmin(pmax(r, -1), 1)

  list(r = r, n = n)
}


#' Two-sided p-value for a correlation coefficient
#'
#' Uses the t-transform \eqn{t = r\sqrt{(n-2)/(1-r^2)}} on \code{n - 2} degrees
#' of freedom. For Pearson under bivariate normality this is the exact test.
#' For Spearman it is the asymptotic approximation -- see the note in
#' \code{\link{correlate_feature_vs_all}}.
#'
#' @param r Numeric vector of correlation coefficients (NA allowed).
#' @param n Integer vector of sample counts, same length as \code{r}.
#' @return Numeric vector of p-values, NA where \code{r} is NA or \code{n < 3}.
#' @keywords internal
.correlation_pvalue <- function(r, n) {
  p <- rep(NA_real_, length(r))
  ok <- is.finite(r) & is.finite(n) & n > 2
  if (!any(ok)) return(p)

  r_ok <- r[ok]
  n_ok <- n[ok]
  out <- numeric(length(r_ok))

  # |r| == 1 makes 1 - r^2 vanish; take the limit directly instead of dividing
  # by zero and hoping pt() copes with an Inf.
  perfect <- abs(r_ok) >= 1 - 1e-12
  out[perfect] <- 0

  if (any(!perfect)) {
    rp <- r_ok[!perfect]
    np <- n_ok[!perfect]
    t_stat <- rp * sqrt((np - 2) / (1 - rp^2))
    out[!perfect] <- 2 * stats::pt(-abs(t_stat), df = np - 2)
  }

  p[ok] <- out
  p
}


#' Pairwise-complete Spearman correlation of one vector against every matrix row
#'
#' Ranks must be recomputed within each pair's shared samples, so a single
#' global rank transform is only valid when the matrix has no missing values.
#' That is the common case (the Shiny payload median-fills NAs), so it gets a
#' vectorised fast path; otherwise we fall back to a loop over features. The
#' loop is over features, not over pairs, so it is a few thousand \code{rank()}
#' calls on short vectors.
#'
#' Because all-tied shared values rank to a constant, the zero-variance guard in
#' \code{.pairwise_pearson()} covers Spearman too, on the same shared subset.
#'
#' @param x Numeric vector, the query profile (may contain NA).
#' @param mat Numeric matrix (features x samples) aligned to \code{x}.
#' @return list with \code{r} (Spearman's rho, NA where undefined) and \code{n}.
#' @keywords internal
.pairwise_spearman <- function(x, mat) {
  x <- as.numeric(x)
  mat <- as.matrix(mat)

  if (!anyNA(mat)) {
    # Every feature shares exactly the same samples with x, so subset once and
    # rank once.
    keep <- !is.na(x)
    if (sum(keep) < 3) {
      return(list(r = rep(NA_real_, nrow(mat)), n = rep(sum(keep), nrow(mat))))
    }
    xs <- x[keep]
    ms <- mat[, keep, drop = FALSE]
    # apply() over rows returns samples-by-features; transpose puts it back.
    mr <- t(apply(ms, 1, rank))
    if (nrow(ms) == 1L) mr <- matrix(mr, nrow = 1L)
    return(.pairwise_pearson(rank(xs), mr))
  }

  n_feat <- nrow(mat)
  r <- rep(NA_real_, n_feat)
  n <- integer(n_feat)

  for (i in seq_len(n_feat)) {
    y <- mat[i, ]
    keep <- !is.na(x) & !is.na(y)
    n[i] <- sum(keep)
    if (n[i] < 3) next

    res <- .pairwise_pearson(rank(x[keep]), matrix(rank(y[keep]), nrow = 1L))
    r[i] <- res$r[1]
  }

  list(r = r, n = n)
}


#' Benjamini-Hochberg adjustment restricted to genuinely testable pairs
#'
#' Adjusting across rows that were never tested inflates every q-value in the
#' table. Untestable rows keep \code{NA} so they stay visible to the reader as
#' "not tested" rather than quietly vanishing or masquerading as non-significant.
#'
#' @param p Numeric vector of raw p-values.
#' @param testable Logical vector marking rows eligible for adjustment.
#' @return Numeric vector of adjusted p-values, NA outside \code{testable}.
#' @keywords internal
.adjust_testable <- function(p, testable) {
  padj <- rep(NA_real_, length(p))
  idx <- which(testable & is.finite(p))
  if (length(idx) > 0L) {
    padj[idx] <- stats::p.adjust(p[idx], method = "BH")
  }
  padj
}


# ---- public: the engine ------------------------------------------------------

#' Rank every feature by how closely it tracks a chosen feature
#'
#' Correlates one feature's profile against every other feature in the matrix
#' and returns them ranked by strength of association. Intended to answer a
#' researcher's question directly: "which other metabolites behave like this
#' one?"
#'
#' Both coefficients are always returned. **Pearson is the primary measure** --
#' the matrix is already log2-transformed, and \code{build_clustering_distance()}
#' also uses Pearson, so this ranking and the pipeline's clustering agree on what
#' "similar" means. **Spearman is a supporting measure**: agreement between the
#' two is a convincing hit, while a large gap usually means either a single
#' outlying sample is driving Pearson, or the relationship is monotonic but not
#' linear.
#'
#' Correlations are computed across whatever samples are given. Pass an
#' already-filtered matrix -- see \code{\link{prepare_correlation_matrix}}, since
#' \code{expr_work} still contains QC/blank/pool samples.
#'
#' @section Missing values:
#' Metabolomics has no imputation stage, so real NAs reach the matrix. Each pair
#' uses its own mutually-observed samples and reports that count in
#' \code{n_used}; pairs sharing fewer than \code{min_n} samples are reported with
#' \code{NA} statistics rather than dropped, so a correlation resting on four
#' points can never quietly top the table. A pair is also \code{NA} when either
#' feature is constant *on the samples that pair shares*, which can happen even
#' when both vary across the full matrix.
#'
#' @section A note on the Spearman p-value:
#' \code{spearman_pvalue} is an **asymptotic approximation**, not an exact
#' permutation test: the Pearson t-transform is applied to the rank correlation,
#' matching \code{stats::cor.test(method = "spearman", exact = FALSE)}. It
#' degrades at small n and in the presence of ties. \code{pearson_pvalue} is the
#' exact t-test under bivariate normality. Read the two accordingly.
#'
#' @param expr_mat Numeric matrix (features x samples), rownames = feature IDs.
#'   Expected to be normalized and log2-scaled, i.e. \code{pre$expr_work} or the
#'   payload's \code{expr_norm}.
#' @param feature_id Single feature ID; must be present in
#'   \code{rownames(expr_mat)}.
#' @param min_n Minimum shared samples for a pair to be tested. Default 5.
#' @param top_n Optionally keep only the strongest \code{top_n} partners. Applied
#'   last -- after the query feature is removed and the table is sorted -- so it
#'   always returns that many genuine partners.
#' @return A data.frame, one row per other feature, ordered by
#'   \code{abs(pearson_r)} descending with ties broken by \code{feature_id} so
#'   the order is fully deterministic. Untested pairs sort last. Columns:
#'   \describe{
#'     \item{feature_id}{The partner feature.}
#'     \item{pearson_r, pearson_pvalue, pearson_padj}{Primary measure; BH-adjusted
#'       across testable pairs only.}
#'     \item{spearman_rho, spearman_pvalue, spearman_padj}{Supporting measure;
#'       p-value is approximate (see above).}
#'     \item{n_used}{Samples shared with the query feature.}
#'     \item{direction}{\code{"positive"}, \code{"negative"}, \code{"none"} for an
#'       exactly zero correlation, or \code{NA} when \code{pearson_r} is \code{NA}.}
#'   }
#' @seealso \code{\link{prepare_correlation_matrix}} to filter out QC samples
#'   first, and \code{\link{plot_feature_correlation_profiles}} to see the
#'   pattern rather than only the number.
#' @examples
#' m <- matrix(c(1, 2, 3, 4, 5,
#'               2, 4, 6, 8, 10,
#'               5, 4, 3, 2, 1),
#'             nrow = 3, byrow = TRUE,
#'             dimnames = list(c("F1", "F2", "F3"), paste0("S", 1:5)))
#' correlate_feature_vs_all(m, "F1", min_n = 3)
correlate_feature_vs_all <- function(expr_mat, feature_id, min_n = 5, top_n = NULL) {
  if (!is.matrix(expr_mat) && !is.data.frame(expr_mat)) {
    stop("correlate_feature_vs_all(): `expr_mat` must be a matrix or data frame ",
         "of features x samples, got ", class(expr_mat)[1], ".", call. = FALSE)
  }
  expr_mat <- as.matrix(expr_mat)

  if (is.null(rownames(expr_mat))) {
    stop("correlate_feature_vs_all(): `expr_mat` needs feature IDs as rownames.",
         call. = FALSE)
  }
  if (!is.character(feature_id) || length(feature_id) != 1L || is.na(feature_id)) {
    stop("correlate_feature_vs_all(): `feature_id` must be a single ",
         "non-missing feature ID.", call. = FALSE)
  }
  if (nrow(expr_mat) < 2L) {
    stop("correlate_feature_vs_all(): need at least 2 features to correlate ",
         "against, `expr_mat` has ", nrow(expr_mat), ".", call. = FALSE)
  }

  feature_ids <- rownames(expr_mat)
  if (!feature_id %in% feature_ids) {
    near <- tryCatch(
      utils::head(agrep(feature_id, feature_ids, max.distance = 0.25,
                        value = TRUE, ignore.case = TRUE), 5),
      error = function(e) character(0)
    )
    stop(
      "correlate_feature_vs_all(): feature '", feature_id,
      "' is not in the matrix (", length(feature_ids), " features available).",
      if (length(near) > 0) {
        paste0("\n  Did you mean: ", paste(near, collapse = ", "), "?")
      } else {
        "\n  Check that you are passing a feature_id, not a display name."
      },
      call. = FALSE
    )
  }
  if (!is.numeric(min_n) || length(min_n) != 1L || is.na(min_n) || min_n < 3) {
    stop("correlate_feature_vs_all(): `min_n` must be a single number >= 3 ",
         "(a correlation needs at least 3 points).", call. = FALSE)
  }

  x <- expr_mat[feature_id, ]

  pea <- .pairwise_pearson(x, expr_mat)
  spe <- .pairwise_spearman(x, expr_mat)

  res <- data.frame(
    feature_id      = feature_ids,
    pearson_r       = pea$r,
    spearman_rho    = spe$r,
    n_used          = as.integer(pea$n),
    stringsAsFactors = FALSE
  )

  # Drop the query feature BEFORE adjusting, so it never contributes to the
  # multiple-testing burden of its own partners.
  res <- res[res$feature_id != feature_id, , drop = FALSE]

  # Below min_n the coefficient may be computable but is not trustworthy; blank
  # it so it cannot be ranked or adjusted on.
  under_powered <- res$n_used < min_n
  res$pearson_r[under_powered] <- NA_real_
  res$spearman_rho[under_powered] <- NA_real_

  res$pearson_pvalue  <- .correlation_pvalue(res$pearson_r, res$n_used)
  res$spearman_pvalue <- .correlation_pvalue(res$spearman_rho, res$n_used)

  # Each coefficient gets its own testable set: a pair can be defined under one
  # and undefined under the other.
  res$pearson_padj <- .adjust_testable(
    res$pearson_pvalue,
    res$n_used >= min_n & is.finite(res$pearson_r) & is.finite(res$pearson_pvalue)
  )
  res$spearman_padj <- .adjust_testable(
    res$spearman_pvalue,
    res$n_used >= min_n & is.finite(res$spearman_rho) & is.finite(res$spearman_pvalue)
  )

  res$direction <- ifelse(
    is.na(res$pearson_r), NA_character_,
    ifelse(res$pearson_r > 0, "positive",
           ifelse(res$pearson_r < 0, "negative", "none"))
  )

  # Explicit sort key: untested pairs last, ties broken by ID so two runs on the
  # same data can never disagree about row order.
  sort_key <- -abs(res$pearson_r)
  sort_key[is.na(sort_key)] <- Inf
  res <- res[order(sort_key, res$feature_id), , drop = FALSE]

  res <- res[, c("feature_id",
                 "pearson_r", "pearson_pvalue", "pearson_padj",
                 "spearman_rho", "spearman_pvalue", "spearman_padj",
                 "n_used", "direction")]
  rownames(res) <- NULL

  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1) {
      stop("correlate_feature_vs_all(): `top_n` must be a single number >= 1, ",
           "or NULL to keep every partner.", call. = FALSE)
    }
    res <- utils::head(res, as.integer(top_n))
  }

  res
}


# ---- public: the filtering boundary ------------------------------------------

#' Drop QC/blank/pool samples before correlating
#'
#' \code{pre$expr_work} -- and therefore the Shiny payload's \code{expr_norm} --
#' still contains technical samples; every consumer in this pipeline filters at
#' its own point of use rather than relying on upstream filtering. Correlating
#' without filtering lets pooled QC samples, which sit at the average of
#' everything, distort every coefficient in the table.
#'
#' This wrapper exists so callers get that right with one line instead of
#' re-implementing the QC/blank/pool naming conventions, while
#' \code{\link{correlate_feature_vs_all}} stays pure and simply trusts the matrix
#' it is handed.
#'
#' @param expr_mat Numeric matrix (features x samples).
#' @param sample_meta Sample metadata; one row per sample. Sample IDs may live in
#'   a column or in the rownames.
#' @param condition_col Column naming the experimental group. From the Shiny
#'   payload this is \code{payload$group}.
#' @param sample_col Column holding sample IDs. Defaults to \code{"sample_id"};
#'   if absent, the rownames of \code{sample_meta} are used.
#' @param qc_flag_column Optional column that flags technical samples, for
#'   projects that set \code{qc.qc_flag_column}. The Shiny payload does not carry
#'   this setting, so such projects must pass it explicitly -- otherwise those
#'   samples are not recognised as technical.
#' @return The filtered matrix, columns restricted to biological samples.
#' @seealso \code{filter_to_biological}, which does the actual detection.
prepare_correlation_matrix <- function(expr_mat, sample_meta, condition_col,
                                       sample_col = "sample_id",
                                       qc_flag_column = NULL) {
  expr_mat <- as.matrix(expr_mat)
  meta <- as.data.frame(sample_meta, stringsAsFactors = FALSE)

  if (!sample_col %in% colnames(meta)) {
    if (is.null(rownames(meta))) {
      stop("prepare_correlation_matrix(): no '", sample_col, "' column in ",
           "`sample_meta` and no rownames to fall back on. Pass `sample_col`.",
           call. = FALSE)
    }
    meta[[sample_col]] <- rownames(meta)
  }
  if (!condition_col %in% colnames(meta)) {
    stop("prepare_correlation_matrix(): condition column '", condition_col,
         "' not found in `sample_meta`. Available: ",
         paste(colnames(meta), collapse = ", "), call. = FALSE)
  }

  # filter_to_biological() subsets the matrix by sample ID, so meta must
  # describe exactly the columns present, in that order.
  idx <- match(colnames(expr_mat), as.character(meta[[sample_col]]))
  if (anyNA(idx)) {
    missing <- colnames(expr_mat)[is.na(idx)]
    stop("prepare_correlation_matrix(): ", length(missing), " matrix column(s) ",
         "have no row in `sample_meta`, e.g. ",
         paste(utils::head(missing, 5), collapse = ", "), ".", call. = FALSE)
  }
  meta <- meta[idx, , drop = FALSE]

  bio <- filter_to_biological(
    mat            = expr_mat,
    meta           = meta,
    condition_col  = condition_col,
    sample_col     = sample_col,
    label          = "feature correlation",
    qc_flag_column = qc_flag_column
  )

  if (ncol(bio$mat) < 3L) {
    stop("prepare_correlation_matrix(): only ", ncol(bio$mat), " biological ",
         "sample(s) remain after removing QC/blank/pool -- not enough to ",
         "correlate. Check `condition_col` ('", condition_col, "') names the ",
         "right column.", call. = FALSE)
  }

  bio$mat
}
