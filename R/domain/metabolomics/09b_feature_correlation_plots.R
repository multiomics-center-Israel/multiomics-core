# =============================================================================
# Feature-correlation plots
# =============================================================================
# Visual companion to R/domain/metabolomics/09_feature_correlation.R. Split out
# per the convention used by 06e_mummichog_plots.R and 07b_multigsea_plots.R:
# the engine computes, the plots render.
# =============================================================================


#' Overlay the group-mean profiles of a feature and its top correlates
#'
#' A correlation coefficient says two features move together; this shows *how*.
#' Profiles are z-scored per feature so a low-abundance metabolite and a
#' high-abundance one can be compared on the same axes -- the shape is the point,
#' not the height.
#'
#' @param expr_mat Numeric matrix (features x samples), already filtered to
#'   biological samples.
#' @param meta Sample metadata, one row per column of \code{expr_mat}.
#' @param feature_id The chosen feature; drawn highlighted.
#' @param partner_ids Feature IDs to overlay, typically the head of
#'   \code{\link{correlate_feature_vs_all}}. IDs not present in the matrix are
#'   dropped with a message.
#' @param group_col Metadata column defining the groups on the x-axis. From the
#'   Shiny payload this is \code{payload$group}.
#' @param sample_col Metadata column holding sample IDs. Default \code{NULL},
#'   meaning take them from \code{rownames(meta)} -- see
#'   \code{\link{prepare_correlation_matrix}} for why \code{"sample_id"} is the
#'   wrong default for the shipped templates.
#' @param label_map Optional named character vector mapping feature IDs to
#'   display names. May be partial or absent; any feature without an entry falls
#'   back to its raw ID.
#' @section Groups with no observations:
#' Once missingness is restored, a feature can have every sample in one group
#' unobserved. The group mean is then undefined, and this draws a \strong{gap}
#' in the line rather than a point. It never interpolates or substitutes: a
#' fabricated midpoint would show the reader a measurement that was never made.
#'
#' @return A \code{ggplot} object. Nothing is written to disk -- saving belongs
#'   in a target, not in a plotting function.
plot_feature_correlation_profiles <- function(expr_mat, meta,
                                              feature_id, partner_ids,
                                              group_col,
                                              sample_col = NULL,
                                              label_map = NULL) {
  expr_mat <- as.matrix(expr_mat)

  if (!feature_id %in% rownames(expr_mat)) {
    stop("plot_feature_correlation_profiles(): feature '", feature_id,
         "' is not in the matrix.", call. = FALSE)
  }

  partner_ids <- setdiff(unique(as.character(partner_ids)), feature_id)
  known <- partner_ids %in% rownames(expr_mat)
  if (any(!known)) {
    message(sprintf(
      "feature correlation plot: dropping %d partner ID(s) not present in the matrix",
      sum(!known)
    ))
  }
  partner_ids <- partner_ids[known]

  # build_group_means_from_effects() reads the group column via
  # get_clustering_group_col(), which insists on cfg$clustering$group_col. This
  # plot has nothing to do with clustering and callers (notably the GUI) only
  # have the column name, so hand it a minimal shim rather than demanding a
  # full clustering config.
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  resolved <- .resolve_sample_ids(meta, expr_mat, sample_col,
                                  context = "plot_feature_correlation_profiles")

  cfg_shim <- list(
    clustering = list(group_col = group_col),
    effects    = list(samples = resolved$sample_col)
  )
  gm <- build_group_means_from_effects(expr_mat, resolved$meta, cfg_shim)

  ids <- c(feature_id, partner_ids)
  gmeans <- gm$group_means[ids, , drop = FALSE]

  # A group in which every sample is unobserved makes rowMeans(na.rm = TRUE)
  # return NaN. Normalise that (and any Inf) to NA deliberately, so geom_line()
  # renders an honest break instead of us relying on ggplot to drop a
  # not-a-number by accident. Never interpolate: an invented midpoint would draw
  # a measurement nobody took.
  gaps <- !is.finite(gmeans)
  if (any(gaps)) {
    gmeans[gaps] <- NA_real_
    message(sprintf(
      "feature correlation plot: %d feature x group cell(s) have no observed values; drawn as gaps",
      sum(gaps)
    ))
  }

  # zscore_rows() already maps a zero-variance row to zeros rather than NaN, so
  # a flat profile draws as a flat line instead of vanishing. It cannot rescue a
  # row that is entirely unobserved though -- its row mean is NaN, which then
  # propagates -- so normalise once more on the way out. After this the plotted
  # values are finite or NA, never NaN/Inf.
  z <- zscore_rows(gmeans)
  z[!is.finite(z)] <- NA_real_

  group_levels <- as.character(gm$group_levels)
  long <- data.frame(
    feature_id = rep(ids, times = length(group_levels)),
    group      = rep(group_levels, each = length(ids)),
    value      = as.vector(z),
    stringsAsFactors = FALSE
  )
  long$group <- factor(long$group, levels = group_levels)
  long$is_query <- long$feature_id == feature_id

  # An unnamed label_map would silently index by position and mislabel every
  # feature, which is worse than no labels at all -- fall back to raw IDs.
  display <- function(id) {
    if (is.null(label_map) || is.null(names(label_map))) return(id)
    lbl <- as.character(label_map[id])
    ifelse(is.na(lbl) | !nzchar(lbl), id, lbl)
  }
  long$label <- display(long$feature_id)

  query_label <- display(feature_id)
  partners <- long[!long$is_query, , drop = FALSE]
  query <- long[long$is_query, , drop = FALSE]

  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = group, y = value))

  if (nrow(partners) > 0) {
    p <- p + ggplot2::geom_line(
      data = partners,
      mapping = ggplot2::aes(group = feature_id, colour = label),
      linewidth = 0.6, alpha = 0.8
    )
  }

  p +
    ggplot2::geom_line(
      data = query,
      mapping = ggplot2::aes(group = feature_id),
      colour = "black", linewidth = 1.6
    ) +
    ggplot2::geom_point(data = query, colour = "black", size = 2.4) +
    ggplot2::labs(
      title    = sprintf("Features tracking %s", query_label),
      subtitle = sprintf("%d correlate(s) overlaid; %s drawn in black",
                         length(partner_ids), query_label),
      x        = group_col,
      y        = "Group mean (z-scored per feature)",
      colour   = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
