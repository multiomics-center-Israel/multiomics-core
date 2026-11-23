# R/domain/metabolomics/06o_mummichog_ma_gsea.R
#
# MetaboAnalyst-compatible mummichog GSEA on the engine's own joint-polarity
# empirical compounds.
#
# MetaboAnalystR's peaks GSEA (.compute.mummichog.RT.fgsea) needs exactly two
# inputs, and both are already produced — and validated — by this engine:
#   * mSetObj$ec_exp: each EC's score, the max of its member features' t-scores
#     (.search.compoundLib; here ecs$ec_exp from mmc_ma_build_ecs()). In mixed
#     mode an EC's members can come from both polarities, so the score is taken
#     over both;
#   * mSetObj$pathways$emp_cpds: pathway -> ECs after the minLib filter
#     (.setup.psea.library; here mmc_ma_pathway_ecs()), i.e. the very pathway
#     sets the ORA tests.
# The statistical part is reused unchanged: mmc_gsea_ranked_inputs() (upstream
# L3470-L3483) and mmc_gsea_metaboanalyst() (util_fgsea.R .run_fgsea_inner),
# both validated against the pinned reference by the GSEA parity fixture. What
# the python adapter run_mummichog_gsea() adds — reading python mummichog's EC
# tables, its m/z merge (mmc_ec_scores) and pathway mapping
# (mmc_pathway_ec_sets) — is not used: MetaboAnalyst builds these inputs itself,
# and this engine reproduces that.
#
# ===========================================================================
# BEHAVIOURAL REFERENCE — xia-lab/MetaboAnalystR @ 398476ae, R/peaks_to_function.R
# .compute.mummichog.RT.fgsea() L3457-L3573:
#   res.mat = Pathway_Total (model compounds), Hits (pathway compounds that
#   matched), P_val, P_adj (BH over the tested pathways), NES; ordered by P_val;
#   signif(, 4).
# Not reproduced, documented: its result CSV binds the EC-hit column in the
# library's pathway order to rows ordered by p-value, so the hits do not belong
# to their rows; the EC hits here are each pathway's own. When exactly one
# pathway is tested, upstream's signif(as.matrix(res.mat[ord.inx, ])) turns the
# row into a column; it is kept a row here.


#' MetaboAnalyst-compatible GSEA over the engine's empirical compounds
#'
#' @param ecs      Output of `mmc_ma_build_ecs()` (uses `ec_exp`, `cpd_exp_dict`).
#' @param pathways Output of `mmc_ma_pathway_ecs()`.
#' @param n_perm,min_size,max_size,seed GSEA settings; defaults are
#'   MetaboAnalyst's (100, 1, Inf, 123). Any other value is reported in
#'   `deviations`.
#' @param metric_label Provenance: how the features' scores were chosen.
#' @return NULL when no pathway passes the size filter; otherwise a list with
#'   `resmat` (MetaboAnalyst's 5-column table, rounded and ordered as upstream),
#'   `fgsea` (the engine's full per-pathway result), `table` (the report table:
#'   one row per tested pathway in `resmat` order, with `Pathway`, model size,
#'   matched compounds, detected ECs, tested size, `ES`, `NES`, `P.Value`,
#'   `padj` (BH), the leading edge, the pathway's ECs and the ES-fallback
#'   flags), `params`, `deviations`, `reference`, `metric_label`, `nes_note`,
#'   `n_ec`, `n_positions`, `n_pathways`.
#' @examples
#' \dontrun{
#' g <- mmc_ma_gsea(ecs, pathways)
#' g$table
#' }
mmc_ma_gsea <- function(ecs, pathways, n_perm = .MMC_GSEA_DEFAULTS$n_perm,
                        min_size = .MMC_GSEA_DEFAULTS$min_size,
                        max_size = .MMC_GSEA_DEFAULTS$max_size,
                        seed = .MMC_GSEA_DEFAULTS$seed,
                        metric_label = NA_character_) {
  sets <- stats::setNames(pathways$emp_cpds, pathways$name)
  ranked <- mmc_gsea_ranked_inputs(ecs$ec_exp)
  res <- mmc_gsea_metaboanalyst(sets, ranked$stats, ranked$ranks,
                                n_perm = n_perm, min_size = min_size,
                                max_size = max_size, gsea_param = 1,
                                seed = seed)
  if (nrow(res) == 0) {
    message("mummichog GSEA (metaboanalyst): no pathway passed the size filter")
    return(NULL)
  }

  # .compute.mummichog.RT.fgsea()'s result table.
  idx <- match(res$pathway, pathways$name)
  total_cpds <- names(ecs$cpd_exp_dict)
  matched_size <- unlist(lapply(lapply(pathways$cpds, function(x)
    intersect(x, total_cpds)), length))[idx]
  path_size <- unlist(lapply(pathways$cpds, length))[idx]
  res_mat <- matrix(0, nrow = length(res$pathway), ncol = 5)
  res_mat[, 1] <- path_size
  res_mat[, 2] <- matched_size
  res_mat[, 3] <- res$pval
  res_mat[, 4] <- res$padj
  res_mat[, 5] <- res$NES
  rownames(res_mat) <- res$pathway
  colnames(res_mat) <- c("Pathway_Total", "Hits", "P_val", "P_adj", "NES")
  ord_inx <- order(res_mat[, 3])
  res_mat <- signif(as.matrix(res_mat[ord_inx, , drop = FALSE]), 4)

  # Each pathway's own ECs, as upstream's pathwaysFiltered holds them.
  ec_hits <- vapply(rownames(res_mat), function(p) {
    hit <- stats::na.omit(match(sets[[p]], names(ranked$ranks)))
    paste(names(ranked$ranks)[hit], collapse = ";")
  }, character(1), USE.NAMES = FALSE)
  r <- res[match(rownames(res_mat), res$pathway), , drop = FALSE]
  table <- data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    "Pathway"                         = rownames(res_mat),
    "Pathway size (model compounds)"  = unname(res_mat[, "Pathway_Total"]),
    "Matched compounds"               = unname(res_mat[, "Hits"]),
    "Detected ECs"                    = unname(lengths(sets[rownames(res_mat)])),
    "Tested size"                     = r$size,
    "ES"                              = round(r$ES, 4),
    "NES"                             = unname(res_mat[, "NES"]),
    "P.Value"                         = unname(res_mat[, "P_val"]),
    "padj"                            = unname(res_mat[, "P_adj"]),
    "Leading-edge EmpiricalCompounds" = r$leadingEdge,
    "EmpiricalCompounds"              = ec_hits,
    "ES defaulted by reference implementation" = r$es_defaulted,
    "ES fallback reason"                       = r$es_reason)

  d <- .MMC_GSEA_DEFAULTS
  deviations <- character(0)
  for (k in c("n_perm", "min_size", "seed")) {
    v <- get(k)
    if (!identical(as.numeric(v), as.numeric(d[[k]]))) {
      deviations <- c(deviations, sprintf("%s = %s (MetaboAnalyst default %s)",
                                          k, v, d[[k]]))
    }
  }
  if (!identical(as.numeric(max_size), as.numeric(d$max_size))) {
    deviations <- c(deviations, sprintf("max_size = %s (MetaboAnalyst default Inf)",
                                        max_size))
  }
  if (isTRUE(attr(res, "mmc_single_pathway_rng"))) {
    deviations <- c(deviations, paste0(
      "one permutable pathway: the permutation loop is seeded here, where the ",
      "reference leaves it to a non-deterministic worker RNG"))
  }
  if (isTRUE(attr(res, "mmc_unpermutable") > 0)) {
    deviations <- c(deviations, sprintf(paste0(
      "%d pathway(s) tested but not permuted (ES = 0, p = 1) where the ",
      "reference implementation aborts"), attr(res, "mmc_unpermutable")))
  }
  if (nrow(res_mat) == 1L) {
    deviations <- c(deviations, paste0(
      "one tested pathway: kept as a table row (upstream's signif(as.matrix()) ",
      "would turn it into a column)"))
  }
  fgsea_note <- .mmc_gsea_check_fgsea(warn = FALSE)
  if (!is.null(fgsea_note)) deviations <- c(deviations, fgsea_note)

  list(
    resmat = res_mat, fgsea = res, table = table,
    params = list(n_perm = n_perm, min_size = min_size, max_size = max_size,
                  gsea_param = 1, seed = seed,
                  fgsea_version = as.character(utils::packageVersion("fgsea"))),
    deviations = deviations,
    reference = list(repo = .MMC_GSEA_REF_REPO, commit = .MMC_GSEA_REF_COMMIT,
                     verified = paste0("EC scores and pathway EC sets (the ",
                                       "engine's own ECs) and the GSEA ",
                                       "calculation from them")),
    metric_label = metric_label, nes_note = mmc_gsea_nes_note(),
    n_ec = length(ecs$ec_exp), n_positions = length(ranked$stats),
    n_pathways = nrow(table))
}
