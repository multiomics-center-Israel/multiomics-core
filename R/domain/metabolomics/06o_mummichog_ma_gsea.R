# R/domain/metabolomics/06o_mummichog_ma_gsea.R
#
# Mummichog GSEA on the MetaboAnalyst-compatible engine's own joint-polarity
# empirical compounds.
#
# Two methods (`gsea_method`, see 06g): the default "fgsea" hands the two inputs
# below to mmc_gsea_fgsea(), i.e. standard fgsea::fgsea() (multilevel) through
# the pipeline's shared run_fgsea_core(), the same function the python engine's
# GSEA uses. "metaboanalyst_legacy" reproduces
# MetaboAnalystR's own peaks GSEA, with its known defects, and is described in
# the rest of this header. The ECs, their scores and the pathway sets are the
# engine's in both cases; only the test differs.
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
# The legacy statistical part is reused unchanged: mmc_gsea_ranked_inputs() (upstream
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


#' GSEA over the MetaboAnalyst-compatible engine's empirical compounds
#'
#' `"fgsea"` (default) is standard `fgsea::fgsea()` (multilevel) through
#' `mmc_gsea_fgsea()`; `"metaboanalyst_legacy"` reproduces MetaboAnalystR's
#' peaks GSEA and warns. Both test the same EC scores (`ecs$ec_exp`) and
#' pathway EC sets (`pathways$emp_cpds`).
#'
#' @param ecs      Output of `mmc_ma_build_ecs()` (uses `ec_exp`, `cpd_exp_dict`).
#' @param pathways Output of `mmc_ma_pathway_ecs()`.
#' @param n_perm Legacy method only: its number of permutations
#'   (MetaboAnalyst's default 100). Must not be given with the default method,
#'   which has no such setting.
#' @param min_size,max_size,seed GSEA settings; defaults are MetaboAnalyst's
#'   (1, Inf, 123). For the legacy method any non-default value, `n_perm`
#'   included, is reported in `deviations`.
#' @param metric_label Provenance: how the features' scores were chosen.
#' @param method `"fgsea"` (default) or `"metaboanalyst_legacy"`; resolve it
#'   from config with `mmc_gsea_method()`.
#' @return NULL when no pathway passes the size filter; otherwise a list with
#'   `fgsea` (the method's full per-pathway result), `table` (the report table:
#'   one row per tested pathway ordered by p-value, with `Pathway`, model size,
#'   matched compounds, detected ECs, tested size, `ES`, `NES`, `P.Value`,
#'   `padj` (BH), the leading edge, the pathway's ECs and the ES-fallback
#'   flags), `method`, `method_label`, `params`, `deviations` (legacy only;
#'   empty for `"fgsea"`), `reference`, `metric_label`, `nes_note`, `ties`
#'   (`mmc_gsea_tie_summary()`; NULL for legacy), `n_p_unavailable`,
#'   `log2err` (fgsea's expected error of each log2 p-value, named by pathway,
#'   as `run_mummichog_gsea()` returns it; NULL for legacy) and
#'   `log2err_summary` (`mmc_gsea_log2err_summary()`; NULL for legacy), `n_ec`,
#'   `n_positions`, `n_pathways`. The legacy method also returns `resmat`
#'   (MetaboAnalyst's 5-column table, rounded and ordered as upstream), and its
#'   `table` carries NES, p and adjusted p at MetaboAnalyst's 4 significant
#'   digits; `"fgsea"` reports them unrounded.
#' @examples
#' \dontrun{
#' g <- mmc_ma_gsea(ecs, pathways)
#' g$table
#' }
mmc_ma_gsea <- function(ecs, pathways, n_perm = .MMC_GSEA_DEFAULTS$n_perm,
                        min_size = .MMC_GSEA_DEFAULTS$min_size,
                        max_size = .MMC_GSEA_DEFAULTS$max_size,
                        seed = .MMC_GSEA_DEFAULTS$seed,
                        metric_label = NA_character_,
                        method = .MMC_GSEA_METHODS[[1]]) {
  method <- mmc_gsea_method(list(gsea_method = method))
  if (identical(method, "fgsea") && !missing(n_perm)) {
    .mmc_stop("n_perm is the legacy method's number of permutations and was ",
              "given with method \"fgsea\", which has no such setting. Drop ",
              "n_perm. (method = \"metaboanalyst_legacy\" takes it, but is a ",
              "different analysis.)")
  }
  sets <- stats::setNames(pathways$emp_cpds, pathways$name)
  fgsea_version <- as.character(utils::packageVersion("fgsea"))

  if (identical(method, "fgsea")) {
    res <- mmc_gsea_fgsea(sets, ecs$ec_exp, min_size = min_size,
                          max_size = max_size, seed = seed)
    if (nrow(res) == 0) {
      message("mummichog GSEA (fgsea): no pathway passed the size filter")
      return(NULL)
    }
    res <- res[order(res$pval, na.last = TRUE), , drop = FALSE]
    idx <- match(res$pathway, pathways$name)
    ranked_ecs <- names(ecs$ec_exp)[is.finite(ecs$ec_exp)]
    table <- data.frame(
      check.names = FALSE, stringsAsFactors = FALSE,
      "Pathway"                         = res$pathway,
      "Pathway size (model compounds)"  = as.numeric(lengths(pathways$cpds)[idx]),
      "Matched compounds"               = as.numeric(vapply(
        pathways$cpds[idx], function(x)
          length(intersect(x, names(ecs$cpd_exp_dict))), integer(1))),
      "Detected ECs"                    = as.numeric(unname(lengths(sets[res$pathway]))),
      "Tested size"                     = res$size,
      "ES"                              = round(res$ES, 4),
      "NES"                             = round(res$NES, 4),
      "P.Value"                         = res$pval,
      "padj"                            = res$padj,
      "Leading-edge EmpiricalCompounds" = res$leadingEdge,
      "EmpiricalCompounds"              = vapply(sets[res$pathway], function(p)
        paste(p[p %in% ranked_ecs], collapse = ";"), character(1),
        USE.NAMES = FALSE),
      "ES defaulted by reference implementation" = res$es_defaulted,
      "ES fallback reason"                       = res$es_reason)
    rownames(table) <- NULL
    return(list(
      fgsea = res, table = table, method = method,
      method_label = mmc_gsea_method_label(method, fgsea_version),
      # `nperm_simple`, not `n_perm_simple`: see run_mummichog_gsea().
      params = list(nperm_simple = .MMC_GSEA_NPERM_SIMPLE, min_size = min_size,
                    max_size = max_size, gsea_param = 1, seed = seed,
                    fgsea_version = fgsea_version),
      deviations = character(0),
      reference = list(method = "fgsea::fgsea (fgseaMultilevel)",
                       package = "fgsea", version = fgsea_version),
      metric_label = metric_label, nes_note = mmc_gsea_nes_note(method),
      ties = attr(res, "mmc_ties"),
      n_p_unavailable = sum(is.na(res$pval)),
      log2err = stats::setNames(as.numeric(res$log2err), res$pathway),
      log2err_summary = mmc_gsea_log2err_summary(res),
      n_ec = length(ecs$ec_exp), n_positions = length(ranked_ecs),
      n_pathways = nrow(table)))
  }

  # ---- legacy method: MetaboAnalystR's peaks GSEA, reproduced ----------------
  warning(.mmc_gsea_legacy_warning(), call. = FALSE)
  ranked <- mmc_gsea_ranked_inputs(ecs$ec_exp)
  res <- mmc_gsea_metaboanalyst(sets, ranked$stats, ranked$ranks,
                                n_perm = n_perm, min_size = min_size,
                                max_size = max_size, gsea_param = 1,
                                seed = seed)
  if (nrow(res) == 0) {
    message("mummichog GSEA (metaboanalyst_legacy): no pathway passed the size filter")
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
    resmat = res_mat, fgsea = res, table = table, method = method,
    method_label = mmc_gsea_method_label(method, fgsea_version),
    params = list(n_perm = n_perm, min_size = min_size, max_size = max_size,
                  gsea_param = 1, seed = seed, fgsea_version = fgsea_version),
    deviations = deviations,
    reference = list(repo = .MMC_GSEA_REF_REPO, commit = .MMC_GSEA_REF_COMMIT,
                     verified = paste0("EC scores and pathway EC sets (the ",
                                       "engine's own ECs) and the GSEA ",
                                       "calculation from them")),
    metric_label = metric_label, nes_note = mmc_gsea_nes_note(method),
    ties = NULL, n_p_unavailable = 0L, log2err = NULL, log2err_summary = NULL,
    n_ec = length(ecs$ec_exp), n_positions = length(ranked$stats),
    n_pathways = nrow(table))
}
