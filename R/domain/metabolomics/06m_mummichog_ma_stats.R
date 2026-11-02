# R/domain/metabolomics/06m_mummichog_ma_stats.R
#
# MetaboAnalyst-compatible mummichog, component 4: pathway statistics and
# permutations.
#
# Everything is counted in empirical-compound (EC) space, over ONE background
# for both polarities:
#   background  all ECs (total_matched_ecpds), whichever polarity their members
#               came from;
#   significant the ECs holding at least one feature with p <= cutoff
#               (input_ecpdlist), from one cutoff over the whole peak list;
#   null        N m/z drawn from the WHOLE peak list (both polarities, matched
#               or not), N = the number of significant features.
#
# ===========================================================================
# BEHAVIOURAL REFERENCE — pinned upstream version
# ===========================================================================
#     repo:   xia-lab/MetaboAnalystR
#     commit: 398476ae2a0c996e390925fa62adc46b5ecde334
#     file:   R/peaks_to_function.R
#
#   mmc_ma_permutations()   .perform.mummichogRTPermutations()  L2884-L2926
#                           ComputeMummichogRTPermPvals()       L2930-L2959
#                           make_ecpdlist()                     L4526
#   mmc_ma_pathway_stats()  .compute.mummichogRTSigPvals()      L3138-L3316
#                           .save.mummichog.restable()          L3339-L3355
#
# Reproduced exactly, including what upstream's own comments misdescribe:
#
#   * Per pathway, with u = significant ECs in the pathway, s = ECs in the
#     pathway, q = significant ECs, T = all ECs, nn = max(0, T + u - s - q):
#         FET  = phyper(u - 1, s, q + nn - u, q, lower.tail = FALSE)
#         EASE = fisher.test(matrix(c(max(0, u - 1), s - u + 1, q - u, nn), 2),
#                            alternative = "greater")$p.value
#   * Permutation i: set.seed(123 + i), then sample(peak list, N); the
#     permutation p-value is phyper(u - 1, s, q + nn, q, lower.tail = FALSE)
#     with nn NOT clamped at 0 and without the "- u" of the observed test.
#   * Gamma: fitdistrplus::fitdist() on |0.9999999999 - p| of ALL permutation
#     p-values pooled (all pathways x all permutations), with parameters named
#     (scale, shape); the adjusted p-value is
#         1 - pgamma(1 - FET, shape = <shape>, rate = <scale estimate>)
#     i.e. the SCALE estimate is passed as the RATE. It is applied to FET (the
#     code comment says EASE). If the fit errors, Gamma is NA.
#     The null is dominated by near-zero values (pathways a permutation does
#     not touch have p = 1), so fitdistrplus often warns while inverting the
#     Hessian for the standard errors (cov2cor(), sqrt(diag(varcovar)): an SE
#     reported as NaN). Those warnings concern the SEs only; the estimates come
#     from optim() and are all the Gamma p-values use. They are not muffled.
#   * Empirical p: sum(perm p <= FET) / n_perm, per pathway.
#   * Pathways with no significant EC are dropped; with one pathway or none
#     left upstream stops.
#   * Row order: by Gamma, else by FET when every Gamma is NA. At that point
#     upstream's table is a CHARACTER matrix (the EC-hit strings were bound
#     into it), so the p-values are ordered as TEXT, not as numbers
#     ("1e-04" sorts after "0.5").
#   * Every numeric column is then rounded with signif(, 5).
#   * Adjusted p-values: p.adjust() with no method — HOLM — on the rounded
#     FET, EASE and Gamma columns (the result table's AdjP.* columns).
#
# Intentional differences:
#   * the permutation seeds are upstream's (123 + i), but drawn through
#     withr::with_seed(), so the caller's RNG state is left untouched;
#   * a missing fitdistrplus stops the run (upstream's tryCatch would quietly
#     report Gamma as NA); an error inside the fit still gives NA, as upstream,
#     with a warning instead of upstream's print();
#   * upstream's JSON, plot and console output are not produced.
#
# Dependency: fitdistrplus (MetaboAnalystR "Suggests" it). Not in renv.lock yet.


#' One permutation's pathway p-values (MetaboAnalystR's ComputeMummichogRTPermPvals)
#'
#' @param input_ecpdlist ECs hit by the permuted m/z draw.
#' @param total_matched_ecpds All ECs (the background).
#' @param pathways Pathway sets from `mmc_ma_pathway_ecs()`.
#' @return `list(<n_pathway x 1 matrix of p-values>, <n_pathway x 1 matrix of
#'   EC overlaps>)`, as upstream returns it.
#' @noRd
.mmc_ma_perm_pvals <- function(input_ecpdlist, total_matched_ecpds, pathways) {
  ora_vec <- input_ecpdlist
  query_set_size <- length(ora_vec)
  current_mset <- pathways$emp_cpds
  total_ecpds <- unique(total_matched_ecpds)
  total_feature_num <- length(total_ecpds)

  ecpds <- lapply(current_mset, function(x) intersect(x, total_ecpds))
  feats <- lapply(current_mset, function(x) intersect(x, ora_vec))
  feat_len <- unlist(lapply(feats, length))
  set_num <- unlist(lapply(ecpds, length))

  negneg <- sizes <- vector(mode = "list", length = length(current_mset))
  for (i in seq_along(current_mset)) {
    sizes[[i]] <- feat_len[i]
    negneg[[i]] <- total_feature_num + sizes[[i]] - set_num[i] - query_set_size
  }
  unsize <- as.integer(unlist(sizes))
  res_mat <- matrix(0, nrow = length(current_mset), ncol = 1)
  fishermatrix <- cbind(unsize - 1, set_num, (query_set_size + unlist(negneg)),
                        query_set_size)
  res_mat[, 1] <- apply(fishermatrix, 1, function(x)
    stats::phyper(x[1], x[2], x[3], x[4], lower.tail = FALSE))
  list(res_mat, as.matrix(unsize))
}

#' Permutation null for the pathway statistics (MetaboAnalystR)
#'
#' Port of `.perform.mummichogRTPermutations()`: permutation `i` draws N m/z
#' from the whole peak list (both polarities, matched or not) under seed
#' `base_seed + i`, maps them to ECs, and records every pathway's permutation
#' p-value and EC overlap.
#'
#' @param prep Peak list from `mmc_ma_prepare_peaks()`.
#' @param ecs ECs from `mmc_ma_build_ecs()`.
#' @param pathways Pathway sets from `mmc_ma_pathway_ecs()`.
#' @param input_mzlist Significant m/z keys (`mmc_ma_input_mzlist()`); only its
#'   length N is used.
#' @param n_perm Number of permutations (upstream `PerformPSEA(permNum = 100)`).
#' @param base_seed Upstream seeds permutation `i` with `123 + i`; any other
#'   value is a deviation from the reference.
#' @return `list(perm_record, perm_hits)`, each a list with one element per
#'   permutation shaped as upstream stores it.
#' @examples
#' \dontrun{
#' perms <- mmc_ma_permutations(prep, ecs, pw, sig, n_perm = 100)
#' }
mmc_ma_permutations <- function(prep, ecs, pathways, input_mzlist,
                                n_perm = 100, base_seed = 123) {
  if (!is.numeric(n_perm) || length(n_perm) != 1L || is.na(n_perm) ||
      n_perm < 1 || n_perm != round(n_perm)) {
    .mmc_stop("n_perm must be a whole number >= 1; got '",
              paste(format(n_perm), collapse = ", "), "'.")
  }
  base_seed <- .mmc_as_seed(base_seed, "base_seed")
  ref_mzlist <- prep$ref_mzlist
  n_sig <- length(input_mzlist)
  record <- hits <- vector("list", n_perm)
  for (i in seq_len(n_perm)) {
    draw <- withr::with_seed(base_seed + i, sample(ref_mzlist, n_sig))
    perm_ecs <- unique(unlist(ecs$mz2ec_dict[draw]))
    perm_ecs <- perm_ecs[!is.null(perm_ecs)]
    perm <- .mmc_ma_perm_pvals(perm_ecs, ecs$total_matched_ecpds, pathways)
    record[[i]] <- perm[1]
    hits[[i]] <- perm[2]
  }
  list(perm_record = record, perm_hits = hits)
}

#' Pathway statistics in EC space (MetaboAnalystR's mummichog ORA)
#'
#' Port of `.compute.mummichogRTSigPvals()` and `.save.mummichog.restable()`:
#' per-pathway FET and EASE, the Gamma-adjusted p-value fitted on the
#' permutation null, the empirical p-value, row ordering, rounding and Holm
#' adjustment, with every quirk listed in the file header.
#'
#' @param ecs ECs from `mmc_ma_build_ecs()`.
#' @param pathways Pathway sets from `mmc_ma_pathway_ecs()`.
#' @param perms Permutations from `mmc_ma_permutations()`.
#' @return A list:
#'   \describe{
#'     \item{pvals}{Unrounded 8-column matrix for every pathway (upstream's
#'       `mSetObj$pvals`): Pathway total, Hits.total, Hits.sig, Expected, FET,
#'       EASE, Total.EC, Sig.EC.}
#'     \item{gamma}{`list(fitted, estimate, error)`: the fitdist estimates
#'       (named scale, shape) or the error that left Gamma NA.}
#'     \item{resmat}{Pathways with a significant EC, ordered and rounded:
#'       the 8 columns + Gamma, Emp.Hits, Empirical (numeric matrix).}
#'     \item{ec_hits}{Significant ECs per row of `resmat`, ";"-joined.}
#'     \item{restable}{Upstream's result table as written to its CSV (a
#'       character matrix with the Holm AdjP.* columns and pathway numbers).}
#'     \item{table}{`resmat` as a data frame with `pathway`, the Holm
#'       `adj_fet`, `adj_ease`, `adj_gamma` and `ec_hits`.}
#'   }
#' @examples
#' \dontrun{
#' res <- mmc_ma_pathway_stats(ecs, pw, perms)
#' res$table
#' }
mmc_ma_pathway_stats <- function(ecs, pathways, perms) {
  if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
    .mmc_stop("the MetaboAnalyst-compatible engine needs the 'fitdistrplus' ",
              "package for its Gamma adjustment (MetaboAnalystR fits it with ",
              "fitdistrplus::fitdist). Without it MetaboAnalyst silently ",
              "reports every Gamma p-value as NA; this engine stops instead.")
  }
  qset <- unique(unlist(ecs$input_ecpdlist))
  query_set_size <- length(qset)
  input_cpd <- unique(unlist(ecs$ecpd_cpd_dict[qset]))
  total_ecpds <- unique(ecs$total_matched_ecpds)
  total_feature_num <- length(total_ecpds)
  total_cpds <- unique(unlist(ecs$ecpd_cpd_dict))

  current_mset <- pathways$emp_cpds
  path_num <- unlist(lapply(current_mset, length))
  ecpds <- lapply(current_mset, function(x) intersect(x, total_ecpds))
  set_num <- unlist(lapply(ecpds, length))
  cpd_num <- unlist(lapply(lapply(pathways$cpds, function(x)
    intersect(x, total_cpds)), length))
  feats <- lapply(current_mset, function(x) intersect(x, qset))
  feat_len <- unlist(lapply(feats, length))
  feat_vec <- sapply(feats, function(x) paste(x, collapse = ";"))
  cpd_len <- lapply(pathways$cpds, function(x) intersect(x, input_cpd))

  negneg <- sizes <- vector(mode = "list", length = length(current_mset))
  for (i in seq_along(current_mset)) {
    sizes[[i]] <- feat_len[i]
    negneg[[i]] <- total_feature_num + sizes[[i]] - set_num[i] - query_set_size
  }
  negneg <- rapply(negneg, function(x) ifelse(x < 0, 0, x), how = "replace")
  unsize <- as.integer(unlist(sizes))
  uniq_count <- length(unique(unlist(current_mset, use.names = FALSE)))

  res_mat <- matrix(0, nrow = length(current_mset), ncol = 8)
  fishermatrix <- cbind(unsize - 1, set_num,
                        (query_set_size + unlist(negneg) - unsize), query_set_size)
  first <- unlist(lapply(sizes, function(x) max(0, x - 1)))
  easematrix <- cbind(first, (set_num - unsize + 1), (query_set_size - unsize),
                      unlist(negneg))
  res_mat[, 1] <- unlist(lapply(pathways$cpds, length))
  res_mat[, 2] <- cpd_num
  res_mat[, 3] <- as.integer(unlist(lapply(cpd_len, length)))
  res_mat[, 4] <- query_set_size * (path_num / uniq_count)
  res_mat[, 5] <- apply(fishermatrix, 1, function(x)
    stats::phyper(x[1], x[2], x[3], x[4], lower.tail = FALSE))
  res_mat[, 6] <- apply(easematrix, 1, function(x)
    stats::fisher.test(matrix(x, nrow = 2), alternative = "greater")$p.value)
  res_mat[, 7] <- set_num
  res_mat[, 8] <- unsize
  colnames(res_mat) <- c("Pathway total", "Hits.total", "Hits.sig", "Expected",
                         "FET", "EASE", "Total.EC", "Sig.EC")
  rownames(res_mat) <- pathways$name
  pvals <- res_mat

  # Gamma on the pooled permutation null.
  sig_hits <- unsize
  sigpvalue <- res_mat[, 5]
  perm_minus <- abs(0.9999999999 - unlist(perms$perm_record))
  gamma <- list(fitted = FALSE, estimate = NULL, error = NULL)
  adjustedp <- rep(NA, length = length(res_mat[, 1]))
  if (!(length(sig_hits[sig_hits != 0]) < round(length(sig_hits) * 0.05))) {
    fit <- tryCatch(
      fitdistrplus::fitdist(perm_minus, distr = "gamma", method = "mle",
                            lower = c(0, 0), start = list(scale = 1, shape = 1)),
      error = function(e) e)
    if (inherits(fit, "error")) {
      gamma$error <- conditionMessage(fit)
      warning("[mummichog] MetaboAnalyst-compatible Gamma fit failed, so every ",
              "Gamma p-value is NA, as in MetaboAnalyst: ", gamma$error,
              call. = FALSE)
    } else {
      gamma$fitted <- TRUE
      gamma$estimate <- fit$estimate
      rawpval <- as.numeric(sigpvalue)
      adjustedp <- 1 - (stats::pgamma(1 - rawpval,
                                      shape = fit$estimate["shape"],
                                      rate = fit$estimate["scale"]))
    }
  }
  res_mat <- cbind(res_mat, Gamma = adjustedp)

  # Empirical p-values: pathways in rows, permutations in columns.
  fisher_p <- as.numeric(res_mat[, 5])
  record_matrix <- do.call(cbind, do.call(cbind, perms$perm_record))
  num_perm <- ncol(record_matrix)
  better_hits <- sapply(seq_along(record_matrix[, 1]), function(i)
    sum(record_matrix[i, ] <= fisher_p[i]))
  emp_p <- sapply(seq_along(record_matrix[, 1]), function(i)
    (sum(record_matrix[i, ] <= fisher_p[i]) / num_perm))
  res_mat <- cbind(res_mat, Emp.Hits = better_hits, Empirical = emp_p,
                   EC.Hits = feat_vec)

  hit_inx <- as.numeric(as.character(res_mat[, 8])) > 0
  res_mat <- res_mat[hit_inx, , drop = FALSE]
  if (nrow(res_mat) <= 1) {
    .mmc_stop("only ", nrow(res_mat), " pathway(s) contain a significant ",
              "empirical compound; MetaboAnalyst stops here as too little ",
              "m/z-to-pathway signal for pathway analysis.")
  }
  path_nms <- pathways$name[hit_inx]

  # res_mat is a character matrix here, so Gamma is ordered as text (upstream).
  if (length(stats::na.omit(res_mat[, 9])) == 0) {
    ord_inx <- order(res_mat[, 5])
  } else {
    ord_inx <- order(res_mat[, 9])
  }
  ec_hits <- res_mat[ord_inx, 12]
  resmat <- signif(apply(as.matrix(res_mat[ord_inx, 1:11]), 2, as.numeric), 5)
  rownames(resmat) <- path_nms[ord_inx]

  # .save.mummichog.restable(): Holm (p.adjust's default) on rounded columns.
  adj_fet   <- stats::p.adjust(resmat[, "FET"])
  adj_ease  <- stats::p.adjust(resmat[, "EASE"])
  adj_gamma <- stats::p.adjust(resmat[, "Gamma"])
  restable <- cbind(resmat, AdjP.Fisher = adj_fet, AdjP.EASE = adj_ease,
                    AdjP.Gamma = adj_gamma)
  restable <- restable[, -c(7, 8)]
  restable <- cbind(restable, paste0("P", seq.int(1, nrow(restable))))
  colnames(restable)[ncol(restable)] <- "Pathway Number"
  colnames(restable)[which(colnames(restable) == "FET")] <- "P(Fisher)"
  colnames(restable)[which(colnames(restable) == "EASE")] <- "P(EASE)"
  colnames(restable)[which(colnames(restable) == "Gamma")] <- "P(Gamma)"
  # upstream's argument is named cpd.hits, which names the column
  restable <- cbind(restable, cpd.hits = ec_hits)

  table <- data.frame(pathway = rownames(resmat), resmat, adj_fet = adj_fet,
                      adj_ease = adj_ease, adj_gamma = adj_gamma,
                      ec_hits = unname(ec_hits), check.names = FALSE,
                      stringsAsFactors = FALSE, row.names = NULL)

  list(pvals = pvals, gamma = gamma, resmat = resmat, ec_hits = ec_hits,
       restable = restable, table = table)
}
