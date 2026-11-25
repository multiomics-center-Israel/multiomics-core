# R/domain/metabolomics/06u_mummichog_across.R
#
# Across-contrast view of the mummichog results: which pathways show evidence
# in more than one contrast.
#
# It is a lookup across results that already exist, not an analysis. Each
# contrast's evidence comes from its overview (06t), which the report and the
# workbook also show, with the evidence level decided by the one shared rule
# (`.mmc_evidence_level()`): ORA by its raw and Holm-adjusted p, GSEA by its raw
# p and BH FDR. Pathways are aligned on the model's pathway id (the name only
# where a pathway has no id). Nothing is combined across contrasts: no
# cross-contrast test, no pooled p-value, no direction. A pathway that recurs
# can do so through different empirical compounds and features, and with
# opposite or no consistent direction; the per-contrast evidence tables keep
# those differences visible.


#' One contrast's evidence for one pathway, as short text
#'
#' @param ora_lvl,gsea_lvl Levels from `.mmc_evidence_level()`.
#' @param gsea_p_unavailable Logical: GSEA tested the pathway but could give no
#'   p-value. Its level is `"absent"` (there is no p to grade), yet it has a
#'   result, so it must not read `"No result"`.
#' @return Character vector, e.g. `"ORA Holm · GSEA nominal"`, `"—"` (results
#'   but no nominal evidence) or `"No result"` (in neither test).
#' @noRd
.mmc_across_cell <- function(ora_lvl, gsea_lvl,
                             gsea_p_unavailable = rep(FALSE, length(gsea_lvl))) {
  ora  <- c(adjusted = "ORA Holm", nominal = "ORA nominal", none = "", absent = "")[ora_lvl]
  gsea <- c(adjusted = "GSEA FDR", nominal = "GSEA nominal", none = "", absent = "")[gsea_lvl]
  out <- ifelse(nzchar(ora) & nzchar(gsea), paste(ora, "·", gsea),
                ifelse(nzchar(ora), ora, gsea))
  out[!nzchar(out)] <- "—"
  out[ora_lvl == "absent" & gsea_lvl == "absent" & !gsea_p_unavailable] <- "No result"
  unname(out)
}

#' Which pathways the report's Across contrasts view lists
#'
#' The report lists a pathway when it has nominal evidence (ORA or GSEA raw p
#' below the cutoff) in at least two contrasts, OR adjusted significance (ORA
#' Holm or GSEA BH FDR below the threshold) in at least one. The workbook's
#' Across Contrasts sheet lists every pathway; this selects only the report's
#' rows.
#'
#' @param wide The `wide` table of `build_mummichog_across_contrasts()`.
#' @return Logical vector, one per row of `wide`.
#' @examples
#' mmc_across_report_rows(data.frame(
#'   "Contrasts with nominal evidence" = c(2, 1, 1),
#'   "Contrasts with adjusted evidence" = c(0, 1, 0), check.names = FALSE))
mmc_across_report_rows <- function(wide) {
  nom <- wide[["Contrasts with nominal evidence"]]
  adj <- wide[["Contrasts with adjusted evidence"]]
  (!is.na(nom) & nom >= 2) | (!is.na(adj) & adj >= 1)
}

#' Pathways across contrasts, from each contrast's overview
#'
#' Recurrence is counted over DE contrasts, not over report sections: the
#' python engine's mixed mode gives one section per polarity of a contrast
#' (labels resolved by `mmc_run_label_parts()` into each section's
#' `de_contrast`), and two polarity runs of one contrast are not two contrasts.
#' The per-section evidence cells are kept as they are.
#'
#' @param sections Report sections (`build_mummichog_report_sections()`), each
#'   carrying `overview`, `p_cutoff` and `de_contrast`. Sections without an
#'   overview are skipped.
#' @param fdr Adjusted-p threshold.
#' @return NULL when no section has an overview; otherwise a list:
#'   \describe{
#'     \item{long}{One row per (contrast, pathway): the raw and adjusted
#'       p-values of both tests as the overview holds them, each test's level
#'       and the short evidence text.}
#'     \item{wide}{One row per pathway: the number of contrasts with adjusted
#'       evidence (ORA Holm or GSEA BH FDR below `fdr`), with nominal evidence
#'       (raw p below the cutoff in either test; includes the adjusted ones)
#'       and with any result; then the evidence text per contrast; then each
#'       contrast's four p-values. Ordered by adjusted count, nominal count,
#'       smallest raw p, name.}
#'     \item{contrasts}{The section labels, in section order (one column of
#'       evidence each).}
#'     \item{de_contrasts}{The distinct DE contrasts behind them.}
#'   }
#' @examples
#' build_mummichog_across_contrasts(list())
build_mummichog_across_contrasts <- function(sections, fdr = .MMC_SUMMARY_FDR) {
  have <- Filter(function(s) is.data.frame(s$overview), sections)
  if (length(have) == 0L) return(NULL)
  contrasts <- names(have)

  long <- do.call(rbind, lapply(contrasts, function(ct) {
    s  <- have[[ct]]
    ov <- s$overview
    if (nrow(ov) == 0L) return(NULL)
    p_cut <- s$p_cutoff %||% 0.05
    ora_lvl  <- .mmc_evidence_level(ov[["ORA raw p"]], ov[["ORA Holm-adjusted p"]], p_cut, fdr)
    gsea_lvl <- .mmc_evidence_level(ov[["GSEA raw p"]], ov[["GSEA BH FDR"]], p_cut, fdr)
    id <- as.character(ov[["Pathway ID"]])
    data.frame(
      check.names = FALSE, stringsAsFactors = FALSE,
      "Contrast"            = ct,
      "DE contrast"         = s$de_contrast %||% ct,
      "Pathway key"         = ifelse(is.na(id) | !nzchar(id),
                                     paste0("name:", ov$Pathway), id),
      "Pathway"             = ov$Pathway,
      "Pathway ID"          = id,
      "ORA raw p"           = ov[["ORA raw p"]],
      "ORA Holm-adjusted p" = ov[["ORA Holm-adjusted p"]],
      "GSEA raw p"          = ov[["GSEA raw p"]],
      "GSEA BH FDR"         = ov[["GSEA BH FDR"]],
      "ORA level"           = ora_lvl,
      "GSEA level"          = gsea_lvl,
      "Evidence"            = .mmc_across_cell(
        ora_lvl, gsea_lvl,
        gsea_p_unavailable = ov[["GSEA status"]] %in% .MMC_GSEA_P_UNAVAILABLE))
  }))
  if (is.null(long) || nrow(long) == 0L) return(NULL)

  keys <- unique(long[["Pathway key"]])
  by_key <- split(long, factor(long[["Pathway key"]], levels = keys))
  rows <- lapply(keys, function(k) {
    d <- by_key[[k]]
    adj <- d[["ORA level"]] == "adjusted" | d[["GSEA level"]] == "adjusted"
    nom <- adj | d[["ORA level"]] == "nominal" | d[["GSEA level"]] == "nominal"
    best <- suppressWarnings(min(c(d[["ORA raw p"]], d[["GSEA raw p"]]), na.rm = TRUE))
    out <- data.frame(
      check.names = FALSE, stringsAsFactors = FALSE,
      "Pathway"    = d$Pathway[[1]],
      "Pathway ID" = d[["Pathway ID"]][[1]],
      "Contrasts with adjusted evidence" = length(unique(d[["DE contrast"]][adj])),
      "Contrasts with nominal evidence"  = length(unique(d[["DE contrast"]][nom])),
      "Contrasts with a result"          = length(unique(d[["DE contrast"]])),
      ".best_p" = if (is.finite(best)) best else NA_real_)
    for (ct in contrasts) {
      r <- d[d$Contrast == ct, , drop = FALSE]
      out[[ct]] <- if (nrow(r) == 0L) "No result" else r$Evidence[[1]]
    }
    for (ct in contrasts) {
      r <- d[d$Contrast == ct, , drop = FALSE]
      pick <- function(col) if (nrow(r) == 0L) NA_real_ else r[[col]][[1]]
      out[[paste0(ct, ": ORA raw p")]]           <- pick("ORA raw p")
      out[[paste0(ct, ": ORA Holm-adjusted p")]] <- pick("ORA Holm-adjusted p")
      out[[paste0(ct, ": GSEA raw p")]]          <- pick("GSEA raw p")
      out[[paste0(ct, ": GSEA BH FDR")]]         <- pick("GSEA BH FDR")
    }
    out
  })
  wide <- do.call(rbind, rows)
  ord <- order(-wide[["Contrasts with adjusted evidence"]],
               -wide[["Contrasts with nominal evidence"]],
               wide[[".best_p"]], wide$Pathway, method = "radix", na.last = TRUE)
  wide <- wide[ord, setdiff(names(wide), ".best_p"), drop = FALSE]
  rownames(wide) <- NULL
  list(long = long, wide = wide, contrasts = contrasts,
       de_contrasts = unique(long[["DE contrast"]]))
}
