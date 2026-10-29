# R/domain/metabolomics/06h_mummichog_summary.R
#
# Interpretation layer for the mummichog report section: one compact Summary
# table per contrast, built from results the pipeline already holds -- the
# mummichog ORA table (06e), the mummichog GSEA table (06g), the pathway
# evidence (06f) and the identified-metabolite ORA / QEA / GSEA tables (06).
# Nothing here computes a statistic. It selects rows, joins them on the shared
# pathway id and labels them, so the report describes these results and does not
# re-decide them.
#
# The cross-analysis join is the exact pathway id both sides share (e.g.
# `mng00010`): the model's own pathway id on the mummichog side, the GMT id on
# the identified side. Never a pathway name, never a normalised id.


# Adjusted-significance threshold this layer reads: a pathway counts as
# supported by an identified-metabolite method, or as surviving BH correction in
# mummichog GSEA, at FDR below this. It mirrors the literal 0.05 the identified
# enrichment functions themselves report against (06_enrichment.R); there is no
# config key for it.
.MMC_SUMMARY_FDR <- 0.05

# Identified-metabolite methods that count as support, in display order, and the
# enrichment_res slot holding each one's per-contrast results. ssGSEA is left
# out on purpose.
.MMC_IDENTIFIED_METHODS <- c(ORA  = "ora_by_contrast",
                             QEA  = "qea_by_contrast",
                             GSEA = "gsea_by_contrast")

# How many pathways each contrast's Summary lists (and exports).
.MMC_SUMMARY_TOP_N <- 20L

# The Summary columns the report shows; the export carries the rest as well.
.MMC_SUMMARY_DISPLAY_COLS <- c(
  "Pathway", "Identified-metabolite support", "Mummichog ORA",
  "Mummichog GSEA", "Mummichog contribution", "Supporting evidence",
  "Annotation agreement")


#' Recover the pathway id from an identified-enrichment pathway label
#'
#' The identified ORA / QEA / GSEA tables label pathways with
#' `make_pathway_labels()`, which writes `"<id> - <description>"`, or the bare id
#' when the GMT has no description. The id is therefore exactly the text before
#' the first `" - "`.
#'
#' @param label Character vector of pathway labels.
#' @return Character vector of ids, same length as `label`.
#' @noRd
.mmc_identified_pathway_id <- function(label) {
  sub(" - .*$", "", as.character(label))
}

#' Pathway ids in the identified-metabolite enrichment library
#'
#' Reads the GMT file(s) the identified enrichment is configured with, resolved
#' the way the enrichment functions resolve them, and returns their pathway ids.
#' This is the shared namespace a mummichog pathway must belong to before it can
#' be compared with identified-metabolite enrichment at all.
#'
#' @param config Full pipeline config.
#' @return Character vector of pathway ids; empty when no GMT is configured or
#'   readable. A GMT that exists but cannot be read is skipped with a warning.
mmc_identified_pathway_library <- function(config) {
  gmt <- tryCatch(
    resolve_input_path(config, config$modes$metabolomics$enrichment$gmt_file),
    error = function(e) NULL
  )
  if (length(gmt) == 0) return(character(0))
  gmt <- gmt[!is.na(gmt) & nzchar(gmt) & file.exists(gmt)]
  # A GMT that exists but cannot be read must not take the whole report down:
  # identified-metabolite support is optional, so that library is skipped.
  ids <- unlist(lapply(gmt, function(f) tryCatch(
    names(read_gmt_list(f)),
    error = function(e) {
      warning("mummichog summary: could not read GMT '", f, "' (",
              conditionMessage(e), ") — its pathways are treated as not ",
              "comparable", call. = FALSE)
      character(0)
    })), use.names = FALSE)
  unique(as.character(ids))
}

#' Which pathways each identified-metabolite method tested, and which it supports
#'
#' Reads one contrast's identified ORA, QEA and GSEA result tables as they are
#' and flags each row's `FDR` against `fdr`. Methods that did not run for the
#' contrast simply contribute no rows.
#'
#' @param enrichment_res Result of `mod_metabolomics_enrichment()`, or NULL.
#' @param contrast       Contrast label (a key of its `*_by_contrast` lists).
#' @param fdr            Adjusted-p threshold for support.
#' @return A data.frame with `pathway_id`, `method` and `supported`, one row per
#'   tested (method, pathway); zero rows when nothing was tested.
mmc_identified_support <- function(enrichment_res, contrast,
                                   fdr = .MMC_SUMMARY_FDR) {
  rows <- lapply(names(.MMC_IDENTIFIED_METHODS), function(m) {
    tbl <- enrichment_res[[.MMC_IDENTIFIED_METHODS[[m]]]][[contrast]]$res$table
    if (!is.data.frame(tbl) || nrow(tbl) == 0 ||
        !all(c("pathway", "FDR") %in% names(tbl))) {
      return(NULL)
    }
    data.frame(pathway_id = .mmc_identified_pathway_id(tbl$pathway),
               method     = m,
               supported  = !is.na(tbl$FDR) & tbl$FDR < fdr,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(pathway_id = character(0), method = character(0),
                      supported = logical(0), stringsAsFactors = FALSE)
  }
  out
}

#' Collapse EC-level agreement counts into one readable value
#'
#' @param pathway_summary Rows of `evidence$pathway_summary` (or any frame with
#'   the `ECs Match` / `ECs Conflict` / `ECs Mixed` / `ECs Not assessed` columns).
#' @return Character vector, one per row, e.g. `"1 Match · 2 Conflict"`; zero
#'   categories are omitted, and a row with no counts gives `NA`.
mmc_format_agreement <- function(pathway_summary) {
  if (!is.data.frame(pathway_summary) || nrow(pathway_summary) == 0) {
    return(character(0))
  }
  vapply(seq_len(nrow(pathway_summary)), function(i) {
    n <- vapply(.MMC_AGREEMENT_STATES, function(st) {
      v <- pathway_summary[[paste("ECs", st)]][i]
      if (is.null(v) || is.na(v)) 0 else as.numeric(v)
    }, numeric(1))
    if (all(n == 0)) return(NA_character_)
    paste(sprintf("%d %s", as.integer(n[n > 0]), names(n)[n > 0]),
          collapse = " · ")
  }, character(1))
}

#' One-line count of nominally significant pathways for a plot tab
#'
#' @param p        Raw pathway p-values.
#' @param p_cutoff Nominal threshold.
#' @param analysis Name of the analysis, e.g. `"Mummichog ORA"`.
#' @param padj     Optional BH-adjusted p-values aligned with `p`; when given,
#'   the sentence also says how many nominal pathways survive correction.
#' @param fdr      Adjusted-p threshold.
#' @return A single sentence.
mmc_nominal_sentence <- function(p, p_cutoff, analysis, padj = NULL,
                                 fdr = .MMC_SUMMARY_FDR) {
  nominal <- !is.na(p) & p < p_cutoff
  n <- sum(nominal)
  lead <- if (n == 0) "No pathway" else if (n == 1) "One pathway" else
    sprintf("%d pathways", n)
  s <- sprintf("%s reached nominal p < %g in %s", lead, p_cutoff, analysis)
  if (is.null(padj) || n == 0) return(paste0(s, "."))
  k <- sum(nominal & !is.na(padj) & padj < fdr)
  sprintf("%s; %s significant after BH correction (FDR < %g).", s,
          if (k == 0) "none remained" else
            sprintf("%d remained", k), fdr)
}

#' Build one contrast's mummichog Summary table
#'
#' Ranks every pathway by the smaller of its raw mummichog ORA `p.value` and
#' GSEA `P.Value` (a missing value counts as `Inf`), keeps the first `top_n`,
#' and places beside each one what the other results already say about it. A
#' row need not be significant: `nominal` (ORA or GSEA below `p_cutoff`) says
#' whether it carries a pathway-level mummichog signal at all.
#'
#' Identified-metabolite support and the mummichog contribution follow three
#' separate questions:
#' \describe{
#'   \item{joinable}{the pathway has a single model id that is also a pathway of
#'     the identified-enrichment library.}
#'   \item{tested}{that id appears in at least one of this contrast's identified
#'     ORA / QEA / GSEA result tables.}
#'   \item{supported}{at least one of those tested methods has FDR < `fdr`.}
#' }
#' \preformatted{
#' Identified-metabolite support:  !joinable -> "Not comparable"
#'                                  !tested   -> "Not tested"
#'                                  supported -> the methods, e.g. "ORA + QEA"
#'                                  otherwise -> "None"
#' Mummichog contribution:          !nominal  -> "No nominal Mummichog signal"
#'                                  !joinable -> "Not comparable"
#'                                  supported -> "Confirms identified signal"
#'                                  otherwise -> "Mummichog-only candidate"
#' }
#'
#' @param ora_table   The section's ORA table (`build_mummichog_pathway_table()`),
#'   or NULL.
#' @param gsea_table  The GSEA result table (`run_mummichog_gsea()$table`), or
#'   NULL when GSEA did not run.
#' @param evidence    The evidence list (`build_mummichog_pathway_evidence()`),
#'   or NULL.
#' @param pathway_ids Model pathway name -> id (`model$pathway_ids`), or NULL.
#' @param identified  Result of `mmc_identified_support()`, or NULL.
#' @param library_ids Result of `mmc_identified_pathway_library()`.
#' @param p_cutoff    Nominal threshold (`mummichog.p_cutoff`).
#' @param fdr         Adjusted-p threshold.
#' @param top_n       Number of top-ranked pathways to keep.
#' @return A data.frame, one row per kept pathway in rank order: the
#'   `.MMC_SUMMARY_DISPLAY_COLS` labels followed by the values behind them
#'   (`Pathway ID`, `Nominal Mummichog signal`, `ORA p.value`, `GSEA P.Value`,
#'   `GSEA padj`, `Supporting ECs`, `Supporting features`,
#'   `Identified methods tested`). The evidence cells, agreement and counts are
#'   filled only for pathways nominal in ORA; a GSEA-only pathway reads
#'   `Not available — GSEA-only signal` even when an ORA evidence row exists for
#'   it, and a pathway nominal in neither reads `—`. Zero rows only when there
#'   are no pathway results at all.
build_mummichog_summary <- function(ora_table, gsea_table = NULL,
                                    evidence = NULL, pathway_ids = NULL,
                                    identified = NULL,
                                    library_ids = character(0),
                                    p_cutoff = 0.05, fdr = .MMC_SUMMARY_FDR,
                                    top_n = .MMC_SUMMARY_TOP_N) {
  lookup <- function(tbl, key, col) {
    if (!is.data.frame(tbl) || !col %in% names(tbl)) {
      return(rep(NA_real_, length(key)))
    }
    suppressWarnings(as.numeric(tbl[[col]][match(key, tbl$Pathway)]))
  }
  all_pw <- as.character(unique(c(
    if (is.data.frame(ora_table)) ora_table$Pathway,
    if (is.data.frame(gsea_table)) gsea_table$Pathway)))
  ora_p  <- lookup(ora_table, all_pw, "p.value")
  gsea_p <- lookup(gsea_table, all_pw, "P.Value")
  rank_p <- pmin(ifelse(is.na(ora_p), Inf, ora_p),
                 ifelse(is.na(gsea_p), Inf, gsea_p))
  keep <- utils::head(order(rank_p, all_pw), top_n)

  pw       <- all_pw[keep]
  ora_p    <- ora_p[keep]
  gsea_p   <- gsea_p[keep]
  gsea_q   <- lookup(gsea_table, pw, "padj")
  ora_nom  <- !is.na(ora_p) & ora_p < p_cutoff
  gsea_nom <- !is.na(gsea_p) & gsea_p < p_cutoff
  nominal  <- ora_nom | gsea_nom
  n <- length(pw)

  # -- the cross-analysis join, on the exact shared pathway id ---------------
  pid <- if (is.null(pathway_ids)) rep(NA_character_, n) else
    unname(pathway_ids[pw])
  joinable <- !is.na(pid) & pid %in% library_ids
  if (is.null(identified)) identified <- mmc_identified_support(NULL, "")
  methods <- lapply(pid, function(id) identified[identified$pathway_id %in% id, ,
                                                 drop = FALSE])
  tested    <- joinable & vapply(methods, nrow, integer(1)) > 0
  sig_label <- vapply(methods, function(m) {
    hit <- intersect(names(.MMC_IDENTIFIED_METHODS), m$method[m$supported])
    paste(hit, collapse = " + ")
  }, character(1))
  supported <- tested & nzchar(sig_label)

  ident <- rep("None", n)
  ident[supported]  <- sig_label[supported]
  ident[!tested]    <- "Not tested"
  ident[!joinable]  <- "Not comparable"
  # A row without a nominal mummichog signal is ranked, not found: it is never
  # called a candidate, whatever the identified side says.
  contrib <- rep("Mummichog-only candidate", n)
  contrib[supported] <- "Confirms identified signal"
  contrib[!joinable] <- "Not comparable"
  contrib[!nominal]  <- "No nominal Mummichog signal"

  # -- evidence (ORA-derived) -------------------------------------------------
  ps <- evidence$pathway_summary
  ev_rows <- if (is.data.frame(ps)) ps[match(pw, ps$Pathway), , drop = FALSE] else
    NULL
  # The evidence cell is support for the ORA signal, so it is shown only for
  # pathways nominal in ORA. A GSEA-only pathway can still have an evidence row
  # (its ORA overlap exists even when ORA p is not nominal); showing it here
  # would read as ORA support it does not have. It stays in the full evidence
  # tables under Advanced.
  has_ev  <- if (is.null(ev_rows)) rep(FALSE, n) else
    !is.na(ev_rows$Pathway) & ora_nom
  n_ec    <- if (is.null(ev_rows)) rep(NA_real_, n) else ev_rows[["Supporting ECs"]]
  n_feat  <- if (is.null(ev_rows)) rep(NA_real_, n) else ev_rows[["Supporting features"]]
  n_ec[!has_ev]   <- NA_real_
  n_feat[!has_ev] <- NA_real_
  plural  <- function(k, word) sprintf("%d %s%s", as.integer(k), word,
                                       ifelse(k == 1, "", "s"))
  ev_label <- rep("Not available", n)
  ev_label[!ora_nom & gsea_nom] <- "Not available — GSEA-only signal"
  ev_label[!nominal] <- "—"
  agree <- rep(NA_character_, n)
  if (any(has_ev)) {
    ev_label[has_ev] <- paste(plural(n_ec[has_ev], "EC"), "/",
                              plural(n_feat[has_ev], "feature"))
    agree[has_ev] <- mmc_format_agreement(ev_rows[has_ev, , drop = FALSE])
  }
  agree[is.na(agree)] <- "—"

  fmt <- function(x) formatC(x, digits = 2, format = "g")
  ora_label <- ifelse(is.na(ora_p), "Not tested", paste("p =", fmt(ora_p)))
  gsea_label <- if (!is.data.frame(gsea_table)) rep("Not run", n) else
    ifelse(is.na(gsea_p), "Not tested",
           sprintf("p = %s; FDR = %s", fmt(gsea_p), fmt(gsea_q)))

  data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    "Pathway"                       = pw,
    "Identified-metabolite support" = ident,
    "Mummichog ORA"                 = as.character(ora_label),
    "Mummichog GSEA"                = as.character(gsea_label),
    "Mummichog contribution"        = contrib,
    "Supporting evidence"           = ev_label,
    "Annotation agreement"          = agree,
    "Pathway ID"                    = pid,
    "Nominal Mummichog signal"      = nominal,
    "ORA p.value"                   = ora_p,
    "GSEA P.Value"                  = gsea_p,
    "GSEA padj"                     = gsea_q,
    "Supporting ECs"                = as.numeric(n_ec),
    "Supporting features"           = as.numeric(n_feat),
    "Identified methods tested"     = vapply(seq_len(n), function(i) {
      if (!joinable[i]) return(NA_character_)
      m <- intersect(names(.MMC_IDENTIFIED_METHODS), methods[[i]]$method)
      if (length(m) == 0) NA_character_ else paste(m, collapse = " + ")
    }, character(1))
  )
}
