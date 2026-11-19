# R/domain/metabolomics/06t_mummichog_presentation.R
#
# Presentation layer for the mummichog report and the mummichog workbook.
#
# The report sections (06e) already hold every result a reader needs: the ORA
# table, the GSEA table, the evidence (06f / 06s) and the Summary (06h). The
# Summary and the ORA table are also the Shiny app's contract
# (07_shiny_export.R, docs/SHINY_PAYLOAD_HANDOFF.md), so their columns and
# values are left exactly as they are. This file builds the reader-facing
# versions BESIDE them, from the same objects:
#
#   build_mummichog_overview()     every pathway of a contrast, one row each,
#                                  with each test's raw and adjusted p-value
#                                  and its evidence level spelled out;
#   build_mummichog_ora_detail()   the complete ORA table with every engine
#                                  column under a name that says what it is;
#   mmc_input_counts()             what entered the engine for a contrast, read
#                                  from its run record and the DE results;
#   .mmc_evidence_level()          the one rule that turns a raw and an
#                                  adjusted p-value into "adjusted" /
#                                  "nominal" / "none" / "absent", shared by the
#                                  overview, the across-contrast view (06u) and
#                                  the workbook (06v).
#
# Nothing here computes a statistic, a significance decision of the engine, or
# an evidence count: the overview is built ON build_mummichog_summary() (same
# ranking, same identified-metabolite join, same evidence rules) and relabels
# it, and the input counts are the engine's own recorded counts.


#' Evidence level of one test from its raw and adjusted p-values
#'
#' The single rule the report, the across-contrast view and the workbook use
#' to say how strong a pathway's evidence is in one test:
#' \describe{
#'   \item{`"adjusted"`}{nominal (below) AND adjusted p (Holm for ORA, BH FDR
#'     for GSEA) < `fdr`;}
#'   \item{`"nominal"`}{raw p < `p_cutoff` but not adjusted-significant;}
#'   \item{`"none"`}{the pathway has a result but is not nominal;}
#'   \item{`"absent"`}{the pathway has no result in this test.}
#' }
#' "Adjusted" is a subset of "nominal" on purpose: the Summary (06h) decides
#' whether a pathway carries a Mummichog signal, and whether its ORA evidence
#' is shown, by raw p < `p_cutoff` alone. With a `p_cutoff` below `fdr`, a
#' pathway could otherwise read "adjusted" here while the Summary calls it
#' "No nominal Mummichog signal". An adjusted value that is NA (e.g. the python
#' engine reports no adjusted ORA p-value) can never give `"adjusted"`.
#'
#' @param p        Raw p-values.
#' @param adj      Adjusted p-values (same length; NA allowed).
#' @param p_cutoff Nominal threshold (`mummichog.p_cutoff`).
#' @param fdr      Adjusted-p threshold.
#' @return Character vector of levels, same length as `p`.
#' @noRd
.mmc_evidence_level <- function(p, adj, p_cutoff, fdr = .MMC_SUMMARY_FDR) {
  p   <- suppressWarnings(as.numeric(p))
  adj <- suppressWarnings(as.numeric(adj))
  if (length(adj) == 0L) adj <- rep(NA_real_, length(p))
  out <- rep("none", length(p))
  nominal <- !is.na(p) & p < p_cutoff
  out[nominal] <- "nominal"
  out[nominal & !is.na(adj) & adj < fdr] <- "adjusted"
  out[is.na(p)] <- "absent"
  out
}

# Summary labels (06h) -> the reader-facing wording. Keyed by 06h's own labels
# so the meaning is decided there once; an unknown label stops rather than
# passing through unexplained.
.MMC_INTERPRETATION_LABELS <- c(
  "No nominal Mummichog signal" = "No nominal Mummichog signal",
  "Not comparable"              = "Not comparable (no shared pathway ID)",
  "Confirms identified signal"  = "Also supported by identified-metabolite enrichment",
  "Mummichog-only candidate"    = "Mummichog-only hypothesis")

#' Map 06h labels to their reader-facing wording, refusing unknown ones
#' @param x   Labels from `build_mummichog_summary()`.
#' @param map Named character vector, 06h label -> wording.
#' @return Character vector.
#' @noRd
.mmc_relabel <- function(x, map) {
  bad <- setdiff(unique(x[!is.na(x)]), names(map))
  if (length(bad) > 0L) {
    .mmc_stop("no reader-facing wording for Summary label(s) ",
              paste0("'", bad, "'", collapse = ", "), ".")
  }
  unname(map[x])
}

#' Reader-facing wording of the identified-metabolite support column
#' @param x `Identified-metabolite support` values from 06h.
#' @param fdr Adjusted-p threshold the support was judged at.
#' @return Character vector.
#' @noRd
.mmc_identified_wording <- function(x, fdr = .MMC_SUMMARY_FDR) {
  out <- x
  out[x == "Not comparable"] <- "Not comparable (no shared pathway ID)"
  out[x == "Not tested"]     <- "Not in this contrast's identified-metabolite results"
  out[x == "None"]           <- sprintf("Tested; none at FDR < %g", fdr)
  methods <- !x %in% c("Not comparable", "Not tested", "None") & !is.na(x)
  out[methods] <- sprintf("%s (FDR < %g)", x[methods], fdr)
  out
}

# Column names of the overview, fixed so the report, the workbook and the
# across-contrast view address the same columns.
.MMC_OVERVIEW_COLS <- c(
  "Pathway", "Pathway ID",
  "ORA status", "ORA raw p", "ORA Holm-adjusted p",
  "GSEA status", "GSEA raw p", "GSEA BH FDR",
  "ORA overlap ECs (Sig.EC)", "Features in overlap ECs",
  "Annotation agreement (ECs)", "ORA evidence",
  "Identified-metabolite enrichment", "Interpretation",
  "Nominal Mummichog signal", "ORA raw p type")

#' Build the reader-facing overview of one contrast's mummichog results
#'
#' Every pathway with an ORA or GSEA result, one row each, in the Summary's
#' ranking (the smaller raw p of ORA and GSEA first). It is built on
#' `build_mummichog_summary()` called for every pathway, so the first
#' `.MMC_SUMMARY_TOP_N` rows are exactly the Summary's rows, and the
#' identified-metabolite join and the evidence rules are the Summary's own.
#' What it adds is wording and the values the Summary does not carry:
#'
#' \itemize{
#'   \item each test's status (`.mmc_evidence_level()`), with ORA's raw p and
#'     its Holm-adjusted p (MetaboAnalyst-compatible engine) and GSEA's raw p
#'     and BH FDR in separate columns;
#'   \item for the MetaboAnalyst-compatible engine, a pathway absent from the
#'     ORA table had no significant EC (the engine reports only pathways with
#'     one), so it reads "No significant EC overlap", never "Not tested";
#'   \item `ORA overlap ECs (Sig.EC)` and `Features in overlap ECs` only for
#'     pathways nominal in ORA, so a GSEA-only pathway is never shown ORA
#'     evidence; the feature count includes features that are not significant
#'     and features that matched no pathway candidate themselves.
#' }
#'
#' @param ora_table   The section's ORA table (`build_mummichog_pathway_table()`).
#' @param gsea_table  The GSEA result table, or NULL when GSEA did not run.
#' @param evidence    The section's evidence, or NULL.
#' @param pathway_ids Model pathway name -> id.
#' @param identified  Result of `mmc_identified_support()`.
#' @param library_ids Result of `mmc_identified_pathway_library()`.
#' @param p_cutoff    Nominal threshold.
#' @param engine      `"metaboanalyst"` or `"python"`.
#' @param evidence_note Why evidence was not built, or NULL.
#' @param fdr         Adjusted-p threshold.
#' @return A data.frame with columns `.MMC_OVERVIEW_COLS`; zero rows when there
#'   are no pathway results.
#' @examples
#' ora <- data.frame(Pathway = c("a", "b"), Overlap = c(2, 1),
#'                   "Pathway size" = c(5, 4), "Enrichment ratio" = c(0.4, 0.25),
#'                   p.value = c(0.01, 0.3), check.names = FALSE)
#' build_mummichog_overview(ora, engine = "python")
build_mummichog_overview <- function(ora_table, gsea_table = NULL,
                                     evidence = NULL, pathway_ids = NULL,
                                     identified = NULL,
                                     library_ids = character(0),
                                     p_cutoff = 0.05, engine = "python",
                                     evidence_note = NULL,
                                     fdr = .MMC_SUMMARY_FDR) {
  base <- build_mummichog_summary(
    ora_table, gsea_table = gsea_table, evidence = evidence,
    pathway_ids = pathway_ids, identified = identified,
    library_ids = library_ids, p_cutoff = p_cutoff, fdr = fdr,
    top_n = .Machine$integer.max)
  is_ma <- identical(engine, "metaboanalyst")
  n <- nrow(base)
  at <- if (is.data.frame(ora_table)) match(base$Pathway, ora_table$Pathway) else
    rep(NA_integer_, n)
  col_at <- function(col) {
    if (!is.data.frame(ora_table) || !col %in% names(ora_table)) {
      return(rep(NA_real_, n))
    }
    suppressWarnings(as.numeric(ora_table[[col]][at]))
  }
  ora_p   <- base[["ORA p.value"]]
  holm    <- if (is_ma) col_at("adj.p Holm FET") else rep(NA_real_, n)
  overlap <- col_at("Overlap")
  gsea_p  <- base[["GSEA P.Value"]]
  gsea_q  <- base[["GSEA padj"]]

  ora_lvl  <- .mmc_evidence_level(ora_p, holm, p_cutoff, fdr)
  gsea_lvl <- .mmc_evidence_level(gsea_p, gsea_q, p_cutoff, fdr)
  ora_status <- c(adjusted = sprintf("Holm-adjusted p < %g", fdr),
                  nominal  = sprintf("Nominal (raw p < %g)", p_cutoff),
                  none     = "Not nominal",
                  absent   = if (is_ma) "No significant EC overlap (not in ORA table)"
                             else "Not in ORA results")[ora_lvl]
  gsea_status <- if (!is.data.frame(gsea_table)) rep("GSEA not run", n) else
    c(adjusted = sprintf("BH FDR < %g", fdr),
      nominal  = sprintf("Nominal (raw p < %g)", p_cutoff),
      none     = "Not nominal",
      absent   = "Not tested by GSEA")[gsea_lvl]

  ora_nom  <- ora_lvl %in% c("adjusted", "nominal")
  gsea_nom <- gsea_lvl %in% c("adjusted", "nominal")
  feats    <- base[["Supporting features"]]
  ev <- rep("—", n)
  ev[ora_nom & !is.na(feats)] <- "See Pathway Evidence"
  ev[ora_nom & is.na(feats)]  <- if (is.null(evidence_note))
    "Not available" else "Not available: evidence not built (see Pathway Evidence)"
  ev[!ora_nom & gsea_nom] <- "GSEA-only signal: no ORA evidence shown"

  agree <- base[["Annotation agreement"]]
  agree[!ora_nom] <- "—"

  data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    "Pathway"                    = base$Pathway,
    "Pathway ID"                 = base[["Pathway ID"]],
    "ORA status"                 = unname(ora_status),
    "ORA raw p"                  = ora_p,
    "ORA Holm-adjusted p"        = holm,
    "GSEA status"                = unname(gsea_status),
    "GSEA raw p"                 = gsea_p,
    "GSEA BH FDR"                = gsea_q,
    "ORA overlap ECs (Sig.EC)"   = ifelse(ora_nom, overlap, NA_real_),
    "Features in overlap ECs"    = ifelse(ora_nom, feats, NA_real_),
    "Annotation agreement (ECs)" = agree,
    "ORA evidence"               = ev,
    "Identified-metabolite enrichment" =
      .mmc_identified_wording(base[["Identified-metabolite support"]], fdr),
    "Interpretation"             = .mmc_relabel(base[["Mummichog contribution"]],
                                                .MMC_INTERPRETATION_LABELS),
    "Nominal Mummichog signal"   = base[["Nominal Mummichog signal"]],
    "ORA raw p type"             = rep(if (is_ma) "FET" else "empirical", n)
  )
}

#' Which columns of a mummichog presentation table are p-values and counts
#'
#' The report's table helper shows p-values in scientific notation only when
#' it recognises them by name, and rounds every other number to 2 decimals.
#' The presentation tables name their columns for what they hold ("ORA raw p
#' (FET)", "Holm-adjusted FET", "GSEA BH FDR"), so this names them for it: a
#' FET of 3e-5 must not display as 0.00, and an EC count as 3.00.
#'
#' @param df A presentation table (overview, detailed ORA table, input counts,
#'   across-contrast view).
#' @return `list(pval, int)`: numeric columns holding p-values, and numeric
#'   columns holding whole numbers only.
#' @examples
#' mmc_dt_column_formats(data.frame("ORA raw p" = 3e-5, "ECs" = 3,
#'                                  check.names = FALSE))
mmc_dt_column_formats <- function(df) {
  if (!is.data.frame(df)) return(list(pval = character(0), int = character(0)))
  num <- names(df)[vapply(df, is.numeric, logical(1))]
  pval <- num[grepl("(raw p|Holm-adjusted|BH FDR|^ORA p \\(|p-value|p\\.value|adj\\.p|^P\\.Value$|^padj$)",
                    num, ignore.case = TRUE)]
  whole <- vapply(num, function(cl) {
    v <- df[[cl]][!is.na(df[[cl]])]
    length(v) > 0 && all(is.finite(v) & v == round(v))
  }, logical(1))
  list(pval = pval, int = setdiff(num[whole], pval))
}

#' Compact display of the overview for the HTML Summary
#'
#' One readable cell per test, so the Summary fits a page: ORA as
#' "<status> · FET p = … · Holm = …", GSEA as "<status> · p = … · FDR = …".
#' The numeric columns stay in the overview (and the workbook).
#'
#' @param overview Rows of `build_mummichog_overview()`.
#' @return A data.frame for `make_dt()`.
#' @noRd
mmc_overview_display <- function(overview) {
  fmt <- function(x) ifelse(is.na(x), NA_character_, formatC(x, digits = 2, format = "g"))
  ora_type <- overview[["ORA raw p type"]]
  ora <- ifelse(is.na(overview[["ORA raw p"]]), overview[["ORA status"]],
                paste0(overview[["ORA status"]], " · ", ora_type, " p = ",
                       fmt(overview[["ORA raw p"]]),
                       ifelse(is.na(overview[["ORA Holm-adjusted p"]]), "",
                              paste0(" · Holm = ", fmt(overview[["ORA Holm-adjusted p"]])))))
  gsea <- ifelse(is.na(overview[["GSEA raw p"]]), overview[["GSEA status"]],
                 paste0(overview[["GSEA status"]], " · p = ", fmt(overview[["GSEA raw p"]]),
                        " · FDR = ", fmt(overview[["GSEA BH FDR"]])))
  data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    "Pathway"                    = overview$Pathway,
    "ORA"                        = ora,
    "GSEA"                       = gsea,
    "ORA overlap ECs"            = overview[["ORA overlap ECs (Sig.EC)"]],
    "Features in overlap ECs"    = overview[["Features in overlap ECs"]],
    "Annotation agreement (ECs)" = overview[["Annotation agreement (ECs)"]],
    "Identified-metabolite enrichment" = overview[["Identified-metabolite enrichment"]],
    "Interpretation"             = overview$Interpretation
  )
}

#' The complete ORA table with every engine column named for what it holds
#'
#' For the MetaboAnalyst-compatible engine, every column of its pathway table
#' (06n), in its own row order: the EC counts the test uses, the compound
#' counts (`Hits.*`, putative compound matches, not identified metabolites),
#' each p-value and its Holm adjustment, and the overlap EC ids. For the python
#' engine, the section's ORA table under the same naming. Values are copied,
#' never recomputed.
#'
#' @param pathways  The contrast's pathway table as read by the stage readers.
#' @param ora_table The section's ORA table (`build_mummichog_pathway_table()`).
#' @return A data.frame, or NULL when there is no ORA table.
#' @examples
#' build_mummichog_ora_detail(NULL, NULL)
build_mummichog_ora_detail <- function(pathways, ora_table) {
  if (!is.data.frame(ora_table) || nrow(ora_table) == 0L) return(NULL)
  if (!.mmc_is_ma_table(pathways)) {
    return(data.frame(
      check.names = FALSE, stringsAsFactors = FALSE,
      "Pathway"                   = ora_table$Pathway,
      "ORA overlap ECs"           = ora_table$Overlap,
      "ECs in pathway"            = ora_table[["Pathway size"]],
      "Enrichment ratio"          = ora_table[["Enrichment ratio"]],
      "ORA raw p (empirical)"     = ora_table[["p.value"]]))
  }
  pw <- pathways[order(pathways$ma_rank), , drop = FALSE]
  num <- function(col) if (col %in% names(pw)) suppressWarnings(as.numeric(pw[[col]])) else
    rep(NA_real_, nrow(pw))
  data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    "Pathway"                                   = as.character(pw$pathway),
    "Pathway ID"                                = as.character(pw$pathway_id),
    "MetaboAnalyst row order"                   = num("ma_rank"),
    "ORA overlap ECs (Sig.EC)"                  = num("Sig.EC"),
    "ECs in pathway (Total.EC)"                 = num("Total.EC"),
    "Expected overlap ECs"                      = num("Expected"),
    "Pathway compounds in model (Pathway total)" = num("Pathway total"),
    "Putative compound matches, all ECs (Hits.total)" = num("Hits.total"),
    "Putative compound matches, overlap ECs (Hits.sig)" = num("Hits.sig"),
    "ORA raw p (FET)"                           = num("FET"),
    "ORA raw p (EASE)"                          = num("EASE"),
    "ORA p (Gamma-adjusted FET)"                = num("Gamma"),
    "ORA p (empirical)"                         = num("Empirical"),
    "Permutations with p <= FET (Emp.Hits)"     = num("Emp.Hits"),
    "Holm-adjusted FET"                         = num("Holm_adj_FET"),
    "Holm-adjusted EASE"                        = num("Holm_adj_EASE"),
    "Holm-adjusted Gamma"                       = num("Holm_adj_Gamma"),
    "Overlap EC IDs"                            = as.character(pw$ec_hits %||% NA_character_)
  )
}

#' What entered the mummichog engine for one contrast
#'
#' Read, never recomputed: the engine's own counts from its run record
#' (`ma_run_info.json`, MetaboAnalyst-compatible engine), the size of the
#' contrast's DE table, and the dataset-wide number of features passing the
#' overall DE selection (`pass_any_contrast`), which is reported so that it is
#' not mistaken for the mummichog input. The python engine records none of
#' the engine counts; they are NA.
#'
#' @param pathways The contrast's pathway table (carries `attr(, "ma_info")`).
#' @param de_table The contrast's DE table, or NULL.
#' @param de_res   The DE results (for `summary_df$pass_any_contrast`), or NULL.
#' @return A one-row data.frame.
#' @examples
#' mmc_input_counts(NULL, NULL, NULL)
mmc_input_counts <- function(pathways, de_table, de_res) {
  info <- attr(pathways, "ma_info")
  ct   <- info$counts %||% list()
  num1 <- function(x) if (is.null(x) || length(x) != 1L) NA_real_ else
    suppressWarnings(as.numeric(x))
  at_de <- if (is.data.frame(de_table)) nrow(de_table) else NA_real_
  given <- num1(ct$features_given)
  sdf <- if (is.list(de_res) && !is.data.frame(de_res)) de_res$summary_df else NULL
  data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    "Features at DE (this contrast)"              = at_de,
    "Passed to the engine"                        = given,
    "Not passed (missing m/z, RT, p or statistic)" = at_de - given,
    "Removed: m/z out of range"                   = num1(ct$removed_out_of_range),
    "Removed: duplicate m/z"                      = num1(ct$removed_duplicate_mz),
    "Eligible input features"                     = num1(ct$features_kept),
    "Significant input features"                  = num1(ct$significant_features),
    "Significance rule"                           = if (is.null(info$params$p_cutoff))
      NA_character_ else sprintf("raw DE p <= %s", info$params$p_cutoff),
    "Empirical compounds"                         = num1(ct$empirical_compounds),
    "Significant ECs"                             = num1(ct$significant_ecs),
    "Pathways tested"                             = num1(ct$pathways_tested),
    "Pathways in ORA table"                       = num1(ct$pathways_reported),
    "Features in DE results (all contrasts)"      = if (is.data.frame(sdf)) nrow(sdf) else NA_real_,
    "Passing overall DE selection (pass_any_contrast; not the Mummichog input)" =
      if (is.data.frame(sdf) && "pass_any_contrast" %in% names(sdf))
        sum(sdf$pass_any_contrast == 1, na.rm = TRUE) else NA_real_
  )
}
