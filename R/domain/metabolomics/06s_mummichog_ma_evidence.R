# R/domain/metabolomics/06s_mummichog_ma_evidence.R
#
# Pathway Evidence for the MetaboAnalyst-compatible mummichog engine (06j-06n).
#
#     Pathway -> supporting EC -> pathway-matching candidate
#             -> measured feature -> original annotation -> agreement
#
# The python engine's evidence (06f) reads that engine's EC tables. This file
# reads what the MetaboAnalyst-compatible engine wrote for a contrast and hands
# the same per-EC shape to the shared assembler, `.mmc_assemble_evidence()`:
#
#   ma_pathway_results.tsv   which ECs support a pathway: `ec_hits`, the
#                            pathway's significant ECs (their count is Sig.EC);
#   ma_matched_features.tsv  each kept (feature, compound, ion form) match with
#                            its final EC, feature id and polarity;
#   ma_input.tsv             the feature's m/z, RT, p-value and the signed
#                            statistic the engine was given;
#   ma_run_info.json         the sha256 of the model the contrast ran on;
# plus the metabolic model (pathway -> compounds, names, KEGG ids), the DE
# table (log2FC) and the dataset's own feature annotations.
#
# Nothing here runs or re-runs the engine, and no statistic, EC grouping or
# pathway membership is recomputed: the ECs a pathway lists are the engine's
# `ec_hits`, and the support counts are checked against its Sig.EC.
#
# Where this differs from the python evidence, and why:
#   * A MetaboAnalyst EC can hold features matched to different compounds
#     (ECs of isobaric matches are merged), so a feature can sit in an EC that
#     reaches a pathway through a compound the feature never matched. Each
#     feature's annotation is therefore compared only with the pathway
#     candidates IT matched; a feature with none reads "Not assessed", and the
#     feature table says so (`Matched via pathway candidate`). EC membership
#     and the support counts are not changed by this.
#   * In mixed mode an EC can hold features of both polarities; each feature
#     keeps its own polarity and ion form.
#   * The engine keeps only the smaller-p feature of a set sharing one m/z
#     (across polarities). A dropped feature was never analysed, so it is in no
#     EC and in no evidence row. Which input features were dropped this way is
#     not recorded per feature by the engine, and is not re-derived here.
#   * Evidence re-reads the model, so it is built only when that model is the
#     one the contrast was analysed with (sha256 from ma_run_info.json).


#' Read one contrast's MetaboAnalyst-compatible matches into per-EC frames
#'
#' Reads `ma_matched_features.tsv` and `ma_input.tsv` (every column as text, so
#' a feature id such as `"0012"` keeps its leading zeros) and checks that every
#' match traces back to exactly one input feature: each match row names a
#' feature id, every id is one of the input features, the input ids are unique,
#' each m/z key belongs to one feature and each feature to one m/z key, and a
#' feature's polarity agrees between the two files. Any failure stops — evidence
#' with an untraceable or ambiguous feature would misattribute support.
#'
#' @param files    One contrast's files (found by name).
#' @param de_table Optional DE table with `feature_id` and `logFC`; supplies
#'   `log2fc` (NA when absent).
#' @return A list:
#'   \describe{
#'     \item{cand}{One row per (EC, compound): `EID`, `compound_id`,
#'       `compound_name` (NA: the engine writes no names, the model supplies
#'       them).}
#'     \item{feat}{One row per (EC, feature): `EID`, `feature_id`, `polarity`,
#'       `mz`, `retention_time`, `p_value`, `statistic`, `log2fc`, `adduct`
#'       (the feature's ion forms in this EC) and `matched_compounds` (the
#'       compounds it matched in this EC, ";"-joined). Within an EC, rows are
#'       ordered by p-value, then feature id.}
#'   }
#' @examples
#' \dontrun{
#' x <- read_mummichog_ma_evidence_inputs(files_for_one_contrast)
#' }
read_mummichog_ma_evidence_inputs <- function(files, de_table = NULL) {
  pick <- function(name) {
    f <- files[basename(files) == name]
    if (length(f) == 0L) {
      .mmc_stop("this contrast has no '", name, "' among its MetaboAnalyst-",
                "compatible result files.")
    }
    f[[1]]
  }
  read_chr <- function(f) as.data.frame(readr::read_tsv(
    f, show_col_types = FALSE, progress = FALSE, na = c("", "NA"),
    col_types = readr::cols(.default = readr::col_character())),
    stringsAsFactors = FALSE)
  m   <- read_chr(pick(.MMC_MA_FEATURES_FILE))
  inp <- read_chr(pick(.MMC_MA_INPUT_FILE))

  need <- function(df, cols, name) {
    miss <- setdiff(cols, names(df))
    if (length(miss) > 0L) {
      .mmc_stop("'", name, "' lacks column(s) ", paste(miss, collapse = ", "),
                ". Present: ", paste(names(df), collapse = ", "), ".")
    }
  }
  need(m, c("feature_id", "polarity", "Query.Mass", "Matched.Compound",
            "Matched.Form", "Empirical.Compound"), .MMC_MA_FEATURES_FILE)
  need(inp, c("feature_id", "mz", "rt", "p_value", "statistic"),
       .MMC_MA_INPUT_FILE)

  # -- traceability ---------------------------------------------------------------
  no_id <- is.na(m$feature_id) | !nzchar(m$feature_id)
  if (any(no_id)) {
    .mmc_stop(sum(no_id), " matched row(s) in '", .MMC_MA_FEATURES_FILE,
              "' carry no feature id (m/z e.g. ",
              paste(utils::head(unique(m$Query.Mass[no_id]), 3), collapse = ", "),
              "), so their support cannot be traced to a measured feature.")
  }
  dup_in <- unique(inp$feature_id[duplicated(inp$feature_id)])
  if (length(dup_in) > 0L) {
    .mmc_stop(length(dup_in), " feature id(s) occur more than once in '",
              .MMC_MA_INPUT_FILE, "' (e.g. ",
              paste(utils::head(dup_in, 3), collapse = ", "),
              "), so a match cannot be traced to one feature.")
  }
  unknown <- setdiff(m$feature_id, inp$feature_id)
  if (length(unknown) > 0L) {
    .mmc_stop(length(unknown), " matched feature id(s) are not in '",
              .MMC_MA_INPUT_FILE, "' (e.g. ",
              paste(utils::head(unknown, 3), collapse = ", "), ").")
  }
  per_mz <- tapply(m$feature_id, m$Query.Mass, function(x) length(unique(x)))
  per_id <- tapply(m$Query.Mass, m$feature_id, function(x) length(unique(x)))
  if (any(per_mz > 1L) || any(per_id > 1L)) {
    .mmc_stop("'", .MMC_MA_FEATURES_FILE, "' maps an m/z to several feature ",
              "ids or a feature id to several m/z (e.g. ",
              paste(utils::head(c(names(per_mz)[per_mz > 1L],
                                  names(per_id)[per_id > 1L]), 3),
                    collapse = ", "),
              "); feature ids are not traceable.")
  }
  idx <- match(m$feature_id, inp$feature_id)
  if ("polarity" %in% names(inp)) {
    ip <- inp$polarity[idx]
    bad_pol <- !is.na(ip) & (is.na(m$polarity) | ip != m$polarity)
    if (any(bad_pol)) {
      .mmc_stop(length(unique(m$feature_id[bad_pol])), " feature(s) have a ",
                "different polarity in '", .MMC_MA_FEATURES_FILE, "' and '",
                .MMC_MA_INPUT_FILE, "' (e.g. ",
                paste(utils::head(unique(m$feature_id[bad_pol]), 3),
                      collapse = ", "), ").")
    }
  }

  # -- per EC ---------------------------------------------------------------------
  cand <- unique(data.frame(EID = m$Empirical.Compound,
                            compound_id = m$Matched.Compound,
                            stringsAsFactors = FALSE))
  cand$compound_name <- rep(NA_character_, nrow(cand))
  rownames(cand) <- NULL

  key <- paste(m$Empirical.Compound, m$feature_id, sep = "\u001f")
  grp <- split(seq_len(nrow(m)), factor(key, levels = unique(key)))
  first <- vapply(grp, `[[`, integer(1), 1L)
  num <- function(x) suppressWarnings(as.numeric(x))
  fid <- m$feature_id[first]
  at  <- match(fid, inp$feature_id)

  log2fc <- rep(NA_real_, length(fid))
  if (is.data.frame(de_table) && all(c("feature_id", "logFC") %in% names(de_table))) {
    log2fc <- num(de_table$logFC[match(fid, as.character(de_table$feature_id))])
  }

  feat <- data.frame(
    EID               = m$Empirical.Compound[first],
    feature_id        = fid,
    polarity          = m$polarity[first],
    mz                = num(inp$mz[at]),
    retention_time    = num(inp$rt[at]),
    p_value           = num(inp$p_value[at]),
    statistic         = num(inp$statistic[at]),
    log2fc            = log2fc,
    adduct            = vapply(grp, function(ix) .mmc_join_unique(m$Matched.Form[ix]),
                               character(1), USE.NAMES = FALSE),
    matched_compounds = vapply(grp, function(ix)
                                 paste(unique(m$Matched.Compound[ix]), collapse = ";"),
                               character(1), USE.NAMES = FALSE),
    stringsAsFactors = FALSE)
  # Explicit, locale-independent order: the file's row order is not a contract.
  feat <- feat[order(feat$EID, feat$p_value, feat$feature_id, method = "radix"), ,
               drop = FALSE]
  rownames(feat) <- NULL
  list(cand = cand, feat = feat)
}

#' Trace the supporting evidence of a MetaboAnalyst-compatible mummichog result
#'
#' The MetaboAnalyst-compatible counterpart of
#' `build_mummichog_pathway_evidence()`, returning the same three frames. A
#' pathway's supporting ECs are the engine's own `ec_hits` (its significant
#' ECs); their candidates are the compounds matched in the EC that belong to
#' the pathway in the model; each EC lists every measured feature the engine
#' grouped into it, significant or not. A feature's agreement compares its
#' annotation only with the pathway candidates it matched itself (see the file
#' header).
#'
#' Refuses (stops) rather than return evidence it cannot vouch for: no model,
#' a model whose sha256 differs from the one the contrast ran on (or that
#' cannot be checked), an untraceable feature (see
#' `read_mummichog_ma_evidence_inputs()`), an `ec_hits` EC absent from the
#' matched features, or support counts that disagree with the pathway's Sig.EC.
#'
#' @param pathways One contrast's pathway table from
#'   `read_mummichog_ma_pathways()` (carries `attr(, "ma_info")`).
#' @param files    That contrast's files.
#' @param model    Model index from `mmc_load_model_pathways()`.
#' @param annot    Normalised annotations from `normalize_metab_annotation()`.
#' @param de_table Optional DE table for the contrast (`feature_id`, `logFC`).
#' @param model_sha256 sha256 of the model file `model` was read from;
#'   defaults to the one `mmc_load_model_pathways()` attaches.
#' @param contrast Contrast label, used only in the log2FC note.
#' @return The list described in `build_mummichog_pathway_evidence()`, with
#'   `attr(, "log2fc_note")` saying why `Feature log2FC (DE)` is empty or
#'   incomplete (NULL when every evidence feature has one); the
#'   feature table carries `Polarity`, `Matched compound(s)`,
#'   `Pathway candidate(s) matched by feature`, `Matched via pathway candidate`,
#'   `Feature statistic (MetaboAnalyst input)` and `Feature log2FC (DE)` in
#'   place of the python `Feature log2FC (mummichog input)`, and the EC table
#'   carries `Polarities`.
#' @examples
#' \dontrun{
#' ev <- build_mummichog_ma_pathway_evidence(pw, files, model, annot, de_tbl)
#' }
build_mummichog_ma_pathway_evidence <- function(pathways, files, model, annot,
                                                de_table = NULL,
                                                model_sha256 = attr(model, "sha256"),
                                                contrast = NULL) {
  if (!.mmc_is_ma_table(pathways) || nrow(pathways) == 0L) {
    .mmc_stop("not a MetaboAnalyst-compatible pathway table.")
  }
  if (!"ec_hits" %in% names(pathways)) {
    .mmc_stop("the pathway table has no 'ec_hits' column, so the ECs ",
              "supporting each pathway are unknown.")
  }
  if (is.null(model) || length(model$pathways) == 0L) {
    .mmc_stop("the metabolic model could not be read, so pathway membership ",
              "is unknown.")
  }
  recorded <- attr(pathways, "ma_info")$model$sha256
  if (is.null(recorded) || length(recorded) != 1L || is.na(recorded) ||
      !nzchar(recorded)) {
    .mmc_stop("the run record (", .MMC_MA_INFO_FILE, ") names no model ",
              "sha256, so the model evidence would be read from cannot be ",
              "checked against the one the analysis used.")
  }
  if (is.null(model_sha256) || length(model_sha256) != 1L ||
      is.na(model_sha256) || !nzchar(model_sha256)) {
    .mmc_stop("the sha256 of the currently configured metabolic model is ",
              "unknown, so it cannot be checked against the one the analysis ",
              "used.")
  }
  if (!identical(tolower(model_sha256), tolower(recorded))) {
    .mmc_stop("the currently configured metabolic model (sha256 ",
              substr(model_sha256, 1, 12), "...) is not the model this contrast ",
              "was analysed with (sha256 ", substr(recorded, 1, 12), "..., in ",
              .MMC_MA_INFO_FILE, "). Evidence would describe pathway membership ",
              "the analysis did not use. Restore that model, or re-run the ",
              "mummichog stage with the current one.")
  }

  inp <- read_mummichog_ma_evidence_inputs(files, de_table = de_table)
  ec_lists <- lapply(as.character(pathways$ec_hits), .mmc_split_cell, sep = ";")
  missing <- setdiff(unlist(ec_lists, use.names = FALSE), inp$cand$EID)
  if (length(missing) > 0L) {
    .mmc_stop(length(missing), " EC(s) listed in the pathway table's ec_hits ",
              "are not in '", .MMC_MA_FEATURES_FILE, "' (e.g. ",
              paste(utils::head(missing, 3), collapse = ", "), ").")
  }

  ev <- .mmc_assemble_evidence(
    pathways, ec_lists,
    cand_by_ec = split(inp$cand, inp$cand$EID),
    feat_by_ec = split(inp$feat, inp$feat$EID),
    model = model, annot = annot, own_matches = TRUE)

  # Each pathway's supporting ECs must be exactly the ECs the engine counted.
  ps  <- if (is.null(ev)) NULL else ev$pathway_summary
  got <- if (is.null(ps)) rep(0, nrow(pathways)) else
    ps[["Supporting ECs"]][match(pathways$pathway, ps$Pathway)]
  got[is.na(got)] <- 0
  want <- suppressWarnings(as.numeric(pathways[["Sig.EC"]]))
  off <- is.na(want) | got != want
  if (any(off)) {
    .mmc_stop("the evidence finds a different number of supporting ECs than ",
              "the engine reported (Sig.EC) for ", sum(off), " pathway(s) (e.g. '",
              paste(utils::head(pathways$pathway[off], 2), collapse = "', '"),
              "'). The model's pathway membership does not match the run.")
  }
  # Why `Feature log2FC (DE)` is empty or incomplete, when it is: carried with
  # the evidence so the report can say so instead of showing bare blanks.
  attr(ev, "log2fc_note") <- .mmc_ma_log2fc_note(de_table, ev$feature_table,
                                                 contrast)
  ev
}

#' Explain an empty or incomplete `Feature log2FC (DE)` column
#'
#' The log2FC is joined from the DE results by feature id, so it can be
#' missing for a reason that has nothing to do with the evidence itself: no DE
#' table could be matched to the contrast, the table has no `feature_id` or
#' `logFC` column, or some evidence features have no logFC in it. Every other
#' evidence column is unaffected.
#'
#' @param de_table      The DE table the evidence was given, or NULL.
#' @param feature_table The evidence's feature table.
#' @param contrast      Contrast label, for the message (may be NULL).
#' @return One sentence, or NULL when every evidence feature has a log2FC.
#' @noRd
.mmc_ma_log2fc_note <- function(de_table, feature_table, contrast = NULL) {
  ct <- if (length(contrast) == 1L && !is.na(contrast) && nzchar(contrast))
    sprintf(" '%s'", contrast) else ""
  tail <- paste0("Every other evidence column, including Feature statistic ",
                 "(MetaboAnalyst input), is unaffected.")
  if (!is.data.frame(de_table)) {
    return(paste0("No differential-abundance table could be matched to ",
                  "contrast", ct, ", so the Feature log2FC (DE) column is ",
                  "empty. ", tail))
  }
  miss <- setdiff(c("feature_id", "logFC"), names(de_table))
  if (length(miss) > 0L) {
    return(paste0("The differential-abundance table for contrast", ct,
                  " has no ", paste(miss, collapse = " or "), " column, so ",
                  "the Feature log2FC (DE) column is empty. ", tail))
  }
  ids <- unique(as.character(feature_table$Feature))
  no_fc <- unique(as.character(
    feature_table$Feature[is.na(feature_table[["Feature log2FC (DE)"]])]))
  if (length(no_fc) == 0L) return(NULL)
  paste0(length(no_fc), " of ", length(ids), " evidence feature(s) have no ",
         "log2FC in the differential-abundance table for contrast", ct,
         "; their Feature log2FC (DE) cells are empty. ", tail)
}

#' Build one MetaboAnalyst-compatible section's evidence, or say why not
#'
#' Wraps `build_mummichog_ma_pathway_evidence()` for the report sections: a
#' refusal or error becomes a warning and a reason the report shows in the
#' contrast's Pathway Evidence tab, so evidence never disappears silently. The
#' ORA and GSEA results are not touched either way.
#'
#' @param pathways One contrast's MetaboAnalyst-compatible pathway table.
#' @param files    That contrast's files, or NULL.
#' @param model    Model index (with its sha256 attribute), or NULL.
#' @param annot    Normalised annotations.
#' @param de_table DE table for the contrast, or NULL.
#' @param contrast Contrast label, for the warning.
#' @return `list(evidence, note, log2fc_note)`: the evidence (NULL when not
#'   built), the reason it was not built (NULL when it was), and why its
#'   `Feature log2FC (DE)` column is empty or incomplete (NULL when it is
#'   complete or there is no evidence). Either note also warns.
#' @examples
#' \dontrun{
#' mmc_ma_section_evidence(pw, files, model, annot, de_tbl, "A_vs_B")
#' }
mmc_ma_section_evidence <- function(pathways, files, model, annot,
                                    de_table = NULL, contrast = "") {
  ev <- if (length(files) == 0L) {
    simpleError("no MetaboAnalyst-compatible result files were found for this contrast.")
  } else {
    tryCatch(build_mummichog_ma_pathway_evidence(pathways, files, model, annot,
                                                 de_table = de_table,
                                                 contrast = contrast),
             error = function(e) e)
  }
  if (!inherits(ev, "error")) {
    fc_note <- attr(ev, "log2fc_note")
    if (!is.null(fc_note)) {
      warning("mummichog evidence: contrast '", contrast, "': ", fc_note,
              call. = FALSE)
    }
    return(list(evidence = ev, note = NULL, log2fc_note = fc_note))
  }
  reason <- sub("^\\[mummichog\\] ", "", conditionMessage(ev))
  warning("mummichog evidence: contrast '", contrast, "': no Pathway Evidence ",
          "built (ORA and GSEA are unaffected): ", reason, call. = FALSE)
  list(evidence = NULL, note = reason, log2fc_note = NULL)
}
