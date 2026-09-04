# R/domain/metabolomics/06f_mummichog_evidence.R
#
# Evidence-tracing layer for the pinned mummichog v2 stage (06c/05b). It answers
# the question the pathway table alone cannot: *which measured features actually
# support this enriched pathway, and does the identity mummichog used to place
# them there agree with our own annotation?*
#
#     Pathway -> EmpiricalCompound -> pathway-matching candidate
#             -> original measured feature -> original annotation -> agreement
#
# This module runs NO analysis. It parses mummichog's own result tables, joins
# them to the metabolic model and to our feature annotations, and returns pure
# report-ready frames. The engine's statistics are never recomputed or filtered.
#
# ---------------------------------------------------------------------------
# Why we do NOT use mummichog's "Best guess" / face_compound
# ---------------------------------------------------------------------------
# Verified against the mummichog 2.7.0 sources:
#
#   * get_user_data.py::EmpiricalCompound.designate_face_cpd() sets
#     `face_compound = chosen_compounds[-1]` — its own docstring says one
#     candidate is "arbitrarily designated" when several are suggested.
#   * functional_analysis.py::collect_hit_Trios() fills `chosen_compounds` from
#     the UNION of every pathway that passed the significance cutoff, not from
#     the pathway currently being inspected.
#   * reporting.py::export_pathway_enrichtest() writes that same
#     `chosen_compounds` into the `overlap_features (id)` column.
#
# So neither `face_compound` nor `overlap_features (id)` is an evidence-ranked,
# per-pathway identification. We derive the pathway-matching candidate(s)
# ourselves: all candidate compounds carried by the EmpiricalCompound,
# intersected with the compounds of the pathway under inspection. Several
# candidates can survive that intersection — we keep every one of them.
#
# ---------------------------------------------------------------------------
# Where pathway membership comes from
# ---------------------------------------------------------------------------
# mummichog does not export pathway -> compound membership (in reporting.py's
# JSON dump, `'all_compounds': P.cpds` is commented out), so it has to be read
# from the metabolic model itself — the very model the stage ran on, resolved
# through the existing mmc_select_model() (06d). The built-in `human_mfn` model
# lives inside the Python package and is not readable from R; with that model
# the evidence layer reports itself unavailable and the report section is simply
# omitted rather than guessed at.
#
# Dependencies (R): jsonlite, readr. Both already in renv.lock.


# ==== internal helpers ======================================================

#' Pick one named result file out of a mummichog output file list
#'
#' The 06c readers each locate their own file inline; new readers here need the
#' same lookup for files 06c does not read, so it is factored into one place
#' instead of being repeated per reader.
#'
#' @param files Character vector of mummichog output files.
#' @param name  Exact basename to find.
#' @return The first matching path, or NULL when absent.
#' @noRd
.mmc_pick_result_file <- function(files, name) {
  if (length(files) == 0) return(NULL)
  hit <- files[basename(files) == name]
  if (length(hit) == 0) NULL else hit[[1]]
}

#' Split a delimited cell into a trimmed, non-empty character vector
#' @param x   A single character scalar (may be NA).
#' @param sep Fixed separator.
#' @return Character vector (possibly empty).
#' @noRd
.mmc_split_cell <- function(x, sep) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return(character(0))
  out <- trimws(strsplit(as.character(x), sep, fixed = TRUE)[[1]])
  out[nzchar(out)]
}

#' Normalise a metabolite name for a conservative string comparison
#'
#' Case, whitespace and punctuation differ freely between vendor software and
#' metabolic models ("3'-AMP" / "3-AMP", "L-Glutamate" / "L-glutamate"), so a
#' raw string comparison would report conflicts that are really the same
#' compound. Lower-casing and dropping every non-alphanumeric character is
#' deliberately the ONLY normalisation applied: nothing here infers identity
#' from mass, formula or m/z.
#'
#' @param x Character vector of names.
#' @return Character vector of normalised keys ("" for unusable input).
#' @noRd
.mmc_norm_name <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[is.na(x)] <- ""
  gsub("[^a-z0-9]", "", x)
}

#' Extract the KEGG compound id from a raw annotation cell
#'
#' Cells arrive as "C00031", "cpd:C00031" or "C00031;C00267"; take the first
#' KEGG compound id present and nothing else.
#'
#' @param x Character vector of raw cells.
#' @return Character vector of KEGG ids, NA where none is present.
#' @noRd
.mmc_extract_kegg <- function(x) {
  x   <- as.character(x)
  hit <- grepl("C[0-9]{5}", x)
  out <- rep(NA_character_, length(x))
  out[hit] <- sub(".*?(C[0-9]{5}).*", "\\1", x[hit], perl = TRUE)
  out
}

#' Extract the HMDB accession from a raw annotation cell
#' @param x Character vector of raw cells.
#' @return Character vector of HMDB ids, NA where none is present.
#' @noRd
.mmc_extract_hmdb <- function(x) {
  x   <- as.character(x)
  hit <- grepl("HMDB[0-9]+", x, ignore.case = TRUE)
  out <- rep(NA_character_, length(x))
  out[hit] <- toupper(sub(".*?(HMDB[0-9]+).*", "\\1", x[hit],
                          perl = TRUE, ignore.case = TRUE))
  out
}


# ==== metabolic model: pathway -> compounds =================================

#' Read pathway membership and compound definitions from a metabolic model JSON
#'
#' Understands both shapes a mummichog `-n <model>` file can take, detected by
#' their keys (never guessed):
#'
#' \describe{
#'   \item{Azimuth}{`list_of_pathways` / `list_of_reactions` /
#'     `list_of_compounds` — the format `model_ref` / `model_json` publish and
#'     mummichog's `models.py::read_user_json_model()` consumes. A pathway's
#'     compounds are the union of its reactions' reactants and products, exactly
#'     as `json_convert_azmuth_mummichog()` derives them.}
#'   \item{mummichog 2 native}{`metabolic_pathways` (each with `cpds`) plus an
#'     optional `dict_cpds_def` name map.}
#' }
#'
#' @param path Path to the model JSON.
#' @return A list with:
#'   \describe{
#'     \item{pathways}{Named list: pathway name -> character vector of compound ids.}
#'     \item{compound_names}{Named character vector: compound id -> name (may be
#'       a `";"`-separated synonym list, as the models themselves store it).}
#'     \item{compound_kegg}{Named character vector: compound id -> KEGG compound
#'       id, `NA` when the model carries none.}
#'     \item{format}{`"azimuth"` or `"mummichog2"`.}
#'   }
#'   Aborts when the file is unreadable or matches neither shape — a malformed
#'   model must never silently degrade into an empty evidence section.
read_mummichog_model_pathways <- function(path) {
  if (is.null(path) || !is.character(path) || length(path) != 1L ||
      !nzchar(path) || !file.exists(path)) {
    .mmc_stop("metabolic model JSON not found: '",
              if (is.null(path)) "<NULL>" else paste(path, collapse = ", "), "'.")
  }
  jm <- tryCatch(
    jsonlite::fromJSON(path, simplifyVector = FALSE),
    error = function(e) .mmc_stop("could not parse the metabolic model '", path,
                                  "': ", conditionMessage(e))
  )

  if (!is.null(jm$list_of_pathways)) {
    return(.mmc_model_from_azimuth(jm, path))
  }
  if (!is.null(jm$metabolic_pathways)) {
    return(.mmc_model_from_mummichog2(jm))
  }
  .mmc_stop("the metabolic model '", path, "' has neither 'list_of_pathways' ",
            "(Azimuth format) nor 'metabolic_pathways' (mummichog 2 format). ",
            "Top-level keys: ", paste(names(jm), collapse = ", "), ".")
}

#' Build the model index from an Azimuth-format model
#' @param jm   Parsed JSON.
#' @param path Source path (for error messages).
#' @return See `read_mummichog_model_pathways()`.
#' @noRd
.mmc_model_from_azimuth <- function(jm, path) {
  rxns <- jm$list_of_reactions %||% list()
  if (length(rxns) == 0) {
    .mmc_stop("the Azimuth model '", path, "' has no 'list_of_reactions', so ",
              "pathway compound membership cannot be derived.")
  }
  # reaction id -> compounds (reactants + products), as mummichog derives it
  rxn_cpds <- lapply(rxns, function(r) {
    unique(c(unlist(r$reactants, use.names = FALSE),
             unlist(r$products,  use.names = FALSE)))
  })
  names(rxn_cpds) <- vapply(rxns, function(r) as.character(r$id %||% NA_character_),
                            character(1))

  pw <- lapply(jm$list_of_pathways, function(p) {
    ids <- as.character(unlist(p$list_of_reactions, use.names = FALSE))
    unique(unlist(rxn_cpds[ids[ids %in% names(rxn_cpds)]], use.names = FALSE))
  })
  names(pw) <- vapply(jm$list_of_pathways,
                      function(p) as.character(p$name %||% p$id %||% NA_character_),
                      character(1))
  pw <- pw[!is.na(names(pw)) & nzchar(names(pw))]

  cpds <- jm$list_of_compounds %||% list()
  cid  <- vapply(cpds, function(c) as.character(c$id %||% NA_character_), character(1))
  cnm  <- vapply(cpds, function(c) as.character(c$name %||% NA_character_), character(1))
  # KEGG: the compound id itself when it is a KEGG accession, else whatever the
  # model's `identifiers` block declares (Azimuth uses "kegg.compound").
  ckegg <- vapply(cpds, function(c) {
    idf <- c$identifiers %||% list()
    raw <- c(idf[["kegg.compound"]], idf[["kegg"]], idf[["KEGG"]], c$id)
    raw <- as.character(unlist(raw, use.names = FALSE))
    k   <- .mmc_extract_kegg(raw)
    k   <- k[!is.na(k)]
    if (length(k) == 0) NA_character_ else k[[1]]
  }, character(1))

  list(
    pathways       = pw,
    compound_names = stats::setNames(cnm, cid),
    compound_kegg  = stats::setNames(ckegg, cid),
    format         = "azimuth"
  )
}

#' Build the model index from a mummichog-2 native model
#' @param jm Parsed JSON.
#' @return See `read_mummichog_model_pathways()`.
#' @noRd
.mmc_model_from_mummichog2 <- function(jm) {
  pw <- lapply(jm$metabolic_pathways, function(p) {
    unique(as.character(unlist(p$cpds, use.names = FALSE)))
  })
  names(pw) <- vapply(jm$metabolic_pathways,
                      function(p) as.character(p$name %||% p$id %||% NA_character_),
                      character(1))
  pw <- pw[!is.na(names(pw)) & nzchar(names(pw))]

  defs <- jm$dict_cpds_def %||% list()
  cnm  <- stats::setNames(
    vapply(defs, function(v) as.character(v)[1], character(1)),
    names(defs)
  )
  all_ids <- unique(c(names(cnm), unlist(pw, use.names = FALSE)))
  list(
    pathways       = pw,
    compound_names = cnm[intersect(all_ids, names(cnm))],
    compound_kegg  = stats::setNames(.mmc_extract_kegg(all_ids), all_ids),
    format         = "mummichog2"
  )
}

#' Load the model pathway index for the model this config runs mummichog on
#'
#' Resolves the model exactly as the stage does — via `mmc_select_model()` (06d),
#' so a `model_ref` is served from its verified cache — and reads its pathway
#' membership. Returns `NULL` (with a message, never an error) when the evidence
#' layer genuinely cannot be built: mummichog's built-in `human_mfn` model ships
#' inside the Python package and is not a file R can read, and a project may not
#' have the model available at report time. The report then omits the evidence
#' section instead of showing guessed identities.
#'
#' @param config    Full pipeline config.
#' @param cache_dir Model cache directory (as used by the stage).
#' @return The list from `read_mummichog_model_pathways()`, or `NULL`.
mmc_load_model_pathways <- function(config,
                                    cache_dir = "envs/mummichog-models") {
  mummi_cfg <- config$modes$metabolomics$enrichment$mummichog %||% list()
  network <- tryCatch(
    mmc_select_model(mummi_cfg,
                     organism  = config$modes$metabolomics$organism,
                     cache_dir = cache_dir),
    error = function(e) {
      message("mummichog evidence: could not resolve the metabolic model (",
              conditionMessage(e), ") — skipping the evidence layer")
      NULL
    }
  )
  if (is.null(network)) return(NULL)
  if (identical(network, "human_mfn") || !file.exists(network)) {
    message("mummichog evidence: pathway membership is unavailable for the ",
            "built-in '", network, "' model (it lives inside the Python ",
            "package, not in a readable JSON) — skipping the evidence layer. ",
            "Configure model_ref / model_json to enable it.")
    return(NULL)
  }
  tryCatch(read_mummichog_model_pathways(network), error = function(e) {
    warning("mummichog evidence: ", conditionMessage(e), call. = FALSE)
    NULL
  })
}


# ==== mummichog result readers (EC level) ===================================

#' Read every candidate compound of every EmpiricalCompound
#'
#' Parses `ListOfEmpiricalCompounds.tsv` into long form: one row per
#' (EmpiricalCompound, candidate compound). Verified against mummichog 2.7.0's
#' `reporting.py::export_EmpiricalCompounds()`: `compounds` is the FULL candidate
#' list of the EC, `";"`-joined, and `compound_names` holds the matching names
#' `"$"`-joined in the same order (each name may itself be a `";"`-separated
#' synonym list, straight from the model's `dict_cpds_def`).
#'
#' @param files Character vector of mummichog output files.
#' @return A data.frame with columns `EID`, `compound_id`, `compound_name`,
#'   `str_row_ion`, `massfeature_rows`; zero rows when the file is absent (the
#'   evidence section then renders empty rather than failing).
read_mummichog_ec_candidates <- function(files) {
  empty <- data.frame(EID = character(0), compound_id = character(0),
                      compound_name = character(0), str_row_ion = character(0),
                      massfeature_rows = character(0),
                      stringsAsFactors = FALSE)
  f <- .mmc_pick_result_file(files, "ListOfEmpiricalCompounds.tsv")
  if (is.null(f)) return(empty)

  raw <- tryCatch(
    readr::read_tsv(f, show_col_types = FALSE,
                    col_types = readr::cols(.default = readr::col_character()),
                    name_repair = "unique_quiet"),
    error = function(e) .mmc_stop("could not read ListOfEmpiricalCompounds.tsv: ",
                                  conditionMessage(e))
  )
  miss <- setdiff(c("EID", "compounds"), names(raw))
  if (length(miss) > 0) {
    .mmc_stop("ListOfEmpiricalCompounds.tsv missing column(s): ",
              paste(miss, collapse = ", "), ". Present: ",
              paste(names(raw), collapse = ", "))
  }
  if (nrow(raw) == 0) return(empty)

  names_col <- if ("compound_names" %in% names(raw)) raw[["compound_names"]] else
    rep(NA_character_, nrow(raw))
  rows_col  <- if ("massfeature_rows" %in% names(raw)) raw[["massfeature_rows"]] else
    rep(NA_character_, nrow(raw))
  ion_col   <- if ("str_row_ion" %in% names(raw)) raw[["str_row_ion"]] else
    rep(NA_character_, nrow(raw))

  parts <- lapply(seq_len(nrow(raw)), function(i) {
    ids <- .mmc_split_cell(raw[["compounds"]][i], ";")
    if (length(ids) == 0) return(NULL)
    nms <- .mmc_split_cell(names_col[i], "$")
    # Positional pairing is mummichog's own contract; pad rather than recycle so
    # a truncated name list can never mislabel a candidate.
    if (length(nms) < length(ids)) {
      nms <- c(nms, rep(NA_character_, length(ids) - length(nms)))
    }
    data.frame(
      EID              = as.character(raw[["EID"]][i]),
      compound_id      = ids,
      compound_name    = nms[seq_along(ids)],
      str_row_ion      = as.character(ion_col[i]),
      massfeature_rows = as.character(rows_col[i]),
      stringsAsFactors = FALSE
    )
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0) return(empty)
  do.call(rbind, parts)
}

#' Read the measured features behind every EmpiricalCompound
#'
#' Parses `userInput_to_EmpiricalCompounds.tsv` into one row per
#' (EmpiricalCompound, input feature), carrying the per-feature statistics
#' mummichog echoes back: m/z, retention time, p-value, statistic and the
#' `feature_id` we sent as the 5th input column. The per-feature adduct is
#' recovered from the EC's `str_row_ion` string (`"<row>_<ion>"` tokens joined
#' by `";"`, built by `EmpiricalCompound.__make_str_row_ion__()`), matched on the
#' feature's own input row — the file has no per-row ion column of its own.
#'
#' Every underlying feature is preserved: nothing is summed, averaged or
#' collapsed into a representative signal.
#'
#' @param files Character vector of mummichog output files.
#' @return A data.frame with columns `EID`, `input_row`, `feature_id`, `mz`,
#'   `retention_time`, `p_value`, `statistic`, `adduct`; zero rows when the file
#'   is absent.
read_mummichog_ec_features <- function(files) {
  empty <- data.frame(EID = character(0), input_row = character(0),
                      feature_id = character(0), mz = numeric(0),
                      retention_time = numeric(0), p_value = numeric(0),
                      statistic = numeric(0), adduct = character(0),
                      stringsAsFactors = FALSE)
  f <- .mmc_pick_result_file(files, "userInput_to_EmpiricalCompounds.tsv")
  if (is.null(f)) return(empty)

  raw <- tryCatch(
    readr::read_tsv(f, show_col_types = FALSE,
                    col_types = readr::cols(.default = readr::col_character()),
                    name_repair = "unique_quiet"),
    error = function(e) .mmc_stop(
      "could not read userInput_to_EmpiricalCompounds.tsv: ", conditionMessage(e))
  )
  miss <- setdiff(c("EID", "CompoundID_from_user"), names(raw))
  if (length(miss) > 0) {
    .mmc_stop("userInput_to_EmpiricalCompounds.tsv missing column(s): ",
              paste(miss, collapse = ", "), ". Present: ",
              paste(names(raw), collapse = ", "))
  }
  if (nrow(raw) == 0) return(empty)

  # The file repeats the `input_row` column (reporting.py writes it once itself
  # and once inside the feature's own output block), so readr renames the
  # duplicate; take the first, which is the EC-side key.
  row_col <- grep("^input_row", names(raw), value = TRUE)
  input_row <- if (length(row_col) > 0) as.character(raw[[row_col[[1]]]]) else
    rep(NA_character_, nrow(raw))

  num <- function(nm) {
    if (!nm %in% names(raw)) return(rep(NA_real_, nrow(raw)))
    suppressWarnings(as.numeric(raw[[nm]]))
  }
  ion_str <- if ("str_row_ion" %in% names(raw)) as.character(raw[["str_row_ion"]]) else
    rep(NA_character_, nrow(raw))

  data.frame(
    EID            = as.character(raw[["EID"]]),
    input_row      = input_row,
    feature_id     = as.character(raw[["CompoundID_from_user"]]),
    mz             = num("m/z"),
    retention_time = num("retention_time"),
    p_value        = num("p_value"),
    statistic      = num("statistic"),
    adduct         = .mmc_adduct_for_rows(input_row, ion_str),
    stringsAsFactors = FALSE
  )
}

#' Recover a feature's adduct from an EmpiricalCompound's str_row_ion string
#'
#' `str_row_ion` is `"<row>_<ion>"` tokens joined by `";"` (e.g.
#' `"row12_M+H[1+];row40_M+Na[1+]"`). The row id is delimited by the first `"_"`,
#' so the ion is everything after it for the token whose row matches.
#'
#' @param input_row Character vector of feature row ids.
#' @param ion_str   Character vector of the EC's `str_row_ion` (same length).
#' @return Character vector of adducts, `NA` when not recoverable.
#' @noRd
.mmc_adduct_for_rows <- function(input_row, ion_str) {
  vapply(seq_along(input_row), function(i) {
    row <- input_row[i]
    if (is.na(row) || !nzchar(row) || is.na(ion_str[i])) return(NA_character_)
    toks <- .mmc_split_cell(ion_str[i], ";")
    hit  <- toks[startsWith(toks, paste0(row, "_"))]
    if (length(hit) == 0) return(NA_character_)
    sub(paste0("^", row, "_"), "", hit[[1]])
  }, character(1), USE.NAMES = FALSE)
}


# ==== original annotation (schema-adaptive) =================================

# Column-name candidates, in preference order. Only columns that unambiguously
# mean "the compound this feature was annotated as" are considered; the feature
# id itself is never treated as an annotation, however name-like it looks.
.MMC_ANNOT_NAME_COLS <- c("Name", "Metabolite", "Molecule", "Compound",
                          "compound_name", "Compound_Name", "metabolite_name",
                          "Annotation", "annotation", "Identification")
.MMC_ANNOT_KEGG_COLS <- c("KEGG", "KEGG_ID", "KEGG ID", "kegg", "kegg_id",
                          "KEGG.ID")
.MMC_ANNOT_HMDB_COLS <- c("HMDB", "HMDB_ID", "HMDB ID", "hmdb_id", "HMDB.ID")
.MMC_ANNOT_CONF_COLS <- c("identification_level", "Identification_level",
                          "id_level", "MSI_level", "Confidence",
                          "confidence_level", "annotation_confidence",
                          "Annotation_confidence", "Level")

# Anchored, case-insensitive fallbacks for the same four meanings, so a dataset
# that merely differs in capitalisation or separator is still understood. They
# are deliberately anchored: a column is used only when its whole name says what
# it holds — never because it happens to contain a keyword.
.MMC_ANNOT_NAME_RX <- "^(name|metabolite|molecule|compound|compound[_. ]?name|metabolite[_. ]?name|annotation|identification)$"
.MMC_ANNOT_KEGG_RX <- "^kegg([_. ]?(id|compound))?$"
.MMC_ANNOT_HMDB_RX <- "^hmdb([_. ]?id)?$"
.MMC_ANNOT_CONF_RX <- "^((identification|id|msi|annotation)[_. ]?(level|confidence)|confidence([_. ]?level)?|level)$"

#' Normalise a dataset's own feature annotations into one generic contract
#'
#' Datasets differ: some carry MSI identification levels, some carry annotations
#' with no level system at all, and some features are simply unannotated. Rather
#' than hard-coding one schema (or "Level 1"), this locates whichever annotation
#' columns `row_data` actually has and flattens them into four fields. Absent
#' information stays `NA` — it is never inferred.
#'
#' Confidence and agreement are deliberately separate concepts: confidence is
#' what the dataset claims about its own annotation, agreement (see
#' `mmc_annotation_agreement()`) is whether that annotation and mummichog's
#' pathway-matching candidate refer to the same compound.
#'
#' @param row_data     Feature annotation table (`pre$row_data`).
#' @param mapping_file Optional HMDB -> KEGG mapping TSV; when given, an HMDB-only
#'   annotation gains a KEGG id through the pipeline's existing
#'   `read_hmdb_kegg_map()` so ID-based comparison stays possible.
#' @return A data.frame with one row per feature and columns `feature_id`,
#'   `original_annotation_name`, `original_annotation_id`,
#'   `original_annotation_id_type`, `original_annotation_kegg`,
#'   `original_annotation_confidence`. Zero rows when `row_data` is unusable.
normalize_metab_annotation <- function(row_data, mapping_file = NULL) {
  empty <- data.frame(feature_id = character(0),
                      original_annotation_name = character(0),
                      original_annotation_id = character(0),
                      original_annotation_id_type = character(0),
                      original_annotation_kegg = character(0),
                      original_annotation_confidence = character(0),
                      stringsAsFactors = FALSE)
  if (is.null(row_data) || !is.data.frame(row_data) || nrow(row_data) == 0) {
    return(empty)
  }
  if (!"feature_id" %in% names(row_data)) {
    if (is.null(rownames(row_data))) return(empty)
    row_data$feature_id <- rownames(row_data)
  }
  n <- nrow(row_data)

  name_col <- .mmc_find_col(row_data, .MMC_ANNOT_NAME_COLS, .MMC_ANNOT_NAME_RX)
  kegg_col <- .mmc_find_col(row_data, .MMC_ANNOT_KEGG_COLS, .MMC_ANNOT_KEGG_RX)
  hmdb_col <- .mmc_find_col(row_data, .MMC_ANNOT_HMDB_COLS, .MMC_ANNOT_HMDB_RX)
  conf_col <- .mmc_find_col(row_data, .MMC_ANNOT_CONF_COLS, .MMC_ANNOT_CONF_RX)

  nm <- if (!is.null(name_col)) trimws(as.character(row_data[[name_col]])) else
    rep(NA_character_, n)
  nm[!is.na(nm) & (!nzchar(nm) | nm %in% c("NA", "-", "unknown", "Unknown"))] <-
    NA_character_

  kegg <- if (!is.null(kegg_col)) .mmc_extract_kegg(row_data[[kegg_col]]) else
    rep(NA_character_, n)
  hmdb <- if (!is.null(hmdb_col)) .mmc_extract_hmdb(row_data[[hmdb_col]]) else
    rep(NA_character_, n)

  # HMDB -> KEGG through the pipeline's existing mapping reader, so an HMDB-only
  # dataset can still be compared on a stable compound id.
  map_vec <- read_hmdb_kegg_map(mapping_file)
  if (!is.null(map_vec)) {
    need <- is.na(kegg) & !is.na(hmdb)
    if (any(need)) {
      got <- .mmc_extract_kegg(map_vec[hmdb[need]])
      kegg[need] <- got
    }
  }

  # Canonical id: a stable compound id, KEGG preferred over HMDB.
  id      <- ifelse(!is.na(kegg), kegg, hmdb)
  id_type <- ifelse(!is.na(kegg), "KEGG", ifelse(!is.na(hmdb), "HMDB",
                                                 NA_character_))

  conf <- if (!is.null(conf_col)) .mmc_format_confidence(row_data[[conf_col]]) else
    rep(NA_character_, n)

  data.frame(
    feature_id                     = as.character(row_data$feature_id),
    original_annotation_name       = nm,
    original_annotation_id         = id,
    original_annotation_id_type    = id_type,
    original_annotation_kegg       = kegg,
    original_annotation_confidence = conf,
    stringsAsFactors = FALSE
  )
}

#' Render an annotation-confidence column as human-readable text
#'
#' A bare numeric column (the pipeline's `identification_level`) becomes
#' `"Level 1"`, `"Level 2"`, ...; anything already textual is passed through
#' verbatim, because a dataset that spells its own confidence scheme knows it
#' better than we do. Datasets with no level system keep `NA`.
#'
#' @param x The raw confidence column.
#' @return Character vector, `NA` where the dataset says nothing.
#' @noRd
.mmc_format_confidence <- function(x) {
  chr <- trimws(as.character(x))
  chr[is.na(x) | !nzchar(chr) | chr == "NA"] <- NA_character_
  num <- suppressWarnings(as.numeric(chr))
  is_num <- !is.na(num) & grepl("^[0-9.]+$", chr)
  out <- chr
  out[is_num] <- paste("Level", format(num[is_num], trim = TRUE,
                                       drop0trailing = TRUE))
  out
}

#' Decide whether an original annotation and pathway-matching candidates agree
#'
#' Three outcomes, kept strictly separate from annotation confidence:
#'
#' \describe{
#'   \item{`"Match"`}{The original annotation and at least one pathway-matching
#'     candidate refer to the same compound.}
#'   \item{`"Conflict"`}{There IS a comparable original annotation, but it
#'     represents a different metabolite.}
#'   \item{`"Not assessed"`}{No usable original annotation to compare against
#'     (unannotated feature, or nothing comparable on either side).}
#' }
#'
#' Comparison prefers stable compound ids — KEGG when both sides have one — and
#' only falls back to a conservatively normalised name/synonym comparison when
#' ids are unavailable. Model compound names are `";"`-separated synonym lists
#' (mummichog's own `dict_cpds_def` convention), and every synonym counts.
#' Identity is NEVER inferred from m/z, molecular formula or mass.
#'
#' A conflict is an annotation, not a veto: nothing here removes an
#' EmpiricalCompound or a pathway from the mummichog result.
#'
#' @param annot_kegg     The feature's original KEGG id (or NA).
#' @param annot_name     The feature's original annotation name (or NA).
#' @param candidate_ids  Character vector of pathway-matching candidate compound ids.
#' @param candidate_kegg Character vector of those candidates' KEGG ids (NA allowed).
#' @param candidate_names Character vector of those candidates' names (may hold
#'   `";"`-separated synonyms).
#' @return One of `"Match"`, `"Conflict"`, `"Not assessed"`.
mmc_annotation_agreement <- function(annot_kegg, annot_name,
                                     candidate_ids, candidate_kegg,
                                     candidate_names) {
  has_annot <- (!is.na(annot_kegg) && nzchar(annot_kegg)) ||
    (!is.na(annot_name) && nzchar(annot_name))
  if (!has_annot || length(candidate_ids) == 0) return("Not assessed")

  # --- 1. stable compound ids (preferred) ---------------------------------
  if (!is.na(annot_kegg) && nzchar(annot_kegg)) {
    # A candidate's own id can itself be a KEGG accession (KEGG-based models),
    # so consider both the declared KEGG id and the id.
    cand_k <- unique(c(candidate_kegg, .mmc_extract_kegg(candidate_ids)))
    cand_k <- cand_k[!is.na(cand_k) & nzchar(cand_k)]
    if (length(cand_k) > 0) {
      return(if (annot_kegg %in% cand_k) "Match" else "Conflict")
    }
  }

  # --- 2. conservative name / synonym comparison --------------------------
  if (!is.na(annot_name) && nzchar(annot_name)) {
    syn <- unique(unlist(lapply(candidate_names, function(x) {
      if (is.na(x)) character(0) else .mmc_split_cell(x, ";")
    }), use.names = FALSE))
    syn <- .mmc_norm_name(syn)
    syn <- syn[nzchar(syn)]
    if (length(syn) > 0) {
      key <- .mmc_norm_name(annot_name)
      if (!nzchar(key)) return("Not assessed")
      return(if (key %in% syn) "Match" else "Conflict")
    }
  }

  # Nothing comparable on the candidate side (no ids, no names).
  "Not assessed"
}


# ==== the evidence layer ====================================================

#' Trace the supporting evidence for every pathway in a mummichog result
#'
#' Builds the `Pathway -> EmpiricalCompound -> pathway-matching candidate ->
#' measured feature -> original annotation -> agreement` chain as three
#' report-ready frames. Pure: it reads mummichog's tables (already on disk), the
#' metabolic model index and the normalised annotations, and computes no
#' statistics of its own — the ORA p-values and overlaps are passed through
#' untouched, conflicts included.
#'
#' Pathway-matching candidates are derived as
#' `intersect(all candidates of the EC, compounds of THIS pathway)` and every
#' surviving candidate is kept; when one EmpiricalCompound contributes two
#' candidates that both belong to the pathway, both are listed and the EC is
#' still counted once.
#'
#' @param pathways      One contrast's mummichog pathway table (from
#'   `read_mummichog_pathways()`).
#' @param files         That contrast's mummichog output files.
#' @param model         Model index from `read_mummichog_model_pathways()` /
#'   `mmc_load_model_pathways()`.
#' @param annot         Normalised annotations from `normalize_metab_annotation()`
#'   (may be zero-row; every feature then reads `"Not assessed"`).
#' @param p_cutoff      Feature significance cutoff — the same
#'   `mummichog.p_cutoff` the engine was run with — used only to flag features as
#'   significant, never to filter them.
#' @param ec_col        Name of the pathway table's EmpiricalCompound column.
#' @return `NULL` when evidence cannot be built at all (no model, no EC tables,
#'   empty pathway table). Otherwise a list with:
#'   \describe{
#'     \item{pathway_summary}{One row per pathway: overlap, detected pathway
#'       size, enrichment ratio, empirical p-value, supporting EC/feature counts
#'       and the Match / Conflict / Not assessed breakdown.}
#'     \item{ec_table}{One row per (pathway, supporting EmpiricalCompound).}
#'     \item{feature_table}{One row per (pathway, EmpiricalCompound, measured
#'       feature) — every underlying signal, nothing collapsed.}
#'   }
build_mummichog_pathway_evidence <- function(pathways, files, model, annot,
                                             p_cutoff = 0.05,
                                             ec_col = "overlap_EmpiricalCompounds (id)") {
  if (is.null(pathways) || !is.data.frame(pathways) || nrow(pathways) == 0) {
    return(NULL)
  }
  if (is.null(model) || length(model$pathways) == 0) return(NULL)
  if (!ec_col %in% names(pathways)) {
    message("mummichog evidence: pathway table has no '", ec_col,
            "' column — skipping the evidence layer")
    return(NULL)
  }

  cand <- read_mummichog_ec_candidates(files)
  feat <- read_mummichog_ec_features(files)
  if (nrow(cand) == 0 || nrow(feat) == 0) {
    message("mummichog evidence: no EmpiricalCompound tables among the ",
            "mummichog outputs — skipping the evidence layer")
    return(NULL)
  }

  p_col <- .mmc_find_col(pathways,
                         c("p-value", "p.value", "pvalue", "p_value", "P.Value"),
                         "^p[._-]?value$")

  cand_by_ec <- split(cand, cand$EID)
  feat_by_ec <- split(feat, feat$EID)
  annot_idx  <- if (nrow(annot) > 0) {
    stats::setNames(seq_len(nrow(annot)), as.character(annot$feature_id))
  } else {
    integer(0)
  }

  ec_rows      <- list()
  feature_rows <- list()
  summary_rows <- list()

  for (i in seq_len(nrow(pathways))) {
    pw_name <- as.character(pathways$pathway[i])
    pw_cpds <- model$pathways[[pw_name]]
    eids    <- .mmc_split_cell(pathways[[ec_col]][i], ",")
    if (length(eids) == 0) next

    ec_acc   <- list()
    feat_acc <- list()
    for (eid in eids) {
      ec_cand <- cand_by_ec[[eid]]
      if (is.null(ec_cand)) next
      # Pathway-matching candidate(s): the EC's candidates that belong to THIS
      # pathway. All of them are kept — never one arbitrary pick.
      keep <- if (is.null(pw_cpds)) logical(0) else
        ec_cand$compound_id %in% pw_cpds
      match_cand <- ec_cand[keep, , drop = FALSE]
      if (nrow(match_cand) == 0) next

      cand_ids   <- match_cand$compound_id
      cand_nms   <- ifelse(is.na(match_cand$compound_name),
                           unname(model$compound_names[cand_ids]),
                           match_cand$compound_name)
      cand_kegg  <- unname(model$compound_kegg[cand_ids])
      if (length(cand_kegg) == 0) cand_kegg <- rep(NA_character_, length(cand_ids))

      ec_feat <- feat_by_ec[[eid]]
      if (is.null(ec_feat) || nrow(ec_feat) == 0) next

      # Per-feature original annotation + agreement. Every feature of the EC is
      # kept as its own row: signals are never summed or reduced to one
      # representative.
      idx <- annot_idx[ec_feat$feature_id]
      a_name <- rep(NA_character_, nrow(ec_feat))
      a_id   <- rep(NA_character_, nrow(ec_feat))
      a_kegg <- rep(NA_character_, nrow(ec_feat))
      a_conf <- rep(NA_character_, nrow(ec_feat))
      ok <- !is.na(idx)
      if (any(ok)) {
        a_name[ok] <- annot$original_annotation_name[idx[ok]]
        a_id[ok]   <- annot$original_annotation_id[idx[ok]]
        a_kegg[ok] <- annot$original_annotation_kegg[idx[ok]]
        a_conf[ok] <- annot$original_annotation_confidence[idx[ok]]
      }
      agree <- vapply(seq_len(nrow(ec_feat)), function(k) {
        mmc_annotation_agreement(a_kegg[k], a_name[k],
                                 cand_ids, cand_kegg, cand_nms)
      }, character(1))

      feat_acc[[length(feat_acc) + 1L]] <- data.frame(
        check.names = FALSE, stringsAsFactors = FALSE,
        "Pathway"                = pw_name,
        "EmpiricalCompound"      = eid,
        "Feature"                = ec_feat$feature_id,
        "m/z"                    = ec_feat$mz,
        "RT"                     = ec_feat$retention_time,
        "Adduct/ion"             = ec_feat$adduct,
        "Feature p-value"        = ec_feat$p_value,
        "Feature statistic"      = ec_feat$statistic,
        "Significant"            = ifelse(is.na(ec_feat$p_value), NA,
                                          ec_feat$p_value < p_cutoff),
        "Original annotation"    = a_name,
        "Annotation ID"          = a_id,
        "Annotation confidence"  = a_conf,
        "Agreement"              = agree
      )

      ec_acc[[length(ec_acc) + 1L]] <- data.frame(
        check.names = FALSE, stringsAsFactors = FALSE,
        "Pathway"                        = pw_name,
        "EmpiricalCompound"              = eid,
        "# Features"                     = nrow(ec_feat),
        "Pathway-matching candidate(s)"  = .mmc_join_unique(cand_ids),
        "Candidate name(s)"              = .mmc_join_unique(
                                             vapply(cand_nms, function(x)
                                               if (is.na(x)) NA_character_ else
                                                 .mmc_split_cell(x, ";")[1],
                                               character(1), USE.NAMES = FALSE)),
        "Candidate KEGG ID(s)"           = .mmc_join_unique(cand_kegg),
        "Original annotation"            = .mmc_join_unique(a_name),
        "Annotation confidence"          = .mmc_join_unique(a_conf),
        "Agreement"                      = .mmc_summarise_agreement(agree)
      )
    }

    if (length(ec_acc) == 0) next
    ec_df   <- do.call(rbind, ec_acc)
    feat_df <- do.call(rbind, feat_acc)
    ec_rows[[length(ec_rows) + 1L]]           <- ec_df
    feature_rows[[length(feature_rows) + 1L]] <- feat_df

    overlap <- suppressWarnings(as.numeric(pathways$overlap_size[i]))
    pw_size <- suppressWarnings(as.numeric(pathways$pathway_size[i]))
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      check.names = FALSE, stringsAsFactors = FALSE,
      "Pathway"                 = pw_name,
      "Overlap"                 = overlap,
      "Detected pathway size"   = pw_size,
      "Enrichment ratio"        = round(overlap / pw_size, 3),
      "p.value"                 = if (is.null(p_col)) NA_real_ else
                                    suppressWarnings(as.numeric(pathways[[p_col]][i])),
      "Supporting ECs"          = nrow(ec_df),
      "Supporting features"     = nrow(feat_df),
      "Match"                   = sum(ec_df$Agreement == "Match"),
      "Conflict"                = sum(ec_df$Agreement == "Conflict"),
      "Not assessed"            = sum(ec_df$Agreement == "Not assessed")
    )
  }

  if (length(ec_rows) == 0) return(NULL)
  summary_df <- do.call(rbind, summary_rows)
  summary_df <- summary_df[order(summary_df[["p.value"]],
                                 na.last = TRUE), , drop = FALSE]
  list(
    pathway_summary = summary_df,
    ec_table        = do.call(rbind, ec_rows),
    feature_table   = do.call(rbind, feature_rows)
  )
}

#' Join unique, non-missing values into one display cell
#' @param x   Character vector.
#' @param sep Separator for the joined output.
#' @return A single character scalar, `NA` when nothing usable remains.
#' @noRd
.mmc_join_unique <- function(x, sep = "; ") {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) NA_character_ else paste(x, collapse = sep)
}

#' Roll per-feature agreements up to one EmpiricalCompound verdict
#'
#' An EmpiricalCompound is supported by several measured signals, each with its
#' own original annotation. Precedence is `Match > Conflict > Not assessed`: if
#' ANY underlying feature was originally annotated as a pathway-matching
#' candidate, the identity mummichog used does agree with our annotation for
#' that EC. The per-feature detail stays visible in `feature_table`.
#'
#' @param x Character vector of per-feature verdicts.
#' @return One of `"Match"`, `"Conflict"`, `"Not assessed"`.
#' @noRd
.mmc_summarise_agreement <- function(x) {
  if (any(x == "Match"))    return("Match")
  if (any(x == "Conflict")) return("Conflict")
  "Not assessed"
}
