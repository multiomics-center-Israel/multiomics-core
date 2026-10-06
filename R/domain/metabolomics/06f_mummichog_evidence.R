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
# the evidence layer reports itself unavailable (returns NULL) rather than
# guessing. Nothing in this file is wired into the pipeline or the report yet.

# Statistical caveat: this layer reports the per-feature p-value and logFC that
# mummichog received as input, but does NOT flag features as significant. That
# cutoff is applied by the mummichog stage itself, and re-deriving it here would
# be a second interpretation of a decision that belongs to the stage.
#
# Dependencies (R): jsonlite, readr. Both already in renv.lock.


# ==== internal helpers ======================================================

#' Pick one named result file out of a mummichog output file list
#'
#' The 06c readers each locate their own file inline and stop when it is
#' missing. The readers here use this lookup instead, to read a file 06c does not
#' read (`ListOfEmpiricalCompounds.tsv`) and to return zero rows, not an error,
#' when an EC table is absent.
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

#' Split a delimited cell keeping every slot, for positional pairing
#'
#' Unlike `.mmc_split_cell()`, empty fields are kept (as `NA`), so element `i`
#' still lines up with element `i` of a parallel list. mummichog 2.7.0 writes an
#' empty name for a candidate the model has no name for
#' (`reporting.py`: `dict_cpds_def.get(x, '')` joined with `"$"`), so `"$Glucose"`
#' means "no name, then Glucose" — dropping the empty field would hand
#' "Glucose" to the first candidate instead of the second.
#'
#' @param x   A single character scalar (may be NA).
#' @param sep Fixed separator.
#' @return Character vector with one element per field, `NA` for empty fields;
#'   `character(0)` for a missing or empty cell.
#' @noRd
.mmc_split_positional <- function(x, sep) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return(character(0))
  x   <- as.character(x)
  out <- trimws(strsplit(x, sep, fixed = TRUE)[[1]])
  # strsplit() drops a trailing empty field ("a$" -> "a"); one field per
  # separator plus one restores it.
  n_fields <- lengths(regmatches(x, gregexpr(sep, x, fixed = TRUE))) + 1L
  length(out) <- n_fields
  out[!is.na(out) & !nzchar(out)] <- NA_character_
  out
}

#' Normalise HMDB accessions to the 7-digit form
#'
#' HMDB ids appear both as the legacy 5-digit form ("HMDB00122") and the current
#' 7-digit form ("HMDB0000122") for the same compound; compare them only after
#' zero-padding to 7 digits.
#'
#' @param x Character vector of HMDB accessions (already extracted).
#' @return Character vector, `NA` where `x` is not an HMDB accession.
#' @noRd
.mmc_norm_hmdb <- function(x) {
  x   <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  ok  <- !is.na(x) & grepl("^HMDB[0-9]+$", x)
  out[ok] <- sprintf("HMDB%07d", as.integer(sub("^HMDB", "", x[ok])))
  out
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
#' have the model available when evidence is built. Callers can then omit the
#' evidence instead of showing guessed identities.
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
#'   `str_row_ion`, `massfeature_rows`; zero rows when the file is absent, so a
#'   caller can treat it as "no evidence" rather than an error.
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
    # Positional pairing is mummichog's own contract, so both lists keep their
    # empty slots until they are paired; pad rather than recycle so a truncated
    # name list can never mislabel a candidate.
    ids <- .mmc_split_positional(raw[["compounds"]][i], ";")
    nms <- .mmc_split_positional(names_col[i], "$")
    if (length(nms) < length(ids)) {
      nms <- c(nms, rep(NA_character_, length(ids) - length(nms)))
    }
    keep <- !is.na(ids)
    if (!any(keep)) return(NULL)
    data.frame(
      EID              = as.character(raw[["EID"]][i]),
      compound_id      = ids[keep],
      compound_name    = nms[seq_along(ids)][keep],
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
#' Shapes `read_mummichog_empirical_map()` (06c) — the one reader of
#' `userInput_to_EmpiricalCompounds.tsv` — into one row per (EmpiricalCompound,
#' input feature) with numeric m/z, retention time, p-value and statistic, and
#' the `feature_id` we sent as the 5th input column. These are the values
#' mummichog echoes back from its input, so `statistic` is the logFC the stage
#' writes there (05b), not a moderated t. The per-feature adduct is recovered
#' from the EC's `str_row_ion` string (`"<row>_<ion>"` tokens joined by `";"`,
#' built by `EmpiricalCompound.__make_str_row_ion__()`), matched on the feature's
#' own input row — the file has no per-row ion column of its own.
#'
#' Every underlying feature is preserved: nothing is summed, averaged or
#' collapsed into a representative signal.
#'
#' Unlike `read_mummichog_empirical_map()`, which stops when the file is absent,
#' this returns zero rows then, so a caller can treat "no EC table" as "no
#' evidence" without catching an error.
#'
#' @param files Character vector of mummichog output files.
#' @return A data.frame with columns `EID`, `input_row`, `feature_id`, `mz`,
#'   `retention_time`, `p_value`, `statistic`, `adduct`; zero rows when the file
#'   is absent or empty.
read_mummichog_ec_features <- function(files) {
  empty <- data.frame(EID = character(0), input_row = character(0),
                      feature_id = character(0), mz = numeric(0),
                      retention_time = numeric(0), p_value = numeric(0),
                      statistic = numeric(0), adduct = character(0),
                      stringsAsFactors = FALSE)
  if (is.null(.mmc_pick_result_file(files, "userInput_to_EmpiricalCompounds.tsv"))) {
    return(empty)
  }

  emap <- read_mummichog_empirical_map(files)
  if (nrow(emap) == 0) return(empty)

  num <- function(x) suppressWarnings(as.numeric(x))
  data.frame(
    EID            = as.character(emap$EID),
    input_row      = emap$input_row,
    feature_id     = as.character(emap$feature_id),
    mz             = num(emap$mz),
    retention_time = num(emap$retention_time),
    p_value        = num(emap$p_value),
    statistic      = num(emap$statistic),
    adduct         = .mmc_adduct_for_rows(emap$input_row, emap$str_row_ion),
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
#' columns `row_data` actually has and flattens them into one set of fields.
#' Absent information stays `NA` — it is never inferred.
#'
#' Several recognised columns can hold the same field: a multi-level dataset is
#' the union of its levels' columns, so one level may fill `Name` and another
#' `Molecule`, each `NA` on the other's rows. Every recognised column is
#' therefore read and the first usable value is taken row by row, in the
#' preference order of the column lists below.
#'
#' Feature ids are resolved with `mmc_feature_ids()` (06c), the same rule the
#' mummichog stage uses for the ids it sends, so the annotations join back to
#' the features mummichog echoes.
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
#'   `original_annotation_hmdb`, `original_annotation_confidence`. Zero rows
#'   when `row_data` is unusable.
normalize_metab_annotation <- function(row_data, mapping_file = NULL) {
  empty <- data.frame(feature_id = character(0),
                      original_annotation_name = character(0),
                      original_annotation_id = character(0),
                      original_annotation_id_type = character(0),
                      original_annotation_kegg = character(0),
                      original_annotation_hmdb = character(0),
                      original_annotation_confidence = character(0),
                      stringsAsFactors = FALSE)
  if (is.null(row_data) || !is.data.frame(row_data) || nrow(row_data) == 0) {
    return(empty)
  }
  if (!"feature_id" %in% names(row_data)) {
    ids <- mmc_feature_ids(row_data)
    if (is.null(ids)) return(empty)
    row_data$feature_id <- ids
  }

  clean_name <- function(x) {
    v <- trimws(as.character(x))
    v[!is.na(v) & (!nzchar(v) | v %in% c("NA", "-", "unknown", "Unknown"))] <-
      NA_character_
    v
  }
  nm <- .mmc_coalesce_annot(row_data, .MMC_ANNOT_NAME_COLS, .MMC_ANNOT_NAME_RX,
                            clean_name)
  kegg <- .mmc_coalesce_annot(row_data, .MMC_ANNOT_KEGG_COLS, .MMC_ANNOT_KEGG_RX,
                              .mmc_extract_kegg)
  hmdb <- .mmc_coalesce_annot(row_data, .MMC_ANNOT_HMDB_COLS, .MMC_ANNOT_HMDB_RX,
                              .mmc_extract_hmdb)

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

  conf <- .mmc_coalesce_annot(row_data, .MMC_ANNOT_CONF_COLS, .MMC_ANNOT_CONF_RX,
                              .mmc_format_confidence)

  data.frame(
    feature_id                     = as.character(row_data$feature_id),
    original_annotation_name       = nm,
    original_annotation_id         = id,
    original_annotation_id_type    = id_type,
    original_annotation_kegg       = kegg,
    original_annotation_hmdb       = hmdb,
    original_annotation_confidence = conf,
    stringsAsFactors = FALSE
  )
}

#' First usable value of one annotation field, row by row across its columns
#'
#' Collects every column of `df` that means the field — the exact names first,
#' in list order, then the anchored case-insensitive regex matches — and fills
#' each row from the first column whose transformed value is not `NA`. The
#' feature id column is never treated as an annotation.
#'
#' @param df        Feature annotation table.
#' @param exact     Character vector of exact column names, in preference order.
#' @param rx        Anchored regex for the same field, matched case-insensitively.
#' @param transform Function turning one raw column into a character vector,
#'   `NA` where the cell holds nothing usable.
#' @return Character vector with one value per row of `df`.
#' @noRd
.mmc_coalesce_annot <- function(df, exact, rx, transform) {
  cols <- setdiff(names(df), "feature_id")
  cols <- unique(c(exact[exact %in% cols],
                   cols[grepl(rx, cols, ignore.case = TRUE)]))
  out <- rep(NA_character_, nrow(df))
  for (cl in cols) {
    v    <- as.character(transform(df[[cl]]))
    fill <- is.na(out) & !is.na(v)
    out[fill] <- v[fill]
  }
  out
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
#'   \item{`"Conflict"`}{At least one candidate could be compared with the
#'     original annotation, and none of them is the same metabolite.}
#'   \item{`"Not assessed"`}{No usable original annotation, or no candidate
#'     that can be compared with it.}
#' }
#'
#' Each candidate is compared on its own, with stable ids taking precedence
#' over names for that candidate: KEGG when both sides have one, else HMDB when
#' both sides have one (5- and 7-digit forms compared after zero-padding), else
#' a conservatively normalised name/synonym comparison. Model compound names are
#' `";"`-separated synonym lists (mummichog's own `dict_cpds_def` convention),
#' and every synonym counts. So a candidate whose KEGG id differs does not stop
#' another candidate, with no KEGG id, from matching by name.
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
#' @param annot_hmdb     The feature's original HMDB id (or NA).
#' @return One of `"Match"`, `"Conflict"`, `"Not assessed"`.
mmc_annotation_agreement <- function(annot_kegg, annot_name,
                                     candidate_ids, candidate_kegg,
                                     candidate_names,
                                     annot_hmdb = NA_character_) {
  usable <- function(x) length(x) == 1 && !is.na(x) && nzchar(x)
  a_kegg <- if (usable(annot_kegg)) annot_kegg else NA_character_
  a_hmdb <- if (usable(annot_hmdb)) .mmc_norm_hmdb(annot_hmdb) else NA_character_
  a_key  <- if (usable(annot_name)) .mmc_norm_name(annot_name) else ""
  if (is.na(a_kegg) && is.na(a_hmdb) && !nzchar(a_key)) return("Not assessed")
  if (length(candidate_ids) == 0) return("Not assessed")

  verdict_for <- function(j) {
    id <- candidate_ids[j]
    # --- 1. KEGG ids (a candidate's own id can itself be a KEGG accession) ---
    ck <- c(candidate_kegg[j], .mmc_extract_kegg(id))
    ck <- unique(ck[!is.na(ck) & nzchar(ck)])
    if (!is.na(a_kegg) && length(ck) > 0) {
      return(if (a_kegg %in% ck) "Match" else "Conflict")
    }
    # --- 2. HMDB ids (custom HMDB-based models) -----------------------------
    ch <- .mmc_norm_hmdb(.mmc_extract_hmdb(id))
    if (!is.na(a_hmdb) && !is.na(ch)) {
      return(if (identical(a_hmdb, ch)) "Match" else "Conflict")
    }
    # --- 3. conservative name / synonym comparison --------------------------
    if (nzchar(a_key)) {
      nm  <- candidate_names[j]
      syn <- if (is.na(nm)) character(0) else .mmc_norm_name(.mmc_split_cell(nm, ";"))
      syn <- syn[nzchar(syn)]
      if (length(syn) > 0) return(if (a_key %in% syn) "Match" else "Conflict")
    }
    "Not assessed"
  }
  verdicts <- vapply(seq_along(candidate_ids), verdict_for, character(1))

  if (any(verdicts == "Match"))    return("Match")
  if (any(verdicts == "Conflict")) return("Conflict")
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
#' @param ec_col        Name of the pathway table's EmpiricalCompound column.
#' @return `NULL` when evidence cannot be built at all (no model, no EC tables,
#'   empty pathway table). Otherwise a list with:
#'   \describe{
#'     \item{pathway_summary}{One row per pathway: overlap, detected pathway
#'       size, enrichment ratio, empirical p-value, `Supporting ECs`,
#'       `Supporting features` (distinct measured features), `Feature-EC links`
#'       (rows of `feature_table`; one feature can sit in several ECs), and the
#'       agreement breakdown at both grains — `ECs Match/Conflict/Mixed/Not
#'       assessed` counts EmpiricalCompounds by their roll-up state, while
#'       `feature-EC links Match/Conflict/Not assessed` counts feature-EC links
#'       by their own verdict and sums to `Feature-EC links`.}
#'     \item{ec_table}{One row per (pathway, supporting EmpiricalCompound), with
#'       the four-state `Agreement` roll-up and the feature-level `n_match`,
#'       `n_conflict`, `n_not_assessed` counts behind it.}
#'     \item{feature_table}{One row per (pathway, EmpiricalCompound, measured
#'       feature) — every underlying signal, nothing collapsed. Its p-value and
#'       `Feature log2FC (mummichog input)` are the values mummichog received as
#'       input (raw DE p-value and logFC). The table does not mark features as
#'       significant: that cutoff is applied by the mummichog stage and is not
#'       re-derived here.}
#'   }
build_mummichog_pathway_evidence <- function(pathways, files, model, annot,
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
      a_hmdb <- rep(NA_character_, nrow(ec_feat))
      a_conf <- rep(NA_character_, nrow(ec_feat))
      ok <- !is.na(idx)
      if (any(ok)) {
        a_name[ok] <- annot$original_annotation_name[idx[ok]]
        a_id[ok]   <- annot$original_annotation_id[idx[ok]]
        a_kegg[ok] <- annot$original_annotation_kegg[idx[ok]]
        if (!is.null(annot$original_annotation_hmdb)) {
          a_hmdb[ok] <- annot$original_annotation_hmdb[idx[ok]]
        }
        a_conf[ok] <- annot$original_annotation_confidence[idx[ok]]
      }
      agree <- vapply(seq_len(nrow(ec_feat)), function(k) {
        mmc_annotation_agreement(a_kegg[k], a_name[k],
                                 cand_ids, cand_kegg, cand_nms,
                                 annot_hmdb = a_hmdb[k])
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
        "Feature log2FC (mummichog input)" = ec_feat$statistic,
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
        "Agreement"                      = .mmc_summarise_agreement(agree),
        # FEATURE-level counts of the evidence behind this EC's roll-up. They
        # count measured features, not ECs, candidates or pathway members —
        # n_match + n_conflict + n_not_assessed == `# Features`.
        "n_match"                        = sum(agree == "Match"),
        "n_conflict"                     = sum(agree == "Conflict"),
        "n_not_assessed"                 = sum(agree == "Not assessed")
      )
    }

    if (length(ec_acc) == 0) next
    ec_df   <- do.call(rbind, ec_acc)
    feat_df <- do.call(rbind, feat_acc)
    ec_rows[[length(ec_rows) + 1L]]           <- ec_df
    feature_rows[[length(feature_rows) + 1L]] <- feat_df

    overlap <- suppressWarnings(as.numeric(pathways$overlap_size[i]))
    pw_size <- suppressWarnings(as.numeric(pathways$pathway_size[i]))
    # Two grains, never mixed: "ECs ..." columns count EmpiricalCompounds by
    # their roll-up state, "feature-EC links ..." columns count rows of the
    # feature table by their own verdict. One measured feature can sit in
    # several EmpiricalCompounds (mummichog writes one row per feature-EC pair),
    # so links are not features: "Supporting features" counts distinct feature
    # ids, and a feature's verdict can differ between its ECs. Column names
    # carry the grain so a reader cannot mistake one for another, and none is
    # the pathway overlap or the candidate count.
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      check.names = FALSE, stringsAsFactors = FALSE,
      "Pathway"                 = pw_name,
      "Overlap"                 = overlap,
      "Detected pathway size"   = pw_size,
      "Enrichment ratio"        = round(overlap / pw_size, 3),
      "p.value"                 = if (is.null(p_col)) NA_real_ else
                                    suppressWarnings(as.numeric(pathways[[p_col]][i])),
      "Supporting ECs"          = nrow(ec_df),
      "Supporting features"     = length(unique(feat_df$Feature)),
      "Feature-EC links"        = nrow(feat_df),
      "ECs Match"               = sum(ec_df$Agreement == "Match"),
      "ECs Conflict"            = sum(ec_df$Agreement == "Conflict"),
      "ECs Mixed"               = sum(ec_df$Agreement == "Mixed"),
      "ECs Not assessed"        = sum(ec_df$Agreement == "Not assessed"),
      "feature-EC links Match"        = sum(feat_df$Agreement == "Match"),
      "feature-EC links Conflict"     = sum(feat_df$Agreement == "Conflict"),
      "feature-EC links Not assessed" = sum(feat_df$Agreement == "Not assessed")
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
#' own original annotation, so its roll-up needs four states rather than a
#' precedence chain — an EC where one feature agrees and another disagrees is
#' genuinely `Mixed`, and collapsing that to `Match` would hide the
#' disagreement:
#'
#' \preformatted{
#' has_match && has_conflict   -> "Mixed"
#' has_match && !has_conflict  -> "Match"
#' !has_match && has_conflict  -> "Conflict"
#' otherwise                   -> "Not assessed"
#' }
#'
#' Features with no usable annotation never override assessed evidence, so
#' `Match + Not assessed` is `Match` and `Conflict + Not assessed` is
#' `Conflict`. Per-feature verdicts stay unchanged and visible in
#' `feature_table`.
#'
#' @param x Character vector of per-feature verdicts.
#' @return One of `"Match"`, `"Conflict"`, `"Mixed"`, `"Not assessed"`.
#' @noRd
.mmc_summarise_agreement <- function(x) {
  has_match    <- any(x == "Match")
  has_conflict <- any(x == "Conflict")
  if (has_match && has_conflict)  return("Mixed")
  if (has_match)                  return("Match")
  if (has_conflict)               return("Conflict")
  "Not assessed"
}

# The four EC-level agreement states, in report order.
.MMC_AGREEMENT_STATES <- c("Match", "Conflict", "Mixed", "Not assessed")
