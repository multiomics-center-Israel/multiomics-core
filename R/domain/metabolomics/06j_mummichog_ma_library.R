# R/domain/metabolomics/06j_mummichog_ma_library.R
#
# MetaboAnalyst-compatible mummichog, component 1: the compound library.
#
# The native MetaboAnalyst-compatible engine (06j-06l) matches m/z values
# against a library in MetaboAnalystR's own `mummichog.lib` shape. MetaboAnalystR
# only ever loads that shape from a pre-built `.qs` file; it has no reader for a
# mummichog/metDataModel JSON model. This file builds the same object from our
# model JSON (e.g. mng_kegg_20260723), so the engine can run on our
# organism-specific models.
#
# ===========================================================================
# BEHAVIOURAL REFERENCE — pinned upstream version
# ===========================================================================
#     repo:   xia-lab/MetaboAnalystR
#     commit: 398476ae2a0c996e390925fa62adc46b5ecde334
#     files:  R/util_psea.R          CreateLibFromKEGG(), Compound_function_mzlist()
#             R/peaks_to_function.R  .setup.psea.library() (how the object is read)
#
# Reproduced verbatim:
#   * adduct m/z per compound — Compound_function_mzlist() for "positive"
#     (25 forms) and "negative" (16 forms): same constants, same column order,
#     same column names;
#   * the m/z lookup tree — CreateLibFromKEGG(): integer-m/z bins 51..1999,
#     each holding the cpd.lib row indices with any adduct m/z flooring into it.
#
# Where our library necessarily differs from a MetaboAnalyst KEGG library:
#   * compound masses are the model's `neutral_mono_mass`, not MetaboAnalyst's
#     KEGG compound table;
#   * pathway membership comes from the model's reactions (reactants + products
#     of every reaction in the pathway), exactly as mummichog derives it from
#     the same model (read_mummichog_model_pathways()), not from KEGG's own
#     pathway-compound lists. A pathway therefore holds the compounds our python
#     engine sees for it, which is what makes the two engines comparable on one
#     model;
#   * the "dpj_positive" adduct set is not built: no analysis mode the engine
#     supports reads it.


# Proton mass used by Compound_function_mzlist() (util_psea.R).
.MMC_MA_PROTON <- 1.00727646677

#' Adduct m/z matrix for a vector of neutral masses (MetaboAnalystR forms)
#'
#' A line-for-line port of MetaboAnalystR's `Compound_function_mzlist()` for the
#' two modes a MetaboAnalyst mummichog library is read in: `"positive"` (25 ion
#' forms) and `"negative"` (16 ion forms). Constants, column order and column
#' names are upstream's, so a library built here matches one built by
#' `CreateLibFromKEGG()` from the same masses.
#'
#' @param mw   Numeric vector of neutral monoisotopic masses.
#' @param mode `"positive"` or `"negative"`.
#' @return Numeric matrix, one row per mass, one named column per ion form.
#' @examples
#' mmc_ma_adduct_matrix(c(180.063388, 146.105528), "negative")[, "M-H[-]"]
mmc_ma_adduct_matrix <- function(mw, mode = c("positive", "negative")) {
  mode <- match.arg(mode)
  if (!is.numeric(mw) || anyNA(mw)) {
    .mmc_stop("adduct masses must be numeric with no missing values.")
  }
  PROTON <- .MMC_MA_PROTON
  if (mode == "positive") {
    m <- cbind(mw, mw + PROTON, mw/2 + PROTON, mw/3 + PROTON,
               mw + 1.0034 + PROTON, mw/2 + 0.5017 + PROTON,
               mw/3 + 0.3344 + PROTON, mw + 1.9958 + PROTON,
               mw + 1.9972 + PROTON, mw + 21.9820 + PROTON,
               mw/2 + 10.991 + PROTON, mw + 37.9555 + PROTON,
               mw + 18.0106 + PROTON, mw - 18.0106 + PROTON,
               mw - 36.0212 + PROTON, mw - 17.0265 + PROTON,
               mw - 27.9950 + PROTON, mw - 43.9898 + PROTON,
               mw - 46.0054 + PROTON, mw + 67.9874 + PROTON,
               mw - 67.9874 + PROTON, mw + 57.9586 + PROTON,
               mw - 72.0211 + PROTON, mw + 83.9613 + PROTON,
               mw - 83.9613 + PROTON)
    colnames(m) <- c("M[1+]", "M+H[1+]", "M+2H[2+]", "M+3H[3+]",
                     "M(C13)+H[1+]", "M(C13)+2H[2+]", "M(C13)+3H[3+]",
                     "M(S34)+H[1+]", "M(Cl37)+H[1+]", "M+Na[1+]",
                     "M+H+Na[2+]", "M+K[1+]", "M+H2O+H[1+]", "M-H2O+H[1+]",
                     "M-H4O2+H[1+]", "M-NH3+H[1+]", "M-CO+H[1+]",
                     "M-CO2+H[1+]", "M-HCOOH+H[1+]", "M+HCOONa[1+]",
                     "M-HCOONa+H[1+]", "M+NaCl[1+]", "M-C3H4O2+H[1+]",
                     "M+HCOOK[1+]", "M-HCOOK+H[1+]")
  } else {
    m <- cbind(mw - PROTON, mw/2 - PROTON, mw + 1.0034 - PROTON,
               mw + 1.9958 - PROTON, mw + 1.9972 - PROTON,
               mw + 21.9820 - 2*PROTON, mw + 37.9555 - 2*PROTON,
               mw - 18.0106 - PROTON, mw + 34.9689, mw + 36.9659,
               mw + 78.9183, mw + 80.9163,
               mw + 2*12 + 3*1.007825 + 14.00307 - PROTON,
               mw + 1.007825 + 12 + 2*15.99491,
               mw + 3*1.007825 + 2*12 + 2*15.99491,
               mw - PROTON + 15.99491)
    colnames(m) <- c("M-H[-]", "M-2H[2-]", "M(C13)-H[-]", "M(S34)-H[-]",
                     "M(Cl37)-H[-]", "M+Na-2H[-]", "M+K-2H[-]", "M-H2O-H[-]",
                     "M+Cl[-]", "M+Cl37[-]", "M+Br[-]", "M+Br81[-]",
                     "M+ACN-H[-]", "M+HCOO[-]", "M+CH3COO[-]", "M-H+O[-]")
  }
  m
}

#' Integer-m/z lookup tree over an adduct matrix (MetaboAnalystR shape)
#'
#' Port of the tree loop in MetaboAnalystR's `CreateLibFromKEGG()`. Element `n`
#' of the result holds the row indices of every compound with at least one
#' adduct m/z whose floor is `n`, for `50 < n < 2000`. Upstream seeds the list
#' with `""` at positions 49 and 2001 (fixing its length at 2001); that quirk is
#' kept so the object is identical to upstream's, and the matcher never reads
#' those two positions.
#'
#' @param mz_mat Numeric adduct matrix from `mmc_ma_adduct_matrix()`.
#' @return A list of length 2001; elements are integer row indices, `NULL`, or
#'   `""` at positions 49 and 2001.
#' @examples
#' tr <- mmc_ma_cpd_tree(mmc_ma_adduct_matrix(180.063388, "positive"))
#' tr[[181]]
mmc_ma_cpd_tree <- function(mz_mat) {
  l2 <- list()
  l2[[49]] <- ""
  l2[[2001]] <- ""
  floor_mzs <- floor(mz_mat)
  for (i in seq_len(nrow(floor_mzs))) {
    for (n in floor_mzs[i, ]) {
      if ((n > 50) & (n < 2000)) {
        l2[[n]] <- append(l2[[n]], i)
      }
    }
  }
  lapply(l2, unique)
}

#' Read the compound table of a model JSON (id, name, neutral mass)
#'
#' Only the Azimuth / metDataModel shape carries `neutral_mono_mass`, which the
#' library needs; mummichog's built-in model format does not, and is refused
#' rather than given invented masses.
#'
#' @param jm   Parsed model JSON (`simplifyVector = FALSE`).
#' @param path Source path, for messages.
#' @return data.frame with `id`, `name`, `mw`, in the model's compound order.
#' @noRd
.mmc_ma_model_compounds <- function(jm, path) {
  cpds <- jm$list_of_compounds
  if (is.null(cpds) || length(cpds) == 0L) {
    .mmc_stop("the model '", path, "' has no 'list_of_compounds'. The ",
              "MetaboAnalyst-compatible engine needs a model in the Azimuth / ",
              "metDataModel shape (as published by multiomics-annotation-prep), ",
              "whose compounds carry 'neutral_mono_mass'.")
  }
  str1 <- function(x) {
    if (is.null(x) || length(x) != 1L) NA_character_ else as.character(x)
  }
  id   <- vapply(cpds, function(c) str1(c$id), character(1))
  name <- vapply(cpds, function(c) str1(c$name), character(1))
  mw   <- vapply(cpds, function(c) {
    v <- c$neutral_mono_mass
    if (is.null(v) || length(v) != 1L || !is.numeric(v)) NA_real_ else as.numeric(v)
  }, numeric(1))

  bad_id <- is.na(id) | !nzchar(id)
  if (any(bad_id)) {
    .mmc_stop(sum(bad_id), " compound(s) in '", path, "' have no 'id'.")
  }
  dup <- unique(id[duplicated(id)])
  if (length(dup) > 0L) {
    .mmc_stop("compound ids must be unique in '", path, "'; duplicated: ",
              paste(utils::head(dup, 5), collapse = ", "), ".")
  }
  bad_mw <- is.na(mw) | !is.finite(mw) | mw <= 0
  if (any(bad_mw)) {
    .mmc_stop(sum(bad_mw), " compound(s) in '", path, "' have no usable ",
              "'neutral_mono_mass' (e.g. ",
              paste(utils::head(id[bad_mw], 5), collapse = ", "), "). Every ",
              "library compound needs its neutral monoisotopic mass; none is ",
              "guessed.")
  }
  name[is.na(name)] <- id[is.na(name)]
  data.frame(id = id, name = name, mw = mw, stringsAsFactors = FALSE)
}

#' Build a MetaboAnalystR mummichog library from a model JSON
#'
#' Produces the object MetaboAnalystR's `.setup.psea.library()` reads from a
#' `<lib>.qs` file, built from one of our metabolic models instead:
#'
#' \describe{
#'   \item{cpd.lib}{`id`, `name`, `mw` (neutral monoisotopic mass), and
#'     `adducts$positive` / `adducts$negative` (adduct m/z matrices, one row per
#'     compound).}
#'   \item{cpd.tree}{`positive` / `negative` integer-m/z lookup trees.}
#'   \item{pathways}{`cpds` (unnamed list of compound-id vectors), `id`, `name`.}
#'   \item{model}{Provenance: the model's `id`, `version` and source path.}
#' }
#'
#' Pathway membership is read with `read_mummichog_model_pathways()`, the same
#' reaction-derived membership the python engine and the evidence layer use.
#' Pathway compounds missing from the compound table are dropped (they can
#' never be matched); upstream's own library builder keeps only compounds of
#' its table the same way. Pathway names must be unique (MetaboAnalyst keys
#' pathways by name); a model with a repeated name stops.
#'
#' @param path Path to a model JSON in the Azimuth / metDataModel shape.
#' @return The library list described above.
#' @examples
#' \dontrun{
#' lib <- mmc_ma_library_from_model("envs/mummichog-models/<sha256>.json")
#' str(lib$cpd.lib$adducts$negative)
#' }
mmc_ma_library_from_model <- function(path) {
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
  cpd <- .mmc_ma_model_compounds(jm, path)
  pw  <- read_mummichog_model_pathways(path)

  pw_cpds <- lapply(unname(pw$pathways), function(x) {
    x <- as.character(x)
    x[x %in% cpd$id]
  })
  pw_names <- names(pw$pathways)
  # MetaboAnalyst keys every pathway by its name (result rows, set lookups), so
  # two pathways sharing a name would be indistinguishable downstream, and a
  # name lookup would hand both the first one's set and statistics. Refuse
  # rather than rename: a renamed pathway would no longer match the model.
  dup <- unique(pw_names[duplicated(pw_names)])
  if (length(dup) > 0L) {
    .mmc_stop(length(dup), " pathway name(s) occur more than once in '", path,
              "' (e.g. '", paste(utils::head(dup, 3), collapse = "', '"),
              "'). The MetaboAnalyst-compatible engine identifies pathways by ",
              "name, as MetaboAnalyst does, so each name must be unique. Give ",
              "the duplicated pathways distinct names in the model, or use ",
              "engine: python.")
  }
  pw_ids   <- unname(pw$pathway_ids[pw_names])

  adducts <- list(positive = mmc_ma_adduct_matrix(cpd$mw, "positive"),
                  negative = mmc_ma_adduct_matrix(cpd$mw, "negative"))
  list(
    pathways = list(cpds = pw_cpds, id = pw_ids, name = pw_names),
    cpd.tree = list(positive = mmc_ma_cpd_tree(adducts$positive),
                    negative = mmc_ma_cpd_tree(adducts$negative)),
    cpd.lib  = list(id = cpd$id, name = cpd$name, mw = cpd$mw,
                    adducts = adducts),
    model    = list(id = as.character(jm$id %||% NA_character_),
                    version = as.character(jm$version %||% NA_character_),
                    path = path)
  )
}
