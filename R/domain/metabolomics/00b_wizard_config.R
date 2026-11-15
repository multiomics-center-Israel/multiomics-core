# R/domain/metabolomics/00b_wizard_config.R
#
# The `preprocessing:` block the config wizards write for metabolomics.
#
# Self-contained on purpose (base R only, no helpers from other files): the
# terminal wizard in run.R sources this one file before any of R/ is loaded.
# wizard.html builds the same block in JavaScript (buildMetabPreprocessingYaml);
# tests/testthat/test-metab-wizard-preprocessing.R keeps the two in step.

#' Sample-normalization methods the metabolomics wizards offer
#'
#' The values are `preprocessing.chosen_norm` values, written to the config as
#' they are. Internal-standard normalization is not among them: the
#' metabolomics pipeline has no such branch.
metab_wizard_norm_methods <- c("pqn", "median", "tss", "eigenms",
                               "eigenms_forced", "bio_factor", "none")

#' Does a normalization method honour `preprocessing.transform`?
#'
#' Only `median` and `none` are built on the configurable transform of
#' `met_log`. Every other method always takes `log2(x + pseudocount)`, so
#' `transform` has no effect on it.
#'
#' @param chosen_norm Character scalar, a `preprocessing.chosen_norm` value.
#' @return `TRUE` when `transform` is applied for that method, else `FALSE`.
metab_norm_honours_transform <- function(chosen_norm) {
  tolower(chosen_norm) %in% c("median", "none")
}

#' Build the metabolomics `preprocessing:` YAML block for a wizard config
#'
#' Writes the keys the pipeline reads, under `preprocessing:`; there is no
#' `normalization:` block, which config validation rejects. A transform other
#' than log2 is refused for the methods that always use log2, so a wizard
#' config never describes a transform that would not run.
#'
#' @param chosen_norm One of `metab_wizard_norm_methods`.
#' @param transform `"log2"`, `"log10"`, `"glog10"` or `"none"`; must be
#'   `"log2"` unless `metab_norm_honours_transform(chosen_norm)`.
#' @param scaling `"none"`, `"center"`, `"auto"`, `"pareto"` or `"range"`.
#' @param qc_only `TRUE` writes `chosen_norm: null` (QC-review mode) and keeps
#'   the picked method as a comment for the later full run.
#' @param biological_factor_col Metadata column for `bio_factor`; required for
#'   that method, ignored otherwise.
#' @param relaxed_missingness `TRUE` adds the looser missingness thresholds
#'   the wizard uses for multi-level input.
#' @return Character scalar: the block, indented for `modes: metabolomics:`,
#'   without a trailing newline.
build_metab_wizard_preprocessing_yaml <- function(chosen_norm,
                                                  transform = "log2",
                                                  scaling = "none",
                                                  qc_only = FALSE,
                                                  biological_factor_col = NULL,
                                                  relaxed_missingness = FALSE) {
  one_of <- function(x, what, choices) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% choices) {
      stop(sprintf("build_metab_wizard_preprocessing_yaml: %s must be one of: %s.",
                   what, paste(choices, collapse = ", ")), call. = FALSE)
    }
  }
  one_of(chosen_norm, "chosen_norm", metab_wizard_norm_methods)
  one_of(transform, "transform", c("log2", "log10", "glog10", "none"))
  one_of(scaling, "scaling", c("none", "center", "auto", "pareto", "range"))

  if (transform != "log2" && !metab_norm_honours_transform(chosen_norm)) {
    stop(sprintf(paste0(
      "build_metab_wizard_preprocessing_yaml: transform \"%s\" cannot be ",
      "combined with chosen_norm \"%s\", which always uses log2. Use ",
      "transform \"log2\", or pick \"median\" or \"none\"."),
      transform, chosen_norm), call. = FALSE)
  }

  lines <- "    preprocessing:"
  if (isTRUE(qc_only)) {
    lines <- c(lines,
      sprintf("      # QC-review mode. For the full analysis set chosen_norm: \"%s\" and rerun.",
              chosen_norm),
      "      chosen_norm: null")
  } else {
    lines <- c(lines, sprintf("      chosen_norm: \"%s\"", chosen_norm))
  }
  lines <- c(lines,
    sprintf("      transform: \"%s\"", transform),
    sprintf("      scaling: \"%s\"", scaling),
    "      pseudocount: 1")

  if (chosen_norm == "bio_factor") {
    col <- biological_factor_col
    if (!is.character(col) || length(col) != 1L || is.na(col) || !nzchar(trimws(col))) {
      stop("build_metab_wizard_preprocessing_yaml: chosen_norm \"bio_factor\" needs ",
           "biological_factor_col, the metadata column holding the per-sample ",
           "amount to divide by (for example \"mg_protein\").", call. = FALSE)
    }
    # Single-quoted YAML scalar: the only character to escape is the quote.
    lines <- c(lines, sprintf("      biological_factor_col: '%s'",
                              gsub("'", "''", trimws(col), fixed = TRUE)))
  }
  if (isTRUE(relaxed_missingness)) {
    lines <- c(lines,
      "      feat_missing_threshold: 0.80",
      "      sample_missing_threshold: 0.99")
  }
  paste(lines, collapse = "\n")
}
