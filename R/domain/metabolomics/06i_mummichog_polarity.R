# R/domain/metabolomics/06i_mummichog_polarity.R
#
# Polarity helpers for mixed-polarity mummichog runs
# (modes.metabolomics.enrichment.mummichog.ionization_mode: mixed).
#
# mummichog 2.7.0 runs one ionization mode per run (`-m`) and matches every
# submitted m/z against that mode's ion forms; it never reads a feature's own
# ion label. A dataset holding positive- and negative-ion features is therefore
# run once per polarity, each run with only that polarity's features as its
# background. These helpers decide each feature's polarity from a configured
# annotation column and name the per-polarity run units. They compute no
# statistic.


#' Read each feature's ion polarity from a configured annotation column
#'
#' Accepts ion labels that end in their charge sign, as Compound Discoverer
#' writes its `Reference Ion` (`[M+H]+1`, `[M-H]-1`, `[M-2H]2-`, `[M]+`), and
#' the plain words positive / negative / pos / neg in any case. Anything else
#' -- missing, blank, or a label without a sign -- stops the run: a feature
#' submitted under the wrong mode is matched against the wrong ion forms, so no
#' polarity is ever assumed.
#'
#' @param row_data Feature annotation table (`pre$row_data`).
#' @param column   Name of the column holding each feature's ion label or polarity.
#' @return Character vector, `"positive"` or `"negative"`, one per row of `row_data`.
#' @examples
#' rd <- data.frame(id = c("F1", "F2", "F3"),
#'                  ion = c("[M+H]+1", "[M-2H]2-", "negative"))
#' mmc_feature_polarity(rd, "ion")
mmc_feature_polarity <- function(row_data, column) {
  if (!is.character(column) || length(column) != 1L || is.na(column) ||
      !nzchar(column)) {
    .mmc_stop("ionization_mode 'mixed' needs polarity_column: the name of the ",
              "annotation column that gives each feature's ion polarity ",
              "(e.g. \"Reference Ion\").")
  }
  if (!column %in% names(row_data)) {
    .mmc_stop("polarity_column '", column, "' is not a column of the feature ",
              "annotations. Available columns: ",
              paste(names(row_data), collapse = ", "))
  }

  x <- trimws(as.character(row_data[[column]]))
  out <- rep(NA_character_, length(x))

  # "[...]" then an optional charge, the sign, an optional charge.
  ion_rx <- "^\\[.+\\]\\s*[0-9]*([+-])[0-9]*$"
  is_ion <- !is.na(x) & grepl(ion_rx, x)
  out[is_ion] <- ifelse(sub(ion_rx, "\\1", x[is_ion]) == "+",
                        "positive", "negative")
  word <- tolower(x)
  out[!is_ion & word %in% c("positive", "pos")] <- "positive"
  out[!is_ion & word %in% c("negative", "neg")] <- "negative"

  bad <- is.na(out)
  if (any(bad)) {
    seen <- x[bad]
    seen[is.na(seen) | !nzchar(seen)] <- "<empty>"
    .mmc_stop(sum(bad), " of ", length(x), " feature(s) have no readable ",
              "polarity in '", column, "' (e.g. ",
              paste(utils::head(unique(seen), 5), collapse = " | "), "). ",
              "Expected an ion label ending in its charge sign, such as [M+H]+1 ",
              "or [M-H]-1, or the word positive/negative. No polarity is ",
              "assumed: fix the column, or point polarity_column at one that ",
              "gives every feature's polarity.")
  }
  out
}

#' Label one per-polarity mummichog run
#'
#' The label is the stage's contrast key in mixed mode: it names the run's
#' output directory, its report tab and its export suffix.
#' `mmc_split_polarity_label()` reverses it.
#'
#' @param contrast DE contrast name.
#' @param polarity `"positive"` or `"negative"`.
#' @return Character label, e.g. `"A_vs_B, positive ions"`.
#' @examples
#' mmc_polarity_label("A_vs_B", "negative")
mmc_polarity_label <- function(contrast, polarity) {
  paste0(contrast, ", ", polarity, " ions")
}

#' Split a mummichog run label into its DE contrast and polarity
#'
#' A label written by `mmc_polarity_label()` gives back its contrast and
#' polarity; any other label (single-polarity runs) is the contrast itself, with
#' polarity `NA`.
#'
#' @param label Character vector of run labels.
#' @return `list(contrast, polarity)`, each a character vector as long as `label`.
#' @examples
#' mmc_split_polarity_label(c("A_vs_B, positive ions", "A_vs_B"))
mmc_split_polarity_label <- function(label) {
  rx  <- "^(.*), (positive|negative) ions$"
  hit <- !is.na(label) & grepl(rx, label)
  list(contrast = ifelse(hit, sub(rx, "\\1", label), label),
       polarity = ifelse(hit, sub(rx, "\\2", label), NA_character_))
}

#' Resolve a mummichog run label to its DE contrast and polarity
#'
#' The one place the report decides whether a run label is a polarity run. Only
#' the python engine's mixed-mode stage writes polarity labels, so the label is
#' split only when `ionization_mode` is `mixed` -- read as
#' `mod_mummichog_pinned()` reads it -- and the engine is python
#' (`mmc_engine()`). The MetaboAnalyst-compatible engine's mixed mode is one
#' joint run per contrast, labelled by the contrast. In any other case the
#' label is the contrast itself, even when a contrast's own name happens to end
#' in ", positive ions".
#'
#' @param label  Character vector of run labels.
#' @param config Full pipeline config.
#' @return `list(contrast, polarity)`, as `mmc_split_polarity_label()`.
#' @examples
#' cfg <- list(modes = list(metabolomics = list(enrichment = list(
#'   mummichog = list(ionization_mode = "mixed")))))
#' mmc_run_label_parts("A_vs_B, positive ions", cfg)
mmc_run_label_parts <- function(label, config) {
  mode <- config$modes$metabolomics$enrichment$mummichog$ionization_mode %||% ""
  if (identical(tolower(as.character(mode)), "mixed") &&
      identical(mmc_engine(config), "python")) {
    return(mmc_split_polarity_label(label))
  }
  list(contrast = label, polarity = rep(NA_character_, length(label)))
}
