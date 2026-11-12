# R/domain/metabolomics/06q_enrichment_reference_files.R
#
# The reference files the metabolomics enrichment reads, and their fingerprint
# for {targets}.
#
# Two files outside the data folder decide enrichment results: the GMT pathway
# libraries (`enrichment.gmt_file`) and the HMDB -> KEGG table
# (`enrichment.mapping_file`). Neither was a dependency of any target, so
# editing one left the cached results in place.
#
# They are tracked by FINGERPRINT, a named vector `path -> md5` taken while the
# plan is being defined (each time _targets.R is sourced) and written into the
# command of a target. {targets} reruns a target whose command changed, so:
#   - a file edited in place changes its md5      -> the target reruns;
#   - a configured file that was missing appears  -> a new entry, it reruns;
#   - nothing changed                             -> same command, nothing runs.
# The third line is why this is not a `format = "file"` target fed by an
# always-running discovery target: a target that always runs is always reported
# as outdated, together with everything downstream of it.
#
# A configured file that does not exist is left out, as the enrichment code
# skips it (with a warning from QEA). It is not an error here either.

#' GMT library files the enrichment stage reads
#'
#' Resolved exactly as the enrichment functions resolve
#' \code{enrichment.gmt_file}: through \code{resolve_input_path()}, so an
#' absolute path is kept and a relative one is located under the raw data
#' folder.
#'
#' @param config Config list (raw YAML or validated; needs \code{project$dir}).
#' @return Character vector of the configured GMT files that exist, without
#'   duplicates, in configured order. \code{character(0)} when none.
find_metab_enrichment_gmt_files <- function(config) {
    enr_cfg <- config$modes$metabolomics$enrichment %||% list()
    .existing_reference_files(resolve_input_path(config, enr_cfg$gmt_file))
}

#' HMDB -> KEGG mapping file, at every path a reader uses
#'
#' Two stages read \code{enrichment.mapping_file}. The enrichment functions
#' resolve it through \code{resolve_input_path()}; \code{mod_met_raw()} passes
#' the configured value on as written. For an absolute path these are one file.
#' For a relative path they can be two, and both are returned so that each
#' reader's file is tracked. Which of them a reader SHOULD use is a separate
#' question that this function does not settle.
#'
#' @param config Config list (raw YAML or validated; needs \code{project$dir}).
#' @return Character vector of the existing mapping file path(s), without
#'   duplicates. \code{character(0)} when none is configured or none exists.
find_metab_kegg_mapping_files <- function(config) {
    enr_cfg  <- config$modes$metabolomics$enrichment %||% list()
    as_given <- unlist(enr_cfg$mapping_file, use.names = FALSE)
    .existing_reference_files(c(resolve_input_path(config, as_given), as_given))
}

#' Keep the paths that are existing files, once each
#'
#' @param paths Character vector of candidate paths (may be NULL).
#' @return Character vector of normalised paths of existing regular files.
#' @noRd
.existing_reference_files <- function(paths) {
    paths <- as.character(unlist(paths, use.names = FALSE))
    paths <- paths[!is.na(paths) & nzchar(paths)]
    paths <- paths[file.exists(paths) & !dir.exists(paths)]
    if (length(paths) == 0) return(character(0))
    unique(normalizePath(paths, winslash = "/", mustWork = TRUE))
}

#' Fingerprint a set of files
#'
#' @param paths Character vector of existing file paths.
#' @return Named character vector, md5 by path, in the order given.
#'   \code{character(0)} for no paths.
fingerprint_reference_files <- function(paths) {
    paths <- as.character(paths)
    if (length(paths) == 0) return(character(0))
    stats::setNames(unname(tools::md5sum(paths)), paths)
}

#' Fingerprints of the metabolomics enrichment reference files
#'
#' Called while the plan is defined, with the raw config, and passed to
#' \code{pipe_metabolomics()}.
#'
#' @param config Config list (raw YAML or validated; needs \code{project$dir}).
#' @return list(gmt, mapping): two named character vectors, md5 by path.
metab_reference_fingerprints <- function(config) {
    list(
        gmt     = fingerprint_reference_files(find_metab_enrichment_gmt_files(config)),
        mapping = fingerprint_reference_files(find_metab_kegg_mapping_files(config))
    )
}
