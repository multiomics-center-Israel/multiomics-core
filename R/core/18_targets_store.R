#' Name the {targets} store that belongs to one project run
#'
#' Every run gets its own store, derived from the config, so a run can never
#' read or destroy another project's cache. The name follows the run's output
#' folder (\code{project$name} + \code{project$analysis_round}), because that is
#' the unit a cache describes. \code{project$targets_store} overrides it, which
#' lets an existing hand-named store keep being used.
#'
#' @param config Config list with \code{project$name} (and optional
#'   \code{project$analysis_round}, \code{project$targets_store}).
#' @return Store directory name, relative to the pipeline root.
targets_store_for <- function(config) {
    override <- config$project$targets_store
    if (!is.null(override) && nzchar(override)) return(override)
    name <- config$project$name
    if (is.null(name) || !nzchar(name)) {
        stop("config$project$name is missing; cannot name the targets store.")
    }
    parts <- c(name, config$project$analysis_round)
    slug <- tolower(gsub("[^A-Za-z0-9]+", "_", paste(parts, collapse = "_")))
    paste0("_targets_", gsub("^_+|_+$", "", slug))
}

#' Owner label written into a store, so the store can say whose it is
#'
#' @param config Config list.
#' @return Single string, e.g. \code{"MyProject A01"}.
targets_store_owner <- function(config) {
    trimws(paste(config$project$name, config$project$analysis_round %||% ""))
}

#' Path of the ownership stamp inside a store
#'
#' Lives under \code{user/}, the folder {targets} reserves for user files, so
#' the pipeline never touches it and \code{tar_destroy()} removes it with the
#' store.
targets_store_stamp_path <- function(store) {
    file.path(store, "user", "multiomics_project")
}

#' Stop if a store is stamped for a different project
#'
#' A store with no stamp (created before stamping existed) is accepted, as is
#' a missing store.
#'
#' @param store Store directory.
#' @param config Config list of the run about to use it.
#' @return TRUE (invisibly) when the run may use the store.
assert_targets_store_owner <- function(store, config) {
    stamp <- targets_store_stamp_path(store)
    if (!file.exists(stamp)) return(invisible(TRUE))
    owner <- readLines(stamp, n = 1L, warn = FALSE)
    expected <- targets_store_owner(config)
    if (!identical(owner, expected)) {
        stop(sprintf(
            "Targets store '%s' belongs to '%s', not '%s'. Refusing to use it. ",
            store, owner, expected),
            "Set project$targets_store to a store of this project.",
            call. = FALSE)
    }
    invisible(TRUE)
}

#' Stamp a store with the project that uses it
#'
#' @param store Store directory.
#' @param config Config list.
#' @return The stamp path (invisibly).
stamp_targets_store <- function(store, config) {
    stamp <- targets_store_stamp_path(store)
    dir.create(dirname(stamp), recursive = TRUE, showWarnings = FALSE)
    writeLines(targets_store_owner(config), stamp)
    invisible(stamp)
}

#' Point {targets} at this run's own store
#'
#' Writes \code{_targets_<run>.yaml} next to \code{_targets.R} and sets
#' \code{TAR_CONFIG} to it. \code{TAR_CONFIG} is an environment variable, so
#' the callr worker behind \code{tar_make()} and any \code{tar_read()} in a
#' report inherit it. The shared \code{_targets.yaml} is no longer consulted.
#'
#' @param config Config list.
#' @param root Pipeline root (where \code{_targets.R} lives).
#' @return Store directory name (invisibly).
use_targets_store <- function(config, root = getwd()) {
    store <- targets_store_for(config)
    yaml_path <- file.path(root, paste0(store, ".yaml"))
    writeLines(c("main:",
                 sprintf("  store: %s", store),
                 "  script: _targets.R"),
               yaml_path)
    Sys.setenv(TAR_CONFIG = yaml_path)
    invisible(store)
}
