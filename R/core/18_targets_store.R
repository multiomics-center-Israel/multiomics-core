#' Exact identity of one project run: its raw name and analysis round
#'
#' The store name, the ownership stamp and the error messages all come from
#' this one reading of the config, so they cannot disagree about which run a
#' store belongs to.
#'
#' @param config Config list with \code{project$name} (and optional
#'   \code{project$analysis_round}).
#' @return List with \code{name} (one string) and \code{analysis_round} (one
#'   string, or NULL when the config gives none).
targets_store_identity <- function(config) {
    one_value <- function(v, what) {
        if (is.null(v)) return(NULL)
        if (!(is.character(v) || is.numeric(v)) || length(v) != 1 || is.na(v)) {
            stop("config$project$", what, " must be a single value to name the ",
                 "targets store.", call. = FALSE)
        }
        enc2utf8(as.character(v))
    }
    name <- one_value(config$project$name, "name")
    if (is.null(name) || !nzchar(trimws(name))) {
        stop("config$project$name is missing; cannot name the targets store.",
             call. = FALSE)
    }
    list(name = name, analysis_round = one_value(config$project$analysis_round,
                                                  "analysis_round"))
}


#' Unambiguous encoding of a run identity
#'
#' Each field is written with its byte length ("<n>:<value>", or "-" when
#' absent), so no two different name/round pairs share an encoding -- unlike
#' \code{paste(name, round)}, under which "A B" + "C" and "A" + "B C" agree.
#'
#' @param identity Output of \code{targets_store_identity()}.
#' @return One string.
targets_store_identity_key <- function(identity) {
    field <- function(v) if (is.null(v)) "-" else paste0(nchar(v, type = "bytes"), ":", v)
    paste0("v1|", field(identity$name), "|", field(identity$analysis_round))
}


#' Name the {targets} store that belongs to one project run
#'
#' Every run gets its own store, derived from the config, so a run can never
#' read or destroy another project's cache. The name is a readable slug of
#' \code{project$name} and \code{project$analysis_round}, followed by a short
#' hash of their exact values: the slug alone would merge "My Project-X" with
#' "My Project_X", or "ProjA" with "proja" on a case-insensitive filesystem.
#' \code{project$targets_store} overrides the name, so an existing hand-named
#' store can keep being used; see \code{validate_targets_store_name()}.
#'
#' @param config Config list with \code{project$name} (and optional
#'   \code{project$analysis_round}, \code{project$targets_store}).
#' @return Store directory name, relative to the pipeline root.
targets_store_for <- function(config) {
    override <- config$project$targets_store
    if (!is.null(override)) return(validate_targets_store_name(override))
    identity <- targets_store_identity(config)
    slug <- tolower(gsub("[^A-Za-z0-9]+", "_",
                         paste(c(identity$name, identity$analysis_round), collapse = "_")))
    slug <- substr(gsub("^_+|_+$", "", slug), 1L, 40L)
    hash <- substr(digest::digest(targets_store_identity_key(identity), algo = "sha256",
                                  serialize = FALSE), 1L, 12L)
    paste0("_targets_", if (nzchar(slug)) paste0(slug, "_") else "", hash)
}


#' Check a configured store name before anything can build in or destroy it
#'
#' \code{run.R --fresh} calls \code{tar_destroy()} on this store, so the name
#' is held to what a store of this pipeline looks like: one directory directly
#' under the pipeline root, named \code{_targets...} from letters, digits,
#' \code{.}, \code{_} and \code{-}. That rules out absolute paths, \code{.} and
#' \code{..}, any traversal, and names such as \code{R} or \code{data} that
#' would turn a project folder into a destroyable store.
#'
#' @param store The configured \code{project$targets_store}.
#' @return \code{store}, when it is acceptable.
validate_targets_store_name <- function(store) {
    if (!is.character(store) || length(store) != 1 || is.na(store) ||
        !nzchar(trimws(store))) {
        stop("config$project$targets_store must be one non-empty store name.",
             call. = FALSE)
    }
    if (!grepl("^_targets[A-Za-z0-9._-]*$", store) || grepl("\\.\\.", store)) {
        stop("config$project$targets_store '", store, "' is not a safe store name: ",
             "it must be a single folder under the pipeline root, named '_targets...' ",
             "with only letters, digits, '.', '_' and '-' (no '/', '\\', '..' or ",
             "absolute path).", call. = FALSE)
    }
    store
}


#' Path of the ownership stamp inside a store
#'
#' Lives under \code{user/}, the folder {targets} reserves for user files, so
#' the pipeline never touches it and \code{tar_destroy()} removes it with the
#' store.
#'
#' @param store_dir Store directory.
#' @return Path of the stamp file.
targets_store_stamp_path <- function(store_dir) {
    file.path(store_dir, "user", "multiomics_project.yaml")
}


#' Read a store's ownership stamp, failing closed on anything malformed
#'
#' @param stamp Path of the stamp file.
#' @return The stamp as a list with \code{identity}, \code{project_name} and
#'   optionally \code{analysis_round}.
read_targets_store_stamp <- function(stamp) {
    bad <- function(why) {
        stop("Targets store ownership stamp '", stamp, "' is unreadable (", why, "). ",
             "Refusing to use the store. Remove the stamp only if you are sure the ",
             "store belongs to this project.", call. = FALSE)
    }
    s <- tryCatch(yaml::read_yaml(stamp), error = function(e) bad(conditionMessage(e)))
    if (!is.list(s)) bad("not a map of fields")
    if (!identical(s$stamp_version, 1L)) bad("unknown stamp_version")
    if (!is.character(s$identity) || length(s$identity) != 1 || is.na(s$identity) ||
        !startsWith(s$identity, "v1|")) {
        bad("no identity")
    }
    s
}


#' Stop unless a store may be used by this run
#'
#' A store stamped for another run is refused, and so is a malformed stamp. An
#' existing directory with no stamp is accepted only if it already looks like a
#' {targets} store (it has \code{meta/meta}) or is empty; anything else is not
#' a store, and pointing {targets} at it would let \code{--fresh} delete it. A
#' missing store is accepted.
#'
#' @param store_dir Store directory (path, not just its name).
#' @param config Config list of the run about to use it.
#' @return TRUE (invisibly) when the run may use the store.
assert_targets_store_owner <- function(store_dir, config) {
    identity <- targets_store_identity(config)
    stamp <- targets_store_stamp_path(store_dir)
    if (file.exists(stamp)) {
        s <- read_targets_store_stamp(stamp)
        if (!identical(s$identity, targets_store_identity_key(identity))) {
            stop(sprintf(
                "Targets store '%s' belongs to project '%s' (round '%s'), not '%s' (round '%s'). Refusing to use it. ",
                store_dir, s$project_name %||% "?", s$analysis_round %||% "none",
                identity$name, identity$analysis_round %||% "none"),
                "Set project$targets_store to a store of this project.",
                call. = FALSE)
        }
        return(invisible(TRUE))
    }
    if (dir.exists(store_dir) && length(list.files(store_dir, all.files = TRUE,
                                                   no.. = TRUE)) > 0 &&
        !file.exists(file.path(store_dir, "meta", "meta"))) {
        stop("'", store_dir, "' exists but is not a targets store (no meta/meta and no ",
             "ownership stamp). Refusing to use it: a --fresh run would delete it.",
             call. = FALSE)
    }
    invisible(TRUE)
}


#' Stamp a store with the run that uses it
#'
#' @param store_dir Store directory.
#' @param config Config list.
#' @return The stamp path (invisibly).
stamp_targets_store <- function(store_dir, config) {
    identity <- targets_store_identity(config)
    stamp <- targets_store_stamp_path(store_dir)
    dir.create(dirname(stamp), recursive = TRUE, showWarnings = FALSE)
    fields <- list(stamp_version = 1L,
                   identity = targets_store_identity_key(identity),
                   project_name = identity$name)
    if (!is.null(identity$analysis_round)) fields$analysis_round <- identity$analysis_round
    yaml::write_yaml(fields, stamp)
    invisible(stamp)
}


#' Point {targets} at this run's own store, after checking it may be used
#'
#' Resolves and validates the store, checks its ownership, and only then writes
#' \code{<store>.yaml} next to \code{_targets.R} and sets \code{TAR_CONFIG} to
#' it -- so a refused store never becomes the active one. \code{TAR_CONFIG} is
#' an environment variable, so the callr worker behind \code{tar_make()} and any
#' \code{tar_read()} in a report inherit it. The caller restores it afterwards
#' with \code{restore_tar_config()}.
#'
#' @param config Config list.
#' @param root Pipeline root (where \code{_targets.R} lives).
#' @return Store directory name (invisibly).
use_targets_store <- function(config, root = getwd()) {
    store <- targets_store_for(config)
    assert_targets_store_owner(file.path(root, store), config)
    yaml_path <- file.path(root, paste0(store, ".yaml"))
    writeLines(c("main:",
                 sprintf("  store: %s", store),
                 "  script: _targets.R"),
               yaml_path)
    Sys.setenv(TAR_CONFIG = yaml_path)
    invisible(store)
}


#' Put TAR_CONFIG back as it was
#'
#' \code{run.R} is also used from an interactive session, which must not stay
#' pointed at the last project's store once a run ends, whether it succeeded or
#' failed.
#'
#' @param old The value \code{Sys.getenv("TAR_CONFIG", unset = NA)} returned
#'   before the run: a string to restore, or NA to unset.
#' @return NULL (invisibly).
restore_tar_config <- function(old) {
    if (is.na(old)) Sys.unsetenv("TAR_CONFIG") else Sys.setenv(TAR_CONFIG = old)
    invisible(NULL)
}
