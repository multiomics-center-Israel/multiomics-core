# Omics bundles: the package one finished single-omics run leaves behind for
# integration. A bundle is a directory
#
#   <root>/bundles/<omics>/<bundle_id>/
#       manifest.json        identity, provenance, and a sha256 per part
#       checksums.sha256     the sha256 of manifest.json
#       parts/<name>.rds     one file per part, exactly as the run produced it
#
# Bundles are immutable. A new run, or any changed part, gets a new ID; an
# existing ID is never written over. A bundle is staged under a hidden
# ".partial_*" folder, its destination is claimed atomically, and its files are
# moved in with checksums.sha256 last; readers ignore any folder without that
# file, so a run that dies half way leaves nothing a reader will pick up.
# Readers verify the manifest and every part they load before returning it.
#
# Trust assumption: bundles are pipeline-generated artefacts read by the same
# team. The checks guard data integrity -- accidental edits, damaged files and
# interrupted writes -- not deliberate tampering. In particular, symbolic links
# inside a bundle are not resolved or checked, so a link under parts/ is
# followed when a part is read.
#
# What goes into the parts is the caller's business (one assembler per omics);
# this file knows nothing about any omics.

.BUNDLE_CONTRACT_VERSION <- "1.0.0"

# Omics and part names become directory and file names, so they are held to
# lowercase letters, digits and "_": nothing that can escape the bundle root.
.BUNDLE_NAME_PATTERN <- "^[a-z][a-z0-9_]*$"

.BUNDLE_TIME_PATTERN <- "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$"


#' sha256 of a file's bytes, lowercase hex
#'
#' @param path File path.
#' @return One 64-character string.
.bundle_sha256_file <- function(path) {
    tolower(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}


#' Readable label for a bundle ID
#'
#' The same slug rule \code{targets_store_for()} uses for store names. The slug
#' only makes an ID readable: the exact project name and round are kept in the
#' manifest, and the content hash in the ID keeps two bundles apart.
#'
#' @param x One string.
#' @return Lowercase slug of at most 40 characters, never empty.
.bundle_slug <- function(x) {
    slug <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
    slug <- substr(gsub("^_+|_+$", "", slug), 1L, 40L)
    if (nzchar(slug)) slug else "x"
}


#' Capture where and with what a bundle is being written
#'
#' Called once by the code that writes a bundle, at write time, so the record
#' describes the run that produced the parts. It is deliberately not copied
#' from \code{execution_info/}: that target is cached, and goes stale when only
#' code changes. Git details are recorded as NA, not an error, when git or the
#' repository is unavailable.
#'
#' @param repo_dir Pipeline root: the git work tree and the folder holding
#'   \code{renv.lock}.
#' @return List with \code{created_at_utc} (ISO 8601, UTC), \code{git_commit},
#'   \code{git_dirty} (uncommitted changes to tracked files),
#'   \code{r_version}, \code{renv_lock_sha256} and \code{packages} (named list
#'   of the versions of every loaded namespace).
bundle_provenance <- function(repo_dir = ".") {
    git <- function(args) {
        out <- tryCatch(
            suppressWarnings(system2("git", c("-C", shQuote(repo_dir), args),
                                     stdout = TRUE, stderr = FALSE)),
            error = function(e) NULL)
        if (!is.null(attr(out, "status"))) NULL else out
    }
    commit <- git(c("rev-parse", "HEAD"))
    commit <- if (length(commit) == 1 && nzchar(commit)) commit else NA_character_
    status <- if (is.na(commit)) NULL else git(c("status", "--porcelain",
                                                 "--untracked-files=no"))
    # NULL means git status failed, so the state is unknown; character(0) means
    # it ran and found no changes.
    dirty <- if (is.null(status)) NA else length(status) > 0

    lock <- file.path(repo_dir, "renv.lock")
    pkgs <- sort(loadedNamespaces())

    list(
        created_at_utc   = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
        git_commit       = commit,
        git_dirty        = dirty,
        r_version        = paste(R.version$major, R.version$minor, sep = "."),
        renv_lock_sha256 = if (file.exists(lock)) .bundle_sha256_file(lock) else NA_character_,
        packages         = stats::setNames(
            lapply(pkgs, function(p) as.character(utils::packageVersion(p))), pkgs)
    )
}


#' Bundle identity of one single-omics run, read from its config
#'
#' Project and analysis round fall back exactly as \code{get_run_out_dir()}
#' does (\code{"Project"}, \code{"Analysis"}), so a bundle always names the run
#' its \code{Results_<name>_<round>} folder belongs to.
#'
#' @param config Validated config of the run.
#' @param omics Omics name for the bundle (e.g. \code{"rnaseq"}).
#' @return List with \code{omics}, \code{project}, \code{analysis_round} and
#'   \code{run_name}, ready to pass to \code{write_omics_bundle()} as
#'   \code{info}.
omics_bundle_info <- function(config, omics) {
    list(omics          = omics,
         project        = as.character(config$project$name %||% "Project"),
         analysis_round = as.character(config$project$analysis_round %||% "Analysis"),
         run_name       = get_run_name(config))
}


#' Write a new, immutable omics bundle
#'
#' Each part is saved as its own RDS file under a hidden \code{.partial_*}
#' directory, hashed, and listed in a manifest; the manifest's own sha256 goes
#' into \code{checksums.sha256}. The ID is
#' \code{<omics>-<project>-<round>-<UTC time>-<content hash>}, where the
#' content hash covers every part file. The destination folder is then claimed
#' with a non-recursive \code{dir.create()}, which fails if it exists at all,
#' even empty, so a bundle is never replaced and of two concurrent writers of
#' one ID only one succeeds. The staged files are moved in, with
#' \code{checksums.sha256} last, which is what marks the bundle complete.
#'
#' If writing fails, the \code{.partial_*} directory (and a destination left
#' without \code{checksums.sha256}) stay where they are and are named in the
#' error: readers never list or load them, and they are removed by hand once
#' looked at.
#'
#' @param parts Named list of R objects, at least one. Names must match
#'   \code{^[a-z][a-z0-9_]*$}; each becomes \code{parts/<name>.rds}.
#' @param info Named list with \code{omics} (same name rule), \code{project}
#'   and \code{analysis_round} (non-empty strings). Any further fields are kept
#'   in the manifest and must each be a single string, number or logical: the
#'   manifest describes a bundle, it does not hold data.
#' @param root Output root; the bundle goes under \code{<root>/bundles/<omics>/}.
#' @param provenance Output of \code{bundle_provenance()}, captured by the
#'   caller at write time.
#' @return Path of the new bundle directory.
write_omics_bundle <- function(parts, info, root, provenance) {
    if (!is.list(parts) || length(parts) == 0 || is.null(names(parts))) {
        stop("Bundle parts must be a non-empty named list.", call. = FALSE)
    }
    bad <- names(parts)[is.na(names(parts)) | !grepl(.BUNDLE_NAME_PATTERN, names(parts))]
    if (length(bad) > 0) {
        stop("Bundle part names must be lowercase letters, digits and '_', starting ",
             "with a letter; got: ", paste(sQuote(bad, FALSE), collapse = ", "), ".",
             call. = FALSE)
    }
    if (anyDuplicated(names(parts))) {
        stop("Bundle part names must be unique; repeated: ",
             paste(unique(names(parts)[duplicated(names(parts))]), collapse = ", "), ".",
             call. = FALSE)
    }

    assert_named_list(info, "info")
    for (f in c("omics", "project", "analysis_round")) assert_scalar_chr(info[[f]], paste0("info$", f))
    if (!grepl(.BUNDLE_NAME_PATTERN, info$omics)) {
        stop("info$omics must be lowercase letters, digits and '_' (e.g. 'rnaseq'); got '",
             info$omics, "'.", call. = FALSE)
    }
    extra <- info[setdiff(names(info), c("omics", "project", "analysis_round"))]
    not_scalar <- names(extra)[!vapply(extra, function(v) {
        (is.character(v) || is.numeric(v) || is.logical(v)) && length(v) == 1
    }, logical(1))]
    if (length(not_scalar) > 0) {
        stop("Extra info fields must each be one string, number or logical; not: ",
             paste(not_scalar, collapse = ", "), ". Put tables and objects in parts.",
             call. = FALSE)
    }

    assert_named_list(provenance, "provenance")
    created <- provenance$created_at_utc
    if (!is.character(created) || length(created) != 1 || !grepl(.BUNDLE_TIME_PATTERN, created)) {
        stop("provenance$created_at_utc must be a UTC time like '2026-01-31T09:30:00Z'; ",
             "pass the output of bundle_provenance().", call. = FALSE)
    }
    assert_scalar_chr(root, "root")

    omics_dir <- file.path(root, "bundles", info$omics)
    ensure_dir(omics_dir)
    if (!dir.exists(omics_dir)) {
        stop("Cannot create the bundle directory '", omics_dir, "'.", call. = FALSE)
    }
    stamp <- gsub("[-:]", "", created)
    # A staging name unique to this call (time, process, and tempfile()'s own
    # suffix), so neither a concurrent writer nor a folder an earlier write
    # left behind can block this one. Claimed with a non-recursive dir.create(),
    # which fails rather than share an existing folder.
    partial <- tempfile(paste0(".partial_", stamp, "_", Sys.getpid(), "_"), tmpdir = omics_dir)
    if (!dir.create(partial, showWarnings = FALSE)) {
        stop("Cannot create the staging folder '", partial, "'.", call. = FALSE)
    }

    # Everything up to the rename happens inside the hidden directory; a failure
    # leaves it there, named, rather than deleting files the user may need to
    # see.
    left_behind <- function(e) {
        stop("Writing the ", info$omics, " bundle failed: ", conditionMessage(e),
             "\nThe incomplete directory '", partial, "' was left in place; readers ",
             "ignore it. Remove it once you have looked at it.", call. = FALSE)
    }
    manifest <- tryCatch({
        dir.create(file.path(partial, "parts"))
        part_entries <- list()
        for (nm in names(parts)) {
            rel <- file.path("parts", paste0(nm, ".rds"))
            saveRDS(parts[[nm]], file.path(partial, rel))
            part_entries[[nm]] <- list(file = rel,
                                       sha256 = .bundle_sha256_file(file.path(partial, rel)),
                                       class = class(parts[[nm]])[1])
        }
        nm_sorted <- sort(names(part_entries))
        content_hash <- digest::digest(
            paste(nm_sorted, vapply(part_entries[nm_sorted], `[[`, "", "sha256"),
                  collapse = "\n"),
            algo = "sha256", serialize = FALSE)
        bundle_id <- paste(info$omics, .bundle_slug(info$project),
                           .bundle_slug(info$analysis_round), stamp,
                           substr(content_hash, 1L, 8L), sep = "-")
        m <- list(
            bundle_contract_version = .BUNDLE_CONTRACT_VERSION,
            bundle_id      = bundle_id,
            status         = "complete",
            omics          = info$omics,
            project        = info$project,
            analysis_round = info$analysis_round,
            created_at_utc = created,
            content_hash   = content_hash,
            info           = if (length(extra) > 0) extra else stats::setNames(list(), character(0)),
            provenance     = provenance,
            parts          = part_entries
        )
        jsonlite::write_json(m, file.path(partial, "manifest.json"), auto_unbox = TRUE,
                             pretty = TRUE, digits = NA, na = "null", null = "null")
        # sha256sum format, so `sha256sum -c checksums.sha256` works outside R too.
        writeLines(paste0(.bundle_sha256_file(file.path(partial, "manifest.json")),
                          "  manifest.json"),
                   file.path(partial, "checksums.sha256"))
        m
    }, error = left_behind)

    # Claim the destination. A non-recursive dir.create() is one atomic step
    # that fails when the folder exists, empty or not, so of two writers with
    # the same ID exactly one gets past here. Renaming the staging folder onto
    # the destination instead would silently replace an empty folder on Unix.
    final <- file.path(omics_dir, manifest$bundle_id)
    if (!dir.create(final, showWarnings = FALSE)) {
        why <- if (file.exists(final)) {
            "already exists; bundles are never overwritten"
        } else {
            "could not be created"
        }
        stop("Bundle '", manifest$bundle_id, "' in '", omics_dir, "' ", why,
             ". The new copy was left at '", partial, "'.", call. = FALSE)
    }

    # Publish into the claimed folder with checksums.sha256 last: readers treat
    # a bundle without it as incomplete, so a crash part way through leaves a
    # folder nobody lists or loads. Each rename stays on one filesystem.
    for (item in c("parts", "manifest.json", "checksums.sha256")) {
        if (!file.rename(file.path(partial, item), file.path(final, item))) {
            stop("Publishing bundle '", manifest$bundle_id, "' failed while moving ", item,
                 ". '", final, "' is incomplete and ignored by readers; remove it and '",
                 partial, "' before writing this bundle again.", call. = FALSE)
        }
    }
    # The staging folder should be empty now; it is removed only if it is.
    # unlink() deletes a folder only with recursive = TRUE, on Windows as on
    # Unix, and the emptiness check keeps that from ever deleting a file.
    if (length(list.files(partial, all.files = TRUE, no.. = TRUE)) == 0) {
        unlink(partial, recursive = TRUE)
    }

    # Best effort: a read-only bundle is harder to change by accident. Checksums
    # catch the changes this does not prevent.
    suppressWarnings(Sys.chmod(list.files(final, recursive = TRUE, full.names = TRUE),
                               mode = "0444"))
    final
}


#' Read a bundle's manifest, or NULL when the directory holds none
#'
#' @param path Bundle directory.
#' @return The manifest as a list, or NULL.
.read_bundle_manifest <- function(path) {
    file <- file.path(path, "manifest.json")
    if (!file.exists(file)) return(NULL)
    tryCatch(jsonlite::read_json(file, simplifyVector = TRUE),
             error = function(e) stop("The manifest of '", path, "' is not readable JSON: ",
                                      conditionMessage(e), call. = FALSE))
}


#' Check that a bundle is complete and its parts are the files it was written with
#'
#' Fails on: a hidden (\code{.partial_*}) or missing directory, no manifest, a
#' missing or malformed \code{checksums.sha256}, a manifest whose sha256 differs
#' from it (checked before the manifest is used), a status other than
#' complete, an unsupported contract version, a directory renamed away from its
#' bundle ID, any manifest part name not matching the part-name rule (checked
#' for all parts, before any path is built), an unknown part, a missing part
#' file, or a part whose sha256 differs from the manifest. Messages name the
#' bundle, the part and the file, never a value from inside a part. These checks
#' assume a trusted, pipeline-generated bundle: symbolic links are not checked
#' (see the note at the top of this file).
#'
#' @param path Bundle directory.
#' @param parts Names of the parts to verify; NULL verifies all of them.
#' @return The manifest (invisibly) when every check passes.
validate_omics_bundle <- function(path, parts = NULL) {
    assert_scalar_chr(path, "path")
    if (!dir.exists(path)) stop("Bundle directory '", path, "' does not exist.", call. = FALSE)
    if (startsWith(basename(path), ".")) {
        stop("'", path, "' is an incomplete bundle left by a failed write, not a bundle.",
             call. = FALSE)
    }
    if (!file.exists(file.path(path, "manifest.json"))) {
        stop("'", path, "' has no manifest.json: it is incomplete or not a bundle.",
             call. = FALSE)
    }
    # The manifest is checked before anything in it is trusted: it names the
    # run, its provenance and the part checksums verified below.
    sums <- file.path(path, "checksums.sha256")
    if (!file.exists(sums)) {
        stop("Bundle '", path, "' has no checksums.sha256: it is incomplete (an ",
             "interrupted write) or its manifest cannot be verified.", call. = FALSE)
    }
    line <- readLines(sums, warn = FALSE)
    # The whole line is checked: a bare hash would pass a hash-only test but
    # break `sha256sum -c`.
    if (length(line) != 1 || !grepl("^[0-9a-f]{64}  manifest\\.json$", line)) {
        stop("Bundle '", path, "': checksums.sha256 is not one '<sha256>  manifest.json' ",
             "line.", call. = FALSE)
    }
    if (!identical(.bundle_sha256_file(file.path(path, "manifest.json")),
                   substr(line, 1L, 64L))) {
        stop("Bundle '", path, "': manifest.json does not match checksums.sha256; the ",
             "manifest was changed after the bundle was written.", call. = FALSE)
    }
    m <- .read_bundle_manifest(path)
    if (!identical(m$status, "complete")) {
        stop("Bundle '", path, "' is not complete (status: ", m$status %||% "missing", ").",
             call. = FALSE)
    }
    major <- strsplit(as.character(m$bundle_contract_version %||% ""), ".", fixed = TRUE)[[1]][1]
    if (!identical(major, strsplit(.BUNDLE_CONTRACT_VERSION, ".", fixed = TRUE)[[1]][1])) {
        stop("Bundle '", path, "' uses contract version '",
             m$bundle_contract_version %||% "missing", "'; this pipeline reads ",
             .BUNDLE_CONTRACT_VERSION, ".", call. = FALSE)
    }
    if (!identical(m$bundle_id, basename(path))) {
        stop("Bundle directory '", basename(path), "' does not match its manifest ID '",
             m$bundle_id %||% "missing", "'. Bundles must keep the name they were written ",
             "with.", call. = FALSE)
    }

    available <- names(m$parts)
    # Every part name in the manifest is checked, not only the requested ones,
    # and before any path is built from one: a name such as "../x" would
    # otherwise produce an "expected" path outside the bundle.
    bad <- available[is.na(available) | !grepl(.BUNDLE_NAME_PATTERN, available)]
    if (length(available) == 0 || length(bad) > 0) {
        stop("Bundle '", m$bundle_id, "': the manifest has no parts or invalid part ",
             "name(s) (", paste(bad, collapse = ", "), "); part names must be lowercase ",
             "letters, digits and '_'.", call. = FALSE)
    }
    if (is.null(parts)) parts <- available
    unknown <- setdiff(parts, available)
    if (length(unknown) > 0) {
        stop("Bundle '", m$bundle_id, "' has no part(s) ", paste(unknown, collapse = ", "),
             ". Available: ", paste(available, collapse = ", "), ".", call. = FALSE)
    }
    for (nm in parts) {
        # Only the path the writer produces is accepted, so an edited manifest
        # cannot point a reader at a file outside the bundle.
        if (!identical(m$parts[[nm]]$file, file.path("parts", paste0(nm, ".rds")))) {
            stop("Bundle '", m$bundle_id, "': the manifest entry for part '", nm,
                 "' does not point to parts/", nm, ".rds.", call. = FALSE)
        }
        file <- file.path(path, m$parts[[nm]]$file)
        if (!file.exists(file)) {
            stop("Bundle '", m$bundle_id, "': part '", nm, "' is missing (", file, ").",
                 call. = FALSE)
        }
        if (!identical(.bundle_sha256_file(file), m$parts[[nm]]$sha256)) {
            stop("Bundle '", m$bundle_id, "': part '", nm, "' does not match its checksum; ",
                 "the file was changed after the bundle was written (", file, ").",
                 call. = FALSE)
        }
    }
    invisible(m)
}


#' Load parts of a bundle, after verifying them
#'
#' @param path Bundle directory.
#' @param parts Names of the parts to load; NULL loads every part.
#' @return List with \code{manifest} and \code{parts} (named list of the
#'   loaded objects).
read_omics_bundle <- function(path, parts = NULL) {
    m <- validate_omics_bundle(path, parts)
    wanted <- if (is.null(parts)) names(m$parts) else parts
    list(manifest = m,
         parts = stats::setNames(lapply(wanted, function(nm) {
             readRDS(file.path(path, m$parts[[nm]]$file))
         }), wanted))
}


#' List the complete bundles under an output root
#'
#' Reads manifests only, so it is cheap; checksums are verified when a bundle
#' is read. Hidden \code{.partial_*} directories, directories without
#' \code{checksums.sha256} (checked before any manifest is parsed),
#' directories without a complete or readable manifest, directories renamed
#' away from their bundle ID (which could never be read), and bundles sitting
#' under another omics' folder are left out. Rows are ordered by omics, then creation time,
#' then ID; the order is for display, and integration still names the bundles
#' it uses.
#'
#' @param root Output root, as passed to \code{write_omics_bundle()}.
#' @param omics Optional omics name to restrict the listing to.
#' @return Data frame with \code{bundle_id}, \code{omics}, \code{project},
#'   \code{analysis_round}, \code{created_at_utc} and \code{path}; zero rows
#'   when there are none.
list_omics_bundles <- function(root, omics = NULL) {
    assert_scalar_chr(root, "root")
    empty <- data.frame(bundle_id = character(0), omics = character(0),
                        project = character(0), analysis_round = character(0),
                        created_at_utc = character(0), path = character(0),
                        stringsAsFactors = FALSE)
    base <- file.path(root, "bundles")
    if (!dir.exists(base)) return(empty)

    omics_dirs <- list.dirs(base, recursive = FALSE)
    if (!is.null(omics)) omics_dirs <- omics_dirs[basename(omics_dirs) %in% omics]
    dirs <- unlist(lapply(omics_dirs, list.dirs, recursive = FALSE), use.names = FALSE)
    dirs <- dirs[!startsWith(basename(dirs), ".")]

    rows <- lapply(dirs, function(d) {
        # The completion marker is checked first, so a manifest still being
        # written is never parsed.
        if (!file.exists(file.path(d, "checksums.sha256"))) return(NULL)
        # One unreadable manifest must not hide every other bundle; it is
        # skipped here, and validate_omics_bundle() still refuses it.
        m <- tryCatch(.read_bundle_manifest(d), error = function(e) NULL)
        # Valid JSON that is not a manifest (a scalar, an array, or fields that
        # are not single strings) is skipped too, before any field is used.
        fields <- c("bundle_id", "omics", "project", "analysis_round",
                    "created_at_utc", "status")
        if (!is.list(m) || !all(vapply(fields, function(f) {
            is.character(m[[f]]) && length(m[[f]]) == 1
        }, logical(1)))) return(NULL)
        # A folder renamed away from its ID would always fail validation, and one
        # moved under another omics folder would be listed under the wrong omics.
        if (!identical(m$status, "complete") ||
            !identical(m$bundle_id, basename(d)) ||
            !identical(m$omics, basename(dirname(d)))) return(NULL)
        data.frame(bundle_id = m$bundle_id, omics = m$omics, project = m$project,
                   analysis_round = m$analysis_round, created_at_utc = m$created_at_utc,
                   path = d, stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, c(list(empty), rows))
    out <- out[order(out$omics, out$created_at_utc, out$bundle_id), , drop = FALSE]
    rownames(out) <- NULL
    out
}
