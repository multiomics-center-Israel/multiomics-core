# R/core/01b_pinned_assets.R
#
# Generic pinned-asset resolver: turn a {url | local_path, sha256} reference
# into a verified local file path. A pinned asset is fetched (or read from a
# pre-provisioned local copy) exactly once, checked against its sha256, and
# cached under <cache_dir>/<sha256><ext>. This is the SINGLE implementation
# shared by the mummichog model layer (domain/metabolomics/06d) and the
# metabolite-network reference (domain/metabolomics/08b).
#
# A checksum is MANDATORY in every path, including local_path: the local
# override is a pre-provisioned copy of the SAME pinned asset, not an arbitrary
# unversioned file. Never a silent "use what we got" -- a mismatch is a hard
# error. No private/authenticated download here (deferred pending the
# distribution decision): only a public URL or a local file.
#
# Dependencies: digest (sha256), curl (download) -- both already in renv.lock.


#' sha256 digest of a file's contents (lowercase hex)
#'
#' @param path Path to the file to hash.
#' @return A 64-character lowercase hex string.
#' @noRd
.sha256_file <- function(path) {
  tolower(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}


#' Resolve a pinned asset reference to a verified, cached local file
#'
#' Accepts a mapping with `sha256` (required) and either `local_path` (a
#' pre-provisioned copy, verified against the checksum and used offline) or
#' `url` (fetched to a scratch file, verified, then published into the cache).
#' The cache is keyed on the expected digest (`<sha256><ext>`), so a
#' present-but-corrupt entry is a miss and is re-fetched rather than trusted.
#'
#' @param ref A list with `sha256` (single 64-char hex string) and either
#'   `local_path` or `url` (both single non-empty strings). `local_path` wins
#'   when both are set.
#' @param cache_dir Directory for the on-disk cache of downloaded assets.
#'   Created if absent.
#' @param ext File extension for cached downloads (e.g. ".json", ".tsv.gz").
#' @param downloader Function `(url, destfile)` that fetches `url` to
#'   `destfile`. Defaults to `curl::curl_download`; injectable so tests need no
#'   network.
#' @return The path to the verified file (normalised, winslash = "/").
resolve_pinned_asset <- function(ref, cache_dir, ext = ".bin",
                                 downloader = curl::curl_download) {
  if (!is.list(ref)) {
    stop("pinned asset reference must be a mapping with 'sha256' and either ",
         "'url' or 'local_path'.", call. = FALSE)
  }

  # -- sha256 is mandatory in every path (including local_path) --------------
  sha256 <- ref$sha256
  if (is.null(sha256) || !is.character(sha256) || length(sha256) != 1L) {
    stop("pinned asset 'sha256' must be a single 64-character hex digest.",
         call. = FALSE)
  }
  sha256 <- tolower(trimws(sha256))
  if (!grepl("^[0-9a-f]{64}$", sha256)) {
    stop("pinned asset 'sha256' must be a 64-character hex digest; got '",
         sha256, "'.", call. = FALSE)
  }

  # -- local override: verify existence + checksum, stay fully offline -------
  local_path <- ref$local_path
  if (!is.null(local_path) && is.character(local_path) &&
      length(local_path) == 1L && nzchar(local_path)) {
    if (!file.exists(local_path)) {
      stop("pinned asset local_path does not exist: '", local_path,
           "'. Provision the checksummed reference there, or configure a 'url'.",
           call. = FALSE)
    }
    got <- .sha256_file(local_path)
    if (!identical(got, sha256)) {
      stop("sha256 mismatch for local_path '", local_path, "': expected ",
           sha256, " but got ", got,
           ". The local copy is not the pinned asset; refusing to use it.",
           call. = FALSE)
    }
    return(normalizePath(local_path, winslash = "/", mustWork = TRUE))
  }

  # -- otherwise a public URL is required ------------------------------------
  url <- ref$url
  if (is.null(url) || !is.character(url) || length(url) != 1L || !nzchar(url)) {
    stop("pinned asset must supply a non-empty 'url' (or a 'local_path').",
         call. = FALSE)
  }

  cache_file <- file.path(cache_dir, paste0(sha256, ext))
  if (file.exists(cache_file) && identical(.sha256_file(cache_file), sha256)) {
    return(normalizePath(cache_file, winslash = "/", mustWork = TRUE))
  }

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  tmp <- tempfile("pinned-asset-", fileext = ext)
  on.exit(unlink(tmp, force = TRUE), add = TRUE)

  tryCatch(
    downloader(url, tmp),
    error = function(e) {
      stop("failed to download pinned asset from '", url, "': ",
           conditionMessage(e), call. = FALSE)
    }
  )
  if (!file.exists(tmp) || file.info(tmp)$size == 0) {
    stop("download of '", url, "' produced no data.", call. = FALSE)
  }

  got <- .sha256_file(tmp)
  if (!identical(got, sha256)) {
    stop("sha256 mismatch for asset downloaded from '", url, "': expected ",
         sha256, " but got ", got, ". Refusing to use an unverified asset.",
         call. = FALSE)
  }

  # Publish atomically where possible; file.rename can fail across filesystems
  # (scratch tempdir vs the cache dir), so fall back to copy.
  if (!file.rename(tmp, cache_file)) {
    if (!file.copy(tmp, cache_file, overwrite = TRUE)) {
      stop("could not write the verified asset into the cache at '",
           cache_file, "'.", call. = FALSE)
    }
  }
  normalizePath(cache_file, winslash = "/", mustWork = TRUE)
}
