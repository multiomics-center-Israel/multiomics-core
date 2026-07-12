# R/domain/metabolomics/06d_mummichog_model.R
#
# Model-selection layer for the pinned mummichog v2 engine (06c). Turns the
# config's mummichog block into the concrete `-n <model>` argument the runner
# needs, honouring this precedence:
#
#     model_ref  (URL + sha256, fetched & verified, cached on disk)
#   > model_json (a local path to a custom metabolic model)
#   > "human_mfn" (mummichog's built-in human model)
#
# `model_ref` lets a project pin an organism-specific model published elsewhere
# (e.g. the cre_kegg_* artifacts) without committing large JSON into this repo:
# the file is downloaded once, verified against its sha256, and cached under
# envs/mummichog-models/<sha256>.json (a gitignored dir). Subsequent runs reuse
# the cache as long as its content still matches the digest.
#
# Dependencies (R): digest (sha256), curl (download). Both are light and already
# in renv.lock — no Bioconductor.


#' sha256 digest of a file's contents (lowercase hex)
#'
#' @param path Path to the file to hash.
#' @return A 64-character lowercase hex string.
#' @noRd
.mmc_sha256_file <- function(path) {
  tolower(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

#' Resolve a URL+sha256 model reference to a verified, cached local file
#'
#' Downloads the model referenced by `model_ref` (a mapping with `url` and
#' `sha256`) unless a verified copy is already cached, checks its sha256, and
#' returns the path to the cached file. The download is fetched to a scratch
#' file and only published into the cache once its digest matches — a mismatch
#' is a hard error, never a silent "use what we got". The cache is keyed on the
#' expected digest (`<sha256>.json`), so a present-but-corrupt entry is treated
#' as a miss and re-fetched rather than trusted.
#'
#' @param model_ref  A list with `url` (single non-empty string) and `sha256`
#'                   (single 64-char hex string).
#' @param cache_dir  Directory for the on-disk model cache. Created if absent.
#' @param downloader Function `(url, destfile)` that fetches `url` to `destfile`.
#'                   Defaults to `curl::curl_download`; injectable so tests need
#'                   no network.
#' @return The path to the verified cached model (normalised, winslash = "/").
#' @examples
#' \dontrun{
#' ref <- list(
#'   url    = "https://example.org/models/cre_kegg_20260711.json",
#'   sha256 = "c403c96fbec8df9ae34b828fec01270c8ea3940acc36e4e5ff770868dc8b912b"
#' )
#' model_path <- mmc_resolve_model(ref)
#' }
mmc_resolve_model <- function(model_ref,
                              cache_dir = "envs/mummichog-models",
                              downloader = curl::curl_download) {
  # -- validate the model_ref block ------------------------------------------
  if (!is.list(model_ref)) {
    .mmc_stop("model_ref must be a mapping with 'url' and 'sha256'.")
  }
  url    <- model_ref$url
  sha256 <- model_ref$sha256
  if (is.null(url) || !is.character(url) || length(url) != 1L || !nzchar(url)) {
    .mmc_stop("model_ref$url must be a single non-empty string.")
  }
  if (is.null(sha256) || !is.character(sha256) || length(sha256) != 1L) {
    .mmc_stop("model_ref$sha256 must be a single 64-character hex sha256 digest.")
  }
  sha256 <- tolower(trimws(sha256))
  if (!grepl("^[0-9a-f]{64}$", sha256)) {
    .mmc_stop("model_ref$sha256 must be a 64-character hex sha256 digest; got '",
              sha256, "'.")
  }

  # -- cache hit: an existing file whose content still hashes as expected -----
  cache_file <- file.path(cache_dir, paste0(sha256, ".json"))
  if (file.exists(cache_file) &&
      identical(.mmc_sha256_file(cache_file), sha256)) {
    return(normalizePath(cache_file, winslash = "/", mustWork = TRUE))
  }

  # -- download to a scratch file, verify, then publish into the cache --------
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  tmp <- tempfile("mmc-model-", fileext = ".json")
  on.exit(unlink(tmp, force = TRUE), add = TRUE)

  tryCatch(
    downloader(url, tmp),
    error = function(e) {
      .mmc_stop("failed to download model from '", url, "': ",
                conditionMessage(e))
    }
  )
  if (!file.exists(tmp) || file.info(tmp)$size == 0) {
    .mmc_stop("download of '", url, "' produced no data.")
  }

  got <- .mmc_sha256_file(tmp)
  if (!identical(got, sha256)) {
    .mmc_stop("sha256 mismatch for model downloaded from '", url, "': expected ",
              sha256, " but got ", got,
              ". Refusing to use an unverified model.")
  }

  # Publish atomically where possible. file.rename can fail across filesystems
  # (scratch tempdir vs the project's cache dir); fall back to copy + cleanup.
  if (!file.rename(tmp, cache_file)) {
    if (!file.copy(tmp, cache_file, overwrite = TRUE)) {
      .mmc_stop("could not write the verified model into the cache at '",
                cache_file, "'.")
    }
  }
  normalizePath(cache_file, winslash = "/", mustWork = TRUE)
}

#' Select the mummichog `-n` model from config, honouring model precedence
#'
#' Resolves the metabolic model for a mummichog run from the config's mummichog
#' block: a URL+sha256 `model_ref` wins, then a local `model_json` path, then
#' the built-in `"human_mfn"`. When falling back to the built-in model, the
#' organism is checked — a non-human organism with no custom model is refused
#' rather than silently analysed against the human network (a species mismatch
#' that would otherwise pass unnoticed).
#'
#' @param mummi_cfg  The `modes.metabolomics.enrichment.mummichog` config list
#'                   (may carry `model_ref`, `model_json`, `organisms`).
#' @param organism   The general `modes.metabolomics.organism` value, used only
#'                   for the human-only guard on the built-in model.
#' @param cache_dir  Directory for the model cache (passed to
#'                   `mmc_resolve_model()`).
#' @return A model spec for `run_mummichog(network = )`: a verified cached path
#'   (model_ref), the configured local path (model_json), or `"human_mfn"`.
mmc_select_model <- function(mummi_cfg, organism = NULL,
                             cache_dir = "envs/mummichog-models") {
  # 1. model_ref: fetch-and-verify a pinned remote model (highest precedence).
  if (!is.null(mummi_cfg$model_ref)) {
    return(mmc_resolve_model(mummi_cfg$model_ref, cache_dir = cache_dir))
  }

  # 2. model_json: a local custom model. Passed through as-is; run_mummichog()
  #    normalises it if it exists (a non-file value is left for mummichog to
  #    reject), matching the pre-model_ref behaviour.
  if (!is.null(mummi_cfg$model_json)) {
    return(mummi_cfg$model_json)
  }

  # 3. built-in human_mfn — refuse to run it on a non-human organism.
  org <- organism %||% mummi_cfg$organisms
  org <- tolower(trimws(paste(unlist(org), collapse = " ")))
  if (nzchar(org) && !grepl("human|homo sapiens|\\bhsa\\b", org)) {
    .mmc_stop(
      "configured organism ('", org, "') is non-human, but the pinned engine's ",
      "built-in model is human-only. Supply a custom model via 'model_ref' ",
      "(URL + sha256) or 'model_json' (a local path), or set the organism to human."
    )
  }
  "human_mfn"
}
