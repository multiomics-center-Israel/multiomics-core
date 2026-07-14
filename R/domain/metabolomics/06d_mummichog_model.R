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
#' Thin alias over the generic `.sha256_file()` (core/01b), kept because the
#' mummichog model tests reference this name directly.
#'
#' @param path Path to the file to hash.
#' @return A 64-character lowercase hex string.
#' @noRd
.mmc_sha256_file <- function(path) .sha256_file(path)

#' Resolve a URL+sha256 model reference to a verified, cached local file
#'
#' Thin wrapper over the generic `resolve_pinned_asset()` (core/01b) — the
#' single URL/cache/sha256 implementation — pinning the cache extension to
#' ".json" for mummichog models. `model_ref` is a mapping with `url` and
#' `sha256`; the download is verified against the checksum and cached under
#' `<cache_dir>/<sha256>.json`, and a mismatch is a hard error.
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
  resolve_pinned_asset(model_ref, cache_dir = cache_dir, ext = ".json",
                       downloader = downloader)
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
  # Resolve the model by precedence:
  #   1. model_ref  — fetch-and-verify a pinned remote model (cached path).
  #   2. model_json — a local custom model, passed through as-is; run_mummichog()
  #      normalises it if it exists (a non-file value is left for mummichog to
  #      reject), matching the pre-model_ref behaviour.
  #   3. built-in "human_mfn".
  network <- if (!is.null(mummi_cfg$model_ref)) {
    mmc_resolve_model(mummi_cfg$model_ref, cache_dir = cache_dir)
  } else if (!is.null(mummi_cfg$model_json)) {
    mummi_cfg$model_json
  } else {
    "human_mfn"
  }

  # The built-in model is human-only. Refuse to run it on a non-human organism
  # whether it was selected by default OR named explicitly (e.g.
  # model_json: human_mfn) — a species mismatch here would otherwise pass with
  # no warning. Guard on the RESOLVED model, not the config source.
  if (identical(network, "human_mfn")) {
    org <- organism %||% mummi_cfg$organisms
    org <- tolower(trimws(paste(unlist(org), collapse = " ")))
    if (nzchar(org) && !grepl("human|homo sapiens|\\bhsa\\b", org)) {
      .mmc_stop(
        "configured organism ('", org, "') is non-human, but the pinned engine's ",
        "built-in model is human-only. Supply a custom model via 'model_ref' ",
        "(URL + sha256) or 'model_json' (a local path), or set the organism to human."
      )
    }
  }
  network
}
