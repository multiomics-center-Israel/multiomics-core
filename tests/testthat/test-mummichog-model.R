# tests/testthat/test-mummichog-model.R
#
# Unit tests for the mummichog model-selection layer (06d): sha256-verified
# fetch + cache (mmc_resolve_model) and the model_ref > model_json > human_mfn
# precedence with the human-only organism guard (mmc_select_model).
#
# All tests are NETWORK-FREE: mmc_resolve_model takes an injectable downloader,
# and the model_ref precedence cases are served from a pre-warmed cache so the
# default curl downloader is never even forced.

# Write `content` with writeLines into a scratch file and return its sha256,
# computed exactly as the code under test computes it. Any downloader that
# writes the same `content` with writeLines yields a file with this digest.
sha_of_written <- function(content) {
  scratch <- withr::local_tempfile(.local_envir = parent.frame())
  writeLines(content, scratch)
  .mmc_sha256_file(scratch)
}

# A downloader that writes fixed content (mimics a successful fetch).
writer <- function(content) function(url, destfile) writeLines(content, destfile)

test_that("mmc_resolve_model verifies sha256 and caches under <sha256>.json", {
  cache   <- withr::local_tempdir()
  content <- '{"model":"ok"}'
  sha     <- sha_of_written(content)

  path <- mmc_resolve_model(list(url = "https://x/m.json", sha256 = sha),
                            cache_dir = cache, downloader = writer(content))

  expect_true(file.exists(path))
  expect_identical(.mmc_sha256_file(path), sha)
  expect_match(path, paste0(sha, "\\.json$"))
})

test_that("mmc_resolve_model errors on a sha256 mismatch and caches nothing", {
  cache     <- withr::local_tempdir()
  wrong_sha <- paste(rep("a", 64), collapse = "")

  expect_error(
    mmc_resolve_model(list(url = "https://x/m.json", sha256 = wrong_sha),
                      cache_dir = cache, downloader = writer("anything")),
    "mismatch"
  )
  # a failed verification must not leave a (mis-named, unverified) cache entry
  expect_length(list.files(cache), 0L)
})

test_that("mmc_resolve_model reuses a verified cache entry without downloading", {
  cache   <- withr::local_tempdir()
  content <- '{"model":"cached"}'
  sha     <- sha_of_written(content)
  writeLines(content, file.path(cache, paste0(sha, ".json")))

  boom <- function(url, destfile) stop("must not download on a cache hit")
  path <- mmc_resolve_model(list(url = "https://x/m.json", sha256 = sha),
                            cache_dir = cache, downloader = boom)

  expect_match(path, paste0(sha, "\\.json$"))
})

test_that("mmc_resolve_model re-fetches when the cached file is corrupt", {
  cache   <- withr::local_tempdir()
  content <- '{"model":"good"}'
  sha     <- sha_of_written(content)
  # a corrupt entry squatting at the expected cache path
  writeLines("corrupt", file.path(cache, paste0(sha, ".json")))

  path <- mmc_resolve_model(list(url = "https://x/m.json", sha256 = sha),
                            cache_dir = cache, downloader = writer(content))

  expect_identical(.mmc_sha256_file(path), sha)   # now holds the good content
})

test_that("mmc_resolve_model rejects a malformed model_ref", {
  cache <- withr::local_tempdir()
  dl    <- writer("x")
  hex64 <- paste(rep("a", 64), collapse = "")

  expect_error(mmc_resolve_model("not-a-list", cache_dir = cache, downloader = dl),
               "mapping")
  expect_error(mmc_resolve_model(list(url = "", sha256 = hex64),
                                 cache_dir = cache, downloader = dl), "url")
  expect_error(mmc_resolve_model(list(url = "https://x/m.json", sha256 = "not-hex"),
                                 cache_dir = cache, downloader = dl), "hex")
})

test_that("mmc_select_model precedence: model_ref > model_json > human_mfn", {
  cache <- withr::local_tempdir()

  # model_json beats the built-in fallback (returned as-is, non-human organism)
  expect_identical(
    mmc_select_model(list(model_json = "/custom/model.json"),
                     organism = "Mus musculus", cache_dir = cache),
    "/custom/model.json"
  )

  # no model + human organism -> the built-in model
  expect_identical(
    mmc_select_model(list(), organism = "Homo sapiens", cache_dir = cache),
    "human_mfn"
  )

  # model_ref beats model_json (served from a pre-warmed cache -> no network)
  content <- '{"model":"ref-wins"}'
  sha     <- sha_of_written(content)
  writeLines(content, file.path(cache, paste0(sha, ".json")))
  out <- mmc_select_model(
    list(model_ref  = list(url = "https://x/m.json", sha256 = sha),
         model_json = "/ignored.json"),
    organism = "Mus musculus", cache_dir = cache
  )
  expect_match(out, paste0(sha, "\\.json$"))
})

test_that("mmc_select_model refuses a non-human organism with no custom model", {
  expect_error(
    mmc_select_model(list(), organism = "Mus musculus",
                     cache_dir = withr::local_tempdir()),
    "non-human"
  )
  # the legacy mummichog.organisms key is honoured too
  expect_error(
    mmc_select_model(list(organisms = "mmu"), organism = NULL,
                     cache_dir = withr::local_tempdir()),
    "non-human"
  )
})

test_that("model_json: human_mfn on a non-human organism is still guarded", {
  # The guard keys on the RESOLVED model, so naming the built-in model
  # explicitly via model_json must not smuggle it past the check.
  expect_error(
    mmc_select_model(list(model_json = "human_mfn"), organism = "Mus musculus",
                     cache_dir = withr::local_tempdir()),
    "non-human"
  )
  # ...but human_mfn on a human organism is fine
  expect_identical(
    mmc_select_model(list(model_json = "human_mfn"), organism = "Homo sapiens",
                     cache_dir = withr::local_tempdir()),
    "human_mfn"
  )
})

test_that("a model_ref supplies a real model, so a non-human organism passes", {
  skip_if_not_installed("yaml")
  cache   <- withr::local_tempdir()
  content <- '{"model":"cre-kegg"}'
  sha     <- sha_of_written(content)
  writeLines(content, file.path(cache, paste0(sha, ".json")))

  # a config block as it would arrive from the YAML loader
  yml <- sprintf(paste(
    "mummichog:",
    "  enabled: true",
    "  model_ref:",
    "    url: https://example.org/models/cre_kegg_20260711.json",
    "    sha256: %s",
    sep = "\n"), sha)
  cfg <- yaml::yaml.load(yml)

  out <- mmc_select_model(cfg$mummichog, organism = "Caenorhabditis elegans",
                          cache_dir = cache)
  expect_match(out, paste0(sha, "\\.json$"))
})
