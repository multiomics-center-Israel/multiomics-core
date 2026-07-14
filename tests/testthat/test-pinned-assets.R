# tests/testthat/test-pinned-assets.R
#
# Unit tests for the generic pinned-asset resolver (core/01b): sha256-verified
# URL download + cache, and the mandatory-checksum local_path override.
# NETWORK-FREE: resolve_pinned_asset takes an injectable downloader.

# sha256 of writeLines(content) — matches how the code hashes file bytes.
sha_of_written <- function(content) {
  scratch <- withr::local_tempfile(.local_envir = parent.frame())
  writeLines(content, scratch)
  .sha256_file(scratch)
}

# A downloader that writes fixed content (mimics a successful fetch).
writer <- function(content) function(url, destfile) writeLines(content, destfile)

test_that("resolve_pinned_asset verifies + caches a URL download under <sha><ext>", {
  cache   <- withr::local_tempdir()
  content <- "reaction_id\tsubstrate_id\tproduct_id"
  sha     <- sha_of_written(content)

  path <- resolve_pinned_asset(
    list(url = "https://x/ref.tsv.gz", sha256 = sha),
    cache_dir = cache, ext = ".tsv.gz", downloader = writer(content))

  expect_true(file.exists(path))
  expect_match(path, paste0(sha, "\\.tsv\\.gz$"))
  expect_identical(.sha256_file(path), sha)
})

test_that("resolve_pinned_asset errors on a URL sha mismatch and caches nothing", {
  cache <- withr::local_tempdir()
  wrong <- paste(rep("a", 64), collapse = "")

  expect_error(
    resolve_pinned_asset(list(url = "https://x/ref", sha256 = wrong),
                         cache_dir = cache, ext = ".tsv.gz",
                         downloader = writer("whatever")),
    "mismatch")
  expect_length(list.files(cache), 0L)
})

test_that("resolve_pinned_asset reuses a verified cache entry without downloading", {
  cache   <- withr::local_tempdir()
  content <- "cached"
  sha     <- sha_of_written(content)
  writeLines(content, file.path(cache, paste0(sha, ".tsv.gz")))

  boom <- function(url, destfile) stop("must not download on a cache hit")
  path <- resolve_pinned_asset(list(url = "https://x/ref", sha256 = sha),
                               cache_dir = cache, ext = ".tsv.gz",
                               downloader = boom)
  expect_match(path, paste0(sha, "\\.tsv\\.gz$"))
})

test_that("resolve_pinned_asset accepts a checksum-valid local_path fully offline", {
  f <- withr::local_tempfile(fileext = ".tsv.gz")
  writeLines("id\tval", f)
  sha  <- .sha256_file(f)
  boom <- function(url, destfile) stop("must not download for a local_path")

  path <- resolve_pinned_asset(list(local_path = f, sha256 = sha),
                               cache_dir = withr::local_tempdir(), ext = ".tsv.gz",
                               downloader = boom)
  expect_identical(normalizePath(path), normalizePath(f))
})

test_that("resolve_pinned_asset rejects a local_path whose checksum is wrong", {
  f <- withr::local_tempfile(fileext = ".tsv.gz")
  writeLines("id\tval", f)
  wrong <- paste(rep("b", 64), collapse = "")

  expect_error(
    resolve_pinned_asset(list(local_path = f, sha256 = wrong),
                         cache_dir = withr::local_tempdir(), ext = ".tsv.gz",
                         downloader = writer("x")),
    "mismatch")
})

test_that("resolve_pinned_asset errors on a missing local_path", {
  wrong <- paste(rep("c", 64), collapse = "")
  expect_error(
    resolve_pinned_asset(list(local_path = "/no/such/ref.tsv.gz", sha256 = wrong),
                         cache_dir = withr::local_tempdir(), ext = ".tsv.gz"),
    "does not exist")
})

test_that("resolve_pinned_asset rejects malformed references", {
  cache <- withr::local_tempdir()
  hex64 <- paste(rep("a", 64), collapse = "")

  expect_error(resolve_pinned_asset("not-a-list", cache_dir = cache), "mapping")
  expect_error(resolve_pinned_asset(list(sha256 = "not-hex"), cache_dir = cache),
               "hex")
  # sha ok but neither url nor local_path -> complains about url
  expect_error(resolve_pinned_asset(list(sha256 = hex64), cache_dir = cache),
               "url")
})
