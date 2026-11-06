# tests/testthat/fixtures/mummichog_ma_parity/generate_library_reference.R
#
# Generates the LIBRARY parity reference for the MetaboAnalyst-compatible
# mummichog engine by running the PINNED MetaboAnalystR library builder itself.
# Run MANUALLY, never in CI — the committed fixture is the reference.
#
#   Rscript tests/testthat/fixtures/mummichog_ma_parity/generate_library_reference.R \
#       /path/to/MetaboAnalystR-checkout
#
# The checkout must be at commit 398476ae2a0c996e390925fa62adc46b5ecde334.
#
# What runs upstream here: R/util_psea.R, sourced verbatim into its own
# environment — CreateLibFromKEGG() and Compound_function_mzlist(). Only two
# side-effect helpers are stubbed: load_stringr() (a package loader the maths
# never uses) and ov_qs_save() (the .qs file write, replaced by capturing the
# object it would have written).
#
# Output: library_reference.rds next to this script.

PINNED_SHA <- "398476ae2a0c996e390925fa62adc46b5ecde334"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: generate_library_reference.R <MetaboAnalystR checkout>")
mar_dir <- normalizePath(args[[1]], mustWork = TRUE)
here <- normalizePath(dirname(sub("^--file=", "", grep("^--file=",
          commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
if (is.na(here) || !nzchar(here)) here <- "."

sha <- system2("git", c("-C", mar_dir, "rev-parse", "HEAD"), stdout = TRUE)
if (!identical(sha, PINNED_SHA)) {
  stop("MetaboAnalystR checkout is at ", sha, ", not the pinned ", PINNED_SHA)
}

# ---- upstream, sourced verbatim ---------------------------------------------
up <- new.env()
up$load_stringr <- function() invisible(NULL)
up$ov_qs_save   <- function(obj, file) assign("captured_lib", obj, envir = up)
sys.source(file.path(mar_dir, "R", "util_psea.R"), envir = up)

# ---- synthetic library inputs -----------------------------------------------
# Masses chosen to exercise the tree edges, not to be real metabolites:
#   C90001 / C90002  ordinary masses;
#   C90003           tiny mass: most adducts are negative or < 50 (never binned);
#   C90004           M+H floors to 50, which upstream EXCLUDES (n > 50);
#   C90005           near 2000: some adducts bin at 1999, the rest are >= 2000;
#   C90006           several adducts floor into one bin (exercises unique()).
cpd.lib <- list(
  id   = sprintf("C9000%d", 1:6),
  name = sprintf("synthetic compound %d", 1:6),
  mw   = c(180.063388, 146.105528, 1.007825, 49.9, 1998.5, 300.9995)
)
pathways <- list(
  cpds = list(c("C90001", "C90002", "C90006"), c("C90003", "C90004"), "C90005"),
  id   = c("syn00001", "syn00002", "syn00003"),
  name = c("Synthetic pathway A", "Synthetic pathway B", "Synthetic pathway C")
)

up$CreateLibFromKEGG(cpd.lib, pathways, "synthetic")
lib <- get("captured_lib", envir = up)

ref <- list(
  metaboanalystr_commit = sha,
  upstream_files        = "R/util_psea.R",
  upstream_functions    = c("CreateLibFromKEGG", "Compound_function_mzlist"),
  r_version             = R.version.string,
  generated_utc         = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  inputs                = list(cpd.lib = cpd.lib, pathways = pathways),
  adducts               = lib$cpd.lib$adducts[c("positive", "negative")],
  cpd.tree              = lib$cpd.tree[c("positive", "negative")]
)
out <- file.path(here, "library_reference.rds")
saveRDS(ref, out, version = 2)
message("wrote ", out)
str(ref$adducts, max.level = 1)
