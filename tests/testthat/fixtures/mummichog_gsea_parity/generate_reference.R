# tests/testthat/fixtures/mummichog_gsea_parity/generate_reference.R
#
# Generates the GSEA parity reference by running the PINNED MetaboAnalystR
# implementation itself. Run MANUALLY, never in CI — the committed fixture is
# the reference, and regenerating it would defeat the point of a parity test.
#
#   Rscript tests/testthat/fixtures/mummichog_gsea_parity/generate_reference.R \
#       /path/to/MetaboAnalystR-checkout
#
# The checkout must be at the pinned commit recorded in REFERENCE.md.
#
# What runs upstream here:
#   * R/util_fgsea.R::.run_fgsea_inner()  — sourced verbatim from the pinned
#     checkout; this is the whole GSEA engine (ES, permutations, NES, p, BH,
#     leading edge). Only its two qs side-effect calls and AddErrMsg are stubbed,
#     because they are I/O for MetaboAnalyst's web session, not maths.
#   * the ranked-input construction of
#     R/peaks_to_function.R::.compute.mummichog.RT.fgsea() (pinned lines
#     3470-3483) — copied verbatim below, and cross-checked against our own
#     mmc_gsea_ranked_inputs() so a drift in either is caught.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: generate_reference.R <MetaboAnalystR checkout>")
mar_dir <- normalizePath(args[[1]], mustWork = TRUE)
here    <- normalizePath(dirname(sub("^--file=", "", grep("^--file=",
             commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
if (is.na(here) || !nzchar(here)) here <- "."
repo_root <- normalizePath(file.path(here, "..", "..", "..", ".."), mustWork = TRUE)

pinned_sha <- system2("git", c("-C", mar_dir, "rev-parse", "HEAD"), stdout = TRUE)

# The implementation depends on fgsea INTERNAL primitives, so the reference must
# be generated under the fgsea series this project locks. renv.lock pins the
# Bioconductor 3.22 release series (1.36.x); generate under that, and record the
# exact build so the parity test can refuse to compare across series.
lock <- jsonlite::fromJSON(file.path(repo_root, "renv.lock"))
locked_fgsea <- lock$Packages$fgsea$Version
fgsea_ver <- as.character(utils::packageVersion("fgsea"))
series <- function(v) paste(strsplit(v, ".", fixed = TRUE)[[1]][1:2], collapse = ".")
if (!identical(series(fgsea_ver), series(locked_fgsea))) {
  stop(sprintf(paste0("installed fgsea %s is not in renv.lock's series (%s). ",
                      "Install the locked series before generating the ",
                      "reference — the engine uses fgsea internals."),
               fgsea_ver, locked_fgsea))
}
message("fgsea ", fgsea_ver, " (renv.lock pins ", locked_fgsea, ")")

# ---- the fixture inputs ----------------------------------------------------
# Signed EC scores: positives and negatives, an exact tie (E02/E03), and a
# +/- pair with identical magnitude (E01/E16) to exercise the ordering.
ec_scores <- c(
  E01 =  4.0, E02 =  3.2, E03 =  3.2, E04 =  2.5, E05 =  1.8, E06 =  1.1,
  E07 =  0.6, E08 =  0.2, E09 = -0.3, E10 = -0.9, E11 = -1.4, E12 = -2.0,
  E13 = -2.6, E14 = -3.2, E15 = -3.9, E16 = -4.0
)
pathways <- list(
  PW_top       = c("E01", "E02", "E03", "E04"),   # incl. two ECs in one tie group
  PW_bottom    = c("E13", "E14", "E15", "E16"),
  PW_mid       = c("E07", "E08", "E09", "E10"),
  PW_overlap_a = c("E01", "E02", "E05", "E06", "E11"),
  PW_overlap_b = c("E02", "E03", "E04", "E05"),
  PW_single    = c("E05"),                        # size 1 (kept at minSize = 1)
  PW_absent    = c("EX1", "EX2")                  # no detected EC -> dropped
)
nperm     <- 100L   # PerformPSEA()'s default permNum
gsea_param <- 1
min_size   <- 1
max_size   <- Inf

# ---- upstream ranked-input construction (verbatim, pinned 3470-3483) -------
total_ecpds <- ec_scores
df.scores <- data.frame(id = names(total_ecpds), scores = total_ecpds)
ag.scores <- aggregate(id ~ scores, data = df.scores, paste, collapse = "; ")
ag.sorted <- ag.scores[order(-ag.scores$scores), ]
row.names(ag.sorted) <- NULL
dt.scores <- data.table::data.table(ag.sorted)
dt.scores.out <- dt.scores[, list(scores = scores,
                                  id = unlist(strsplit(id, "; ", fixed = TRUE))),
                           by = 1:nrow(dt.scores)]
rank.vec <- as.numeric(dt.scores.out$nrow)
names(rank.vec) <- as.character(dt.scores.out$id)
scores.vec <- as.numeric(ag.sorted$scores)
names(scores.vec) <- as.character(ag.sorted$id)

# ---- cross-check our own construction against it ---------------------------
for (f in list.files(file.path(repo_root, "R"), pattern = "\\.[Rr]$",
                     recursive = TRUE, full.names = TRUE)) {
  try(source(f), silent = TRUE)
}
ours <- mmc_gsea_ranked_inputs(ec_scores)
stopifnot(identical(names(ours$stats), names(scores.vec)),
          isTRUE(all.equal(unname(ours$stats), unname(scores.vec))),
          identical(names(ours$ranks), names(rank.vec)),
          isTRUE(all.equal(unname(ours$ranks), unname(rank.vec))))
message("ranked-input construction matches upstream verbatim block")

# ---- run the PINNED upstream engine ---------------------------------------
ov_qs_save <- function(...) invisible(NULL)
ov_qs_read <- function(...) data.frame()
AddErrMsg  <- function(...) invisible(NULL)
suppressPackageStartupMessages({
  library(fgsea); library(BiocParallel); library(fastmatch); library(data.table)
})
source(file.path(mar_dir, "R", "util_fgsea.R"))

mSetObj <- list(dataSet = list(paramSet = list(mumRT = TRUE)))
upstream_warnings <- character(0)
ref <- withCallingHandlers(
  .run_fgsea_inner(mSetObj, pathways, scores.vec, rank.vec, nperm,
                   minSize = min_size, maxSize = max_size,
                   gseaParam = gsea_param),
  warning = function(w) {
    upstream_warnings <<- c(upstream_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

ref_df <- data.frame(
  pathway            = as.character(ref$pathway),
  pval               = as.numeric(ref$pval),
  padj               = as.numeric(ref$padj),
  ES                 = as.numeric(ref$ES),
  NES                = as.numeric(ref$NES),
  nMoreExtreme       = as.numeric(ref$nMoreExtreme),
  size               = as.numeric(ref$size),
  leadingEdgeMatched = as.character(ref$leadingEdgeMatched),
  stringsAsFactors   = FALSE
)

# fastmatch::fmatch caches a hash table as a `.match.hash` attribute on
# names(ranks) while the engine runs. That is an implementation artefact (an
# external pointer), not part of the reference, so strip it before storing.
strip_attrs <- function(x) {
  nm <- names(x)
  attributes(nm) <- NULL
  stats::setNames(as.numeric(unname(x)), nm)
}
scores.vec <- strip_attrs(scores.vec)
rank.vec   <- strip_attrs(rank.vec)

out <- list(
  metaboanalystr_commit = pinned_sha,
  metaboanalystr_files  = c("R/util_fgsea.R", "R/peaks_to_function.R"),
  upstream_entry        = ".compute.mummichog.RT.fgsea -> fgsea2 -> my.fgsea -> .run_fgsea_inner",
  fgsea_version         = fgsea_ver,
  fgsea_series          = series(fgsea_ver),
  fgsea_locked          = locked_fgsea,
  fgsea_branch          = if (utils::packageVersion("fgsea") > "1.24.0")
                            "post-1.24.0 (stats re-sorted decreasing after abs())"
                          else "pre-1.24.0 (no re-sort)",
  r_version             = R.version.string,
  bpparam               = class(BiocParallel::bpparam())[1],
  nperm = nperm, gsea_param = gsea_param,
  min_size = min_size, max_size = max_size,
  set_seed = 123L,
  ec_scores = ec_scores, pathways = pathways,
  stats = scores.vec, ranks = rank.vec,
  reference = ref_df,
  upstream_warnings = upstream_warnings,
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
saveRDS(out, file.path(here, "reference.rds"), version = 2)
readr::write_tsv(ref_df, file.path(here, "reference.tsv"))
cat("\n=== pinned MetaboAnalystR reference ===\n")
cat("commit:", pinned_sha, "\nfgsea:", out$fgsea_version, "-", out$fgsea_branch, "\n")
cat("upstream warnings:", length(upstream_warnings), "\n\n")
print(ref_df)
