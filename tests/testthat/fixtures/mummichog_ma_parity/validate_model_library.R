# tests/testthat/fixtures/mummichog_ma_parity/validate_model_library.R
#
# One-off validation (run MANUALLY, never in CI; writes nothing): does our
# converter, mmc_ma_library_from_model(), build from a real model JSON the very
# library the PINNED MetaboAnalystR builder builds from the same compounds?
#
#   Rscript tests/testthat/fixtures/mummichog_ma_parity/validate_model_library.R \
#       /path/to/MetaboAnalystR-checkout /path/to/mng_kegg_20260723.json
#
# Checks, each reported PASS/FAIL (exit status 1 if any fails):
#   1. compound ids / masses carried over unchanged;
#   2. adduct matrices (positive, negative) identical to upstream's
#      Compound_function_mzlist() on the same masses;
#   3. m/z trees (positive, negative) identical to upstream's CreateLibFromKEGG();
#   4. known compounds: neutral mass and [M+H]+ / [M-H]- against textbook
#      values computed here from the elemental masses, independent of both
#      implementations;
#   5. pathway membership equal to an independent re-derivation from the
#      model's reactions (reactants + products), pathway by pathway.
# It also prints how many pathways survive upstream's minLib = 3 size filter.

PINNED_SHA <- "398476ae2a0c996e390925fa62adc46b5ecde334"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: validate_model_library.R <MetaboAnalystR checkout> <model.json>")
}
mar_dir    <- normalizePath(args[[1]], mustWork = TRUE)
model_path <- normalizePath(args[[2]], mustWork = TRUE)
here <- normalizePath(dirname(sub("^--file=", "", grep("^--file=",
          commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
if (is.na(here) || !nzchar(here)) here <- "."
repo_root <- normalizePath(file.path(here, "..", "..", "..", ".."), mustWork = TRUE)

sha <- system2("git", c("-C", mar_dir, "rev-parse", "HEAD"), stdout = TRUE)
if (!identical(sha, PINNED_SHA)) {
  stop("MetaboAnalystR checkout is at ", sha, ", not the pinned ", PINNED_SHA)
}

# ---- ours -------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x
for (f in c("06c_mummichog_pinned.R", "06f_mummichog_evidence.R",
            "06j_mummichog_ma_library.R")) {
  source(file.path(repo_root, "R", "domain", "metabolomics", f))
}
ours <- mmc_ma_library_from_model(model_path)

# ---- upstream, sourced verbatim, on the same compounds ----------------------
up <- new.env()
up$load_stringr <- function() invisible(NULL)
up$ov_qs_save   <- function(obj, file) assign("captured_lib", obj, envir = up)
sys.source(file.path(mar_dir, "R", "util_psea.R"), envir = up)

jm <- jsonlite::fromJSON(model_path, simplifyVector = FALSE)
cpd.lib <- list(
  id   = vapply(jm$list_of_compounds, function(c) c$id, character(1)),
  name = vapply(jm$list_of_compounds, function(c) c$name, character(1)),
  mw   = vapply(jm$list_of_compounds, function(c) as.numeric(c$neutral_mono_mass),
                numeric(1))
)
invisible(capture.output(up$CreateLibFromKEGG(cpd.lib, ours$pathways, "model")))
upl <- get("captured_lib", envir = up)

# ---- checks -----------------------------------------------------------------
results <- list()
check <- function(name, ok, detail = "") {
  results[[name]] <<- isTRUE(ok)
  cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", name,
              if (nzchar(detail)) paste0(" -- ", detail) else ""))
}

check("1a compound ids unchanged", identical(ours$cpd.lib$id, cpd.lib$id))
check("1b compound masses unchanged", identical(ours$cpd.lib$mw, cpd.lib$mw))
for (m in c("positive", "negative")) {
  check(paste0("2 adduct matrix identical: ", m),
        identical(ours$cpd.lib$adducts[[m]], upl$cpd.lib$adducts[[m]]))
  check(paste0("3 m/z tree identical: ", m),
        identical(ours$cpd.tree[[m]], upl$cpd.tree[[m]]))
}

# Textbook monoisotopic masses from elemental masses (IUPAC), independent of
# both implementations. Electron mass is accounted for via the proton mass.
el <- c(C = 12, H = 1.00782503207, N = 14.0030740048, O = 15.99491461956)
proton <- 1.007276466812
known <- list(
  C00031 = c(C = 6, H = 12, O = 6),       # D-glucose
  C00022 = c(C = 3, H = 4, O = 3),        # pyruvate
  C00025 = c(C = 5, H = 9, N = 1, O = 4)  # L-glutamate
)
for (cid in names(known)) {
  i <- match(cid, ours$cpd.lib$id)
  if (is.na(i)) { cat(sprintf("[SKIP] 4 %s not in this model\n", cid)); next }
  f <- known[[cid]]
  mono <- sum(el[names(f)] * f)
  mh  <- ours$cpd.lib$adducts$positive[i, "M+H[1+]"]
  mmh <- ours$cpd.lib$adducts$negative[i, "M-H[-]"]
  ok <- abs(ours$cpd.lib$mw[i] - mono) < 1e-4 &&
        abs(mh - (mono + proton)) < 1e-4 && abs(mmh - (mono - proton)) < 1e-4
  check(paste0("4 known compound ", cid), ok,
        sprintf("mono %.6f (expected %.6f), M+H %.6f, M-H %.6f",
                ours$cpd.lib$mw[i], mono, mh, mmh))
}

# Independent pathway membership: reactants + products of each reaction.
rx <- lapply(jm$list_of_reactions, function(r)
  unique(c(unlist(r$reactants), unlist(r$products))))
names(rx) <- vapply(jm$list_of_reactions, function(r) r$id, character(1))
indep <- lapply(jm$list_of_pathways, function(p) {
  ids <- unlist(p$list_of_reactions)
  intersect(unique(unlist(rx[ids])), cpd.lib$id)
})
names(indep) <- vapply(jm$list_of_pathways, function(p) p$name, character(1))
same <- vapply(seq_along(ours$pathways$name), function(k) {
  nm <- ours$pathways$name[[k]]
  setequal(ours$pathways$cpds[[k]], indep[[nm]])
}, logical(1))
check("5 pathway membership equals reaction-derived sets",
      length(same) == length(indep) && all(same),
      sprintf("%d of %d pathways equal", sum(same), length(indep)))

cat(sprintf("pathways with >= 3 compounds (upstream minLib = 3): %d of %d\n",
            sum(lengths(ours$pathways$cpds) >= 3), length(ours$pathways$cpds)))
cat("model:", ours$model$id, "version", ours$model$version, "\n")
cat("R:", R.version.string, "| MetaboAnalystR:", sha, "\n")
if (!all(unlist(results))) quit(status = 1)
