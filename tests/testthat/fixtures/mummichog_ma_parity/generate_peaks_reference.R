# tests/testthat/fixtures/mummichog_ma_parity/generate_peaks_reference.R
#
# Generates the PEAKS parity reference (component 2: peak-list preparation,
# significant-peak selection and polarity-specific m/z matching) by running the
# PINNED MetaboAnalystR code itself. Run MANUALLY, never in CI.
#
#   Rscript tests/testthat/fixtures/mummichog_ma_parity/generate_peaks_reference.R \
#       /path/to/MetaboAnalystR-checkout
#
# The checkout must be at commit 398476ae2a0c996e390925fa62adc46b5ecde334.
#
# What runs upstream, verbatim (R/peaks_to_function.R and R/util_psea.R sourced
# into one environment):
#   * SanityCheckMummichogData()  — ordering, range trim, de-duplication, m/z
#                                   keys, RT tolerance;
#   * SetMummichogPval()          — significant peaks;
#   * CreateLibFromKEGG()         — the library (validated in component 1), then
#                                   the cpd.lib list exactly as
#                                   .setup.psea.library() builds it;
#   * .search.compoundLib()       — its source lines from the function start up
#                                   to (not including) the line
#                                   "# now create empirical compounds if
#                                   necessary!", i.e. the matching and currency
#                                   removal, closed with a return of the match
#                                   table and the RT ranks. Nothing else in the
#                                   function is altered.
# Stubbed (session I/O only, no maths): .get.mSet/.set.mSet (identity),
# AddErrMsg (stop), ov_qs_save/ov_qs_read (an in-memory store), save.image,
# load_stringr.
#
# Not run upstream: Read.PeakListData(), which only reads a file. The mSet it
# would produce is built directly below with its own lines
# (mummi.orig <- cbind(p.value, m.z, t.score, r.t); pos_inx <- mode == "positive").
#
# Runs the fixture in all three modes: "mixed", "positive", "negative".
# Output: peaks_reference.rds next to this script.

PINNED_SHA <- "398476ae2a0c996e390925fa62adc46b5ecde334"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: generate_peaks_reference.R <MetaboAnalystR checkout>")
mar_dir <- normalizePath(args[[1]], mustWork = TRUE)
here <- normalizePath(dirname(sub("^--file=", "", grep("^--file=",
          commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
if (is.na(here) || !nzchar(here)) here <- "."

sha <- system2("git", c("-C", mar_dir, "rev-parse", "HEAD"), stdout = TRUE)
if (!identical(sha, PINNED_SHA)) {
  stop("MetaboAnalystR checkout is at ", sha, ", not the pinned ", PINNED_SHA)
}
if (!identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection"))) {
  stop("generate under R's default RNG kinds; got ", paste(RNGkind(), collapse = "/"))
}

source(file.path(here, "fixture_inputs.R"))
inp <- ma_fixture_inputs()

# ---- upstream environment -----------------------------------------------------
new_upstream <- function() {
  up <- new.env()
  up$.on.public.web <- FALSE
  up$.get.mSet    <- function(obj = NA) obj
  up$.set.mSet    <- function(obj = NA) obj
  up$AddErrMsg    <- function(msg) stop("upstream AddErrMsg: ", msg, call. = FALSE)
  up$qs_store     <- list()
  up$ov_qs_save   <- function(obj, file, ...) up$qs_store[[file]] <- obj
  up$ov_qs_read   <- function(file, ...) up$qs_store[[file]]
  up$save.image   <- function(...) invisible(NULL)
  up$load_stringr <- function() invisible(NULL)
  sys.source(file.path(mar_dir, "R", "peaks_to_function.R"), envir = up)
  sys.source(file.path(mar_dir, "R", "util_psea.R"), envir = up)

  # .search.compoundLib() up to the empirical-compound step, verbatim.
  src <- readLines(file.path(mar_dir, "R", "peaks_to_function.R"))
  start  <- grep("^\\.search\\.compoundLib <- function\\(", src)
  anchor <- grep("# now create empirical compounds if necessary!", src, fixed = TRUE)
  anchor <- anchor[anchor > start][1]
  if (length(start) != 1L || is.na(anchor)) stop("anchors not found in upstream source")
  body <- c(src[start:(anchor - 1L)],
            "  return(list(matched_res = matched_res, ord_rt = if (exists('ord_rt')) ord_rt else NULL))",
            "}")
  eval(parse(text = sub("^\\.search\\.compoundLib <- ", ".search_matching_only <- ",
                        body)), envir = up)
  up
}

run_mode <- function(mode) {
  up <- new_upstream()
  pk <- inp$peaks
  input <- data.frame(m.z = pk$mz, p.value = pk$p_value, t.score = pk$statistic,
                      r.t = pk$rt)
  up$qs_store[["mum_raw.qs"]] <- input

  mSet <- list(
    dataSet = list(mum.type = "list", instrument = inp$ppm, mode = mode,
                   primary_ion = "yes", rt_tol = NA, rt_frac = inp$rt_frac),
    paramSet = list(mumRT = TRUE, version = "v2", mumDataContainsPval = 1,
                    peakFormat = "mprt", anal.type = "mummichog"),
    msgSet = list(read.msg = "synthetic")
  )
  # Read.PeakListData()'s own lines for these two objects:
  mSet$dataSet$mummi.orig <- cbind(input$p.value, input$m.z, input$t.score, input$r.t)
  colnames(mSet$dataSet$mummi.orig) <- c("p.value", "m.z", "t.score", "r.t")
  mSet$dataSet$pos_inx <- switch(mode,
    positive = rep(TRUE, nrow(input)),
    negative = rep(FALSE, nrow(input)),
    mixed    = pk$polarity == "positive")

  mSet <- up$SanityCheckMummichogData(mSet)
  mSet <- up$SetMummichogPval(mSet, inp$cutoff)

  invisible(capture.output(up$CreateLibFromKEGG(
    list(id = inp$compounds$id, name = inp$compounds$name, mw = inp$compounds$mw),
    list(cpds = unname(inp$pathways), id = names(inp$pathways),
         name = inp$pathway_names), "synthetic")))
  mlib <- up$qs_store[["synthetic_kegg.qs"]]
  # .setup.psea.library()'s non-MS2, non-custom-adduct branch:
  cpd.lib <- list(mz.matp = mlib$cpd.lib$adducts[["positive"]],
                  mz.matn = mlib$cpd.lib$adducts[["negative"]],
                  mw = mlib$cpd.lib$mw, id = mlib$cpd.lib$id,
                  name = mlib$cpd.lib$name)
  invisible(capture.output(m <- up$.search_matching_only(
    mSet, cpd.lib, mlib$cpd.tree[["positive"]], mlib$cpd.tree[["negative"]])))

  ds <- mSet$dataSet
  list(mode = mode, mummi.proc = ds$mummi.proc, ref_mzlist = ds$ref_mzlist,
       pos_inx = ds$pos_inx, ret_time = ds$ret_time, expr_dic = ds$expr_dic,
       rt_tol = ds$rt_tol, input_mzlist = ds$input_mzlist, N = ds$N,
       ord_rt = m$ord_rt, matched_res = m$matched_res)
}

ref <- list(
  metaboanalystr_commit = sha,
  upstream_files        = c("R/peaks_to_function.R", "R/util_psea.R"),
  upstream_functions    = c("SanityCheckMummichogData", "SetMummichogPval",
                            "CreateLibFromKEGG",
                            ".search.compoundLib (matching part)"),
  r_version             = R.version.string,
  rng_kind              = RNGkind(),
  generated_utc         = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  inputs                = inp,
  modes                 = lapply(stats::setNames(nm = c("mixed", "positive",
                                                        "negative")), run_mode)
)
out <- file.path(here, "peaks_reference.rds")
saveRDS(ref, out, version = 2)
message("wrote ", out)
for (md in names(ref$modes)) {
  r <- ref$modes[[md]]
  message(sprintf("%-8s features kept %d | significant %d | matches %d",
                  md, length(r$ref_mzlist), length(r$input_mzlist),
                  nrow(r$matched_res)))
}
