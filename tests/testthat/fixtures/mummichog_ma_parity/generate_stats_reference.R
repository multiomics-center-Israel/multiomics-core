# tests/testthat/fixtures/mummichog_ma_parity/generate_stats_reference.R
#
# Generates the PATHWAY-STATISTICS parity reference (component 4: permutations,
# FET / EASE, Gamma adjustment, empirical p-values, ordering, rounding and the
# Holm-adjusted result table) by running the PINNED MetaboAnalystR code itself.
# Run MANUALLY, never in CI.
#
#   Rscript tests/testthat/fixtures/mummichog_ma_parity/generate_stats_reference.R \
#       /path/to/MetaboAnalystR-checkout
#
# The checkout must be at commit 398476ae2a0c996e390925fa62adc46b5ecde334.
# Needs fitdistrplus installed: without it upstream's own tryCatch would quietly
# make every Gamma NA, and the reference would not test the Gamma path.
#
# What runs upstream, verbatim: everything component 3's generator runs
# (SanityCheckMummichogData, SetMummichogPval, CreateLibFromKEGG,
# .setup.psea.library), then — as .init.RT.Permutations() does for
# anal.type "mummichog" —
#   * .perform.mummichogRTPermutations(mSet, permNum = 100), whole;
#   * .compute.mummichogRTSigPvals(mSet): its source lines from the function
#     start up to (not including) the line
#     "json.mat <- RJSONIO::toJSON(json.res);" — i.e. every statistic, the
#     ordering, rounding and .save.mummichog.restable() — closed with
#     return(mSetObj). The cut drops only the JSON written for upstream's web
#     network view (which would need RJSONIO) and the web-table bookkeeping.
# Stubbed (session I/O only, no maths): as component 3, plus fast.write.csv,
# which stores the table it is given (the result CSV's content).
#
# Runs the four component 3 cases. Output: stats_reference.rds.

PINNED_SHA <- "398476ae2a0c996e390925fa62adc46b5ecde334"
PERM_NUM   <- 100

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: generate_stats_reference.R <MetaboAnalystR checkout>")
mar_dir <- normalizePath(args[[1]], mustWork = TRUE)
here <- normalizePath(dirname(sub("^--file=", "", grep("^--file=",
          commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
if (is.na(here) || !nzchar(here)) here <- "."
here <- normalizePath(here, mustWork = TRUE)

sha <- system2("git", c("-C", mar_dir, "rev-parse", "HEAD"), stdout = TRUE)
if (!identical(sha, PINNED_SHA)) {
  stop("MetaboAnalystR checkout is at ", sha, ", not the pinned ", PINNED_SHA)
}
if (!identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection"))) {
  stop("generate under R's default RNG kinds; got ", paste(RNGkind(), collapse = "/"))
}
if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
  stop("fitdistrplus is not installed: upstream's Gamma fit would fail ",
       "silently and the reference would carry no Gamma values.")
}

source(file.path(here, "fixture_inputs.R"))
inp <- ma_fixture_inputs()
LIB <- "synthetic"

new_upstream <- function() {
  up <- new.env()
  up$.on.public.web <- FALSE
  up$.get.mSet      <- function(obj = NA) obj
  up$.set.mSet      <- function(obj = NA) obj
  up$AddErrMsg      <- function(msg) stop("upstream AddErrMsg: ", msg, call. = FALSE)
  up$AddMsg         <- function(msg) invisible(NULL)
  up$qs_store       <- list()
  up$csv_store      <- list()
  up$ov_qs_save     <- function(obj, file, ...) up$qs_store[[file]] <- obj
  up$ov_qs_read     <- function(file, ...) up$qs_store[[file]]
  up$fast.write.csv <- function(dat, file, ...) up$csv_store[[file]] <- dat
  up$save.image     <- function(...) invisible(NULL)
  up$load_stringr   <- function() invisible(NULL)
  if (!exists("%||%", envir = baseenv(), inherits = FALSE)) {
    up$`%||%` <- function(x, y) if (is.null(x)) y else x
  }
  for (f in c("peaks_to_function.R", "util_psea.R", "depends_refs.R")) {
    sys.source(file.path(mar_dir, "R", f), envir = up)
  }

  # .compute.mummichogRTSigPvals() up to its JSON export, verbatim.
  src <- readLines(file.path(mar_dir, "R", "peaks_to_function.R"))
  start  <- grep("^\\.compute\\.mummichogRTSigPvals <- function\\(", src)
  anchor <- grep("json.mat <- RJSONIO::toJSON(json.res);", src, fixed = TRUE)
  anchor <- anchor[anchor > start][1]
  if (length(start) != 1L || is.na(anchor)) stop("anchors not found in upstream source")
  body <- c(src[start:(anchor - 1L)], "  return(mSetObj)", "}")
  eval(parse(text = sub("^\\.compute\\.mummichogRTSigPvals <- ",
                        ".compute_sig_pvals_no_json <- ", body)), envir = up)
  up
}

run_case <- function(mode, primary_ion) {
  up <- new_upstream()
  pk <- inp$peaks
  input <- data.frame(m.z = pk$mz, p.value = pk$p_value, t.score = pk$statistic,
                      r.t = pk$rt)
  up$qs_store[["mum_raw.qs"]] <- input
  mSet <- list(
    dataSet = list(mum.type = "list", instrument = inp$ppm, mode = mode,
                   primary_ion = primary_ion, rt_tol = NA, rt_frac = inp$rt_frac),
    paramSet = list(mumRT = TRUE, version = "v2", mumDataContainsPval = 1,
                    peakFormat = "mprt", anal.type = "mummichog"),
    msgSet = list(read.msg = "synthetic")
  )
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
         name = inp$pathway_names), LIB)))
  up$qs_store[[paste0(LIB, ".qs")]] <- up$qs_store[[paste0(LIB, "_kegg.qs")]]

  wd <- tempfile("ma_stats_")
  dir.create(wd)
  file.create(file.path(wd, paste0(LIB, ".qs")))
  withr::with_dir(wd, invisible(capture.output({
    mSet <- up$.setup.psea.library(mSet, LIB, "current", minLib = 3)
    mSet$dataSet$paramSet <- mSet$paramSet               # as PerformPSEA()
    mSet <- up$.perform.mummichogRTPermutations(mSet, PERM_NUM)
    mSet <- up$.compute_sig_pvals_no_json(mSet)
  })))
  if (!is.list(mSet)) stop("upstream statistics failed for ", mode, "/", primary_ion)

  list(mode = mode, primary_ion = primary_ion,
       perm_record = mSet$perm_record, perm_hits = mSet$perm_hits,
       pvals = mSet$pvals, mummi.resmat = mSet$mummi.resmat,
       path.nms = mSet$path.nms, path.pval = mSet$path.pval,
       restable = up$csv_store[[mSet$mum_nm_csv]])
}

cases <- list(
  mixed            = run_case("mixed", "yes"),
  positive         = run_case("positive", "yes"),
  negative         = run_case("negative", "yes"),
  mixed_no_primary = run_case("mixed", "no")
)

ref <- list(
  metaboanalystr_commit = sha,
  upstream_files        = c("R/peaks_to_function.R", "R/util_psea.R",
                            "R/depends_refs.R"),
  upstream_functions    = c(".perform.mummichogRTPermutations",
                            "ComputeMummichogRTPermPvals", "make_ecpdlist",
                            ".compute.mummichogRTSigPvals (up to its JSON export)",
                            ".save.mummichog.restable"),
  r_version             = R.version.string,
  rng_kind              = RNGkind(),
  locale_collate        = Sys.getlocale("LC_COLLATE"),
  fitdistrplus_version  = as.character(utils::packageVersion("fitdistrplus")),
  igraph_version        = as.character(utils::packageVersion("igraph")),
  generated_utc         = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  perm_num              = PERM_NUM,
  min_lib               = 3,
  inputs                = inp,
  cases                 = cases
)
out <- file.path(here, "stats_reference.rds")
saveRDS(ref, out, version = 2)
message("wrote ", out, " (LC_COLLATE = ", ref$locale_collate,
        ", fitdistrplus ", ref$fitdistrplus_version, ")")
for (nm in names(cases)) {
  r <- cases[[nm]]
  g <- r$mummi.resmat[, "Gamma"]
  message(sprintf("%-17s pathways reported %d | Gamma %s", nm,
                  nrow(r$mummi.resmat),
                  if (all(is.na(g))) "ALL NA (fit failed or skipped)" else "fitted"))
}
