# tests/testthat/fixtures/mummichog_ma_parity/generate_ecs_reference.R
#
# Generates the EMPIRICAL-COMPOUND parity reference (component 3: EC
# construction, the dictionaries, and pathway -> EC sets) by running the PINNED
# MetaboAnalystR code itself. Run MANUALLY, never in CI.
#
#   Rscript tests/testthat/fixtures/mummichog_ma_parity/generate_ecs_reference.R \
#       /path/to/MetaboAnalystR-checkout
#
# The checkout must be at commit 398476ae2a0c996e390925fa62adc46b5ecde334.
#
# What runs upstream, verbatim (R/peaks_to_function.R, R/util_psea.R and
# R/depends_refs.R sourced into one environment): SanityCheckMummichogData(),
# SetMummichogPval(), CreateLibFromKEGG(), and then the WHOLE
# .setup.psea.library() — library read, minLib filter, the complete
# .search.compoundLib() (matching, ECs, dictionaries) and the pathway -> EC
# mapping — exactly as PerformPSEA() calls it.
#
# Stubbed (session I/O only, no maths): .get.mSet/.set.mSet (identity),
# AddErrMsg (stop), AddMsg (no-op), ov_qs_save/ov_qs_read (an in-memory store;
# the library is served from it under "<lib>.qs"), fast.write.csv, save.image,
# load_stringr; `%||%` is defined only when this R lacks base::`%||%`.
# .setup.psea.library() checks file.exists("<lib>.qs"), so it runs in a
# temporary working directory holding an empty file of that name.
#
# EC labels follow the collation of the session locale (see 06l); the locale
# used is recorded and the parity test compares under the same locale only.
#
# Runs: mixed, positive and negative mode with primary_ion = "yes"; mixed with
# primary_ion = "no". Output: ecs_reference.rds next to this script.

PINNED_SHA <- "398476ae2a0c996e390925fa62adc46b5ecde334"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: generate_ecs_reference.R <MetaboAnalystR checkout>")
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
  up$ov_qs_save     <- function(obj, file, ...) up$qs_store[[file]] <- obj
  up$ov_qs_read     <- function(file, ...) up$qs_store[[file]]
  up$fast.write.csv <- function(...) invisible(NULL)
  up$save.image     <- function(...) invisible(NULL)
  up$load_stringr   <- function() invisible(NULL)
  if (!exists("%||%", envir = baseenv(), inherits = FALSE)) {
    up$`%||%` <- function(x, y) if (is.null(x)) y else x
  }
  for (f in c("peaks_to_function.R", "util_psea.R", "depends_refs.R")) {
    sys.source(file.path(mar_dir, "R", f), envir = up)
  }
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
         name = inp$pathway_names), LIB)))
  up$qs_store[[paste0(LIB, ".qs")]] <- up$qs_store[[paste0(LIB, "_kegg.qs")]]

  wd <- tempfile("ma_ecs_")
  dir.create(wd)
  file.create(file.path(wd, paste0(LIB, ".qs")))
  withr::with_dir(wd, invisible(capture.output(
    mSet <- up$.setup.psea.library(mSet, LIB, "current", minLib = 3))))
  if (!is.list(mSet)) stop("upstream .setup.psea.library() failed")

  keep <- c("ec_exp_dict", "cpd_exp_dict", "ec_exp", "mz2cpd_dict",
            "mz2ec_dict", "ec2mz_dict", "ecpd_cpd_dict", "cpd_ecpd_dict",
            "cpd_ecpd_counts", "input_ecpdlist", "total_matched_ecpds",
            "cpd_form_dict", "pathways")
  c(list(mode = mode, primary_ion = primary_ion,
         initial_ecs = up$qs_store[["initial_ecs.qs"]],
         matched_res = up$qs_store[["mum_res.qs"]]),
    mSet[keep])
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
  upstream_functions    = c("SanityCheckMummichogData", "SetMummichogPval",
                            "CreateLibFromKEGG", ".setup.psea.library",
                            ".search.compoundLib", "Convert2Dictionary",
                            "melt.list"),
  r_version             = R.version.string,
  rng_kind              = RNGkind(),
  locale_collate        = Sys.getlocale("LC_COLLATE"),
  igraph_version        = as.character(utils::packageVersion("igraph")),
  generated_utc         = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  min_lib               = 3,
  inputs                = inp,
  cases                 = cases
)
out <- file.path(here, "ecs_reference.rds")
saveRDS(ref, out, version = 2)
message("wrote ", out, " (LC_COLLATE = ", ref$locale_collate, ")")
for (nm in names(cases)) {
  r <- cases[[nm]]
  message(sprintf("%-17s ECs %d | significant ECs %d | pathways %d | rows %d",
                  nm, length(r$total_matched_ecpds), length(r$input_ecpdlist),
                  length(r$pathways$name), nrow(r$matched_res)))
}
