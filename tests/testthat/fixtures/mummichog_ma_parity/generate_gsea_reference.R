# tests/testthat/fixtures/mummichog_ma_parity/generate_gsea_reference.R
#
# Generates the MetaboAnalyst-compatible GSEA parity reference: the pinned
# MetaboAnalystR peaks GSEA run END TO END on the shared synthetic fixture, on
# the ECs its own library search builds — the joint mixed-polarity ECs in mixed
# mode. Run MANUALLY, never in CI.
#
#   Rscript tests/testthat/fixtures/mummichog_ma_parity/generate_gsea_reference.R \
#       /path/to/MetaboAnalystR-checkout
#
# What runs upstream, verbatim: SanityCheckMummichogData, SetMummichogPval,
# CreateLibFromKEGG, the whole .setup.psea.library (matching, ECs, ec_exp,
# pathway -> ECs), and .compute.mummichog.RT.fgsea() up to its JSON export
# (cut at "json.mat <- RJSONIO::toJSON(json.res);", as the other generators
# cut theirs) — which calls fgsea2() -> my.fgsea() -> .run_fgsea_inner().
# my.fgsea() is replaced by a direct call of .run_fgsea_inner() with the
# arguments it passes through its subprocess bridge: the bridge is session
# plumbing (qs files, a microservice), not statistics. The result is captured.
# Other stubs as in generate_stats_reference.R.
#
# Runs the four cases of the EC fixture. Output: gsea_reference.rds.

PINNED_SHA <- "398476ae2a0c996e390925fa62adc46b5ecde334"
PERM_NUM   <- 100

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: generate_gsea_reference.R <MetaboAnalystR checkout>")
mar_dir <- normalizePath(args[[1]], mustWork = TRUE)
here <- normalizePath(dirname(sub("^--file=", "", grep("^--file=",
          commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
if (is.na(here) || !nzchar(here)) here <- "."
here <- normalizePath(here, mustWork = TRUE)

sha <- system2("git", c("-C", mar_dir, "rev-parse", "HEAD"), stdout = TRUE)
if (!identical(sha, PINNED_SHA)) stop("MetaboAnalystR checkout is at ", sha)
if (!identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection"))) {
  stop("generate under R's default RNG kinds")
}
for (pkg in c("fgsea", "fastmatch", "BiocParallel", "data.table")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop(pkg, " is required")
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
  for (f in c("peaks_to_function.R", "util_psea.R", "depends_refs.R", "util_fgsea.R")) {
    sys.source(file.path(mar_dir, "R", f), envir = up)
  }
  # my.fgsea's subprocess bridge, without the subprocess: the same call its
  # child session makes, on the same arguments.
  up$my.fgsea <- function(mSetObj, pathways, stats, ranks, nperm, minSize = 1,
                          maxSize = Inf, nproc = 0, gseaParam = 1, BPPARAM = NULL) {
    res <- up$.run_fgsea_inner(mSetObj, pathways, stats, ranks, nperm, minSize,
                               maxSize, gseaParam)
    up$fgsea_res <- res
    res
  }
  src <- readLines(file.path(mar_dir, "R", "peaks_to_function.R"))
  start  <- grep("^\\.compute\\.mummichog\\.RT\\.fgsea <- function\\(", src)
  anchor <- grep("json.mat <- RJSONIO::toJSON(json.res);", src, fixed = TRUE)
  anchor <- anchor[anchor > start][1]
  if (length(start) != 1L || is.na(anchor)) stop("anchors not found in upstream source")
  eval(parse(text = sub("^\\.compute\\.mummichog\\.RT\\.fgsea <- ",
                        ".compute_rt_fgsea_no_json <- ",
                        c(src[start:(anchor - 1L)], "  return(mSetObj)", "}"))),
       envir = up)
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
                    peakFormat = "mprt", anal.type = "gsea_peaks"),
    msgSet = list(read.msg = "synthetic"),
    mum_nm_csv = "mummichog_pathway_enrichment_gsea.csv")
  mSet$dataSet$mummi.orig <- cbind(input$p.value, input$m.z, input$t.score, input$r.t)
  colnames(mSet$dataSet$mummi.orig) <- c("p.value", "m.z", "t.score", "r.t")
  mSet$dataSet$pos_inx <- switch(mode,
    positive = rep(TRUE, nrow(input)), negative = rep(FALSE, nrow(input)),
    mixed = pk$polarity == "positive")

  mSet <- up$SanityCheckMummichogData(mSet)
  mSet <- up$SetMummichogPval(mSet, inp$cutoff)
  invisible(capture.output(up$CreateLibFromKEGG(
    list(id = inp$compounds$id, name = inp$compounds$name, mw = inp$compounds$mw),
    list(cpds = unname(inp$pathways), id = names(inp$pathways),
         name = inp$pathway_names), LIB)))
  up$qs_store[[paste0(LIB, ".qs")]] <- up$qs_store[[paste0(LIB, "_kegg.qs")]]
  wd <- tempfile("ma_gsea_"); dir.create(wd)
  file.create(file.path(wd, paste0(LIB, ".qs")))
  withr::with_dir(wd, invisible(capture.output(suppressWarnings({
    mSet <- up$.setup.psea.library(mSet, LIB, "current", minLib = 3)
    mSet$dataSet$paramSet <- mSet$paramSet
    mSet <- up$.compute_rt_fgsea_no_json(mSet, PERM_NUM)
  }))))
  if (!is.list(mSet)) stop("upstream GSEA failed for ", mode, "/", primary_ion)
  fr <- up$fgsea_res
  list(mode = mode, primary_ion = primary_ion,
       ec_exp = mSet$ec_exp, emp_cpds = mSet$pathways$emp_cpds,
       pathway_names = mSet$pathways$name,
       fgsea = data.frame(
         pathway = as.character(fr$pathway), pval = as.numeric(fr$pval),
         padj = as.numeric(fr$padj), ES = as.numeric(fr$ES), NES = as.numeric(fr$NES),
         nMoreExtreme = as.numeric(fr$nMoreExtreme), size = as.numeric(fr$size),
         leadingEdgeMatched = as.character(fr$leadingEdgeMatched),
         stringsAsFactors = FALSE),
       resmat = mSet$mummi.gsea.resmat,
       pathways_filtered = up$qs_store[["pathwaysFiltered.qs"]],
       csv = up$csv_store[["mummichog_pathway_enrichment_gsea.csv"]])
}

cases <- list(
  mixed            = run_case("mixed", "yes"),
  positive         = run_case("positive", "yes"),
  negative         = run_case("negative", "yes"),
  mixed_no_primary = run_case("mixed", "no"))

ref <- list(
  metaboanalystr_commit = sha,
  upstream_entry = ".compute.mummichog.RT.fgsea -> fgsea2 -> my.fgsea -> .run_fgsea_inner",
  r_version = R.version.string, rng_kind = RNGkind(),
  locale_collate = Sys.getlocale("LC_COLLATE"),
  fgsea_version = as.character(utils::packageVersion("fgsea")),
  generated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  perm_num = PERM_NUM, min_lib = 3, inputs = inp, cases = cases)
out <- file.path(here, "gsea_reference.rds")
saveRDS(ref, out, version = 2)
message("wrote ", out, " (fgsea ", ref$fgsea_version, ", LC_COLLATE ", ref$locale_collate, ")")
for (nm in names(cases)) {
  message(sprintf("%-17s tested pathways %d | ECs scored %d", nm,
                  nrow(cases[[nm]]$fgsea), length(cases[[nm]]$ec_exp)))
}
