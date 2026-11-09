# tests/testthat/fixtures/mummichog_ma_parity/validate_real_parity.R
#
# Real-data parity: runs ONE contrast of a MetaboAnalyst-compatible run through
# the PINNED MetaboAnalystR itself and through our engine, from the identical
# input, and compares every intermediate. Run MANUALLY; writes only to --out.
#
#   Rscript tests/testthat/fixtures/mummichog_ma_parity/validate_real_parity.R \
#       --mar /path/to/MetaboAnalystR-checkout \
#       --ma-dir <run>/metabolomics/mummichog_pinned/<contrast>/metaboanalyst \
#       --model /path/to/model.json \
#       --out /path/to/parity_out
#
# Inputs are taken from the run itself, so nothing can drift between the two
# implementations: ma_input.tsv (the exact peak list the engine was given, in
# its order) and ma_run_info.json (mode, ppm, cutoff, permutations, primary-ion
# rule, rt_frac, min pathway size, model sha256). --model must be that model
# (its sha256 is checked).
#
# Upstream runs verbatim, as in the fixture generators: SanityCheckMummichogData,
# SetMummichogPval, CreateLibFromKEGG (from the model's compounds and the
# pathway membership our converter derives — validated identical to upstream's
# builder in component 1), the whole .setup.psea.library (matching, ECs,
# dictionaries, pathway -> ECs), .perform.mummichogRTPermutations, and
# .compute.mummichogRTSigPvals up to its JSON export. Only session I/O is
# stubbed. Neither implementation is modified.
#
# Output: PASS/FAIL per compared object, with the number of differing elements
# (no data values) on the console; full objects and a difference report in
# --out for local inspection.

PINNED_SHA <- "398476ae2a0c996e390925fa62adc46b5ecde334"

args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing ", flag)
  args[[i + 1L]]
}
mar_dir <- normalizePath(opt("--mar"), mustWork = TRUE)
ma_dir  <- normalizePath(opt("--ma-dir"), mustWork = TRUE)
model   <- normalizePath(opt("--model"), mustWork = TRUE)
out     <- opt("--out")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!file.exists("_targets.R") || !dir.exists("R")) stop("run from the multiomics-core root")

for (layer in c("core", "services", "domain", "modules")) {
  for (f in sort(list.files(file.path("R", layer), pattern = "\\.R$",
                            full.names = TRUE, recursive = TRUE))) source(f)
}

results <- list()
check <- function(name, ok, detail = "") {
  results[[name]] <<- isTRUE(ok)
  cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", name,
              if (nzchar(detail)) paste0(" -- ", detail) else ""))
}

# ---- environment and inputs ---------------------------------------------------
sha <- system2("git", c("-C", mar_dir, "rev-parse", "HEAD"), stdout = TRUE)
check("MetaboAnalystR at the pinned commit", identical(sha, PINNED_SHA), sha)
check("default RNG kinds", identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection")))
check("fitdistrplus 1.2.6", requireNamespace("fitdistrplus", quietly = TRUE) &&
        identical(as.character(utils::packageVersion("fitdistrplus")), "1.2.6"))
info  <- jsonlite::fromJSON(file.path(ma_dir, "ma_run_info.json"), simplifyVector = TRUE)
p     <- info$params
peaks <- as.data.frame(readr::read_tsv(file.path(ma_dir, "ma_input.tsv"),
                                       show_col_types = FALSE,
                                       col_types = readr::cols(
                                         feature_id = "c", polarity = "c",
                                         .default = "d")))
check("model is the one the run used", identical(.mmc_sha256_file(model), info$model$sha256))
if (!all(unlist(results))) stop("environment/input checks failed; nothing compared")
fpi <- isTRUE(p$force_primary_ion)

# ---- ours -----------------------------------------------------------------------
lib <- mmc_ma_library_from_model(model)
gp <- info$gsea$params
gsea_settings <- if (isTRUE(info$gsea$ran)) list(
  n_perm = gp$n_perm, min_size = gp$min_size,
  max_size = suppressWarnings(as.numeric(gp$max_size)), seed = gp$seed,
  metric_label = info$gsea$metric_label) else NULL
ours <- suppressWarnings(mmc_ma_run_contrast(
  peaks, lib, p$ionization_mode, ppm = p$tolerance_ppm, cutoff = p$p_cutoff,
  n_perm = p$n_permutations, force_primary_ion = fpi, rt_frac = p$rt_frac,
  min_lib = p$min_pathway_size, gsea = gsea_settings))

# ---- upstream -------------------------------------------------------------------
up <- new.env()
up$.on.public.web <- FALSE
up$.get.mSet <- function(obj = NA) obj
up$.set.mSet <- function(obj = NA) obj
up$AddErrMsg <- function(msg) stop("upstream AddErrMsg: ", msg, call. = FALSE)
up$AddMsg <- function(msg) invisible(NULL)
up$qs_store <- list(); up$csv_store <- list()
up$ov_qs_save <- function(obj, file, ...) up$qs_store[[file]] <- obj
up$ov_qs_read <- function(file, ...) up$qs_store[[file]]
up$fast.write.csv <- function(dat, file, ...) up$csv_store[[file]] <- dat
up$save.image <- function(...) invisible(NULL)
up$load_stringr <- function() invisible(NULL)
if (!exists("%||%", envir = baseenv(), inherits = FALSE)) {
  up$`%||%` <- function(x, y) if (is.null(x)) y else x
}
for (f in c("peaks_to_function.R", "util_psea.R", "depends_refs.R", "util_fgsea.R")) {
  sys.source(file.path(mar_dir, "R", f), envir = up)
}
# my.fgsea's subprocess bridge without the subprocess: the same call, same args.
up$my.fgsea <- function(mSetObj, pathways, stats, ranks, nperm, minSize = 1,
                        maxSize = Inf, nproc = 0, gseaParam = 1, BPPARAM = NULL) {
  res <- up$.run_fgsea_inner(mSetObj, pathways, stats, ranks, nperm, minSize,
                             maxSize, gseaParam)
  up$fgsea_res <- res
  res
}
gstart <- grep("^\\.compute\\.mummichog\\.RT\\.fgsea <- function\\(", src_all <-
                 readLines(file.path(mar_dir, "R", "peaks_to_function.R")))
ganchor <- grep("json.mat <- RJSONIO::toJSON(json.res);", src_all, fixed = TRUE)
ganchor <- ganchor[ganchor > gstart][1]
eval(parse(text = sub("^\\.compute\\.mummichog\\.RT\\.fgsea <- ",
                      ".compute_rt_fgsea_no_json <- ",
                      c(src_all[gstart:(ganchor - 1L)], "  return(mSetObj)", "}"))),
     envir = up)
src <- readLines(file.path(mar_dir, "R", "peaks_to_function.R"))
start  <- grep("^\\.compute\\.mummichogRTSigPvals <- function\\(", src)
anchor <- grep("json.mat <- RJSONIO::toJSON(json.res);", src, fixed = TRUE)
anchor <- anchor[anchor > start][1]
eval(parse(text = sub("^\\.compute\\.mummichogRTSigPvals <- ",
                      ".compute_sig_pvals_no_json <- ",
                      c(src[start:(anchor - 1L)], "  return(mSetObj)", "}"))), envir = up)

input <- data.frame(m.z = peaks$mz, p.value = peaks$p_value,
                    t.score = peaks$statistic, r.t = peaks$rt)
up$qs_store[["mum_raw.qs"]] <- input
mSet <- list(
  dataSet = list(mum.type = "list", instrument = p$tolerance_ppm,
                 mode = p$ionization_mode, primary_ion = if (fpi) "yes" else "no",
                 rt_tol = NA, rt_frac = p$rt_frac),
  paramSet = list(mumRT = TRUE, version = "v2", mumDataContainsPval = 1,
                  peakFormat = "mprt", anal.type = "mummichog"),
  msgSet = list(read.msg = "real-data parity"))
mSet$dataSet$mummi.orig <- cbind(input$p.value, input$m.z, input$t.score, input$r.t)
colnames(mSet$dataSet$mummi.orig) <- c("p.value", "m.z", "t.score", "r.t")
mSet$dataSet$pos_inx <- switch(p$ionization_mode,
  positive = rep(TRUE, nrow(input)), negative = rep(FALSE, nrow(input)),
  mixed = peaks$polarity == "positive")

t0 <- Sys.time()
mSet <- up$SanityCheckMummichogData(mSet)
mSet <- up$SetMummichogPval(mSet, p$p_cutoff)
invisible(capture.output(up$CreateLibFromKEGG(
  list(id = lib$cpd.lib$id, name = lib$cpd.lib$name, mw = lib$cpd.lib$mw),
  lib$pathways[c("cpds", "id", "name")], "realdata")))
up$qs_store[["realdata.qs"]] <- up$qs_store[["realdata_kegg.qs"]]
wd <- tempfile("ma_parity_"); dir.create(wd)
file.create(file.path(wd, "realdata.qs"))
up_warn <- character()
withr::with_dir(wd, invisible(capture.output(withCallingHandlers({
  mSet <- up$.setup.psea.library(mSet, "realdata", "current",
                                 minLib = p$min_pathway_size)
  mSet$dataSet$paramSet <- mSet$paramSet
  mSet <- up$.perform.mummichogRTPermutations(mSet, p$n_permutations)
  mSet <- up$.compute_sig_pvals_no_json(mSet)
  if (!is.null(gsea_settings)) {
    gset <- mSet
    gset$mum_nm_csv <- "gsea.csv"
    gset <- up$.compute_rt_fgsea_no_json(gset, gsea_settings$n_perm)
    up$gsea_resmat <- gset$mummi.gsea.resmat
  }
}, warning = function(w) { up_warn <<- c(up_warn, conditionMessage(w))
                           invokeRestart("muffleWarning") }))))
cat(sprintf("upstream ran in %.1f min (%d warning(s))\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), length(up_warn)))

# ---- compare --------------------------------------------------------------------
n_diff <- function(a, b) {
  if (identical(a, b)) return(0L)
  if (length(a) != length(b)) return(NA_integer_)
  sum(mapply(function(x, y) !identical(x, y), as.list(a), as.list(b)))
}
cmp <- function(name, a, b) {
  same <- identical(a, b)
  d <- if (same) "" else {
    k <- n_diff(a, b)
    if (is.na(k)) sprintf("lengths %d vs %d", length(a), length(b)) else
      sprintf("%d element(s) differ", k)
  }
  check(name, same, d)
}
ds <- mSet$dataSet
cmp("peak list (m/z keys, order)",          ours$prep$ref_mzlist, ds$ref_mzlist)
cmp("polarity per feature",                 ours$prep$pos_inx, ds$pos_inx)
cmp("retention times",                      ours$prep$ret_time, ds$ret_time)
cmp("RT tolerance",                         ours$prep$rt_tol, ds$rt_tol)
cmp("significant m/z",                      ours$input_mzlist, ds$input_mzlist)
cmp("initial empirical compounds",          ours$ecs$initial_ecs, up$qs_store[["initial_ecs.qs"]])
cmp("matched features + final ECs",         ours$ecs$matched_res, up$qs_store[["mum_res.qs"]])
for (d in c("cpd_ecpd_dict", "ecpd_cpd_dict", "mz2ec_dict", "ec_exp",
            "cpd_ecpd_counts", "input_ecpdlist", "total_matched_ecpds")) {
  cmp(paste("EC dictionary", d), ours$ecs[[d]], mSet[[d]])
}
cmp("pathway compound sets",                ours$pathways$cpds, mSet$pathways$cpds)
cmp("pathway names",                        ours$pathways$name, mSet$pathways$name)
cmp("pathway EC sets",                      ours$pathways$emp_cpds, mSet$pathways$emp_cpds)
cmp("permutation records",                  ours$perms$perm_record, mSet$perm_record)
cmp("permutation hits",                     ours$perms$perm_hits, mSet$perm_hits)
cmp("raw pathway statistics (FET, EASE, counts)", ours$stats$pvals, mSet$pvals)
cmp("reported table (Gamma, empirical, order)", ours$stats$resmat[, -11, drop = FALSE],
    mSet$mummi.resmat)
cmp("result CSV incl. Holm columns",        ours$stats$restable,
    up$csv_store[[mSet$mum_nm_csv]])

# ---- GSEA ---------------------------------------------------------------------
if (!is.null(gsea_settings)) {
  check("GSEA ran with MetaboAnalyst's defaults (permutations, size, seed)",
        identical(as.numeric(c(gsea_settings$n_perm, gsea_settings$min_size, gsea_settings$seed)),
                  c(100, 1, 123)) && is.infinite(gsea_settings$max_size))
  cmp("GSEA EC scores (ec_exp)", ours$ecs$ec_exp, mSet$ec_exp)
  f <- ours$gsea$fgsea; e <- up$fgsea_res
  cmp("GSEA tested pathways and order", f$pathway, as.character(e$pathway))
  cmp("GSEA tested sizes", as.numeric(f$size), as.numeric(e$size))
  cmp("GSEA leading edges", ifelse(is.na(f$leadingEdge), "", f$leadingEdge),
      as.character(e$leadingEdgeMatched))
  for (col in c("ES", "NES", "pval", "padj", "nMoreExtreme")) {
    d <- max(abs(as.numeric(f[[col]]) - as.numeric(e[[col]])), na.rm = TRUE)
    check(paste("GSEA", col, "(|diff| <= 1e-8)"), is.finite(d) && d <= 1e-8,
          sprintf("max |diff| %.3g", d))
  }
  cmp("GSEA result table (rounded, ordered)", ours$gsea$resmat, up$gsea_resmat)
  own <- vapply(rownames(ours$gsea$resmat), function(pw)
    paste(names(up$qs_store[["pathwaysFiltered.qs"]][[pw]]), collapse = ";"),
    character(1), USE.NAMES = FALSE)
  cmp("GSEA EC hits per pathway", ours$gsea$table$EmpiricalCompounds, own)
  gnom <- function(m) rownames(m)[m[, "P_val"] < p$p_cutoff]
  gbh  <- function(m) rownames(m)[m[, "P_adj"] < p$p_cutoff]
  check("same pathways nominal in GSEA (p < cutoff)", setequal(gnom(ours$gsea$resmat), gnom(up$gsea_resmat)))
  check("same pathways passing GSEA BH FDR < cutoff", setequal(gbh(ours$gsea$resmat), gbh(up$gsea_resmat)))
  gw <- file.path(ma_dir, "ma_gsea_results.tsv")
  if (file.exists(gw)) {
    wr <- readr::read_tsv(gw, show_col_types = FALSE)
    check("the run's written GSEA table equals this re-run",
          identical(wr$Pathway, ours$gsea$table$Pathway) &&
            isTRUE(all.equal(wr$P.Value, ours$gsea$table$P.Value)) &&
            isTRUE(all.equal(wr$padj, ours$gsea$table$padj)) &&
            isTRUE(all.equal(wr$NES, ours$gsea$table$NES)))
  }
} else {
  check("the run recorded a GSEA", FALSE, "ma_run_info.json says GSEA did not run")
}

# the run's own written CSV is this engine's output
written <- utils::read.csv(file.path(ma_dir, "mummichog_pathway_enrichment_mummichog.csv"),
                           check.names = FALSE, row.names = 1, colClasses = "character")
check("the run's written CSV equals this re-run",
      identical(unname(as.matrix(written)), unname(ours$stats$restable)))

# ---- materiality, if anything differs ---------------------------------------------
cut <- p$p_cutoff
a <- ours$stats$resmat; b <- mSet$mummi.resmat
nom <- function(m) rownames(m)[m[, "FET"] < cut]
holm <- function(m) { x <- stats::p.adjust(m[, "FET"]); rownames(m)[x < cut] }
check("same pathways nominal at FET < cutoff", setequal(nom(a), nom(b)),
      sprintf("ours only %d, upstream only %d", length(setdiff(nom(a), nom(b))),
              length(setdiff(nom(b), nom(a)))))
check("same pathways passing Holm < cutoff", setequal(holm(a), holm(b)))

saveRDS(list(ours = ours, upstream_gsea = list(fgsea = up$fgsea_res, resmat = up$gsea_resmat),
             upstream = mSet[c("pvals", "mummi.resmat", "perm_record",
                                            "pathways", "input_ecpdlist",
                                            "total_matched_ecpds")],
             upstream_initial_ecs = up$qs_store[["initial_ecs.qs"]],
             upstream_mum_res = up$qs_store[["mum_res.qs"]],
             upstream_restable = up$csv_store[[mSet$mum_nm_csv]],
             upstream_warnings = up_warn, run_info = info,
             compared = unlist(results)),
        file.path(out, "real_parity_objects.rds"))
writeLines(c(sprintf("MetaboAnalystR %s | R %s | fitdistrplus %s | LC_COLLATE %s",
                     sha, R.version.string,
                     as.character(utils::packageVersion("fitdistrplus")),
                     Sys.getlocale("LC_COLLATE")),
             sprintf("%s: %s", names(results), ifelse(unlist(results), "PASS", "FAIL"))),
           file.path(out, "real_parity_summary.txt"))
cat("\nwrote ", normalizePath(out), "\n", sep = "")
if (!all(unlist(results))) quit(status = 1)
