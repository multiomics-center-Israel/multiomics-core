# tests/testthat/fixtures/mummichog_ma_parity/generate_stats_targeted_reference.R
#
# Generates the TARGETED pathway-statistics reference: a hand-built EC universe
# whose pathway p-values are informative, so the reference can tell apart
#   * Gamma with the scale estimate passed as `rate` (upstream) from `scale`;
#   * Holm (upstream's p.adjust default) from BH.
# The end-to-end fixture (generate_stats_reference.R) cannot: its significant
# set is so large relative to its EC universe that every reported FET is 1.
# Run MANUALLY, never in CI.
#
#   Rscript tests/testthat/fixtures/mummichog_ma_parity/generate_stats_targeted_reference.R \
#       /path/to/MetaboAnalystR-checkout
#
# What runs upstream, verbatim: .perform.mummichogRTPermutations(permNum = 100)
# whole, and .compute.mummichogRTSigPvals() up to its JSON export (exactly as
# generate_stats_reference.R cuts it). Their inputs — the EC universe, pathway
# EC sets, the m/z -> EC dictionary and the significant ECs — are given here
# directly instead of being built from peaks: components 2-3 are covered by
# their own references, and this fixture isolates the statistics.
#
# Design (all synthetic): 1000 m/z keys, the first 200 mapping one-to-one to
# ECs E001-E200 (compound C90001-C90200); 6 disjoint pathways of 20, 20, 25,
# 20, 15, 10 ECs holding 10, 8, 8, 7, 2, 0 of the 35 significant ECs; 40
# significant m/z (the 35 + 5 unmatched). Expected FET roughly 4e-4, 0.011,
# 0.045, 0.038, 0.78 (and no hit: dropped).
#
# Also stored, as documented DIAGNOSTICS (not upstream output): the same
# fitdistrplus::fitdist() call on the same permutation null, with its
# estimates, standard errors, optimiser convergence code and the warnings it
# raised, to show which of them the Gamma p-values use.
#
# Output: stats_targeted_reference.rds next to this script.

PINNED_SHA <- "398476ae2a0c996e390925fa62adc46b5ecde334"
PERM_NUM   <- 100

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: generate_stats_targeted_reference.R <MetaboAnalystR checkout>")
}
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
  stop("fitdistrplus is not installed: upstream's Gamma fit would fail silently.")
}

# ---- the targeted inputs ----------------------------------------------------
ec  <- sprintf("E%03d", 1:200)
cpd <- sprintf("C9%04d", 1:200)
ref_mzlist <- as.character(round(100 + (1:1000) * 0.731, 4))
sizes <- c(20, 20, 25, 20, 15, 10)
n_sig <- c(10, 8, 8, 7, 2, 0)
ends <- cumsum(sizes)
starts <- ends - sizes + 1
emp_cpds <- lapply(seq_along(sizes), function(k) ec[starts[k]:ends[k]])
sig_ec <- unlist(lapply(seq_along(sizes), function(k)
  ec[starts[k] + seq_len(n_sig[k]) - 1]))
inputs <- list(
  ref_mzlist = ref_mzlist,
  input_mzlist = c(ref_mzlist[match(sig_ec, ec)], ref_mzlist[201:205]),
  pathways = list(
    cpds = lapply(emp_cpds, function(e) cpd[match(e, ec)]),
    name = sprintf("Targeted pathway %d", seq_along(sizes)),
    emp_cpds = emp_cpds),
  input_ecpdlist = sig_ec,
  total_matched_ecpds = ec,
  mz2ec_dict = stats::setNames(as.list(ec), ref_mzlist[1:200]),
  ecpd_cpd_dict = stats::setNames(as.list(cpd), ec),
  cpd_ecpd_dict = stats::setNames(as.list(ec), cpd),
  ec_exp = stats::setNames(rep(1, 200), ec),
  cpd_form_dict = stats::setNames(as.list(rep("M+H[1+]", 200)), cpd)
)

# ---- upstream -----------------------------------------------------------------
up <- new.env()
up$.on.public.web <- FALSE
up$AddErrMsg      <- function(msg) stop("upstream AddErrMsg: ", msg, call. = FALSE)
up$AddMsg         <- function(msg) invisible(NULL)
up$qs_store       <- list(mum_res.qs = data.frame(Query.Mass = character()))
up$csv_store      <- list()
up$ov_qs_save     <- function(obj, file, ...) up$qs_store[[file]] <- obj
up$ov_qs_read     <- function(file, ...) up$qs_store[[file]]
up$fast.write.csv <- function(dat, file, ...) up$csv_store[[file]] <- dat
sys.source(file.path(mar_dir, "R", "peaks_to_function.R"), envir = up)
src <- readLines(file.path(mar_dir, "R", "peaks_to_function.R"))
start  <- grep("^\\.compute\\.mummichogRTSigPvals <- function\\(", src)
anchor <- grep("json.mat <- RJSONIO::toJSON(json.res);", src, fixed = TRUE)
anchor <- anchor[anchor > start][1]
if (length(start) != 1L || is.na(anchor)) stop("anchors not found in upstream source")
eval(parse(text = sub("^\\.compute\\.mummichogRTSigPvals <- ",
                      ".compute_sig_pvals_no_json <- ",
                      c(src[start:(anchor - 1L)], "  return(mSetObj)", "}"))),
     envir = up)

mSet <- c(list(
  dataSet = list(ref_mzlist = inputs$ref_mzlist, N = length(inputs$input_mzlist),
                 input_mzlist = inputs$input_mzlist, mum.type = "list"),
  paramSet = list(mumRT = TRUE, version = "v2", anal.type = "mummichog"),
  mum_nm_csv = "targeted.csv", lib.organism = "targeted"),
  inputs[c("pathways", "input_ecpdlist", "total_matched_ecpds", "mz2ec_dict",
           "ecpd_cpd_dict", "cpd_ecpd_dict", "ec_exp", "cpd_form_dict")])

# A warning as "<call>: <message>", so its origin (e.g. cov2cor) is recorded.
warning_text <- function(w) {
  paste0(paste(deparse(conditionCall(w)), collapse = " "), ": ",
         conditionMessage(w))
}

upstream_warnings <- character()
withr::with_dir(tempdir(), invisible(capture.output(withCallingHandlers({
  mSet <- up$.perform.mummichogRTPermutations(mSet, PERM_NUM)
  mSet <- up$.compute_sig_pvals_no_json(mSet)
}, warning = function(w) {
  upstream_warnings <<- c(upstream_warnings, warning_text(w))
  invokeRestart("muffleWarning")
}))))
if (!is.list(mSet)) stop("upstream statistics failed on the targeted fixture")

# ---- diagnostics: the same fit, inspected (not upstream output) ----------------
perm_minus <- abs(0.9999999999 - unlist(mSet$perm_record))
fit_warnings <- character()
fit <- withCallingHandlers(
  fitdistrplus::fitdist(perm_minus, distr = "gamma", method = "mle",
                        lower = c(0, 0), start = list(scale = 1, shape = 1)),
  warning = function(w) {
    fit_warnings <<- c(fit_warnings, warning_text(w))
    invokeRestart("muffleWarning")
  })

ref <- list(
  metaboanalystr_commit = sha,
  upstream_functions    = c(".perform.mummichogRTPermutations",
                            "ComputeMummichogRTPermPvals", "make_ecpdlist",
                            ".compute.mummichogRTSigPvals (up to its JSON export)",
                            ".save.mummichog.restable"),
  r_version             = R.version.string,
  rng_kind              = RNGkind(),
  locale_collate        = Sys.getlocale("LC_COLLATE"),
  fitdistrplus_version  = as.character(utils::packageVersion("fitdistrplus")),
  generated_utc         = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  perm_num              = PERM_NUM,
  inputs                = inputs,
  upstream = list(perm_record = mSet$perm_record, perm_hits = mSet$perm_hits,
                  pvals = mSet$pvals, mummi.resmat = mSet$mummi.resmat,
                  path.nms = mSet$path.nms, path.pval = mSet$path.pval,
                  restable = up$csv_store[["targeted.csv"]],
                  warnings = upstream_warnings),
  diagnostics = list(estimate = fit$estimate, sd = fit$sd,
                     convergence = fit$convergence, loglik = fit$loglik,
                     warnings = fit_warnings)
)
out <- file.path(here, "stats_targeted_reference.rds")
saveRDS(ref, out, version = 2)
message("wrote ", out)
res <- ref$upstream$mummi.resmat
est <- ref$diagnostics$estimate
fet <- ref$upstream$pvals[rownames(res), "FET"]
print(signif(cbind(FET = fet, Gamma = res[, "Gamma"],
  as_scale = 1 - stats::pgamma(1 - fet, shape = est["shape"], scale = est["scale"]),
  holm = stats::p.adjust(res[, "FET"], "holm"),
  bh = stats::p.adjust(res[, "FET"], "BH")), 5))
message("fit: shape ", signif(est["shape"], 6), ", scale ", signif(est["scale"], 6),
        ", convergence ", ref$diagnostics$convergence, ", sd ",
        paste(signif(ref$diagnostics$sd, 4), collapse = "/"), ", ",
        length(fit_warnings), " fit warning(s), ", length(upstream_warnings),
        " upstream warning(s)")
