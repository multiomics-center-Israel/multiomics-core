# tests/testthat/fixtures/mummichog_ma_parity/validate_real_run.R
#
# Real-data validation of the MetaboAnalyst-compatible mummichog engine
# (mummichog.engine: metaboanalyst) on an EXISTING pipeline run, without
# re-running anything upstream and without touching that run.
#
#   Rscript tests/testthat/fixtures/mummichog_ma_parity/validate_real_run.R \
#       --config /path/to/mum09_config.yaml \
#       --out    /path/to/separate/validation_dir \
#       [--store _targets_<...>]  [--polarity-column "Reference Ion"]  [--no-report]
#
# Run from the multiomics-core root. What it does:
#   1. reads the run's config and its targets store READ-ONLY (targets::tar_read):
#      the preprocessed data (metab_pre), the DE results (metab_de_res) and the
#      other inputs of the report and Shiny payload, exactly as the pipeline
#      stored them;
#   2. overrides, IN MEMORY only, mummichog.engine = metaboanalyst and
#      ionization_mode = mixed (preprocessing, normalization and DE keys are not
#      touched; the config file on disk is not modified); a copy of the
#      effective mummichog block is written to <out>/mummichog_ma_config.yaml;
#   3. preflights the inputs and prints PASS / FAIL lines (no data values);
#   4. runs the engine, the report sections, their exports, the Shiny payload
#      and (unless --no-report) the HTML report, all into <out>/metabolomics/
#      and <out>/, never into the original run directory.
#
# It writes nothing into the store or the original results, and it never calls
# tar_make(). The DE results it uses are the stored object itself, so they are
# by construction the ones the existing python-engine results were built on;
# the digest printed at the end lets you confirm that against the store later.

args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(flag, " needs a value")
  args[[i + 1L]]
}
config_path <- opt("--config")
out_root    <- opt("--out")
if (is.null(config_path) || is.null(out_root)) {
  stop("usage: Rscript tests/testthat/fixtures/mummichog_ma_parity/validate_real_run.R --config <config.yaml> ",
       "--out <separate dir> [--store <store>] [--polarity-column <col>] [--no-report]")
}
do_report <- !("--no-report" %in% args)
if (!file.exists("_targets.R") || !dir.exists("R")) {
  stop("run this from the multiomics-core root")
}

# Source R/ in the pipeline's load order (core -> services -> domain -> modules).
for (layer in c("core", "services", "domain", "modules")) {
  for (f in sort(list.files(file.path("R", layer), pattern = "\\.R$",
                            full.names = TRUE, recursive = TRUE))) {
    source(f)
  }
}

results <- list()
check <- function(name, ok, detail = "") {
  results[[name]] <<- isTRUE(ok)
  cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", name,
              if (nzchar(detail)) paste0(" -- ", detail) else ""))
  invisible(isTRUE(ok))
}
stop_if_failed <- function() {
  if (!all(unlist(results))) stop("preflight failed; nothing was run", call. = FALSE)
}

# ---- config and store, read-only ---------------------------------------------
config <- validate_config(load_config(config_path))
store <- opt("--store", targets_store_for(config))
check("targets store exists", dir.exists(store), store)
stop_if_failed()
read_t <- function(name) targets::tar_read_raw(name, store = store)
# tar_meta(), not tar_objects(): the latter omits format = "file" targets.
has_t <- function(name) name %in% targets::tar_meta(store = store, fields = "name")$name

need <- c("metab_pre", "metab_de_res", "metab_inputs", "metab_enrichment_res")
for (n in need) check(paste("stored target", n), has_t(n))
stop_if_failed()

cfg_ma <- config
mm <- cfg_ma$modes$metabolomics$enrichment$mummichog %||% list()
pol_col <- opt("--polarity-column", mm$polarity_column)
mm$enabled <- TRUE
mm$engine <- "metaboanalyst"
mm$ionization_mode <- "mixed"
mm$polarity_column <- pol_col
cfg_ma$modes$metabolomics$enrichment$mummichog <- mm
check("only the mummichog block differs from the run's config", identical(
  { a <- config; a$modes$metabolomics$enrichment$mummichog <- NULL; a },
  { b <- cfg_ma; b$modes$metabolomics$enrichment$mummichog <- NULL; b }))

pre    <- read_t("metab_pre")
de_res <- read_t("metab_de_res")
rd <- pre$row_data
check("feature annotations have an m/z column",
      !is.null(.mmc_find_col(rd, c("m/z", "mz", "m.z", "MZ", "Mass", "m.z."), "^m[./]z")))
check("feature annotations have a retention-time column",
      !is.null(.mmc_find_col(rd, c("RT [min]", "rt", "RT", "retention_time",
                                   "RetentionTime", "RT..min.", "RT.min."),
                             "^RT|retention")))
pol <- tryCatch(mmc_feature_polarity(rd, pol_col), error = function(e) e)
check("every feature has a readable polarity", !inherits(pol, "error"),
      if (inherits(pol, "error")) conditionMessage(pol) else
        sprintf("column '%s'", pol_col))
if (!inherits(pol, "error")) {
  check("both polarities are present", all(c("positive", "negative") %in% pol))
}
de_tables <- de_res$de_tables %||% list()
check("DE results have per-contrast tables", length(de_tables) > 0,
      paste(names(de_tables), collapse = ", "))
check("every DE table has feature_id, P.Value and logFC",
      length(de_tables) > 0 && all(vapply(de_tables, function(t)
        all(c("feature_id", "P.Value", "logFC") %in% names(t)), logical(1))))
model <- tryCatch(mmc_select_model(mm, organism = cfg_ma$modes$metabolomics$organism,
                                   cache_dir = "envs/mummichog-models"),
                  error = function(e) e)
check("the configured model resolves to a readable JSON",
      !inherits(model, "error") && !identical(model, "human_mfn") && file.exists(model),
      if (inherits(model, "error")) conditionMessage(model) else basename(model))
check("fitdistrplus 1.2.6 installed (the parity fixtures' version)",
      requireNamespace("fitdistrplus", quietly = TRUE) &&
        identical(as.character(utils::packageVersion("fitdistrplus")), "1.2.6"))
stop_if_failed()

# ---- run, into the separate directory only -------------------------------------
out_dir <- file.path(out_root, "metabolomics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
yaml::write_yaml(list(mummichog = mm), file.path(out_root, "mummichog_ma_config.yaml"))
de_digest <- digest::digest(de_res)

files <- mod_mummichog_metaboanalyst(pre, de_res, cfg_ma, out_dir)
pathways <- read_mummichog_pathways_by_contrast(files)
enrichment_res <- read_t("metab_enrichment_res")
sections <- mod_mummichog_sections(pathways, files, pre, de_res, cfg_ma,
                                   enrichment_res = enrichment_res)
exports <- mod_mummichog_exports(sections, files, out_dir)

opt_t <- function(name) if (has_t(name)) read_t(name) else NULL
fs <- opt_t("metab_feature_sel_res")
shiny_file <- file.path(out_dir, sprintf("shiny_payload_metabolomics_%s_ma.rds",
                                         get_run_name(cfg_ma)))
save_shiny_payload_metabolomics(
  pre = pre, de_res = de_res, inputs = read_t("metab_inputs"), config = cfg_ma,
  pca_res = opt_t("metab_qc_pre_obj"), clustering_res = opt_t("metab_clustering_obj"),
  rf_res = if (!is.null(fs)) fs$rf else NULL,
  plsda_res = if (!is.null(fs)) fs$plsda else NULL,
  enrichment_res = enrichment_res, include_legacy = TRUE,
  xlsx_files = opt_t("metab_final_results"), mummichog_pathways = pathways,
  out_file = shiny_file)

report_file <- NULL
if (do_report) {
  report_file <- tryCatch(mod_metabolomics_report(
    pre = pre, qc_res = opt_t("metab_qc_pre_obj"), de_res = de_res,
    clustering_res = opt_t("metab_clustering_obj"), feature_sel_res = fs,
    enrichment_res = enrichment_res, network_res = opt_t("metab_network"),
    config = cfg_ma, out_dir = out_dir, qc_comparison_file = NULL,
    qc_suite_files = NULL, commentary_file = opt_t("metab_commentary"),
    mummichog_sections = sections),
    error = function(e) { message("report: ", conditionMessage(e)); NULL })
}

# ---- post-run checks (structure only) -------------------------------------------
cat("\n")
check("every DE contrast produced a pathway table",
      setequal(names(pathways), names(de_tables)),
      paste(names(pathways), collapse = ", "))
for (ct in names(pathways)) {
  tb <- sections[[ct]]$table
  check(sprintf("[%s] section built with the metaboanalyst engine", ct),
        identical(sections[[ct]]$engine, "metaboanalyst"))
  check(sprintf("[%s] table has FET p.value and the Holm columns", ct),
        !is.null(tb) && all(c("p.value", "p.value EASE", "p.value Gamma",
                              "p.value empirical", "adj.p Holm FET",
                              "adj.p Holm EASE", "adj.p Holm Gamma") %in% names(tb)) &&
          identical(tb$p.value, pathways[[ct]]$FET))
  check(sprintf("[%s] no BH/FDR-labelled column", ct),
        !any(grepl("BH|FDR|fdr|padj", names(tb))))
  check(sprintf("[%s] GSEA and evidence skipped", ct),
        is.null(sections[[ct]]$gsea) && is.null(sections[[ct]]$evidence))
  g <- attr(pathways[[ct]], "ma_info")$gamma
  cat(sprintf("       Gamma fit: %s\n", if (isTRUE(g$fitted)) "fitted" else
    paste0("NOT fitted (Gamma column NA): ", g$error %||% "too few hits")))
}
payload <- readRDS(shiny_file)
check("Shiny payload has a mummichog entry per contrast",
      setequal(names(payload$mummichog), names(pathways)))
check("Shiny mummichog keys unchanged",
      all(vapply(payload$mummichog, function(s)
        identical(names(s), .METAB_SHINY_MUMMICHOG_KEYS), logical(1))))
if (do_report) check("HTML report rendered", !is.null(report_file) &&
                       all(file.exists(report_file)))
check("DE results unchanged by the run", identical(digest::digest(de_res), de_digest))

cat("\nOutputs under: ", normalizePath(out_root), "\n", sep = "")
cat("DE results digest (store object metab_de_res): ", de_digest, "\n", sep = "")
if (!all(unlist(results))) quit(status = 1)
