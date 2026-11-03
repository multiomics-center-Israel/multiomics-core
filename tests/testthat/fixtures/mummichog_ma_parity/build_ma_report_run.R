# tests/testthat/fixtures/mummichog_ma_parity/build_ma_report_run.R
#
# Builds a COMPLETE report deliverable for an existing run with the Mummichog
# section replaced by the MetaboAnalyst-compatible joint-polarity engine, in a
# NEW run directory, leaving the original run untouched. Run MANUALLY, from the
# multiomics-core root:
#
#   Rscript tests/testthat/fixtures/mummichog_ma_parity/build_ma_report_run.R \
#       --config /path/to/run_config.yaml --dest /path/to/Results_<run>_MA \
#       [--store _targets_<...>] [--polarity-column "Reference Ion"]
#
# What it does:
#   1. snapshots every file of the original run directory (size + mtime);
#   2. copies the run into --dest, EXCEPT the old Mummichog outputs
#      (metabolomics/mummichog_pinned/), the old report (report_metabolomics.*)
#      and the old Shiny payload (metabolomics/shiny_payload_metabolomics_*.rds),
#      which are rebuilt; everything else — preprocessing, QC, normalization,
#      DE, clustering, enrichment, final Excel files — is carried over as is;
#   3. reads the stored pipeline objects read-only and points their file paths
#      at the copy, so the report reads and writes only inside --dest (the
#      report embeds figures from those paths and writes heatmaps beside them);
#   4. runs the MetaboAnalyst-compatible engine (engine / ionization_mode
#      overridden in memory only), its sections and exports, the Shiny payload
#      and the HTML report into --dest;
#   5. writes MUMMICHOG_UPDATE_NOTE.md (methods / provenance);
#   6. checks: original unchanged; no old Mummichog artefact in --dest; the new
#      report states the new method and none of the old per-polarity wording;
#      every report heading of the original report is still present, except the
#      Mummichog sub-tabs this engine does not produce; Excel files present.
#
# --dest should be a SIBLING of the original run directory: the report finds
# its interactive cutoff-panel helpers by walking up from where it is rendered.

args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(flag, " needs a value")
  args[[i + 1L]]
}
config_path <- opt("--config"); dest <- opt("--dest")
if (is.null(config_path) || is.null(dest)) {
  stop("usage: build_ma_report_run.R --config <config.yaml> --dest <new run dir> ",
       "[--store <store>] [--polarity-column <col>] [--rerender]")
}
# --rerender: --dest is an earlier build of this script. Nothing is copied and
# the Mummichog engine is not re-run: its results already in --dest are read
# back, and only the sections, exports, report, Shiny payload and note are
# rebuilt there.
rerender <- "--rerender" %in% args
# --rerun-mummichog (with --rerender): also re-run the MetaboAnalyst-compatible
# stage (ORA + GSEA) in --dest. The ORA result tables it had are kept aside and
# must come back identical.
rerun_mm <- "--rerun-mummichog" %in% args
if (rerun_mm && !rerender) stop("--rerun-mummichog needs --rerender")
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
halt_if_failed <- function() if (!all(unlist(results))) stop("checks failed; stopping", call. = FALSE)

config <- validate_config(load_config(config_path))
orig_given <- get_run_out_dir(config)
orig <- normalizePath(orig_given, mustWork = FALSE)
store <- opt("--store", targets_store_for(config))
check("original run directory exists", dir.exists(orig), orig)
check("targets store exists", dir.exists(store), store)
if (rerender) {
  check("--dest holds an earlier build (MetaboAnalyst-compatible results present)",
        file.exists(file.path(dest, "metabolomics", "mummichog_pinned", "contrasts.tsv")), dest)
} else {
  check("--dest does not exist yet", !file.exists(dest), dest)
}
halt_if_failed()
if (!rerender) dir.create(dest, recursive = TRUE)
dest <- normalizePath(dest)
check("--dest is outside the original run", !startsWith(paste0(dest, "/"), paste0(orig, "/")))
check("--dest is a sibling of the original run (report helpers resolve the same)",
      identical(dirname(dest), dirname(orig)))
halt_if_failed()

# ---- 1. snapshot of the original ------------------------------------------------
snapshot <- function(root) {
  f <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  i <- file.info(file.path(root, f))
  data.frame(file = f, size = i$size, mtime = as.numeric(i$mtime),
             stringsAsFactors = FALSE)
}
before <- snapshot(orig)

# ---- 2. copy, without the old Mummichog / report / Shiny payload ----------------
excluded <- "^metabolomics/mummichog_pinned/|^report_metabolomics\\.(html|Rmd)$|^metabolomics/shiny_payload_metabolomics_.*\\.rds$"
if (!rerender) {
  to_copy <- before$file[!grepl(excluded, before$file)]
  for (d in unique(dirname(file.path(dest, to_copy)))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(file.path(orig, to_copy), file.path(dest, to_copy), copy.date = TRUE)
  check("run copied (old Mummichog, report and Shiny payload left out)", all(ok),
        sprintf("%d file(s) not copied", sum(!ok)))
  halt_if_failed()
}

# ---- 3. stored objects, paths pointed at the copy --------------------------------
prefixes <- unique(c(orig, orig_given, sub("/$", "", orig_given)))
# Graphics objects are never descended into. A compound grid `unit` is a list
# whose `[[` returns another unit, so walking one recursed until the node stack
# overflowed; plots and grobs carry no output paths anyway. A depth cap guards
# against any other self-similar object.
.no_descend <- c("unit", "grob", "gList", "gTree", "gtable", "gg", "ggplot",
                 "theme", "patchwork", "Heatmap", "HeatmapList", "pheatmap",
                 "htmlwidget", "recordedplot", "trellis")
walk_paths <- function(x, f, depth = 0L) {
  if (is.character(x)) return(f(x))
  if (is.list(x) && depth < 40L && !inherits(x, .no_descend)) {
    for (i in seq_along(x)) {
      if (!is.null(x[[i]])) x[[i]] <- walk_paths(x[[i]], f, depth + 1L)
    }
  }
  x
}
in_orig <- function(x) {
  hit <- rep(FALSE, length(x))
  for (pre in prefixes) hit <- hit | (!is.na(x) & startsWith(x, paste0(pre, "/")))
  hit
}
relocate <- function(x) walk_paths(x, function(s) {
  for (pre in prefixes) {
    hit <- !is.na(s) & startsWith(s, paste0(pre, "/"))
    s[hit] <- paste0(dest, substring(s[hit], nchar(pre) + 1L))
  }
  s
})
left_in_orig <- function(x) {
  n <- 0L
  walk_paths(x, function(s) { n <<- n + sum(in_orig(s)); s })
  n
}
read_t <- function(n) targets::tar_read_raw(n, store = store)
# Read and relocated EAGERLY, one target at a time, so a failure names its target
# instead of surfacing inside whichever function first forces a lazy argument.
# Existence from the store's metadata: tar_objects() lists only targets stored
# as R objects, so format = "file" targets (the final Excel files, the AI
# commentary) were silently read as absent.
stored_names <- targets::tar_meta(store = store, fields = "name")$name
opt_t <- function(n) {
  if (!n %in% stored_names) return(NULL)
  tryCatch(relocate(read_t(n)), error = function(e)
    stop("relocating stored target '", n, "': ", conditionMessage(e), call. = FALSE))
}
pre <- relocate(read_t("metab_pre"))
de_raw <- read_t("metab_de_res")
de_res <- relocate(de_raw)
check("DE tables unchanged by path relocation", identical(de_res$de_tables, de_raw$de_tables))
inputs <- relocate(read_t("metab_inputs"))
enrichment_res <- opt_t("metab_enrichment_res")
fs <- opt_t("metab_feature_sel_res")
qc_obj <- opt_t("metab_qc_pre_obj")
clustering_obj <- opt_t("metab_clustering_obj")
final_results <- opt_t("metab_final_results")
network_obj <- opt_t("metab_network")
commentary <- opt_t("metab_commentary")
check("final-results Excel files found in the store",
      !is.null(final_results) && any(grepl("\\.xlsx$", final_results)) &&
        all(file.exists(final_results)))
check("AI commentary file found in the store (if the original run has one)",
      !("metab_commentary" %in% stored_names) ||
        (!is.null(commentary) && all(file.exists(commentary))))
stored <- list(pre = pre, de_res = de_res, inputs = inputs,
               enrichment = enrichment_res, feature_sel = fs, qc = qc_obj,
               clustering = clustering_obj, final_results = final_results,
               network = network_obj, commentary = commentary)
left <- vapply(stored, left_in_orig, integer(1))
check("no stored path still points into the original run", all(left == 0L),
      if (any(left > 0L)) paste(names(left)[left > 0L], collapse = ", ") else "")
halt_if_failed()

cfg_ma <- config
mm <- cfg_ma$modes$metabolomics$enrichment$mummichog %||% list()
mm$enabled <- TRUE; mm$engine <- "metaboanalyst"; mm$ionization_mode <- "mixed"
mm$polarity_column <- opt("--polarity-column", mm$polarity_column)
cfg_ma$modes$metabolomics$enrichment$mummichog <- mm

# ---- 4. engine, sections, exports, Shiny payload, report --------------------------
out_dir <- file.path(dest, "metabolomics")
ora_files <- function() {
  f <- list.files(file.path(out_dir, "mummichog_pinned"), recursive = TRUE,
                  full.names = TRUE)
  f[basename(f) %in% c("ma_pathway_results.tsv",
                       "mummichog_pathway_enrichment_mummichog.csv")]
}
ora_before <- if (rerun_mm) lapply(stats::setNames(nm = ora_files()), readLines) else NULL
files <- if (rerender && !rerun_mm) {
  list.files(file.path(out_dir, "mummichog_pinned"), recursive = TRUE, full.names = TRUE)
} else {
  mod_mummichog_metaboanalyst(pre, de_res, cfg_ma, out_dir)
}
if (rerun_mm) {
  ora_after <- lapply(stats::setNames(nm = ora_files()), readLines)
  check("ORA result tables identical to the previous build (GSEA added, ORA unchanged)",
        identical(ora_before, ora_after),
        sprintf("%d file(s) before, %d after", length(ora_before), length(ora_after)))
}
check("GSEA results written for every contrast",
      length(files[basename(files) == "ma_gsea_results.tsv"]) == length(de_res$de_tables))
pathways <- read_mummichog_pathways_by_contrast(files)
sections <- mod_mummichog_sections(pathways, files, pre, de_res, cfg_ma,
                                   enrichment_res = enrichment_res)
exports <- mod_mummichog_exports(sections, files, out_dir)
shiny_file <- file.path(out_dir, sprintf("shiny_payload_metabolomics_%s.rds",
                                         get_run_name(config)))
# The report does not read the Shiny payload, so it is rendered first and a
# payload failure cannot cost it.
report <- mod_metabolomics_report(
  pre = pre, qc_res = qc_obj, de_res = de_res,
  clustering_res = clustering_obj, feature_sel_res = fs,
  enrichment_res = enrichment_res, network_res = network_obj,
  config = cfg_ma, out_dir = out_dir, qc_comparison_file = NULL,
  qc_suite_files = NULL, commentary_file = commentary,
  mummichog_sections = sections)
shiny_ok <- tryCatch({
  save_shiny_payload_metabolomics(
    pre = pre, de_res = de_res, inputs = inputs, config = cfg_ma,
    pca_res = qc_obj, clustering_res = clustering_obj,
    rf_res = if (!is.null(fs)) fs$rf else NULL,
    plsda_res = if (!is.null(fs)) fs$plsda else NULL,
    enrichment_res = enrichment_res, include_legacy = TRUE,
    xlsx_files = final_results, mummichog_pathways = pathways,
    out_file = shiny_file)
  TRUE
}, error = function(e) { message("Shiny payload: ", conditionMessage(e)); FALSE })
check("Shiny payload written", shiny_ok)

# ---- 5. methods / provenance note --------------------------------------------------
info1 <- attr(pathways[[1]], "ma_info")
core_sha <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                     error = function(e) "unknown")
note <- c(
  "# Mummichog update — methods and provenance", "",
  sprintf("Run: %s (copied from `%s`; the original is unchanged).", get_run_name(config), orig), "",
  "Only the Mummichog (m/z-based pathway) analysis was redone. Preprocessing, QC, normalization, differential abundance, clustering, feature selection, identified-metabolite enrichment and the final Excel files are the original run's, carried over unchanged; the new analysis reads the stored differential-abundance results of that run.", "",
  "## What changed", "",
  "- **Before:** python mummichog 2.7.0, mixed-polarity data run **separately per ion polarity** (each polarity its own background, significant list and permutation null; separate result tabs; p-values not combined).",
  "- **Now:** ORA and GSEA from a native R reimplementation of MetaboAnalyst's MS Peaks-to-Pathways mummichog (MetaboAnalystR commit 398476ae2a0c996e390925fa62adc46b5ecde334), run **jointly over both polarities** per contrast: one peak list, each feature matched only to its own polarity's ion forms, a compound's positive- and negative-ion matches at similar retention times forming one empirical compound, one background and one permutation null. Each step was validated against that MetaboAnalystR version's own output.", "",
  "## How to read the results", "",
  "- `p.value` is the **unadjusted** Fisher exact test (FET) p-value per pathway; *nominal* means FET p < the cutoff.",
  "- `adj.p Holm *` columns are **Holm** (family-wise) adjustments over the reported pathways, as MetaboAnalyst computes them — **not** Benjamini–Hochberg FDR.",
  "- A pathway nominal by FET is a pathway-level hypothesis, **not** a result that is significant after multiple-testing correction.",
  "- GSEA (MetaboAnalyst's peaks GSEA) runs on the same joint-polarity empirical compounds and pathway sets as the ORA, over the complete ranked list; `P.Value` is its raw permutation p-value and `padj` its **Benjamini–Hochberg FDR** (ORA: Holm; GSEA: BH — as in MetaboAnalyst). NES sign is not a biological direction.",
  "- The Pathway Evidence section (python engine only) is not produced.", "",
  "## Provenance", "",
  sprintf("- multiomics-core commit: %s", core_sha),
  sprintf("- MetaboAnalystR reference: %s @ %s", info1$reference$repo, info1$reference$commit),
  sprintf("- fitdistrplus: %s", info1$fitdistrplus_version),
  sprintf("- Metabolic model: %s (version %s, sha256 %s)", info1$model$id, info1$model$version, info1$model$sha256),
  sprintf("- Parameters: ionization_mode %s (polarity column `%s`), %s ppm, significant features p <= %s, %s permutations (seed 123 + i), primary-ion rule %s, rt_frac %s, min pathway size %s",
          info1$params$ionization_mode, info1$params$polarity_column, info1$params$tolerance_ppm,
          info1$params$p_cutoff, info1$params$n_permutations,
          if (isTRUE(info1$params$force_primary_ion)) "on" else "off",
          info1$params$rt_frac, info1$params$min_pathway_size),
  sprintf("- Built: %s", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  "- Per-contrast outputs: `metabolomics/mummichog_pinned/<contrast>/metaboanalyst/` (incl. MetaboAnalyst's own result table, `mummichog_pathway_enrichment_mummichog.csv`).")
writeLines(note, file.path(dest, "MUMMICHOG_UPDATE_NOTE.md"))

# ---- 6. checks ------------------------------------------------------------------
cat("\n")
after <- snapshot(orig)
check("original run unchanged (same files, sizes, mtimes)", identical(before, after))

dest_files <- list.files(dest, recursive = TRUE)
mm_files <- dest_files[grepl("mummichog|mcg_", dest_files, ignore.case = TRUE)]
# Python-engine artefacts. GSEA exports are legitimate now (the new engine's
# own GSEA); evidence exports exist only for the python engine.
old <- mm_files[grepl("/v2/|mcg_|_positive_ions|_negative_ions|mummichog_evidence_",
                      mm_files)]
check("no old Mummichog artefact in the new run", length(old) == 0,
      if (length(old)) paste(utils::head(old, 3), collapse = "; ") else "")
check("every Mummichog file outside the new engine's tree is a new export",
      all(grepl("^metabolomics/mummichog_pinned/", mm_files) |
            grepl("^metabolomics/shiny_payload", mm_files) |
            grepl("^MUMMICHOG_UPDATE_NOTE", mm_files)))
check("final Excel files present", any(grepl("\\.xlsx$", dest_files)))

# Whitespace-normalised: pandoc wraps HTML text at 72 columns, so a phrase can
# be split across a line break in the file.
html_new <- gsub("\\s+", " ", paste(readLines(file.path(dest, "report_metabolomics.html"),
                                             warn = FALSE), collapse = "\n"))
check("report names the new engine", grepl("MetaboAnalyst-compatible engine", html_new, fixed = TRUE))
check("report states joint mixed-polarity analysis", grepl("analysed jointly", html_new, fixed = TRUE))
check("report states Holm, not BH/FDR", grepl("Holm, not Benjamini", html_new, fixed = TRUE))
check("report carries none of the old per-polarity Mummichog wording",
      !grepl("run separately per ion polarity", html_new, fixed = TRUE) &&
        !grepl("Mummichog ORA reports empirical raw pathway p-values", html_new, fixed = TRUE))
check("report never calls a Mummichog pathway significant after correction",
      !grepl("significant after (multiple-testing )?correction", html_new, perl = TRUE))

source(file.path("tests", "testthat", "fixtures", "mummichog_ma_parity", "check_ma_report.R"))
check_ma_report(orig, dest, check)

cat("\nNew run: ", dest, "\n", sep = "")
if (!all(unlist(results))) quit(status = 1)
