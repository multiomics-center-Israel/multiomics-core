#!/usr/bin/env Rscript
#' Validate the pinned mummichog engine against a MetaboAnalyst run
#'
#' One-off analysis (NOT a pipeline stage): run our version-pinned mummichog
#' 2.7.0 engine (R/domain/metabolomics/06c) on the SAME peak input + model a
#' technician submitted to MetaboAnalyst for the 24h High-Light vs Low-Light
#' contrast, then compare pathway p-values to her MetaboAnalyst output.
#'
#' Apples-to-apples by design: same features + stats + metabolic model (cre), so
#' any difference is engine-vs-engine (standalone mummichog 2.7.0 vs
#' MetaboAnalyst's reimplementation), not input drift.
#'
#' It does NOT modify the engine — it only calls the existing engine/model
#' functions and reads the two result tables.
#'
#' ---------------------------------------------------------------------------
#' USAGE
#'   1. Build the pinned venv once:  make setup   (sets MUMMICHOG_PYTHON)
#'   2. Edit the CONFIG block below (or pass the same names as VAR=value env
#'      vars), then:
#'        Rscript scripts/validate_mummichog_vs_metaboanalyst.R
#'   Outputs (aggregate pathway stats only — no raw sample data) land in OUT_DIR:
#'     - comparison_pathways.csv   merged per-pathway p-values + ranks
#'     - comparison_summary.md     overlap + correlation metrics + caveats
#' ---------------------------------------------------------------------------

# ==== CONFIG (edit these, or override with same-named env vars) ==============
cfg <- list(
  # Technician's MetaboAnalyst PEAK INPUT for 24h HL-LL (the table she
  # submitted): columns m.z, r.t, mode, t.score, p.value. Set the path to wherever
  # the MetaboAnalyst files live on your machine.
  ma_input   = Sys.getenv("MA_INPUT",  "MetaboAnalyst Files/Mummichog/24h HL-LL/FOR UPLOAD HL vs LL 24h.csv"),
  # Technician's MetaboAnalyst mummichog pathway OUTPUT to compare against.
  ma_output  = Sys.getenv("MA_OUTPUT", "MetaboAnalyst Files/Mummichog/24h HL-LL/mummichog_pathway_enrichment_mummichog.csv"),

  # cre model: the published model_ref we use (fetched + sha256-verified +
  # cached). Leave model_json empty to use model_ref; set model_json to a local
  # path to skip the download.
  model_url    = Sys.getenv("MODEL_URL",    "https://github.com/multiomics-center-Israel/multiomics-annotation-prep/releases/download/cre_kegg_20260711/cre_kegg_20260711.json"),
  model_sha256 = Sys.getenv("MODEL_SHA256", "c403c96fbec8df9ae34b828fec01270c8ea3940acc36e4e5ff770868dc8b912b"),
  model_json   = Sys.getenv("MODEL_JSON",   ""),

  # Match the parameters the technician used in MetaboAnalyst.
  mode         = Sys.getenv("MODE",         "pos_default"),  # pos_default|positive|negative
  ppm          = as.numeric(Sys.getenv("PPM",          "10")),
  permutations = as.numeric(Sys.getenv("PERMUTATIONS", "100")),
  cutoff       = suppressWarnings(as.numeric(Sys.getenv("CUTOFF", ""))),  # NA -> mummichog's auto cutoff

  # MetaboAnalyst can run MIXED ionization (a per-feature `mode` column), but the
  # pinned engine (standalone mummichog 2.7.0) applies a SINGLE global mode. Set
  # MODE_FILTER (e.g. "Positive") to run on that subset only — and compare against
  # a MATCHING single-mode MetaboAnalyst reference. Leave "" for single-mode input.
  mode_filter = Sys.getenv("MODE_FILTER", ""),

  # If you already have our pathway table from a prior run, point here to SKIP
  # running the engine and just compare (e.g. .../mcg_pathwayanalysis_*.tsv).
  our_pathways = Sys.getenv("OUR_PATHWAYS", ""),

  # Reference p-value / pathway columns. REF_PCOL defaults to "P(Fisher)" — the
  # confirmed MetaboAnalyst mummichog p-value for this run (== integ
  # Mummichog_Pvals). NB our engine reports mummichog's PERMUTATION p, not Fisher,
  # so expect rank agreement rather than identical values. Override if needed.
  ref_pcol    = Sys.getenv("REF_PCOL",    "P(Fisher)"),
  ref_pathcol = Sys.getenv("REF_PATHCOL", ""),

  out_dir = Sys.getenv("OUT_DIR", "mummichog_validation_24h_HL_LL")
)

# ==== load the engine (no changes to it) ====================================
# Resolve repo root relative to this script so it runs from anywhere.
this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo <- if (!is.na(this_file) && nzchar(this_file)) {
  normalizePath(file.path(dirname(this_file), ".."))
} else {
  normalizePath(".")
}
for (f in c("R/core/00_paths.R",
            "R/domain/metabolomics/06c_mummichog_pinned.R",
            "R/domain/metabolomics/06d_mummichog_model.R")) {
  source(file.path(repo, f))
}

# ==== small helpers =========================================================
# Normalise a pathway name for joining (case/punctuation/spacing-insensitive).
norm_name <- function(x) {
  x <- tolower(trimws(as.character(x)))
  trimws(gsub("\\s+", " ", gsub("[^a-z0-9]+", " ", x)))
}

# Read a delimited table, choosing CSV vs TSV by file extension.
read_table_any <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    readr::read_csv(path, show_col_types = FALSE)
  } else {
    readr::read_tsv(path, show_col_types = FALSE)
  }
}

# Find the first present column among candidates (then a regex fallback).
find_col <- function(df, candidates, regex = NULL) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) return(hit[[1]])
  if (!is.null(regex)) {
    i <- grep(regex, names(df), ignore.case = TRUE)
    if (length(i)) return(names(df)[i[[1]]])
  }
  NULL
}

# ==== 1. our pinned run (or read an existing result) ========================
if (nzchar(cfg$our_pathways)) {
  message("Reading existing pinned pathway table: ", cfg$our_pathways)
  our <- read_table_any(cfg$our_pathways)
} else {
  stopifnot("MetaboAnalyst input not found" = file.exists(cfg$ma_input))
  peaks <- read_table_any(cfg$ma_input)

  # Ionization mode: the pinned engine runs ONE global mode and ignores any
  # per-feature `mode` column. If the input is mixed-mode, report the split and
  # either subset to MODE_FILTER (honest single-mode run) or warn loudly.
  mode_col <- find_col(peaks, c("mode", "Mode", "ion_mode", "ionization"),
                       "^mode$|ion")
  if (!is.null(mode_col)) {
    modes <- tolower(trimws(as.character(peaks[[mode_col]])))
    tab   <- table(modes)
    message("Per-feature mode column '", mode_col, "': ",
            paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", "))
    if (nzchar(cfg$mode_filter)) {
      keep <- modes == tolower(trimws(cfg$mode_filter))
      message("MODE_FILTER='", cfg$mode_filter, "' -> keeping ", sum(keep), "/",
              length(keep), " features.")
      peaks <- peaks[keep, , drop = FALSE]
      if (nrow(peaks) == 0L) {
        stop("No features match MODE_FILTER='", cfg$mode_filter, "'.")
      }
    } else if (length(tab) > 1L) {
      warning("Input is MIXED-mode but MODE_FILTER is unset: running all features ",
              "under a single '", cfg$mode, "' mode is NOT apples-to-apples with a ",
              "mixed MetaboAnalyst run. Set MODE_FILTER (e.g. Positive) and compare ",
              "against a matching single-mode reference, or add mixed-mode support (B3).",
              call. = FALSE)
    }
  }

  mz_c <- find_col(peaks, c("m.z", "m/z", "mz", "mass"), "^m[./]?z$")
  rt_c <- find_col(peaks, c("r.t", "rt", "retention_time", "RT"), "^r[._]?t$|retention")
  p_c  <- find_col(peaks, c("p.value", "pvalue", "p_value", "p-value", "P.Value"), "^p[._-]?val")
  s_c  <- find_col(peaks, c("t.score", "tscore", "t.stat", "statistic", "logFC", "stat"), "score|stat|fc")
  if (is.null(mz_c) || is.null(p_c) || is.null(s_c)) {
    stop("Could not locate m/z, p-value and statistic columns in ", cfg$ma_input,
         " (found: ", paste(names(peaks), collapse = ", "), ")")
  }
  if (is.null(rt_c)) {
    stop("No retention-time column in the MetaboAnalyst input; the pinned v2 ",
         "engine needs RT. Re-export the peak list with r.t, or supply OUR_PATHWAYS.")
  }

  stats_tbl <- data.frame(
    feature_id     = if ("CompoundID" %in% names(peaks)) as.character(peaks$CompoundID)
                     else paste0("f", seq_len(nrow(peaks))),
    mz             = as.numeric(peaks[[mz_c]]),
    retention_time = as.numeric(peaks[[rt_c]]),
    p_value        = as.numeric(peaks[[p_c]]),
    statistic      = as.numeric(peaks[[s_c]]),
    stringsAsFactors = FALSE
  )

  # cre model: model_ref (fetch + verify + cache) unless a local model_json is set.
  network <- if (nzchar(cfg$model_json)) {
    cfg$model_json
  } else {
    mmc_resolve_model(list(url = cfg$model_url, sha256 = cfg$model_sha256))
  }

  ensure_dir(cfg$out_dir)
  run_dir <- file.path(normalizePath(cfg$out_dir), "our_run")
  infile  <- write_mummichog_input(
    stats_tbl = stats_tbl,
    path      = file.path(run_dir, "input.tsv"),
    mz_col = "mz", rt_col = "retention_time",
    p_col = "p_value", stat_col = "statistic", id_col = "feature_id"
  )
  files <- run_mummichog_v2(
    infile         = infile,
    out_dir        = file.path(run_dir, "v2"),
    project        = "mcg_validation",
    network        = network,
    mode           = cfg$mode,
    instrument_ppm = cfg$ppm,
    permutations   = cfg$permutations,
    cutoff         = if (is.na(cfg$cutoff)) NULL else cfg$cutoff
  )
  our <- read_mummichog_pathways(files)
}

# our columns arrive verbatim: pathway, overlap_size, pathway_size, `p-value`
our_p_c <- find_col(our, c("p-value", "p.value", "pvalue", "p_value"), "^p[._-]?val")
our_df  <- data.frame(
  key      = norm_name(our$pathway),
  pathway  = as.character(our$pathway),
  our_p    = suppressWarnings(as.numeric(our[[our_p_c]])),
  stringsAsFactors = FALSE
)
our_df <- our_df[is.finite(our_df$our_p), , drop = FALSE]

# ==== 2. the technician's MetaboAnalyst output ==============================
stopifnot("MetaboAnalyst output not found" = file.exists(cfg$ma_output))
ref <- read_table_any(cfg$ma_output)

# Pathway name: usually the first (row-name) column; MetaboAnalyst writes it as
# an unnamed/`...1`/`X` column.
ref_path_c <- if (nzchar(cfg$ref_pathcol)) cfg$ref_pathcol else
  find_col(ref, c("...1", "X", "Pathway", "pathway", "Name", "name"), "path|name")
if (is.null(ref_path_c)) ref_path_c <- names(ref)[1]

# mummichog permutation p-value: MetaboAnalyst names it Gamma / Mummichog.P.value
# (NOT the Fisher FET/EASE columns). Auto-detect, overridable via REF_PCOL.
ref_p_c <- if (nzchar(cfg$ref_pcol)) cfg$ref_pcol else
  find_col(ref,
           c("Mummichog.P.value", "mummichog.pvalue", "Mummichog_Pvalue",
             "Gamma", "gamma", "Pval_Mummichog", "P.value", "Pvalue"),
           "mummich|gamma|p[._-]?val")
if (is.null(ref_p_c)) {
  stop("Could not find a mummichog p-value column in ", cfg$ma_output,
       " (columns: ", paste(names(ref), collapse = ", "),
       "). Set REF_PCOL to the right one.")
}
message("Reference pathway column: '", ref_path_c, "'  p-value column: '", ref_p_c, "'")

ref_df <- data.frame(
  key     = norm_name(ref[[ref_path_c]]),
  their_p = suppressWarnings(as.numeric(ref[[ref_p_c]])),
  stringsAsFactors = FALSE
)
ref_df <- ref_df[is.finite(ref_df$their_p), , drop = FALSE]

# ==== 3. compare ============================================================
merged <- merge(our_df, ref_df, by = "key")
merged$our_rank   <- rank(merged$our_p,   ties.method = "average")
merged$their_rank <- rank(merged$their_p, ties.method = "average")
merged <- merged[order(merged$our_p), c("pathway", "our_p", "their_p",
                                        "our_rank", "their_rank")]

top_overlap <- function(n) {
  a <- head(our_df$key[order(our_df$our_p)], n)
  b <- head(ref_df$key[order(ref_df$their_p)], n)
  inter <- length(intersect(a, b)); uni <- length(union(a, b))
  list(n = n, shared = inter, jaccard = if (uni) inter / uni else NA_real_)
}
n_matched <- nrow(merged)
spearman  <- if (n_matched >= 3) suppressWarnings(
  cor(merged$our_p, merged$their_p, method = "spearman")) else NA_real_
pearson_nlp <- if (n_matched >= 3) suppressWarnings(
  cor(-log10(merged$our_p), -log10(merged$their_p))) else NA_real_
t10 <- top_overlap(10); t20 <- top_overlap(20)

# ==== 4. write outputs (aggregate pathway stats only) =======================
ensure_dir(cfg$out_dir)
readr::write_csv(merged, file.path(cfg$out_dir, "comparison_pathways.csv"))

only_ours   <- setdiff(our_df$key, ref_df$key)
only_theirs <- setdiff(ref_df$key, our_df$key)
fmt <- function(x) if (is.na(x)) "NA" else formatC(x, digits = 3, format = "f")

summary_md <- c(
  "# Mummichog validation — 24h HL vs LL (ours vs MetaboAnalyst)",
  "",
  "Engine-vs-engine on the technician's MetaboAnalyst peak input + the cre model.",
  "",
  "## Coverage",
  sprintf("- Pathways — ours: %d, MetaboAnalyst: %d, matched by name: %d",
          nrow(our_df), nrow(ref_df), n_matched),
  sprintf("- Only in ours: %d; only in MetaboAnalyst: %d",
          length(only_ours), length(only_theirs)),
  "",
  "## Agreement (on matched pathways)",
  sprintf("- Spearman rho (p-values): %s", fmt(spearman)),
  sprintf("- Pearson r (-log10 p): %s", fmt(pearson_nlp)),
  sprintf("- Top-10 overlap: %d shared (Jaccard %s)", t10$shared, fmt(t10$jaccard)),
  sprintf("- Top-20 overlap: %d shared (Jaccard %s)", t20$shared, fmt(t20$jaccard)),
  "",
  "## Notes / expected differences",
  "- mummichog 2.7.0 estimates significance by UNSEEDED permutation, so p-values",
  "  vary run-to-run; judge by rank agreement + top-N overlap, not exact equality.",
  "- Standalone mummichog 2.7.0 != MetaboAnalyst's reimplementation (permutation",
  "  handling, EASE/gamma, and pathway-library versioning differ).",
  "- Compared against the MetaboAnalyst mummichog/Gamma p-value column",
  sprintf("  ('%s'); set REF_PCOL to change.", ref_p_c),
  ""
)
writeLines(summary_md, file.path(cfg$out_dir, "comparison_summary.md"))

cat(paste(summary_md, collapse = "\n"), "\n")
cat("\nWrote:\n  ", file.path(cfg$out_dir, "comparison_pathways.csv"),
    "\n  ", file.path(cfg$out_dir, "comparison_summary.md"), "\n")
