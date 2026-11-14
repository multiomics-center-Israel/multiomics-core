# R/domain/metabolomics/06n_mummichog_ma_results.R
#
# MetaboAnalyst-compatible mummichog, component 5 (domain side): the engine
# switch, the per-contrast driver that chains the validated components
# (06j library -> 06k peaks -> 06l ECs -> 06m statistics), and the writer and
# reader of the engine's result files.
#
# The stage module (05d) runs one contrast at a time through
# mmc_ma_run_contrast() and writes its results with write_mummichog_ma_results()
# under mummichog_pinned/<contrast>/metaboanalyst/ — beside where the python
# engine writes its v2/ tree, so the contrast map, the exports directory and
# group_mummichog_files_by_contrast() are shared. The report reads the results
# back with read_mummichog_ma_pathways(), which adds the columns the existing
# plot / table builders key on, keeping their contract.
#
# Nothing here changes a statistic: every number comes from the 06j-06m ports
# validated against MetaboAnalystR @ 398476ae.


#' Which mummichog engine a config selects
#'
#' The single resolver for `modes.metabolomics.enrichment.mummichog.engine`:
#' `"metaboanalyst"` (the MetaboAnalyst-compatible engine) or `"python"` (the
#' pinned mummichog 2.7.0 engine, kept as an explicit option). Every consumer
#' — the DAG, the stage, the report text, the polarity-label resolver, the
#' engine record — asks this rather than reading the key itself.
#'
#' A config that names no engine gets `"metaboanalyst"`. That is a default,
#' not a fallback: if the MetaboAnalyst-compatible engine's prerequisites are
#' missing, its stage stops; nothing reroutes the run to python.
#'
#' The resolver does not look at `enabled`. A disabled config runs no engine,
#' so its callers do not ask (see `mmc_engine_record()`).
#'
#' @param config Full pipeline config.
#' @return `"metaboanalyst"` or `"python"`. An unknown engine stops.
#' @examples
#' mmc_engine(list(modes = list(metabolomics = list(enrichment = list(
#'   mummichog = list(enabled = TRUE, engine = "python"))))))
mmc_engine <- function(config) {
  mm <- config$modes$metabolomics$enrichment$mummichog %||% list()
  if (is.null(mm$engine)) return(.MMC_DEFAULT_ENGINE)
  eng <- tolower(trimws(as.character(mm$engine)))
  if (length(eng) != 1L || is.na(eng) || !eng %in% c("python", "metaboanalyst")) {
    .mmc_stop("mummichog.engine must be \"metaboanalyst\" (the MetaboAnalyst-",
              "compatible engine, the default) or \"python\" (the pinned ",
              "mummichog 2.7.0 engine); got '",
              paste(mm$engine, collapse = ", "), "'.")
  }
  eng
}

# The engine an enabled config gets when it names none.
.MMC_DEFAULT_ENGINE <- "metaboanalyst"

# The engine record, at the top of mummichog_pinned/ (beside the contrast
# directories, which the stages own and clear).
.MMC_ENGINE_INFO_FILE <- "mummichog_engine.json"

#' The engine an enabled mummichog config selects, and where the choice came from
#'
#' What `write_mummichog_engine_info()` writes, and what the stage reads to say
#' that the default was applied. A disabled config selects nothing: the engine
#' key is neither resolved nor validated, and the result is NULL.
#'
#' @param config Full pipeline config.
#' @return NULL when mummichog is not enabled; otherwise a list with `engine`
#'   (as `mmc_engine()`), `engine_source` (`"config"` when the config names the
#'   engine, `"default"` when it does not) and `configured_engine` (the config
#'   value as written, NULL when absent).
#' @examples
#' mmc_engine_record(list(modes = list(metabolomics = list(enrichment = list(
#'   mummichog = list(enabled = TRUE))))))
mmc_engine_record <- function(config) {
  mm <- config$modes$metabolomics$enrichment$mummichog %||% list()
  if (!isTRUE(mm$enabled)) return(NULL)
  list(engine            = mmc_engine(config),
       engine_source     = if (is.null(mm$engine)) "default" else "config",
       configured_engine = mm$engine)
}

#' Record the selected mummichog engine beside the run's mummichog outputs
#'
#' Writes `mummichog_pinned/mummichog_engine.json`. The config snapshot in
#' `execution_info/` holds the config as written, so an engine that was
#' defaulted is not in it; this file names the engine either way and says
#' whether the config chose it. It records the selection for the run, not that
#' the stage completed — the engine's own result files are the evidence of
#' that. No timestamp, so an unchanged selection rewrites an identical file.
#'
#' @param config  Full pipeline config.
#' @param out_dir Metabolomics output dir; the file goes in
#'   `<out_dir>/mummichog_pinned/`.
#' @return The path of the record, for a `format = "file"` target; NULL
#'   (nothing written) when mummichog is not enabled.
#' @examples
#' \dontrun{
#' write_mummichog_engine_info(config, "out/metabolomics")
#' }
write_mummichog_engine_info <- function(config, out_dir) {
  rec <- mmc_engine_record(config)
  if (is.null(rec)) return(NULL)
  dir <- file.path(out_dir, "mummichog_pinned")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, .MMC_ENGINE_INFO_FILE)
  jsonlite::write_json(rec, path, auto_unbox = TRUE, pretty = TRUE,
                       null = "null")
  path
}

# The engine's own result files, per contrast.
.MMC_MA_SUBDIR        <- "metaboanalyst"
.MMC_MA_RESULTS_FILE  <- "ma_pathway_results.tsv"
.MMC_MA_UPSTREAM_CSV  <- "mummichog_pathway_enrichment_mummichog.csv"
.MMC_MA_FEATURES_FILE <- "ma_matched_features.tsv"
.MMC_MA_INPUT_FILE    <- "ma_input.tsv"
.MMC_MA_INFO_FILE     <- "ma_run_info.json"
.MMC_MA_GSEA_FILE     <- "ma_gsea_results.tsv"
.MMC_MA_FAILED_FILE   <- "ma_failed.json"
.MMC_MA_ALL_FILES     <- c(.MMC_MA_RESULTS_FILE, .MMC_MA_UPSTREAM_CSV,
                           .MMC_MA_FEATURES_FILE, .MMC_MA_INPUT_FILE,
                           .MMC_MA_INFO_FILE, .MMC_MA_GSEA_FILE,
                           .MMC_MA_FAILED_FILE)

#' Remove one contrast's MetaboAnalyst-compatible result files
#'
#' Removes exactly the engine's own file names (`.MMC_MA_ALL_FILES`) from
#' `dir`, never a pattern, so nothing else there is touched. Every write of a
#' contrast — a result or a failure record — starts here, so a rerun cannot
#' leave the previous run's files reading as current.
#'
#' @param dir The contrast's `metaboanalyst/` directory.
#' @return Character vector of the files removed (may be empty).
#' @noRd
.mmc_ma_clear_results <- function(dir) {
  old <- file.path(dir, .MMC_MA_ALL_FILES)
  old <- old[file.exists(old)]
  unlink(old)
  old
}

#' Record that one contrast failed in the MetaboAnalyst-compatible engine
#'
#' Clears the contrast's previous result files and writes `ma_failed.json`
#' (contrast label and error message) in their place. The record keeps the
#' contrast in the stage's tracked files, so the report layer still sees it:
#' it gets no section (there is no result table), and
#' `clear_orphan_mummichog_gsea_exports()` removes its earlier presentation
#' exports instead of leaving them on disk looking current.
#'
#' @param dir      The contrast's `metaboanalyst/` directory (created if absent).
#' @param contrast The DE contrast label.
#' @param error    The error message.
#' @return The path of the failure record.
#' @examples
#' \dontrun{
#' write_mummichog_ma_failure(tempfile(), "A_vs_B", "no usable statistic")
#' }
write_mummichog_ma_failure <- function(dir, contrast, error) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  .mmc_ma_clear_results(dir)
  path <- file.path(dir, .MMC_MA_FAILED_FILE)
  jsonlite::write_json(list(engine = "metaboanalyst", contrast = contrast,
                            failed = TRUE, error = error),
                       path, auto_unbox = TRUE, pretty = TRUE)
  path
}

#' Does a contrast's file list carry a failure record?
#'
#' @param files One contrast's files.
#' @return TRUE when `ma_failed.json` is among them.
#' @noRd
.mmc_ma_contrast_failed <- function(files) {
  any(basename(files) == .MMC_MA_FAILED_FILE)
}

#' Run the MetaboAnalyst-compatible mummichog analysis for one peak list
#'
#' Chains the validated components in MetaboAnalystR's order: prepare the peak
#' list, select the significant m/z, match against the library by polarity,
#' build empirical compounds, map pathways to ECs, permute, and compute the
#' pathway statistics. In `"mixed"` mode both polarities form ONE analysis.
#'
#' @param peaks Data frame: `feature_id`, `mz`, `rt`, `p_value`, `statistic`,
#'   and `polarity` in mixed mode (see `mmc_ma_prepare_peaks()`).
#' @param lib Library from `mmc_ma_library_from_model()`.
#' @param mode `"mixed"`, `"positive"` or `"negative"`.
#' @param ppm Mass tolerance in ppm.
#' @param cutoff Significance cutoff (inclusive, `p <= cutoff`).
#' @param n_perm Permutations (MetaboAnalyst default 100).
#' @param force_primary_ion Apply the primary-ion rule (default TRUE).
#' @param rt_frac Retention-time tolerance fraction (default 0.02).
#' @param min_lib Smallest pathway kept, in compounds (default 3).
#' @param gsea NULL to skip GSEA, or a list of `mmc_ma_gsea()` settings
#'   (`n_perm`, `min_size`, `max_size`, `seed`, `metric_label`). It runs on the
#'   same ECs and pathway sets as the ORA, after it, and changes nothing in it.
#'   GSEA is optional: an error in it is caught and returned as `gsea_error`,
#'   so the ORA result is kept.
#' @return A list: `prep`, `input_mzlist`, `matched`, `ecs`, `pathways`,
#'   `perms`, `stats`, `gsea` (NULL when skipped, failed, or nothing passed its
#'   size filter) and `gsea_error` (the error message when GSEA failed, else
#'   NULL).
#' @examples
#' \dontrun{
#' run <- mmc_ma_run_contrast(peaks, lib, "mixed", ppm = 10, cutoff = 0.05)
#' run$stats$table
#' }
mmc_ma_run_contrast <- function(peaks, lib, mode, ppm, cutoff, n_perm = 100,
                                force_primary_ion = TRUE, rt_frac = 0.02,
                                min_lib = 3, gsea = NULL) {
  prep <- mmc_ma_prepare_peaks(peaks, mode, rt_frac = rt_frac)
  sig <- mmc_ma_input_mzlist(prep, cutoff)
  matched <- mmc_ma_match_peaks(prep, lib, ppm = ppm)
  ecs <- mmc_ma_build_ecs(matched, prep, sig,
                          force_primary_ion = force_primary_ion)
  pw <- mmc_ma_pathway_ecs(lib, ecs, min_lib = min_lib)
  perms <- mmc_ma_permutations(prep, ecs, pw, sig, n_perm = n_perm)
  stats <- mmc_ma_pathway_stats(ecs, pw, perms)
  # GSEA is the optional half (as when fgsea is missing): its failure must not
  # take the already-computed ORA down with it.
  gsea_error <- NULL
  g <- if (is.null(gsea)) NULL else tryCatch(
    do.call(mmc_ma_gsea, c(list(ecs, pw), gsea)),
    error = function(e) {
      gsea_error <<- conditionMessage(e)
      NULL
    })
  list(prep = prep, input_mzlist = sig, matched = matched, ecs = ecs,
       pathways = pw, perms = perms, stats = stats, gsea = g,
       gsea_error = gsea_error)
}

#' Write one contrast's MetaboAnalyst-compatible results
#'
#' Writes, into `dir`:
#' \describe{
#'   \item{ma_pathway_results.tsv}{The pathway table in MetaboAnalyst's order
#'     (`ma_rank`), with the model's pathway id, every statistic, the Holm
#'     columns named `Holm_adj_FET` / `Holm_adj_EASE` / `Holm_adj_Gamma`, and
#'     the significant ECs.}
#'   \item{mummichog_pathway_enrichment_mummichog.csv}{The table MetaboAnalystR
#'     writes under that name, verbatim (`stats$restable`).}
#'   \item{ma_matched_features.tsv}{Every kept (feature, compound, ion form)
#'     match with its EC, feature id and polarity.}
#'   \item{ma_input.tsv}{The peak list the analysis was given.}
#'   \item{ma_run_info.json}{Engine, reference, parameters, model, Gamma fit
#'     and counts.}
#' }
#' and, when GSEA ran, `ma_gsea_results.tsv` (the GSEA table in
#' MetaboAnalyst's order; its settings and provenance go into the run record;
#' when GSEA failed, the record says so with its error). Only the engine's own
#' file names, including a previous failure record, are removed before writing,
#' so nothing else in `dir` is touched.
#'
#' @param run  Output of `mmc_ma_run_contrast()`.
#' @param peaks The peak list given to it.
#' @param dir  Destination directory (created if absent).
#' @param info Named list of run information to record (parameters, model).
#' @param pathway_ids Named character vector, pathway name -> model pathway id.
#' @return Character vector of the written paths.
#' @examples
#' \dontrun{
#' write_mummichog_ma_results(run, peaks, "out/mummichog_pinned/A_vs_B/metaboanalyst",
#'                            info = list(), pathway_ids = c())
#' }
write_mummichog_ma_results <- function(run, peaks, dir, info, pathway_ids) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(dir, c(.MMC_MA_RESULTS_FILE, .MMC_MA_UPSTREAM_CSV,
                            .MMC_MA_FEATURES_FILE, .MMC_MA_INPUT_FILE,
                            .MMC_MA_INFO_FILE))
  gsea_path <- file.path(dir, .MMC_MA_GSEA_FILE)
  .mmc_ma_clear_results(dir)

  st <- run$stats
  tbl <- st$table
  names(tbl)[names(tbl) == "adj_fet"]   <- "Holm_adj_FET"
  names(tbl)[names(tbl) == "adj_ease"]  <- "Holm_adj_EASE"
  names(tbl)[names(tbl) == "adj_gamma"] <- "Holm_adj_Gamma"
  pid <- if (is.null(pathway_ids)) rep(NA_character_, nrow(tbl)) else
    unname(pathway_ids[tbl$pathway])
  tbl <- data.frame(ma_rank = seq_len(nrow(tbl)), pathway_id = pid, tbl,
                    check.names = FALSE, stringsAsFactors = FALSE)
  readr::write_tsv(tbl, paths[1], na = "NA")
  utils::write.csv(st$restable, paths[2], row.names = TRUE)

  key <- run$prep$table
  m <- run$ecs$matched_res
  hit <- match(m$Query.Mass, key$mz)
  feats <- data.frame(feature_id = key$feature_id[hit],
                      polarity = ifelse(key$positive[hit], "positive", "negative"),
                      m, stringsAsFactors = FALSE, check.names = FALSE)
  readr::write_tsv(feats, paths[3], na = "NA")
  readr::write_tsv(peaks, paths[4], na = "NA")

  g <- st$gamma
  info <- c(info, list(
    engine = "metaboanalyst",
    reference = list(repo = "xia-lab/MetaboAnalystR",
                     commit = "398476ae2a0c996e390925fa62adc46b5ecde334"),
    fitdistrplus_version = as.character(utils::packageVersion("fitdistrplus")),
    rt_tol = run$prep$rt_tol,
    gamma = list(fitted = isTRUE(g$fitted),
                 shape = if (isTRUE(g$fitted)) unname(g$estimate["shape"]) else NULL,
                 scale = if (isTRUE(g$fitted)) unname(g$estimate["scale"]) else NULL,
                 error = g$error),
    counts = list(features_given = nrow(peaks),
                  features_kept = length(run$prep$ref_mzlist),
                  removed_out_of_range = unname(run$prep$removed[["out_of_range"]]),
                  removed_duplicate_mz = unname(run$prep$removed[["duplicate_mz"]]),
                  significant_features = length(run$input_mzlist),
                  empirical_compounds = length(run$ecs$total_matched_ecpds),
                  significant_ecs = length(run$ecs$input_ecpdlist),
                  pathways_tested = length(run$pathways$name),
                  pathways_reported = nrow(st$resmat))))
  gs <- run$gsea
  if (!is.null(gs)) {
    readr::write_tsv(gs$table, gsea_path, na = "NA")
    info$gsea <- list(ran = TRUE, metric_label = gs$metric_label,
                      params = gs$params, deviations = gs$deviations,
                      reference = gs$reference, adjustment = "BH (p.adjust 'fdr')",
                      n_ec = gs$n_ec, n_positions = gs$n_positions,
                      n_pathways = gs$n_pathways)
    paths <- c(paths, gsea_path)
  } else {
    info$gsea <- list(ran = FALSE, error = run$gsea_error)
  }
  jsonlite::write_json(info, file.path(dir, .MMC_MA_INFO_FILE), auto_unbox = TRUE,
                       pretty = TRUE, digits = NA, null = "null")
  paths
}

#' Read one contrast's MetaboAnalyst-compatible pathway table for the report
#'
#' Reads `ma_pathway_results.tsv` and adds the three columns the existing plot
#' and table builders key on — `overlap_size` (significant ECs in the pathway),
#' `pathway_size` (ECs in the pathway) and `p-value` (FET, MetaboAnalyst's own
#' plotted p-value) — without changing any value. The run information is
#' attached as `attr(, "ma_info")`.
#'
#' @param files The contrast's files (the result table is found by name).
#' @return A data frame in MetaboAnalyst's row order, or NULL when the contrast
#'   has no result table.
#' @examples
#' \dontrun{
#' read_mummichog_ma_pathways(list.files("out/mummichog_pinned/A_vs_B",
#'                                       recursive = TRUE, full.names = TRUE))
#' }
read_mummichog_ma_pathways <- function(files) {
  res <- files[basename(files) == .MMC_MA_RESULTS_FILE]
  if (length(res) == 0L) return(NULL)
  tbl <- as.data.frame(readr::read_tsv(res[[1]], show_col_types = FALSE,
                                       na = "NA"), check.names = FALSE)
  # A statistic that is NA throughout (e.g. Gamma when the fit failed) would be
  # read back as logical; every statistic column is numeric.
  num_cols <- c("ma_rank", "Pathway total", "Hits.total", "Hits.sig", "Expected",
                "FET", "EASE", "Total.EC", "Sig.EC", "Gamma", "Emp.Hits",
                "Empirical", "Holm_adj_FET", "Holm_adj_EASE", "Holm_adj_Gamma")
  for (col in intersect(num_cols, names(tbl))) {
    tbl[[col]] <- as.numeric(tbl[[col]])
  }
  tbl$pathway_id <- as.character(tbl$pathway_id)
  tbl$overlap_size <- tbl[["Sig.EC"]]
  tbl$pathway_size <- tbl[["Total.EC"]]
  tbl[["p-value"]] <- tbl[["FET"]]
  info_file <- files[basename(files) == .MMC_MA_INFO_FILE]
  info <- NULL
  if (length(info_file) > 0L) {
    info <- jsonlite::fromJSON(info_file[[1]], simplifyVector = TRUE)
    attr(tbl, "ma_info") <- info
  }
  gsea_file <- files[basename(files) == .MMC_MA_GSEA_FILE]
  if (length(gsea_file) > 0L) {
    g <- as.data.frame(readr::read_tsv(
      gsea_file[[1]], show_col_types = FALSE, na = "NA",
      col_types = readr::cols(
        "Pathway" = "c", "Leading-edge EmpiricalCompounds" = "c",
        "EmpiricalCompounds" = "c", "ES fallback reason" = "c",
        "ES defaulted by reference implementation" = "l", .default = "d")),
      check.names = FALSE)
    gi <- info$gsea %||% list()
    # The shape the report's GSEA tab and Summary read (as for python GSEA).
    attr(tbl, "ma_gsea") <- list(
      table = g, metric_label = gi$metric_label, params = gi$params,
      deviations = as.character(unlist(gi$deviations)),
      reference = gi$reference, nes_note = mmc_gsea_nes_note(),
      n_ec = gi$n_ec, n_positions = gi$n_positions, n_pathways = nrow(g),
      adjustment = gi$adjustment, engine = "metaboanalyst")
  }
  tbl
}

#' Is a pathway table the MetaboAnalyst-compatible engine's?
#'
#' @param pathways A pathway data frame.
#' @return TRUE when it carries the engine's Holm and rank columns.
#' @noRd
.mmc_is_ma_table <- function(pathways) {
  is.data.frame(pathways) &&
    all(c("ma_rank", "Holm_adj_FET", "FET", "EASE", "Gamma") %in% names(pathways))
}
