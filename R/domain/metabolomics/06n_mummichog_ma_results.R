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
#' `"python"` (the pinned mummichog 2.7.0 engine, the default when unset) or
#' `"metaboanalyst"` (the native MetaboAnalyst-compatible engine). Every
#' consumer — the DAG, the stage, the report sections, the report text, the
#' polarity-label resolver — asks this rather than reading the key itself.
#'
#' @param config Full pipeline config.
#' @return `"python"` or `"metaboanalyst"`. Any other value stops.
#' @examples
#' mmc_engine(list(modes = list(metabolomics = list(enrichment = list(
#'   mummichog = list(engine = "metaboanalyst"))))))
mmc_engine <- function(config) {
  eng <- config$modes$metabolomics$enrichment$mummichog$engine %||% "python"
  eng <- tolower(trimws(as.character(eng)))
  if (length(eng) != 1L || is.na(eng) || !eng %in% c("python", "metaboanalyst")) {
    .mmc_stop("mummichog.engine must be \"python\" (the pinned mummichog 2.7.0 ",
              "engine, default) or \"metaboanalyst\" (the MetaboAnalyst-",
              "compatible engine); got '",
              paste(config$modes$metabolomics$enrichment$mummichog$engine,
                    collapse = ", "), "'.")
  }
  eng
}

# The engine's own result files, per contrast.
.MMC_MA_SUBDIR        <- "metaboanalyst"
.MMC_MA_RESULTS_FILE  <- "ma_pathway_results.tsv"
.MMC_MA_UPSTREAM_CSV  <- "mummichog_pathway_enrichment_mummichog.csv"
.MMC_MA_FEATURES_FILE <- "ma_matched_features.tsv"
.MMC_MA_INPUT_FILE    <- "ma_input.tsv"
.MMC_MA_INFO_FILE     <- "ma_run_info.json"
.MMC_MA_GSEA_FILE     <- "ma_gsea_results.tsv"

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
#' @return A list: `prep`, `input_mzlist`, `matched`, `ecs`, `pathways`,
#'   `perms`, `stats`, `gsea` (NULL when skipped or nothing passed its size
#'   filter).
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
  g <- if (is.null(gsea)) NULL else do.call(mmc_ma_gsea, c(list(ecs, pw), gsea))
  list(prep = prep, input_mzlist = sig, matched = matched, ecs = ecs,
       pathways = pw, perms = perms, stats = stats, gsea = g)
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
#' MetaboAnalyst's order; its settings and provenance go into the run record).
#' Only these six names are removed before writing, so nothing else in `dir`
#' is touched.
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
  stale <- c(paths, gsea_path)
  unlink(stale[file.exists(stale)])

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
    info$gsea <- list(ran = FALSE)
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
