# R/modules/metabolomics/05d_mod_mummichog_ma.R
#
# Target-ready stage for the MetaboAnalyst-compatible mummichog engine
# (06j-06n), selected with modes.metabolomics.enrichment.mummichog.engine:
# metaboanalyst. It fills the same target as the python stage (05b),
# metab_mummichog_pinned_files, and writes beside it — per contrast under
# mummichog_pinned/<contrast>/metaboanalyst/, with the same contrasts.tsv map —
# so the report, exports and Shiny payload read it through the same targets.
#
# Differences from the python stage that a reader of the results must know:
#   * ionization_mode: mixed is ONE joint analysis of both polarities (one
#     peak list, one background, one permutation null), as MetaboAnalyst does —
#     not one run per polarity. positive / negative read every feature as that
#     polarity, with MetaboAnalyst's ion lists (25 positive, 16 negative forms).
#   * the cutoff is inclusive (p <= p_cutoff), MetaboAnalyst's SetMummichogPval;
#   * permutation i is seeded 123 + i (MetaboAnalyst's own seeds), not from
#     params.seed;
#   * the model must be a readable JSON (model_ref / model_json): the built-in
#     human_mfn model lives inside the python package and has no masses here;
#   * GSEA runs here, in the stage, on the same ECs and pathway sets as the
#     ORA (06o), not in the report layer from python EC tables.


#' Run the MetaboAnalyst-compatible mummichog engine for every contrast
#'
#' Builds the library once from the configured model, then for each DE contrast
#' assembles the peak list (DE p-value and logFC joined to the m/z, retention
#' time and — in mixed mode — polarity of every feature), runs
#' `mmc_ma_run_contrast()` and writes its results. A contrast that fails warns,
#' and its previous result files are replaced by a failure record
#' (`ma_failed.json`); if every contrast fails the stage stops. A GSEA failure
#' only warns: the contrast's ORA is still written.
#'
#' @param pre     Preprocessing results; uses `pre$row_data`.
#' @param de_res  DE results; each table in `de_res$de_tables` (or
#'   `summary_df`, or a bare data frame) needs `feature_id`, `P.Value`, `logFC`.
#' @param config  Full pipeline config.
#' @param out_dir Metabolomics output directory.
#' @return Character vector of every file written, for a `format = "file"`
#'   target; NULL when mummichog is disabled.
#' @examples
#' \dontrun{
#' files <- mod_mummichog_metaboanalyst(pre, de_res, config, "out/metabolomics")
#' }
mod_mummichog_metaboanalyst <- function(pre, de_res, config, out_dir) {
  mummi_cfg <- config$modes$metabolomics$enrichment$mummichog %||% list()
  if (!isTRUE(mummi_cfg$enabled)) {
    message("mummichog (metaboanalyst): disabled in config — skipping")
    return(NULL)
  }
  if (!identical(mmc_engine(config), "metaboanalyst")) {
    .mmc_stop("mod_mummichog_metaboanalyst() called with mummichog.engine '",
              mmc_engine(config), "'.")
  }
  # Said once per stage run, ahead of the prerequisite checks below: a config
  # written for the python engine that names no engine lands here, and should
  # read why before it reads "needs a model".
  if (identical(mmc_engine_record(config)$engine_source, "default")) {
    message("mummichog: the config names no engine, so the default applies: ",
            "engine: metaboanalyst (MetaboAnalyst-compatible). To run the ",
            "pinned python mummichog instead, or to reproduce an earlier ",
            "analysis that used it, set engine: python. See docs/mummichog.md.")
  }

  # -- parameters ---------------------------------------------------------------
  mode <- tolower(as.character(mummi_cfg$ionization_mode %||% ""))
  if (length(mode) != 1L || !mode %in% c("mixed", "positive", "negative")) {
    .mmc_stop("engine 'metaboanalyst' needs ionization_mode \"mixed\" (one joint ",
              "analysis of both polarities), \"positive\" or \"negative\"; got '",
              paste(mummi_cfg$ionization_mode, collapse = ", "), "'. ",
              "(\"pos_default\" is a python mummichog ion list with no ",
              "MetaboAnalyst equivalent.)")
  }
  num1 <- function(x, key, default) {
    v <- x %||% default
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v <= 0) {
      .mmc_stop("mummichog.", key, " must be a single positive number; got '",
                paste(format(v), collapse = ", "), "'.")
    }
    v
  }
  p_cutoff <- num1(mummi_cfg$p_cutoff, "p_cutoff", 0.05)
  ppm      <- num1(mummi_cfg$tolerance_ppm, "tolerance_ppm", 10)
  n_perm   <- num1(mummi_cfg$n_permutations, "n_permutations", 100)
  rt_frac  <- num1(mummi_cfg$rt_frac, "rt_frac", 0.02)
  min_lib  <- num1(mummi_cfg$min_pathway_size, "min_pathway_size", 3)
  force_primary_ion <- .mmc_as_flag(mummi_cfg$force_primary_ion,
                                    "force_primary_ion") %||% TRUE
  # GSEA: same switch and overrides as the python engine's; defaults are
  # MetaboAnalyst's. It needs fgsea.
  gsea_on <- isTRUE(mummi_cfg$gsea_enabled %||% TRUE)
  if (gsea_on && !requireNamespace("fgsea", quietly = TRUE)) {
    message("mummichog (metaboanalyst): fgsea not available — GSEA skipped")
    gsea_on <- FALSE
  }
  gsea_settings <- list(
    n_perm   = mummi_cfg$gsea_permutations %||% .MMC_GSEA_DEFAULTS$n_perm,
    min_size = mummi_cfg$gsea_min_size     %||% .MMC_GSEA_DEFAULTS$min_size,
    max_size = mummi_cfg$gsea_max_size     %||% .MMC_GSEA_DEFAULTS$max_size,
    seed     = mummi_cfg$gsea_seed         %||% .MMC_GSEA_DEFAULTS$seed)

  # -- model -> library, once ---------------------------------------------------
  network <- mmc_select_model(mummi_cfg,
                              organism  = config$modes$metabolomics$organism,
                              cache_dir = "envs/mummichog-models")
  if (identical(network, "human_mfn") || !file.exists(network)) {
    .mmc_stop("engine 'metaboanalyst' needs a readable metabolic model JSON ",
              "(model_ref or model_json); the built-in human_mfn model lives ",
              "inside the python mummichog package and carries no compound ",
              "masses here.")
  }
  lib <- mmc_ma_library_from_model(network)
  pathway_ids <- stats::setNames(lib$pathways$id, lib$pathways$name)

  # -- features -----------------------------------------------------------------
  row_data <- pre$row_data
  mz_col <- .mmc_find_col(row_data, c("m/z", "mz", "m.z", "MZ", "Mass", "m.z."),
                          "^m[./]z")
  rt_col <- .mmc_find_col(row_data, c("RT [min]", "rt", "RT", "retention_time",
                                      "RetentionTime", "RT..min.", "RT.min."),
                          "^RT|retention")
  if (is.null(mz_col) || is.null(rt_col)) {
    .mmc_stop("engine 'metaboanalyst' needs m/z and retention-time columns in ",
              "the feature annotations. Columns: ",
              paste(names(row_data), collapse = ", "))
  }
  if (!"feature_id" %in% names(row_data)) {
    row_data$feature_id <- mmc_feature_ids(row_data)
  }
  row_data$.polarity <- if (mode == "mixed") {
    mmc_feature_polarity(row_data, mummi_cfg$polarity_column)
  } else {
    mode
  }

  de_tables <- if (!is.null(de_res$de_tables) && length(de_res$de_tables) > 0) {
    de_res$de_tables
  } else if (!is.null(de_res$summary_df)) {
    list(contrast = de_res$summary_df)
  } else if (is.data.frame(de_res)) {
    list(contrast = de_res)
  } else {
    .mmc_stop("could not extract any DE table (expected de_res$de_tables / ",
              "$summary_df / a data.frame).")
  }
  labels <- names(de_tables)
  if (is.null(labels) || any(!nzchar(labels))) {
    labels <- paste0("contrast_", seq_along(de_tables))
  }

  info <- list(params = list(
    ionization_mode = mode, polarity_column = mummi_cfg$polarity_column,
    tolerance_ppm = ppm, p_cutoff = p_cutoff, cutoff_rule = "p <= p_cutoff",
    n_permutations = n_perm, permutation_seeds = "123 + i",
    force_primary_ion = force_primary_ion, rt_frac = rt_frac,
    min_pathway_size = min_lib),
    model = list(id = lib$model$id, version = lib$model$version,
                 sha256 = .mmc_sha256_file(network)))

  # -- per contrast ---------------------------------------------------------------
  mummi_dir <- file.path(out_dir, "mummichog_pinned")
  dirs <- make.unique(gsub("[^A-Za-z0-9]+", "_", labels), sep = "_")
  written <- vector("list", length(de_tables))
  errors  <- vector("list", length(de_tables))
  failed  <- logical(length(de_tables))
  for (i in seq_along(de_tables)) {
    contrast <- labels[[i]]
    written[[i]] <- tryCatch({
      de <- de_tables[[i]]
      for (col in c("feature_id", "P.Value", "logFC")) {
        if (!col %in% names(de)) {
          .mmc_stop("contrast '", contrast, "': DE table missing column '",
                    col, "'.")
        }
      }
      merged <- merge(de, row_data[, c("feature_id", mz_col, rt_col, ".polarity")],
                      by = "feature_id")
      # The peak list's score is MetaboAnalyst's t.score: the signed DE
      # statistic the GSEA ranks on, chosen from the DE method exactly as for
      # the python engine's GSEA (mmc_gsea_rank_metric). The ORA does not use
      # it; it becomes each EC's GSEA score (the max over its members).
      rank <- mmc_gsea_rank_metric(
        merged, de_method = if (is.data.frame(de_res)) NULL else de_res$method)
      if (is.null(rank)) {
        .mmc_stop("contrast '", contrast, "': no usable signed statistic or ",
                  "logFC in the DE table to score features.")
      }
      peaks <- data.frame(feature_id = as.character(merged$feature_id),
                          mz = as.numeric(merged[[mz_col]]),
                          rt = as.numeric(merged[[rt_col]]),
                          p_value = as.numeric(merged$P.Value),
                          statistic = unname(rank$values[as.character(merged$feature_id)]),
                          polarity = merged$.polarity,
                          stringsAsFactors = FALSE)
      # MetaboAnalyst drops rows with a missing value before anything else
      # (SanityCheckMummichogData); the python stage does too. Done here,
      # visibly, because mmc_ma_prepare_peaks() itself refuses gaps.
      complete <- stats::complete.cases(peaks)
      if (!all(complete)) {
        message("mummichog (metaboanalyst): contrast '", contrast, "': ",
                sum(!complete), " feature(s) without a p-value, statistic, m/z ",
                "or retention time left out, as MetaboAnalyst does")
        peaks <- peaks[complete, , drop = FALSE]
      }
      # Feature-table order, as a table uploaded to MetaboAnalyst would arrive.
      # It matters: MetaboAnalyst's order() keeps tied p-values in input order,
      # and its permutations sample the list in that order. merge() re-sorts.
      peaks <- peaks[order(match(peaks$feature_id, row_data$feature_id)), ,
                     drop = FALSE]
      rownames(peaks) <- NULL
      run <- mmc_ma_run_contrast(peaks, lib, mode, ppm = ppm, cutoff = p_cutoff,
                                 n_perm = n_perm,
                                 force_primary_ion = force_primary_ion,
                                 rt_frac = rt_frac, min_lib = min_lib,
                                 gsea = if (gsea_on) c(gsea_settings,
                                   list(metric_label = rank$label)) else NULL)
      if (!is.null(run$gsea_error)) {
        warning("mummichog (metaboanalyst): contrast '", contrast, "': GSEA ",
                "failed, ORA kept: ", run$gsea_error, call. = FALSE)
      }
      write_mummichog_ma_results(
        run, peaks, file.path(mummi_dir, dirs[[i]], .MMC_MA_SUBDIR),
        info = c(list(contrast = contrast), info), pathway_ids = pathway_ids)
    }, error = function(e) {
      warning("mummichog (metaboanalyst): contrast '", contrast, "' failed: ",
              conditionMessage(e), call. = FALSE)
      errors[[i]] <<- conditionMessage(e)
      failed[[i]] <<- TRUE
      character(0)
    })
  }
  # A failed contrast's earlier results must not survive a rerun as if current:
  # its result files are replaced by a failure record, which stays tracked so
  # the export cleanup downstream still sees the contrast.
  for (i in which(failed)) {
    written[[i]] <- write_mummichog_ma_failure(
      file.path(mummi_dir, dirs[[i]], .MMC_MA_SUBDIR), labels[[i]], errors[[i]])
  }
  produced <- unlist(written, use.names = FALSE)
  if (all(failed)) {
    # Stopping here means the downstream sections and exports targets never
    # run, so their cleanup cannot catch this case: clear every contrast's
    # earlier presentation exports now. The set of export slugs a full run
    # writes is the set of contrast directory names (same sanitise and
    # make.unique over every label), and every one of them is stale here.
    old <- unlist(lapply(dirs, function(d) .mmc_contrast_export_paths(
      mummi_dir, d, all = TRUE)), use.names = FALSE)
    unlink(old[file.exists(old)])
    first <- Filter(Negate(is.null), errors)
    .mmc_stop("every contrast failed in the MetaboAnalyst-compatible engine. ",
              "First error: ", if (length(first)) first[[1]] else "unknown")
  }

  map_file <- file.path(mummi_dir, "contrasts.tsv")
  readr::write_tsv(data.frame(dir = dirs, contrast = labels,
                              stringsAsFactors = FALSE), map_file)
  sort(unique(c(produced, map_file)))
}
