# R/modules/metabolomics/05b_mod_mummichog_pinned.R
#
# Target-ready wrapper around the pinned mummichog v2 engine (06c). It assembles
# the input from DE results joined to feature annotations — DE table, m/z and
# retention-time columns from row_data (RT passed through unchanged, in minutes),
# logFC as the statistic.
#
# Wired into 00_pipe_metabolomics.R (analysis outputs) as the metab_mummichog_*
# targets, added only when modes.metabolomics.enrichment.mummichog.enabled is
# true. Returns a character vector of produced files for a format = "file" target.


#' Run pinned mummichog v2 enrichment for metabolomics
#'
#' Builds the mummichog input from DE results joined to feature annotations,
#' then runs the isolated, version-pinned mummichog v2 subprocess (06c). Knobs
#' are read from config$modes$metabolomics$enrichment$mummichog.
#'
#' @param pre     Preprocessing results; uses `pre$row_data` for m/z and RT.
#' @param de_res  DE results from mod_metabolomics_de(); uses de_tables[[1]]
#'                (falls back to summary_df, or a bare data.frame), with columns
#'                feature_id, P.Value, logFC.
#' @param config  Full pipeline config.
#' @param out_dir Output directory for the metabolomics mode.
#' @param python  Path to the pinned venv python (defaults to $MUMMICHOG_PYTHON,
#'                else envs/mummichog/bin/python).
#' @return Character vector of produced file paths (input, id-map, result tree,
#'   manifest), suitable for a format = "file" target.
mod_mummichog_pinned <- function(pre, de_res, config, out_dir,
                                 python = Sys.getenv("MUMMICHOG_PYTHON",
                                                     .mmc_default_python())) {
  mummi_cfg <- config$modes$metabolomics$enrichment$mummichog %||% list()

  # Enable gate: do nothing — no input written, no venv invoked — unless the
  # config explicitly enables mummichog. The DAG also only adds this target when
  # enabled; this guard additionally covers direct callers.
  if (!isTRUE(mummi_cfg$enabled)) {
    message("mummichog (pinned): disabled in config — skipping")
    return(NULL)
  }

  # Analysis knobs from config (sensible defaults when unset).
  p_cutoff <- mummi_cfg$p_cutoff       %||% 0.05
  n_perm   <- mummi_cfg$n_permutations %||% 100
  ppm      <- mummi_cfg$tolerance_ppm  %||% 10
  ion_mode <- mummi_cfg$ionization_mode %||% "pos_default"
  mode_flag <- switch(
    tolower(as.character(ion_mode)),
    "positive" = "positive",
    "negative" = "negative",
    "pos_default"
  )
  # An optional local model JSON (model_json) selects a custom metabolic network;
  # otherwise mummichog's built-in human model is used. (URL+sha256 model_ref is
  # a follow-up.)
  network <- mummi_cfg$model_json %||% "human_mfn"

  # The built-in model is human. Refuse to silently run human pathways on a
  # non-human organism — that would produce species-mismatched results with no
  # warning. Require a custom model (model_json now; model_ref in a follow-up) or
  # a human organism. Organism comes from the general metabolomics key or the
  # legacy mummichog.organisms.
  if (identical(network, "human_mfn")) {
    org <- config$modes$metabolomics$organism %||% mummi_cfg$organisms
    org <- tolower(trimws(paste(unlist(org), collapse = " ")))
    if (nzchar(org) && !grepl("human|homo sapiens|\\bhsa\\b", org)) {
      .mmc_stop(
        "configured organism ('", org, "') is non-human, but the pinned engine's ",
        "built-in model is human-only. Supply a custom model via 'model_json' ",
        "(URL+sha256 'model_ref' support is coming), or set the organism to human."
      )
    }
  }

  # -- Extract the DE table ---------------------------------------------------
  de_table <- if (!is.null(de_res$de_tables) && length(de_res$de_tables) > 0) {
    de_res$de_tables[[1]]
  } else if (!is.null(de_res$summary_df)) {
    de_res$summary_df
  } else if (is.data.frame(de_res)) {
    de_res
  } else {
    NULL
  }
  if (is.null(de_table)) {
    .mmc_stop("could not extract a DE table (expected de_res$de_tables / $summary_df / a data.frame).")
  }
  for (col in c("feature_id", "P.Value", "logFC")) {
    if (!col %in% names(de_table)) {
      .mmc_stop("DE table missing required column '", col, "'.")
    }
  }

  # -- Locate m/z and RT columns in the feature annotations -------------------
  row_data <- pre$row_data
  mz_col <- .mmc_find_col(row_data,
                          c("m/z", "mz", "m.z", "MZ", "Mass", "m.z."),
                          "^m[./]z")
  if (is.null(mz_col)) {
    .mmc_stop("No m/z column found in row_data. Columns: ",
              paste(names(row_data), collapse = ", "))
  }
  rt_col <- .mmc_find_col(row_data,
                          c("RT [min]", "rt", "RT", "retention_time",
                            "RetentionTime", "RT..min.", "RT.min."),
                          "^RT|retention")
  if (is.null(rt_col)) {
    .mmc_stop("No retention time column found in row_data. Columns: ",
              paste(names(row_data), collapse = ", "))
  }

  if (!"feature_id" %in% names(row_data)) {
    row_data$feature_id <- if ("Metabolite" %in% names(row_data)) {
      row_data$Metabolite
    } else {
      rownames(row_data)
    }
  }

  # -- Merge DE + annotations, build the canonical stats table ----------------
  merged <- merge(de_table, row_data[, c("feature_id", mz_col, rt_col)],
                  by = "feature_id")
  stats_tbl <- data.frame(
    feature_id     = as.character(merged$feature_id),
    mz             = as.numeric(merged[[mz_col]]),
    # Retention time passed through unchanged (minutes). mummichog only uses RT
    # for relative coelution grouping, so the unit must simply be consistent.
    retention_time = as.numeric(merged[[rt_col]]),
    p_value        = merged$P.Value,
    statistic      = merged$logFC,
    stringsAsFactors = FALSE
  )
  stats_tbl <- stats_tbl[stats::complete.cases(stats_tbl), , drop = FALSE]
  if (nrow(stats_tbl) == 0L) {
    .mmc_stop("no complete features after joining DE results with m/z / RT annotations.")
  }

  # -- Write input, run the pinned v2 engine ----------------------------------
  mummi_dir  <- file.path(out_dir, "mummichog_pinned")
  input_file <- write_mummichog_input(
    stats_tbl = stats_tbl,
    path      = file.path(mummi_dir, "input.tsv"),
    mz_col    = "mz",
    rt_col    = "retention_time",
    p_col     = "p_value",
    stat_col  = "statistic",
    id_col    = "feature_id"
  )

  result_files <- run_mummichog_v2(
    infile         = input_file,
    out_dir        = file.path(mummi_dir, "v2"),
    project        = "mcg_pinned",
    python         = python,
    network        = network,
    mode           = mode_flag,
    instrument_ppm = ppm,
    permutations   = n_perm,
    cutoff         = p_cutoff
  )

  sort(unique(c(input_file, paste0(input_file, ".idmap.tsv"), result_files)))
}
