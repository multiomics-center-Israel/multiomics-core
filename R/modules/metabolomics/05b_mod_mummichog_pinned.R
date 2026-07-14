# R/modules/metabolomics/05b_mod_mummichog_pinned.R
#
# Target-ready wrapper around the pinned mummichog v2 engine (06c). It assembles
# the input from DE results joined to feature annotations — DE table, m/z and
# retention-time columns from row_data (RT passed through unchanged, in minutes),
# logFC as the statistic.
#
# Wired into 00_pipe_metabolomics.R as the metab_mummichog_pinned_files target
# (format = "file"), added only when
# modes.metabolomics.enrichment.mummichog.enabled is true. Runs mummichog once
# PER CONTRAST and returns the flat list of every file produced across all
# contrasts; the report reads each contrast's pathway table back from those files.


#' Run the pinned mummichog engine for a single contrast
#'
#' Assembles the mummichog input from one contrast's DE table joined to the
#' shared m/z / RT annotations and runs the v2 subprocess into `contrast_dir`.
#' Independent of other contrasts. Returns only the produced file paths — the
#' per-contrast pathway tables are read back downstream from the tracked files
#' (see `read_mummichog_pathways_by_contrast()`), keeping the side-effecting run
#' inside the `format = "file"` target so it self-heals on output deletion.
#'
#' @param de_table One contrast's DE table (feature_id, P.Value, logFC).
#' @param row_data Feature annotations (shared); m/z, RT, feature_id.
#' @param mz_col,rt_col Column names located once by the caller.
#' @param contrast_label Human-readable contrast name (for messages).
#' @param contrast_dir Per-contrast output dir (mummichog_pinned/<contrast>/).
#' @param network,mode,instrument_ppm,permutations,cutoff,force_primary_ion,python
#'   Run parameters (a single global model + params for every contrast).
#' @return Character vector of files this contrast produced (input, id-map,
#'   result tree, manifest).
#' @noRd
.mmc_run_one_contrast <- function(de_table, row_data, mz_col, rt_col,
                                  contrast_label, contrast_dir,
                                  network, mode, instrument_ppm, permutations,
                                  cutoff, force_primary_ion, python) {
  for (col in c("feature_id", "P.Value", "logFC")) {
    if (!col %in% names(de_table)) {
      .mmc_stop("contrast '", contrast_label,
                "': DE table missing required column '", col, "'.")
    }
  }
  merged <- merge(de_table, row_data[, c("feature_id", mz_col, rt_col)],
                  by = "feature_id")
  stats_tbl <- data.frame(
    feature_id     = as.character(merged$feature_id),
    mz             = as.numeric(merged[[mz_col]]),
    # Retention time passed through unchanged (minutes); mummichog only uses RT
    # for relative coelution grouping, so the unit must simply be consistent.
    retention_time = as.numeric(merged[[rt_col]]),
    p_value        = merged$P.Value,
    statistic      = merged$logFC,
    stringsAsFactors = FALSE
  )
  stats_tbl <- stats_tbl[stats::complete.cases(stats_tbl), , drop = FALSE]
  if (nrow(stats_tbl) == 0L) {
    .mmc_stop("contrast '", contrast_label,
              "': no complete features after joining DE with m/z / RT.")
  }

  input_file <- write_mummichog_input(
    stats_tbl = stats_tbl,
    path      = file.path(contrast_dir, "input.tsv"),
    mz_col = "mz", rt_col = "retention_time",
    p_col = "p_value", stat_col = "statistic", id_col = "feature_id"
  )
  result_files <- run_mummichog_v2(
    infile         = input_file,
    out_dir        = file.path(contrast_dir, "v2"),
    project        = "mcg_pinned",
    python         = python,
    network        = network,
    mode           = mode,
    instrument_ppm = instrument_ppm,
    permutations   = permutations,
    cutoff         = cutoff,
    force_primary_ion = force_primary_ion
  )
  sort(unique(c(input_file, paste0(input_file, ".idmap.tsv"), result_files)))
}

#' Run pinned mummichog v2 enrichment for metabolomics, per contrast
#'
#' Runs the isolated, version-pinned mummichog v2 subprocess (06c) once for each
#' DE contrast independently — each contrast's own p-values define its
#' significant set against all features, using a single global model + params.
#' Knobs are read from config$modes$metabolomics$enrichment$mummichog.
#'
#' @param pre     Preprocessing results; uses `pre$row_data` for m/z and RT.
#' @param de_res  DE results from mod_metabolomics_de(); each contrast in
#'                `de_tables` (falls back to summary_df, or a bare data.frame) is
#'                run independently. DE tables need columns feature_id, P.Value,
#'                logFC.
#' @param config  Full pipeline config.
#' @param out_dir Output directory for the metabolomics mode.
#' @param python  Path to the pinned venv python (defaults to $MUMMICHOG_PYTHON,
#'                else envs/mummichog/bin/python).
#' @return Character vector of every file produced across all contrasts
#'   (per-contrast inputs, id-maps, result trees, manifests), suitable for a
#'   `format = "file"` target; NULL when mummichog is disabled in config.
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
  # Optional mummichog -z (force primary ion). NULL when unset -> mummichog
  # 2.7.0's default (require a primary ion); set false to allow non-primary
  # adducts. Passed through to run_mummichog_v2().
  force_primary_ion <- mummi_cfg$force_primary_ion
  # Select the metabolic model (06d): a URL+sha256 `model_ref` (fetched, verified
  # and cached) wins, then a local `model_json` path, then the built-in human
  # model. Falling back to human on a non-human organism is refused there rather
  # than silently producing species-mismatched results.
  network <- mmc_select_model(
    mummi_cfg,
    organism  = config$modes$metabolomics$organism,
    cache_dir = "envs/mummichog-models"
  )

  # -- Locate m/z and RT columns (shared across contrasts) --------------------
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

  # -- Resolve the set of contrasts to run ------------------------------------
  # Each contrast runs mummichog INDEPENDENTLY: its own DE p-values define the
  # significant set (Lsig) against Lref = all features. Prefer the named
  # per-contrast de_tables; fall back to a single table for non-standard de_res.
  de_tables <- if (!is.null(de_res$de_tables) && length(de_res$de_tables) > 0) {
    de_res$de_tables
  } else if (!is.null(de_res$summary_df)) {
    list(contrast = de_res$summary_df)
  } else if (is.data.frame(de_res)) {
    list(contrast = de_res)
  } else {
    .mmc_stop("could not extract any DE table (expected de_res$de_tables / $summary_df / a data.frame).")
  }
  contrast_names <- names(de_tables)
  if (is.null(contrast_names) || any(!nzchar(contrast_names))) {
    contrast_names <- paste0("contrast_", seq_along(de_tables))
  }

  # -- Preflight the shared runner ONCE ---------------------------------------
  # The interpreter + pinned-version check is the same for every contrast, so a
  # broken/stale venv is a stage-level setup error, not a per-contrast no-result.
  # Fail loud here, BEFORE the loop, so it isn't masked by the per-contrast
  # tryCatch below (which is meant only for per-contrast data failures).
  python <- .mmc_preflight_runner(python)

  # -- Run mummichog per contrast, into its own subdir ------------------------
  # A single contrast failing (bad DE table, empty overlap, ...) yields a NULL
  # result for that contrast + a warning; the rest still run and the report
  # section hides gracefully for the failed one.
  mummi_dir <- file.path(out_dir, "mummichog_pinned")
  # Sanitise contrast names into subdir names, de-duplicating so two labels that
  # collapse to the same token (e.g. "A-B" and "A_B") get distinct dirs instead
  # of clobbering each other's on-disk results. sep = "_" keeps names in charset.
  contrast_dirs <- make.unique(gsub("[^A-Za-z0-9]+", "_", contrast_names),
                               sep = "_")
  files  <- vector("list", length(de_tables))
  errors <- vector("list", length(de_tables))
  for (i in seq_along(de_tables)) {
    contrast     <- contrast_names[[i]]
    contrast_dir <- file.path(mummi_dir, contrast_dirs[[i]])
    files[[i]] <- tryCatch(
      .mmc_run_one_contrast(
        de_table          = de_tables[[i]],
        row_data          = row_data,
        mz_col            = mz_col,
        rt_col            = rt_col,
        contrast_label    = contrast,
        contrast_dir      = contrast_dir,
        network           = network,
        mode              = mode_flag,
        instrument_ppm    = ppm,
        permutations      = n_perm,
        cutoff            = p_cutoff,
        force_primary_ion = force_primary_ion,
        python            = python
      ),
      error = function(e) {
        warning("mummichog (pinned): contrast '", contrast, "' failed: ",
                conditionMessage(e), call. = FALSE)
        errors[[i]] <<- conditionMessage(e)
        character(0)
      }
    )
  }
  # Flat, de-duplicated file list across all contrasts. This is the value of the
  # format = "file" target, so deleting any tracked output re-runs the whole
  # stage (regenerating every contrast). The report reads per-contrast pathway
  # tables back from these paths via read_mummichog_pathways_by_contrast().
  produced <- sort(unique(unlist(files, use.names = FALSE)))

  # If NOTHING was produced across every contrast, the failure is systemic — a
  # shared runner/config problem (e.g. an invalid knob validated inside the
  # engine), not one contrast's data — so fail loud rather than let the report
  # silently hide an enabled stage. Isolated per-contrast data failures are still
  # skipped as long as at least one contrast succeeds. (The python preflight
  # above already catches interpreter/version issues before any work.)
  if (length(produced) == 0L) {
    first <- Filter(Negate(is.null), errors)
    .mmc_stop("every contrast failed with no mummichog output — likely a shared ",
              "config/runner problem, not per-contrast data. First error: ",
              if (length(first)) first[[1]] else "unknown")
  }

  # Record the sanitised-dir -> original-contrast map so the report can label
  # tabs, subtitles and exports with the real contrast names — the flat file list
  # otherwise only carries the sanitised directory names.
  dir.create(mummi_dir, recursive = TRUE, showWarnings = FALSE)
  map_file <- file.path(mummi_dir, "contrasts.tsv")
  readr::write_tsv(
    data.frame(dir = contrast_dirs, contrast = contrast_names,
               stringsAsFactors = FALSE),
    map_file
  )
  sort(unique(c(produced, map_file)))
}
