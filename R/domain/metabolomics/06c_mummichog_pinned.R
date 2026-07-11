# R/domain/metabolomics/06c_mummichog_pinned.R
#
# Version-pinned mummichog v2 engine, run as an ISOLATED Python subprocess via
# processx (never imported into R). Companion to 06b_mummichog.R, which takes a
# different route (reticulate + KEGGREST organism models). The two coexist:
#   * 06b  — reticulate, builds custom KEGG organism models (bee/E. coli/human).
#   * 06c  — processx, pinned venv (mummichog==2.7.0), light deps only
#            (readr/processx/jsonlite), built-in metabolic models by default.
#
# Design principles:
#   * Do NOT import mummichog internals. Call `python -m mummichog.main` only.
#   * mummichog lives in its OWN pinned venv (see setup-mummichog-venv.sh); the
#     main R pipeline never depends on Python packages.
#   * Fail loudly and early on bad input and on missing/unexpected output.
#   * Return the full set of produced files (+ a manifest) so {targets} can hash
#     them with format = "file". v2 emits a whole result tree, not one file.
#   * v2 output is STOCHASTIC (permutation-based, no seed control). Do not expect
#     bit-identical reruns; {targets} only reruns on input change, so that's fine.
#
# Verified against mummichog==2.7.0 (get_user_data.py / reporting.py + a real
# run):
#   * Input (`-f`): tab-delimited, POSITIONAL over the first four columns
#     [m/z, retention_time, p-value, statistic]; an optional 5th column is
#     carried through as CompoundID_from_user. The first line is ALWAYS consumed
#     as a header and skipped, so a header row is required.
#   * m/z is hard-filtered to 50 < m/z < 2000 (MASS_RANGE); out-of-range features
#     are silently dropped by mummichog — we pre-filter and warn instead.
#   * Duplicate input lines are silently de-duplicated before row numbers are
#     assigned; writing feature_id as the 5th column keeps every line unique and
#     lets us join results back on the carried id rather than fragile row order.
#   * Output tree: <out_dir>/<timestamp>.<project>/ with result.html, a tables/
#     subdir (mcg_pathwayanalysis_*.tsv/.xlsx, mcg_modularanalysis_*.tsv/.xlsx,
#     ListOfEmpiricalCompounds.tsv, userInputData.txt,
#     userInput_to_EmpiricalCompounds.tsv), figures/ and js/. Tables are
#     .tsv/.xlsx — never .csv.
#
# Dependencies (R): processx, readr, jsonlite. (All light; no Bioconductor.)


# ==== internal helpers ======================================================

#' Signal a mummichog wrapper error with a consistent prefix
#' @param ... Parts of the message, pasted together.
#' @return Never returns; raises an error.
#' @noRd
.mmc_stop <- function(...) stop(paste0("[mummichog] ", ...), call. = FALSE)

#' Find the first matching column name in a data.frame
#'
#' Tries exact candidates first, then a case-insensitive regex fallback. Used to
#' locate messy m/z and retention-time columns (e.g. "m/z", "RT [min]"). Mirrors
#' the detection used in 06b_mummichog.R so both engines pick the same columns.
#'
#' @param df         A data.frame.
#' @param candidates Character vector of exact column names to try, in order.
#' @param regex      Fallback regex matched case-insensitively against names.
#' @return The matched column name, or NULL if none matched.
#' @noRd
.mmc_find_col <- function(df, candidates, regex = NULL) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0) return(hit[[1]])
  if (!is.null(regex)) {
    idx <- grep(regex, names(df), ignore.case = TRUE)
    if (length(idx) > 0) return(names(df)[idx[1]])
  }
  NULL
}

#' Validate a statistics table for mummichog input
#'
#' Coerces the requested columns to numeric and enforces mummichog's sanity
#' rules, returning a clean numeric frame in canonical order (row order
#' preserved, no rows dropped here — range filtering happens in the writer).
#'
#' @param stats_tbl  A data.frame / tibble.
#' @param mz_col,rt_col,p_col,stat_col Column names to read.
#' @param require_rt Whether retention time is required (TRUE for v2).
#' @return A data.frame with columns mz, rtime, pvalue, statistic.
#' @noRd
.mmc_validate_stats <- function(stats_tbl,
                                mz_col, rt_col, p_col, stat_col,
                                require_rt = TRUE) {
  if (!is.data.frame(stats_tbl)) {
    .mmc_stop("stats_tbl must be a data.frame / tibble.")
  }
  if (nrow(stats_tbl) == 0L) {
    .mmc_stop("stats_tbl has zero rows.")
  }

  required <- if (require_rt) c(mz_col, rt_col, p_col, stat_col) else c(mz_col, p_col, stat_col)
  missing <- setdiff(required, names(stats_tbl))
  if (length(missing) > 0) {
    .mmc_stop("Missing required columns: ", paste(missing, collapse = ", "),
              ". Present: ", paste(names(stats_tbl), collapse = ", "))
  }

  num <- function(x, nm) {
    v <- suppressWarnings(as.numeric(as.character(x)))
    if (any(is.na(v) & !is.na(x))) {
      .mmc_stop("Column '", nm, "' contains non-numeric values.")
    }
    v
  }

  mz   <- num(stats_tbl[[mz_col]], mz_col)
  pval <- num(stats_tbl[[p_col]], p_col)
  stat <- num(stats_tbl[[stat_col]], stat_col)
  rt   <- if (require_rt) num(stats_tbl[[rt_col]], rt_col) else rep(NA_real_, length(mz))

  if (any(!is.finite(mz)) || any(mz <= 0)) {
    .mmc_stop("Invalid m/z values (must be finite and > 0).")
  }
  if (require_rt && (any(!is.finite(rt)) || any(rt < 0))) {
    .mmc_stop("Invalid retention time values (must be finite and >= 0).")
  }
  if (any(!is.finite(pval)) || any(pval < 0 | pval > 1)) {
    .mmc_stop("Invalid p-values (must be finite and within [0, 1]).")
  }
  if (any(!is.finite(stat))) {
    .mmc_stop("Invalid statistic values (must be finite).")
  }

  data.frame(
    mz = mz, rtime = rt, pvalue = pval, statistic = stat,
    stringsAsFactors = FALSE
  )
}

# mummichog's MASS_RANGE: features outside are silently excluded (strict bounds).
.MMC_MZ_MIN <- 50
.MMC_MZ_MAX <- 2000


# ==== input writer ==========================================================

#' Write a mummichog v2 input table (tab-delimited)
#'
#' Emits the four positional columns mummichog v2 reads — m/z, retention time,
#' p-value, statistic — plus the feature id as a 5th column (CompoundID_from_user).
#' The 5th column serves two purposes verified against 2.7.0: it keeps every line
#' unique (so mummichog's silent line de-duplication can't collapse rows and shift
#' row numbers), and mummichog echoes it into userInputData.txt and
#' userInput_to_EmpiricalCompounds.tsv, giving a robust join key back to your
#' original features. A header row is written because v2 always skips line 1.
#'
#' Features with m/z outside (50, 2000) are pre-filtered with a warning, because
#' mummichog would otherwise drop them silently. Row order is preserved. An
#' id-map sidecar of exactly what was sent is written to paste0(path, ".idmap.tsv").
#'
#' @param stats_tbl A data.frame with the feature statistics.
#' @param path      Output path for the input .tsv.
#' @param mz_col,rt_col,p_col,stat_col Column names in stats_tbl.
#' @param id_col    Column holding feature ids; if NULL or absent, row numbers
#'                  are used.
#' @return Normalized path to the written input file.
write_mummichog_input <- function(stats_tbl, path,
                                   mz_col = "mz",
                                   rt_col = "retention_time",
                                   p_col  = "p_value",
                                   stat_col = "statistic",
                                   id_col = NULL) {
  clean <- .mmc_validate_stats(stats_tbl, mz_col, rt_col, p_col, stat_col,
                               require_rt = TRUE)

  feature_id <- if (!is.null(id_col) && id_col %in% names(stats_tbl)) {
    as.character(stats_tbl[[id_col]])
  } else {
    as.character(seq_len(nrow(clean)))
  }

  # Pre-filter to mummichog's accepted m/z window and report what we drop, so a
  # silent exclusion never surprises anyone comparing input vs output counts.
  in_range <- clean$mz > .MMC_MZ_MIN & clean$mz < .MMC_MZ_MAX
  if (any(!in_range)) {
    warning(sprintf(
      "[mummichog] Dropping %d/%d feature(s) with m/z outside (%g, %g); mummichog would exclude these silently.",
      sum(!in_range), length(in_range), .MMC_MZ_MIN, .MMC_MZ_MAX
    ), call. = FALSE)
  }
  clean      <- clean[in_range, , drop = FALSE]
  feature_id <- feature_id[in_range]
  if (nrow(clean) == 0L) {
    .mmc_stop("No features remain after the m/z (", .MMC_MZ_MIN, ", ", .MMC_MZ_MAX,
              ") filter — nothing to send to mummichog.")
  }

  out <- data.frame(
    `mz`         = clean$mz,
    `rtime`      = clean$rtime,
    `p-value`    = clean$pvalue,
    `t-score`    = clean$statistic,
    `CompoundID` = feature_id,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(out, path)
  .mmc_write_idmap(feature_id, paste0(path, ".idmap.tsv"))
  normalizePath(path, mustWork = TRUE)
}

#' Write the id-map sidecar (what we actually sent to mummichog)
#'
#' Provenance only: maps the written row order to feature_id. The authoritative
#' join back to features uses the id mummichog echoes into its own tables (see
#' join_features_to_results), not this positional map.
#'
#' @param feature_id Character vector of feature ids, in written order.
#' @param path       Output path for the id-map .tsv.
#' @return Normalized path to the written id-map.
#' @noRd
.mmc_write_idmap <- function(feature_id, path) {
  idmap <- data.frame(
    written_row = seq_along(feature_id),
    feature_id  = as.character(feature_id),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(idmap, path)
  normalizePath(path, mustWork = TRUE)
}


# ==== output discovery ======================================================

#' List every file (not directory) under a mummichog output directory
#' @param out_dir Directory to scan recursively.
#' @return Character vector of file paths.
list_mummichog_files <- function(out_dir) {
  files <- list.files(out_dir, recursive = TRUE, full.names = TRUE,
                      all.files = FALSE, no.. = TRUE)
  files[file.exists(files) & !dir.exists(files)]
}

#' Classify a mummichog output file by its basename
#'
#' Kinds reflect the real v2.7.0 tree: pathway/module analyses are .tsv and
#' .xlsx (never .csv), tables live under tables/, figures are .pdf.
#'
#' @param path A file path.
#' @return A single character label.
classify_mummichog_file <- function(path) {
  b <- basename(path)
  if (b == "result.html")                              return("v2_result_html")
  if (grepl("^mcg_pathwayanalysis.*\\.tsv$",  b))      return("pathway_tsv")
  if (grepl("^mcg_pathwayanalysis.*\\.xlsx$", b))      return("pathway_xlsx")
  if (grepl("^mcg_modularanalysis.*\\.tsv$",  b))      return("module_tsv")
  if (grepl("^mcg_modularanalysis.*\\.xlsx$", b))      return("module_xlsx")
  if (b == "ListOfEmpiricalCompounds.tsv")             return("empirical_compounds")
  if (b == "userInput_to_EmpiricalCompounds.tsv")      return("user_input_to_ec")
  if (b == "userInputData.txt")                        return("user_input_copy")
  if (b == "exported_Compounds.json")                  return("exported_compounds_json")
  if (grepl("\\.pdf$", b))                             return("figure_pdf")
  if (grepl("\\.js$",  b))                             return("js")
  if (b %in% c("mummichog.log", "runner.log"))         return("log")
  "other"
}

#' Write a manifest describing every produced file
#' @param files         Character vector of file paths.
#' @param manifest_file Output path for the manifest .tsv.
#' @return Normalized path to the manifest.
write_mummichog_manifest <- function(files, manifest_file) {
  np <- normalizePath(files, mustWork = TRUE)
  manifest <- data.frame(
    path      = np,
    file      = basename(np),
    directory = dirname(np),
    kind      = vapply(np, classify_mummichog_file, character(1)),
    bytes     = file.info(np)$size,
    stringsAsFactors = FALSE
  )
  readr::write_tsv(manifest, manifest_file)
  normalizePath(manifest_file, mustWork = TRUE)
}


# ==== runner ================================================================

#' Run mummichog v2 as an isolated subprocess
#'
#' Invokes `python -m mummichog.main` in the pinned venv with the working
#' directory set to out_dir (so v2's default workdir resolves there). out_dir is
#' WIPED and recreated on each run for a clean, hashable result tree.
#'
#' @param infile         Path to the input table (from write_mummichog_input).
#' @param out_dir        Directory for results (wiped and recreated).
#' @param project        Output identifier; simple charset only (validated).
#' @param python         Path to the venv python that has mummichog installed.
#' @param network        Metabolic model: "human_mfn" / "worm", or a path to a
#'                        custom model JSON (e.g. one produced by 06b) for a
#'                        like-for-like organism comparison.
#' @param mode           Ionization mode: "pos_default" / "positive" / "negative".
#' @param instrument_ppm Instrument accuracy in ppm (integer).
#' @param permutations   Number of permutations for the null distribution.
#' @param cutoff         Optional significance p-value cutoff; NULL = mummichog's
#'                        automated cutoff.
#' @param timeout        Seconds before the process is killed.
#' @param extra_args     Extra CLI args passed through verbatim.
#' @return Sorted character vector of all produced files plus the manifest.
run_mummichog <- function(infile, out_dir, project = "mummichog_run",
                          python = Sys.getenv("MUMMICHOG_PYTHON",
                                              "envs/mummichog/bin/python"),
                          network = "human_mfn",
                          mode = "pos_default",
                          instrument_ppm = 10,
                          permutations = 100,
                          cutoff = NULL,
                          timeout = 3600,
                          extra_args = character()) {
  if (!nzchar(python) || !file.exists(python)) {
    .mmc_stop("Python executable not found: '", python, "'. ",
              "Set MUMMICHOG_PYTHON or run setup-mummichog-venv.sh.")
  }
  # Resolve to an absolute path: processx changes to `wd` (out_dir) before exec,
  # so a relative command (e.g. the default envs/mummichog/bin/python) would be
  # looked up under out_dir and fail to start. Same reasoning for a file-backed
  # model passed via `network`/-n; built-in names (human_mfn, worm) are not files.
  python <- normalizePath(python, mustWork = TRUE)
  if (file.exists(network)) network <- normalizePath(network, mustWork = TRUE)
  if (!grepl("^[A-Za-z0-9._-]+$", project)) {
    .mmc_stop("Use a simple project name (letters, digits, dot, underscore, hyphen).")
  }
  infile <- normalizePath(infile, mustWork = TRUE)

  unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  log_file <- file.path(out_dir, "runner.log")

  # Short flags only: mummichog's getopt long options for --cutoff and
  # --force_primary_ion are declared WITHOUT an "=", so they would not accept a
  # value. The first `-m mummichog.main` is Python's module flag; the later
  # `-m <mode>` is mummichog's ionization mode (Python passes it through).
  args <- c("-m", "mummichog.main",
            "-f", infile,
            "-o", project,
            "-n", network,
            "-m", mode,
            "-u", as.character(instrument_ppm),
            "-p", as.character(permutations))
  if (!is.null(cutoff)) args <- c(args, "-c", as.character(cutoff))
  args <- c(args, extra_args)

  # Force a headless matplotlib backend: mummichog imports matplotlib.pyplot and
  # writes figures, and a non-framework venv python on macOS must not reach for a
  # GUI backend. "current" inherits the caller's environment.
  result <- processx::run(
    command = python, args = args, wd = out_dir,
    env = c("current", MPLBACKEND = "Agg"),
    echo = TRUE, error_on_status = FALSE, timeout = timeout
  )

  writeLines(
    c(paste("Version:", "v2 (mummichog==2.7.0)"),
      paste("Timeout (s):", timeout),
      paste("Working directory:", out_dir),
      paste("Command:", python, paste(args, collapse = " ")),
      paste("Exit status:", result$status),
      paste("Timed out:", isTRUE(result$timeout)),
      "", "STDOUT:", result$stdout, "", "STDERR:", result$stderr),
    con = log_file
  )

  if (isTRUE(result$timeout)) {
    .mmc_stop("mummichog timed out after ", timeout, "s. See: ", log_file)
  }
  if (result$status != 0) {
    .mmc_stop("mummichog exited with status ", result$status, ". See: ", log_file)
  }

  files <- unique(c(list_mummichog_files(out_dir), log_file))
  base  <- basename(files)
  has_v2 <- any(base == "result.html") ||
    any(grepl("^mcg_(pathway|modular)analysis.*\\.(tsv|xlsx)$", base))
  if (!has_v2) {
    .mmc_stop("mummichog finished but no expected v2 outputs were found. See: ", log_file)
  }

  manifest_file <- write_mummichog_manifest(
    files, file.path(out_dir, "mummichog_manifest.tsv")
  )
  sort(unique(c(files, manifest_file)))
}

#' Convenience alias matching the v2 target name
#' @inheritParams run_mummichog
#' @param ... Passed to run_mummichog().
#' @return Sorted character vector of produced files plus the manifest.
run_mummichog_v2 <- function(infile, out_dir, project = "mummichog_run", ...) {
  run_mummichog(infile, out_dir, project, ...)
}


# ==== result readers ========================================================

#' Read the mummichog pathway-analysis table
#' @param files Character vector of mummichog output files.
#' @return A tibble of pathway results.
read_mummichog_pathways <- function(files) {
  tsv <- files[grepl("^mcg_pathwayanalysis.*\\.tsv$", basename(files))]
  if (length(tsv) == 0) .mmc_stop("No pathway analysis .tsv found among mummichog outputs.")
  readr::read_tsv(tsv[[1]], show_col_types = FALSE)
}

#' Read the mummichog module-analysis table
#' @param files Character vector of mummichog output files.
#' @return A tibble of module results.
read_mummichog_modules <- function(files) {
  tsv <- files[grepl("^mcg_modularanalysis.*\\.tsv$", basename(files))]
  if (length(tsv) == 0) .mmc_stop("No modular analysis .tsv found among mummichog outputs.")
  readr::read_tsv(tsv[[1]], show_col_types = FALSE)
}

#' Read the echoed user input table (userInputData.txt)
#'
#' This is mummichog's authoritative record of the features it kept, carrying the
#' feature id we sent (as CompoundID_from_user) alongside its own row id.
#'
#' @param files Character vector of mummichog output files.
#' @return A tibble with all columns as character.
read_mummichog_user_input <- function(files) {
  f <- files[basename(files) == "userInputData.txt"]
  if (length(f) == 0) .mmc_stop("No userInputData.txt found among mummichog outputs.")
  readr::read_tsv(f[[1]], show_col_types = FALSE,
                  col_types = readr::cols(.default = readr::col_character()),
                  name_repair = "minimal")
}

#' Read the input-to-empirical-compound map (carries feature ids)
#'
#' userInput_to_EmpiricalCompounds.tsv links each empirical compound (EID, the
#' unit mummichog reports in pathway/module tables) to the input feature and the
#' feature id we sent (CompoundID_from_user).
#'
#' @param files Character vector of mummichog output files.
#' @return A tibble with columns EID, feature_id, compounds, compound_names.
read_mummichog_empirical_map <- function(files) {
  f <- files[basename(files) == "userInput_to_EmpiricalCompounds.tsv"]
  if (length(f) == 0) .mmc_stop("No userInput_to_EmpiricalCompounds.tsv found among mummichog outputs.")
  raw <- readr::read_tsv(f[[1]], show_col_types = FALSE,
                         col_types = readr::cols(.default = readr::col_character()),
                         name_repair = "unique_quiet")
  needed <- c("EID", "CompoundID_from_user")
  miss <- setdiff(needed, names(raw))
  if (length(miss) > 0) {
    .mmc_stop("userInput_to_EmpiricalCompounds.tsv missing column(s): ",
              paste(miss, collapse = ", "))
  }
  data.frame(
    EID           = raw[["EID"]],
    feature_id    = raw[["CompoundID_from_user"]],
    compounds     = if ("compounds" %in% names(raw)) raw[["compounds"]] else NA_character_,
    compound_names = if ("compound_names" %in% names(raw)) raw[["compound_names"]] else NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Attach original feature ids to a mummichog results table
#'
#' Reworked to key on the feature id mummichog carries through its own output
#' (via the 5th input column), NOT the post-de-duplication "rowN" strings. Each
#' results row lists the overlapping empirical compounds (EIDs); this expands
#' those to the feature ids we sent and adds them as a ";"-joined column.
#'
#' @param results A results table (e.g. from read_mummichog_pathways) containing
#'                a column of comma-separated empirical-compound ids.
#' @param files   Character vector of mummichog output files (for the id map).
#' @param ec_col  Name of the column holding the empirical-compound ids.
#' @param sep     Separator between ids within ec_col.
#' @return `results` with an added `feature_ids` column.
join_features_to_results <- function(results, files,
                                     ec_col = "overlap_EmpiricalCompounds (id)",
                                     sep = ",") {
  if (!ec_col %in% names(results)) {
    .mmc_stop("results has no column '", ec_col, "' to map empirical compounds from.")
  }
  emap <- read_mummichog_empirical_map(files)
  eid_to_features <- split(emap$feature_id, emap$EID)

  results$feature_ids <- vapply(as.character(results[[ec_col]]), function(cell) {
    if (is.na(cell) || !nzchar(cell)) return(NA_character_)
    eids <- trimws(strsplit(cell, sep, fixed = TRUE)[[1]])
    eids <- eids[nzchar(eids)]
    feats <- unique(unlist(eid_to_features[eids], use.names = FALSE))
    feats <- feats[!is.na(feats) & nzchar(feats)]
    if (length(feats) == 0) NA_character_ else paste(feats, collapse = ";")
  }, character(1), USE.NAMES = FALSE)

  results
}
