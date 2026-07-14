# R/domain/metabolomics/06c_mummichog_pinned.R
#
# Version-pinned mummichog v2 engine, run as an ISOLATED Python subprocess via
# processx (never imported into R): pinned venv (mummichog==2.7.0), light deps
# only (readr/processx/jsonlite), built-in metabolic model by default. This is
# the canonical mummichog engine (it replaced an earlier reticulate/KEGGREST
# implementation). Wired into the metabolomics DAG via mod_mummichog_pinned()
# (05b), gated by modes.metabolomics.enrichment.mummichog.enabled.
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

#' Platform-appropriate default path to the pinned venv's Python
#'
#' `python -m venv` puts the interpreter under bin/ on Unix and Scripts/ on
#' Windows, so the documented default must branch on the platform.
#' @return Relative path to the venv Python for the current OS.
#' @noRd
.mmc_default_python <- function() {
  if (.Platform$OS.type == "windows") {
    "envs/mummichog/Scripts/python.exe"
  } else {
    "envs/mummichog/bin/python"
  }
}

#' Does a path exist, treating a symlink as present even if unresolved?
#'
#' `file.exists()` follows symlinks. We only want to confirm the configured
#' interpreter path is a real directory entry, without depending on where a
#' symlink points — see `.mmc_abs_keep_symlink()` for why we must not resolve it.
#'
#' @param path A file path.
#' @return TRUE if `path` exists as a file/dir or is a symlink (even a dangling one).
#' @noRd
.mmc_exists_nofollow <- function(path) {
  if (!nzchar(path)) return(FALSE)
  path <- path.expand(path)           # expand ~ / ~user; does NOT resolve symlinks
  link <- Sys.readlink(path)          # "" if not a link, target if a link, NA if absent
  file.exists(path) || (!is.na(link) && nzchar(link))
}

#' Make a path absolute WITHOUT resolving symlinks
#'
#' processx runs the engine with `wd = out_dir`, so a relative command (the
#' default `envs/mummichog/bin/python`) would be looked up under out_dir and
#' fail — the interpreter path must be absolute. But `normalizePath()` also
#' *canonicalises* it, following symlinks, and that BREAKS a venv: a venv's
#' `bin/python` is typically a symlink to the base interpreter, and running the
#' resolved target (e.g. /Library/Frameworks/.../python3.11) starts the BASE
#' environment with no venv site-packages — so mummichog isn't importable and
#' the 2.7.0 check fails with "found <none>". Invoke the venv symlink itself so
#' Python detects the venv (via pyvenv.cfg). Only prefix a relative path with
#' the working directory; leave the final (possibly symlinked) component alone.
#'
#' @param path A file path (may be relative, absolute, or home-relative).
#' @return An absolute path with any symlink in the final component preserved.
#' @noRd
.mmc_abs_keep_symlink <- function(path) {
  # Expand ~ / ~user first (path.expand does NOT resolve symlinks), so a
  # home-relative MUMMICHOG_PYTHON becomes absolute instead of <cwd>/~/....
  path <- path.expand(path)
  # Already absolute? POSIX "/...", Windows drive "C:...", or UNC "\\\\host".
  if (grepl("^(/|[A-Za-z]:|\\\\\\\\)", path)) return(path)
  file.path(getwd(), path)
}

#' Fail unless the given interpreter has the exact pinned mummichog version
#'
#' Guards reproducibility: if MUMMICHOG_PYTHON points at a stale or shared
#' environment with a different mummichog, results would carry misleading
#' "2.7.0" provenance. Query the actual version and abort on any mismatch.
#'
#' @param python   Path to the interpreter to probe.
#' @param expected Exact version string required.
#' @return Invisibly `expected` on success; aborts otherwise.
#' @noRd
.mmc_check_mummichog_version <- function(python, expected = "2.7.0") {
  probe <- processx::run(
    python, c("-c", "from mummichog.config import VERSION; print(VERSION)"),
    error_on_status = FALSE
  )
  got <- trimws(probe$stdout)
  if (probe$status != 0 || !identical(got, expected)) {
    .mmc_stop("expected mummichog ", expected, " at '", python, "' but found '",
              if (nzchar(got)) got else "<none>",
              "'. Rebuild the pinned venv with `make setup` (or `make mummichog-lock`) ",
              "or fix MUMMICHOG_PYTHON.")
  }
  invisible(expected)
}

#' Preflight the shared mummichog runner
#'
#' Validates the environment every run shares — the venv interpreter exists and
#' carries the pinned mummichog version — and returns its absolute path. Split
#' out so a caller that runs the engine many times (e.g. once per contrast) can
#' fail loud ONCE up front on a broken or stale venv, instead of letting a
#' per-run `tryCatch` downgrade that setup error into an empty result.
#'
#' @param python Path to the interpreter (may be relative or a venv symlink).
#' @return The absolute interpreter path (symlink preserved); aborts otherwise.
#' @noRd
.mmc_preflight_runner <- function(python) {
  if (!nzchar(python) || !.mmc_exists_nofollow(python)) {
    .mmc_stop("Python executable not found: '", python, "'. ",
              "Run `make setup` once per machine to build the pinned venv and print ",
              "the MUMMICHOG_PYTHON line for your .Renviron, or set MUMMICHOG_PYTHON yourself.")
  }
  # Make the interpreter path absolute (processx runs with wd = out_dir, so a
  # relative command would be looked up there) but do NOT canonicalise it:
  # normalizePath() follows symlinks, and a venv bin/python is usually a symlink
  # to the base interpreter — running the resolved target loses the venv (and its
  # mummichog). See .mmc_abs_keep_symlink().
  python <- .mmc_abs_keep_symlink(python)
  # Reproducibility guard: confirm this interpreter actually has mummichog 2.7.0
  # before we record 2.7.0 provenance and run it.
  .mmc_check_mummichog_version(python, "2.7.0")
  python
}

#' Find the first matching column name in a data.frame
#'
#' Tries exact candidates first, then a case-insensitive regex fallback. Used to
#' locate messy m/z and retention-time columns (e.g. "m/z", "RT [min]").
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

#' Stage the input file outside a directory that is about to be wiped
#'
#' run_mummichog() wipes out_dir at the start of each run. If the caller placed
#' the input inside out_dir (a natural layout, e.g.
#' write_mummichog_input(file.path(out_dir, "input.tsv"))), copy it to a tempfile
#' outside out_dir and return that path so the input survives the wipe and
#' mummichog can still read it; otherwise return the normalized path unchanged.
#'
#' @param infile  Path to the input table (must exist).
#' @param out_dir Directory that will be wiped (need not exist yet).
#' @return Absolute path to the input to actually read from.
#' @noRd
.mmc_stage_infile <- function(infile, out_dir) {
  # winslash = "/" so the prefix comparison is separator-consistent on Windows:
  # normalizePath() there returns "\" while .Platform$file.sep is "/", so mixing
  # them made startsWith() miss the inside-out_dir case (input then got wiped).
  infile <- normalizePath(infile, winslash = "/", mustWork = TRUE)
  out_dir_abs <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  inside <- infile == out_dir_abs || startsWith(infile, paste0(out_dir_abs, "/"))
  if (inside) {
    copy <- tempfile(fileext = ".tsv")
    if (!file.copy(infile, copy, overwrite = TRUE)) {
      .mmc_stop("failed to copy the input out of out_dir before wiping it.")
    }
    return(normalizePath(copy, winslash = "/", mustWork = TRUE))
  }
  infile
}


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

#' Coerce a config flag to a single logical (TRUE/FALSE), or fail loudly
#'
#' NULL passes through (unset). Accepts a real logical scalar, 0/1, and the usual
#' boolean spellings (true/false, t/f, yes/no, on/off) — YAML normally yields a
#' logical, but a QUOTED value or a typo arrives as a string. Anything
#' unrecognised is a hard error rather than a silent default, so a bad
#' `force_primary_ion` can never quietly relax mummichog's primary-ion filter.
#'
#' @param x   The raw config value.
#' @param key Config key name, for the error message.
#' @return TRUE/FALSE, or NULL when `x` is NULL.
#' @noRd
.mmc_as_flag <- function(x, key = "value") {
  if (is.null(x)) return(NULL)
  if (is.logical(x) && length(x) == 1L && !is.na(x)) return(x)
  v <- tolower(trimws(as.character(x)))
  if (length(v) == 1L && !is.na(v)) {
    if (v %in% c("true", "t", "yes", "y", "on", "1"))  return(TRUE)
    if (v %in% c("false", "f", "no", "n", "off", "0")) return(FALSE)
  }
  .mmc_stop("'", key, "' must be true or false; got '",
            paste(as.character(x), collapse = ", "), "'.")
}

#' Assemble the `python -m mummichog.main` argument vector
#'
#' Pure builder split out of run_mummichog() so the CLI flags are unit-testable
#' without a live subprocess. SHORT flags only: mummichog 2.7.0's getopt long
#' options for --cutoff and --force_primary_ion are declared without a trailing
#' "=", so they cannot take a value — the short forms (-c, -z) do.
#'
#' @param infile,project,network,mode Resolved run inputs.
#' @param instrument_ppm,permutations Run parameters.
#' @param cutoff Optional significance cutoff; NULL = mummichog's auto cutoff.
#' @param force_primary_ion Optional logical. NULL leaves mummichog 2.7.0's
#'   default (force_primary_ion = TRUE, i.e. require a primary ion). TRUE/FALSE
#'   emit `-z True` / `-z False`; `-z` takes a VALUE in 2.7.0 (getopt "z:"), so a
#'   bare -z is invalid, and FALSE is what allows non-primary adducts through.
#' @param extra_args Extra CLI args appended verbatim.
#' @return Character vector of arguments following the interpreter.
#' @noRd
.mmc_build_cli_args <- function(infile, project, network, mode,
                                instrument_ppm, permutations,
                                cutoff = NULL, force_primary_ion = NULL,
                                extra_args = character()) {
  # `-m mummichog.main` is Python's module flag; the later `-m <mode>` is
  # mummichog's ionization mode (Python passes it through).
  args <- c("-m", "mummichog.main",
            "-f", infile,
            "-o", project,
            "-n", network,
            "-m", mode,
            "-u", as.character(instrument_ppm),
            "-p", as.character(permutations))
  if (!is.null(cutoff)) args <- c(args, "-c", as.character(cutoff))
  # Only emit -z when explicitly set, so an unset config leaves mummichog's
  # default (require primary ion) exactly as-is. Coerce to a strict logical first
  # so a quoted/typo'd config value errors instead of silently relaxing to False.
  fpi <- .mmc_as_flag(force_primary_ion, "force_primary_ion")
  if (!is.null(fpi)) {
    args <- c(args, "-z", if (fpi) "True" else "False")
  }
  c(args, extra_args)
}

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
#'                        custom model JSON for an organism-specific analysis.
#' @param mode           Ionization mode: "pos_default" / "positive" / "negative".
#' @param instrument_ppm Instrument accuracy in ppm (integer).
#' @param permutations   Number of permutations for the null distribution.
#' @param cutoff         Optional significance p-value cutoff; NULL = mummichog's
#'                        automated cutoff.
#' @param force_primary_ion Optional logical passed to mummichog's `-z`. NULL
#'                        keeps 2.7.0's default (require a primary ion); FALSE
#'                        emits `-z False` to allow non-primary adducts.
#' @param timeout        Seconds before the process is killed.
#' @param extra_args     Extra CLI args passed through verbatim.
#' @return Sorted character vector of all produced files plus the manifest.
run_mummichog <- function(infile, out_dir, project = "mummichog_run",
                          python = Sys.getenv("MUMMICHOG_PYTHON",
                                              .mmc_default_python()),
                          network = "human_mfn",
                          mode = "pos_default",
                          instrument_ppm = 10,
                          permutations = 100,
                          cutoff = NULL,
                          force_primary_ion = NULL,
                          timeout = 3600,
                          extra_args = character()) {
  # Validate the shared runner (interpreter exists, pinned version) and take its
  # absolute path. Kept here too — defensive when run_mummichog() is called
  # directly — even though mod_mummichog_pinned() also preflights before its loop.
  python <- .mmc_preflight_runner(python)
  if (file.exists(network)) network <- normalizePath(network, mustWork = TRUE)
  if (!grepl("^[A-Za-z0-9._-]+$", project)) {
    .mmc_stop("Use a simple project name (letters, digits, dot, underscore, hyphen).")
  }
  # out_dir is wiped below; if the input lives inside it, stage a copy outside
  # first so unlink() cannot delete the input before mummichog reads it.
  infile <- .mmc_stage_infile(infile, out_dir)

  unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  log_file <- file.path(out_dir, "runner.log")

  args <- .mmc_build_cli_args(infile, project, network, mode, instrument_ppm,
                              permutations, cutoff, force_primary_ion, extra_args)

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

#' Read per-contrast mummichog pathway tables from a flat file list
#'
#' The per-contrast stage writes each contrast's outputs under
#' `mummichog_pinned/<contrast>/...` and tracks every file in one flat
#' `format = "file"` target. This regroups that flat list by contrast directory
#' — the path segment directly under `mummichog_pinned/` — and reads each
#' contrast's pathway table, so the report can render one section per contrast
#' without the stage having to hold the tables in memory.
#'
#' @param files Character vector of mummichog output files across all contrasts.
#' @return A named list keyed by the original contrast label (recovered from the
#'   stage's `contrasts.tsv` map; falls back to the sanitised subdir name when
#'   the map is absent), each a pathway tibble. Contrasts that produced no
#'   pathway table are omitted, so the report simply shows no section for them;
#'   an empty list results when `files` is empty or none live under
#'   `mummichog_pinned/`.
read_mummichog_pathways_by_contrast <- function(files) {
  if (length(files) == 0) return(list())
  # Tolerate Windows backslash separators (some paths come from normalizePath(),
  # which uses "\\" there): match on a "/"-normalised copy, read from originals.
  norm <- gsub("\\\\", "/", files)

  # Recover original contrast labels from the sanitised-dir -> contrast map the
  # stage writes, so report tabs/subtitles show the real DE contrast name rather
  # than the sanitised directory token. Absent map -> fall back to the dir name.
  label_of <- function(dir) dir
  is_map   <- grepl("/mummichog_pinned/contrasts\\.tsv$", norm)
  if (any(is_map)) {
    m <- tryCatch(readr::read_tsv(files[is_map][[1]], show_col_types = FALSE),
                  error = function(e) NULL)
    if (!is.null(m) && all(c("dir", "contrast") %in% names(m))) {
      lut <- stats::setNames(as.character(m$contrast), as.character(m$dir))
      label_of <- function(dir) if (!is.na(lut[dir])) unname(lut[dir]) else dir
    }
  }

  # Contrast dir = the path segment immediately under mummichog_pinned/.
  pat   <- ".*/mummichog_pinned/([^/]+)/.*"
  under <- grepl(pat, norm)
  if (!any(under)) return(list())
  keep     <- files[under]                 # original paths, for reading
  contrast <- sub(pat, "\\1", norm[under])

  out <- list()
  for (dir in unique(contrast)) {
    grp <- keep[contrast == dir]
    # A contrast with no pathway table -> read errors -> NULL; only assign
    # non-NULL results (an `out[[x]] <- NULL` would drop the key anyway).
    pw <- tryCatch(read_mummichog_pathways(grp), error = function(e) NULL)
    if (!is.null(pw)) out[[label_of(dir)]] <- pw
  }
  out
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
