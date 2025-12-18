# ============================================================
# multiomics-core — run execution info utilities
#
# Goal:
# - Snapshot the *exact* inputs needed to reproduce a run:
#   config, git commit, timestamp, and session info.
#
# Notes:
# - This function is safe to call both interactively and via {targets}.
# - It writes files under <run_dir>/run_metadata/.
# ============================================================

write_execution_info <- function(config, run_dir) {
  if (!is.list(config)) {
    stop("`config` must be a list (parsed YAML config).")
  }
  if (!is.character(run_dir) || length(run_dir) != 1L || !nzchar(run_dir)) {
    stop("`run_dir` must be a single, non-empty character path.")
  }
  if (!dir.exists(run_dir)) {
    stop(sprintf("`run_dir` does not exist: %s", run_dir))
  }
  
  info_dir <- file.path(run_dir, "execution_info")
  dir.create(info_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---- 1) config snapshot ----
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required for write_run_metadata(). Please add it to renv and tar_option_set().")
  }
  yaml::write_yaml(config, file.path(info_dir, "config_used.yaml"))
  
  # ---- 2) git commit ----
  git_commit <- tryCatch(
    system("git rev-parse HEAD", intern = TRUE),
    error = function(e) "UNKNOWN",
    warning = function(w) "UNKNOWN"
  )
  if (length(git_commit) < 1L) git_commit <- "UNKNOWN"
  writeLines(git_commit[1], file.path(info_dir, "git_commit.txt"))
  
  # ---- 3) timestamp ----
  writeLines(as.character(Sys.time()), file.path(info_dir, "timestamp.txt"))
  
  # ---- 4) session info ----
  capture.output(sessionInfo(), file = file.path(info_dir, "sessionInfo.txt"))
  
  # Return the directory path so {targets} can track this as a file target.
  return(info_dir)
}
