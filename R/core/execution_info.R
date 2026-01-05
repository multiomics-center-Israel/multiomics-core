# ============================================================
# multiomics-core — run execution info utilities
#
# Goal:
# Snapshot the exact inputs needed to reproduce a run:
# - config snapshot
# - git commit + timestamp + session info (existing code)
# - optionally copy _targets.R
# - optionally save project workspace (project.Rdata) using save.image()
#
# Notes:
# - safe both interactively and via {targets}
# - writes under <run_dir>/execution_info/
# - returns a character vector of file paths (targets-friendly)
# ============================================================

write_execution_info <- function(config,
                                 run_dir,
                                 config_path = NULL,
                                 manifest_df = NULL) {
  
  exec_dir <- file.path(run_dir, "execution_info")
  dir.create(exec_dir, recursive = TRUE, showWarnings = FALSE)
  
  written_files <- character(0)
  
  # --- config snapshot (full resolved list) ---
  cfg_path <- file.path(exec_dir, "config_used.yaml")
  yaml::write_yaml(config, cfg_path)
  written_files <- c(written_files, cfg_path)
  
  # --- record which config file was used (path) ---
  if (!is.null(config_path)) {
    cfgp_path <- file.path(exec_dir, "config_path.txt")
    writeLines(normalizePath(config_path, winslash = "/", mustWork = FALSE), cfgp_path)
    written_files <- c(written_files, cfgp_path)
  }
  
  # --- git commit ---
  git_path <- file.path(exec_dir, "git_commit.txt")
  git_commit <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) "UNKNOWN")
  writeLines(git_commit, git_path)
  written_files <- c(written_files, git_path)
  
  # --- session info ---
  sess_path <- file.path(exec_dir, "sessionInfo.txt")
  writeLines(capture.output(sessionInfo()), sess_path)
  written_files <- c(written_files, sess_path)
  
  # --- timestamp ---
  time_path <- file.path(exec_dir, "timestamp.txt")
  writeLines(as.character(Sys.time()), time_path)
  written_files <- c(written_files, time_path)
  
  # --- (A) copy renv.lock ---
  renv_lock_src <- "renv.lock"
  if (file.exists(renv_lock_src)) {
    renv_lock_dst <- file.path(exec_dir, "renv.lock")
    ok <- file.copy(renv_lock_src, renv_lock_dst, overwrite = TRUE)
    if (isTRUE(ok)) written_files <- c(written_files, renv_lock_dst)
  } else {
    warning("renv.lock not found in project root; run may be harder to reproduce.")
  }
  
  # --- (B) copy _targets.R ---
  save_targets <- isTRUE(config$execution_info$save_targets_file %||% TRUE)
  if (save_targets) {
    src_targets <- "_targets.R"
    if (!file.exists(src_targets)) stop("Cannot find _targets.R in project root.")
    dst_targets <- file.path(exec_dir, "_targets.R")
    ok <- file.copy(src_targets, dst_targets, overwrite = TRUE)
    if (!isTRUE(ok)) stop("Failed to copy _targets.R into execution_info.")
    written_files <- c(written_files, dst_targets)
  }

  
  # --- (C) save targets manifest (recommended to pass in from a target) ---
  if (!is.null(manifest_df)) {
    man_path <- file.path(exec_dir, "targets_manifest.csv")
    utils::write.csv(manifest_df, man_path, row.names = FALSE)
    written_files <- c(written_files, man_path)
  }
  
  # --- (D) optional workspace snapshot (project.Rdata) ---
  save_ws <- isTRUE(config$execution_info$save_workspace %||% FALSE)
  if (save_ws) {
    ws_name <- config$execution_info$workspace_filename %||% "project.Rdata"
    ws_path <- file.path(exec_dir, ws_name)
    save.image(file = ws_path)
    if (!file.exists(ws_path)) stop("save.image() did not create project workspace file.")
    written_files <- c(written_files, ws_path)
  }
  
  written_files
}

# Optional helper if you want an explicit function name to create project.Rdata
# (still uses save.image(), but keeps callsite clearer).
write_project_rdata <- function(run_dir, ..., filename = "project.Rdata") {
  exec_dir <- file.path(run_dir, "execution_info")
  dir.create(exec_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(exec_dir, filename)
  
  objects_list <- list(...)
  save(list = names(objects_list), envir = list2env(objects_list), file = out)
  
  if (!file.exists(out)) stop("save() did not create project.Rdata.")
  out
}
