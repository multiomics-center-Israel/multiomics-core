ensure_dir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  invisible(dir)
}

p_tag <- function(config, default = "NA") {
  p <- config$modes$proteomics$de$p_cutoff
  if (is.null(p) || is.na(p)) return(default)
  format(p, trim = TRUE, scientific = FALSE)
}

save_tsv <- function(x, dir, filename) {
  ensure_dir(dir)
  path <- file.path(dir, filename)
  readr::write_tsv(x, path)
  path
}

save_tsv_path <- function(x, path) {
  ensure_dir(dirname(path))
  readr::write_tsv(x, path)
  path
}



