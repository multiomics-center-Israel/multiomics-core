ensure_dir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  invisible(dir)
}

write_tsv <- function(x, dir, filename) {
  ensure_dir(dir)
  path <- file.path(dir, filename)
  readr::write_tsv(x, path)
  path
}

p_tag <- function(config, default = "NA") {
  p <- config$modes$proteomics$limma$p_cutoff
  if (is.null(p) || is.na(p)) return(default)
  format(p, trim = TRUE, scientific = FALSE)
}
write_tsv <- function(x, dir, filename) {
  ensure_dir(dir)
  path <- file.path(dir, filename)
  readr::write_tsv(x, path)
  path
}

# Helper for safe saving (returns path)
do_write <- function(dat, out_dir, fname) {
  fpath <- file.path(out_dir, fname)
  # Using readr if available, otherwise base write.table could be used
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_tsv(dat, fpath)
  } else {
    utils::write.table(dat, fpath, sep="\t", quote=FALSE, row.names=FALSE)
  }
  fpath
}