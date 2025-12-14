write_proteomics_limma_outputs <- function(de_tables,
                                           out_dir = "limma",
                                           prefix  = "limma_",
                                           ext     = ".csv") {
  stopifnot(is.list(de_tables))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  out_files <- vapply(names(de_tables), function(cn) {
    df <- de_tables[[cn]]
    stopifnot(is.data.frame(df) || inherits(df, "matrix"))
    f <- file.path(out_dir, paste0(prefix, cn, ext))
    utils::write.csv(df, f, row.names = FALSE, quote = FALSE)
    f
  }, character(1))
  
  unname(out_files)
}
