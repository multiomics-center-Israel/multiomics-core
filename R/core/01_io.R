#' Save a data frame as TSV (creating parent dir if needed)
#'
#' @param x Data frame.
#' @param dir Directory path.
#' @param filename Filename.
#' @return The full path.
save_tsv <- function(x, dir, filename) {
    ensure_dir(dir)
    path <- file.path(dir, filename)
    readr::write_tsv(x, path)
    path
}

#' Save a data frame as TSV to a full path
save_tsv_path <- function(x, path) {
    ensure_dir(dirname(path))
    readr::write_tsv(x, path)
    path
}

#' Read a table automatically detecting TSV vs CSV by extension
read_table_auto <- function(path) {
    ext <- tolower(tools::file_ext(path))
    df <- if (ext %in% c("tsv", "txt")) {
        readr::read_tsv(path, show_col_types = FALSE)
    } else {
        readr::read_csv(path, show_col_types = FALSE)
    }
    # Convert tibble to data.frame to support rownames and proper subsetting
    as.data.frame(df)
}
