#' Save a plotting expression to PNG (for base plots / pheatmap)
#'
#' @param out_file path to png file (if NULL -> draw to screen)
#' @param expr plotting expression (e.g. pheatmap::pheatmap(...))
#' @param width,height in pixels
#' @param res resolution (dpi)
with_png <- function(out_file, expr, width = 1600, height = 1200, res = 150) {
  # NOTE:
  # We must close the PNG device by device ID, not by dev.off(),
  # because some plotting functions (e.g. pheatmap) may change
  # the active graphics device internally.
  
  # If no output file is provided, just evaluate the expression
  # without opening any graphics device
  if (is.null(out_file)) {
    force(expr)
    return(invisible(NULL))
  }
  
  # Ensure the output file has a .png extension
  if (!grepl("\\.png$", out_file, ignore.case = TRUE)) {
    out_file <- paste0(out_file, ".png")
  }
  
  # Create output directory if it does not exist
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  
  # Open a PNG graphics device
  grDevices::png(
    filename = out_file,
    width    = width,
    height   = height,
    res      = res
  )
  
  # Store the device ID of the PNG device we just opened
  # (important: plotting functions may change the active device internally)
  png_dev <- grDevices::dev.cur()
  
  # Ensure that we close exactly this PNG device on exit,
  # even if the active device changes inside `expr`
  on.exit({
    if (png_dev %in% grDevices::dev.list()) {
      try(grDevices::dev.off(which = png_dev), silent = TRUE)
    }
  }, add = TRUE)
  
  # Evaluate the plotting expression
  force(expr)
  
  invisible(out_file)
  
}


