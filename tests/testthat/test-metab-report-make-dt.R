# make_dt() in the metabolomics report: small non-p values must stay readable.
# With fixed 2 decimals, importances such as 0.006 and 0.0053 both showed as
# 0.01, so rows that differed looked tied. Values below are synthetic.

METAB_ROOT <- normalizePath(if (dir.exists("R")) "." else "../..",
                            mustWork = FALSE)
METAB_TEMPLATE <- file.path(METAB_ROOT, "R", "pipeline", "metabolomics",
                            "templates", "report_metabolomics.Rmd")

# Parse the first complete `make_dt <- function(...)` expression out of the Rmd.
load_metab_make_dt <- function(env) {
  lines <- readLines(METAB_TEMPLATE, warn = FALSE)
  start <- grep("^make_dt <- function", lines)[1]
  if (is.na(start)) stop("make_dt() not found in ", METAB_TEMPLATE)
  for (end in start:length(lines)) {
    ex <- tryCatch(parse(text = lines[start:end]), error = function(e) NULL)
    if (!is.null(ex) && length(ex) == 1L) {
      eval(ex[[1]], envir = env)
      return(invisible(env$make_dt))
    }
  }
  stop("could not parse make_dt() from ", METAB_TEMPLATE)
}

test_that("make_dt shows non-p numbers with 3 significant digits", {
  skip_if_not_installed("DT")
  env <- new.env()
  env$`%||%` <- function(x, y) if (is.null(x)) y else x
  load_metab_make_dt(env)

  df <- data.frame(feature_id = c("F1", "F2"),
                   importance = c(0.006, 0.0053),
                   P.Value    = c(1e-4, 0.2))
  # make_dt() returns htmltools::tagList(<datatable>).
  dt <- env$make_dt(df, caption = "test")[[1]]
  expect_s3_class(dt, "datatables")
  render <- vapply(dt$x$options$columnDefs,
                   function(def) as.character(def$render %||% ""), character(1))
  targets <- lapply(dt$x$options$columnDefs, function(def) unlist(def$targets))

  imp_i <- which(names(df) == "importance") - 1L
  on_imp <- vapply(targets, function(t) imp_i %in% t, logical(1))
  expect_true(any(grepl("formatSignif(data, 3,", render[on_imp], fixed = TRUE)))
  expect_false(any(grepl("formatRound(data, 2,", render[on_imp], fixed = TRUE)))
  # The stored values are untouched; only their display is formatted.
  expect_equal(dt$x$data$importance, c(0.006, 0.0053))
})
