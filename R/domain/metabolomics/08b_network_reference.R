# R/domain/metabolomics/08b_network_reference.R
#
# Read + validate the pinned KEGG reaction-pair reference (built by
# multiomic-annotation-prep) for the metabolite-network stage. This layer knows
# nothing about KEGG or downloading: it is handed a local file path (already
# resolved and checksum-verified by resolve_pinned_asset()) and returns a
# validated data.frame of reaction-keyed compound pairs.
#
# The .tsv.gz carries a commented "# key: value" header that is the CANONICAL
# runtime metadata (schema_version, pair_definition_method, snapshot). The
# sidecar .manifest.json is release documentation and is NOT consumed here.
#
# Dependencies: readr (gz-aware TSV reader; already used elsewhere).


#' Parse the "# key: value" header block from a reference file
#'
#' @param path Path to the (optionally gzipped) reference file.
#' @return Named list of header values (character), one per "# key: value" line.
#' @noRd
.read_reference_header <- function(path) {
  lines <- readr::read_lines(path, n_max = 200L, progress = FALSE)
  hdr <- list()
  for (ln in lines) {
    if (!startsWith(ln, "#")) break
    kv <- sub("^#\\s*", "", ln)
    if (!grepl(":", kv)) next
    key <- trimws(sub(":.*$", "", kv))
    val <- trimws(sub("^[^:]*:\\s*", "", kv))
    if (nzchar(key)) hdr[[key]] <- val
  }
  hdr
}


#' Validate a loaded KEGG reaction-pair reference
#'
#' Enforces the schema-v1 contract. Every failure is a hard error with an
#' actionable message — there is deliberately no lenient/repair mode and no
#' fallback to live KEGG.
#'
#' @param df data.frame read from the reference file.
#' @param header Named list from `.read_reference_header()`.
#' @param expected_schema_version Integer schema this pipeline supports.
#' @param expected_method Pair-definition method this pipeline expects.
#' @return `invisible(TRUE)` on success; otherwise stops.
validate_kegg_reaction_pairs <- function(df, header,
                                         expected_schema_version = 1L,
                                         expected_method = "equation_side_cartesian_product") {
  # -- header contract (canonical runtime metadata) --------------------------
  sv <- suppressWarnings(as.integer(header$schema_version))
  if (is.null(header$schema_version) || is.na(sv)) {
    stop("KEGG reaction-pair reference: missing or non-integer 'schema_version' ",
         "header line.", call. = FALSE)
  }
  if (!identical(sv, as.integer(expected_schema_version))) {
    stop("KEGG reaction-pair reference: unsupported schema_version ", sv,
         " (this pipeline expects ", expected_schema_version,
         "). Repin a matching reference asset.", call. = FALSE)
  }
  if (is.null(header$pair_definition_method) ||
      !identical(header$pair_definition_method, expected_method)) {
    stop("KEGG reaction-pair reference: pair_definition_method '",
         header$pair_definition_method %||% "<missing>",
         "' does not match the expected '", expected_method, "'.", call. = FALSE)
  }

  # -- table contract --------------------------------------------------------
  required <- c("reaction_id", "substrate_id", "product_id", "equation",
                "equation_arrow")
  missing_cols <- setdiff(required, colnames(df))
  if (length(missing_cols) > 0) {
    stop("KEGG reaction-pair reference is missing required column(s): ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  if (nrow(df) == 0L) {
    stop("KEGG reaction-pair reference is empty (no data rows).", call. = FALSE)
  }
  if (anyNA(df$reaction_id) || anyNA(df$substrate_id) || anyNA(df$product_id)) {
    stop("KEGG reaction-pair reference has NA reaction/substrate/product ids.",
         call. = FALSE)
  }
  if (!all(grepl("^R[0-9]{5}$", df$reaction_id))) {
    bad <- unique(df$reaction_id[!grepl("^R[0-9]{5}$", df$reaction_id)])
    stop("KEGG reaction-pair reference has malformed reaction id(s): ",
         paste(utils::head(bad, 5L), collapse = ", "), ".", call. = FALSE)
  }
  if (!all(grepl("^C[0-9]{5}$", df$substrate_id) &
           grepl("^C[0-9]{5}$", df$product_id))) {
    stop("KEGG reaction-pair reference has malformed compound id(s) ",
         "(expected C#####).", call. = FALSE)
  }
  # Self-pairs must be removed + reported upstream by the builder. A valid
  # reference has none; reject rather than silently drop them again (corr. #6).
  if (any(df$substrate_id == df$product_id)) {
    bad <- unique(df$reaction_id[df$substrate_id == df$product_id])
    stop("KEGG reaction-pair reference contains self-pair(s) ",
         "(substrate == product), e.g. reaction ", utils::head(bad, 1L),
         ". A valid reference removes these; refusing to use it.", call. = FALSE)
  }
  key <- paste(df$reaction_id, df$substrate_id, df$product_id, sep = "\r")
  if (anyDuplicated(key) > 0L) {
    stop("KEGG reaction-pair reference has duplicate rows under ",
         "(reaction_id, substrate_id, product_id).", call. = FALSE)
  }
  invisible(TRUE)
}


#' Read + validate the pinned KEGG reaction-pair reference
#'
#' @param path Local path to the checksummed .tsv.gz (from
#'   `resolve_pinned_asset()`).
#' @param expected_schema_version Integer schema this pipeline supports.
#' @param expected_method Pair-definition method this pipeline expects.
#' @return A validated data.frame with columns `reaction_id`, `substrate_id`,
#'   `product_id`, `equation`, `equation_arrow`; the parsed header is attached as
#'   the "reference_header" attribute.
read_kegg_reaction_pairs <- function(path, expected_schema_version = 1L,
                                     expected_method = "equation_side_cartesian_product") {
  if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
    stop("KEGG reaction-pair reference file not found: '", path, "'.",
         call. = FALSE)
  }
  header <- .read_reference_header(path)
  df <- readr::read_tsv(
    path, comment = "#", progress = FALSE, show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  validate_kegg_reaction_pairs(
    df, header,
    expected_schema_version = expected_schema_version,
    expected_method = expected_method
  )
  attr(df, "reference_header") <- header
  df
}
