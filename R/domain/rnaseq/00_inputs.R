#' Load rna inputs
#'
#' @param config list as returned by load_config()
#' @return list (counts, sample_map, meta, contrasts, engine, ...)
load_rna_inputs <- function(config) {
    load_omics_inputs(config, mode = "rna")
}

#' Validate RNA inputs
#'
#' @param inputs List returned by \code{load_rna_inputs()}
#' @param cfg    The \code{config$modes$rna} sub-list
validate_rna_inputs <- function(inputs, cfg) {
    gene_id_col <- cfg$id_columns$gene_id
    sample_col <- cfg$id_columns$sample_col %||% "SampleID"

    if (!gene_id_col %in% names(inputs$counts)) {
        stop("RNA counts missing gene id column: ", gene_id_col)
    }

    meta <- inputs$metadata
    if (!sample_col %in% names(meta)) {
        stop("metadata missing sample column: ", sample_col)
    }

    # Basic checks for counts (numeric columns)
    sample_cols <- setdiff(names(inputs$counts), gene_id_col)
    if (length(sample_cols) == 0) stop("Counts table has no sample columns.")

    # Check alignment consistency
    meta_samples <- unique(as.character(meta[[sample_col]]))
    count_samples <- sample_cols

    missing_in_counts <- setdiff(meta_samples, count_samples)
    if (length(missing_in_counts) > 0) {
        stop("metadata contains samples not present in counts: ", paste(missing_in_counts, collapse = ", "))
    }

    invisible(TRUE)
}

#' Generic loader helper (likely available in domain/proteomics/00_inputs too, or should be in core if strictly dry)
#' Keeping local copy or assuming sourced if not in core.
# For safety in this refactor step, I will rely on the fact that we might have moved load_omics_inputs to core?
# No, I put it in R/domain/proteomics/00_inputs.R. I should properly put it in R/core/01_io.R or duplicating it.
# Duplicating small helpers avoids cross-domain fragility if core isn't strict.
load_omics_inputs <- function(config, mode = c("proteomics", "rna")) {
    mode <- match.arg(mode)
    cfg <- config$modes[[mode]]
    if (is.null(cfg)) stop("No config for mode ", mode)

    files <- cfg$files
    inputs <- list()

    for (nm in names(files)) {
        rel <- files[[nm]]
        if (!nzchar(rel)) stop("Empty file path for ", nm)
        abs <- resolve_raw_path(config, rel)
        if (!file.exists(abs)) stop("File not found: ", abs)
        inputs[[nm]] <- read_table_auto(abs)
    }

    if (!is.null(cfg$engine)) inputs$engine <- cfg$engine
    inputs
}
