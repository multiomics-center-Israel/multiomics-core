#' Validate RNA Configuration
validate_rna_config <- function(cfg) {
    assert_named_list(cfg, "modes$rna")

    # Detect if tximport is being used (txi file configured instead of counts)
    uses_tximport <- !is.null(cfg$files$txi) && nzchar(cfg$files$txi %||% "")

    # 1. ID Columns
    assert_named_list(cfg$id_columns, "rna$id_columns")

    # gene_id is required for raw counts, but optional for tximport
    # (tximport uses rownames for gene IDs)
    if (!uses_tximport) {
        assert_scalar_chr(cfg$id_columns$gene_id, "id_columns$gene_id", allow_null = TRUE)
    }
    # sample_col is optional (defaults to SampleID), map_from/map_to optional too

    # 2. Normalization
    if (!is.null(cfg$normalization)) {
        n <- cfg$normalization
        assert_one_of(n$method, "normalization$method", c("TMMlogCPM", "VST"), allow_null = TRUE)
        assert_scalar_num(n$prior.count, "normalization$prior.count", allow_null = TRUE, min_val = 0)

        # tximport requires VST (TMMlogCPM not supported)
        if (uses_tximport && identical(n$method, "TMMlogCPM")) {
            stop(
                "[config] normalization$method = 'TMMlogCPM' is not compatible with tximport input. ",
                "Use 'VST' instead.",
                call. = FALSE
            )
        }
    }

    # 3. Filtering
    if (!is.null(cfg$filtering)) {
        f <- cfg$filtering
        assert_scalar_chr(f$group_col, "filtering$group_col", allow_null = TRUE)
        assert_scalar_num(f$threshold, "filtering$threshold", allow_null = TRUE, min_val = 0)
    }

    # 4. DE
    if (!is.null(cfg$de)) {
        d <- cfg$de
        assert_scalar_num(d$padj_cutoff, "de$padj_cutoff", allow_null = TRUE, min_val = 0, max_val = 1)
        assert_scalar_num(d$linear_fc_cutoff, "de$linear_fc_cutoff", allow_null = TRUE, min_val = 1)
    }

    invisible(TRUE)
}
