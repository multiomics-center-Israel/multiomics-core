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

    # 2. Files — counts and metadata are required
    if (!is.null(cfg$files)) {
        assert_scalar_chr(cfg$files$counts, "files$counts")
        assert_scalar_chr(cfg$files$metadata, "files$metadata")
        # sample_map and contrasts are optional
    }

    # 3. Normalization
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

    # 4. Filtering
    if (!is.null(cfg$filtering)) {
        f <- cfg$filtering
        assert_scalar_chr(f$group_col, "filtering$group_col", allow_null = TRUE)
        assert_scalar_num(f$threshold, "filtering$threshold", allow_null = TRUE, min_val = 0)
    }

    # 5. DE
    if (!is.null(cfg$de)) {
        d <- cfg$de
        assert_scalar_num(d$padj_cutoff, "de$padj_cutoff", allow_null = TRUE, min_val = 0, max_val = 1)
        assert_scalar_num(d$linear_fc_cutoff, "de$linear_fc_cutoff", allow_null = TRUE, min_val = 1)
    }

    # 6. Annotation (optional)
    if (!is.null(cfg$annotation)) {
        a <- cfg$annotation
        if (!is.null(a$source)) {
            assert_one_of(
                a$source,
                "annotation$source",
                c("Ensembl", "NCBI_GTF", "Generic"),
                allow_null = TRUE
            )
        }
        assert_scalar_chr(a$organism, "annotation$organism", allow_null = TRUE)
        assert_scalar_chr(a$id_type, "annotation$id_type", allow_null = TRUE)
        assert_scalar_bool(a$skip_annotation, "annotation$skip_annotation", allow_null = TRUE)
    }

    # 7. Batch Correction
    if (!is.null(cfg$batch_correction)) {
        bc <- cfg$batch_correction
        assert_scalar_bool(bc$enabled, "batch_correction$enabled")
        if (isTRUE(bc$enabled)) {
            assert_one_of(bc$method, "batch_correction$method",
                          c("combat_seq", "sva", "ruv"), allow_null = TRUE)
        }
    }

    # 8. Deconvolution
    if (!is.null(cfg$deconvolution)) {
        dc <- cfg$deconvolution
        assert_scalar_bool(dc$enabled, "deconvolution$enabled")
    }

    # 9. Pathway
    if (!is.null(cfg$pathway)) {
        p <- cfg$pathway
        assert_scalar_bool(p$enabled, "pathway$enabled", allow_null = TRUE)
        if (isTRUE(p$enabled)) {
            assert_one_of(p$method, "pathway$method",
                          c("fgsea", "ora", "gsea"), allow_null = TRUE)
            if (!is.null(p$min_size) && !is.null(p$max_size)) {
                assert_scalar_num(p$min_size, "pathway$min_size", min_val = 1)
                assert_scalar_num(p$max_size, "pathway$max_size", min_val = 1)
                if (p$min_size >= p$max_size) {
                    stop("pathway$min_size must be less than pathway$max_size")
                }
            }
        }
    }

    invisible(TRUE)
}
