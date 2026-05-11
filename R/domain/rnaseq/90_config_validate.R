#' Validate RNA Configuration
validate_rna_config <- function(cfg) {
    assert_named_list(cfg, "modes$rna")

    # Detect input format
    uses_tximport <- !is.null(cfg$files$txi) && nzchar(cfg$files$txi %||% "")
    is_preprocessed <- identical(cfg$input$format, "preprocessed")

    # Validate input$format if provided
    if (!is.null(cfg$input$format)) {
        assert_one_of(cfg$input$format, "input$format", c("preprocessed"), allow_null = TRUE)
    }

    # 1. ID Columns
    assert_named_list(cfg$id_columns, "rna$id_columns")

    # gene_id is required for raw counts and preprocessed, but optional for tximport.
    # For preprocessed input, preprocess_rna() stops if gene_id is missing, so
    # enforce it here at validation time for a clearer error message.
    if (is_preprocessed) {
        assert_scalar_chr(cfg$id_columns$gene_id, "id_columns$gene_id", allow_null = FALSE)
    } else if (!uses_tximport) {
        assert_scalar_chr(cfg$id_columns$gene_id, "id_columns$gene_id", allow_null = TRUE)
    }

    # 2. Files — counts OR txi required; metadata always required
    if (!is.null(cfg$files)) {
        has_counts <- !is.null(cfg$files$counts) && nzchar(cfg$files$counts %||% "")
        has_txi <- uses_tximport
        has_preprocessed <- !is.null(cfg$files$preprocessed_counts) &&
                            nzchar(cfg$files$preprocessed_counts %||% "")
        if (!has_counts && !has_txi && !has_preprocessed) {
            stop("One of 'files$counts', 'files$txi', or 'files$preprocessed_counts' must be specified",
                 call. = FALSE)
        }
        assert_scalar_chr(cfg$files$metadata, "files$metadata")
    }

    # 3. Normalization (not required for preprocessed)
    if (!is.null(cfg$normalization) && !is_preprocessed) {
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
                          c("fgsea", "ora", "gsea", "both"), allow_null = TRUE)
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
