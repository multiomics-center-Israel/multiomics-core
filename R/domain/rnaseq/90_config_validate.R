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
        assert_scalar_chr(d$block_col, "de$block_col", allow_null = TRUE)
        assert_scalar_bool(d$block_pca, "de$block_pca", allow_null = TRUE)
        if (!is.null(d$block_model)) {
            assert_one_of(tolower(d$block_model), "de$block_model", c("fixed", "random"))
            if (identical(resolve_rna_block_model(d), "random") &&
                !identical(tolower(d$method %||% "deseq2"), "limma_voom")) {
                stop("de$block_model is 'random', but de$method is '", d$method %||% "deseq2",
                     "'. Only limma_voom fits a random block; DESeq2 and edgeR have no random effect.\n",
                     "  Use de$method 'limma_voom' for a random block, or de$block_model 'fixed'.")
            }
        }
        # Only "limma_voom" and "edger_ql" change the count-path fit; any other
        # value runs DESeq2, as it always has, so older configs keep working.
        de_method <- tolower(d$method %||% "deseq2")
        if (de_method %in% c("limma_voom", "edger_ql") && (is_preprocessed || uses_tximport)) {
            stop("de$method is '", de_method, "', which needs a raw count matrix.\n",
                 "  This config uses ",
                 if (is_preprocessed) "preprocessed input (input$format: preprocessed)" else
                     "tximport input (files$txi)",
                 ". Use de$method 'deseq2', or supply a gene-level count matrix.")
        }
        if (!de_method %in% c("deseq2", "limma_voom", "edger_ql") && !is_preprocessed) {
            warning("de$method '", d$method, "' is not an RNA-seq count method; ",
                    "DESeq2 will be fitted. Use 'deseq2', 'limma_voom' or 'edger_ql'.", call. = FALSE)
        }
        # Blocking is wired into the count paths only (DESeq2 and edgeR fixed
        # term, limma-voom fixed or duplicateCorrelation); the preprocessed (limma) and
        # pre-computed paths would silently fit or load an unpaired result,
        # hiding a design error.
        if (!is.null(d$block_col) && nzchar(d$block_col)) {
            has_precomputed <- length(unlist(cfg$files$de_table)) > 0 &&
                all(nzchar(unlist(cfg$files$de_table)))
            if (is_preprocessed || has_precomputed) {
                stop("de$block_col is set, but blocking applies only when DESeq2, ",
                     "limma-voom or edgeR is fitted on counts.\n",
                     "  This config uses ",
                     if (has_precomputed) "pre-computed DE tables (files$de_table)" else
                         "preprocessed input (input$format: preprocessed)",
                     ", so the pairing would be ignored. Remove de$block_col, or supply raw counts.")
            }
        }
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

    # 7. Pathway
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

    # 8. Enrichment (the `enrichment` block — offline/online RNA pathway enrichment)
    if (!is.null(cfg$enrichment)) {
        e <- cfg$enrichment
        assert_scalar_bool(e$enabled, "enrichment$enabled", allow_null = TRUE)
        assert_scalar_chr(e$annotation_dir, "enrichment$annotation_dir", allow_null = TRUE)
        assert_scalar_num(e$workers, "enrichment$workers", allow_null = TRUE, min_val = 1)
        assert_scalar_num(e$pvalue_cutoff, "enrichment$pvalue_cutoff",
                          allow_null = TRUE, min_val = 0, max_val = 1)
        assert_scalar_num(e$gsea_pvalue_cutoff, "enrichment$gsea_pvalue_cutoff",
                          allow_null = TRUE, min_val = 0, max_val = 1)
        assert_scalar_num(e$max_terms_in_dotplot, "enrichment$max_terms_in_dotplot",
                          allow_null = TRUE, min_val = 1)
        assert_one_of(e$padj_method, "enrichment$padj_method",
                      stats::p.adjust.methods, allow_null = TRUE)
        assert_one_of(e$gsea_padj_method, "enrichment$gsea_padj_method",
                      stats::p.adjust.methods, allow_null = TRUE)

        if (!is.null(e$databases)) {
            valid_db <- c("KEGG", "GO_BP", "GO_MF", "GO_CC")
            if (!is.character(e$databases) || length(e$databases) == 0) {
                stop("enrichment$databases must be a non-empty character vector, e.g. ",
                     "[\"KEGG\", \"GO_BP\", \"GO_MF\", \"GO_CC\"].")
            }
            bad_db <- setdiff(e$databases, valid_db)
            if (length(bad_db) > 0) {
                stop("enrichment$databases has unknown value(s): ",
                     paste(bad_db, collapse = ", "),
                     ". Valid: ", paste(valid_db, collapse = ", "), ".")
            }
        }

        if (!is.null(e$plots)) {
            if (!is.list(e$plots)) {
                stop("enrichment$plots must be a mapping of <toggle>: true/false.")
            }
            for (nm in c("dotplot", "ridgeplot", "ridgeplot_all_genes",
                         "pathway_heatmaps", "shared_genes")) {
                assert_scalar_bool(e$plots[[nm]], paste0("enrichment$plots$", nm),
                                   allow_null = TRUE)
            }
        }

        if (!is.null(e$gsea)) {
            if (!is.list(e$gsea)) stop("enrichment$gsea must be a mapping.")
            assert_scalar_bool(e$gsea$per_pathway_artifacts,
                               "enrichment$gsea$per_pathway_artifacts", allow_null = TRUE)
        }
        assert_scalar_bool(e$gsea_per_pathway_artifacts,
                           "enrichment$gsea_per_pathway_artifacts", allow_null = TRUE)
    }

    invisible(TRUE)
}
