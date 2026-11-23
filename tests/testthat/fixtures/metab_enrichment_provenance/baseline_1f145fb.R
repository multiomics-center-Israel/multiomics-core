# tests/testthat/fixtures/metab_enrichment_provenance/baseline_1f145fb.R
#
# FROZEN BASELINE. Do not edit.
#
# Verbatim copies of four functions from R/domain/metabolomics/06_enrichment.R
# as they were on main at 1f145fb, before provenance was added:
#   read_gmt_list(), map_compounds_for_enrichment(),
#   run_metabolomics_ora(), run_metabolomics_gsea().
# The only change is the `baseline_` prefix on their names and on the calls
# between them, so that they run beside the current functions. Every other
# helper they call is the live one; none of those is changed by the provenance
# work.
#
# test-metab-enrichment-provenance.R runs both versions on the same synthetic
# input and requires the same table, ranks, contrast and method: provenance
# must not move a single statistic.

baseline_read_gmt_list <- function(gmt_file, include_descriptions = FALSE) {
    # Return empty containers when the file does not exist
    if (!file.exists(gmt_file)) {
        if (include_descriptions) return(list(sets = list(), descriptions = character(0)))
        return(list())
    }

    lines    <- readLines(gmt_file, warn = FALSE)
    gmt_list <- list()
    desc_map <- character(0)

    # Each GMT line is: pathway_id <TAB> description <TAB> member1 <TAB> ...
    for (line in lines) {
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) < 3) next

        pathway_id   <- parts[1]
        pathway_desc <- parts[2]
        members      <- parts[3:length(parts)]
        # Drop empty strings that can appear from trailing tabs
        members      <- members[nzchar(members)]

        if (length(members) > 0) {
            gmt_list[[pathway_id]] <- members
            if (nzchar(pathway_desc)) {
                desc_map[[pathway_id]] <- pathway_desc
            }
        }
    }

    if (include_descriptions) {
        return(list(sets = gmt_list, descriptions = desc_map))
    }
    gmt_list
}

baseline_map_compounds_for_enrichment <- function(row_data, expr_mat, mapping_file = NULL,
                                         annotated_only = TRUE) {
    # Align the expression matrix to the annotation rows. row_data can be a
    # filtered subset of expr_mat (e.g. missingness filtering drops features, so
    # expr_raw ends up with more rows than row_data); matching by feature id
    # keeps compound_names (taken from row_data) and the matrix rows in lockstep.
    # Without this, the logical filtering below indexes a mismatched matrix and
    # errors, which the module's tryCatch turns into a silent NULL for every
    # enrichment method.
    if (!is.null(row_data) && "feature_id" %in% colnames(row_data) &&
        !is.null(rownames(expr_mat))) {
        common <- intersect(as.character(row_data$feature_id), rownames(expr_mat))
        row_data <- row_data[match(common, as.character(row_data$feature_id)), ,
                             drop = FALSE]
        expr_mat <- expr_mat[common, , drop = FALSE]
    }

    # Extract compound names from row_data
    if ("Name" %in% colnames(row_data)) {
        compound_names <- as.character(row_data$Name)
    } else if ("Molecule" %in% colnames(row_data)) {
        compound_names <- as.character(row_data$Molecule)
    } else {
        compound_names <- rownames(expr_mat)
    }

    # HMDB-to-KEGG mapping: replace compound names with KEGG IDs when a
    # mapping file is provided, so that IDs match KEGG-based GMT pathways.
    map_vec <- read_hmdb_kegg_map(mapping_file)
    if (!is.null(map_vec)) {
        # Locate HMDB IDs — try common column-name variants, then parse feature_id
        hmdb_ids <- NULL
        hmdb_col <- intersect(c("HMDB", "HMDB_ID", "HMDB ID", "hmdb_id"),
                              colnames(row_data))[1]
        if (!is.na(hmdb_col)) {
            hmdb_ids <- as.character(row_data[[hmdb_col]])
        } else if ("feature_id" %in% colnames(row_data)) {
            # Try extracting HMDB from feature_id (HMDB|Name format)
            hmdb_ids <- sub("\\|.*$", "", as.character(row_data$feature_id))
        }

        if (!is.null(hmdb_ids)) {
            mapped_kegg <- map_vec[hmdb_ids]
            has_map <- !is.na(mapped_kegg) & mapped_kegg != ""
            compound_names[has_map] <- mapped_kegg[has_map]
            message("enrichment: mapped ", sum(has_map), " features to KEGG IDs")
        }
    }

    # Prefer explicit KEGG compound IDs when the annotation carries them: a
    # KEGG column is unambiguous, whereas compound names are spelled
    # inconsistently across software. Applied AFTER the HMDB->KEGG step so a
    # direct KEGG ID wins over a lookup. Cells may hold "C00031", "cpd:C00031"
    # or "C00031;C00267" — take the first KEGG compound ID present.
    # mod_met_raw() resolves every accepted spelling (KEGG_ID, "KEGG ID", ...)
    # into `KEGG` with canonicalize_kegg_column(), so `KEGG` is the only column
    # read here.
    if ("KEGG" %in% colnames(row_data)) {
        kegg_raw <- as.character(row_data$KEGG)
        has_kegg <- grepl("C[0-9]{5}", kegg_raw)
        kegg_id  <- sub(".*?(C[0-9]{5}).*", "\\1", kegg_raw, perl = TRUE)
        compound_names[has_kegg] <- kegg_id[has_kegg]
        message("enrichment: using KEGG IDs from 'KEGG' for ",
                sum(has_kegg), " features")
    }

    # Optionally restrict the universe to features carrying a KEGG compound ID.
    # With a KEGG-based GMT, name-only features can never match a pathway, so
    # keeping them only inflates the background; MetaboAnalyst-style ORA uses the
    # KEGG-mappable set as the reference. Detect KEGG IDs by their Cxxxxx shape
    # (set above from the KEGG column or the HMDB->KEGG mapping).
    if (isTRUE(annotated_only)) {
        is_kegg <- grepl("^C[0-9]{5}$", trimws(compound_names))
        # No-op when the data carries no KEGG ids at all (e.g. a name-only
        # dataset with a name-based GMT), so the background is never emptied.
        if (any(is_kegg)) {
            n_drop <- sum(!is_kegg & !is.na(compound_names))
            compound_names[!is_kegg] <- NA_character_
            message("enrichment: restricting background to ", sum(is_kegg),
                    " KEGG-annotated features (dropped ", n_drop, " name-only).")
        }
    }

    # Per-feature id -> compound map, keyed by the feature ids callers hold
    # (ORA/GSEA look up by feature_id). compound_names is row-aligned to
    # row_data, so key on row_data$feature_id when available (falling back to
    # the matrix rownames). ORA/GSEA use this to translate their feature ids to
    # compounds; deriving that positionally after the filtering/dedup below
    # would misalign, so build it here beforehand.
    feat_key <- if ("feature_id" %in% colnames(row_data) &&
                    length(row_data$feature_id) == length(compound_names)) {
        as.character(row_data$feature_id)
    } else {
        rownames(expr_mat)
    }
    feature_map <- stats::setNames(compound_names, feat_key)

    # Filter valid compound names
    valid <- !is.na(compound_names) & nzchar(compound_names) & compound_names != "NA"
    mat_use   <- expr_mat[valid, , drop = FALSE]
    names_use <- compound_names[valid]

    # Deduplicate (keep first occurrence)
    dup <- duplicated(names_use)
    mat_use <- mat_use[!dup, , drop = FALSE]
    rownames(mat_use) <- names_use[!dup]

    list(
        expr_mapped    = mat_use,
        compound_names = names_use[!dup],
        feature_map    = feature_map
    )
}

baseline_run_metabolomics_ora <- function(pre, de_res, config, contrast = NULL) {
    cfg     <- config$modes$metabolomics
    enr_cfg <- cfg$enrichment %||% list()

    # Resolve gmt_file / mapping_file like the data files: absolute paths are
    # kept as-is, relative paths are located under paths$raw.
    enr_cfg$gmt_file     <- resolve_input_path(config, enr_cfg$gmt_file)
    enr_cfg$mapping_file <- resolve_input_path(config, enr_cfg$mapping_file)

    if (!isTRUE(enr_cfg$run_enrichment)) {
        message("metabolomics ORA: disabled — skipping")
        return(NULL)
    }

    # ---- Determine significance thresholds ----
    de_cfg       <- cfg$de %||% list()
    p_cutoff     <- de_cfg$p_cutoff %||% 0.05
    if (!is.null(de_cfg$logfc_cutoff)) {
        lfc_cutoff <- de_cfg$logfc_cutoff
    } else {
        lfc_cutoff <- log2(de_cfg$linear_fc_cutoff %||% 1.5)
    }
    mapping_file <- enr_cfg$mapping_file

    # ---- Build background (all measured features mapped to IDs) ----
    mapped <- baseline_map_compounds_for_enrichment(pre$row_data, pre$expr_raw, mapping_file)
    bg_ids <- mapped$compound_names
    if (length(bg_ids) < 5) {
        message("metabolomics ORA: too few background compounds — skipping")
        return(NULL)
    }

    # ---- Build foreground (significant features) ----
    if (is.null(de_res) || is.null(de_res$de_tables) ||
        length(de_res$de_tables) == 0) {
        message("metabolomics ORA: no DE results — skipping")
        return(NULL)
    }

    contrast_name <- .resolve_enrichment_contrast(contrast, de_res, "ORA")
    de_tbl <- de_res$de_tables[[contrast_name]]

    # Significance filtering: choose p-value column via config.
    # Defaults to "P.Value" (raw) — MetaboAnalyst-style ORA convention.
    # Set enr_cfg$pvalue_column = "adj.P.Val" for FDR-controlled filtering.
    pval_col <- enr_cfg$pvalue_column %||% "P.Value"
    if (!pval_col %in% colnames(de_tbl)) {
      stop("metabolomics ORA: pvalue_column '", pval_col,
           "' not found in DE table. Available columns: ",
           paste(colnames(de_tbl), collapse = ", "))
    }
    
    sig_mask <- !is.na(de_tbl[[pval_col]]) & de_tbl[[pval_col]] < p_cutoff &
      !is.na(de_tbl$logFC) & abs(de_tbl$logFC) >= lfc_cutoff
    sig_features <- de_tbl$feature_id[sig_mask]
    
    if (length(sig_features) == 0) {
      message("metabolomics ORA: no significant features at ", pval_col, "<", p_cutoff,
              " |logFC|>=", round(lfc_cutoff, 3), " — skipping")
      return(NULL)
    }

    # Map significant feature IDs to the compound namespace via the aligned
    # per-feature map (keyed by feature id), then keep those present in bg.
    fg_ids <- unique(na.omit(mapped$feature_map[sig_features]))
    fg_ids <- fg_ids[fg_ids %in% bg_ids]

    message("ORA: ", length(fg_ids), " foreground / ", length(bg_ids),
            " background compounds (contrast: ", contrast_name, ")")

    if (length(fg_ids) < 2) {
        message("metabolomics ORA: fewer than 2 mapped foreground compounds — skipping")
        return(NULL)
    }

    # ---- Load GMT sets ----
    gene_sets <- list()
    gmt_files <- unlist(enr_cfg$gmt_file)
    if (is.null(gmt_files)) {
        message("metabolomics ORA: no GMT files — skipping")
        return(NULL)
    }
    for (gf in gmt_files) {
        if (!file.exists(gf)) next
        gmt_parsed <- baseline_read_gmt_list(gf, include_descriptions = TRUE)
        gmt <- gmt_parsed$sets
        desc_map <- gmt_parsed$descriptions
        gmt <- translate_gmt_hmdb_to_kegg(gmt, mapping_file)
        names(gmt) <- make_pathway_labels(names(gmt), desc_map)
        lib_label <- tools::file_path_sans_ext(basename(gf))
        for (nm in names(gmt)) {
            gene_sets[[nm]] <- gmt[[nm]]
            attr(gene_sets[[nm]], "library") <- lib_label
        }
        message("ORA: loaded ", length(gmt), " sets from ", basename(gf))
    }

    if (length(gene_sets) == 0) {
        message("metabolomics ORA: no gene sets — skipping")
        return(NULL)
    }

    # Filter to sets with >= 2 members in background (case-insensitive)
    bg_map <- stats::setNames(bg_ids, tolower(bg_ids))
    gene_sets_filt <- lapply(gene_sets, function(cpds) {
        m <- bg_map[tolower(cpds)]
        unname(m[!is.na(m)])
    })
    keep <- vapply(gene_sets_filt, length, integer(1)) >= 2L
    gene_sets_filt <- gene_sets_filt[keep]
    gene_sets      <- gene_sets[keep]

    if (length(gene_sets_filt) == 0) {
        message("ORA: no pathways with >= 2 matching background compounds")
        return(NULL)
    }
    message("ORA: testing ", length(gene_sets_filt), " pathways")

    # ---- Fisher's exact test per pathway ----
    N <- length(bg_ids)   # total background
    K <- length(fg_ids)   # total foreground (significant)

    ora_rows <- lapply(names(gene_sets_filt), function(pw) {
        pw_members <- gene_sets_filt[[pw]]
        n  <- length(pw_members)                          # pathway size in bg
        k  <- sum(fg_ids %in% pw_members)                 # overlap fg ∩ pathway

        # 2x2 contingency: fg-in-pw, fg-not-in-pw, bg-in-pw-not-fg, bg-not-in-pw-not-fg
        mat <- matrix(c(k, K - k, n - k, N - K - (n - k)), nrow = 2)
        mat[mat < 0] <- 0  # safety

        pval <- tryCatch(
            stats::fisher.test(mat, alternative = "greater")$p.value,
            error = function(e) NA_real_
        )

        lib <- attr(gene_sets[[pw]], "library") %||% NA_character_

        data.frame(
            pathway     = pw,
            overlap     = k,
            pathway_size = n,
            fg_size     = K,
            bg_size     = N,
            fold_enrichment = (k / K) / (n / N),
            raw_p       = pval,
            library     = lib,
            stringsAsFactors = FALSE
        )
    })

    result_df <- do.call(rbind, ora_rows)
    result_df$FDR <- stats::p.adjust(result_df$raw_p, method = "fdr")
    result_df$overlap_genes <- vapply(names(gene_sets_filt), function(pw) {
        paste(fg_ids[fg_ids %in% gene_sets_filt[[pw]]], collapse = ";")
    }, character(1))
    result_df <- result_df[order(result_df$FDR), ]

    n_sig <- sum(result_df$FDR < 0.05, na.rm = TRUE)
    message("ORA complete: ", nrow(result_df), " pathways, ",
            n_sig, " with FDR < 0.05")

    list(
        table    = result_df,
        contrast = contrast_name,
        method   = "fisher_ora"
    )
}

baseline_run_metabolomics_gsea <- function(pre, de_res, config, contrast = NULL) {
    cfg     <- config$modes$metabolomics
    enr_cfg <- cfg$enrichment %||% list()

    # Resolve gmt_file / mapping_file like the data files: absolute paths are
    # kept as-is, relative paths are located under paths$raw.
    enr_cfg$gmt_file     <- resolve_input_path(config, enr_cfg$gmt_file)
    enr_cfg$mapping_file <- resolve_input_path(config, enr_cfg$mapping_file)

    if (!isTRUE(enr_cfg$run_enrichment)) {
        message("metabolomics GSEA: disabled — skipping")
        return(NULL)
    }

    if (!requireNamespace("fgsea", quietly = TRUE)) {
        message("metabolomics GSEA: fgsea not available — skipping")
        return(NULL)
    }

    if (is.null(de_res) || is.null(de_res$de_tables) ||
        length(de_res$de_tables) == 0) {
        message("metabolomics GSEA: no DE results — skipping")
        return(NULL)
    }

    mapping_file <- enr_cfg$mapping_file

    # ---- Build ranked list from DE log2FC ----
    contrast_name <- .resolve_enrichment_contrast(contrast, de_res, "GSEA")
    de_tbl <- de_res$de_tables[[contrast_name]]

    # Map feature IDs to the compound namespace via the aligned per-feature map.
    mapped <- baseline_map_compounds_for_enrichment(pre$row_data, pre$expr_raw, mapping_file)
    de_tbl$compound_id <- mapped$feature_map[de_tbl$feature_id]
    de_tbl <- de_tbl[!is.na(de_tbl$compound_id) & nzchar(de_tbl$compound_id), ]
    de_tbl <- de_tbl[!is.na(de_tbl$logFC), ]

    if (nrow(de_tbl) < 5) {
        message("metabolomics GSEA: too few mapped compounds — skipping")
        return(NULL)
    }

    # Ranking metric: log2 fold-change
    ranks <- stats::setNames(de_tbl$logFC, de_tbl$compound_id)

    # Deduplicate: keep entry with largest absolute logFC
    if (any(duplicated(names(ranks)))) {
        ord <- order(abs(ranks), decreasing = TRUE)
        ranks <- ranks[ord]
        ranks <- ranks[!duplicated(names(ranks))]
    }

    ranks <- sort(ranks, decreasing = TRUE)
    message("GSEA: ranked list of ", length(ranks),
            " compounds by log2FC (contrast: ", contrast_name, ")")

    # ---- Load GMT sets ----
    gene_sets <- list()
    gmt_files <- unlist(enr_cfg$gmt_file)
    if (is.null(gmt_files)) {
        message("metabolomics GSEA: no GMT files — skipping")
        return(NULL)
    }
    for (gf in gmt_files) {
        if (!file.exists(gf)) next
        gmt_parsed <- baseline_read_gmt_list(gf, include_descriptions = TRUE)
        gmt <- gmt_parsed$sets
        desc_map <- gmt_parsed$descriptions
        gmt <- translate_gmt_hmdb_to_kegg(gmt, mapping_file)
        names(gmt) <- make_pathway_labels(names(gmt), desc_map)
        lib_label <- tools::file_path_sans_ext(basename(gf))
        for (nm in names(gmt)) {
            gene_sets[[nm]] <- gmt[[nm]]
            attr(gene_sets[[nm]], "library") <- lib_label
        }
        message("GSEA: loaded ", length(gmt), " sets from ", basename(gf))
    }

    if (length(gene_sets) == 0) {
        message("metabolomics GSEA: no gene sets — skipping")
        return(NULL)
    }

    # Filter to sets with >= 2 members present in ranked list
    # (case-insensitive match — map GMT compound names to actual ranks keys)
    avail_map <- stats::setNames(names(ranks), tolower(names(ranks)))
    gene_sets_filt <- lapply(gene_sets, function(cpds) {
        m <- avail_map[tolower(cpds)]
        unname(m[!is.na(m)])
    })
    keep <- vapply(gene_sets_filt, length, integer(1)) >= 2L
    gene_sets_filt <- gene_sets_filt[keep]
    gene_sets      <- gene_sets[keep]

    if (length(gene_sets_filt) == 0) {
        message("GSEA: no pathways with >= 2 matching compounds in ranked list")
        return(NULL)
    }
    message("GSEA: testing ", length(gene_sets_filt), " pathways")

    # ---- Run fgsea ----
    fgsea_res <- tryCatch(
        fgsea::fgsea(pathways = gene_sets_filt, stats = ranks,
                     minSize = 2, maxSize = 500, nPermSimple = 10000),
        error = function(e) {
            warning("fgsea failed: ", e$message)
            NULL
        }
    )

    if (is.null(fgsea_res) || nrow(fgsea_res) == 0) {
        message("GSEA: no results from fgsea")
        return(NULL)
    }

    result_df <- data.frame(
        pathway        = fgsea_res$pathway,
        pval           = fgsea_res$pval,
        FDR            = fgsea_res$padj,
        NES            = fgsea_res$NES,
        ES             = fgsea_res$ES,
        pathway_size   = fgsea_res$size,
        stringsAsFactors = FALSE
    )

    result_df$library <- vapply(result_df$pathway, function(pw) {
        attr(gene_sets[[pw]], "library") %||% NA_character_
    }, character(1))

    result_df$leading_edge <- vapply(seq_len(nrow(fgsea_res)), function(i) {
        le <- fgsea_res$leadingEdge[[i]]
        if (is.null(le)) return("")
        paste(le, collapse = ";")
    }, character(1))

    result_df <- result_df[order(result_df$FDR), ]

    n_sig <- sum(result_df$FDR < 0.05, na.rm = TRUE)
    message("GSEA complete: ", nrow(result_df), " pathways, ",
            n_sig, " with FDR < 0.05")

    list(
        table    = result_df,
        ranks    = ranks,
        contrast = contrast_name,
        method   = "fgsea"
    )
}
