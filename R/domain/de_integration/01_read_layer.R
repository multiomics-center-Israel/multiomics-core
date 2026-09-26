# R/domain/de_integration/01_read_layer.R
#
# Turn one layer's DE table into standardized per-contrast frames: signed
# log2FC, p-values, hit flags, observed counts and a matching key per feature.
# Column names and contrasts are resolved in 00_inputs.R, observed counts in
# 02_observed_counts.R.

#' A numeric column of a layer's table, checked
#'
#' \code{as.numeric()} turns text into NA with only a warning, and nothing
#' stops a p-value of -0.1 from passing \code{p <= cutoff}; either would then
#' be reported as a result. A column that is not numeric, or holds values
#' outside its range, stops the run instead. Missing values stay missing.
#'
#' @param df The layer's table.
#' @param col Column name.
#' @param where Which layer and file, for the error message.
#' @param lower,upper Allowed range of the non-missing values.
#' @param finite Whether non-missing values must also be finite.
#' @param what What the column holds, for the error message.
#' @return Numeric vector.
.dei_numeric_col <- function(df, col, where, lower = -Inf, upper = Inf, finite = FALSE,
                             what = "numbers") {
    x <- df[[col]]
    if (!is.numeric(x) && !all(is.na(x))) {
        stop(where, ": column '", col, "' must hold ", what, " but is not numeric.",
             call. = FALSE)
    }
    x <- as.numeric(x)
    ok <- is.na(x) | (x >= lower & x <= upper & (!finite | is.finite(x)))
    if (!all(ok)) {
        stop(where, ": column '", col, "' must hold ", what, "; it has values ",
             if (finite) "that are not finite or " else "", "outside [", lower, ", ",
             upper, "].", call. = FALSE)
    }
    x
}


#' Signed log2 fold changes for one contrast, and where they came from
#'
#' In order of preference: a stored log2FC column; the unrounded
#' \code{linearRatio} with its sign taken from the signed linear fold change;
#' the signed linear fold change itself, which is rounded to three significant
#' digits. The middle step matters because the pre-computed proteomics loader
#' writes \code{linearRatio = 2^|log2FC|}, unsigned
#' (\code{load_precomputed_proteomics_de()}), so a plain \code{log2()} of it
#' would turn every decrease into an increase. \code{sign(linearFC) *
#' |log2(linearRatio)|} is exact for that export and for the main one alike.
#'
#' @param df The layer's table.
#' @param cols Output of \code{de_layer_columns()}.
#' @param where Which layer and file, for error messages.
#' @return List with \code{values} (numeric, one per row) and \code{source}
#'   (text naming the column or rule used).
layer_log2fc <- function(df, cols, where = "The layer's table") {
    num <- function(nm) .dei_numeric_col(df, nm, where, what = "fold changes")
    for (nm in cols$log2fc) {
        if (nm %in% names(df)) return(list(values = num(nm), source = nm))
    }
    has <- function(nm) !is.null(nm) && nm %in% names(df)
    if (has(cols$linear_ratio) && has(cols$linear_fc)) {
        ratio <- num(cols$linear_ratio)
        sgn <- sign(num(cols$linear_fc))
        return(list(values = sgn * abs(log2(ratio)),
                    source = sprintf("sign(%s) * |log2(%s)|", cols$linear_fc,
                                     cols$linear_ratio)))
    }
    if (has(cols$linear_fc)) {
        return(list(values = signed_fc_to_log2(num(cols$linear_fc)),
                    source = paste0(cols$linear_fc, " (3 significant digits)")))
    }
    stop("No fold-change column found; looked for: ",
         paste(c(cols$log2fc, cols$linear_ratio, cols$linear_fc), collapse = ", "),
         ".", call. = FALSE)
}


#' TRUE for the ways a table writes "yes"
#'
#' @param x Vector from a hit-flag column.
#' @return Logical vector, FALSE where \code{x} is missing.
.dei_truthy <- function(x) {
    if (is.logical(x)) return(!is.na(x) & x)
    if (is.numeric(x)) return(!is.na(x) & x != 0)
    v <- tolower(trimws(as.character(x)))
    !is.na(v) & v %in% c("1", "true", "t", "yes", "y", "pass")
}


#' Which features are hits for one contrast, and by what rule
#'
#' The table's own per-contrast flag when the config allows it and the table
#' has one: for the proteomics export that flag is the multi-imputation vote,
#' which a p-value and fold-change cutoff cannot reproduce. Otherwise the
#' cutoffs: (adjusted) p <= \code{p_cutoff} and |log2FC| >=
#' log2(\code{linear_fc_cutoff}). \code{pass_any_contrast} is never used: it
#' marks a hit in any contrast, not in this one.
#'
#' @param df The layer's table.
#' @param cols Output of \code{de_layer_columns()}.
#' @param log2fc Numeric vector from \code{layer_log2fc()}.
#' @param pvalue,padj Numeric vectors for this contrast.
#' @param hits_cfg The hits settings (global, overridden per layer).
#' @return List with \code{hit} (logical, never NA) and \code{source}.
layer_hit_flags <- function(df, cols, log2fc, pvalue, padj, hits_cfg) {
    if (isTRUE(hits_cfg$use_table_flag) && !is.null(cols$hit) &&
        cols$hit %in% names(df)) {
        v <- df[[cols$hit]]
        hit <- if (identical(cols$hit_kind, "updown")) {
            !is.na(v) & nzchar(trimws(as.character(v)))
        } else {
            .dei_truthy(v)
        }
        return(list(hit = hit, source = cols$hit))
    }
    p <- if (isTRUE(hits_cfg$use_adjusted)) padj else pvalue
    hit <- !is.na(p) & p <= hits_cfg$p_cutoff &
        !is.na(log2fc) & abs(log2fc) >= log2(hits_cfg$linear_fc_cutoff)
    list(hit = hit,
         source = sprintf("%s <= %s and |log2FC| >= log2(%s)",
                          if (isTRUE(hits_cfg$use_adjusted)) "padj" else "p",
                          hits_cfg$p_cutoff, hits_cfg$linear_fc_cutoff))
}


#' Matching key for each feature: a gene symbol, else a gene id
#'
#' Symbols are what lets a protein be paired with its transcript, or with the
#' same protein in another run whose protein groups were built differently.
#' Where a layer carries no symbol, its gene id stands in, and
#' \code{symbol_source} says so -- a gene id matches only the same id in
#' another layer, so for a non-model organism the ids have to share one
#' namespace (or a mapping file bridges them).
#'
#' Proteomics exports keep gene names in \code{Genes}, several for a group
#' ("A;B"); the first is taken, the convention of
#' \code{extract_protein_symbols()}. RNA exports carry gene ids only, so their
#' symbols come from \code{annotation_file} (the pipeline's
#' \code{Enrichment/gene_annotation.csv}: \code{gene_id}, \code{symbol}) when
#' given, and otherwise the gene id is the key.
#'
#' @param ly The layer's config.
#' @param df The layer's table.
#' @param cols Output of \code{de_layer_columns()}.
#' @param ids Feature ids, one per row of \code{df}.
#' @param annotation Annotation table read from \code{annotation_file}, or NULL.
#' @return Data frame with \code{symbol} (NA where nothing is known) and
#'   \code{symbol_source} ("symbol", "gene_id" or "none"), one row per row of
#'   \code{df}.
layer_symbols <- function(ly, df, cols, ids, annotation = NULL) {
    first_of <- function(x) {
        x <- trimws(vapply(strsplit(as.character(x), ";", fixed = TRUE),
                           function(p) if (length(p)) p[1] else NA_character_,
                           character(1)))
        x[!is.na(x) & !nzchar(x)] <- NA_character_
        x
    }
    symbol <- rep(NA_character_, nrow(df))
    source <- rep("none", nrow(df))

    sym_col <- cols$symbol[cols$symbol %in% names(df)][1]
    if (length(sym_col) == 1 && !is.na(sym_col)) {
        symbol <- first_of(df[[sym_col]])
        source[!is.na(symbol)] <- "symbol"
    }

    if (!is.null(annotation)) {
        ann <- first_of(annotation$symbol)[match(ids, as.character(annotation$gene_id))]
        fill <- is.na(symbol) & !is.na(ann)
        symbol[fill] <- ann[fill]
        source[fill] <- "symbol"
    }

    gid_col <- cols$gene_id[cols$gene_id %in% names(df)][1]
    if (length(gid_col) == 1 && !is.na(gid_col)) {
        gid <- first_of(df[[gid_col]])
        fill <- is.na(symbol) & !is.na(gid)
        symbol[fill] <- gid[fill]
        source[fill] <- "gene_id"
    }

    data.frame(symbol = symbol, symbol_source = source, stringsAsFactors = FALSE)
}


#' Read a layer's gene annotation, and check it
#'
#' A configured \code{annotation_file} is how an RNA layer gets symbols, so a
#' file without its \code{gene_id} and \code{symbol} columns stops the run
#' rather than leaving every feature keyed on its gene id with nothing said.
#'
#' @param ly The layer's config.
#' @param config Full config, for path resolution.
#' @return The annotation table, or NULL when the layer names none.
read_layer_annotation <- function(ly, config) {
    if (is.null(ly$annotation_file)) return(NULL)
    path <- resolve_input_path(config, ly$annotation_file)
    ann <- read_table_auto(path)
    gap <- setdiff(c("gene_id", "symbol"), names(ann))
    if (length(gap) > 0) {
        stop("Layer '", ly$name, "': annotation_file (", path, ") lacks column(s): ",
             paste(gap, collapse = ", "), ".", call. = FALSE)
    }
    # match() would silently take the first of two symbols for one gene id.
    pairs <- unique(data.frame(gene_id = as.character(ann$gene_id),
                               symbol = trimws(as.character(ann$symbol)),
                               stringsAsFactors = FALSE))
    pairs <- pairs[!is.na(pairs$gene_id), , drop = FALSE]
    clash <- unique(pairs$gene_id[duplicated(pairs$gene_id)])
    if (length(clash) > 0) {
        stop("Layer '", ly$name, "': annotation_file (", path, ") gives more than one ",
             "symbol for gene id(s) ", paste(utils::head(clash, 5), collapse = ", "),
             if (length(clash) > 5) ", ..." else "", ".", call. = FALSE)
    }
    ann
}


#' Read one layer's DE table into standardized per-contrast frames
#'
#' @param ly The layer's validated config.
#' @param contrasts Contrast labels to extract, as named in the comparisons.
#' @param config Full config, for path resolution.
#' @param hits_default The global \code{hits} settings; \code{ly$hits} overrides.
#' @param well_observed_min Observed values each group needs for a fold change
#'   to count as well observed.
#' @param df The table, when the caller has already read it; read from
#'   \code{ly$path} otherwise.
#' @return List with \code{name}, \code{label}, \code{omics_type},
#'   \code{format}, \code{path}, \code{tables} (named by requested contrast;
#'   each a data frame with \code{feature_id}, \code{symbol},
#'   \code{symbol_source}, \code{description}, \code{log2fc}, \code{pvalue},
#'   \code{padj}, \code{hit}, \code{n_obs_num}, \code{n_obs_den},
#'   \code{well_observed}) and \code{provenance} (one row per contrast).
read_de_layer <- function(ly, contrasts, config, hits_default,
                          well_observed_min = 2, df = NULL) {
    path <- resolve_input_path(config, ly$path)
    if (is.null(df)) df <- read_table_auto(path)
    cn <- names(df)
    available <- list_layer_contrasts(cn, ly$format, ly$contrast %||% "contrast")
    hits_cfg <- utils::modifyList(hits_default, ly$hits %||% list())
    observed <- read_observed_inputs(ly, config)
    annotation <- read_layer_annotation(ly, config)

    # A column the map names is a column the reader is told to use: absent, it
    # would quietly become BH, cutoffs or no symbols instead. The required ones
    # (id, pvalue, fold change) are reported below with their own messages.
    if (identical(ly$format, "generic")) {
        optional <- c("symbol", "gene_id", "description", "padj", "hit",
                      "n_obs_num", "n_obs_den")
        mapped <- unlist(ly$columns[intersect(optional, names(ly$columns))])
        gap <- mapped[!mapped %in% cn]
        if (length(gap) > 0) {
            stop("Layer '", ly$name, "': ",
                 paste0("columns.", names(gap), " names '", gap, "'", collapse = ", "),
                 ", not found in ", basename(path), ".", call. = FALSE)
        }
    }

    tables <- list()
    prov <- list()
    for (requested in unique(contrasts)) {
        contrast <- resolve_layer_contrast(requested, available, ly$name)
        cols <- de_layer_columns(ly$format, contrast, ly$columns, ly$id_col)

        id_col <- cols$id[cols$id %in% cn][1]
        if (is.na(id_col)) {
            stop("Layer '", ly$name, "': no feature id column (looked for ",
                 paste(cols$id, collapse = ", "), ") in ", basename(path), ".",
                 call. = FALSE)
        }
        if (!cols$pvalue %in% cn) {
            stop("Layer '", ly$name, "': column '", cols$pvalue, "' not found in ",
                 basename(path), " for contrast '", contrast, "'.", call. = FALSE)
        }
        where <- sprintf("Layer '%s' (%s)", ly$name, basename(path))
        ids <- as.character(df[[id_col]])
        lfc <- layer_log2fc(df, cols, where)
        pvalue <- .dei_numeric_col(df, cols$pvalue, where, 0, 1, finite = TRUE,
                                   what = "p-values")
        # Which it was goes into the provenance: a BH computed here, over the
        # rows this table kept, is not the adjusted p-value its producer reported.
        has_padj <- !is.null(cols$padj) && cols$padj %in% cn
        padj <- if (has_padj) {
            .dei_numeric_col(df, cols$padj, where, 0, 1, finite = TRUE,
                             what = "adjusted p-values")
        } else {
            stats::p.adjust(pvalue, method = "BH")
        }
        padj_source <- if (has_padj) cols$padj else paste0("BH of ", cols$pvalue)
        hit <- layer_hit_flags(df, cols, lfc$values, pvalue, padj, hits_cfg)
        obs <- layer_observed_counts(ly, df, ids, contrast, config, observed)
        sym <- layer_symbols(ly, df, cols, ids, annotation)
        desc <- if (!is.null(cols$description) && cols$description %in% cn) {
            as.character(df[[cols$description]])
        } else rep(NA_character_, nrow(df))

        tab <- data.frame(
            feature_id = ids, symbol = sym$symbol, symbol_source = sym$symbol_source,
            description = desc, log2fc = lfc$values, pvalue = pvalue, padj = padj,
            hit = hit$hit, n_obs_num = obs$num, n_obs_den = obs$den,
            stringsAsFactors = FALSE)
        tab$well_observed <- ifelse(is.na(tab$n_obs_num) | is.na(tab$n_obs_den), NA,
                                    pmin(tab$n_obs_num, tab$n_obs_den) >= well_observed_min)

        # One row per feature: a duplicated id would pair twice later on. The
        # smallest p-value is kept, then the first in id order, so reruns agree.
        # Radix sorts ids bytewise, so the order is the same under any locale.
        tab <- tab[!is.na(tab$feature_id) & nzchar(tab$feature_id), , drop = FALSE]
        tab <- tab[order(tab$feature_id, tab$pvalue, na.last = TRUE, method = "radix"), ,
                   drop = FALSE]
        n_dup <- sum(duplicated(tab$feature_id))
        if (n_dup > 0) {
            warning("Layer '", ly$name, "', contrast '", contrast, "': ", n_dup,
                    " duplicated feature id(s); kept the row with the smallest p-value.",
                    call. = FALSE)
            tab <- tab[!duplicated(tab$feature_id), , drop = FALSE]
        }
        rownames(tab) <- NULL

        tables[[requested]] <- tab
        prov[[length(prov) + 1]] <- data.frame(
            layer = ly$name, requested_contrast = requested, contrast = contrast,
            file = path, n_features = nrow(tab), n_hits = sum(tab$hit),
            log2fc_source = lfc$source, padj_source = padj_source,
            hit_source = hit$source,
            n_obs_source = obs$source,
            symbol_sources = paste(sprintf("%s: %d", names(table(tab$symbol_source)),
                                           as.integer(table(tab$symbol_source))),
                                   collapse = "; "),
            n_duplicates_dropped = n_dup, stringsAsFactors = FALSE)
    }

    list(name = ly$name, label = ly$label %||% ly$name, omics_type = ly$omics_type,
         format = ly$format, path = path, tables = tables,
         provenance = do.call(rbind, prov))
}
