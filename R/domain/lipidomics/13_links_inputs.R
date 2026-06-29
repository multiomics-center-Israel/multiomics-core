# R/domain/lipidomics/13_links_inputs.R
#
# Lipid-links extension: load and harmonize DE tables from the OTHER omics
# layers (RNAseq, proteomics, metabolomics) into a long-format tibble that
# downstream linkage functions consume.
#
# Output schema (one row per (layer, contrast, feature)):
#   layer        chr  one of: rnaseq | proteomics | metabolomics
#   contrast     chr  canonical lipid-side contrast label
#                     (e.g. "ppm1.56_vs_ppm0", "ppm12.5_vs_ppm0")
#   feature_id   chr  primary stable id (WBGene / UniProt / HMDB)
#   feature_name chr  display name (gene symbol / protein name / metabolite)
#   logFC        num  log2 fold change
#   p_value      num  raw p-value
#   p_adj        num  FDR / adjusted p-value
#   is_sig       lgl  significance flag (per upstream criterion)
#   entrez_id    chr  optional NCBI Entrez gene ID
#   wormbase_id  chr  optional WBGene id (for protein rows)
#   kegg_cpd     chr  optional KEGG compound id (for metabolite rows)
#   hmdb_id      chr  optional HMDB id (for metabolite rows)


# ---------------------------------------------------------------------------
# RNAseq DE table (Excel with multi-row header)
# ---------------------------------------------------------------------------
load_external_de_rnaseq <- function(path, contrasts_map, sheet = "DEG",
                                     header_row = 9) {
    stopifnot(file.exists(path))
    raw <- readxl::read_excel(path, sheet = sheet, col_names = FALSE,
                              skip = header_row - 1)
    header <- as.character(raw[1, ])
    body <- raw[-1, , drop = FALSE]
    colnames(body) <- header
    body <- body[!is.na(body[[1]]) & nzchar(as.character(body[[1]])), ,
                 drop = FALSE]

    out_rows <- list()
    for (lipid_ctr in names(contrasts_map)) {
        rna_ctr <- contrasts_map[[lipid_ctr]]
        fc <- as.numeric(body[[paste0(rna_ctr, "_log2FoldChange")]])
        p  <- as.numeric(body[[paste0(rna_ctr, "_pvalue")]])
        q  <- as.numeric(body[[paste0(rna_ctr, "_padj")]])
        sig_col <- body[[paste0(rna_ctr, "_is_sig")]]
        sig <- !is.na(sig_col) & tolower(sig_col) %in% c("yes", "true", "1")

        out_rows[[lipid_ctr]] <- data.frame(
            layer        = "rnaseq",
            contrast     = lipid_ctr,
            feature_id   = as.character(body[["gene_id"]]),
            feature_name = as.character(body[["gene_name"]]),
            logFC        = fc,
            p_value      = p,
            p_adj        = q,
            is_sig       = sig,
            entrez_id    = as.character(body[["Entrz_id"]]),
            wormbase_id  = NA_character_,
            kegg_cpd     = NA_character_,
            hmdb_id      = NA_character_,
            stringsAsFactors = FALSE
        )
    }
    do.call(rbind, out_rows)
}


# ---------------------------------------------------------------------------
# Proteomics DE table (Excel; header rows 1 + 2, use row 2)
# ---------------------------------------------------------------------------
load_external_de_proteomics <- function(path, contrasts_map, sheet = "DE",
                                         header_row = 2) {
    stopifnot(file.exists(path))
    raw <- readxl::read_excel(path, sheet = sheet, col_names = FALSE,
                              skip = header_row - 1)
    header <- as.character(raw[1, ])
    # Some columns have repeated/blank names — make unique
    header <- make.unique(ifelse(is.na(header) | header == "",
                                  paste0("V", seq_along(header)), header))
    body <- raw[-1, , drop = FALSE]
    colnames(body) <- header
    body <- body[!is.na(body[[1]]) & nzchar(as.character(body[[1]])), ,
                 drop = FALSE]

    out_rows <- list()
    for (lipid_ctr in names(contrasts_map)) {
        prot_ctr <- contrasts_map[[lipid_ctr]]
        fc <- as.numeric(body[[paste0(prot_ctr, "_ratio")]])
        p  <- as.numeric(body[[paste0(prot_ctr, "_p.val")]])
        q  <- as.numeric(body[[paste0(prot_ctr, "_p.adj")]])
        sig_col <- body[[paste0(prot_ctr, "_significant")]]
        sig <- !is.na(sig_col) & tolower(sig_col) %in% c("true", "yes", "1")

        uniprot <- as.character(body[["ID"]])
        wb <- as.character(body[["Wormbase_id"]])
        # Protein name: prefer the human-readable token before "_CAEEL"
        name_raw <- as.character(body[["Name"]])
        name_pretty <- sub("_CAEEL$", "", name_raw)

        out_rows[[lipid_ctr]] <- data.frame(
            layer        = "proteomics",
            contrast     = lipid_ctr,
            feature_id   = uniprot,
            feature_name = name_pretty,
            logFC        = fc,
            p_value      = p,
            p_adj        = q,
            is_sig       = sig,
            entrez_id    = NA_character_,
            wormbase_id  = wb,
            kegg_cpd     = NA_character_,
            hmdb_id      = NA_character_,
            stringsAsFactors = FALSE
        )
    }
    do.call(rbind, out_rows)
}


# ---------------------------------------------------------------------------
# Metabolomics DE table (tab-separated)
# ---------------------------------------------------------------------------
load_external_de_metabolomics <- function(path, contrasts_map) {
    stopifnot(file.exists(path))
    body <- utils::read.table(path, header = TRUE, sep = "\t",
                              quote = "", comment.char = "",
                              check.names = FALSE, stringsAsFactors = FALSE,
                              fill = TRUE, na.strings = c("", "NA"))

    # Parse id_orig: "Name|HMDB|CAS|KEGG|SMILES|level" — extract HMDB + KEGG
    id_orig <- as.character(body[["id_orig"]])
    parts <- strsplit(id_orig, "|", fixed = TRUE)
    pick <- function(i) vapply(parts, function(x) if (length(x) >= i) x[i]
                                                  else NA_character_,
                                character(1))
    kegg_cpd <- pick(4)
    kegg_cpd[!grepl("^C\\d{5}$", kegg_cpd)] <- NA_character_

    out_rows <- list()
    for (lipid_ctr in names(contrasts_map)) {
        met_ctr <- contrasts_map[[lipid_ctr]]
        fc <- as.numeric(body[[paste0(met_ctr, "_log2FC")]])
        p  <- as.numeric(body[[paste0(met_ctr, "_p_value")]])
        q  <- as.numeric(body[[paste0(met_ctr, "_fdr")]])
        sig <- !is.na(q) & q < 0.05

        out_rows[[lipid_ctr]] <- data.frame(
            layer        = "metabolomics",
            contrast     = lipid_ctr,
            feature_id   = ifelse(!is.na(body[["HMDB_ID"]]) &
                                  nzchar(body[["HMDB_ID"]]),
                                  as.character(body[["HMDB_ID"]]),
                                  as.character(body[["name_orig"]])),
            feature_name = as.character(body[["name_orig"]]),
            logFC        = fc,
            p_value      = p,
            p_adj        = q,
            is_sig       = sig,
            entrez_id    = NA_character_,
            wormbase_id  = NA_character_,
            kegg_cpd     = kegg_cpd,
            hmdb_id      = as.character(body[["HMDB_ID"]]),
            stringsAsFactors = FALSE
        )
    }
    do.call(rbind, out_rows)
}


# ---------------------------------------------------------------------------
# Top-level: load all configured external omics
# ---------------------------------------------------------------------------
load_external_omics_for_links <- function(links_cfg) {
    ext <- links_cfg$external_omics %||% list()
    out <- list()

    if (!is.null(ext$rnaseq)) {
        message("[links] loading external RNAseq DE: ", ext$rnaseq$path)
        out$rnaseq <- load_external_de_rnaseq(
            path = ext$rnaseq$path,
            contrasts_map = ext$rnaseq$contrasts,
            sheet = ext$rnaseq$sheet %||% "DEG",
            header_row = ext$rnaseq$header_row %||% 9
        )
    }
    if (!is.null(ext$proteomics)) {
        message("[links] loading external proteomics DE: ", ext$proteomics$path)
        out$proteomics <- load_external_de_proteomics(
            path = ext$proteomics$path,
            contrasts_map = ext$proteomics$contrasts,
            sheet = ext$proteomics$sheet %||% "DE",
            header_row = ext$proteomics$header_row %||% 2
        )
    }
    if (!is.null(ext$metabolomics)) {
        message("[links] loading external metabolomics DE: ",
                ext$metabolomics$path)
        out$metabolomics <- load_external_de_metabolomics(
            path = ext$metabolomics$path,
            contrasts_map = ext$metabolomics$contrasts
        )
    }

    # Enrich proteomics feature_name with the gene symbol from RNAseq when
    # available (UniProt accessions alone are useless for PubMed search).
    if (!is.null(out$proteomics) && !is.null(out$rnaseq)) {
        wb2sym <- unique(out$rnaseq[, c("feature_id", "feature_name"),
                                     drop = FALSE])
        wb2sym <- wb2sym[!is.na(wb2sym$feature_id) &
                          nzchar(wb2sym$feature_id) &
                          !is.na(wb2sym$feature_name) &
                          nzchar(wb2sym$feature_name), , drop = FALSE]
        sym_lut <- stats::setNames(wb2sym$feature_name, wb2sym$feature_id)
        sym_for_prot <- unname(sym_lut[out$proteomics$wormbase_id])
        # Keep the original UniProt-stripped name as a fallback when no
        # mapping exists.
        better <- !is.na(sym_for_prot) & nzchar(sym_for_prot)
        out$proteomics$feature_name[better] <- sym_for_prot[better]
    }

    out
}


# ---------------------------------------------------------------------------
# Convert lipid_de_res$de_tables into the same harmonized schema, then bind
# ---------------------------------------------------------------------------
harmonize_lipid_de <- function(lipid_de_res, row_data, contrast_labels = NULL) {
    de_tables <- lipid_de_res$de_tables %||% list()
    if (length(de_tables) == 0) return(NULL)

    if (is.null(contrast_labels)) contrast_labels <- names(de_tables)

    rd <- row_data
    if (is.null(rd) || !"feature_id" %in% colnames(rd)) {
        rd <- data.frame(feature_id = character(0), ClassKey = character(0),
                         stringsAsFactors = FALSE)
    }
    # Translated-CSV and processed-wide inputs carry lipid_class but not
    # ClassKey; fall back to it so lipids_to_pathways() (which drops rows with
    # missing ClassKey) still finds the parsed class annotations.
    if (!"ClassKey" %in% colnames(rd)) {
        rd$ClassKey <- if ("lipid_class" %in% colnames(rd)) {
            as.character(rd$lipid_class)
        } else {
            NA_character_
        }
    }

    rows <- list()
    for (ctr in contrast_labels) {
        tbl <- de_tables[[ctr]]
        if (is.null(tbl) || nrow(tbl) == 0) next
        # DE table columns: feature_id, Name, logFC, AveExpr, P.Value, adj.P.Val
        rows[[ctr]] <- data.frame(
            layer        = "lipidomics",
            contrast     = ctr,
            feature_id   = as.character(tbl$feature_id),
            feature_name = as.character(tbl$Name %||% tbl$feature_id),
            logFC        = as.numeric(tbl$logFC),
            p_value      = as.numeric(tbl$P.Value),
            p_adj        = as.numeric(tbl$adj.P.Val),
            is_sig       = !is.na(tbl$adj.P.Val) & tbl$adj.P.Val < 0.05,
            entrez_id    = NA_character_,
            wormbase_id  = NA_character_,
            kegg_cpd     = NA_character_,
            hmdb_id      = NA_character_,
            stringsAsFactors = FALSE
        )
    }
    if (length(rows) == 0) return(NULL)
    out <- do.call(rbind, rows)

    # Attach ClassKey via row_data
    cls <- rd[, c("feature_id", "ClassKey"), drop = FALSE]
    cls <- cls[!duplicated(cls$feature_id), , drop = FALSE]
    out <- merge(out, cls, by = "feature_id", all.x = TRUE, sort = FALSE)
    out
}
