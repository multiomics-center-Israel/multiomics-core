# R/domain/lipidomics/14_links_pathway_map.R
#
# Map lipid feature -> ClassKey -> KEGG pathway IDs using the curated
# data/lipidomics/lipid_class_kegg_pathways.tsv table.


load_lipid_class_pathway_map <- function(file = NULL) {
    if (is.null(file)) {
        # Try resolution against the running project's data/ dir.
        candidates <- c(
            file.path(getwd(), "data", "lipidomics",
                      "lipid_class_kegg_pathways.tsv"),
            file.path("/home/ozsol/multiomics-core", "data", "lipidomics",
                      "lipid_class_kegg_pathways.tsv")
        )
        file <- candidates[file.exists(candidates)][1]
    }
    if (is.na(file) || !nzchar(file) || !file.exists(file)) {
        stop("lipid_class_kegg_pathways.tsv not found")
    }
    map <- utils::read.table(file, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, check.names = FALSE,
                             quote = "")
    map
}


# DE-significant lipid -> pathway expansion.
# Returns long table: feature_id, ClassKey, kegg_pathway_id, pathway_name,
# contrast, lipid_logFC, lipid_padj.
lipids_to_pathways <- function(lipid_long, pathway_map,
                                fdr_cutoff = 0.05, organism = "cel") {
    sig <- lipid_long[!is.na(lipid_long$p_adj) &
                       lipid_long$p_adj < fdr_cutoff &
                       !is.na(lipid_long$ClassKey), , drop = FALSE]
    if (nrow(sig) == 0) {
        message("[links] no DE-significant lipids at FDR < ", fdr_cutoff)
        return(NULL)
    }
    pm <- pathway_map[pathway_map$organism == organism, , drop = FALSE]
    merged <- merge(
        sig[, c("feature_id", "feature_name", "ClassKey", "contrast",
                "logFC", "p_adj"), drop = FALSE],
        pm[, c("class_key", "kegg_pathway_id", "pathway_name")],
        by.x = "ClassKey", by.y = "class_key", all.x = FALSE
    )
    names(merged)[names(merged) == "logFC"] <- "lipid_logFC"
    names(merged)[names(merged) == "p_adj"] <- "lipid_padj"
    merged
}
