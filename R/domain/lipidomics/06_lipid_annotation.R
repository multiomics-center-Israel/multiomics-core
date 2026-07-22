# R/domain/lipidomics/06_lipid_annotation.R
#
# Lipid class dictionary and annotation enrichment utilities.
#   1. lipid_class_dictionary()           - Reference table of lipid classes
#   2. annotate_lipid_row_data()          - Enrich row_data with class metadata
#   3. format_lipid_display_name()        - Human-readable lipid names for reports
#   4. compute_lipid_category_composition() - Category-level composition (GP/SP/GL/ST/FA)
#   5. plot_lipid_category_donut()        - Nested donut: categories (inner) + classes (outer)


# ==== LIPID CLASS DICTIONARY =================================================

#' Return a comprehensive lipid class dictionary
#'
#' Maps abbreviations to full names, categories, and LIPID MAPS two-letter
#' category codes.
#'
#' @return data.frame with columns: abbreviation, full_name, category,
#'   lipid_maps_category
lipid_class_dictionary <- function() {
    rows <- list(
        # ---- Glycerophospholipids (GP) ----
        c("PC",      "Phosphatidylcholine",                   "Glycerophospholipid", "GP"),
        c("PE",      "Phosphatidylethanolamine",              "Glycerophospholipid", "GP"),
        c("PS",      "Phosphatidylserine",                    "Glycerophospholipid", "GP"),
        c("PI",      "Phosphatidylinositol",                  "Glycerophospholipid", "GP"),
        c("PG",      "Phosphatidylglycerol",                  "Glycerophospholipid", "GP"),
        c("PA",      "Phosphatidic acid",                     "Glycerophospholipid", "GP"),
        c("CL",      "Cardiolipin",                           "Glycerophospholipid", "GP"),
        c("BMP",     "Bis(monoacylglycero)phosphate",         "Glycerophospholipid", "GP"),
        c("LPC",     "Lysophosphatidylcholine",               "Glycerophospholipid", "GP"),
        c("LPE",     "Lysophosphatidylethanolamine",          "Glycerophospholipid", "GP"),
        c("LPG",     "Lysophosphatidylglycerol",              "Glycerophospholipid", "GP"),
        c("LPI",     "Lysophosphatidylinositol",              "Glycerophospholipid", "GP"),
        c("LPS",     "Lysophosphatidylserine",                "Glycerophospholipid", "GP"),
        c("LPA",     "Lysophosphatidic acid",                 "Glycerophospholipid", "GP"),
        c("PIP",     "Phosphatidylinositol phosphate",        "Glycerophospholipid", "GP"),
        c("PIP2",    "Phosphatidylinositol bisphosphate",     "Glycerophospholipid", "GP"),
        c("PIP3",    "Phosphatidylinositol trisphosphate",    "Glycerophospholipid", "GP"),
        c("PEth",    "Phosphatidylethanol",                   "Glycerophospholipid", "GP"),
        c("PC-O",    "Plasmanyl-phosphatidylcholine",         "Glycerophospholipid", "GP"),
        c("PC-P",    "Plasmenyl-phosphatidylcholine",         "Glycerophospholipid", "GP"),
        c("PE-O",    "Plasmanyl-phosphatidylethanolamine",    "Glycerophospholipid", "GP"),
        c("PE-P",    "Plasmenyl-phosphatidylethanolamine",    "Glycerophospholipid", "GP"),
        c("ePC",     "Ether-phosphatidylcholine",             "Glycerophospholipid", "GP"),
        c("ePE",     "Ether-phosphatidylethanolamine",        "Glycerophospholipid", "GP"),

        # ---- Sphingolipids (SP) ----
        c("SM",      "Sphingomyelin",                         "Sphingolipid",        "SP"),
        c("Cer",     "Ceramide",                              "Sphingolipid",        "SP"),
        c("CerP",    "Ceramide-1-phosphate",                  "Sphingolipid",        "SP"),
        c("HexCer",  "Hexosylceramide",                       "Sphingolipid",        "SP"),
        c("Hex2Cer", "Dihexosylceramide (Lactosylceramide)",  "Sphingolipid",        "SP"),
        c("Hex3Cer", "Trihexosylceramide",                    "Sphingolipid",        "SP"),
        c("GlcCer",  "Glucosylceramide",                      "Sphingolipid",        "SP"),
        c("GalCer",  "Galactosylceramide",                    "Sphingolipid",        "SP"),
        c("LacCer",  "Lactosylceramide",                      "Sphingolipid",        "SP"),
        c("SulfoHexCer", "Sulfatide",                         "Sphingolipid",        "SP"),
        c("PE-Cer",  "Ceramide phosphoethanolamine",          "Sphingolipid",        "SP"),
        c("Sph",     "Sphingosine",                           "Sphingolipid",        "SP"),
        c("SphP",    "Sphingosine-1-phosphate",               "Sphingolipid",        "SP"),
        c("dhSph",   "Dihydrosphingosine (Sphinganine)",      "Sphingolipid",        "SP"),
        c("dhCer",   "Dihydroceramide",                       "Sphingolipid",        "SP"),
        c("phSM",    "Phytosphingomyelin",                    "Sphingolipid",        "SP"),
        c("GM1",     "Ganglioside GM1",                       "Sphingolipid",        "SP"),
        c("GM2",     "Ganglioside GM2",                       "Sphingolipid",        "SP"),
        c("GM3",     "Ganglioside GM3",                       "Sphingolipid",        "SP"),
        c("GD1a",    "Ganglioside GD1a",                      "Sphingolipid",        "SP"),
        c("GD1b",    "Ganglioside GD1b",                      "Sphingolipid",        "SP"),
        c("GD2",     "Ganglioside GD2",                       "Sphingolipid",        "SP"),
        c("GD3",     "Ganglioside GD3",                       "Sphingolipid",        "SP"),
        c("GT1b",    "Ganglioside GT1b",                      "Sphingolipid",        "SP"),
        c("GQ1b",    "Ganglioside GQ1b",                      "Sphingolipid",        "SP"),

        # ---- Glycerolipids (GL) ----
        c("TG",      "Triacylglycerol",                       "Glycerolipid",        "GL"),
        c("DG",      "Diacylglycerol",                        "Glycerolipid",        "GL"),
        c("MG",      "Monoacylglycerol",                      "Glycerolipid",        "GL"),
        c("MGDG",    "Monogalactosyldiacylglycerol",          "Glycerolipid",        "GL"),
        c("DGDG",    "Digalactosyldiacylglycerol",            "Glycerolipid",        "GL"),
        c("SQDG",    "Sulfoquinovosyldiacylglycerol",         "Glycerolipid",        "GL"),
        c("DGTS",    "Diacylglyceryl trimethylhomoserine",    "Glycerolipid",        "GL"),
        c("LDGTS",   "Lysodiacylglyceryl trimethylhomoserine","Glycerolipid",        "GL"),

        # ---- Sterol lipids (ST) ----
        c("CE",      "Cholesteryl ester",                     "Sterol lipid",        "ST"),
        c("Chol",    "Cholesterol",                           "Sterol lipid",        "ST"),
        c("ST",      "Sterol",                                "Sterol lipid",        "ST"),
        c("SE",      "Steryl ester",                          "Sterol lipid",        "ST"),
        c("BA",      "Bile acid",                             "Sterol lipid",        "ST"),
        c("SHex",    "Steryl hexoside",                       "Sterol lipid",        "ST"),
        c("BASulfate","Bile acid sulfate",                    "Sterol lipid",        "ST"),

        # ---- Fatty acyls (FA) ----
        c("CAR",     "Acylcarnitine",                         "Fatty acyl",          "FA"),
        c("FA",      "Fatty acid",                            "Fatty acyl",          "FA"),
        c("FAHFA",   "Fatty acid ester of hydroxyl fatty acid","Fatty acyl",         "FA"),
        c("WE",      "Wax ester",                             "Fatty acyl",          "FA"),
        c("NAE",     "N-acylethanolamine",                    "Fatty acyl",          "FA"),
        c("NAT",     "N-acyltaurine",                         "Fatty acyl",          "FA"),
        c("CoA",     "Acyl-CoA",                              "Fatty acyl",          "FA"),
        c("OAHFA",   "O-acyl hydroxy fatty acid",             "Fatty acyl",          "FA"),

        # ---- Prenol lipids (PR) ----
        c("CoQ",     "Coenzyme Q (Ubiquinone)",               "Prenol lipid",        "PR"),
        c("VAE",     "Vitamin A ester (Retinyl ester)",       "Prenol lipid",        "PR"),
        c("VDE",     "Vitamin D ester",                       "Prenol lipid",        "PR"),
        c("VEE",     "Vitamin E ester (Tocopherol ester)",    "Prenol lipid",        "PR"),

        # ---- Saccharolipids (SL) ----
        c("LipidA",  "Lipid A",                               "Saccharolipid",       "SL"),

        # ---- Polyketides (PK) ----
        c("AC",      "Acylphloroglucinol",                    "Polyketide",          "PK")
    )

    df <- data.frame(
        abbreviation       = vapply(rows, `[`, character(1), 1),
        full_name          = vapply(rows, `[`, character(1), 2),
        category           = vapply(rows, `[`, character(1), 3),
        lipid_maps_category = vapply(rows, `[`, character(1), 4),
        stringsAsFactors   = FALSE
    )
    rownames(df) <- NULL
    df
}


# ==== ANNOTATION ENRICHMENT ==================================================

#' Annotate row_data with class dictionary metadata
#'
#' Joins the lipid class dictionary onto an existing row_data data.frame,
#' adding full name, category, and LIPID MAPS category columns.
#'
#' @param row_data data.frame with at least feature_id and lipid_class columns.
#' @return data.frame with added columns: class_full_name, lipid_category,
#'   lipid_maps_category.
annotate_lipid_row_data <- function(row_data) {
    if (is.null(row_data) || nrow(row_data) == 0) return(row_data)

    dict <- lipid_class_dictionary()
    lookup <- stats::setNames(dict$full_name, dict$abbreviation)
    cat_lookup <- stats::setNames(dict$category, dict$abbreviation)
    maps_lookup <- stats::setNames(dict$lipid_maps_category, dict$abbreviation)

    cls <- as.character(row_data$lipid_class)

    row_data$class_full_name     <- unname(lookup[cls])
    row_data$lipid_category      <- unname(cat_lookup[cls])
    row_data$lipid_maps_category <- unname(maps_lookup[cls])

    # Replace NAs for unrecognised classes
    row_data$class_full_name[is.na(row_data$class_full_name)]         <- "Unknown"
    row_data$lipid_category[is.na(row_data$lipid_category)]           <- "Unknown"
    row_data$lipid_maps_category[is.na(row_data$lipid_maps_category)] <- "Unknown"

    # Add human-readable display name: "Phosphatidylcholine (16:0/18:1)" etc.
    row_data$display_name <- format_lipid_display_name(as.character(row_data$feature_id))

    # Also update the Name column so existing plotting/report code uses full names
    row_data$Name <- row_data$display_name

    row_data
}


# ==== DISPLAY NAME FORMATTER =================================================

#' Format lipid shorthand into a readable display name
#'
#' Converts e.g. "PC 16:0/18:1" -> "Phosphatidylcholine (16:0/18:1)".
#' Unknown classes are returned as-is.
#'
#' @param feature_id Character vector of lipid shorthand names.
#' @return Character vector of formatted display names.
format_lipid_display_name <- function(feature_id) {
    dict <- lipid_class_dictionary()
    lookup <- stats::setNames(dict$full_name, dict$abbreviation)

    vapply(feature_id, function(fid) {
        # Extract class abbreviation (everything before first space or parenthesis)
        abbr <- sub("^([A-Za-z0-9_-]+)[[:space:](].*", "\\1", fid)
        abbr <- trimws(abbr)

        full <- lookup[abbr]
        if (is.na(full)) return(fid)

        # Extract the chain/species part (everything after the class token)
        chain_part <- sub("^[A-Za-z0-9_-]+[[:space:]]*", "", fid)
        chain_part <- trimws(chain_part)

        if (nchar(chain_part) == 0) return(full)

        # Wrap in parentheses if not already
        if (!grepl("^\\(", chain_part)) {
            chain_part <- paste0("(", chain_part, ")")
        }

        paste(full, chain_part)
    }, character(1), USE.NAMES = FALSE)
}


# ==== CATEGORY-LEVEL COMPOSITION =============================================

#' Compute lipid category composition per sample group
#'
#' Like compute_lipid_class_composition() in 05_lipid_class_analysis.R but
#' aggregates at the LIPID MAPS category level (GP, SP, GL, ST, FA, etc.).
#'
#' @param pre    List from preprocess_lipidomics() (pre contract).
#' @param config Full pipeline config.
#' @return list(cat_abs, cat_norm, cat_count, per_sample_abs, per_sample_norm,
#'   categories)
compute_lipid_category_composition <- function(pre, config) {
    cfg <- config$modes$lipidomics

    condition_col <- cfg$de$condition_column %||% cfg$effects$color %||% "sample_type"
    sample_col    <- cfg$effects$samples %||% "sample_id"

    mat      <- pre$expr_filt
    meta     <- pre$meta
    row_data <- pre$row_data

    # Ensure category annotation is present
    if (!"lipid_maps_category" %in% colnames(row_data)) {
        row_data <- annotate_lipid_row_data(row_data)
    }

    meta      <- meta[match(colnames(mat), meta[[sample_col]]), , drop = FALSE]
    condition <- as.character(meta[[condition_col]])
    categories <- row_data$lipid_maps_category

    unique_cats <- sort(unique(categories))

    # Per-sample category totals
    cat_mat <- matrix(0, nrow = length(unique_cats), ncol = ncol(mat),
                      dimnames = list(unique_cats, colnames(mat)))

    for (cat in unique_cats) {
        idx <- which(categories == cat)
        if (length(idx) == 1) {
            cat_mat[cat, ] <- mat[idx, ]
        } else {
            cat_mat[cat, ] <- colSums(mat[idx, , drop = FALSE], na.rm = TRUE)
        }
    }

    # Normalized: each sample sums to 100%
    col_totals <- colSums(cat_mat, na.rm = TRUE)
    col_totals[col_totals == 0] <- 1
    cat_mat_norm <- sweep(cat_mat, 2, col_totals, FUN = "/") * 100

    # Per-group means
    groups <- unique(condition)
    cat_abs_mean <- sapply(groups, function(g) {
        idx <- which(condition == g)
        if (length(idx) == 1) cat_mat[, idx] else rowMeans(cat_mat[, idx, drop = FALSE], na.rm = TRUE)
    })
    colnames(cat_abs_mean) <- groups

    cat_norm_mean <- sapply(groups, function(g) {
        idx <- which(condition == g)
        if (length(idx) == 1) cat_mat_norm[, idx] else rowMeans(cat_mat_norm[, idx, drop = FALSE], na.rm = TRUE)
    })
    colnames(cat_norm_mean) <- groups

    cat_count <- table(categories)

    list(
        cat_abs         = cat_abs_mean,
        cat_norm        = cat_norm_mean,
        cat_count       = cat_count,
        per_sample_abs  = cat_mat,
        per_sample_norm = cat_mat_norm,
        categories      = unique_cats
    )
}


# ==== NESTED DONUT PLOT ======================================================

#' Nested donut (sunburst) plot: categories (inner ring) and classes (outer ring)
#'
#' If \code{group} is NULL the plot uses the overall mean across all samples.
#' If a group name is given, only that group's samples are used.
#'
#' @param cat_data List produced by compute_lipid_category_composition()
#'   **extended** with class-level detail.  In practice the caller should also
#'   pass the class composition data.  This function accepts a single list
#'   containing both \code{cat_norm} (category-level) and a supplementary
#'   \code{class_norm} matrix (from compute_lipid_class_composition()).
#'   When \code{class_norm} is absent, only the single-ring category donut is
#'   drawn.
#' @param group Character scalar. If NULL, averages across all groups.
#' @return ggplot object.
plot_lipid_category_donut <- function(cat_data, group = NULL) {

    # ---- category (inner) ring data ----
    cat_mat <- cat_data$cat_norm
    if (is.null(cat_mat) || length(cat_mat) == 0) {
        return(NULL)
    }

    if (!is.null(group) && group %in% colnames(cat_mat)) {
        cat_vals <- cat_mat[, group]
    } else {
        cat_vals <- if (ncol(cat_mat) == 1) {
            cat_mat[, 1]
        } else {
            rowMeans(cat_mat, na.rm = TRUE)
        }
    }

    cat_df <- data.frame(
        label = names(cat_vals),
        value = as.numeric(cat_vals),
        stringsAsFactors = FALSE
    )
    cat_df <- cat_df[cat_df$value > 0, , drop = FALSE]
    if (nrow(cat_df) == 0) return(NULL)

    # Compute positions for labels
    cat_df <- cat_df[order(-cat_df$value), , drop = FALSE]
    cat_df$fraction <- cat_df$value / sum(cat_df$value)
    cat_df$ymax <- cumsum(cat_df$fraction)
    cat_df$ymin <- c(0, utils::head(cat_df$ymax, -1))
    cat_df$ymid <- (cat_df$ymin + cat_df$ymax) / 2
    cat_df$ring <- "inner"

    # ---- class (outer) ring data ----
    class_mat <- cat_data$class_norm
    has_outer <- !is.null(class_mat) && length(class_mat) > 0

    if (has_outer) {
        if (!is.null(group) && group %in% colnames(class_mat)) {
            cls_vals <- class_mat[, group]
        } else {
            cls_vals <- if (ncol(class_mat) == 1) {
                class_mat[, 1]
            } else {
                rowMeans(class_mat, na.rm = TRUE)
            }
        }

        cls_df <- data.frame(
            label = names(cls_vals),
            value = as.numeric(cls_vals),
            stringsAsFactors = FALSE
        )
        cls_df <- cls_df[cls_df$value > 0, , drop = FALSE]

        if (nrow(cls_df) > 0) {
            # Map each class to its category
            dict <- lipid_class_dictionary()
            cls_to_cat <- stats::setNames(dict$lipid_maps_category, dict$abbreviation)
            cls_df$parent_cat <- unname(cls_to_cat[cls_df$label])
            cls_df$parent_cat[is.na(cls_df$parent_cat)] <- "Unknown"

            # Order: group by parent category to align with inner ring
            cat_order <- cat_df$label
            cls_df$cat_factor <- factor(cls_df$parent_cat, levels = cat_order)
            cls_df <- cls_df[order(cls_df$cat_factor, -cls_df$value), , drop = FALSE]

            cls_df$fraction <- cls_df$value / sum(cls_df$value)
            cls_df$ymax <- cumsum(cls_df$fraction)
            cls_df$ymin <- c(0, utils::head(cls_df$ymax, -1))
            cls_df$ymid <- (cls_df$ymin + cls_df$ymax) / 2
            cls_df$ring <- "outer"
        } else {
            has_outer <- FALSE
        }
    }

    # ---- category colour palette ----
    cat_colours <- c(
        GP = "#1f77b4", SP = "#ff7f0e", GL = "#2ca02c", ST = "#d62728",
        FA = "#9467bd", PR = "#8c564b", SL = "#e377c2", PK = "#7f7f7f",
        Unknown = "#bcbd22"
    )

    # ---- build plot ----
    title_suffix <- if (!is.null(group)) paste0(" - ", group) else ""
    plot_title <- paste0("Lipid Category & Class Composition", title_suffix)

    if (has_outer) {
        # Inner ring: x = 2 (width 1), outer ring: x = 3.5 (width 1)
        inner_df <- cat_df
        inner_df$xmin <- 1.5
        inner_df$xmax <- 2.5

        outer_df <- cls_df
        outer_df$xmin <- 3.0
        outer_df$xmax <- 4.0

        # Assign colours: inner by category, outer inherits lighter tint
        inner_df$fill <- cat_colours[inner_df$label]
        inner_df$fill[is.na(inner_df$fill)] <- "#bcbd22"

        # For outer ring: generate lighter shades per category
        outer_df$fill <- vapply(seq_len(nrow(outer_df)), function(i) {
            base <- cat_colours[outer_df$parent_cat[i]]
            if (is.na(base)) base <- "#bcbd22"
            # Lighten by mixing with white
            grDevices::adjustcolor(base, alpha.f = 0.7)
        }, character(1))

        combined <- rbind(
            data.frame(label = inner_df$label, ymin = inner_df$ymin,
                       ymax = inner_df$ymax, xmin = inner_df$xmin,
                       xmax = inner_df$xmax, ymid = inner_df$ymid,
                       fill = inner_df$fill, ring = "inner",
                       stringsAsFactors = FALSE),
            data.frame(label = outer_df$label, ymin = outer_df$ymin,
                       ymax = outer_df$ymax, xmin = outer_df$xmin,
                       xmax = outer_df$xmax, ymid = outer_df$ymid,
                       fill = outer_df$fill, ring = "outer",
                       stringsAsFactors = FALSE)
        )

        # Build named colour vector for legend (categories + classes)
        legend_fills <- c(
            stats::setNames(inner_df$fill, paste0(inner_df$label, " (", round(inner_df$fraction * 100, 1), "%)")),
            stats::setNames(outer_df$fill, outer_df$label)
        )
        combined$legend_label <- ifelse(
            combined$ring == "inner",
            paste0(combined$label, " (", round(combined$ymax * 100 - combined$ymin * 100, 1), "%)"),
            combined$label
        )

        p <- ggplot2::ggplot(combined) +
            ggplot2::geom_rect(
                ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                             fill = legend_label),
                colour = "white", linewidth = 0.3
            ) +
            ggplot2::scale_fill_manual(
                values = legend_fills,
                name = NULL
            ) +
            ggplot2::geom_text(
                data = combined[combined$ring == "inner", , drop = FALSE],
                ggplot2::aes(x = (xmin + xmax) / 2, y = ymid, label = label),
                size = 3.5, fontface = "bold"
            ) +
            ggplot2::geom_text(
                data = combined[combined$ring == "outer" &
                                combined$ymax - combined$ymin > 0.03, , drop = FALSE],
                ggplot2::aes(x = (xmin + xmax) / 2, y = ymid, label = label),
                size = 2.5
            ) +
            ggplot2::coord_polar(theta = "y") +
            ggplot2::xlim(0, 4.5) +
            ggplot2::theme_void() +
            ggplot2::labs(title = plot_title) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(face = "bold", size = 13,
                                                    hjust = 0.5),
                legend.position = "bottom",
                legend.text = ggplot2::element_text(size = 8)
            ) +
            ggplot2::guides(fill = ggplot2::guide_legend(ncol = 4))
    } else {
        # Single ring (category only)
        cat_df$xmin <- 2
        cat_df$xmax <- 4

        cat_df$fill <- cat_colours[cat_df$label]
        cat_df$fill[is.na(cat_df$fill)] <- "#bcbd22"

        cat_df$legend_label <- paste0(cat_df$label, " (", round(cat_df$fraction * 100, 1), "%)")
        legend_fills_single <- stats::setNames(cat_df$fill, cat_df$legend_label)

        p <- ggplot2::ggplot(cat_df) +
            ggplot2::geom_rect(
                ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                             fill = legend_label),
                colour = "white", linewidth = 0.5
            ) +
            ggplot2::scale_fill_manual(
                values = legend_fills_single,
                name = NULL
            ) +
            ggplot2::geom_text(
                ggplot2::aes(x = 3, y = ymid, label = label),
                size = 3.5, fontface = "bold"
            ) +
            ggplot2::coord_polar(theta = "y") +
            ggplot2::xlim(0, 4.5) +
            ggplot2::theme_void() +
            ggplot2::labs(title = plot_title) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(face = "bold", size = 13,
                                                    hjust = 0.5),
                legend.position = "bottom",
                legend.text = ggplot2::element_text(size = 8)
            ) +
            ggplot2::guides(fill = ggplot2::guide_legend(ncol = 4))
    }

    p
}
