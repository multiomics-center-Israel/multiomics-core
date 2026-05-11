# R/core/16_pptx_helpers.R
#
# Shared helpers for PowerPoint generation across all omics domains.
# Provides text formatting, slide templates, and AI bottom-line generation.
#
# Requires: officer

# ==============================================================================
# AI BOTTOM LINE — Claude Code CLI
# ==============================================================================

#' Truncate bottom-line text to fit in the PPTX text box
#'
#' Limits text to max_chars, cutting at the last sentence boundary that fits.
#' @param txt Character string.
#' @param max_chars Maximum allowed characters (default 350).
#' @return Truncated string.
truncate_bl <- function(txt, max_chars = 350) {
    if (nchar(txt) <= max_chars) return(txt)
    # Try to cut at last sentence boundary within limit
    short <- substr(txt, 1, max_chars)
    last_period <- max(gregexpr("\\.", short)[[1]])
    if (last_period > 0 && last_period > max_chars * 0.5) {
        return(substr(txt, 1, last_period))
    }
    paste0(substr(txt, 1, max_chars - 3), "...")
}


#' Generate an AI bottom-line summary for a slide
#'
#' @param slide_context Character string describing what the slide shows.
#' @param stats_json    Optional JSON string with relevant statistics.
#' @param image_path    Optional path to a PNG image for vision-based analysis.
#' @param config        Pipeline config (for model selection).
#' @return Character string with the AI bottom line.
generate_slide_bottom_line <- function(slide_context, stats_json = NULL,
                                       image_path = NULL, config = NULL,
                                       mode = NULL) {

    # Resolve mode-scoped commentary settings (same pattern as 12_commentary.R)
    if (!is.null(mode) && !is.null(config$modes[[mode]]$commentary)) {
        comment_cfg <- config$modes[[mode]]$commentary
    } else {
        comment_cfg <- config$commentary %||% list()
    }

    # Respect the enabled flag
    if (!isTRUE(comment_cfg$enabled)) return("")

    # Respect the backend setting
    backend <- comment_cfg$backend %||% "none"
    if (backend != "claude-code") return("")

    # Respect commentary.enabled setting
    comment_cfg <- config$commentary %||% list()
    if (!isTRUE(comment_cfg$enabled)) return("")

    if (Sys.which("claude") == "") return("")

    model <- comment_cfg$claude_code_model %||% "sonnet"

    prompt_parts <- paste0(
        "You are an expert bioinformatics scientist writing a bottom-line summary ",
        "for a PowerPoint slide in an omics analysis presentation. ",
        "Write EXACTLY 1-2 short sentences (MAXIMUM 250 characters total). ",
        "No markdown, no bullet points. Be direct and specific. ",
        "Do NOT exceed 250 characters.\n\n",
        "SLIDE CONTENT:\n", slide_context
    )

    if (!is.null(stats_json) && nzchar(stats_json))
        prompt_parts <- paste0(prompt_parts, "\n\nSTATISTICS:\n", stats_json)
    if (!is.null(image_path) && file.exists(image_path))
        prompt_parts <- paste0(prompt_parts,
            "\n\nAlso read the image at '", image_path, "' to inform your summary.")

    cmd <- sprintf("unset CLAUDECODE; claude --print --model %s --no-session-persistence %s",
                    model, shQuote(prompt_parts))

    tryCatch({
        raw <- system(cmd, intern = TRUE, timeout = 120)
        txt <- trimws(paste(raw, collapse = " "))
        if (nzchar(txt)) truncate_bl(txt) else ""
    }, error = function(e) { message("  Claude CLI error: ", e$message); "" })
}


# ==============================================================================
# TEXT FORMATTING CONSTANTS & HELPERS
# ==============================================================================

PPTX_COL_TITLE  <- "#1B2A4A"
PPTX_COL_BODY   <- "#333333"
PPTX_COL_ACCENT <- "#2563EB"
PPTX_COL_SUBTLE <- "#6B7280"
PPTX_COL_UP     <- "#16A34A"
PPTX_COL_DOWN   <- "#DC2626"
PPTX_COL_BL_BG  <- "#F0F4FF"

PPTX_SLIDE_W <- 10
PPTX_SLIDE_H <- 7.5

pptx_fp_title <- function() {
    officer::fp_text(font.size = 20, bold = TRUE, color = PPTX_COL_TITLE,
                      font.family = "Calibri")
}

pptx_fp_body <- function(sz = 11) {
    officer::fp_text(font.size = sz, color = PPTX_COL_BODY, font.family = "Calibri")
}

pptx_fp_body_bold <- function(sz = 11) {
    officer::fp_text(font.size = sz, bold = TRUE, color = PPTX_COL_BODY,
                      font.family = "Calibri")
}

pptx_fp_label <- function() {
    officer::fp_text(font.size = 9, color = PPTX_COL_SUBTLE, font.family = "Calibri")
}

pptx_fp_stat_value <- function(color = PPTX_COL_ACCENT) {
    officer::fp_text(font.size = 22, bold = TRUE, color = color,
                      font.family = "Calibri")
}

pptx_fp_stat_label <- function() {
    officer::fp_text(font.size = 9, color = PPTX_COL_SUBTLE, font.family = "Calibri")
}

pptx_fp_bottom_line <- function() {
    officer::fp_text(font.size = 9, italic = TRUE, color = PPTX_COL_TITLE,
                      font.family = "Calibri")
}

pptx_fp_bl_prefix <- function() {
    officer::fp_text(font.size = 9, bold = TRUE, italic = TRUE, color = PPTX_COL_ACCENT,
                      font.family = "Calibri")
}

pptx_fp_footer <- function() {
    officer::fp_text(font.size = 7, color = PPTX_COL_SUBTLE, font.family = "Calibri")
}


#' Add a figure slide with title, image, subtitle, AI bottom line, and footer
#'
#' @param pptx     officer pptx object.
#' @param title    Slide title.
#' @param img      Path to PNG file.
#' @param subtitle Optional subtitle text.
#' @param bl       AI bottom line text.
#' @param footer   Footer text.
#' @return Updated pptx object.
pptx_add_figure_slide <- function(pptx, title, img, subtitle = NULL, bl = "",
                                   footer = "") {
    pptx <- officer::add_slide(pptx, layout = "Blank", master = "Office Theme")

    # Title
    pptx <- officer::ph_with(pptx,
        value = officer::fpar(officer::ftext(title, pptx_fp_title())),
        location = officer::ph_location(left = 0.5, top = 0.3, width = 9, height = 0.5))

    # Subtitle
    top_img <- 0.95
    if (!is.null(subtitle) && nzchar(subtitle)) {
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(officer::ftext(subtitle, pptx_fp_label())),
            location = officer::ph_location(left = 0.5, top = 0.85, width = 9, height = 0.3))
        top_img <- 1.15
    }

    # Image
    if (file.exists(img)) {
        img_w <- 7.5; img_h <- 4.3
        pptx <- officer::ph_with(pptx,
            value = officer::external_img(img, width = img_w, height = img_h),
            location = officer::ph_location(left = (10 - img_w) / 2, top = top_img,
                                             width = img_w, height = img_h))
    }

    # AI bottom line
    if (nzchar(bl)) {
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(
                officer::ftext("AI Summary: ", pptx_fp_bl_prefix()),
                officer::ftext(bl, pptx_fp_bottom_line())),
            location = officer::ph_location(left = 0.5, top = 5.7, width = 9, height = 1.1))
    }

    # Footer
    if (nzchar(footer)) {
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(officer::ftext(footer, pptx_fp_footer())),
            location = officer::ph_location(left = 0.3, top = 7.1, width = 9.4, height = 0.3))
    }

    pptx
}


#' Add a title slide
#'
#' @param pptx         officer pptx object.
#' @param main_title   Main title text.
#' @param project_name Project name (subtitle).
#' @param analyst      Analyst name.
#' @param date_str     Date string.
#' @param stats_text   One-line stats summary.
#' @param bl           AI bottom line.
#' @param footer       Footer text.
#' @return Updated pptx object.
pptx_add_title_slide <- function(pptx, main_title, project_name, analyst,
                                  date_str, stats_text, bl = "", footer = "") {
    pptx <- officer::add_slide(pptx, layout = "Blank", master = "Office Theme")

    pptx <- officer::ph_with(pptx,
        value = officer::fpar(officer::ftext(main_title,
            officer::fp_text(font.size = 28, bold = TRUE, color = PPTX_COL_TITLE,
                              font.family = "Calibri"))),
        location = officer::ph_location(left = 1, top = 1.5, width = 8, height = 0.8))

    pptx <- officer::ph_with(pptx,
        value = officer::fpar(officer::ftext(project_name,
            officer::fp_text(font.size = 18, color = PPTX_COL_ACCENT,
                              font.family = "Calibri"))),
        location = officer::ph_location(left = 1, top = 2.4, width = 8, height = 0.5))

    pptx <- officer::ph_with(pptx,
        value = officer::fpar(officer::ftext(paste0(analyst, "  |  ", date_str),
            officer::fp_text(font.size = 12, color = PPTX_COL_SUBTLE,
                              font.family = "Calibri"))),
        location = officer::ph_location(left = 1, top = 3.1, width = 8, height = 0.4))

    pptx <- officer::ph_with(pptx,
        value = officer::fpar(officer::ftext(stats_text, pptx_fp_body(11))),
        location = officer::ph_location(left = 1, top = 4.0, width = 8, height = 0.4))

    if (nzchar(bl)) {
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(
                officer::ftext("AI Summary: ", pptx_fp_bl_prefix()),
                officer::ftext(bl, pptx_fp_bottom_line())),
            location = officer::ph_location(left = 1, top = 5.0, width = 8, height = 1.0))
    }

    if (nzchar(footer)) {
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(officer::ftext(footer, pptx_fp_footer())),
            location = officer::ph_location(left = 0.3, top = 7.1, width = 9.4, height = 0.3))
    }

    pptx
}


#' Add a study design slide
#'
#' @param pptx    officer pptx object.
#' @param items   Named list of design items (label = value).
#' @param bl      AI bottom line.
#' @param footer  Footer text.
#' @return Updated pptx object.
pptx_add_design_slide <- function(pptx, items, bl = "", footer = "") {
    pptx <- officer::add_slide(pptx, layout = "Blank", master = "Office Theme")

    pptx <- officer::ph_with(pptx,
        value = officer::fpar(officer::ftext("Study Design", pptx_fp_title())),
        location = officer::ph_location(left = 0.5, top = 0.3, width = 9, height = 0.5))

    design_fpars <- lapply(items, function(item) {
        officer::fpar(
            officer::ftext(paste0(item$label, ":  "), pptx_fp_body_bold(12)),
            officer::ftext(item$value, pptx_fp_body(12)))
    })
    pptx <- officer::ph_with(pptx,
        value = do.call(officer::block_list, design_fpars),
        location = officer::ph_location(left = 1.5, top = 1.2, width = 7, height = 4.2))

    if (nzchar(bl)) {
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(
                officer::ftext("AI Summary: ", pptx_fp_bl_prefix()),
                officer::ftext(bl, pptx_fp_bottom_line())),
            location = officer::ph_location(left = 0.5, top = 5.7, width = 9, height = 1.1))
    }

    if (nzchar(footer)) {
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(officer::ftext(footer, pptx_fp_footer())),
            location = officer::ph_location(left = 0.3, top = 7.1, width = 9.4, height = 0.3))
    }

    pptx
}


#' Add a table slide using flextable
#'
#' @param pptx       officer pptx object.
#' @param title      Slide title.
#' @param display_df Data frame to display.
#' @param subtitle   Optional subtitle.
#' @param bl         AI bottom line.
#' @param footer     Footer text.
#' @return Updated pptx object.
pptx_add_table_slide <- function(pptx, title, display_df, subtitle = NULL,
                                  bl = "", footer = "") {
    pptx <- officer::add_slide(pptx, layout = "Blank", master = "Office Theme")

    pptx <- officer::ph_with(pptx,
        value = officer::fpar(officer::ftext(title, pptx_fp_title())),
        location = officer::ph_location(left = 0.5, top = 0.3, width = 9, height = 0.5))

    if (!is.null(subtitle) && nzchar(subtitle)) {
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(officer::ftext(subtitle, pptx_fp_label())),
            location = officer::ph_location(left = 0.5, top = 0.85, width = 9, height = 0.3))
    }

    # Table on left, AI sidebar on right to avoid overlap
    tbl_w <- if (nzchar(bl)) 6.0 else 9.0
    if (requireNamespace("flextable", quietly = TRUE)) {
        ft <- flextable::flextable(display_df)
        ft <- flextable::fontsize(ft, size = 8, part = "body")
        ft <- flextable::fontsize(ft, size = 9, part = "header")
        ft <- flextable::bold(ft, part = "header")
        ft <- flextable::color(ft, color = "#FFFFFF", part = "header")
        ft <- flextable::bg(ft, bg = PPTX_COL_TITLE, part = "header")
        ft <- flextable::autofit(ft)
        ft <- flextable::padding(ft, padding = 3, part = "all")
        pptx <- officer::ph_with(pptx, value = ft,
            location = officer::ph_location(left = 0.5, top = 1.2, width = tbl_w, height = 5.5))
    } else {
        tbl_text <- paste(
            apply(display_df, 1, function(r) paste(r, collapse = "    ")),
            collapse = "\n")
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(officer::ftext(tbl_text, pptx_fp_body(9))),
            location = officer::ph_location(left = 0.5, top = 1.2, width = tbl_w, height = 5.5))
    }

    if (nzchar(bl)) {
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(
                officer::ftext("AI Summary: ", pptx_fp_bl_prefix()),
                officer::ftext(bl, pptx_fp_bottom_line())),
            location = officer::ph_location(left = 6.8, top = 1.2, width = 2.9, height = 5.5))
    }

    if (nzchar(footer)) {
        pptx <- officer::ph_with(pptx,
            value = officer::fpar(officer::ftext(footer, pptx_fp_footer())),
            location = officer::ph_location(left = 0.3, top = 7.1, width = 9.4, height = 0.3))
    }

    pptx
}


#' Check if PowerPoint generation is enabled in config
#'
#' @param config Full pipeline config.
#' @param mode   Mode name ("metabolomics", "rna", "proteomics").
#' @return Logical.
pptx_enabled <- function(config, mode) {
    mode_cfg <- config$modes[[mode]]
    if (is.null(mode_cfg)) return(FALSE)
    # Check mode-specific setting
    gen <- mode_cfg$outputs$generate_pptx
    if (!is.null(gen)) return(isTRUE(gen))
    # Check top-level setting
    gen <- config$outputs$generate_pptx
    if (!is.null(gen)) return(isTRUE(gen))
    # Default: TRUE (always generate)
    TRUE
}
