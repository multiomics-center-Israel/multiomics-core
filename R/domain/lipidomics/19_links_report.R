# R/domain/lipidomics/19_links_report.R
#
# Render the lipid-links HTML report. Builds an inline Rmd from R objects
# the module already prepared (combined per-layer tables, raw pathway
# tables, literature paper rows). Knit with rmarkdown::render() if
# available; otherwise emit a static HTML.


render_lipid_links_report <- function(combined_list, pathway_links_list,
                                        literature_res, lipid_long,
                                        ext_long_list, links_cfg,
                                        out_dir,
                                        enzyme_lipid_figs = character(0)) {
    ensure_dir(out_dir)
    rmd_path  <- file.path(out_dir, "lipid_links_report.Rmd")
    html_path <- file.path(out_dir, "lipid_links_report.html")
    bundle_path <- file.path(out_dir, "lipid_links_report.bundle.rds")

    enz_pngs <- enzyme_lipid_figs[grepl(
        "enzyme_lipid_network\\.[^/]+\\.png$", enzyme_lipid_figs)]
    enz_pngs_rel <- vapply(enz_pngs, function(p) {
        # Path relative to the Rmd location for self-contained knitting
        f <- normalizePath(p, mustWork = FALSE)
        d <- normalizePath(out_dir, mustWork = FALSE)
        if (startsWith(f, d)) substr(f, nchar(d) + 2, nchar(f)) else f
    }, character(1))

    bundle <- list(
        combined_list      = combined_list,
        pathway_links_list = pathway_links_list,
        literature_res     = literature_res,
        lipid_long         = lipid_long,
        ext_long_list      = ext_long_list,
        links_cfg          = links_cfg,
        enzyme_lipid_pngs  = enz_pngs_rel,
        generated_at       = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )
    saveRDS(bundle, bundle_path)

    rmd <- .lipid_links_rmd_template(basename(bundle_path))
    writeLines(rmd, rmd_path)

    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        message("[links/report] rmarkdown not available — emitting Rmd only")
        return(c(rmd_path, bundle_path))
    }
    tryCatch({
        rmarkdown::render(rmd_path, output_file = basename(html_path),
                          quiet = TRUE, envir = new.env())
    }, error = function(e) {
        message("[links/report] rmarkdown::render failed: ",
                conditionMessage(e))
    })
    files <- c(rmd_path, bundle_path)
    if (file.exists(html_path)) files <- c(files, html_path)
    files
}


.lipid_links_rmd_template <- function(bundle_filename) {
    paste0(
"---
title: \"Lipid ↔ Genes / Proteins / Metabolites: Pathway- and Literature-based Links\"
output:
  html_document:
    toc: true
    toc_float: true
    theme: cosmo
    df_print: paged
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
bundle <- readRDS('", bundle_filename, "')
combined_list      <- bundle$combined_list
pathway_links_list <- bundle$pathway_links_list
literature_res     <- bundle$literature_res
lipid_long         <- bundle$lipid_long
ext_long_list      <- bundle$ext_long_list
links_cfg          <- bundle$links_cfg
n_top              <- links_cfg$output$top_n_per_layer %||% 50
fdr_lipid          <- links_cfg$output$fdr_cutoff_lipid %||% 0.05

fmt_table <- function(df, n = 30) {
  if (is.null(df) || nrow(df) == 0) return(htmltools::tags$em('No rows.'))
  head(df, n)
}
```

# Overview

Cross-omics linkage pipeline that connects significantly changed **lipids**
to candidate **genes / proteins / metabolites** in the same study via two
independent evidence channels:

* **Pathway co-membership** — lipid class → KEGG pathway → gene/compound
  members → layer-side DE features.
* **Literature** — PubMed query restricted to journals with IF ≥
  `r links_cfg$literature$min_journal_if` from `r links_cfg$literature$year_min`
  onward; co-mention required in the same sentence of the abstract.

Generated `r bundle$generated_at`.

# DE input summaries

## Lipidomics
```{r}
if (!is.null(lipid_long)) {
  sig <- lipid_long[!is.na(lipid_long$p_adj) &
                     lipid_long$p_adj < fdr_lipid, , drop = FALSE]
  cat('Significant lipids (FDR <', fdr_lipid, '):',
      nrow(sig), 'rows over', length(unique(sig$feature_id)),
      'unique features\\n')
  if (nrow(sig) > 0 && 'ClassKey' %in% colnames(sig)) {
    knitr::kable(table(sig$ClassKey, sig$contrast))
  }
}
```

## External omics layers
```{r}
if (!is.null(ext_long_list)) {
  for (lay in names(ext_long_list)) {
    e <- ext_long_list[[lay]]
    cat('* **', lay, '**: ',
        nrow(e), ' rows, ',
        sum(e$is_sig, na.rm = TRUE), ' significant\\n', sep = '')
  }
}
```

# Enzyme ↔ Lipid network (LipidONE-style)

Lipid classes (nodes) coloured by the mean log2FC of their species; enzymatic
reaction edges coloured by the best-matched **C. elegans ortholog**
log2FC from the external RNAseq / proteomics DE table. Two variants per
contrast:

* `any` — keeps reaction edges even when no DE data was found
  (greyed-out arrows) to show the structural backbone.
* `with_data` — keeps only reactions whose enzyme family had a matching
  DE feature in this study.

```{r enzyme-lipid-figs, results=\"asis\"}
pngs <- bundle$enzyme_lipid_pngs
if (length(pngs) == 0) {
  cat(\"_No enzyme-lipid figures were produced._\\n\")
} else {
  for (p in pngs) {
    cap <- sub(\"^figures/enzyme_lipid_network\\\\.\", \"\", p)
    cap <- sub(\"\\\\.png$\", \"\", cap)
    cat(sprintf('\\n![%s](%s)\\n\\n', cap, p))
  }
}
```

# Top combined-evidence hits per layer

```{r, results='asis'}
for (lay in names(combined_list)) {
  cat('\\n## ', lay, '\\n\\n', sep = '')
  df <- combined_list[[lay]]
  if (is.null(df) || nrow(df) == 0) {
    cat('No links.\\n'); next
  }
  show_cols <- c('contrast','lipid_class','lipid_name','lipid_logFC',
                 'feature_name','feature_logFC','sign_concordant',
                 'n_pathways','pathway_ids','n_papers','top_pmid',
                 'combined_score')
  show_cols <- intersect(show_cols, colnames(df))
  print(knitr::kable(head(df[, show_cols, drop = FALSE], n_top),
                     digits = 3,
                     caption = paste('Top', n_top, 'lipid ↔', lay,
                                      'links (combined score)')))
  cat('\\n')
}
```

# Pathway-only candidate tables

```{r, results='asis'}
for (lay in names(pathway_links_list)) {
  cat('\\n## ', lay, '\\n\\n', sep = '')
  df <- pathway_links_list[[lay]]
  if (is.null(df) || nrow(df) == 0) {
    cat('No pathway-based hits.\\n'); next
  }
  show_cols <- c('contrast','lipid_class','lipid_name','pathway_name',
                 'feature_name','feature_kegg_member','feature_logFC',
                 'sign_concordant','pathway_link_score')
  show_cols <- intersect(show_cols, colnames(df))
  print(knitr::kable(head(df[, show_cols, drop = FALSE], n_top),
                     digits = 3,
                     caption = paste('Top', n_top, lay,
                                      'pathway-link hits')))
  cat('\\n')
}
```

# Literature evidence

```{r}
papers <- literature_res$papers
if (is.null(papers) || nrow(papers) == 0) {
  cat('No same-sentence literature co-mentions found in journals with IF ≥',
      links_cfg$literature$min_journal_if, 'from',
      links_cfg$literature$year_min, 'onward.\\n')
} else {
  show_cols <- c('layer','contrast','lipid_class','feature_name','pmid',
                 'year','journal','journal_if','title','sentence')
  show_cols <- intersect(show_cols, colnames(papers))
  papers <- papers[order(-papers$journal_if, -papers$year), show_cols,
                    drop = FALSE]
  knitr::kable(head(papers, 100), digits = 1)
}
```

# Session info {.unnumbered}

```{r session-info-meta}
cat('Report generated: ', bundle$generated_at, '\\n', sep='')
cat('Working directory: ', normalizePath(getwd()), '\\n', sep='')
sys_cmd <- Sys.getenv('MULTIOMICS_CONFIG', unset = '')
if (nzchar(sys_cmd)) cat('MULTIOMICS_CONFIG: ', sys_cmd, '\\n', sep='')
cat('Hostname: ', Sys.info()[['nodename']], '\\n', sep='')
cat('User: ', Sys.info()[['user']], '\\n', sep='')
cat('System: ', Sys.info()[['sysname']], ' ', Sys.info()[['release']],
    '\\n', sep='')
```

<details>
<summary>**R session details** (click to expand)</summary>

```{r session-info}
sessionInfo()
```

</details>
"
    )
}
