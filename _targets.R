# _targets.R
library(targets)

tar_option_set(packages = c("readr","dplyr","tidyr","ggplot2"))

# ---- Pipeline ----
list(
  # Config (edit paths later)
  tar_target(cfg_counts_path, "data/example_counts.csv", format = "file"),
  tar_target(cfg_meta_path,   "data/example_meta.csv",   format = "file"),
  
  # Load data
  tar_target(counts_raw, readr::read_csv(cfg_counts_path, show_col_types = FALSE)),
  tar_target(meta,       readr::read_csv(cfg_meta_path,   show_col_types = FALSE)),
  
  # Convert to matrix/data.frame as needed by your functions
  tar_target(counts_mat, {
    m <- as.matrix(counts_raw[,-1, drop = FALSE])
    rownames(m) <- counts_raw[[1]]
    m
  }),
  
  # Normalization step
  tar_target(
    counts_norm,
    {
      source("R/normalize.R")   # ensures the function is loaded
      normalize_counts(counts_mat, method = "TMMlogCPM")
    }),
  
  tar_target(
    counts_norm_csv,
    {
      dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
      out <- "outputs/counts_norm_TMMlogCPM.csv"
      readr::write_csv(
        tibble::tibble(gene = rownames(counts_norm)) |>
          dplyr::bind_cols(as.data.frame(counts_norm)),
        out
      )
      out
    },
    format = "file"
  ),
  
  tar_target(
    counts_norm_density_png,
    {
      dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
      out <- "outputs/counts_norm_density.png"
      
      # Build a long table for plotting
      df <- as.data.frame(counts_norm)
      df$gene <- rownames(counts_norm)
      long <- tidyr::pivot_longer(
        df, cols = -gene, names_to = "sample", values_to = "value"
      )
      long$log_value <- long$value  # already logCPM
      
      p <- ggplot2::ggplot(long, ggplot2::aes(x = log_value)) +
        ggplot2::geom_density(alpha = 0.4) +
        ggplot2::theme_bw() +
        ggplot2::labs(
          title = "Density of normalized expression (TMM logCPM)",
          x = "logCPM", y = "Density"
        )
      
      ggplot2::ggsave(out, p, width = 7, height = 4, dpi = 300)
      out
    },
    format = "file"
  ) ,
  tar_target(
    counts_pca_png,
    {
      dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
      out <- "outputs/counts_pca.png"
      
      # counts_norm is genes x samples (logCPM). PCA is on samples:
      mat <- t(scale(t(as.matrix(counts_norm))))   # z-score per gene
      pca <- prcomp(t(mat), center = FALSE, scale. = FALSE)
      pca_df <- as.data.frame(pca$x[, 1:2, drop = FALSE])
      pca_df$sample <- rownames(pca_df)
      
      p <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, label = sample)) +
        ggplot2::geom_point(size = 2.5, alpha = 0.9) +
        ggplot2::geom_text(vjust = -0.8, size = 3) +
        ggplot2::theme_bw() +
        ggplot2::labs(title = "PCA of normalized counts (z-scored by gene)",
                      x = "PC1", y = "PC2")
      
      ggplot2::ggsave(out, p, width = 6.5, height = 4.5, dpi = 300)
      out
    },
    format = "file"
  ),
      tar_target(
        counts_pca_outliers_csv,
        {
          dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
          
          mat <- t(scale(t(as.matrix(counts_norm))))
          pca <- prcomp(t(mat), center = FALSE, scale. = FALSE)
          
          df <- as.data.frame(pca$x[, 1:min(5, ncol(pca$x)), drop = FALSE])
          
          mu <- colMeans(df)
          S  <- stats::cov(df)
          eps <- 1e-6 * mean(diag(S))
          Sr  <- S + diag(eps, ncol(S))
          d2 <- stats::mahalanobis(df, mu, Sr)
          
          out <- data.frame(
            sample = rownames(df),
            mahalanobis_d2 = d2,
            row.names = NULL
          )
          
          readr::write_csv(out, "outputs/pca_outliers.csv")
          "outputs/pca_outliers.csv"
        },
        format = "file"
      ), tar_target(
        counts_pca_plot_png,
        {
          dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
          
          # use the normalized matrix from counts_norm
          mat <- t(scale(t(as.matrix(counts_norm))))
          pca <- prcomp(t(mat), center = FALSE, scale. = FALSE)
          
          df <- as.data.frame(pca$x[, 1:2, drop = FALSE])
          df$sample <- rownames(df)
          
          # if you have a metadata table called `meta` with a column "label",
          # join it to color points by label (optional but nice)
          if (exists("meta")) {
            df <- dplyr::left_join(df, meta, by = c("sample" = "SampleID"))
          }
          
          var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
          xlab_txt <- paste0("PC1 (", scales::percent(var_expl[1]), ")")
          ylab_txt <- paste0("PC2 (", scales::percent(var_expl[2]), ")")
          
          p <- ggplot2::ggplot(
            df,
            ggplot2::aes(x = PC1, y = PC2, label = sample,
                         color = if ("label" %in% names(df)) .data$label else NULL)
          ) +
            ggplot2::geom_point(size = 2.4, alpha = 0.9) +
            ggplot2::labs(title = "PCA (normalized counts)", x = xlab_txt, y = ylab_txt, color = "label") +
            ggplot2::theme_bw(base_size = 11)
          
          out_file <- "outputs/pca_scatter.png"
          ggplot2::ggsave(out_file, p, width = 7, height = 5, dpi = 300)
          out_file
        },
        format = "file"
      )
  
  
  
  
)
