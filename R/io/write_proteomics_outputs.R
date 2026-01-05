write_proteomics_multimpute_outputs <- function(pre, de_res, inputs, config, run_dir, write_runs = FALSE) {
  
  files <- character(0)
  dirs  <- create_legacy_output_dirs(run_dir)
  
  # 1) datasets
  runs_for_datasets <- if (isTRUE(write_runs)) de_res$runs else NULL
  files <- c(files, write_proteomics_datasets_legacy(pre, runs_for_datasets, config, run_dir))
  
  # 2) summary
  if (!is.null(de_res$summary_df)) {
    files <- c(files, write_limma_multimp_summary_legacy(de_res$summary_df, config, run_dir))
  }
  
  # 3) wide limma per contrast (optional)
  if (!is.null(de_res$runs_de_tables) && length(de_res$runs_de_tables) > 0) {
    contrast_names <- names(de_res$runs_de_tables[[1]])
    for (cn in contrast_names) {
      files <- c(files, write_limma_results_multimp_legacy(de_res = de_res, contrast_name = cn, config = config, run_dir = run_dir))
    }
  }
  
  # 4) final results TSV (optional but useful)
  if (!is.null(inputs$contrasts) && !is.null(de_res$summary_df)) {
    final_results <- build_final_results_proteomics(
      pre          = pre,
      summary_df   = de_res$summary_df,
      contrasts_df = inputs$contrasts,
      row_data     = pre$row_data
    )
    files <- c(files, write_tsv(final_results, dirs$datasets, "final_results.tsv"))
    
    # 5) Excel outputs (delegated to dedicated file)
    files <- c(files, write_final_results_excels_legacy(final_results, pre, config, run_dir))
  }
  
  unique(files)
}

#' Write cluster data in exact legacy format
#' Columns: Name (Sample), Group, Exp (Absolute Expression)
#' Summary File Columns: Cluster, Group, Mean, SE, Mean_SE.y, Mean_SE.ymin, Mean_SE.ymax
#'
#' @param expr_mat Original expression matrix (log-intensities), NOT z-scores
#' @param meta Metadata dataframe
#' @param clusters Named integer vector of clusters
#' @param cfg Config list
#' @param out_dir Output directory
write_legacy_cluster_data <- function(expr_mat, meta, clusters, cfg, out_dir) {
  
  require(dplyr)
  require(tidyr)
  require(tibble)
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # 1. Prepare Metadata Map (Sample -> Group)
  group_col  <- cfg$effects$color
  sample_col <- cfg$effects$samples
  
  meta_map <- meta %>% 
    dplyr::select(Name = all_of(sample_col), Group = all_of(group_col))
  
  # 2. Convert Expression Matrix to Long Format
  # Rows = Genes, Cols = Samples -> Melt
  df_long <- as.data.frame(expr_mat) %>%
    tibble::rownames_to_column("Gene") %>%
    tidyr::pivot_longer(cols = -Gene, names_to = "Name", values_to = "Exp")
  
  # 3. Join with Metadata
  df_annotated <- df_long %>%
    dplyr::inner_join(meta_map, by = "Name")
  
  # 4. Map Genes to Clusters
  cluster_map <- data.frame(
    Gene = names(clusters), 
    Cluster = as.integer(clusters), 
    stringsAsFactors = FALSE
  )
  
  # Final Join: Only keep genes that are in a cluster
  df_final <- df_annotated %>%
    dplyr::inner_join(cluster_map, by = "Gene")
  
  files_written <- character(0)
  
  # 5. Write Per-Cluster Files (Raw Data)
  # Keeping raw data as is (or rounding if you want, usually raw is kept precise)
  unique_clusters <- sort(unique(df_final$Cluster))
  
  for (k in unique_clusters) {
    clus_data <- df_final %>% 
      dplyr::filter(Cluster == k) %>%
      dplyr::select(Name, Group, Exp) 
    
    fname <- file.path(out_dir, sprintf("cluster_profiles_cluster%s_data.txt", k))
    
    write.table(clus_data, fname, sep = "\t", quote = FALSE, row.names = FALSE)
    files_written <- c(files_written, fname)
  }
  
  # 6. Write Summary File (Calculated Stats)
  fname_all <- file.path(out_dir, "cluster_profiles_data.txt")
  
  summary_df <- df_final %>%
    dplyr::group_by(Cluster, Group) %>%
    dplyr::summarise(
      Mean = mean(Exp, na.rm = TRUE),
      SE   = sd(Exp, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Mean_SE.y    = Mean,
      Mean_SE.ymin = Mean - SE,
      Mean_SE.ymax = Mean + SE
    ) %>%
    # --- Rounding to 4 decimal places ---
    dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
  
  write.table(summary_df, fname_all, sep = "\t", quote = FALSE, row.names = FALSE)
  files_written <- c(files_written, fname_all)
  
  return(files_written)
}