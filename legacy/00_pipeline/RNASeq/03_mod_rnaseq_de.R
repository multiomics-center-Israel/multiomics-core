# 03_mod_rnaseq_de
# .build_contrasts_matrix <- function(contrasts_df, design, factor_col) {
#   # contrasts_df: Contrast_name, Factor, Numerator, Denominator
#   stopifnot(all(c("Contrast_name", "Factor", "Numerator", "Denominator") %in% names(contrasts_df)))
#   
#   contrasts_df <- contrasts_df[contrasts_df$Factor == factor_col, , drop = FALSE]
#   if (nrow(contrasts_df) == 0) stop("No contrasts match Factor='", factor_col, "'.")
#   
#   # design columns look like: factor_colLevel (because ~0 + factor)
#   levs <- sub(paste0("^", factor_col), "", colnames(design))
#   
#   make_one <- function(num, den) {
#     num_col <- paste0(factor_col, num)
#     den_col <- paste0(factor_col, den)
#     if (!num_col %in% colnames(design)) stop("Numerator level not in design: ", num)
#     if (!den_col %in% colnames(design)) stop("Denominator level not in design: ", den)
#     v <- rep(0, ncol(design)); names(v) <- colnames(design)
#     v[num_col] <- 1
#     v[den_col] <- -1
#     v
#   }
#   
#   cm <- vapply(seq_len(nrow(contrasts_df)), function(i) {
#     make_one(contrasts_df$Numerator[i], contrasts_df$Denominator[i])
#   }, FUN.VALUE = numeric(ncol(design)))
#   
#   cm <- t(cm)
#   rownames(cm) <- contrasts_df$Contrast_name
#   cm
# }

# mod_rnaseq_de
mod_rnaseq_de <- function(pre, inputs, config, verbose = FALSE) {
  assert_pre_contract(pre, stage = "rna")
  
  cfg <- config$modes$rna
  contrasts_df <- inputs$contrasts
  
  counts <- pre$expr_filt
  meta   <- pre$meta
  storage.mode(counts) <- "integer"
  
  factor_col <- unique(contrasts_df$Factor)
  if (length(factor_col) != 1) stop("RNA contrasts must use a single Factor.")
  factor_col <- factor_col[[1]]
  if (!factor_col %in% colnames(meta)) stop("meta missing Factor column: ", factor_col)
  
  meta[[factor_col]] <- as.factor(meta[[factor_col]])
  
  dds0 <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData   = as.data.frame(meta),
    design    = stats::as.formula(paste0("~ ", factor_col))
  )
  
  # run DESeq
  dds <- DESeq2::DESeq(dds0)
  
  # tables per contrast (keep minimal cols to reduce size)
  tables <- list()
  for (i in seq_len(nrow(contrasts_df))) {
    cn  <- contrasts_df$Contrast_name[i]
    num <- contrasts_df$Numerator[i]
    den <- contrasts_df$Denominator[i]
    
    res <- DESeq2::results(
      dds,
      contrast = c(factor_col, num, den),
      independentFiltering = FALSE,
      cooksCutoff = FALSE
    )
    
    tab <- as.data.frame(res)
    tab$FeatureID <- rownames(tab)
    tab$Contrast  <- cn
    
    tab <- tab[, intersect(c("FeatureID","Contrast","log2FoldChange","pvalue","padj"), colnames(tab)), drop = FALSE]
    tables[[cn]] <- tab
  }
  
  list(
    method    = "DESeq2",
    tables    = tables,
    contrasts = contrasts_df,
    info      = list(
      design    = paste0("~", factor_col),
      n_genes   = nrow(dds),
      n_samples = ncol(dds)
    )
  )
}
