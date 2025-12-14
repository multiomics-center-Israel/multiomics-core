#' Normalize RNA-seq counts
#'
#' Options:
#' - "TMMlogCPM" – edgeR TMM normalization + logCPM
#' - "VST"       – DESeq2 variance stabilizing transform
#'
#' @param counts matrix/data.frame of raw counts (genes x samples)
#' @param meta   data.frame of sample metadata (required for VST);
#'               rownames(meta) must match colnames(counts)
#' @param method one of c("TMMlogCPM","VST")
#' @param prior.count numeric added before log in logCPM (default 1)
#' @return numeric matrix; attr(., "method") indicates method used
normalize_counts <- function(counts,
                             meta = NULL,
                             method = c("TMMlogCPM", "VST"),
                             prior.count = 1) {
  stopifnot(is.matrix(counts) || is.data.frame(counts))
  counts <- as.matrix(counts)
  
  # match and normalize legacy alias
  method <- match.arg(method)
 
  if (method == "TMMlogCPM") {
    if (!requireNamespace("edgeR", quietly = TRUE))
      stop("Package 'edgeR' is required. Install via BiocManager::install('edgeR').")
    
    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    mat <- edgeR::cpm(dge, log = TRUE, prior.count = prior.count)
    attr(mat, "method") <- "TMMlogCPM"
    return(mat)
  }
  
  # VST branch
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("Package 'DESeq2' is required. Install via BiocManager::install('DESeq2').")
  
  if (is.null(meta))
    stop("For method = 'VST', 'meta' must be provided.")
  
  if (is.null(rownames(meta)) || !identical(colnames(counts), rownames(meta))) {
    stop("rownames(meta) must exactly match colnames(counts) in order and identity for VST.")
  }
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData   = as.data.frame(meta),
    design    = ~ 1
  )
  keep_nonzero <- rowSums(DESeq2::counts(dds)) > 0
  dds <- dds[keep_nonzero, , drop = FALSE]
  if (nrow(dds) == 0)
    stop("No rows with nonzero counts available for VST.")
  
  mat <- tryCatch({
    # Small matrices: use varianceStabilizingTransformation directly
    if (nrow(dds) < 50) {
      vt <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE, fitType = "mean")
      SummarizedExperiment::assay(vt)
    } else {
      # Larger: use vst with sub-sampling guard
      nsub <- min(1000L, nrow(dds))
      vt <- DESeq2::vst(dds, blind = TRUE, nsub = nsub, fitType = "mean")
      SummarizedExperiment::assay(vt)
    }
  }, error = function(e1) {
    message("[VST] Fallback to TMMlogCPM due to error: ", conditionMessage(e1))
    if (!requireNamespace("edgeR", quietly = TRUE))
      stop("Fallback requires 'edgeR'. Install via BiocManager::install('edgeR').")
    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    fb  <- edgeR::cpm(dge, log = TRUE, prior.count = prior.count)
    attr(fb, "method") <- "TMMlogCPM_fallback_from_VST"
    fb
  })
  
  # If we actually got VST values (not fallback), tag accordingly
  if (is.null(attr(mat, "method"))) attr(mat, "method") <- "VST"
  return(mat)
}
