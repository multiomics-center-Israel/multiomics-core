#' Run SNF (Similarity Network Fusion) integration analysis
#'
#' SNF integrates multiple omics layers by fusing patient similarity networks.
#' It's an unsupervised method that identifies patient clusters.
#'
#' @param mae MultiAssayExperiment object with aligned samples
#' @param config Full config object
#' @param out_dir Output directory for plots
#' @return List with: fused_network, clusters, silhouette, plots
run_snf_integration <- function(mae, config, out_dir = NULL) {

    if (!requireNamespace("SNFtool", quietly = TRUE)) {
        stop(
            "Package 'SNFtool' is required for SNF integration. ",
            "Install with: install.packages('SNFtool')",
            call. = FALSE
        )
    }

    cfg <- config$modes$multiomics$integration$snf %||% list()
    K <- cfg$K %||% 20              # Number of neighbors
    alpha <- cfg$alpha %||% 0.5     # Hyperparameter (0-1)
    T <- cfg$T %||% 20              # Number of iterations
    n_clusters <- cfg$n_clusters %||% 2

    # Cap K at n_samples - 1 (SNFtool requires K < n_samples)
    n_samples <- ncol(mae[[1]])
    if (K >= n_samples) {
        K <- max(2, n_samples - 1)
        message(sprintf("  SNF: K reduced to %d (n_samples = %d)", K, n_samples))
    }

    message("Running SNF integration...")

    # Extract data from MAE (transpose to samples x features)
    omics <- names(mae@ExperimentList)
    data_list <- lapply(omics, function(om) {
        expr_mat <- tryCatch(
            SummarizedExperiment::assay(mae[[om]], "expr"),
            error = function(e) {
                avail <- SummarizedExperiment::assayNames(mae[[om]])
                if (length(avail) > 0) {
                    message("  SNF: '", om, "' has no 'expr' assay, using '", avail[1], "'")
                    SummarizedExperiment::assay(mae[[om]], avail[1])
                } else {
                    stop("No assays available for ", om)
                }
            }
        )
        t(expr_mat)
    })
    names(data_list) <- omics

    # Standardize each omics layer (mean=0, sd=1)
    data_list <- lapply(data_list, function(x) {
        scale(x, center = TRUE, scale = TRUE)
    })

    # Build distance matrices (Euclidean)
    dist_list <- lapply(data_list, function(x) {
        as.matrix(dist(x, method = "euclidean"))
    })

    # Convert distance matrices to affinity matrices (patient similarity networks)
    affinity_list <- lapply(dist_list, function(dist_mat) {
        SNFtool::affinityMatrix(dist_mat, K = K, sigma = alpha)
    })

    # Fuse networks using SNF
    if (length(affinity_list) < 2) {
        stop("SNF requires at least 2 omics layers")
    }

    fused_network <- tryCatch({
        if (length(affinity_list) == 2) {
            SNFtool::SNF(list(affinity_list[[1]], affinity_list[[2]]), K = K, t = T)
        } else if (length(affinity_list) == 3) {
            SNFtool::SNF(list(affinity_list[[1]], affinity_list[[2]], affinity_list[[3]]),
                         K = K, t = T)
        } else {
            # For >3 omics, fuse pairwise iteratively
            fused <- affinity_list[[1]]
            for (i in 2:length(affinity_list)) {
                fused <- SNFtool::SNF(list(fused, affinity_list[[i]]), K = K, t = T)
            }
            fused
        }
    }, error = function(e) {
        stop("SNF fusion failed: ", e$message, call. = FALSE)
    })

    # Spectral clustering on fused network
    clusters <- SNFtool::spectralClustering(fused_network, K = n_clusters)
    names(clusters) <- rownames(data_list[[1]])

    # Compute silhouette score for cluster quality
    fused_dist <- as.dist(1 - fused_network)
    sil <- cluster::silhouette(clusters, fused_dist)
    sil_avg <- mean(sil[, "sil_width"])

    message(sprintf(
        "SNF integration complete: %d clusters, silhouette = %.3f",
        n_clusters, sil_avg
    ))

    # Generate plots
    plots <- list()
    if (!is.null(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

        # Heatmap of fused network
        plots$network_heatmap <- file.path(out_dir, "snf_network_heatmap.png")
        png(plots$network_heatmap, width = 800, height = 800, res = 120)
        tryCatch({
            SNFtool::displayClustersWithHeatmap(fused_network, clusters)
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("Heatmap failed:", e$message), cex = 1.2)
        })
        dev.off()

        # Silhouette plot
        plots$silhouette_plot <- file.path(out_dir, "snf_silhouette.png")
        png(plots$silhouette_plot, width = 800, height = 600, res = 120)
        tryCatch({
            plot(sil, col = clusters, border = NA,
                 main = sprintf("SNF Clusters (avg silhouette = %.3f)", sil_avg))
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("Silhouette plot failed:", e$message), cex = 1.2)
        })
        dev.off()

        message("  SNF plots saved to: ", out_dir)
    }

    list(
        fused_network = fused_network,
        clusters = clusters,
        silhouette = sil,
        silhouette_avg = sil_avg,
        affinity_list = affinity_list,
        plots = plots,
        config = cfg
    )
}


#' Estimate optimal number of clusters for SNF
#'
#' Uses eigengap heuristic from SNFtool to suggest optimal cluster count.
#'
#' @param fused_network Fused similarity network from SNF
#' @param K_range Range of cluster numbers to test
#' @return List with: K_opt, eigengaps
estimate_snf_clusters <- function(fused_network, K_range = 2:10) {

    if (!requireNamespace("SNFtool", quietly = TRUE)) {
        stop("Package 'SNFtool' is required")
    }

    K_opt <- SNFtool::estimateNumberOfClustersGivenGraph(fused_network, K_range)

    list(
        K_opt = K_opt[[1]],  # First estimate
        K_range = K_range
    )
}


#' Write SNF results to CSV
#'
#' @param snf_results Output from run_snf_integration()
#' @param out_dir Output directory
write_snf_results <- function(snf_results, out_dir) {

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Cluster assignments
    clusters_df <- data.frame(
        sample_id = names(snf_results$clusters),
        cluster = snf_results$clusters,
        row.names = NULL,
        stringsAsFactors = FALSE
    )
    write.csv(clusters_df, file.path(out_dir, "snf_clusters.csv"), row.names = FALSE)

    # Silhouette scores
    sil_df <- as.data.frame(snf_results$silhouette[, c("cluster", "sil_width")])
    sil_df$sample_id <- rownames(snf_results$silhouette)
    write.csv(sil_df, file.path(out_dir, "snf_silhouette_scores.csv"), row.names = FALSE)

    # Fused network (adjacency matrix)
    write.csv(snf_results$fused_network, file.path(out_dir, "snf_fused_network.csv"),
              row.names = TRUE)

    message("SNF results written to: ", out_dir)

    invisible(NULL)
}
