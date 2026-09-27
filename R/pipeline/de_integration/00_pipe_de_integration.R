# R/pipeline/de_integration/00_pipe_de_integration.R

#' DE-integration pipeline
#'
#' Compares finished DE tables across layers (proteomics or RNA-seq, from this
#' pipeline or elsewhere) without the sample-level matrices the multiomics mode
#' needs. It runs on its own, from \code{modes$de_integration} alone.
#'
#' @return List of target objects.
pipe_de_integration <- function() {
    list(
        tar_target(
            dei_out_dir,
            {
                d <- get_mode_out_dir(run_dir, "de_integration")
                ensure_dir(d)
                d
            },
            format = "file"
        ),

        # Every table and side file the layers name, so editing one reruns what
        # was built from it.
        tar_target(
            dei_input_files,
            de_integration_input_files(config),
            format = "file"
        ),

        tar_target(
            dei_layers,
            {
                dei_input_files
                mod_dei_load_layers(config)
            }
        ),

        tar_target(
            dei_layer_summary_file,
            write_dei_layer_summary(dei_layers, dei_out_dir),
            format = "file"
        )
    )
}
