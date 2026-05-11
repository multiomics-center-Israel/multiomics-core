#!/usr/bin/env Rscript
# Regenerate synthetic lipidomics example data used by the wizard.
# Run from the repo root: Rscript data/example_lipidomics/lipidomics/generate_example.R

set.seed(1234)

out_dir       <- "data/example_lipidomics/lipidomics"
n_lipids      <- 120
n_de          <- 20
samples_ctrl  <- c("CTRL_01", "CTRL_02", "CTRL_03", "CTRL_04")
samples_trt   <- c("TRT_01",  "TRT_02",  "TRT_03",  "TRT_04")
all_samples   <- c(samples_ctrl, samples_trt)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Lipid species: class prefix + acyl-chain notation. Composition reflects
# typical plasma lipidome distribution (phospholipids dominant, plus SM,
# Cer, DG, TG, CE, LPC).
make_species <- function(class, total_C, db, n = 1) {
    if (class %in% c("TG", "DG")) {
        sprintf("%s %d:%d", class, total_C, db)
    } else {
        sprintf("%s %d:%d", class, total_C, db)
    }
}

species <- c(
    paste0("PC ",  c("32:0","32:1","32:2","34:1","34:2","34:3","34:4","36:1","36:2","36:3",
                      "36:4","36:5","38:3","38:4","38:5","38:6","40:5","40:6","40:7")),
    paste0("PE ",  c("34:1","34:2","36:1","36:2","36:3","36:4","38:4","38:5","38:6","40:6","40:7")),
    paste0("PS ",  c("36:1","36:2","38:4","40:6")),
    paste0("PI ",  c("36:2","36:3","38:3","38:4","38:5")),
    paste0("PG ",  c("34:1","34:2","36:2")),
    paste0("PA ",  c("36:2","38:4")),
    paste0("LPC ", c("16:0","16:1","18:0","18:1","18:2","20:4","22:6")),
    paste0("LPE ", c("16:0","18:0","18:1","18:2","20:4")),
    paste0("SM ",  c("d18:1/16:0","d18:1/18:0","d18:1/22:0","d18:1/24:0","d18:1/24:1",
                      "d18:1/16:1","d18:1/18:1")),
    paste0("Cer ", c("d18:1/16:0","d18:1/18:0","d18:1/22:0","d18:1/24:0","d18:1/24:1")),
    paste0("HexCer ", c("d18:1/16:0","d18:1/22:0","d18:1/24:1")),
    paste0("DG ",  c("32:0","34:1","36:2","36:3","38:4")),
    paste0("TG ",  c("48:1","50:2","50:3","52:2","52:3","52:4","54:3","54:4","54:5","54:6",
                      "56:4","56:5","56:6","58:6","58:7")),
    paste0("CE ",  c("16:0","18:0","18:1","18:2","20:4","22:6")),
    paste0("FA ",  c("16:0","18:0","18:1","18:2","20:4","22:6")),
    paste0("MG ",  c("16:0","18:1"))
)
# Ensure we have enough
if (length(species) < n_lipids) {
    # Pad with extra TG variants
    extra <- paste0("TG ", sprintf("%d:%d", sample(48:58, n_lipids - length(species), TRUE),
                                            sample(0:8, n_lipids - length(species), TRUE)))
    species <- c(species, extra)
}
species      <- species[seq_len(n_lipids)]
feature_ids  <- sprintf("LIPID_%03d", seq_len(n_lipids))

# Intensities on log2 scale (typical LipidSearch ranges)
base_log2 <- rnorm(n_lipids, mean = 20, sd = 2.2)
mat_log2  <- matrix(rep(base_log2, each = length(all_samples)),
                    nrow = n_lipids, byrow = TRUE)
mat_log2  <- mat_log2 + matrix(rnorm(length(mat_log2), 0, 0.25),
                                nrow = nrow(mat_log2))

# Inject DE signal — TG and Cer often elevated in metabolic perturbations
de_idx <- seq_len(n_de)
lfc    <- sample(c(-1, 1), n_de, replace = TRUE) * runif(n_de, 1.2, 3.0)
trt_cols <- match(samples_trt, all_samples)
mat_log2[de_idx, trt_cols] <- mat_log2[de_idx, trt_cols] + lfc

# Light missingness (~3%)
na_mask <- matrix(runif(length(mat_log2)) < 0.03, nrow = nrow(mat_log2))
na_mask[de_idx, ] <- FALSE
mat_linear <- 2 ^ mat_log2
mat_linear[na_mask] <- NA

# Write data file — translated_csv format with a group-label row
header <- c("Metabolite", all_samples)
group_row <- c("Label", rep("Control", length(samples_ctrl)), rep("Treatment", length(samples_trt)))
data_rows <- cbind(species, mat_linear)
colnames(data_rows) <- header

data_lines <- c(
    paste(header, collapse = ","),
    paste(group_row, collapse = ","),
    apply(data_rows, 1, function(r) paste(r, collapse = ","))
)
writeLines(data_lines, file.path(out_dir, "data.csv"))

# Metadata
meta_df <- data.frame(
    sample_id = all_samples,
    condition = c(rep("Control", length(samples_ctrl)),
                  rep("Treatment", length(samples_trt))),
    stringsAsFactors = FALSE
)
write.csv(meta_df, file.path(out_dir, "metadata.csv"), row.names = FALSE)

# Contrasts
contrasts_df <- data.frame(
    Contrast_name = "Treatment_vs_Control",
    Factor        = "condition",
    Numerator     = "Treatment",
    Denominator   = "Control",
    stringsAsFactors = FALSE
)
write.csv(contrasts_df, file.path(out_dir, "contrasts.csv"), row.names = FALSE)

cat(sprintf("Wrote %d lipid species x %d samples to %s\n",
            n_lipids, length(all_samples), out_dir))
cat(sprintf("  DE lipids: %d (rows 1..%d)\n", n_de, n_de))
cat("  Format: translated_csv with group-label row 1 (has_group_row=true)\n")
