#!/usr/bin/env Rscript
# Regenerate synthetic metabolomics example data used by the wizard.
# Run from the repo root: Rscript data/example_metabolomics/metabolomics/generate_example.R

set.seed(1234)

out_dir       <- "data/example_metabolomics/metabolomics"
n_metabolites <- 150
n_de          <- 25
samples_ctrl  <- c("CTRL_01", "CTRL_02", "CTRL_03")
samples_trt   <- c("TRT_01",  "TRT_02",  "TRT_03")
all_samples   <- c(samples_ctrl, samples_trt)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Metabolite names pulled from real KEGG / HMDB compounds so enrichment hits
# central carbon, lipid, amino acid, nucleotide metabolism.
pool <- c(
    "Glucose", "Glucose-6-phosphate", "Fructose-6-phosphate", "Fructose-1,6-bisphosphate",
    "Dihydroxyacetone phosphate", "Glyceraldehyde-3-phosphate", "3-Phosphoglycerate",
    "Phosphoenolpyruvate", "Pyruvate", "Lactate", "Acetyl-CoA", "Citrate", "Isocitrate",
    "alpha-Ketoglutarate", "Succinyl-CoA", "Succinate", "Fumarate", "Malate", "Oxaloacetate",
    "Glutamine", "Glutamate", "Aspartate", "Asparagine", "Alanine", "Serine", "Glycine",
    "Threonine", "Cysteine", "Methionine", "Valine", "Leucine", "Isoleucine", "Lysine",
    "Arginine", "Histidine", "Phenylalanine", "Tyrosine", "Tryptophan", "Proline",
    "ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "CTP", "CDP", "CMP", "UTP", "UDP", "UMP",
    "NAD+", "NADH", "NADP+", "NADPH", "FAD", "FADH2", "Coenzyme A",
    "S-adenosylmethionine", "S-adenosylhomocysteine", "Homocysteine", "Cystathionine",
    "Choline", "Phosphocholine", "Glycerophosphocholine", "Betaine", "Sarcosine",
    "Creatine", "Creatinine", "Urea", "Uric acid", "Hypoxanthine", "Xanthine", "Inosine",
    "Adenosine", "Guanosine", "Cytidine", "Uridine", "Thymidine",
    "Palmitate", "Stearate", "Oleate", "Linoleate", "Arachidonate", "DHA", "EPA",
    "Glycerol", "Glycerol-3-phosphate", "Sphinganine", "Sphingosine", "Ceramide",
    "Cholesterol", "Cholesteryl ester", "Bile acid",
    "Riboflavin", "FMN", "Pyridoxal", "Pyridoxine", "Thiamine", "Thiamine pyrophosphate",
    "Pantothenate", "Biotin", "Folate", "Tetrahydrofolate", "Niacin", "Nicotinamide",
    "Ascorbate", "Dehydroascorbate", "Tocopherol", "Retinal", "Retinol",
    "Carnitine", "Acetylcarnitine", "Propionylcarnitine", "Octanoylcarnitine",
    "Putrescine", "Spermidine", "Spermine", "Ornithine", "Citrulline",
    "GABA", "Dopamine", "Serotonin", "Histamine", "Tyramine", "Adrenaline",
    "Acetylcholine", "Cortisol", "Estradiol", "Testosterone", "Progesterone",
    "Glutathione", "Glutathione disulfide", "Taurine", "Hypotaurine",
    "Malonyl-CoA", "Succinyl-CoA", "Methylmalonyl-CoA", "Propionyl-CoA",
    "Hydroxyproline", "Beta-alanine", "Aminoisobutyrate", "Aminoadipate",
    "Indole", "Indole-3-acetate", "Kynurenine", "Anthranilate", "Quinolinate",
    "Cadaverine", "Trimethylamine", "TMAO", "Choline sulfate", "Pyroglutamate"
)
metab_names <- sample(unique(pool), n_metabolites, replace = FALSE)
feature_ids <- sprintf("METAB_%03d", seq_len(n_metabolites))

# m/z and RT plausible for HILIC + RP combined
mz <- round(runif(n_metabolites, 80, 900), 4)
rt <- round(runif(n_metabolites, 0.3, 18), 2)

# Intensities on log2 scale (typical UHPLC-MS ranges)
base_log2 <- rnorm(n_metabolites, mean = 18, sd = 2.5)
mat_log2 <- matrix(rep(base_log2, each = length(all_samples)),
                   nrow = n_metabolites, byrow = TRUE)
mat_log2 <- mat_log2 + matrix(rnorm(length(mat_log2), 0, 0.3),
                              nrow = nrow(mat_log2))

# Inject DE signal: first n_de metabolites up/down in TRT vs CTRL
de_idx <- seq_len(n_de)
lfc    <- sample(c(-1, 1), n_de, replace = TRUE) * runif(n_de, 1.2, 3.0)
trt_cols <- match(samples_trt, all_samples)
mat_log2[de_idx, trt_cols] <- mat_log2[de_idx, trt_cols] + lfc

# Realistic ~5% missingness on non-DE metabolites
na_mask <- matrix(runif(length(mat_log2)) < 0.05, nrow = nrow(mat_log2))
na_mask[de_idx, ] <- FALSE
mat_linear <- 2 ^ mat_log2
mat_linear[na_mask] <- NA

# Write data file (processed_wide format)
data_df <- data.frame(
    feature_id = feature_ids,
    Name       = metab_names,
    `m/z`      = mz,
    `RT [min]` = rt,
    mat_linear,
    check.names = FALSE,
    stringsAsFactors = FALSE
)
colnames(data_df)[5:ncol(data_df)] <- all_samples
write.csv(data_df, file.path(out_dir, "data.csv"), row.names = FALSE)

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

cat(sprintf("Wrote %d metabolites x %d samples to %s\n",
            n_metabolites, length(all_samples), out_dir))
cat(sprintf("  DE metabolites: %d (rows 1..%d)\n", n_de, n_de))
