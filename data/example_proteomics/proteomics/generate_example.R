#!/usr/bin/env Rscript
# Regenerate synthetic proteomics example data used by the wizard.
# Run this file from the repo root if you ever want to reproduce the fixtures.

set.seed(1234)

out_dir     <- "data/example_proteomics/proteomics"
n_proteins  <- 80
n_de        <- 20
n_contam    <- 3
samples_ctrl <- c("CTRL_01", "CTRL_02", "CTRL_03")
samples_trt  <- c("TRT_01",  "TRT_02",  "TRT_03")
all_samples  <- c(samples_ctrl, samples_trt)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Gene/protein names pulled from real human symbols so pathway enrichment has
# something to latch on (mix of metabolic + proteostasis + housekeeping).
pool <- c(
  "ACTB","GAPDH","HSPA8","HSP90AA1","HSP90AB1","EEF1A1","EEF2","RPS3","RPS6",
  "RPL4","RPL7","ENO1","PKM","LDHA","ALDOA","TPI1","PGK1","PGAM1","PFKP","HK1",
  "CS","IDH1","IDH2","SDHA","SDHB","FH","MDH2","ACO2","OGDH","SUCLA2",
  "ATP5F1A","ATP5F1B","ATP5F1C","NDUFA1","NDUFA2","UQCRC1","UQCRC2","COX4I1",
  "SOD1","SOD2","CAT","GPX1","GPX4","PRDX1","PRDX2","PRDX3","TXN","GSR","GCLC",
  "HMOX1","NRF2","KEAP1","NFE2L2","ATF4","ATF6","XBP1","HSPA5","DDIT3","EIF2AK3",
  "TUBB","TUBA1A","TUBA1B","VIM","KRT8","KRT18","ACTG1","TLN1","VCL","CAPZA1",
  "YWHAB","YWHAE","YWHAG","YWHAH","YWHAQ","YWHAZ","CALM1","CALM2","CALR","CANX"
)
genes <- sample(pool, n_proteins - n_contam, replace = FALSE)
uniprots <- sprintf("P%05d", seq_len(n_proteins - n_contam) + 10000)
prot_names <- paste0(genes, "_HUMAN")

# Build intensity matrix on log2 scale, then convert to linear.
base_log2 <- rnorm(n_proteins - n_contam, mean = 22, sd = 2.5)
mat_log2  <- matrix(rep(base_log2, each = length(all_samples)),
                    nrow = n_proteins - n_contam, byrow = TRUE)
# Sample-level noise
mat_log2 <- mat_log2 + matrix(rnorm(length(mat_log2), 0, 0.25),
                               nrow = nrow(mat_log2))

# Inject DE signal: first n_de proteins up/down in TRT vs CTRL.
de_idx <- seq_len(n_de)
lfc    <- sample(c(-1, 1), n_de, replace = TRUE) * runif(n_de, 1.5, 3.5)
trt_cols <- match(samples_trt, all_samples)
mat_log2[de_idx, trt_cols] <- mat_log2[de_idx, trt_cols] + lfc

# Introduce realistic missingness: ~8% of non-DE values NA.
na_mask <- matrix(runif(length(mat_log2)) < 0.08, nrow = nrow(mat_log2))
na_mask[de_idx, ] <- FALSE
mat_linear <- 2^mat_log2
mat_linear[na_mask] <- NA

# Add contaminants (cRAP-*) with high steady intensities.
contam_rows <- data.frame(
  Protein.Group  = c("cRAP-P00761", "cRAP-P02769", "cRAP-P00883"),
  Protein.Names  = c("cRAP-TRYP_PIG", "cRAP-ALBU_BOVIN", "cRAP-KRT1_HUMAN"),
  Genes          = c("", "", "KRT1"),
  First.Protein.Description = c("Trypsin", "Serum albumin BSA", "Keratin type II")
)
contam_mat <- matrix(runif(n_contam * length(all_samples), 1.2e10, 1.8e10),
                      nrow = n_contam, dimnames = list(NULL, all_samples))
contam_full <- cbind(contam_rows, as.data.frame(contam_mat))

real_rows <- data.frame(
  Protein.Group  = uniprots,
  Protein.Names  = prot_names,
  Genes          = genes,
  First.Protein.Description = paste(genes, "synthetic example protein")
)
real_full <- cbind(real_rows, as.data.frame(mat_linear))
names(real_full)[5:ncol(real_full)] <- all_samples

pg_matrix <- rbind(contam_full, real_full)

write.csv(pg_matrix,
          file.path(out_dir, "protein_matrix.csv"),
          row.names = FALSE, na = "")

# Metadata: SampleID + Condition
metadata <- data.frame(
  SampleID  = all_samples,
  Condition = c(rep("Control", 3), rep("Treatment", 3))
)
write.csv(metadata, file.path(out_dir, "metadata.csv"), row.names = FALSE)

# Sample map: orig -> SampleID (identity here because columns already match)
sample_map <- data.frame(
  sample_orig_name = all_samples,
  SampleID         = all_samples
)
write.csv(sample_map, file.path(out_dir, "sample_map.csv"), row.names = FALSE)

# Contrasts
contrasts <- data.frame(
  Contrast_name = "Treatment_vs_Control",
  Factor        = "Condition",
  Numerator     = "Treatment",
  Denominator   = "Control"
)
write.csv(contrasts, file.path(out_dir, "contrasts.csv"), row.names = FALSE)

cat(sprintf("Wrote %d proteins x %d samples to %s\n",
            nrow(pg_matrix), length(all_samples), out_dir))
