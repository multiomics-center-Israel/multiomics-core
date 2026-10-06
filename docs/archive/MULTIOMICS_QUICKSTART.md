# Multi-Omics Pipeline Quick Start Guide

## Prerequisites

### 1. Install Required R Packages

```r
# Install if not already installed
install.packages(c("renv", "targets", "tarchetypes"))

# Restore project environment
renv::restore()

# Multi-omics specific packages
BiocManager::install(c("MultiAssayExperiment", "SummarizedExperiment"))
install.packages(c("SNFtool"))
BiocManager::install("mixOmics")  # For DIABLO
```

### 2. Prepare Your Config

Create a YAML config with at least 2 omics modes + multiomics section:

```yaml
modes:
  rna: { ... }
  proteomics: { ... }  # or metabolomics

  multiomics:
    feature_selection:
      method: "variance"
      top_n: 500
    integration:
      methods: ["DIABLO", "SNF"]
    condition_column: "condition"
```

See: `config/templates/multiomics_config.yaml` for full example
Test config: `config/multiomics_GT15_test.yaml`

## Running the Pipeline

### Option 1: From R Console

```r
# Set config path
Sys.setenv(MULTIOMICS_CONFIG = "config/multiomics_GT15_test.yaml")

# Load targets
library(targets)

# Visualize pipeline DAG (optional)
tar_visnetwork()

# Run the full pipeline
tar_make()

# Or run specific targets
tar_make(multiomics_harmonization)
tar_make(multiomics_integration)
```

### Option 2: From Command Line

```bash
# Set config
export MULTIOMICS_CONFIG="config/multiomics_GT15_test.yaml"

# Run pipeline
Rscript -e "library(targets); tar_make()"

# Check status
Rscript -e "library(targets); tar_progress()"
```

### Option 3: Using run.R Script

```bash
Rscript run.R --config config/multiomics_GT15_test.yaml
```

## Pipeline Execution Flow

```
┌─────────────────────────────────────────────────┐
│ Single-Omics Pipelines (run first)             │
├─────────────────────────────────────────────────┤
│ rna_pre       → rna_de_res       → rna_pathway │
│ prot_pre      → prot_de_res      → prot_pathway│
│ metab_pre     → metab_de_res     → metab_enrich│
└─────────────┬───────────────────────────────────┘
              │
              ↓
┌─────────────────────────────────────────────────┐
│ Multi-Omics Integration                         │
├─────────────────────────────────────────────────┤
│ 1. multiomics_harmonization                     │
│    - Load preprocessed data                     │
│    - Build gene-protein mapping                 │
│    - Create MultiAssayExperiment (MAE)          │
│                                                  │
│ 2. multiomics_integration                       │
│    - Feature selection (top 500 per omics)      │
│    - DIABLO: supervised PLS-DA                  │
│    - SNF: unsupervised clustering               │
│                                                  │
│ 3. multiomics_concordance                       │
│    - RNA-protein logFC correlation              │
│    - Classify concordant/discordant changes     │
│                                                  │
│ 4. multiomics_cross_enrichment                  │
│    - Fisher's method for pathway meta-analysis  │
│    - Identify consistent pathway changes        │
│                                                  │
│ 5. multiomics_report                            │
│    - Generate HTML report                       │
└─────────────────────────────────────────────────┘
```

## Accessing Results

```r
# Read targets after pipeline completes
library(targets)

# Harmonization results
harm <- tar_read(multiomics_harmonization)
mae <- harm$mae                           # MultiAssayExperiment
mapping <- harm$gene_protein_mapping      # Gene-protein pairs

# Integration results
integ <- tar_read(multiomics_integration)
diablo <- integ$diablo_results            # DIABLO model + plots
snf <- integ$snf_results                  # SNF clusters + network

# Concordance
conc <- tar_read(multiomics_concordance)
conc$concordance_table                    # RNA-protein agreement

# Cross-omics enrichment
enrich <- tar_read(multiomics_cross_enrichment)
enrich$meta_analysis                      # Combined pathway p-values
```

## Output Directory Structure

```
outputs/{project_name}/multiomics/
├── gene_protein_mapping.csv
├── integration/
│   ├── diablo/
│   │   ├── diablo_sample_plot.png        # PCA-like scores
│   │   ├── diablo_circos_plot.png        # Cross-omics correlations
│   │   ├── diablo_variable_plot.png      # Feature loadings
│   │   ├── diablo_scores_*.csv
│   │   └── diablo_top_features_*.csv
│   └── snf/
│       ├── snf_network_heatmap.png       # Fused similarity network
│       ├── snf_silhouette.png            # Cluster quality
│       ├── snf_clusters.csv
│       └── snf_fused_network.csv
├── concordance/
│   ├── concordance_rna_protein.png       # Scatter plot
│   └── concordance_*.csv
├── cross_enrichment/
│   ├── cross_omics_pathway_heatmap.png   # Pathway significance
│   └── cross_omics_pathways_meta_analysis.csv
└── report_multiomics.html                # Final report
```

## Troubleshooting

### "No common samples found"
- Check that sample IDs match exactly across omics metadata
- Use `sample_col` to specify the correct ID column
- Try `sample_alignment: lenient` in RNA config

### "Cannot find logFC columns"
- Ensure single-omics DE analysis completed successfully
- Check that DE results have `logFC` or `log2FoldChange` columns

### "mixOmics package not installed"
```r
BiocManager::install("mixOmics")
```

### "SNFtool package not installed"
```r
install.packages("SNFtool")
```

### Pipeline stuck or slow
```r
# Check what's running
tar_progress()

# View outdated targets
tar_outdated()

# Invalidate specific target to re-run
tar_invalidate(multiomics_integration)
```

## Customization

### Use only SNF (skip DIABLO if samples don't match exactly)

```yaml
multiomics:
  integration:
    methods: ["SNF"]  # Remove DIABLO
```

### Adjust feature selection

```yaml
multiomics:
  feature_selection:
    method: "de"      # Use DE features instead of variance
    top_n: 1000       # More features
```

### Custom gene-protein mapping

```yaml
multiomics:
  gene_protein_mapping_file: "data/my_gene_protein_map.csv"
```

File format:
```csv
gene_id,protein_id,mapping_source
GENE001,PROT001,custom
GENE002,PROT002,custom
```

## Next Steps

1. Review the integration results in the HTML report
2. Examine top features from DIABLO loadings
3. Validate concordant features in the lab
4. Explore pathway enrichment for biological insights
5. Consider stability analysis (bootstrap validation)

## Support

- Documentation: `R/domain/multiomics/README.md`
- Config templates: `config/templates/multiomics_config.yaml`
- Issues: Check function comments in R files for details
