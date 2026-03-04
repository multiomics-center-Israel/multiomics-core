# Multi-Omics Integration Domain

This directory contains the multi-omics integration functionality for the `multiomics-core` pipeline. It enables joint analysis of 2+ omics layers (RNA-seq, proteomics, metabolomics) to identify coordinated molecular changes across biological scales.

## Architecture

The multi-omics domain follows the same layered architecture as single-omics modes:

```
R/domain/multiomics/   → Domain logic (integration algorithms)
R/modules/multiomics/  → Target-ready wrappers
R/pipeline/multiomics/ → {targets} orchestration
```

## Domain Files (R/domain/multiomics/)

| File | Purpose |
|------|---------|
| `00_inputs.R` | Load preprocessed data from single-omics pipelines |
| `01_id_mapping.R` | Gene-protein ID mapping via symbols or custom file |
| `02_mae.R` | MultiAssayExperiment construction & feature harmonization |
| `03_feature_selection.R` | Select informative features for integration (DE/variance-based) |
| `04_integration_diablo.R` | DIABLO (mixOmics) supervised integration |
| `05_integration_snf.R` | SNF (Similarity Network Fusion) unsupervised clustering |
| `06_concordance.R` | RNA-protein concordance analysis |
| `07_enrichment.R` | Cross-omics pathway enrichment (Fisher's method) |
| `08_report.R` | Multi-omics HTML report generation |
| `90_config_validate.R` | Config validation for multi-omics mode |

## Module Files (R/modules/multiomics/)

| File | Purpose |
|------|---------|
| `01_mod_harmonization.R` | Harmonize data → build MAE → gene-protein mapping |
| `02_mod_integration.R` | Run DIABLO/SNF integration with feature selection |
| `03_mod_concordance.R` | Concordance analysis across omics |
| `04_mod_enrichment.R` | Cross-omics pathway meta-analysis |
| `05_mod_report.R` | Generate final integration report |

## Pipeline Orchestration (R/pipeline/multiomics/)

`00_pipe_multiomics.R` defines the {targets} pipeline:

```
harmonization → de_results → enrichment_results
    ↓               ↓              ↓
    └──> integration ──────────────┘
         (DIABLO/SNF)
              ↓
         concordance + cross_enrichment
              ↓
           report
```

The pipeline automatically detects which single-omics pipelines have run and integrates only available data.

## Data Flow

### 1. Harmonization
- Loads preprocessed data from `rna_pre`, `prot_pre`, `metab_pre` targets
- Validates sample overlap (requires ≥3 common samples)
- Builds gene-protein mapping (symbol-based or custom file)
- Constructs `MultiAssayExperiment` (MAE) with aligned samples

### 2. Feature Selection
- **DE-based**: Select top differentially expressed features per omics
- **Variance-based**: Select high-variance features (default)
- **Combined**: Union or intersection of DE + variance

### 3. Integration Methods

#### DIABLO (mixOmics)
- **Type**: Supervised (requires condition labels)
- **Method**: Partial Least Squares Discriminant Analysis (PLS-DA) for multi-block data
- **Output**:
  - Sample scores (latent components)
  - Feature loadings (contribution to each component)
  - Circos plot (cross-omics correlations)
  - Cross-validation performance

#### SNF (Similarity Network Fusion)
- **Type**: Unsupervised
- **Method**: Fuses patient similarity networks from each omics
- **Output**:
  - Fused similarity network
  - Spectral clustering assignments
  - Silhouette scores (cluster quality)

### 4. Concordance Analysis
- Compares logFC patterns between RNA and protein
- Classifies features as: concordant, discordant, RNA-only, protein-only
- Computes Pearson correlation of logFC values

### 5. Cross-Omics Enrichment
- Identifies pathways consistently dysregulated across omics
- Combines p-values using Fisher's method
- Generates heatmap of pathway significance across layers

## Configuration

Add a `multiomics` mode to your config YAML:

```yaml
modes:
  # ... rna, proteomics, metabolomics configs ...

  multiomics:
    # Gene-protein mapping
    gene_protein_mapping_file: null  # or "path/to/mapping.csv"
    require_one_to_one_mapping: false

    # Feature selection
    feature_selection:
      method: "variance"  # "de", "variance", "combined"
      top_n: 500

    # Integration methods
    integration:
      methods: ["DIABLO", "SNF"]

      diablo:
        ncomp: 2
        design_matrix: "full"
        cv_folds: 5

      snf:
        K: 20
        alpha: 0.5
        n_clusters: 2

    # Condition for supervised methods
    condition_column: "condition"

    # Cross-omics enrichment
    enrichment:
      run_enrichment: true
```

See `config/templates/multiomics_config.yaml` for a full example.

## Requirements

### R Packages
- **Core**: `MultiAssayExperiment`, `SummarizedExperiment`
- **Integration**: `mixOmics` (DIABLO), `SNFtool` (SNF)
- **Utilities**: `dplyr`, `ggplot2`, `pheatmap`, `cluster`

### Prerequisites
- At least 2 single-omics modes configured (`rna`, `proteomics`, `metabolomics`)
- Individual omics pipelines must run first to generate preprocessed data
- ≥3 common samples across omics layers

## Running the Pipeline

```r
# Set config with multi-omics mode enabled
Sys.setenv(MULTIOMICS_CONFIG = "config/multiomics_example.yaml")

# Run pipeline
library(targets)
tar_make()

# Visualize DAG
tar_visnetwork()

# Read results
mae <- tar_read(multiomics_harmonization)$mae
diablo_res <- tar_read(multiomics_integration)$diablo_results
snf_res <- tar_read(multiomics_integration)$snf_results
```

## Outputs

All results are saved to `outputs/{project_name}/multiomics/`:

```
multiomics/
├── gene_protein_mapping.csv
├── integration/
│   ├── diablo/
│   │   ├── diablo_sample_plot.png
│   │   ├── diablo_circos_plot.png
│   │   ├── diablo_scores_*.csv
│   │   └── diablo_top_features_*.csv
│   └── snf/
│       ├── snf_network_heatmap.png
│       ├── snf_clusters.csv
│       └── snf_silhouette_scores.csv
├── concordance/
│   ├── concordance_rna_protein.png
│   └── concordance_transcriptomics_vs_proteomics.csv
├── cross_enrichment/
│   ├── cross_omics_pathway_heatmap.png
│   └── cross_omics_pathways_meta_analysis.csv
└── report_multiomics.html
```

## Future Enhancements

- **MOFA2 integration**: Deferred due to complex Python/reticulate setup
- **Advanced report template**: R Markdown with detailed method descriptions
- **Stability analysis**: Bootstrap validation of integration results
- **Network analysis**: Cross-omics regulatory networks
- **Executive summary**: AI-powered interpretation of integration results

## Troubleshooting

### "No common samples found"
- Check sample ID column names match across omics metadata
- Verify sample naming is consistent (e.g., underscores vs hyphens)

### "Cannot find logFC columns"
- Ensure single-omics DE analysis completed successfully
- Check that DE results have `logFC` or `log2FoldChange` columns

### "mixOmics package not installed"
- Install with: `BiocManager::install("mixOmics")`

### "SNFtool package not installed"
- Install with: `install.packages("SNFtool")`

## References

1. **DIABLO**: Singh et al. (2019). DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays. *Bioinformatics*.
2. **SNF**: Wang et al. (2014). Similarity network fusion for aggregating data types on a genomic scale. *Nature Methods*.
3. **MultiAssayExperiment**: Ramos et al. (2017). Software for the integration of multi-omics experiments in Bioconductor. *Cancer Research*.
