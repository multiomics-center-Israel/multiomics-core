# RNA-seq QC Configuration Guide

This guide explains how to configure QC plotting for your RNA-seq pipeline.

## Quick Start

The pipeline generates QC plots based on your configuration in `config.yaml`. The key sections are:

1. **effects**: Controls plot aesthetics and annotations
2. **qc**: Controls which plots are generated based on sample counts

---

## Effects Configuration

### Single Color (Default)

For simple experiments with one grouping variable:

```yaml
effects:
  color: "treatment"
  shape: "batch"
  samples: "SampleID"
```

**Generates:**
- PCA plots colored by `treatment`
- Heatmaps annotated with `treatment` column

---

### Multiple Colors (Multi-PCA)

For experiments where you want to visualize multiple grouping variables:

```yaml
effects:
  color:
    - "cell_line"
    - "time_point"
    - "treatment"
  shape: "batch"
  samples: "SampleID"
```

**Generates:**
- `PCA_PC1.vs.PC2.png` (colored by `cell_line`)
- `PCA_PC1.vs.PC2_by_time_point.png` (colored by `time_point`)
- `PCA_PC1.vs.PC2_by_treatment.png` (colored by `treatment`)
- Same for PC1 vs PC3 (if `generate_pca_1v3: true`)

**Important:**
- If `color` is an array, `shape` must be a single value
- `color` and `shape` cannot contain the same column

---

### Custom Heatmap Annotations

By default, heatmaps use the primary color variable for annotations. To add multiple annotation columns:

```yaml
effects:
  color: "treatment"
  shape: "batch"
  samples: "SampleID"

  heatmap_annotations:
    - "treatment"
    - "batch"
    - "time_point"
    - "replicate"
```

**Result:** Distance, correlation, and expression heatmaps will show all 4 annotation columns as colored bars.

---

## QC Plot Configuration

### Adaptive Plotting (Threshold-Based)

Control which plots are generated based on sample counts:

```yaml
qc:
  adaptive_plots: true  # Default
  generate_pca_1v3: false  # Generate PC1 vs PC3 plots

  thresholds:
    min_samples_for_pca3d: 10         # Skip 3D PCA if < 10 samples
    max_samples_for_heatmaps: 120     # Skip distance/corr heatmaps if > 120 samples
    max_samples_for_expr_heatmap: 60  # Skip expression heatmap if > 60 samples
```

**Example behaviors:**

| Samples | 3D PCA | Distance/Corr Heatmaps | Expr Heatmap |
|---------|--------|------------------------|--------------|
| 8       | ❌ Skip | ✅ Generate | ✅ Generate |
| 50      | ✅ Generate | ✅ Generate | ✅ Generate |
| 100     | ✅ Generate | ✅ Generate | ❌ Skip |
| 200     | ✅ Generate | ❌ Skip | ❌ Skip |

---

### Disable Adaptive Plotting

To always generate all plots regardless of sample count:

```yaml
qc:
  adaptive_plots: false
```

**Use case:** You want full QC output even for very large or very small datasets.

---

## Validation Rules

The pipeline validates your configuration and will error if:

1. **Both color AND shape are arrays:**
   ```yaml
   # ❌ INVALID
   effects:
     color: ["treatment", "batch"]
     shape: ["time_point", "replicate"]
   ```

2. **Color and shape overlap:**
   ```yaml
   # ❌ INVALID
   effects:
     color: "treatment"
     shape: "treatment"  # Same column!
   ```

---

## Complete Examples

### Example 1: Simple 2-Group Comparison

```yaml
effects:
  color: "condition"
  shape: null
  samples: "SampleID"

qc:
  adaptive_plots: true
  generate_pca_1v3: false
```

**Output:**
- PCA PC1 vs PC2 (colored by condition)
- 3D PCA (if ≥10 samples)
- Density, boxplot
- Heatmaps (if ≤120 samples)

---

### Example 2: Multi-Factor Time Course

```yaml
effects:
  color:
    - "time_point"
    - "treatment"
    - "cell_line"
  shape: "batch"
  samples: "SampleID"

  heatmap_annotations:
    - "time_point"
    - "treatment"
    - "cell_line"
    - "batch"

qc:
  adaptive_plots: true
  generate_pca_1v3: true  # Want to see PC1 vs PC3 too

  thresholds:
    max_samples_for_heatmaps: 150  # Higher threshold
```

**Output:**
- 3 PCA PC1 vs PC2 plots (one per color variable)
- 3 PCA PC1 vs PC3 plots (one per color variable)
- 3D PCA
- Heatmaps with 4 annotation columns (if ≤150 samples)

---

### Example 3: Large Dataset (200 samples)

```yaml
effects:
  color: "tissue"
  samples: "SampleID"

qc:
  adaptive_plots: true

  thresholds:
    max_samples_for_heatmaps: 120  # Default
    max_samples_for_expr_heatmap: 60  # Default
```

**Output:**
- PCA plots (generated)
- 3D PCA (generated)
- Density, boxplot (generated)
- Distance/correlation heatmaps (❌ skipped - 200 > 120)
- Expression heatmap (❌ skipped - 200 > 60)

**Console messages:**
```
Skipping distance heatmap (200 samples > 120 threshold)
Skipping correlation heatmap (200 samples > 120 threshold)
Skipping expression heatmap (200 samples > 60 threshold)
```

---

## Files Modified

This configuration system was implemented with changes to:

| File | What Changed |
|------|--------------|
| [R/modules/rnaseq/01_mod_qc_pre.R](R/modules/rnaseq/01_mod_qc_pre.R) | Multi-color PCA generation, config-based annotations, validation |
| [R/domain/rnaseq/03_preprocess.R](R/domain/rnaseq/03_preprocess.R) | Removed profiler integration |
| [R/core/07_qc.R](R/core/07_qc.R) | Multi-column heatmap annotations (unchanged) |
| [config/templates/rna_config.yaml](config/templates/rna_config.yaml) | Updated effects and qc sections |

---

## Troubleshooting

### Config validation errors

**Error:** "Both 'color' and 'shape' cannot be arrays"

**Fix:** Choose only one to be an array (usually `color`)

---

**Error:** "'color' and 'shape' cannot overlap"

**Fix:** Use different columns for color and shape

---

### Plots not appearing

If expected plots are missing, check:

1. **Adaptive plotting is enabled** - check console for skip messages
2. **Sample count exceeds threshold** - adjust thresholds in config
3. **To force all plots:** set `adaptive_plots: false`

---

### Heatmap annotations not showing

If heatmaps only show one annotation column:

1. **Add `heatmap_annotations`** in config
2. **Ensure columns exist** in your metadata
3. **Check column names match** (case-sensitive)

---

## Migration from Profiler

If you previously used the metadata profiler system:

**Old system:**
- Profiler auto-detected complexity
- Auto-selected annotation columns
- Required `assess_metadata_complexity()`

**New system:**
- Explicit config-based control
- User specifies annotation columns
- Simpler, more predictable

**Migration:**
1. Remove any profiler references from your config
2. Add `heatmap_annotations` if you want multi-column heatmaps
3. Add multiple `color` values if you want multiple PCA plots
4. Adjust `qc$thresholds` to match your dataset size

---

## Summary

- **Single color:** Simple, one PCA plot per PC pair
- **Multiple colors:** One PCA plot per color variable
- **Heatmap annotations:** Explicitly configured via `heatmap_annotations`
- **Adaptive plotting:** Skip expensive plots for large datasets
- **Validation:** Prevents invalid color/shape configurations

Your pipeline now adapts based on sample counts while giving you full control over plot aesthetics!
