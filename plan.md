# Plan: Multi-Level File Merge for Metabolomics Pipeline

## Overview

Add support for accepting a **folder** of multiple metabolomics data files (each representing a different confidence level), generating unique IDs in `RT{rt}_MZ{mz}` format, merging all files into a single table with a `level` column, and validating global ID uniqueness with `stop()` on failure.

**Key design decision:** The level is derived from the filename (sans extension). All files share the same sample columns (same experiment, different feature lists).

---

## Files to Modify (4 files) + 1 New Test File

| # | File | Change Type |
|---|------|-------------|
| 1 | `R/domain/metabolomics/00_inputs.R` | **Primary** — add 3 new functions, modify 2 existing |
| 2 | `R/domain/metabolomics/02_preprocess.R` | **1-line change** — add `multi_level` branch to switch |
| 3 | `R/pipeline/metabolomics/00_pipe_metabolomics.R` | **Modify** `metab_input_files` target for folder tracking |
| 4 | `config/templates/metabolomics_template.yaml` | **Add** `multi_level` config section + example |
| 5 | `tests/testthat/test-multi-level-metabolomics.R` | **New** — unit tests |

**No changes needed** in: `01_normalization.R`, `05_outputs_legacy.R`, `07_shiny_export.R`, `02_validation.R`, QC modules, or Shiny export. The `level` column flows through automatically via `row_data`.

---

## Step-by-Step Implementation

### Step 1: Add `build_feature_ids_mzrt()` to `00_inputs.R`

New function alongside existing `build_feature_ids()` (which remains unchanged for backward compatibility).

- Takes `data_df` and `id_cfg`, extracts `mz_col` and `rt_col` values
- Produces IDs in format `RT8.93_MZ308.09092` (RT rounded to 2 decimals, m/z at full precision)
- `stop()` if m/z or RT columns are missing or contain NA values
- Does **NOT** call `make.unique()` — uniqueness is validated at merge time

### Step 2: Add `load_multi_level_files()` to `00_inputs.R`

New function that discovers and reads all files from a folder.

- Resolves `files$data` as a folder path via `resolve_raw_path()`
- Lists all `.xlsx/.xls/.csv/.tsv/.txt` files in the folder
- Derives the level label from each filename (sans extension)
- Reads each file using existing `read_metab_file()`
- Returns a list of `{data, level, file_path}` entries
- `stop()` if folder doesn't exist or contains no data files

### Step 3: Add `merge_multi_level_data()` to `00_inputs.R`

Core merge function. For each file:

1. Parse through existing `parse_cd_raw()` or `parse_processed_wide()` (reuses current parsers)
2. Generate IDs via `build_feature_ids_mzrt()` (overrides parser-generated IDs)
3. Validate within-file ID uniqueness — `stop()` on duplicates
4. Add `level` column to `row_data`
5. Validate that sample columns match across all files — `stop()` if they differ

After processing all files:

6. Validate **cross-file** global ID uniqueness — `stop()` if any ID appears in more than one file
7. `rbind` all expression matrices (features stacked vertically, same sample columns)
8. `rbind` all `row_data` frames (with `level` column)
9. Merge sample maps (if CD format)
10. Return `{expr_raw, row_data, sample_map, sample_ids}` — same structure as single-file parsers

### Step 4: Modify `load_metabolomics_inputs()` in `00_inputs.R`

Add a branch for `format == "multi_level"`:

- Calls `load_multi_level_files()` instead of reading a single file
- Still loads metadata if provided
- Returns `list(data = NULL, file_list = [...], metadata, format = "multi_level")`
- Existing single-file path remains completely unchanged

### Step 5: Update `validate_metabolomics_config()` in `00_inputs.R`

- Add `"multi_level"` to the allowed format values
- Add validation for the new `multi_level` config section (`file_format`, `level_source`)

### Step 6: Add `multi_level` branch in `preprocess_metabolomics()` — `02_preprocess.R`

One-line addition to the `switch()` at line ~24:

```r
multi_level = merge_multi_level_data(inputs$file_list, cfg, inputs$metadata),
```

Everything downstream (metadata alignment, filtering, normalization, QC) operates on the merged data unchanged.

### Step 7: Update `metab_input_files` target — `00_pipe_metabolomics.R`

When format is `multi_level`, enumerate all files in the folder so `{targets}` tracks each one as a dependency. If any file is added/removed/modified, the pipeline invalidates and re-runs.

### Step 8: Update config template — `metabolomics_template.yaml`

Add a new `multi_level` config section:

```yaml
    multi_level:
      # Format of each individual file in the folder
      file_format: "cd_raw"    # "cd_raw" | "processed_wide"
```

Add a commented example showing how to use `format: "multi_level"` with `files.data` pointing to a folder.

### Step 9: Add unit tests — `test-multi-level-metabolomics.R`

Test cases:
1. `build_feature_ids_mzrt()` produces correct `RT_MZ` format
2. `build_feature_ids_mzrt()` stops on missing/NA m/z or RT
3. `merge_multi_level_data()` merges correctly with distinct IDs
4. `merge_multi_level_data()` stops on cross-file duplicate IDs
5. `merge_multi_level_data()` stops on within-file duplicate IDs
6. `merge_multi_level_data()` stops when sample columns differ between files
7. Level column is correctly populated in merged `row_data`
8. Backward compatibility: `cd_raw` and `processed_wide` formats still work identically

---

## Data Flow (multi_level path)

```
Folder of files (Level_1.xlsx, Level_2.xlsx, ...)
    │
    ▼
load_metabolomics_inputs()  ──── format == "multi_level"
    │ reads each file via read_metab_file()
    │ returns file_list + metadata
    │
    ▼
preprocess_metabolomics()
    │ dispatches to merge_multi_level_data()
    │   ├── per file: parse_cd_raw() → build_feature_ids_mzrt() → add level
    │   ├── validate within-file uniqueness (stop on fail)
    │   ├── validate cross-file uniqueness (stop on fail)
    │   ├── validate sample column consistency (stop on fail)
    │   └── rbind expr_raw + row_data → unified output
    │
    │ metadata alignment (unchanged)
    │ sample filtering (unchanged)
    │ feature filtering (unchanged)
    │ normalization pipeline (unchanged)
    │
    ▼
Outputs: feature_annotations.tsv now includes "level" column
         expr_normalized.tsv has features from all files
         QC / Shiny export work unchanged
```

---

## What Stays the Same

- `build_feature_ids()` — unchanged, continues producing `Name_mzXXXX_rtXX` for cd_raw/processed_wide
- All normalization logic — unchanged
- All QC modules — unchanged
- Output writers — unchanged (level column flows through `row_data` automatically)
- Shiny export — unchanged
- Existing `cd_raw` and `processed_wide` formats — fully backward compatible
