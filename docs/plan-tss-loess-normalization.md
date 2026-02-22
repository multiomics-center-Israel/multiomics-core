# Architectural and Statistical Implementation Plan: TSS and LOESS Normalization

**Scope:** Add TSS (Total Sum Scaling), CLR (Centered Log-Ratio) transformation, cyclic LOESS (intensity-dependent bias correction), and QC-RLSC (signal drift correction) to the metabolomics normalization pipeline.

**Target files (primary changes):**
- `R/domain/metabolomics/01_normalization.R` — normalization functions and pipeline wrapper
- `R/domain/metabolomics/02_preprocess.R` — preprocessing orchestrator
- `R/domain/metabolomics/00_inputs.R` — config validation (`validate_metabolomics_config`)
- `config/templates/metabolomics_template.yaml` — config schema and defaults
- `R/core/06_plots.R` — new diagnostic plot functions
- `R/modules/metabolomics/01_mod_qc_pre.R` — QC diagnostic wrappers
- `tests/testthat/` — new test file(s) for metabolomics normalization

---

## 1. Conceptual Design

### 1.1 Where Normalization Sits in the Pipeline

The current pipeline in `02_preprocess.R` follows this order:

```
parse → metadata align → sample filter → feature filter → NA policy → normalization → method evaluation → missingness summary
```

Normalization itself is a 3-stage chain in `apply_normalization_pipeline()`:

```
Stage 1: sample_norm (column-wise: sum, median, pqn, is)
Stage 2: transform   (element-wise: log2, log10)
Stage 3: scaling     (row-wise: auto, pareto, range)
```

The new methods occupy distinct positions:

- **TSS** is a sample normalization method (column-wise division by column total). It belongs in **Stage 1** alongside the existing methods. It is dispatched as a new option `"tss"` in `normalize_samples()`.

- **CLR (Centered Log-Ratio)** is a transformation. It belongs in **Stage 2** alongside `log2` and `log10`. It is the statistically recommended companion to TSS for downstream analysis. Dispatched as a new option `"clr"` in `transform_metab()`.

- **Cyclic LOESS** corrects intensity-dependent bias between samples by fitting LOESS curves in M-A space (log-ratio vs. average log-intensity). It requires log-transformed data as input. It belongs in a **new Stage 3** (between current transformation and scaling stages). Dispatched via a new `correct_bias()` function.

- **QC-RLSC** corrects per-feature signal drift over injection order using pooled QC samples. It operates on raw linear intensities and must precede all normalization. It belongs in a **new Stage 0** (before sample normalization). Dispatched via a new `correct_drift()` function. Requires metadata columns `injection_order` and QC sample flag.

### 1.2 Extended Pipeline Architecture: 5-Stage Design

```
Stage 0: drift_correction      (NEW, optional) — QC-RLSC
Stage 1: sample_norm           (extended)       — none | sum | median | pqn | is | tss
Stage 2: transform             (extended)       — none | log2 | log10 | clr
Stage 3: bias_correction       (NEW, optional) — none | cyclic_loess
Stage 4: scaling               (unchanged)      — none | auto | pareto | range
```

The pipeline wrapper `apply_normalization_pipeline()` is extended from 3 calls to 5. The return value's `applied` record is extended to include all 5 stages plus their parameters.

### 1.3 Composability

TSS, LOESS, and QC-RLSC are **composable, not mutually exclusive**. They address orthogonal problems:

| Method | Problem Addressed | Scale | Scope |
|--------|-------------------|-------|-------|
| QC-RLSC | Instrument signal drift over time | Raw linear | Per-feature, per-batch |
| TSS | Unequal sample loading / total biomass | Raw linear | Per-sample (column-wise) |
| CLR | Map compositional data to Euclidean space | Log-ratio | Per-sample (column-wise) |
| Cyclic LOESS | Intensity-dependent systematic bias | Log | Global pairwise |

Valid compositions include:
- QC-RLSC → PQN → log2 → cyclic LOESS → pareto (full correction chain)
- QC-RLSC → TSS → CLR → none → none (compositional analysis with drift correction)
- none → TSS → log2 → cyclic LOESS → pareto (no drift data, compositional + bias correction)

Mutually exclusive only within Stage 1 (TSS vs PQN vs median — they occupy the same slot) and within Stage 2 (CLR vs log2 vs log10). This is enforced naturally by the single-method-per-stage config.

### 1.4 Configuration Exposure

The YAML config schema extends the existing `normalization:` block:

```yaml
normalization:
  # Stage 0: Drift correction (NEW)
  drift_correction: "none"              # "none" | "qc_rlsc"
  drift_correction_params:
    injection_order_col: "injection_order"
    qc_flag_col: "is_QC"
    batch_col: null                     # null = single batch; string = batch column name
    loess_span: "auto"                  # numeric (0,1] or "auto"
    min_qc_samples: 5                   # hard minimum per batch

  # Stage 1: Sample normalization (extended)
  sample_norm: "pqn"                    # ... | "tss"

  # Stage 2: Transformation (extended)
  transform: "log2"                     # ... | "clr"

  # Stage 3: Bias correction (NEW)
  bias_correction: "none"              # "none" | "cyclic_loess"
  bias_correction_params:
    loess_span: 0.4                    # LOESS span for MA fitting
    iterations: 3                      # cyclic iteration count

  # Stage 4: Scaling (unchanged)
  scaling: "none"
```

The `evaluate_methods` list is extended to optionally include `drift_correction` and `bias_correction` fields in each evaluation combo. Combos without these fields default to `"none"`.

### 1.5 Dispatch Pattern

Following the existing codebase convention, all new methods use `switch()` with `assert_one_of()` validation. Each new method is implemented as a standalone function with a consistent signature, dispatched from a stage-level dispatcher function:

- `normalize_samples()` — extended switch (add `"tss"`)
- `transform_metab()` — extended switch (add `"clr"`)
- `correct_bias(mat, method, params)` — new dispatcher with switch
- `correct_drift(mat, meta, method, params)` — new dispatcher with switch

This preserves consistency with all three omics modalities (`normalize_counts()`, `transform_metab()`, `scale_metab()` all use this pattern). The functions remain registry-ready: each method is a standalone function that could later be registered in a lookup table without changing its implementation.

---

## 2. Statistical Considerations

### 2.1 TSS Assumptions in Metabolomics

TSS converts absolute intensities to relative abundances by dividing each sample (column) by its total intensity. The output is compositional data constrained to the simplex (all values in [0,1], summing to 1.0 per sample).

**Core assumptions:**
- Total metabolite pool size is not biologically informative (or is a nuisance variable to be removed).
- Relative abundances are the quantities of interest for downstream comparison.
- The feature set represents a sufficiently complete snapshot of the metabolome such that relative proportions are stable representations.

**Compositional constraints (Aitchison geometry):**
- Closure induces spurious negative correlations between features. If one feature's proportion increases, others must decrease even if their absolute abundance is unchanged.
- Standard Pearson/Spearman correlations, Euclidean distances, and PCA on raw proportions are biased. The mathematically correct space for compositional data is the Aitchison simplex, accessed via log-ratio transformations (CLR, ALR, ILR).
- This is why the plan includes CLR as a companion transformation.

**Implication for limma-based downstream analysis:**
Limma operates on log-transformed data and assumes approximate normality. TSS + log2 is compatible with limma, but the user must understand that differential analysis results indicate changes in *relative abundance* (compositional shifts), not absolute intensity changes. TSS + CLR is statistically more principled, as CLR maps compositional data into unconstrained Euclidean space where standard linear models are valid.

### 2.2 When TSS Distorts Biological Signal

TSS is problematic in specific scenarios:

1. **Dominant metabolites:** If one or a few metabolites constitute >20% of total intensity, all other metabolites' relative abundances become inversely coupled to the dominant feature. A change in a single dominant metabolite artificially shifts all others. The pipeline should compute a dominance metric (maximum feature contribution to total per sample) and warn if it exceeds a threshold.

2. **Differential total biomass:** If treatment genuinely increases total metabolite pool size (e.g., metabolically active vs. quiescent cells), TSS removes this information by forcing all samples to equal total. This masks biologically meaningful global shifts.

3. **Sparse data with many zeros:** In sparse matrices, TSS distributes the total across fewer detected features, inflating their apparent relative abundance. Features with many zeros across samples get disproportionately high proportions in samples where they are detected.

4. **Rare/low-abundance features:** Small absolute changes in low-abundance features produce large relative changes after TSS, inflating apparent variability. This interacts with the pseudocount in log transformation.

### 2.3 LOESS Use Cases

**Use Case A: Intensity-dependent bias correction (cyclic LOESS)**

Cyclic LOESS normalizes pairwise differences between samples by fitting LOESS curves through M-A plots:
- M = log2(sample_i) - log2(sample_j) for each feature
- A = 0.5 * (log2(sample_i) + log2(sample_j))

The fitted M~A curve captures systematic intensity-dependent bias: if low-intensity features are shifted relative to high-intensity features, the LOESS curve deviates from M=0. The correction subtracts the fitted M values.

"Cyclic" means this pairwise normalization is iterated across all sample pairs until convergence (typically 3 iterations). `limma::normalizeCyclicLoess()` implements this directly and is already available (limma 3.66.0 is in renv.lock).

**Core assumption:** Most features are NOT differentially abundant between any pair of samples. The LOESS curve captures technical bias because the majority of features define a "null" background trend. If a large fraction of features are truly differential (e.g., >30%), the LOESS fit is biased toward removing real biological signal.

**Use Case B: QC-RLSC signal drift correction**

QC-RLSC uses pooled QC samples (identical aliquots injected at regular intervals) to estimate and correct signal drift over injection order. For each feature independently:
1. Extract QC sample intensities at their injection positions.
2. Fit a LOESS curve: QC_intensity ~ injection_order.
3. Predict expected drift at each biological sample's injection position.
4. Normalize each biological sample's intensity relative to the LOESS-predicted value (typically by dividing by the predicted value and multiplying by the median QC level).

This corrects for:
- Gradual sensitivity loss or gain in the mass spectrometer.
- Carry-over effects that accumulate over a run.
- Column degradation in LC-MS.

### 2.4 LOESS Operating Scope

**Cyclic LOESS:** Operates globally across all samples. It does not require batch or group information. For the initial implementation, apply it globally. Reserve an optional `batch_col` parameter for future within-batch application if batch-specific intensity biases exist.

**QC-RLSC:** Operates per-batch if a `batch_col` is specified (each batch has its own QC trend and drift characteristics). If `batch_col` is null (single run), applies globally. Within each batch, the LOESS fit uses only that batch's QC samples. The critical requirement is that QC samples are distributed throughout the injection order range within each batch.

### 2.5 Overfitting Risks

**Cyclic LOESS span:**
- Default in `limma::normalizeCyclicLoess` is 0.7. For metabolomics with hundreds to thousands of features, this is generally appropriate.
- If the number of features is very small (<50), a smaller span could overfit to noise. The pipeline should warn if features < 100.
- The span must be between 0.1 and 1.0; validated at config time.

**QC-RLSC span and QC sample count:**
- With fewer than 5 QC samples, the LOESS fit is underdetermined. Hard minimum: 5 QC samples per batch.
- With 3-4 QC samples, downgrade to degree=1 (locally linear) with a warning.
- Fewer than 3 QC samples: refuse to run with an informative error.
- `"auto"` span mode: set span to `max(0.5, 3 / n_qc)`, ensuring at least 3 QC points per local fit window.
- QC samples should span the injection order range. Warn if QC injection orders cover less than 50% of the full range.

### 2.6 CLR Statistical Properties

CLR (Centered Log-Ratio) for a sample vector x:
```
CLR(x_i) = log(x_i / geometric_mean(x))
```

This maps compositional data from the simplex to unconstrained Euclidean space where:
- Standard Euclidean distances are equivalent to Aitchison distances on the simplex.
- PCA on CLR-transformed data yields a principal components analysis in Aitchison geometry.
- Linear models (limma) operate in the correct mathematical space.

**Key consideration:** CLR requires all values to be strictly positive. Zeros must be handled before CLR (either by prior imputation, pseudocount addition, or a multiplicative replacement strategy). When `transform = "clr"`, the pseudocount parameter from the config is added before computing the geometric mean and log-ratios.

**Relationship to log2:** CLR produces values centered per sample (mean of CLR values = 0 for each sample). Log2 does not have this centering property. If `sample_norm = "tss"`, CLR is preferred over log2 because it respects the compositional geometry. If `sample_norm` is not TSS (e.g., PQN), log2 remains the standard choice.

---

## 3. Data Integrity & Edge Cases

### 3.1 Zeros and Missing Values

The pipeline handles NAs via `na_policy` ("keep" or "zero") before normalization (lines 73-86 of `02_preprocess.R`). Existing normalization functions use `na.rm = TRUE` throughout. New methods must maintain this convention.

**TSS with NAs:** `colSums(mat, na.rm = TRUE)` computes the total of non-NA values. Division by this sum means relative abundances sum to 1.0 *for observed values only*. This is consistent with existing TIC behavior. Document explicitly.

**TSS with zeros:** Zero values remain zero after TSS (0 / total = 0). With subsequent log transformation, pseudocount ensures `log2(0 + 1) = 0`. This is the expected behavior.

**CLR with zeros:** CLR requires strictly positive values because it computes a geometric mean (which involves a product of values). If any value is zero, the geometric mean is zero and the log-ratio is undefined. The pseudocount parameter handles this: add pseudocount before computing CLR. The pipeline should validate that pseudocount > 0 when `transform = "clr"`.

**Cyclic LOESS with NAs:** `limma::normalizeCyclicLoess` does not natively handle NAs in the input matrix. Two approaches: (a) require that NAs have been resolved before reaching Stage 3 (recommend `na_policy: "zero"` when using cyclic LOESS), or (b) temporarily fill NAs with row medians for the LOESS fit, then restore NAs in the output. Approach (a) is cleaner and should be the primary recommendation, with approach (b) as a documented fallback.

**QC-RLSC with NAs in QC samples:** The per-feature LOESS fit proceeds with available (non-NA) QC points. If all QC values for a feature are NA, that feature cannot be drift-corrected; leave unchanged and flag it in the return metadata.

### 3.2 Near-Zero Total Sum for TSS

If a sample's total intensity is zero or near-zero (all features are zero or NA), TSS produces division by near-zero/zero. Handling:
- Compute column sums with `na.rm = TRUE`.
- Define epsilon as `1e-10 * median(col_sums[col_sums > 0])`.
- If a sample's total < epsilon: emit a warning naming the sample, leave that sample's values unchanged (do not divide), and record the sample ID in the return metadata (`tss_near_zero_samples`).
- This is consistent with the existing `norm_total_sum()` approach of setting scale factors to 1 for zero-sum samples (line 46 of `01_normalization.R`).

### 3.3 TSS with Negative Values

Negative values can arise from baseline correction in mass spectrometry. TSS requires non-negative input (negative proportions are meaningless).

Handling: Before computing column sums, clamp negative values to zero with a warning reporting the count and percentage clamped. This follows MetaboAnalyst convention. Record the count in return metadata (`tss_negative_clamped`).

### 3.4 Interaction with Log Transformation

- **TSS must occur before log transformation (Stage 1 before Stage 2).** TSS operates on raw linear intensities. The "total" is meaningful on a linear scale. log(TSS(x)) = log(x/total) = log(x) - log(total), which correctly yields log-ratios relative to a sample constant.

- **CLR replaces the separate TSS + log2 operation when used.** CLR = log(x / geometric_mean(x)) is both a compositional normalization and a log transformation. When `transform = "clr"`, it implicitly assumes the input may be compositional (TSS-normalized) or not. If `sample_norm = "tss"` and `transform = "clr"`, the pipeline first produces proportions, then applies CLR (log of proportion divided by geometric mean of proportions). This is mathematically valid and is the standard approach for compositional data analysis.

- **Cyclic LOESS requires log-transformed input (Stage 3 after Stage 2).** M-A statistics are only meaningful in log-space. The pipeline must validate that `transform != "none"` when `bias_correction = "cyclic_loess"`, or emit an error at config validation time.

- **QC-RLSC operates on raw linear intensities (Stage 0 before Stage 1).** Drift is multiplicative on the linear scale. Correcting in log-space would convert multiplicative drift to additive, which also works but deviates from the standard QC-RLSC literature.

### 3.5 LOESS with Sparse QC Samples

Minimum viability per batch:
- LOESS degree=2 (locally quadratic, default): minimum 4 QC samples per fit window.
- LOESS degree=1 (locally linear): minimum 3 QC samples per fit window.
- With span=0.75 and 6 QC samples → window covers 4-5 points → barely sufficient for degree=2.

Pipeline enforcement:
- >= 5 QC samples per batch: proceed normally (configurable threshold via `min_qc_samples`).
- 3-4 QC samples: downgrade to degree=1 with a warning.
- < 3 QC samples: refuse with an error.
- QC distribution check: warn if QC injection orders span less than 50% of the full injection order range within a batch.

---

## 4. Order of Operations

### 4.1 Full Extended Pipeline Order

The complete `preprocess_metabolomics()` flow becomes:

```
 1. Parse input                              (existing, unchanged)
 2. Metadata alignment                       (existing, unchanged)
 3. Sample filtering (blanks, QC exclusion)  (existing, unchanged)
 4. Feature filtering (all-NA / all-zero)    (existing, unchanged)
 5. Missing value policy (keep / zero)       (existing, unchanged)
─── normalization chain begins ───
 6. Stage 0: Drift correction (QC-RLSC)     (NEW, optional)
 7. Stage 1: Sample normalization            (extended: + tss)
 8. Stage 2: Transformation                  (extended: + clr)
 9. Stage 3: Bias correction (cyclic LOESS)  (NEW, optional)
10. Stage 4: Scaling                         (unchanged)
─── normalization chain ends ───
11. Method evaluation                        (existing, extended for new stages)
12. Missingness summary                      (existing, unchanged)
```

### 4.2 Rationale for Stage Ordering

**QC-RLSC before sample normalization (Stage 0 → Stage 1):**
QC-RLSC is feature-specific — different metabolites drift at different rates. If sample normalization (PQN, median, TSS) is applied first, the per-sample scale factors absorb some drift in a feature-nonspecific way, destroying the feature-level drift information that QC-RLSC needs. Correct order: remove feature-specific drift first, then adjust residual sample-level differences.

**TSS before log transformation (Stage 1 → Stage 2):**
TSS divides by column totals on the linear scale. This total is only meaningful for raw (or drift-corrected) intensities. Applying log first would make the "total" a sum of log-values, which has no clear biological interpretation.

**CLR in the transformation stage (Stage 2):**
CLR = log(x / geometric_mean(x)). It is both a compositional mapping and a log transformation. It replaces log2/log10 in the chain, not supplements it. The pipeline should not apply log2 *after* CLR.

**Cyclic LOESS after log transformation (Stage 2 → Stage 3):**
Cyclic LOESS operates in M-A space where M = log2(sample_i) - log2(sample_j). These are only meaningful for log-transformed data. Additionally, the symmetry assumption (most features have M near 0) holds better in log-space where multiplicative biases become additive.

**Scaling after bias correction (Stage 3 → Stage 4):**
Scaling (auto, pareto, range) adjusts per-feature variance. It should be the last step to ensure that all upstream corrections are in place before variance standardization.

### 4.3 Config Validation Constraints

The following cross-stage constraints should be enforced at config validation time (in `validate_metabolomics_config()`):

1. If `bias_correction = "cyclic_loess"`, then `transform` must not be `"none"` (LOESS requires log-transformed input).
2. If `transform = "clr"`, then `pseudocount` must be > 0 (CLR requires strictly positive values).
3. If `drift_correction = "qc_rlsc"`, then `drift_correction_params` must specify `injection_order_col` and `qc_flag_col`.
4. If `sample_norm = "tss"` and `transform = "log2"`, emit an informational message recommending CLR instead (not an error — log2 is valid, just suboptimal for compositional data).

---

## 5. Reproducibility & Auditability

### 5.1 Metadata Stored in Output Objects

The `applied` field in the return value of `apply_normalization_pipeline()` (currently at lines 241-248 of `01_normalization.R`) is extended:

**Current record:**
```
sample_norm, transform, scaling, pseudocount
```

**Extended record:**

```
pipeline_version: "2.0"

# Stage 0
drift_correction: "qc_rlsc" | "none"
drift_params:
  injection_order_col: <string>
  qc_flag_col: <string>
  batch_col: <string | null>
  loess_span_used: <numeric>         # actual span (may differ from "auto" input)
  qc_samples_used: <character vector> # per batch
  features_skipped: <character vector> # features with insufficient QC data
  median_correction_factor: <numeric vector> # per-feature summary of drift magnitude

# Stage 1
sample_norm: "tss" | "pqn" | ...
sample_norm_params:
  tss_total_sums: <named numeric vector>  # per-sample totals used as denominators
  tss_near_zero_samples: <character vector>
  tss_negative_clamped: <integer>         # count of negatives clamped to 0

# Stage 2
transform: "clr" | "log2" | ...
pseudocount: <numeric>
clr_geometric_means: <named numeric vector>  # per-sample geometric means (if CLR)

# Stage 3
bias_correction: "cyclic_loess" | "none"
bias_params:
  loess_span: <numeric>
  iterations: <integer>
  converged: <logical>

# Stage 4
scaling: "pareto" | "auto" | ...
```

### 5.2 Recording Parameters for Downstream Reproducibility

The preprocessing info record (`pre$info$normalization` at lines 135-143 of `02_preprocess.R`) should capture the full `applied` record. This is already serialized as `normalization_applied.tsv` by the QC module. For the extended pipeline:

- The TSV output should flatten the nested record for human readability (one row of key-value pairs, or a multi-row "parameter: value" table).
- The Shiny payload (built in `07_shiny_export.R`) constructs `config_norm_method` as a slash-delimited label `"pqn/log2/none"`. Extend to 5 fields: `"none/tss/clr/cyclic_loess/none"` (drift/sample/transform/bias/scaling).
- The full nested list should also be saved as an RDS for programmatic access downstream.

### 5.3 Diagnostics That Should Be Auto-Generated

New plots added to the QC module (`mod_metabolomics_qc_pre`), following the pure `plot_*` function pattern in `R/core/06_plots.R`:

**For TSS:**
- **Total intensity barplot:** Per-sample total intensities before TSS. Visualizes sample loading variation. Flag samples below the 5th percentile or above the 95th percentile.
- **Dominance barplot:** For each sample, fraction contributed by the top 1, top 5, and top 10 features. Warns if any single feature exceeds 20% of total.

**For CLR:**
- **CLR distribution density:** Per-sample density overlay of CLR-transformed values (compare to log2 density overlay to show the centering effect).

**For Cyclic LOESS:**
- **MA plot before/after:** For a representative sample pair (the pair with the highest median absolute M before correction), show M vs A scatter with LOESS curve overlaid (before) and the corrected M vs A (after). The before-plot should show intensity-dependent bias; the after-plot should show M centered at zero.
- **M-value density before/after:** Overlay of M-value densities across all sample pairs, before and after correction. Should shift from skewed to symmetric around zero.
- **Boxplot of |M| per sample:** Before and after correction, showing reduction in systematic bias magnitude.

**For QC-RLSC:**
- **Drift plot per representative feature:** Top 3-5 features by drift magnitude. Each shows: intensity (y-axis) vs injection order (x-axis), with QC samples highlighted, LOESS fit curve overlaid, and corrected values shown. Before and after panels.
- **QC RSD comparison:** Per-feature RSD of QC samples before and after drift correction, shown as a paired scatter or box plot. Should decrease.
- **PCA of QC samples before/after:** QC samples should cluster more tightly after correction.
- **Correction factor heatmap:** Features (rows) vs samples (columns), colored by correction factor magnitude. Visualizes drift pattern across the run.

---

## 6. Testing Strategy

### 6.1 Unit Testing Logic

New tests in `tests/testthat/test-metab-normalization.R`, following the existing `test-normalize.R` pattern.

**TSS unit tests:**
- Dimension preservation: output has same nrow/ncol as input.
- Compositional closure: each column sums to 1.0 (within floating-point tolerance, e.g., `abs(colSum - 1.0) < 1e-12`).
- Non-negativity: all output values >= 0 given non-negative input.
- Identity case: matrix with uniform columns (all features equal) produces uniform proportions (1/n_features per value).
- NA preservation: NAs remain NA; non-NA values in the same column still sum to 1.0.
- Zero-sum column: column of all zeros returns all zeros with a warning.
- Near-zero column: column with total < epsilon triggers warning and leaves values unchanged.
- Negative clamping: negative input values are clamped to zero with a warning; remaining values are normalized correctly.
- Row/column names preserved.
- Single-feature matrix: TSS produces all 1.0s.
- Single-sample matrix: TSS produces correct proportions.

**CLR unit tests:**
- Dimension preservation.
- Per-sample mean of CLR values is zero (within tolerance): `all(abs(colMeans(clr_mat)) < 1e-10)`.
- Known-answer test: for a 3-feature sample with values (2, 8, 32) + pseudocount=0, CLR should yield `log2(c(2,8,32)) - mean(log2(c(2,8,32)))`.
- Pseudocount handling: with zeros in input, pseudocount > 0 ensures finite output.
- Error when pseudocount = 0 and zeros exist (produces -Inf).
- Row/column names preserved.

**Cyclic LOESS unit tests:**
- Dimension preservation.
- Row/column names preserved.
- No-bias case: matrix of identical samples → output approximately equals input (within floating-point tolerance).
- Config validation: error when `transform = "none"` and `bias_correction = "cyclic_loess"`.
- Span boundary: span=1.0 and span=0.1 both produce valid output without errors.
- Warning on very few features (<50).

**QC-RLSC unit tests:**
- Dimension preservation (all samples, including QCs, present in output).
- Row/column names preserved.
- No-drift case: constant QC values → output approximately equals input.
- Error when fewer than 3 QC samples.
- Warning when 3-4 QC samples (downgrade to degree=1).
- Error when `injection_order_col` missing from metadata.
- Error when `qc_flag_col` missing from metadata.
- Features with all-NA QC values are left unchanged and reported.
- Correction factors are all positive (no sign inversions).

### 6.2 Statistical Validation Tests

**TSS closure test:**
- Create a matrix with known, unequal column sums (e.g., column totals of 100, 500, 1000). Verify TSS equalizes them to 1.0 while preserving within-sample proportions.

**CLR isometry test:**
- Create two known compositional vectors. Compute Aitchison distance in the simplex and Euclidean distance in CLR space. Verify they are equal.

**LOESS bias recovery test:**
- Construct a synthetic log-matrix where M_biased = M_true + f(A) for a known function f (e.g., f(A) = 0.5 * sin(A) or f(A) = 0.3 * A).
- Apply cyclic LOESS.
- Verify that residual bias (mean |M| across A bins) is reduced by at least 80% compared to pre-correction.

**QC-RLSC drift recovery test:**
- Construct a dataset with known linear drift: intensity_observed = intensity_true * (1 + rate * injection_order), with per-feature rates drawn from U(0.005, 0.02).
- Insert QC samples at known positions with the same drift.
- Apply QC-RLSC.
- Verify corrected intensities are within 5% of true intensities for all features and samples.

### 6.3 Simulation Ideas

**Intensity-dependent bias simulation (end-to-end):**
1. Generate a reference matrix: 500 features, 20 samples, log-normal intensities.
2. For each pair of samples, inject intensity-dependent bias: M_biased = M_true + beta * A, with beta = 0.3.
3. Reconstruct biased raw intensities.
4. Run full pipeline with cyclic LOESS enabled.
5. Compare DE results (limma) with and without cyclic LOESS against the known truth (which features were spiked as differential).
6. Measure: false positive rate, sensitivity, and bias curve slope before/after.

**Drift simulation (end-to-end):**
1. Generate a reference matrix: 200 features, 50 biological samples.
2. Assign injection orders 1-50. Insert QC at positions 1, 10, 20, 30, 40, 50.
3. Apply per-feature random drift (rate from U(0.005, 0.02) per injection step).
4. Run full pipeline with QC-RLSC enabled.
5. Measure: QC RSD before/after, PCA separation of biological groups before/after, mean absolute deviation from truth.

**Compositional artifact simulation:**
1. Create a matrix with one dominant metabolite (50% of total) that increases 2-fold between conditions.
2. Apply TSS + log2 and TSS + CLR.
3. Run limma DE on both.
4. Verify that TSS + log2 shows spurious changes in non-differential features (the compositional artifact), and that TSS + CLR mitigates this.
5. Compare with PQN + log2 as a non-compositional baseline.

---

## 7. Future Extensibility

### 7.1 Normalization as a Plug-In Module

The current `switch()` + `assert_one_of()` dispatch pattern works well for the current scale (6-7 methods at Stage 1, 3-4 at Stage 2). The new methods bring the total to 8 for Stage 1 and 4 for Stage 2, with 2 each for the new Stages 0 and 3. This is still manageable with `switch()`.

**Registry pattern for the future:** When any stage exceeds approximately 8-10 methods, or when external contributors want to add methods without modifying core files, introduce a formal registry:

- A `NormRegistry` environment (or named list) with `register(stage, name, fn, validator)` and `dispatch(stage, name, mat, ...)`.
- Registration occurs at source time (when `tar_source` loads all R files).
- The validator function checks stage-specific config parameters at planning time.
- Each method function maintains the same signature as today.

**For now:** Keep `switch()` dispatch, but structure all new methods as standalone functions with consistent signatures. A future migration to a registry would only change the dispatch mechanism, not the method implementations themselves.

### 7.2 Scaling for Future Methods

| Method | Stage | Notes |
|--------|-------|-------|
| Quantile normalization | Stage 1 or new Stage 1.5 | `limma::normalizeQuantiles`. Forces all samples to have identical distributions. More aggressive than PQN. Less common in metabolomics due to lower feature count. Add as `sample_norm: "quantile"`. |
| Median normalization | Stage 1 | Already exists as `"median"`. No changes needed. |
| VSN | Combined Stage 1+2 | Variance-stabilizing normalization combines normalization and transformation. Would replace both Stage 1 and Stage 2. Requires special handling: if `sample_norm: "vsn"`, skip Stage 2 (set `transform: "none"` implicitly). Uses `vsn::justvsn()`. |
| Cyclic LOESS | Stage 3 | Already designed. Future sub-variants: fast LOESS (for >10k features), weighted LOESS, per-batch LOESS. Add as sub-options in `bias_correction_params`. |

The 5-stage pipeline architecture naturally accommodates these additions:
- New methods at existing stages: add to the `switch()` and `assert_one_of()` allowed values.
- Combined methods (like VSN): add conditional logic to skip downstream stages, with config validation enforcing compatible combinations.
- New sub-variants: use the `_params` sub-block for method-specific tuning without expanding the top-level config.

### 7.3 Config Evolution

The config schema is designed to be forward-compatible:
- Each stage has a method selector string and an optional `_params` sub-block.
- New methods only require extending the allowed values in `assert_one_of()` and adding their `_params` validation.
- The `evaluate_methods` list supports arbitrary combinations of all 5 stages, so new methods are automatically available for comparative evaluation.
- No breaking changes to existing configs: the new Stage 0 and Stage 3 default to `"none"`, making them backward-compatible.
