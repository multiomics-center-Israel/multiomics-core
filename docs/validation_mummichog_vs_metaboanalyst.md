# Validation: pinned mummichog vs MetaboAnalyst (24h HL vs LL)

**Goal.** Confirm our version-pinned mummichog 2.7.0 engine
(`R/domain/metabolomics/06c`) agrees with a technician's MetaboAnalyst
mummichog run for the **24h High-Light vs Low-Light** contrast — apples-to-apples.

This is a **one-off validation**, not a pipeline feature: a single script
(`utils/validate_mummichog_vs_metaboanalyst.R`) plus this note. The engine is
**not modified**.

## Design

To isolate *engine vs engine* (rather than input drift), we feed our engine the
**exact peak input the technician submitted to MetaboAnalyst** (the
m/z · r.t · p.value · t.score table) and the **same metabolic model** (`cre`,
the published `model_ref` we already use), with the **same** ionization mode,
mass tolerance (ppm), permutations, and cutoff. Any remaining difference is then
attributable to the two implementations, not to different features or
thresholds.

(Running our own pipeline DE table for this contrast end-to-end — where the
feature set and significance cutoff differ from the MetaboAnalyst input — is a
separate, later check.)

### Ionization mode (mixed-mode limitation)

The technician's input carries a **per-feature `mode` column** — MetaboAnalyst
ran **mixed** (positive + negative) ionization. The pinned engine (standalone
mummichog 2.7.0) applies a **single global** mode (`-m`) and has no per-feature
mode; a `mode` column in the input is ignored. Per-feature mixed handling is
**roadmap item B3** (engine work, out of scope for this validation).

So the honest single-mode comparison here is **positive-only**: set
`MODE_FILTER=Positive` (subsets the input to positive features, reported at run
time) and `MODE=positive`, and compare against a **positive-only MetaboAnalyst
reference**. Comparing a positive-only run to the *mixed* reference is a
feature-set mismatch, not a clean validation. Full mixed-mode validation waits
on B3.

## How to run

```bash
make setup      # once: build the pinned venv, set MUMMICHOG_PYTHON in .Renviron
```

Point the script at the two MetaboAnalyst files (peak input + pathway output)
via the CONFIG block at the top of the script, or as env vars, then:

```bash
MA_INPUT="MetaboAnalyst Files/Mummichog/24h HL-LL/FOR UPLOAD HL vs LL 24h.csv" \
MA_OUTPUT="MetaboAnalyst Files/Mummichog/24h HL-LL/mummichog_pathway_enrichment_mummichog.csv" \
REF_PCOL="P(Fisher)" \
MODE_FILTER=Positive MODE=positive PPM=10 PERMUTATIONS=100 \
Rscript utils/validate_mummichog_vs_metaboanalyst.R
```

Match `MODE`/`PPM`/`PERMUTATIONS`/`CUTOFF` to what the technician used in
MetaboAnalyst. `MODE_FILTER=Positive` runs the positive-only subset (see the
mixed-mode note above) — pair it with a positive-only reference. `REF_PCOL`
defaults to `P(Fisher)` (the confirmed MetaboAnalyst mummichog p-value column).
If you already have our pathway table from a pipeline run, set
`OUR_PATHWAYS=.../tables/mcg_pathwayanalysis_*.tsv` to skip re-running the engine
and just compare.

### Outputs (written to `OUT_DIR`, default `mummichog_validation_24h_HL_LL/`)

- `comparison_pathways.csv` — matched pathways with `our_p`, `their_p`, and ranks.
- `comparison_summary.md` — coverage, Spearman ρ of p-values, Pearson r of
  −log10(p), and top-10 / top-20 overlap (Jaccard).

Both hold **aggregate pathway-level statistics only** (KEGG pathway names +
enrichment p-values) — no raw sample or feature values.

## Reading the comparison

- **Rank agreement (Spearman ρ) and top-N overlap are the primary signals.**
  Standalone mummichog 2.7.0 estimates significance by *unseeded* permutation, so
  absolute p-values vary between runs (even our own reruns) — do **not** expect
  exact equality.
- Pathways are matched by a normalized name (case/punctuation/spacing-insensitive).
  Unmatched pathways on either side are counted and usually reflect
  pathway-library version differences.
- We compare against MetaboAnalyst's **mummichog / Gamma** p-value column (the
  permutation-based one), not the Fisher (FET/EASE) columns. If auto-detection
  picks the wrong column, set `REF_PCOL` to the exact name.

## Expected differences (not bugs)

- **Stochasticity:** unseeded permutations → p-value jitter; rankings are stable,
  exact values are not.
- **Implementation:** standalone mummichog 2.7.0 ≠ MetaboAnalyst's
  reimplementation (permutation/gamma/EASE handling).
- **Library version:** the `cre` model snapshot vs MetaboAnalyst's pathway
  library edition can differ, changing pathway sizes and membership slightly.
- **Statistic type:** we compare our mummichog *permutation* p-value against
  MetaboAnalyst's `P(Fisher)` (its reported mummichog p for this run). These are
  related but not identical, so judge by rank agreement, not equal values.
- **Ionization mode:** single-mode (ours, positive-only) vs mixed (MetaboAnalyst)
  differ unless compared against a positive-only reference — see above.

## Results

_To be filled once the script is run (paste `comparison_summary.md` here):
Spearman ρ, top-N overlap, coverage, and a one-line verdict on agreement._
