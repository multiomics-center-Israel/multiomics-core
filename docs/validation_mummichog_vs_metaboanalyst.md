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

## How to run

```bash
make setup      # once: build the pinned venv, set MUMMICHOG_PYTHON in .Renviron
```

Point the script at the two MetaboAnalyst files (peak input + pathway output)
via the CONFIG block at the top of the script, or as env vars, then:

```bash
MA_INPUT="MetaboAnalyst Files/Mummichog/24h HL-LL/mummichog_input.txt" \
MA_OUTPUT="MetaboAnalyst Files/Mummichog/24h HL-LL/mummichog_pathway_enrichment_mummichog.csv" \
MODE=pos_default PPM=10 PERMUTATIONS=100 \
Rscript utils/validate_mummichog_vs_metaboanalyst.R
```

Match `MODE`/`PPM`/`PERMUTATIONS`/`CUTOFF` to what the technician used in
MetaboAnalyst. If you already have our pathway table from a pipeline run, set
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

## Results

_To be filled once the script is run (paste `comparison_summary.md` here):
Spearman ρ, top-N overlap, coverage, and a one-line verdict on agreement._
