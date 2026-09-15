---
name: results-deck
description: Build a slide deck or figure set from a pipeline run in this repo, where every number is verified against the artifact it came from and every claim carries its source. Use whenever asked for a deck, slides, a summary presentation, or a figure set describing RNA-seq, proteomics, metabolomics or lipidomics results. Also use when asked to check numbers that appear in an existing deck or report.
---

# Results deck

Turning a run into slides. The failure mode is not ugly slides, it is a number
on a slide that nobody can trace, or a claim the data does not support.

## 1. Read the run before designing anything

Never take figures from a report, an executive summary, or an earlier deck.
Those are downstream of the data and may be stale. Go to the run directory.

```bash
ls outputs/<run>/            # report, pipeline summary, per-mode folders
ls outputs/<run>/<mode>/Datasets/ outputs/<mode>/Enrichment/
cat outputs/<run>/<mode>/results_fact_sheet.tsv     # start here if present
```

`results_fact_sheet.tsv` already pairs the headline numbers with their sources.
When it exists, it is the fastest correct start; still confirm anything you plan
to put on a slide. Check `execution_info/git_commit.txt` and `timestamp.txt` — a
results directory is often older or newer than you assume, and a run may have
been repeated part-way through a session.

`docs/CONTRACT_rnaseq_final_results.md` explains what each column means.

## 2. Verify every figure and every plot

For each number: recompute it from the file it supposedly comes from. For each
plot: open it and confirm it shows what you are about to claim. Things that have
actually gone wrong here:

- A fold change of `2.97e+07` that meant "not detected in the other group",
  because the denominator was zero.
- A "top two components" claim where the groups did not separate at all; the
  spread came from two samples in one group.
- KEGG `Parkinson disease` and `Huntington disease` read as separate findings
  when they share 73% and 61% of their leading edge with oxidative
  phosphorylation.
- A quoted enrichment result from ranked `fgsea` presented alongside a claim
  about the differentially expressed gene list, which the over-representation
  files showed to be null.

Back the verified values into a TSV of `claim`, `value`, `source_file` beside
the deck. That table is the deliverable's audit trail and takes a minute.

## 3. Tone

The user's standing preference, and it is not negotiable:

- Do not oversell. "scores higher", "associated with", "one modest term" —
  never "switched on", "dominant axis", "strongly".
- No bold, no underline, no capitalised words inside a sentence. Small uppercase
  kicker labels above a title are fine.
- No speculation about mechanism or biological meaning that the run does not
  show. Report what was measured.
- Keep caveats explicit: n per group, single contrast, hypothesis-generating.
- Cite the exact analysis and file for every claim, in small muted type at the
  foot of the slide. Note that functional enrichment from different tools lives
  in different files and they are not interchangeable.
- Include a slide on what the data does not support. It is usually the most
  useful one.

## 4. Build

`pptxgenjs` is in the repo's `node_modules`. The shell working directory resets
between commands, so use an absolute `NODE_PATH`:

```bash
NODE_PATH=/home/ozsol/multiomics-core/node_modules node build_deck.js
```

Check image aspect ratios before sizing (`PIL.Image.open(p).size`) so plots are
not stretched. Keep 0.5 inch margins, 0.3 inch between blocks.

## 5. QA — this machine cannot render pptx normally

There is no LibreOffice, so `pptx -> pdf` fails. `pdftoppm` exists but has
nothing to convert. Do not stop at content checks.

```bash
python -m markitdown deck.pptx                                  # text and order
python scripts/render_pptx_slide.py deck.pptx --all out_dir/    # look at it
```

The renderer reconstructs each slide from the file's own geometry, substituting
DejaVu for Georgia and Calibri, so line breaks shift slightly. Say that when
showing the output.

**Then actually look at the images.** A geometry checker over positions,
overlaps and estimated overflow passed a slide whose text column ended at 60% of
the height while its figure ran to 85%. Only the rendered view caught it. Expect
to find problems on the first pass; if you find none, look harder.

When measuring gaps programmatically, compare rendered text height, not box
height — an oversized box reads as a collision when nothing is wrong. Kicker and
title, or a statistic and its label, are meant to sit together.

## 6. Deliver

Write the deck next to the run in `outputs/<run>/`, matching the naming already
used there (for example `Elad_Chiel_project_RNAseq_analysis_report_Aug2026.html`
implies `..._results_deck_Aug2026.pptx`). Ship the fact TSV alongside it. Say
plainly that visual rendering was a reconstruction and suggest one look in
PowerPoint before it is presented.
