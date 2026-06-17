# Project conventions for Claude Code

This is a `{targets}`-based bioinformatics pipeline in R, maintained by a small team of 5. **Read this file before doing anything in this repo.**

The rules below protect our architecture and our data. The last three sections (creativity, breaking rules, when to ask) protect *you* from being too cautious in the wrong places.

---

## Before you write any code

1. **Read the existing structure first.** Look at `_targets.R` and the files in `R/` before adding anything. Notice the patterns — function naming, file organization, target naming. Match what's there.

2. **Search before adding.** If you're about to write a helper, first search the repo for something similar:
   ```bash
   grep -rn "<keyword>" R/
   ```
   If a similar function exists, extend it or call it. Don't create a duplicate. Tell the user what you found.

3. **For non-trivial tasks, propose a plan first — then wait for approval.** If your change will touch more than one file, introduce a new pattern, or take more than a few minutes — *don't dive in*. Write a short plan (3–6 bullets: what files, what functions, what tests) and let the user approve it before you start coding. This single habit prevents most rework.

4. **Don't create new top-level directories.** The structure is intentional. New directories require a conversation with the user first.

5. **Don't change `renv.lock`** unless the user explicitly asked for a dependency change. If your work seems to need a new package, stop and ask.

---

## Where things live

```
_targets.R         pipeline definition only — NO helper functions go here
R/                 source of truth for every function, in a strict 5-layer load order:
  core/              generic utilities (io, validation, plotting, Excel export, clustering)
  services/          external integrations (AI commentary)
  domain/<omics>/    omics-specific logic (rnaseq, proteomics, metabolomics, multiomics)
  modules/<omics>/   target-ready wrappers (mod_*) around domain logic
  pipeline/<omics>/  {targets} orchestration (sourced LAST)
tests/             testthat tests
config/            YAML configs; templates/ is gitted, per-project configs are gitignored
renv.lock          dependency lockfile — gitted, do NOT casually change
_targets/          pipeline cache — gitignored, NEVER touch
```

> **Repo map:** see `PROJECT_STRUCTURE.md` for the full scanned map (layers, per-file roles, config layout, Shiny contract) — read it after this file. If the layout here drifts from reality, fix it — *but only after asking the user.*

---

## Metabolomics specifics

A few conventions specific to the metabolomics mode:

- **`preprocessing` is the single source of truth for normalization.** `preprocessing.chosen_norm` selects the method (`tss`/`median`/`pqn`/`eigenms`/`eigenms_forced`; `null` = QC-review mode); `scaling`, `pseudocount`, `sample_norm` and `transform` live alongside it. There is **no separate `normalization:` block** and **no imputation stage** — only missingness *filtering*.
- **Say "features"/"metabolites", not "genes".** User-facing plot titles and report text use "features"/"metabolites" and "differential abundance" — not "genes"/"differential expression". But **don't rename internal identifiers** (the `de:` config key, `de_res`/`de_tables`, DE targets, Shiny payload keys like `de_stats`): those are shared contracts, and the DE→DA rename of internal names is a deliberate, separate task.
- **Final Excel: `original_id` comes immediately after `feature_id`.**

---

## Naming and documentation

**Naming:**

- **Functions:** `snake_case` (`load_counts`, not `loadCounts` or `LoadCounts`).
- **One concept per R file.** Split files when they grow past ~300 lines.
- **Targets:** `<verb>_<noun>` (e.g., `compute_de`, `plot_volcano`, `load_samplesheet`).
- **Branches:** `feat/<short>`, `fix/<short>`, `refactor/<short>`.
- **Commits:** Conventional Commits — `feat(scope): summary` in imperative voice, ≤72 chars.

**Documentation:**

- **Every new function in `R/` gets a `roxygen2` docstring.** At minimum: a one-line title, `@param` for each argument, `@return`. Examples (`@examples`) are appreciated but not required.
   ```r
   #' Load and validate a sample sheet
   #'
   #' @param path Path to the sample sheet (TSV).
   #' @param required_cols Character vector of columns that must be present.
   #' @return A tibble with one row per sample.
   load_samplesheet <- function(path, required_cols) { ... }
   ```
- **Inline comments explain *why*, not *what*.** The code already shows what it does. Comments earn their keep by explaining intent or non-obvious tradeoffs.
- **Don't leave commented-out code.** If it's not used, delete it. Git remembers.

---

## Function discipline

- **Pure where possible** — same inputs → same outputs.
- **No `setwd()`**, no hard-coded absolute paths.
- **No `library()` calls inside functions.** Use `pkg::fn()` for external calls, or declare imports in `DESCRIPTION` if this is a package.
- **Side effects (writing files, plotting to disk) go in their own targets**, separated from computation.
- **Match the dependency family already used in the codebase.** If the existing code uses `dplyr`, don't introduce `data.table` for a new helper, and vice versa. Same for `purrr` vs `lapply`, `cli` vs `message()`, etc. Consistency beats your personal preference. If the codebase is mixed, ask the user which family to use for new code.

---

## Reproducibility

This is a scientific pipeline. Results must be reproducible across machines and across time.

- **Every function with a random component must take a seed (or use `withr::with_seed()`).** This includes: sampling, permutations, stochastic ML (`set.seed()` before model training), bootstrap, random projections, anything in `sample()`, `runif()`, `rnorm()`, `rbinom()`.
   ```r
   # ✅ Good
   compute_perm_pval <- function(x, n_perm = 1000, seed = 42) {
     withr::with_seed(seed, { ... })
   }

   # ❌ Bad
   compute_perm_pval <- function(x, n_perm = 1000) {
     # uses global RNG state — fragile
   }
   ```
- **Tests with randomness must be seeded.** A test that passes 95% of the time is a test that wastes everyone's morning.
- **Don't rely on row ordering** unless the function explicitly guarantees it. Always order explicitly when order matters (e.g., `arrange()` / `setkey()`).
- **Don't use `Sys.time()` or `Sys.Date()`** as default inputs to functions — they make the same call return different results tomorrow.

---

## Do NOT touch

These are non-negotiable. If your task seems to require touching one of these, stop and ask the user.

**Files & folders:**
- **The Shiny app** (e.g. `app/`, `shiny/`, `inst/shiny/`, or wherever it lives). The Shiny code has its own owner and conventions. Never edit it autonomously.
- **`_targets/` folder** — the pipeline cache. It's machine-specific and gitignored.
- **`.Renviron`, `.env`, `*.key`, `*.pem`** — these contain secrets. Never read them, never echo their contents, never write to them.
- **`renv.lock`** — only the user decides when dependencies change.
- **`.github/CODEOWNERS`, `.github/workflows/`** — CI and review rules. Changes require explicit ask.

**Destructive operations — never run without explicit user confirmation in the chat:**
- `git push --force` / `git push -f` (on any branch)
- `git reset --hard` on `main`
- `git rebase -i` that rewrites already-pushed commits
- `rm -rf` on anything outside the current working tree
- `targets::tar_destroy()` (deletes the entire pipeline cache)
- `renv::clean()` / `renv::deactivate()`
- `unlink()` on a directory recursively
- Any command that drops a database table or schema

If you think one of these is needed, *propose* it in chat and wait for an explicit "yes, do that."

---

## Data safety

We work with biological data. Some of it may be unpublished, patient-associated, or regulated.

- **Never copy data file contents into commit messages, PR descriptions, comments, issue threads, or any text you produce.** This includes sample IDs that could re-identify a subject, raw expression values, sequencing reads, anything from `data/` or `data-raw/`.
- **In examples and docstrings, use synthetic data only.** Don't paste real values from the user's files. If you need to illustrate behavior, generate something like `data.frame(sample = c("S1", "S2"), counts = c(100, 200))`.
- **Don't `cat()` / `print()` / `head()` data into the chat to "understand its structure."** Use `str()`, `glimpse()`, `colnames()`, `dim()` — they describe shape without exposing values. If you really need to see values, ask the user to share a redacted excerpt.
- **Don't read large raw files** (FASTQ, BAM, large H5) "to understand the format." Read 5 lines, or read the docs of the package that parses them.

---

## Working with `renv`

This project uses `{renv}` for reproducible dependencies. The rules:

- **After pulling changes from `main`**, the user runs `renv::restore()` to sync their local library with `renv.lock`. You don't need to do this — but be aware their environment may be in an inconsistent state until they do.
- **Never run `renv::snapshot()` autonomously.** That updates `renv.lock`, which is a deliberate human decision.
- **If you need a new package**, say so clearly. Don't just `install.packages()` it and move on.

---

## Working with `{targets}`

- **`_targets.R` contains targets only** — no helper functions, no business logic. Helpers live in `R/`.
- **Don't run `tar_make()` autonomously** unless the user asked. The pipeline may take a long time or rebuild expensive computations.
- Before declaring a change done, suggest the user run:
  - `targets::tar_validate()` — does the pipeline parse?
  - `targets::tar_outdated()` — what would re-run? Surprises mean you touched something deeper than intended.

---

## When to stop and ask — and when to just do

Asking too much is annoying. Asking too little causes drift. Use these thresholds.

**Stop and ask before:**

- Touching more than 3 files in a single change
- Adding a new dependency (any package not already in `renv.lock`)
- Introducing a new design pattern not seen elsewhere in the repo
- Renaming or moving an existing function that other code may call
- Doing anything the user didn't explicitly ask for (even if it seems helpful)
- Any item from the "Do NOT touch" list
- Any destructive operation

**Don't ask — just do — for:**

- Typos, formatting fixes, obvious bugs in the code you're already editing
- Adding a test for a function you just wrote
- Matching established patterns in adjacent code
- Adding `roxygen2` docstrings to functions you're touching
- Removing your own debug `print()` statements before declaring done

When in doubt — err on asking. A 30-second clarifying question saves 30 minutes of cleanup.

---

## Where you SHOULD be creative

The rules above protect our structure. But there are places we *want* fresh thinking — please bring it.

- **Algorithm choice and statistical approach.** If you can suggest a better method for what we're trying to compute, say so before writing code. We'd rather hear "have you considered X?" than discover we used a suboptimal approach a month later.
- **Refactoring proposals.** If you notice existing code is awkward, duplicative, or fragile — *propose* a refactor. Don't just do it. A short note in chat ("I noticed `compute_de` and `compute_dea` share 80% of their logic — want me to consolidate?") is welcome.
- **Test design.** What edge cases should we cover? NA handling, empty inputs, single-row inputs, very large inputs, malformed sample sheets? You think about this better than we do — bring suggestions.
- **Naming.** Within our conventions (snake_case, verb_noun), you have full freedom. Pick names that read well in context.
- **Plot design.** Within `ggplot2`, you decide aesthetics. Just match the existing look of plots in the repo if there are any.
- **Error messages.** When you write a `stop()` or `cli::cli_abort()`, write the message you'd want to read at 11 PM debugging — actionable, specific, kind.
- **Documentation prose.** When you write a docstring or comment, you have full creative freedom. Be clear, be a little human.

**Mode of creativity:** *propose, don't impose.* Surface ideas to the user; let them decide. If they say "go ahead," go ahead.

---

## When the user asks you to break a rule in this file

Sometimes the user will ask you to do something this file says not to do — touch the Shiny app, bump `renv.lock`, create a new directory, run a destructive operation. **Their call to make.** When that happens:

1. **Acknowledge** that you're stepping outside the convention: *"This file usually says don't touch X."*
2. **Ask once** for confirmation: *"Are you sure you want me to do this?"*
3. If they confirm, **proceed**.
4. **Note in the commit message that this is intentional**: e.g., `chore: bump renv.lock — intentional, see issue #42`.

This protects the user from autopilot mistakes while not blocking them when they genuinely want to override.

---

## Before you tell the user you're done

Mentally walk through:

- [ ] Did you add to existing files where possible, instead of creating new ones?
- [ ] Did you invent any function signatures, or do they match what exists in `R/`?
- [ ] Did you avoid changing `renv.lock`?
- [ ] Did you avoid touching the Shiny app?
- [ ] Did you remove debug `print()` calls and commented-out code?
- [ ] If you added a function, did you also add a `roxygen2` docstring and a test?
- [ ] If anything is stochastic, is it seeded?
- [ ] Is the diff small and focused on one thing?

---

## Writing the PR description

When you've finished a change and the user asks for a PR description, use this format:

```markdown
## What
<1–3 sentences: what does this PR do?>

## Why
Closes #<issue number>

## How to test
\```r
# commands the reviewer can paste
\```

## Checklist
- [ ] `targets::tar_validate()` passes
- [ ] Self-reviewed in `git diff`
- [ ] No `renv.lock` changes
- [ ] No Shiny app changes
- [ ] New functions have docstrings + tests
- [ ] Random components are seeded
```

---

## How the team works

- Small, single-concern PRs. If your PR description has the words "and also," split the work.
- The architect (Michal) reviews changes to `_targets.R`, `R/`, and `renv.lock`.
- We use **`git switch`** for branches, not `git checkout`.
- We sync with main using **`git fetch origin` + `git rebase origin/main`** (not `git pull`).

---

## The bottom line

You are a teammate, not a tool. Protect the structure we've built; bring fresh ideas where they're welcome; ask before you wander outside the lines. When the user says "go," go.
