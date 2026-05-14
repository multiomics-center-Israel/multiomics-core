# Our Git Cookbook

A practical, command-by-command guide for our 5-person team. No theory. Every section answers three questions: **when** to do something, **why** we do it, and **what to type**.

> **Audience.** You know `git clone`, `git commit`, and `git push`. Branches, PRs, and rebase still feel a bit fuzzy. That's exactly who this is written for.

---

## 1. Why we need this

Right now we have four problems hurting us:

- **Architecture drift** — features land in places they don't belong.
- **Duplicated and orphan code** — the same helper exists 3 times; some functions are never called.
- **Hard-to-review PRs** — too many files and too many concerns at once.
- **Claude Code surprises** — generated code that ignores the structure we already have.

The goal of this guide is simple: **`main` is always trustworthy, reviews take minutes not hours, and nobody is afraid to push.** Follow the recipes below and we get there.

---

## 2. Words you need to know (90 seconds)

| Word | What it means |
|---|---|
| **`main`** | The "official" branch. The clean, working version of the pipeline. **We never edit it directly.** |
| **branch** | A separate line of work. Your sandbox where you make changes without touching `main`. |
| **`origin`** | A nickname for "the version on GitHub" (as opposed to your local copy). |
| **commit** | A saved snapshot of your changes, with a message. Like "Save As" with a note. |
| **push** | Send your local commits up to GitHub. |
| **pull** | Bring down what's on GitHub into your local copy. |
| **PR (Pull Request)** | A GitHub page where you propose: "merge my branch into `main`." Where review happens. |
| **rebase** | "Replay my work on top of the latest `main`." Keeps your branch up-to-date with the team's progress. |
| **conflict** | When git can't auto-merge two changes to the same lines. You decide what wins. |

That's the whole vocabulary. If a word below is fuzzy, come back here.

---

## 3. Recipe: Starting a new branch

**When:** You have a new feature, a bug fix, or any change to make. Always start here. Never edit `main` directly.

**Why:** `main` must always work. We isolate every change on a branch. That way, if something goes wrong, `main` is untouched and the rest of the team is not blocked.

**The commands:**

```bash
# 1. Move to the main branch
git switch main

# 2. Pull the latest version from GitHub
git pull origin main

# 3. Create your new branch and move to it
#    Use a clear name: feat/<short-thing>, fix/<short-thing>, refactor/<short-thing>
git switch -c feat/add-deseq2-target
```

**What just happened:**
- `switch main` → you're now standing on `main` locally.
- `pull origin main` → you downloaded any commits the team merged today. Your `main` is current.
- `switch -c feat/...` → you created a new branch from `main` and switched to it. The `-c` means "create." Anything you do now stays on this branch.

> **Note on `git switch` vs `git checkout`:** You may see `git checkout` in older tutorials. They do the same thing for branches, but `git switch` is newer and clearer (it only switches branches; `checkout` does many other things too). We use `git switch`.

**Branch naming:** start with a prefix.

| Prefix | Use it for | Example |
|---|---|---|
| `feat/` | new feature, new target, new function | `feat/add-deseq2-target` |
| `fix/` | a bug fix | `fix/sample-sheet-na` |
| `refactor/` | restructuring without changing behavior | `refactor/split-utils` |
| `docs/` | docs / comments only | `docs/onboarding` |
| `chore/` | deps, CI, formatting | `chore/bump-renv` |

---

## 4. Recipe: Continuing a branch the next day (fetch + rebase)

**When:** You opened a branch yesterday. Today you want to keep working on it.

**Why:** While you were away, the team probably merged things into `main`. If you ignore that and keep working, you'll pile up a huge mess later. We sync **your branch** with the latest `main` *now*, so the only conflicts you ever face are tiny ones from one day of changes.

**The commands:**

```bash
# 1. Move to your branch
git switch feat/add-deseq2-target

# 2. Download the latest main from GitHub (this does NOT change your files yet)
git fetch origin

# 3. Replay your work on top of the latest main
git rebase origin/main
```

**What just happened (in plain English):**
- `git fetch origin` → goes to GitHub, downloads any new commits, and stores them under the name `origin/main`. Your actual branch and your files are untouched. Think of it as "let me see what's new, without applying anything."
- `git rebase origin/main` → takes your branch's commits and **replays them on top** of what you just fetched. Now your branch sits on the latest `main`, as if you'd started it just now.

Separating these two steps is on purpose: `fetch` is *safe* (just downloading), `rebase` is the one that actually changes your branch. Once you've done this for a week, you'll understand exactly what each command does.

> **The shortcut you'll see later:** `git pull --rebase origin main` does steps 2 and 3 in one command. It's the same thing — just fused together. Use it once you feel comfortable with the two-step version.

**Visual:**

```
Before pull --rebase:
  main:        A───B───C───D     (D = newest)
                    \
  your branch:       E───F        (your commits)

After pull --rebase origin main:
  main:        A───B───C───D
                            \
  your branch:               E'──F'   (same changes, now on top of D)
```

**If git stops and says "CONFLICT":** see Recipe 7. Don't panic — there's an undo button.

**Do this every morning when you sit down to work, before doing anything else.** It takes 5 seconds when there are no conflicts, and saves hours later.

---

## 5. Recipe: Saving your work (commit & push)

**When:** You made changes you want to save. Do this often — small commits, several times a day.

**Why:** A commit is a checkpoint. If you break things later, you can return to the last good commit. Pushing puts your commits on GitHub so they're backed up and visible.

**The commands:**

```bash
# 1. See what changed
git status

# 2. Stage the files you want in this commit
git add R/de.R R/io.R         # specific files
# — or —
git add .                      # everything that changed

# 3. Save the snapshot with a clear message
git commit -m "feat(targets): add deseq2 target"

# 4. Send to GitHub
git push
```

**The first time you push a new branch**, git will ask for a slightly different command. Just paste what it tells you:

```bash
git push --set-upstream origin feat/add-deseq2-target
```

After the first push, plain `git push` works.

**Commit message format we use:**

```
<type>(<scope>): <short summary in imperative voice>
```

Types: `feat`, `fix`, `refactor`, `docs`, `chore`.
Scope is the area of the code: `targets`, `qc`, `de`, `R`, etc.

Good examples:
```
feat(targets): add deseq2 target
fix(qc): handle NA values in sample sheet parser
refactor(R): split utils.R into io.R and stats.R
docs(readme): document input data layout
```

**Rule of thumb:** if the diff is over ~300 lines or your commit message says "and also…", split it into multiple commits.

---

## 6. Recipe: Opening a Pull Request and responding to review

**When:** Your branch is ready for someone to look at it (or you want early feedback even though it's not done).

**Why:** PRs are how every change gets into `main`. Nothing reaches `main` without going through a PR.

**The steps:**

```bash
# 1. Make sure your branch is up to date with main, then push
git fetch origin
git rebase origin/main
git push
```

Then on GitHub:

1. Open the repo. You'll see a yellow banner: **"Compare & pull request"**. Click it.
2. **Title:** use the same format as a commit message — `feat(targets): add deseq2 target`.
3. **Description:** fill in this template (we keep this in the repo as `.github/pull_request_template.md`):
   ```markdown
   ## What
   <one or two sentences: what does this PR do?>

   ## Why
   Closes #<issue number>

   ## How to test
   <commands the reviewer can paste>
   ```
4. **If it's not ready for review yet:** click the dropdown next to "Create pull request" and choose **Create draft pull request**. Drafts let you push commits and run CI without bothering reviewers.
5. **When ready:** click "Ready for review" and request review from Michal.

**While the PR is open, when comments come in:**

```bash
# 1. Make the requested changes locally
# 2. Commit them
git add .
git commit -m "fix: address review comments — handle NA case"

# 3. Push (this updates the PR automatically)
git push
```

Then on GitHub:
- Reply to every comment. Even "done" or 👍 closes the loop.
- Click **Re-request review** so the reviewer knows it's their turn again.

**After your PR is merged:**

```bash
git switch main
git pull origin main
git branch -d feat/add-deseq2-target   # delete your local branch — it's done
```

---

## 7. Recipe: Working with Claude Code on shared code

**When:** Before you ask Claude to write or change anything in the repo.

**Why:** If Claude does something we don't like, we want a clean undo. And Claude has no idea how *our* codebase is structured unless we tell it — that's why we end up with duplicate helpers and orphan files.

**The commands and habits:**

**Step 1 — Always commit before a Claude session.**

```bash
git status                                  # check what's uncommitted
git add .
git commit -m "wip: save before claude session"
```

If Claude makes a mess, you can throw it away cleanly:

```bash
git restore .                  # discard unstaged changes
git reset --hard HEAD          # also discard staged changes
```

**Step 2 — Tell Claude what to read first.** Don't just say *"add a target that runs DESeq2."* Say:

> *"Read `R/io.R` and `R/stats.R` first. Notice the helper pattern. Then add a target that runs DESeq2, following the same pattern, in the right file."*

**Step 3 — Make Claude search before writing.** Before adding new code, ask:

> *"Search the repo for an existing function that already does this. List what you found. Justify why a new one is needed."*

This single prompt eliminates most duplication.

**Step 4 — Read every line Claude writes.**

```bash
git diff           # show what changed since the last commit
```

Look for:
- Did it create new files when it should have edited existing ones?
- Did it add packages to `renv.lock` we didn't ask for?
- Did it leave debug `print()` calls or commented-out code?

If you don't like what you see, tell Claude what to fix. If it's beyond saving:

```bash
git restore .
```

You're back where you started.

**Step 5 — Have Claude write the PR description.** *"Write a PR description for this branch in our template format."* It does this well, and forces it to articulate what changed.

---

## 8. Recipe: Resolving conflicts during rebase

**When:** You ran `git rebase origin/main` and got `CONFLICT (content): Merge conflict in path/to/file.R`. This means someone else changed the same lines you did.

**Why:** Git can't auto-decide which version wins. You have to tell it.

**The commands:**

```bash
# 1. See which files conflict
git status
```

Open each conflicting file. You'll see this in the file:

```r
<<<<<<< HEAD
# the version from main
result <- compute_de(counts, design)
=======
# your version
result <- run_deseq(counts, design, formula)
>>>>>>> your-branch
```

```bash
# 2. Edit the file:
#    - Pick one version, OR combine them
#    - Delete all three marker lines (<<<<<<<, =======, >>>>>>>)
#    - Save the file

# 3. Tell git this file is resolved
git add path/to/file.R

# 4. Continue the rebase
git rebase --continue
```

If there are more conflicts, git will stop again. Repeat steps 2–4 for each file.

**The escape hatch:** if you're stuck or scared, you can always undo:

```bash
git rebase --abort
```

This puts you back **exactly** where you were before the rebase started. No harm done. Then come find Michal.

---

## 9. Recipe: Common mistakes and how to undo them

### "I committed to `main` by mistake"

```bash
# Save your work to a new branch
git branch feat/my-fix

# Reset main back to what's on GitHub (your commit moves with the branch)
git reset --hard origin/main

# Switch to your new branch
git switch feat/my-fix
```

### "Claude rewrote things and I want to undo"

```bash
git diff                # see what changed
git restore .           # throw away unstaged changes
git reset --hard HEAD   # also throw away staged changes
```

This is why we always commit *before* a Claude session.

### "I committed a secret (token, password, key)"

**Stop. Do not push.** Then:

1. **Rotate the secret immediately** — assume it's already public.
2. Tell the team in the channel.
3. Come find Michal — we'll scrub it from history together.

> Removing it from your latest commit is **not enough**. The secret is still in git history and is still leaked. Always rotate first, scrub second.

### "I lost my work"

You probably didn't. Try:

```bash
git reflog              # shows every state your repo was recently in
```

You'll see a list. Find the line that looks like the moment before things went wrong. Then:

```bash
git reset --hard HEAD@{3}    # use the number from reflog
```

If you're not sure — **stop and ask before running this**. `--hard` cannot be undone easily.

### "My push was rejected"

Probably someone pushed to your branch (or to `main`) since you last pulled. Sync first, then push:

```bash
git fetch origin
git rebase origin/<your-branch>
git push
```

---

## 10. The everyday loop (your daily flow)

This is what your day looks like:

```bash
# ── MORNING (sit down to work) ──────────────────────────────
git switch feat/my-thing                      # if continuing a branch
git fetch origin                              # download latest main from GitHub
git rebase origin/main                        # replay your work on top

# ── DURING THE DAY (every hour or two) ──────────────────────
git status                                    # what changed?
git add .
git commit -m "feat(scope): clear message"

# ── BEFORE A CLAUDE SESSION ────────────────────────────────
git add . && git commit -m "wip: save before claude"

# ── END OF DAY ─────────────────────────────────────────────
git push                                      # back up to GitHub

# ── WHEN READY FOR REVIEW ──────────────────────────────────
git fetch origin                              # sync one more time
git rebase origin/main
git push
# → open PR on GitHub, request Michal's review

# ── AFTER MERGE ────────────────────────────────────────────
git switch main
git pull origin main
git branch -d feat/my-thing                   # delete the local branch
```

**Print this section.** Stick it next to your screen for the first two weeks.

---

## 11. Quick reference card

```bash
# Where am I, and what's changed?
git status
git branch                    # which branch am I on?
git log --oneline -10         # last 10 commits

# Sync with the team (two-step version)
git fetch origin              # download latest from GitHub
git rebase origin/main        # replay your work on top

# Save work
git add <file>     OR     git add .
git commit -m "type(scope): message"
git push

# Start something new
git switch main
git pull origin main
git switch -c feat/my-thing

# Look at what changed
git diff                      # unstaged changes
git diff --staged             # staged but not committed
git diff main                 # everything different from main

# Undo
git restore <file>            # discard changes to a file
git restore .                 # discard ALL unstaged changes
git reset --hard HEAD         # nuclear undo of uncommitted work
git rebase --abort            # bail out of a rebase

# In trouble?
git reflog                    # see recent history of your repo
# — and ask Michal
```

---

*If a recipe here is unclear, raise it in the next standup. We'll fix this guide.*
