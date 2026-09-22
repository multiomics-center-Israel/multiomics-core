# tests/testthat/test-report-explorer-missing.R
#
# Two defects and one robustness change in the proteomics report, all in code
# that only exists inside report_template_proteomics.Rmd.
#
# The explorer reads the measured `<sample>` block of final_results.tsv (not the
# imputed `<sample>.norm` block), so its cells are NA in R wherever a protein
# was not measured. Those reach the browser as the string "NA" -- that is how
# jsonlite writes NA in a numeric vector -- which matters for what each site
# had to do about them.
#
# 1. A group with no measured values disappeared. Its trace was built from every
#    sample in the group, missing ones included, so it carried only "NA"
#    placeholders. Plotly does not plot non-numeric values, and a box trace with
#    no plottable point contributes no category, so the group left the x axis
#    instead of being shown as unmeasured. Note that stating categoryarray is
#    not a fix: it orders the categories the traces created, and an empty group
#    creates none. The axis is labelled from an explicit tick array instead,
#    which does not consult the data at all.
#
# 2. ORA section headings stripped `_up`/`_down` along with the file suffix, so
#    a database with both files produced two sections headed identically.
#
# 3. What counts as a measured value was decided separately at each read site.
#    The summary line's `!isNaN(v)` did reject "NA" correctly, so its mean and
#    sample count were right; but it is not the same test the trace loops now
#    need, and it accepts null, which Number() turns into a measured zero. One
#    shared bpIsMeasured() covers the representation actually emitted and the
#    ones a serialisation change could produce.
#
# These are assertions against the template source, not behavioural tests. The
# first two live in browser JavaScript that this suite cannot execute, and the
# third runs only inside a knitr chunk during a render. Source-level pinning is
# the same idiom used by test-imputation-seeding.R for the preprocessing seed.
# The render checklist in the PR is what actually verifies the output.

# Locate a repo file from either the repo root or tests/testthat.
repo_file <- function(...) {
    candidates <- c(testthat::test_path("..", "..", ...), file.path(...))
    candidates[file.exists(candidates)][1]
}

report_src <- function() {
    f <- repo_file("R", "domain", "proteomics", "report_template_proteomics.Rmd")
    skip_if(is.na(f) || !file.exists(f), "proteomics report template not found")
    paste(readLines(f, warn = FALSE), collapse = "\n")
}


# =============================================================================
# 1. One shared test for what counts as a measured value
# =============================================================================

test_that("the protein explorer has one guard for what counts as measured", {
    src <- report_src()

    expect_match(src, "function bpIsMeasured(v)", fixed = TRUE)
    # isNaN() is what rejects the "NA" string the template actually emits.
    expect_match(src, "!isNaN(v)", fixed = TRUE)
    # null is the case isNaN() would accept, since Number(null) is 0.
    expect_match(src, "v !== null", fixed = TRUE)
    expect_match(src, "v !== undefined", fixed = TRUE)
})

test_that("the summary line uses that guard rather than its own test", {
    src <- report_src()

    # The old inline filter rejected "NA" correctly, so its numbers were right.
    # It is replaced so that every read site decides this the same way, not
    # because it was producing a wrong mean.
    expect_false(grepl("filter(function(v) { return !isNaN(v); })", src, fixed = TRUE))
    expect_match(src, "values.filter(bpIsMeasured)", fixed = TRUE)
    # Say "measured samples", so a reader can tell that a protein with gaps is
    # summarised over fewer samples than the run has.
    expect_match(src, "measured samples", fixed = TRUE)
    # A protein with nothing measured reads as such rather than as a NaN mean.
    expect_match(src, "no measured values", fixed = TRUE)
})

test_that("both explorer traces drop missing values with that same guard", {
    src <- report_src()

    # This is the defect: the trace loops pushed every sample in the group,
    # missing ones included, so the trace carried unplottable placeholders.
    # Value and sample name stay in lockstep, because the hover text is
    # positional and would otherwise mislabel points.
    expect_match(src, "if (groups[i] === g && bpIsMeasured(values[i])) {", fixed = TRUE)
    # Every use: the definition, the two trace loops, the empty-group scan, and
    # the info-line filter -- which names the function without calling it, so
    # this counts the bare name. Catches a guard being dropped again without
    # blocking a new one being added.
    n_uses <- lengths(regmatches(src, gregexpr("bpIsMeasured", src, fixed = TRUE)))
    expect_gte(n_uses, 5L)
})

test_that("the peptide explorer keeps the guard this fix was modelled on", {
    # The protein explorer now matches what the peptide explorer already did.
    # If that reference is ever weakened, this fix loses its precedent.
    expect_match(report_src(), "pep.expr[i] !== null && !isNaN(pep.expr[i])", fixed = TRUE)
})


# =============================================================================
# 2. A group with nothing measured stays visible
# =============================================================================

test_that("every group keeps a labelled slot on the axis", {
    src <- report_src()

    # Ticks, not categories. A category axis only knows the categories its own
    # traces supply, so a group with nothing measured never creates one --
    # categoryarray orders the categories that exist, it does not bring a
    # missing one into being, which a rendered report established the hard way.
    # tickmode "array" draws exactly the ticks listed, with no dependency on
    # what the data happens to contain.
    expect_match(src, "tickmode: \"array\"", fixed = TRUE)
    expect_match(src, "ticktext: uniqueGroups", fixed = TRUE)
    # Stated for the same reason: an empty group at either end must not be
    # autoranged off the plot.
    expect_match(src, "range: [-0.5, uniqueGroups.length - 0.5]", fixed = TRUE)
    # The mechanism that did not work must not creep back as the one that does.
    expect_false(grepl("categoryarray: uniqueGroups", src, fixed = TRUE))
})

test_that("each box sits in its own group's slot", {
    src <- report_src()

    # The slot is the group's index in uniqueGroups, which is both what keeps
    # the order and what stops the groups after an unmeasured one from shifting
    # left into the gap it leaves.
    expect_match(src, "x: gValues.map(function() { return gi; })", fixed = TRUE)
    expect_match(src, "tickvals: uniqueGroups.map(function(g, i) { return i; })", fixed = TRUE)
})

test_that("the multi-protein comparison keeps a slot for every selection", {
    # The same defect one branch down: a protein measured in no sample
    # contributes no point, so naming its x position left it without one and it
    # disappeared from the comparison the moment a second protein was selected.
    # Reachable because filtering$min_count has no lower bound -- 0 keeps every
    # feature, all-missing included -- so the two branches must agree.
    src <- report_src()

    expect_match(src, "tickvals: protNames.map(function(p, i) { return i; })", fixed = TRUE)
    expect_match(src, "ticktext: protNames", fixed = TRUE)
    expect_match(src, "range: [-0.5, protNames.length - 0.5]", fixed = TRUE)
    # The position is the protein's index in the selection, so the slots stay
    # put whether or not each protein has anything to show.
    expect_match(src, "selectedGenes.forEach(function(sg, pi) {", fixed = TRUE)
    expect_match(src, "xVals.push(pi);", fixed = TRUE)
    # The superseded mechanism must not creep back here either.
    expect_false(grepl("categoryarray: protNames", src, fixed = TRUE))
})

test_that("the multi-protein hover still names the protein", {
    # x carries a position now, so %{x} would render a number. The name travels
    # in customdata instead, alongside the sample name in text, both pushed in
    # lockstep with the value so the hover cannot drift off its point.
    src <- report_src()

    expect_match(src, "customdata: protTexts", fixed = TRUE)
    expect_match(src, "protTexts.push(sg.name);", fixed = TRUE)
    expect_match(src, "hovertemplate: \"%{customdata}<br>Sample: %{text}", fixed = TRUE)
    expect_false(grepl("hovertemplate: \"%{x}<br>Sample: %{text}", src, fixed = TRUE))
})

test_that("all R chunks in the proteomics report parse", {
    # The durable guard behind the apostrophe scanner below. CI can otherwise
    # pass on a template whose R does not parse, because nothing in this suite
    # renders it -- which is exactly how two heads on this branch went green
    # while the report could not have knitted.
    #
    # Deliberately not knitr::purl(): tangling resolves the `purl`, `eval` and
    # `child` chunk options, and a chunk whose option cannot be determined is
    # assigned purl = FALSE and dropped from the tangled script. The explorer
    # chunks here are guarded by `eval=show_protein_explorer && explorer_ready`,
    # so the very chunks worth checking are the ones that would be omitted, and
    # a broken template could tangle to something that parses cleanly. Reading
    # the fences directly evaluates nothing at all.
    f <- repo_file("R", "domain", "proteomics", "report_template_proteomics.Rmd")
    skip_if(is.na(f) || !file.exists(f), "proteomics report template not found")

    lines <- readLines(f, warn = FALSE)
    i <- 1L
    n_chunks <- 0L

    while (i <= length(lines)) {
        if (!grepl("^```\\{r(?:[ ,}]|$)", lines[i], perl = TRUE)) {
            i <- i + 1L
            next
        }

        n_chunks <- n_chunks + 1L
        start <- i

        rest <- if (start < length(lines)) {
            lines[(start + 1L):length(lines)]
        } else {
            character(0)
        }
        rel_end <- which(grepl("^```\\s*$", rest, perl = TRUE))[1L]

        if (is.na(rel_end)) {
            fail(sprintf("Unclosed R chunk at line %d: %s", start, lines[start]))
            break
        }

        end <- start + rel_end
        body <- if (end > start + 1L) lines[(start + 1L):(end - 1L)] else character(0)

        err <- tryCatch({ parse(text = body); NULL }, error = function(e) e)
        if (!is.null(err)) {
            fail(sprintf("R chunk starting at line %d does not parse (%s): %s",
                         start, lines[start], conditionMessage(err)))
        }

        i <- end + 1L
    }

    # An order-of-magnitude floor, not an exact count: a walker that stopped
    # after the first few chunks would otherwise pass by checking almost
    # nothing, while adding or removing a chunk should not fail the suite.
    expect_gt(n_chunks, 100L)
})

test_that("the explorer script blocks carry no unescaped apostrophe", {
    # Every one of these JS blocks is emitted from inside cat('...'), a
    # single-quoted R string, so one bare apostrophe in the JavaScript closes
    # that string early and breaks the chunk -- a comment reading "the group's
    # slot" is enough, which is exactly how this got in. Nothing in this suite
    # renders the template, so CI cannot see it; it surfaces only when someone
    # knits the report, which is the most expensive place to find out.
    f <- repo_file("R", "domain", "proteomics", "report_template_proteomics.Rmd")
    skip_if(is.na(f) || !file.exists(f), "proteomics report template not found")
    lines <- readLines(f, warn = FALSE)

    opens <- which(trimws(lines) == "cat('")
    closes <- which(trimws(lines) == "')")
    expect_gt(length(opens), 0L)

    offenders <- character(0)
    for (a in opens) {
        b <- closes[closes > a][1]
        if (is.na(b) || b <= a + 1L) next
        body <- lines[(a + 1L):(b - 1L)]
        bare <- grepl("(^|[^\\\\])'", body)
        if (any(bare)) offenders <- c(offenders, trimws(body[bare]))
    }
    expect_equal(offenders, character(0))
})

test_that("the multi-protein grouping and colours are untouched", {
    src <- report_src()

    # boxmode "group" is what offsets the groups within each protein's slot;
    # without it the fix would stack them.
    expect_match(src, "boxmode: \"group\"", fixed = TRUE)
    expect_match(src, "boxgroupgap: 0.1", fixed = TRUE)
})

test_that("a group with nothing measured is named rather than left blank", {
    src <- report_src()

    expect_match(src, "no measured values in: ", fixed = TRUE)
    expect_match(src, "emptyGroups", fixed = TRUE)
})

test_that("labels are escaped before they reach innerHTML", {
    src <- report_src()

    # The info line is assembled as HTML, and both the group labels and the
    # protein name come from the sample sheet and the annotation columns. They
    # are data: a label carrying angle brackets must be shown, not parsed.
    expect_match(src, "function bpEscapeHtml(s)", fixed = TRUE)
    expect_match(src, "bpEscapeHtml(sg.name)", fixed = TRUE)
    expect_match(src, "return bpEscapeHtml(g);", fixed = TRUE)
    # Ampersand first, or the escapes would re-escape each other's output.
    expect_match(src, ".replace(/&/g, \"&amp;\")", fixed = TRUE)
})


# =============================================================================
# 3. ORA headings keep their direction
# =============================================================================

test_that("ORA sections name their direction", {
    src <- report_src()

    expect_match(src, "_ora_up\\\\.csv$", fixed = TRUE)
    expect_match(src, "_ora_down\\\\.csv$", fixed = TRUE)
    expect_match(src, "(up-regulated)", fixed = TRUE)
    expect_match(src, "(down-regulated)", fixed = TRUE)
})

test_that("ORA file discovery and statistics are untouched", {
    src <- report_src()

    # The fix renames sections only: same files read, same significance rule.
    expect_match(src, "pattern = \"pathway_.*_ora.*\\\\.csv$\"", fixed = TRUE)
    expect_match(src, "ora_df[[padj_col]] < 0.05", fixed = TRUE)
})
