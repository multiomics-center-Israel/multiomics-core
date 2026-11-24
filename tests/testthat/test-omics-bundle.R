# Omics bundles (R/core/20_omics_bundle.R): written once, never overwritten,
# verified on read, and invisible to readers until complete.

fixed_provenance <- function(time = "2026-01-31T09:30:00Z") {
    list(created_at_utc = time,
         git_commit = "0123456789abcdef0123456789abcdef01234567",
         git_dirty = FALSE,
         r_version = "4.5.1",
         renv_lock_sha256 = strrep("a", 64),
         packages = list(base = "4.5.1"))
}

example_parts <- function() {
    list(expr = matrix(c(1.5, 2, 3, 4.25), nrow = 2,
                       dimnames = list(c("F1", "F2"), c("S1", "S2"))),
         meta = data.frame(sample = c("S1", "S2"), group = c("A", "B")))
}

example_info <- function(...) {
    utils::modifyList(list(omics = "proteomics", project = "Demo Project",
                           analysis_round = "R1"), list(...))
}

# Bundle files are made read-only; a test that tampers with one has to undo that.
make_writable <- function(path) Sys.chmod(path, mode = "0644")

# Edits manifest.json. With refresh_checksum = TRUE the edit is made consistent,
# which lets a test reach the checks behind the manifest checksum; FALSE leaves
# the edit for the manifest checksum to catch.
rewrite_manifest <- function(bundle, edit, refresh_checksum = TRUE) {
    file <- file.path(bundle, "manifest.json")
    m <- edit(jsonlite::read_json(file, simplifyVector = TRUE))
    make_writable(file)
    jsonlite::write_json(m, file, auto_unbox = TRUE, pretty = TRUE, digits = NA)
    if (refresh_checksum) {
        sums <- file.path(bundle, "checksums.sha256")
        make_writable(sums)
        writeLines(paste0(.bundle_sha256_file(file), "  manifest.json"), sums)
    }
}


test_that("a bundle round-trips its parts exactly", {
    root <- withr::local_tempdir()
    path <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())

    expect_true(dir.exists(path))
    expect_match(basename(path),
                 "^proteomics-demo_project-r1-20260131T093000Z-[0-9a-f]{8}$")
    b <- read_omics_bundle(path)
    expect_identical(b$parts, example_parts())
})

test_that("the manifest keeps the exact identity and the provenance", {
    root <- withr::local_tempdir()
    path <- write_omics_bundle(example_parts(), example_info(instrument = "Orbitrap"),
                               root, fixed_provenance())
    m <- read_omics_bundle(path, parts = "meta")$manifest

    expect_identical(m$status, "complete")
    expect_identical(m$project, "Demo Project")
    expect_identical(m$analysis_round, "R1")
    expect_identical(m$bundle_id, basename(path))
    expect_identical(m$info$instrument, "Orbitrap")
    expect_identical(m$provenance$git_commit, fixed_provenance()$git_commit)
    expect_setequal(names(m$parts), c("expr", "meta"))
    expect_match(m$parts$expr$sha256, "^[0-9a-f]{64}$")
})

test_that("only the requested parts are loaded, and unknown parts are named", {
    root <- withr::local_tempdir()
    path <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())

    b <- read_omics_bundle(path, parts = "meta")
    expect_named(b$parts, "meta")
    expect_error(read_omics_bundle(path, parts = "nope"),
                 "has no part\\(s\\) nope\\. Available: expr, meta")
})

test_that("an existing bundle ID is refused, never overwritten", {
    root <- withr::local_tempdir()
    path <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())
    before <- read_omics_bundle(path)$manifest$content_hash

    # Same content written at the same recorded time resolves to the same ID.
    expect_error(write_omics_bundle(example_parts(), example_info(), root, fixed_provenance()),
                 "already exists.*never overwritten")
    expect_identical(read_omics_bundle(path)$manifest$content_hash, before)
    expect_equal(nrow(list_omics_bundles(root)), 1L)
})

test_that("an existing empty folder at the destination is refused, not replaced", {
    # The ID depends only on content, identity and recorded time, so a bundle
    # written elsewhere first tells the test which folder to occupy.
    id <- basename(write_omics_bundle(example_parts(), example_info(),
                                      withr::local_tempdir(), fixed_provenance()))
    root <- withr::local_tempdir()
    occupied <- file.path(root, "bundles", "proteomics", id)
    dir.create(occupied, recursive = TRUE)

    expect_error(write_omics_bundle(example_parts(), example_info(), root, fixed_provenance()),
                 "already exists; bundles are never overwritten")
    expect_length(list.files(occupied, all.files = TRUE, no.. = TRUE), 0L)
    expect_equal(nrow(list_omics_bundles(root)), 0L)
})

test_that("of two concurrent writers of one bundle, exactly one succeeds", {
    skip_on_os("windows")  # mclapply needs fork()
    root <- withr::local_tempdir()

    results <- parallel::mclapply(1:2, function(i) {
        tryCatch(list(ok = TRUE, value = write_omics_bundle(example_parts(), example_info(),
                                                            root, fixed_provenance())),
                 error = function(e) list(ok = FALSE, value = conditionMessage(e)))
    }, mc.cores = 2)
    ok <- vapply(results, `[[`, logical(1), "ok")

    expect_equal(sum(ok), 1L)
    expect_match(results[[which(!ok)]]$value, "already exists")
    expect_identical(read_omics_bundle(results[[which(ok)]]$value)$parts, example_parts())
    expect_equal(nrow(list_omics_bundles(root)), 1L)
})

test_that("a successful write leaves no staging folder behind", {
    root <- withr::local_tempdir()
    path <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())

    expect_identical(list.files(dirname(path), all.files = TRUE, no.. = TRUE),
                     basename(path))
})

test_that("consecutive writes with one recorded time both succeed", {
    # Same process, same second, different content: each write stages under
    # its own folder, so the first cannot block the second.
    root <- withr::local_tempdir()
    changed <- example_parts()
    changed$expr[1, 1] <- 9

    a <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())
    b <- write_omics_bundle(changed, example_info(), root, fixed_provenance())

    expect_false(identical(a, b))
    expect_identical(read_omics_bundle(b)$parts, changed)
    expect_setequal(list.files(dirname(a), all.files = TRUE, no.. = TRUE),
                    basename(c(a, b)))
})

test_that("new content or a new time gives a new bundle beside the old one", {
    root <- withr::local_tempdir()
    a <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())
    changed <- example_parts()
    changed$expr[1, 1] <- 9
    b <- write_omics_bundle(changed, example_info(), root, fixed_provenance())
    c <- write_omics_bundle(example_parts(), example_info(), root,
                            fixed_provenance("2026-02-01T10:00:00Z"))

    expect_length(unique(basename(c(a, b, c))), 3L)
    expect_equal(nrow(list_omics_bundles(root)), 3L)
    expect_identical(read_omics_bundle(a)$parts, example_parts())
})

test_that("a part changed after writing fails validation without showing its content", {
    root <- withr::local_tempdir()
    parts <- c(example_parts(), list(note = "SECRET-VALUE-42"))
    path <- write_omics_bundle(parts, example_info(), root, fixed_provenance())

    file <- file.path(path, "parts", "note.rds")
    make_writable(file)
    saveRDS("SECRET-VALUE-43", file)

    err <- expect_error(validate_omics_bundle(path),
                        "part 'note' does not match its checksum")
    expect_false(grepl("SECRET", conditionMessage(err)))
    expect_error(read_omics_bundle(path, parts = "note"), "does not match its checksum")
    # Untouched parts still read.
    expect_identical(read_omics_bundle(path, parts = "meta")$parts$meta, parts$meta)
})

test_that("an unchanged bundle passes the manifest check", {
    root <- withr::local_tempdir()
    path <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())

    sums <- readLines(file.path(path, "checksums.sha256"))
    expect_identical(sums, paste0(.bundle_sha256_file(file.path(path, "manifest.json")),
                                  "  manifest.json"))
    expect_identical(validate_omics_bundle(path)$bundle_id, basename(path))
})

test_that("an edited project name or provenance fails the manifest check", {
    root <- withr::local_tempdir()

    p1 <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())
    rewrite_manifest(p1, function(m) { m$project <- "Other Project"; m },
                     refresh_checksum = FALSE)
    expect_error(validate_omics_bundle(p1), "manifest.json does not match checksums.sha256")
    expect_error(read_omics_bundle(p1, parts = "meta"), "does not match checksums.sha256")

    p2 <- write_omics_bundle(example_parts(), example_info(), root,
                             fixed_provenance("2026-02-01T10:00:00Z"))
    rewrite_manifest(p2, function(m) { m$provenance$git_commit <- strrep("f", 40); m },
                     refresh_checksum = FALSE)
    expect_error(validate_omics_bundle(p2), "manifest.json does not match checksums.sha256")
})

test_that("a missing or corrupted checksums.sha256 fails validation", {
    root <- withr::local_tempdir()

    p1 <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())
    sums <- file.path(p1, "checksums.sha256")
    make_writable(sums)
    file.remove(sums)
    expect_error(validate_omics_bundle(p1), "has no checksums.sha256")

    p2 <- write_omics_bundle(example_parts(), example_info(), root,
                             fixed_provenance("2026-02-01T10:00:00Z"))
    sums <- file.path(p2, "checksums.sha256")
    make_writable(sums)
    writeLines("not a checksum line", sums)
    expect_error(validate_omics_bundle(p2), "is not one '<sha256>  manifest.json' line")

    # Well-formed, but the hash of something else.
    writeLines(paste0(strrep("0", 64), "  manifest.json"), sums)
    expect_error(validate_omics_bundle(p2), "manifest.json does not match checksums.sha256")

    # The right hash without the file name breaks `sha256sum -c`, so it is refused.
    writeLines(.bundle_sha256_file(file.path(p2, "manifest.json")), sums)
    expect_error(validate_omics_bundle(p2), "is not one '<sha256>  manifest.json' line")
})

test_that("a missing part file is reported by name", {
    root <- withr::local_tempdir()
    path <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())
    file <- file.path(path, "parts", "expr.rds")
    make_writable(file)
    file.remove(file)

    expect_error(read_omics_bundle(path), "part 'expr' is missing")
})

test_that("an invalid part name in the manifest is refused before any path is built", {
    root <- withr::local_tempdir()
    path <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())
    # Self-consistent: the file matches the name, and the manifest checksum is
    # refreshed, so only the part-name check can stop it.
    rewrite_manifest(path, function(m) {
        m$parts[["../../outside"]] <- list(file = "parts/../../outside.rds",
                                           sha256 = strrep("0", 64), class = "character")
        m
    })

    expect_error(validate_omics_bundle(path), "invalid part name")
    # Checked for every part, even when the caller asks only for a valid one.
    expect_error(read_omics_bundle(path, parts = "meta"), "invalid part name")
})

test_that("a manifest cannot point a reader outside the bundle", {
    root <- withr::local_tempdir()
    path <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())
    rewrite_manifest(path, function(m) {
        m$parts$expr$file <- "../../outside.rds"
        m
    })

    expect_error(read_omics_bundle(path, parts = "expr"), "does not point to parts/expr.rds")
})

test_that("a bundle that is not complete, renamed, or of another contract is refused", {
    root <- withr::local_tempdir()

    p1 <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())
    rewrite_manifest(p1, function(m) { m$status <- "writing"; m })
    expect_error(validate_omics_bundle(p1), "is not complete")
    expect_equal(nrow(list_omics_bundles(root)), 0L)

    p2 <- write_omics_bundle(example_parts(), example_info(), root,
                             fixed_provenance("2026-02-01T10:00:00Z"))
    rewrite_manifest(p2, function(m) { m$bundle_contract_version <- "2.0.0"; m })
    expect_error(validate_omics_bundle(p2), "contract version '2.0.0'")

    p3 <- write_omics_bundle(example_parts(), example_info(), root,
                             fixed_provenance("2026-02-02T10:00:00Z"))
    moved <- file.path(dirname(p3), "renamed")
    expect_true(file.rename(p3, moved))
    expect_error(validate_omics_bundle(moved), "does not match its manifest ID")
    # A renamed bundle cannot be read, so it is not offered either.
    expect_false(moved %in% list_omics_bundles(root)$path)
})

test_that("half-written directories are invisible to readers", {
    root <- withr::local_tempdir()
    omics_dir <- file.path(root, "bundles", "proteomics")
    partial <- file.path(omics_dir, ".partial_20260131T093000Z_123")
    dir.create(file.path(partial, "parts"), recursive = TRUE)
    saveRDS(1, file.path(partial, "parts", "expr.rds"))
    no_manifest <- file.path(omics_dir, "proteomics-demo-r1-20260131T093000Z-deadbeef")
    dir.create(no_manifest)

    expect_equal(nrow(list_omics_bundles(root)), 0L)
    expect_error(validate_omics_bundle(partial), "incomplete bundle left by a failed write")
    expect_error(validate_omics_bundle(no_manifest), "has no manifest.json")
})

test_that("a bundle whose publish was interrupted is not listed or loaded", {
    root <- withr::local_tempdir()

    # Moved parts and manifest, but checksums.sha256 never arrived.
    p1 <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())
    sums <- file.path(p1, "checksums.sha256")
    make_writable(sums)
    file.remove(sums)
    expect_equal(nrow(list_omics_bundles(root)), 0L)
    expect_error(read_omics_bundle(p1), "has no checksums.sha256: it is incomplete")

    # Parts only, and a manifest cut off mid-write: listing must not parse it.
    omics_dir <- file.path(root, "bundles", "proteomics")
    cut <- file.path(omics_dir, "proteomics-demo-r1-20260201T100000Z-0badc0de")
    dir.create(file.path(cut, "parts"), recursive = TRUE)
    writeLines('{"bundle_contract_version": "1.0', file.path(cut, "manifest.json"))
    parts_only <- file.path(omics_dir, "proteomics-demo-r1-20260202T100000Z-0badc0de")
    dir.create(file.path(parts_only, "parts"), recursive = TRUE)

    expect_equal(nrow(list_omics_bundles(root)), 0L)
    expect_error(validate_omics_bundle(cut), "has no checksums.sha256")
    expect_error(validate_omics_bundle(parts_only), "has no manifest.json")
})

test_that("a failed write leaves no bundle behind", {
    root <- withr::local_tempdir()
    not_a_dir <- file.path(root, "file_not_dir")
    writeLines("x", not_a_dir)

    expect_error(write_omics_bundle(example_parts(), example_info(), not_a_dir,
                                    fixed_provenance()),
                 "Cannot create the bundle directory")
    expect_equal(nrow(list_omics_bundles(root)), 0L)
})

test_that("an unreadable manifest is skipped by the listing but still refused on read", {
    root <- withr::local_tempdir()
    good <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())
    bad <- write_omics_bundle(example_parts(), example_info(), root,
                              fixed_provenance("2026-02-01T10:00:00Z"))
    # Malformed JSON with a matching checksum, so only the parse can fail.
    manifest <- file.path(bad, "manifest.json")
    make_writable(manifest)
    writeLines('{"bundle_contract_version": "1.0', manifest)
    sums <- file.path(bad, "checksums.sha256")
    make_writable(sums)
    writeLines(paste0(.bundle_sha256_file(manifest), "  manifest.json"), sums)

    expect_identical(list_omics_bundles(root)$bundle_id, basename(good))
    expect_error(validate_omics_bundle(bad), "not readable JSON")
    expect_error(read_omics_bundle(bad), "not readable JSON")
})

test_that("a bundle moved under another omics folder is not listed", {
    root <- withr::local_tempdir()
    path <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())
    wrong <- file.path(root, "bundles", "rnaseq")
    dir.create(wrong)
    expect_true(file.rename(path, file.path(wrong, basename(path))))

    expect_equal(nrow(list_omics_bundles(root)), 0L)
    expect_equal(nrow(list_omics_bundles(root, omics = "rnaseq")), 0L)
})

test_that("listing filters by omics and is empty for a root with no bundles", {
    root <- withr::local_tempdir()
    empty <- list_omics_bundles(root)
    expect_equal(nrow(empty), 0L)
    expect_named(empty, c("bundle_id", "omics", "project", "analysis_round",
                          "created_at_utc", "path"))

    write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())
    write_omics_bundle(example_parts(), example_info(omics = "rnaseq"), root,
                       fixed_provenance())
    expect_equal(nrow(list_omics_bundles(root)), 2L)
    expect_identical(list_omics_bundles(root, omics = "rnaseq")$omics, "rnaseq")
})

test_that("invalid inputs are refused before anything is written", {
    root <- withr::local_tempdir()
    prov <- fixed_provenance()

    expect_error(write_omics_bundle(list(), example_info(), root, prov), "non-empty named list")
    expect_error(write_omics_bundle(list(1), example_info(), root, prov), "non-empty named list")
    expect_error(write_omics_bundle(list(Bad = 1), example_info(), root, prov), "part names")
    expect_error(write_omics_bundle(list(`../x` = 1), example_info(), root, prov), "part names")
    expect_error(write_omics_bundle(list(a = 1, a = 2), example_info(), root, prov), "unique")
    expect_error(write_omics_bundle(example_parts(), list(project = "P", analysis_round = "R"),
                                    root, prov), "info\\$omics")
    expect_error(write_omics_bundle(example_parts(), example_info(omics = "RNA-seq"), root, prov),
                 "info\\$omics must be lowercase")
    expect_error(write_omics_bundle(example_parts(), example_info(table = data.frame(x = 1)),
                                    root, prov), "Extra info fields")
    expect_error(write_omics_bundle(example_parts(), example_info(), root,
                                    fixed_provenance("2026-01-31 09:30")), "created_at_utc")
    expect_false(dir.exists(file.path(root, "bundles")))
})

test_that("bundle files are made read-only", {
    skip_on_os("windows")
    root <- withr::local_tempdir()
    path <- write_omics_bundle(example_parts(), example_info(), root, fixed_provenance())
    files <- list.files(path, recursive = TRUE, full.names = TRUE)
    expect_true(all(file.mode(files) == as.octmode("0444")))
})

test_that("provenance is captured without failing outside a git repository", {
    dir <- withr::local_tempdir()
    p <- bundle_provenance(dir)

    expect_match(p$created_at_utc, "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$")
    expect_true(is.na(p$git_commit))
    expect_true(is.na(p$git_dirty))
    expect_true(is.na(p$renv_lock_sha256))
    expect_true("base" %in% names(p$packages))
    # What bundle_provenance() returns is accepted by the writer as is.
    path <- write_omics_bundle(example_parts(), example_info(), dir, p)
    # An unknown commit is JSON null in the manifest; read back it is NULL or NA.
    commit <- read_omics_bundle(path, parts = "meta")$manifest$provenance$git_commit
    expect_true(is.null(commit) || isTRUE(is.na(commit)))
})

test_that("git_dirty is FALSE or TRUE when git answers, and NA when git status fails", {
    skip_if(!nzchar(Sys.which("git")), "git is not installed")
    # Everything happens in a throw-away repository; identity and signing are
    # set per command with -c, so no global git configuration is read or changed.
    repo <- withr::local_tempdir()
    git <- function(...) {
        system2("git", c("-C", shQuote(repo), "-c", "user.name=test",
                         "-c", "user.email=test@example.com",
                         "-c", "commit.gpgsign=false", ...),
                stdout = FALSE, stderr = FALSE)
    }
    expect_equal(git("init", "-q"), 0L)
    writeLines("a", file.path(repo, "tracked.txt"))
    expect_equal(git("add", "tracked.txt"), 0L)
    expect_equal(git("commit", "-q", "-m", "init"), 0L)

    expect_false(bundle_provenance(repo)$git_dirty)
    writeLines("b", file.path(repo, "tracked.txt"))
    expect_true(bundle_provenance(repo)$git_dirty)

    # A corrupt index: HEAD still resolves, but git status cannot run, so the
    # state is unknown rather than clean.
    writeLines("not an index", file.path(repo, ".git", "index"))
    p <- bundle_provenance(repo)
    expect_false(is.na(p$git_commit))
    expect_true(is.na(p$git_dirty))
})
