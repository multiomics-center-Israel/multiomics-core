#!/usr/bin/env Rscript
# =============================================================================
# install_pandoc.R
#
# Idempotent pandoc installer for multiomics-core. Called by install.bat so a
# fresh machine (Technion proteomics center, student laptop, etc.) can render
# the HTML report without requiring RStudio Desktop.
#
# Behaviour:
#   - If pandoc >= 2.0 is already discoverable (on PATH, RStudio bundled,
#     tools/pandoc/), exit 0 with no download.
#   - Otherwise on Windows, download a portable pandoc zip from GitHub,
#     extract pandoc.exe into tools/pandoc/, and validate.
#   - On macOS / Linux, print instructions (we don't auto-install on those
#     platforms).
# =============================================================================

PANDOC_MIN <- "2.0"
PANDOC_DL  <- "3.1.12.3"  # pinned release; update when bumping

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    install.packages("rmarkdown", repos = "https://cloud.r-project.org")
}

probe_dirs <- function() {
    if (.Platform$OS.type == "windows") {
        pf   <- Sys.getenv("ProgramFiles",      "C:/Program Files")
        pf86 <- Sys.getenv("ProgramFiles(x86)", "C:/Program Files (x86)")
        lad  <- Sys.getenv("LOCALAPPDATA",
                           file.path(Sys.getenv("USERPROFILE"), "AppData/Local"))
        d <- c(
            file.path(pf,   "RStudio", "resources", "app", "bin", "quarto", "bin", "tools"),
            file.path(pf,   "RStudio", "bin", "quarto", "bin", "tools"),
            file.path(pf,   "RStudio", "bin", "pandoc"),
            file.path(pf86, "RStudio", "bin", "pandoc"),
            file.path(pf,   "Pandoc"),
            file.path(pf86, "Pandoc"),
            file.path(lad,  "Pandoc"),
            file.path(getwd(), "tools", "pandoc")
        )
    } else {
        d <- c("/usr/bin", "/usr/local/bin", "/opt/homebrew/bin",
               file.path(getwd(), "tools", "pandoc"))
    }
    d[dir.exists(d)]
}

try_register <- function(dir) {
    tryCatch({
        rmarkdown::find_pandoc(dir = dir)
        isTRUE(rmarkdown::pandoc_available(version = PANDOC_MIN))
    }, error = function(e) FALSE)
}

# 1. Is pandoc already reachable (PATH / RSTUDIO_PANDOC)?
if (isTRUE(rmarkdown::pandoc_available(version = PANDOC_MIN))) {
    message("  pandoc already available (version ",
            rmarkdown::pandoc_version(), ")")
    quit(save = "no", status = 0)
}

# 2. Try each common install location.
exe <- if (.Platform$OS.type == "windows") "pandoc.exe" else "pandoc"
for (d in probe_dirs()) {
    if (file.exists(file.path(d, exe)) && try_register(d)) {
        message("  Found pandoc at: ", d,
                " (version ", rmarkdown::pandoc_version(), ")")
        quit(save = "no", status = 0)
    }
}

# 3. Not found.
if (.Platform$OS.type != "windows") {
    message("  pandoc not found. Install it via your system package manager:")
    message("    Ubuntu/Debian:  sudo apt-get install pandoc")
    message("    macOS:          brew install pandoc")
    quit(save = "no", status = 1)
}

# 4. Windows: auto-download a portable pandoc.
url <- sprintf(
    "https://github.com/jgm/pandoc/releases/download/%s/pandoc-%s-windows-x86_64.zip",
    PANDOC_DL, PANDOC_DL
)
dest_dir <- file.path(getwd(), "tools", "pandoc")
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
zip_path <- file.path(dest_dir, "pandoc.zip")

message("  Downloading pandoc ", PANDOC_DL, " ...")
ok <- tryCatch({
    # mode = "wb" is critical on Windows for binary zip files.
    utils::download.file(url, zip_path, mode = "wb", quiet = FALSE)
    TRUE
}, error = function(e) {
    message("  Download failed: ", e$message)
    FALSE
})

if (!ok) {
    message("  Please install pandoc manually: https://pandoc.org/installing.html")
    quit(save = "no", status = 1)
}

tmp_ex <- file.path(dest_dir, "_tmp")
dir.create(tmp_ex, showWarnings = FALSE)
utils::unzip(zip_path, exdir = tmp_ex)

exe_paths <- list.files(tmp_ex, pattern = "^pandoc\\.exe$",
                         recursive = TRUE, full.names = TRUE)
if (length(exe_paths) == 0) {
    message("  Could not locate pandoc.exe inside the downloaded archive.")
    unlink(zip_path); unlink(tmp_ex, recursive = TRUE)
    quit(save = "no", status = 1)
}
file.copy(exe_paths[1], file.path(dest_dir, "pandoc.exe"), overwrite = TRUE)
unlink(zip_path)
unlink(tmp_ex, recursive = TRUE)

if (try_register(dest_dir)) {
    message("  pandoc ", rmarkdown::pandoc_version(), " installed to: ", dest_dir)
    quit(save = "no", status = 0)
}

message("  Downloaded pandoc but rmarkdown still can't detect it at: ", dest_dir)
quit(save = "no", status = 1)
