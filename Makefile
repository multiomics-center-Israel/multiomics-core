# Makefile — project helper targets.
#
# Only the mummichog venv helpers live here for now; the R pipeline is driven by
# {targets} (targets::tar_make()), not make.

SHELL := /bin/bash

.PHONY: help setup mummichog-venv mummichog-lock

help:
	@echo "Targets:"
	@echo "  setup           One-time per machine: build the mummichog venv AND record"
	@echo "                  MUMMICHOG_PYTHON in .Renviron so tar_make() just works."
	@echo "  mummichog-venv  Create the isolated mummichog v2 venv (envs/mummichog)"
	@echo "                  and freeze requirements-mummichog.lock."
	@echo "  mummichog-lock  Recreate the venv from the committed lockfile exactly."

# One-time machine setup: build the venv, then pin its interpreter in .Renviron
# (idempotent — updates the line if present) so R picks up MUMMICHOG_PYTHON on
# start and the pipeline needs no manual `export`. .Renviron is gitignored.
setup: mummichog-venv
	@Rscript -e 'source("R/core/18_renviron.R"); cands <- c("envs/mummichog/bin/python", "envs/mummichog/Scripts/python.exe"); py <- cands[file.exists(cands)][1]; if (is.na(py)) stop("venv python not found under envs/mummichog after build", call. = FALSE); py <- normalizePath(py); act <- set_renviron_var("MUMMICHOG_PYTHON", py); message(sprintf("[setup] %s MUMMICHOG_PYTHON=%s in .Renviron", act, py))'
	@echo "[setup] Done. Start a fresh R session so it reads .Renviron; then targets::tar_make() picks up MUMMICHOG_PYTHON automatically."

# Build the pinned venv from requirements-mummichog.txt and refresh the lockfile.
mummichog-venv:
	bash setup-mummichog-venv.sh
	@echo "Now run 'make setup' to record MUMMICHOG_PYTHON in .Renviron (or export it yourself)."

# Reproduce the exact locked environment (no re-freeze).
mummichog-lock:
	USE_LOCK=1 bash setup-mummichog-venv.sh
