# Makefile — project helper targets.
#
# Only the mummichog venv helpers live here for now; the R pipeline is driven by
# {targets} (targets::tar_make()), not make.

SHELL := /bin/bash

.PHONY: help mummichog-venv mummichog-lock

help:
	@echo "Targets:"
	@echo "  mummichog-venv  Create the isolated mummichog v2 venv (envs/mummichog)"
	@echo "                  and freeze requirements-mummichog.lock."
	@echo "  mummichog-lock  Recreate the venv from the committed lockfile exactly."

# Build the pinned venv from requirements-mummichog.txt and refresh the lockfile.
mummichog-venv:
	bash setup-mummichog-venv.sh
	@echo "Now: export MUMMICHOG_PYTHON=\"$$(pwd)/envs/mummichog/bin/python\""

# Reproduce the exact locked environment (no re-freeze).
mummichog-lock:
	USE_LOCK=1 bash setup-mummichog-venv.sh
