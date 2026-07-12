# Makefile — project helper targets.
#
# Only the mummichog venv helpers live here for now; the R pipeline is driven by
# {targets} (targets::tar_make()), not make.

SHELL := /bin/bash

.PHONY: help setup mummichog-venv mummichog-lock

help:
	@echo "Targets:"
	@echo "  setup           One-time per machine: build the mummichog venv and print"
	@echo "                  the MUMMICHOG_PYTHON line to add to your .Renviron."
	@echo "  mummichog-venv  Create the isolated mummichog v2 venv (envs/mummichog)"
	@echo "                  and freeze requirements-mummichog.lock."
	@echo "  mummichog-lock  Recreate the venv from the committed lockfile exactly."

# One-time machine setup: build the venv from the COMMITTED lockfile (via
# mummichog-lock, so the install is reproducible and never re-freezes the
# tracked requirements-mummichog.lock), then print the exact MUMMICHOG_PYTHON
# line for this checkout. We do NOT write .Renviron (a protected, machine-
# specific, gitignored file) — the user adds the line themselves, e.g. by
# copying .Renviron.example.
setup: mummichog-lock
	@venv="$${VENV_DIR:-envs/mummichog}"; \
	case "$$venv" in /*) ;; *) venv="$$(pwd)/$$venv";; esac; \
	py="$$venv/bin/python"; \
	[ -x "$$py" ] || py="$$venv/Scripts/python.exe"; \
	if command -v cygpath >/dev/null 2>&1; then py="$$(cygpath -ma "$$py")"; fi; \
	echo ""; \
	echo "[setup] Venv ready. Add this line to your .Renviron (create the file if you"; \
	echo "        don't have one — it lives in the project root and is gitignored):"; \
	echo ""; \
	echo "    MUMMICHOG_PYTHON=$$py"; \
	echo ""; \
	echo "[setup] (.Renviron.example shows the format. Add just the line above — don't"; \
	echo "        copy the whole template unless you also set its MULTIOMICS_CONFIG.)"; \
	echo "[setup] Then start a fresh R session; targets::tar_make() picks it up."

# Build the pinned venv from requirements-mummichog.txt and refresh the lockfile.
mummichog-venv:
	bash setup-mummichog-venv.sh
	@echo "Venv ready. Run 'make setup' to print the MUMMICHOG_PYTHON line for your .Renviron."

# Reproduce the exact locked environment (no re-freeze).
mummichog-lock:
	USE_LOCK=1 bash setup-mummichog-venv.sh
