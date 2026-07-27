#!/usr/bin/env bash
#
# Create an isolated Python venv for the TabPFN classifier-benchmark POC and
# freeze the fully-resolved dependency tree to a lockfile for reproducibility.
#
# Usage:
#   bash setup-tabpfn-venv.sh            # install from the top-level pin
#   USE_LOCK=1 bash setup-tabpfn-venv.sh # install the exact locked tree
#
# Commit BOTH requirements-tabpfn.txt (top-level pin) and
# requirements-tabpfn.lock (full resolved tree) to git.
#
# Point the R harness at the resulting interpreter:
#   export TABPFN_PYTHON="$(pwd)/envs/tabpfn/bin/python"
#
# Note: TabPFN downloads its pretrained weights on first use, so the first
# benchmark run needs one-time network access.

set -euo pipefail

VENV_DIR="${VENV_DIR:-envs/tabpfn}"
REQ_FILE="${REQ_FILE:-requirements-tabpfn.txt}"
LOCK_FILE="${LOCK_FILE:-requirements-tabpfn.lock}"
USE_LOCK="${USE_LOCK:-0}"

if [[ ! -f "${REQ_FILE}" ]]; then
  echo "ERROR: ${REQ_FILE} not found (run from the project root)." >&2
  exit 1
fi

echo ">> Creating venv at ${VENV_DIR}"
rm -rf "${VENV_DIR}"          # start clean so the lockfile is an exact reproduction
# --copies: real interpreter copy instead of a symlink. A symlinked bin/python can
# be canonicalised back to the base interpreter by anything that resolves paths,
# which then runs WITHOUT the venv (no tabpfn). A copy is immune.
python3 -m venv --copies "${VENV_DIR}"

# venv Python lives under bin/ on Unix and Scripts/ on Windows (e.g. Git-Bash).
VENV_PY="${VENV_DIR}/bin/python"
[ -x "${VENV_PY}" ] || VENV_PY="${VENV_DIR}/Scripts/python.exe"

echo ">> Upgrading pip"
"${VENV_PY}" -m pip install --upgrade pip

if [[ "${USE_LOCK}" == "1" && -f "${LOCK_FILE}" ]]; then
  echo ">> Installing pinned tree from ${LOCK_FILE}"
  "${VENV_PY}" -m pip install -r "${LOCK_FILE}"
else
  echo ">> Installing ${REQ_FILE}"
  "${VENV_PY}" -m pip install -r "${REQ_FILE}"
  echo ">> Freezing resolved tree to ${LOCK_FILE}"
  "${VENV_PY}" -m pip freeze > "${LOCK_FILE}"
fi

echo ">> Smoke test: tabpfn + sklearn import"
"${VENV_PY}" -c "import tabpfn, sklearn; print('tabpfn', getattr(tabpfn, '__version__', '?'), '| sklearn', sklearn.__version__)"

echo
echo "Done. Point the R harness at this interpreter:"
echo "  export TABPFN_PYTHON=\"\$(pwd)/${VENV_PY}\""
