#!/usr/bin/env bash
#
# Create an isolated Python venv for the mummichog v2 pipeline stage and freeze
# the fully-resolved dependency tree to a lockfile for reproducibility.
#
# Usage:
#   bash setup-mummichog-venv.sh            # install from the top-level pin
#   USE_LOCK=1 bash setup-mummichog-venv.sh # install the exact locked tree
#
# Commit BOTH requirements-mummichog.txt (top-level pin) and
# requirements-mummichog.lock (full resolved tree) to git.
#
# Point the R wrapper at the resulting interpreter:
#   export MUMMICHOG_PYTHON="$(pwd)/envs/mummichog/bin/python"

set -euo pipefail

VENV_DIR="${VENV_DIR:-envs/mummichog}"
REQ_FILE="${REQ_FILE:-requirements-mummichog.txt}"
LOCK_FILE="${LOCK_FILE:-requirements-mummichog.lock}"
USE_LOCK="${USE_LOCK:-0}"

if [[ ! -f "${REQ_FILE}" ]]; then
  echo "ERROR: ${REQ_FILE} not found (run from the project root)." >&2
  exit 1
fi

echo ">> Creating venv at ${VENV_DIR}"
python3 -m venv "${VENV_DIR}"

echo ">> Upgrading pip"
"${VENV_DIR}/bin/python" -m pip install --upgrade pip

if [[ "${USE_LOCK}" == "1" && -f "${LOCK_FILE}" ]]; then
  echo ">> Installing pinned tree from ${LOCK_FILE}"
  "${VENV_DIR}/bin/python" -m pip install -r "${LOCK_FILE}"
else
  echo ">> Installing ${REQ_FILE}"
  "${VENV_DIR}/bin/python" -m pip install -r "${REQ_FILE}"
  echo ">> Freezing resolved tree to ${LOCK_FILE}"
  "${VENV_DIR}/bin/python" -m pip freeze > "${LOCK_FILE}"
fi

echo ">> Smoke test: mummichog imports and reports its version"
"${VENV_DIR}/bin/python" -c "from mummichog.config import VERSION; print('mummichog', VERSION)"

echo
echo "Done. Point the R wrapper at this interpreter:"
echo "  export MUMMICHOG_PYTHON=\"\$(pwd)/${VENV_DIR}/bin/python\""
