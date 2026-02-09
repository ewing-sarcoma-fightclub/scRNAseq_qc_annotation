#!/bin/bash
set -euo pipefail

# Run full pipeline: Cell Ranger -> QC pipeline
# Usage:
#   run_all.sh [--config /path/to/env/config.env] [--cellranger-root PATH] [--qc-out PATH] [--skip-cellranger]

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PIPELINE_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
if [[ -f "${PIPELINE_ROOT}/env/config.local.env" ]]; then
  CONFIG_FILE="${PIPELINE_ROOT}/env/config.local.env"
else
  CONFIG_FILE="${PIPELINE_ROOT}/env/config.env"
fi
CELLRANGER_ROOT=""
QC_OUT=""
SKIP_CELLRANGER="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_FILE="$2"
      shift 2
      ;;
    --cellranger-root)
      CELLRANGER_ROOT="$2"
      shift 2
      ;;
    --qc-out)
      QC_OUT="$2"
      shift 2
      ;;
    --skip-cellranger)
      SKIP_CELLRANGER="true"
      shift 1
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -f "${CONFIG_FILE}" ]]; then
  CONFIG_FILE="$(cd -- "$(dirname -- "${CONFIG_FILE}")" && pwd)/$(basename -- "${CONFIG_FILE}")"
fi

if [[ -f "${CONFIG_FILE}" ]]; then
  # shellcheck source=/dev/null
  source "${CONFIG_FILE}"
else
  echo "[WARN] Config file not found: ${CONFIG_FILE} (using defaults)" >&2
fi

# Defaults (only if not set in config or CLI)
CELLRANGER_ROOT="${CELLRANGER_ROOT:-${CELLRANGER_OUT_ROOT:-"${PIPELINE_ROOT}/outputs/cellranger"}}"
QC_OUT="${QC_OUT:-${QC_OUT_ROOT:-"${PIPELINE_ROOT}/outputs/qc"}}"

mkdir -p "${QC_OUT}"

# Run Cell Ranger (unless skipped)
if [[ "${SKIP_CELLRANGER}" != "true" ]]; then
  mkdir -p "${CELLRANGER_ROOT}"
  bash "${PIPELINE_ROOT}/bin/CellRanger" --config "${CONFIG_FILE}" --out-root "${CELLRANGER_ROOT}"
else
  if [[ ! -d "${CELLRANGER_ROOT}" ]]; then
    echo "[ERROR] --skip-cellranger set but CELLRANGER_ROOT does not exist: ${CELLRANGER_ROOT}" >&2
    echo "Provide --cellranger-root to an existing Cell Ranger output directory." >&2
    exit 1
  fi
fi

# Run QC pipeline
bash "${PIPELINE_ROOT}/bin/pipeline_QC_after_cellranger.sh" --config "${CONFIG_FILE}" --root "${CELLRANGER_ROOT}" --out "${QC_OUT}"

echo "All done."
