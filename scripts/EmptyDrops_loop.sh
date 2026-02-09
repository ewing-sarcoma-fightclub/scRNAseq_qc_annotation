#!/bin/bash
set -euo pipefail

# Loop EmptyDrops per sample (annotation-only)
# Usage:
#   EmptyDrops_loop.sh <root_dir> [out_root]
#
# root_dir: directory containing one subfolder per sample
# out_root: optional base output directory (default: per-sample EmptyDrops/)

ROOT_DIR=${1:?"Provide root directory with per-sample folders"}
OUT_ROOT=${2:-""}

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PIPELINE_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
ED_SCRIPT="${PIPELINE_ROOT}/r/EmptyDrops_per_sample.R"
R_BIN="${R_BIN:-Rscript}"
RESUME="${RESUME:-true}"

if [[ ! -f "${ED_SCRIPT}" ]]; then
  echo "EmptyDrops script not found: ${ED_SCRIPT}" >&2
  exit 1
fi

while IFS= read -r -d '' sample_dir; do
  sample_name=$(basename "${sample_dir}")

  if [[ -n "${OUT_ROOT}" ]]; then
    out_dir="${OUT_ROOT}/${sample_name}"
    out_csv="${out_dir}/${sample_name}_emptydrops_results.csv"
    if [[ "${RESUME}" == "true" ]] && [[ -f "${out_csv}" ]]; then
      echo "[INFO] Skipping EmptyDrops (resume): ${out_csv} exists" >&2
      continue
    fi
    "${R_BIN}" "${ED_SCRIPT}" "${sample_dir}" "${sample_name}" "${out_dir}"
  else
    out_dir="${sample_dir}/EmptyDrops"
    out_csv="${out_dir}/${sample_name}_emptydrops_results.csv"
    if [[ "${RESUME}" == "true" ]] && [[ -f "${out_csv}" ]]; then
      echo "[INFO] Skipping EmptyDrops (resume): ${out_csv} exists" >&2
      continue
    fi
    "${R_BIN}" "${ED_SCRIPT}" "${sample_dir}" "${sample_name}"
  fi

done < <(find "${ROOT_DIR}" -mindepth 1 -maxdepth 1 -type d -print0)

echo "Done."
