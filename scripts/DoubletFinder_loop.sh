#!/bin/bash
set -euo pipefail

# Loop DoubletFinder per sample (post-SoupX)
# Usage:
#   DoubletFinder_loop.sh <root_dir> [doublet_rate] [npcs] [out_root] [soupx_root] [metadata_root] [min_features] [min_cells] [max_percent_mt] [seed]
#
# root_dir: directory containing one subfolder per sample
# doublet_rate: expected doublet rate (default: 0.075)
# npcs: number of PCs (default: 30)
# out_root: optional base output directory (default: per-sample DoubletFinder/)
# soupx_root: optional SoupX output root (default: per-sample SoupX/)
# metadata_root: optional metadata output root (default: per-sample metadata.csv)
# min_features/min_cells/max_percent_mt: basic Seurat QC filters prior to DoubletFinder

ROOT_DIR=${1:?"Provide root directory with per-sample folders"}
DOUBLEt_RATE=${2:-0.075}
NPCS=${3:-30}
OUT_ROOT=${4:-""}
SOUPX_ROOT=${5:-""}
METADATA_ROOT=${6:-""}
MIN_FEATURES=${7:-200}
MIN_CELLS=${8:-3}
MAX_PERCENT_MT=${9:-""}
SEED=${10:-1}

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PIPELINE_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
DF_SCRIPT="${PIPELINE_ROOT}/r/DoubletFinder_after_SoupX.R"
R_BIN="${R_BIN:-Rscript}"
RESUME="${RESUME:-true}"

if [[ ! -f "${DF_SCRIPT}" ]]; then
  echo "DoubletFinder script not found: ${DF_SCRIPT}" >&2
  exit 1
fi

while IFS= read -r -d '' sample_dir; do
  sample_name=$(basename "${sample_dir}")

  if [[ -n "${SOUPX_ROOT}" ]]; then
    soupx_mtx="${SOUPX_ROOT}/${sample_name}/${sample_name}_soupx.mtx"
  else
    soupx_mtx=""
  fi

  if [[ -n "${METADATA_ROOT}" ]]; then
    metadata_csv="${METADATA_ROOT}/${sample_name}/metadata.csv"
  else
    metadata_csv=""
  fi

  if [[ -n "${OUT_ROOT}" ]]; then
    out_dir="${OUT_ROOT}/${sample_name}"
    out_csv="${out_dir}/${sample_name}_metadata_with_DoubletFinder.csv"
    if [[ "${RESUME}" == "true" ]] && [[ -f "${out_csv}" ]]; then
      echo "[INFO] Skipping DoubletFinder (resume): ${out_csv} exists" >&2
      continue
    fi
    "${R_BIN}" "${DF_SCRIPT}" "${sample_dir}" "${sample_name}" "${soupx_mtx}" "${out_dir}" \
      "${DOUBLEt_RATE}" "${NPCS}" "${metadata_csv}" "${MIN_FEATURES}" "${MIN_CELLS}" "${MAX_PERCENT_MT}" "${SEED}"
  else
    out_dir="${sample_dir}/DoubletFinder"
    out_csv="${out_dir}/${sample_name}_metadata_with_DoubletFinder.csv"
    if [[ "${RESUME}" == "true" ]] && [[ -f "${out_csv}" ]]; then
      echo "[INFO] Skipping DoubletFinder (resume): ${out_csv} exists" >&2
      continue
    fi
    "${R_BIN}" "${DF_SCRIPT}" "${sample_dir}" "${sample_name}" "${soupx_mtx}" "" \
      "${DOUBLEt_RATE}" "${NPCS}" "${metadata_csv}" "${MIN_FEATURES}" "${MIN_CELLS}" "${MAX_PERCENT_MT}" "${SEED}"
  fi

done < <(find "${ROOT_DIR}" -mindepth 1 -maxdepth 1 -type d -print0)

echo "Done."
