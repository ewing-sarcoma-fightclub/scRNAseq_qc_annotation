#!/bin/bash
set -euo pipefail

# Loop DropletQC per sample (post-SoupX)
# Usage:
#   DropletQC_loop.sh <root_dir> [out_root] [tiles] [cores] [metadata_root] [metadata_col]

ROOT_DIR=${1:?"Provide root directory with per-sample folders"}
OUT_ROOT=${2:-""}
TILES=${3:-100}
CORES=${4:-""}
METADATA_ROOT=${5:-""}
METADATA_COL=${6:-""}

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PIPELINE_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
DQC_SCRIPT="${PIPELINE_ROOT}/r/DropletQC_per_sample.R"
R_BIN="${R_BIN:-Rscript}"
RESUME="${RESUME:-true}"

if [[ ! -f "${DQC_SCRIPT}" ]]; then
  echo "DropletQC script not found: ${DQC_SCRIPT}" >&2
  exit 1
fi

while IFS= read -r -d '' sample_dir; do
  sample_name=$(basename "${sample_dir}")

  if [[ -n "${OUT_ROOT}" ]]; then
    out_dir="${OUT_ROOT}/${sample_name}"
  else
    out_dir="${sample_dir}/DropletQC"
  fi
  out_csv="${out_dir}/${sample_name}_dropletqc_results.csv"
  if [[ "${RESUME}" == "true" ]] && [[ -f "${out_csv}" ]]; then
    echo "[INFO] Skipping DropletQC (resume): ${out_csv} exists" >&2
    continue
  fi

  bam_path=""
  for candidate in \
    "${sample_dir}/possorted_genome_bam.bam" \
    "${sample_dir}/outs/possorted_genome_bam.bam"; do
    if [[ -f "${candidate}" ]]; then
      bam_path="${candidate}"
      break
    fi
  done

  if [[ -z "${bam_path}" ]]; then
    echo "[WARN] DropletQC skipped for ${sample_name}: possorted_genome_bam.bam not found (Cell Ranger may have been run with --no-bam)." >&2
    continue
  fi

  if [[ ! -f "${bam_path}.bai" && ! -f "${bam_path%.bam}.bai" ]]; then
    if command -v samtools >/dev/null 2>&1; then
      echo "[INFO] Indexing BAM for ${sample_name}" >&2
      if ! samtools index "${bam_path}"; then
        echo "[WARN] DropletQC skipped for ${sample_name}: failed to index BAM with samtools." >&2
        continue
      fi
    else
      echo "[WARN] DropletQC skipped for ${sample_name}: BAM index missing and samtools not found." >&2
      continue
    fi
  fi

  metadata_csv=""
  if [[ -n "${METADATA_ROOT}" ]]; then
    metadata_csv="${METADATA_ROOT}/${sample_name}/metadata.csv"
  fi

  if [[ -n "${CORES}" ]]; then
    "${R_BIN}" "${DQC_SCRIPT}" "${sample_dir}" "${sample_name}" "${out_dir}" "${TILES}" "${CORES}" "${metadata_csv}" "${METADATA_COL}"
  else
    "${R_BIN}" "${DQC_SCRIPT}" "${sample_dir}" "${sample_name}" "${out_dir}" "${TILES}" "" "${metadata_csv}" "${METADATA_COL}"
  fi

done < <(find "${ROOT_DIR}" -mindepth 1 -maxdepth 1 -type d -print0)

echo "Done."
