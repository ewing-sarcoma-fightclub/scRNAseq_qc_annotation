#!/bin/bash
set -euo pipefail

# QC pipeline after Cell Ranger per-sample outs
# Steps:
# 0) AmbiQuant (pre-SoupX)
# 1) EmptyDrops (annotation only)
# 2) Seurat clustering metadata for SoupX (no filtering)
# 3) SoupX per sample
# 4) DoubletFinder (post-SoupX)
# 5) DropletQC (post-SoupX)
# 6) AmbiQuant (post-SoupX)
# 7) Final subsetting + normalization/clustering in Seurat
# 8) Merge (non-integrated) + AddModuleScore + optional integration

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PIPELINE_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
if [[ -f "${PIPELINE_ROOT}/env/config.local.env" ]]; then
  CONFIG_FILE="${PIPELINE_ROOT}/env/config.local.env"
else
  CONFIG_FILE="${PIPELINE_ROOT}/env/config.env"
fi
ROOT_DIR=""
OUT_QC=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_FILE="$2"
      shift 2
      ;;
    --root)
      ROOT_DIR="$2"
      shift 2
      ;;
    --out)
      OUT_QC="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "${ROOT_DIR}" || -z "${OUT_QC}" ]]; then
  echo "Usage: pipeline_QC_after_cellranger.sh --root ROOT_DIR --out OUT_QC [--config CONFIG]" >&2
  exit 1
fi

if [[ -f "${CONFIG_FILE}" ]]; then
  CONFIG_FILE="$(cd -- "$(dirname -- "${CONFIG_FILE}")" && pwd)/$(basename -- "${CONFIG_FILE}")"
fi
if [[ -f "${CONFIG_FILE}" ]]; then
  # shellcheck source=/dev/null
  source "${CONFIG_FILE}"
fi

# ---- configure ----
EMPTYDROPS_LOOP="${PIPELINE_ROOT}/scripts/EmptyDrops_loop.sh"
AMBIQUANT_SCRIPT="${PIPELINE_ROOT}/python/AmbiQuant_pre_post.py"
SOUPX_SCRIPT="${PIPELINE_ROOT}/r/SoupX_per_sample.R"
DOUBLETFINDER_LOOP="${PIPELINE_ROOT}/scripts/DoubletFinder_loop.sh"
DROPLETQC_LOOP="${PIPELINE_ROOT}/scripts/DropletQC_loop.sh"
SEURAT_FINAL_SCRIPT="${PIPELINE_ROOT}/r/Seurat_final_subsetting.R"
SEURAT_ANNOTATE_SCRIPT="${PIPELINE_ROOT}/r/Seurat_merge_annotate_integrate.R"

MICROMAMBA_BIN="${MICROMAMBA_BIN:-micromamba}"
AMBIQUANT_ENV="${AMBIQUANT_ENV:-}"
AMBIQUANT_REPO="${AMBIQUANT_REPO:-}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
R_BIN="${R_BIN:-Rscript}"
R_LIBS_USER="${R_LIBS_USER:-NULL}"

export PYTHON_BIN
export R_BIN
export R_LIBS_USER
unset R_LIBS || true

EMPTYDROPS_OUT="${OUT_QC}/emptydrops"
SEURAT_META_OUT="${OUT_QC}/seurat_metadata"
SOUPX_OUT="${OUT_QC}/soupx"
DOUBLETFINDER_OUT="${OUT_QC}/doubletfinder"
DROPLETQC_OUT="${OUT_QC}/dropletqc"
DOUBLETFINDER_RATE="${DOUBLETFINDER_RATE:-0.075}"
DOUBLETFINDER_NPCS="${DOUBLETFINDER_NPCS:-30}"
DOUBLETFINDER_MIN_FEATURES="${DOUBLETFINDER_MIN_FEATURES:-200}"
DOUBLETFINDER_MIN_CELLS="${DOUBLETFINDER_MIN_CELLS:-3}"
DOUBLETFINDER_MAX_PERCENT_MT="${DOUBLETFINDER_MAX_PERCENT_MT:-}"
SEURAT_NPCS="${SEURAT_NPCS:-30}"
SEURAT_RESOLUTION="${SEURAT_RESOLUTION:-0.5}"
SEURAT_NFEATURES="${SEURAT_NFEATURES:-2000}"
SEURAT_SEED="${SEURAT_SEED:-1}"
SEURAT_FINAL_OUT="${OUT_QC}/seurat_qc"
SEURAT_FINAL_USE_SOUPX="${SEURAT_FINAL_USE_SOUPX:-true}"
SEURAT_FINAL_REQUIRE_EMPTYDROPS="${SEURAT_FINAL_REQUIRE_EMPTYDROPS:-true}"
SEURAT_FINAL_REQUIRE_DROPLETQC="${SEURAT_FINAL_REQUIRE_DROPLETQC:-true}"
SEURAT_FINAL_NPCS="${SEURAT_FINAL_NPCS:-30}"
SEURAT_FINAL_RESOLUTION="${SEURAT_FINAL_RESOLUTION:-0.5}"
SEURAT_FINAL_SEED="${SEURAT_FINAL_SEED:-1}"
SEURAT_FINAL_MIN_CELLS="${SEURAT_FINAL_MIN_CELLS:-3}"
SEURAT_FINAL_MIN_FEATURES="${SEURAT_FINAL_MIN_FEATURES:-200}"
SEURAT_FINAL_REGRESS="${SEURAT_FINAL_REGRESS:-false}"
SEURAT_FINAL_REGRESS_CELL_CYCLE="${SEURAT_FINAL_REGRESS_CELL_CYCLE:-false}"
SEURAT_FINAL_MAX_PERCENT_MT="${SEURAT_FINAL_MAX_PERCENT_MT:-${DOUBLETFINDER_MAX_PERCENT_MT:-}}"
SEURAT_FINAL_MAX_PERCENT_RIBO="${SEURAT_FINAL_MAX_PERCENT_RIBO:-}"
SEURAT_FINAL_BASIC_MIN_FEATURES="${SEURAT_FINAL_BASIC_MIN_FEATURES:-${DOUBLETFINDER_MIN_FEATURES:-}}"
SEURAT_ANNOTATE_OUT="${OUT_QC}/../annotation_integration/seurat_annotation_integration"
DOUBLETFINDER_SEED="${DOUBLETFINDER_SEED:-1}"
CELL_MARKER_XLSX="${CELL_MARKER_XLSX:-${PIPELINE_ROOT}/resources/Cell_marker_Seq_human.xlsx}"
CELL_MARKER_TISSUES="${CELL_MARKER_TISSUES:-}"
CELL_MARKER_MIN_GENES="${CELL_MARKER_MIN_GENES:-5}"
CELL_MARKER_CTRL="${CELL_MARKER_CTRL:-20}"
CELL_MARKER_SEED="${CELL_MARKER_SEED:-1}"
SEURAT_DO_INTEGRATION="${SEURAT_DO_INTEGRATION:-true}"
SEURAT_INTEGRATION_NPCS="${SEURAT_INTEGRATION_NPCS:-30}"
SEURAT_INTEGRATION_RESOLUTION="${SEURAT_INTEGRATION_RESOLUTION:-0.5}"
SEURAT_INTEGRATION_KWEIGHT="${SEURAT_INTEGRATION_KWEIGHT:-80}"
SEURAT_INTEGRATION_KANCHOR="${SEURAT_INTEGRATION_KANCHOR:-5}"
EWING_SIGNATURE_CSV="${EWING_SIGNATURE_CSV:-${PIPELINE_ROOT}/resources/Aynaud.csv}"
RESUME="${RESUME:-true}"
CELLTYPE_GROUPING_FILE="${CELLTYPE_GROUPING_FILE:-${PIPELINE_ROOT}/resources/celltype_grouping.txt}"
CELLTYPE_GROUPING_PERCENTILE="${CELLTYPE_GROUPING_PERCENTILE:-0.85}"
CELLTYPE_GROUPING_AZ_LINEAGE_MIN="${CELLTYPE_GROUPING_AZ_LINEAGE_MIN:-0.60}"
CELLTYPE_GROUPING_AZ_LINEAGE_DELTA="${CELLTYPE_GROUPING_AZ_LINEAGE_DELTA:-0.07}"

export PIPELINE_ROOT
export CELLTYPE_GROUPING_FILE
export CELLTYPE_GROUPING_PERCENTILE
export CELLTYPE_GROUPING_AZ_LINEAGE_MIN
export CELLTYPE_GROUPING_AZ_LINEAGE_DELTA
export SEURAT_INTEGRATION_KWEIGHT
export SEURAT_INTEGRATION_KANCHOR

mkdir -p "${OUT_QC}" \
  "${EMPTYDROPS_OUT}" \
  "${SEURAT_META_OUT}" \
  "${SOUPX_OUT}" \
  "${DOUBLETFINDER_OUT}" \
  "${DROPLETQC_OUT}" \
  "${SEURAT_FINAL_OUT}" \
  "${SEURAT_ANNOTATE_OUT}"

# ---- run ----

SKIP_AMBIQUANT="false"
if [[ -z "${AMBIQUANT_REPO}" ]]; then
  echo "[WARN] AMBIQUANT_REPO not set; skipping AmbiQuant steps" >&2
  SKIP_AMBIQUANT="true"
fi

# 0) AmbiQuant (pre-SoupX)
if [[ "${SKIP_AMBIQUANT}" != "true" ]]; then
  if [[ "${RESUME}" == "true" ]] && [[ -d "${OUT_QC}/ambiquant_pre" ]] && [[ -n "$(find "${OUT_QC}/ambiquant_pre" -type f -print -quit)" ]]; then
    echo "[INFO] Skipping AmbiQuant pre (resume): ${OUT_QC}/ambiquant_pre already exists" >&2
  else
  if [[ -n "${AMBIQUANT_ENV}" ]]; then
    "${MICROMAMBA_BIN}" run -p "${AMBIQUANT_ENV}" "${PYTHON_BIN}" "${AMBIQUANT_SCRIPT}" \
      --root "${ROOT_DIR}" --out "${OUT_QC}/ambiquant_pre" --skip-post --ambiquant-root "${AMBIQUANT_REPO}"
  else
    "${PYTHON_BIN}" "${AMBIQUANT_SCRIPT}" \
      --root "${ROOT_DIR}" --out "${OUT_QC}/ambiquant_pre" --skip-post --ambiquant-root "${AMBIQUANT_REPO}"
  fi
  fi
fi

# 1) EmptyDrops (outputs under OUT_QC)
bash "${EMPTYDROPS_LOOP}" "${ROOT_DIR}" "${EMPTYDROPS_OUT}"

# 2) Seurat clustering metadata for SoupX (no filtering) + SoupX per sample
while IFS= read -r -d '' sample_dir; do
  sample_name=$(basename "${sample_dir}")
  metadata_csv="${SEURAT_META_OUT}/${sample_name}/metadata.csv"
  out_dir="${SOUPX_OUT}/${sample_name}"
  soupx_mtx="${out_dir}/${sample_name}_soupx.mtx"
  if [[ "${RESUME}" == "true" ]] && [[ -f "${soupx_mtx}" ]]; then
    echo "[INFO] Skipping SoupX (resume): ${soupx_mtx} exists" >&2
    continue
  fi
  "${R_BIN}" "${SOUPX_SCRIPT}" "${sample_dir}" "${sample_name}" "${metadata_csv}" "${out_dir}" \
    "${SEURAT_NPCS}" "${SEURAT_RESOLUTION}" "${SEURAT_NFEATURES}" "${SEURAT_SEED}"
done < <(find "${ROOT_DIR}" -mindepth 1 -maxdepth 1 -type d -print0)

# 4) DoubletFinder per sample
bash "${DOUBLETFINDER_LOOP}" "${ROOT_DIR}" "${DOUBLETFINDER_RATE}" "${DOUBLETFINDER_NPCS}" \
  "${DOUBLETFINDER_OUT}" "${SOUPX_OUT}" "${SEURAT_META_OUT}" \
  "${DOUBLETFINDER_MIN_FEATURES}" "${DOUBLETFINDER_MIN_CELLS}" "${DOUBLETFINDER_MAX_PERCENT_MT}" "${DOUBLETFINDER_SEED}"

# 5) DropletQC per sample (use SoupX metadata Leiden clusters if available)
bash "${DROPLETQC_LOOP}" "${ROOT_DIR}" "${DROPLETQC_OUT}" "100" "" "${SEURAT_META_OUT}"

# 6) AmbiQuant (post-SoupX)
if [[ "${SKIP_AMBIQUANT}" != "true" ]]; then
  if [[ "${RESUME}" == "true" ]] && [[ -d "${OUT_QC}/ambiquant_post" ]] && [[ -n "$(find "${OUT_QC}/ambiquant_post" -type f -print -quit)" ]]; then
    echo "[INFO] Skipping AmbiQuant post (resume): ${OUT_QC}/ambiquant_post already exists" >&2
  else
  if [[ -n "${AMBIQUANT_ENV}" ]]; then
    "${MICROMAMBA_BIN}" run -p "${AMBIQUANT_ENV}" "${PYTHON_BIN}" "${AMBIQUANT_SCRIPT}" \
      --root "${ROOT_DIR}" --out "${OUT_QC}/ambiquant_post" --skip-pre --soupx-root "${SOUPX_OUT}" \
      --ambiquant-root "${AMBIQUANT_REPO}"
  else
    "${PYTHON_BIN}" "${AMBIQUANT_SCRIPT}" \
      --root "${ROOT_DIR}" --out "${OUT_QC}/ambiquant_post" --skip-pre --soupx-root "${SOUPX_OUT}" \
      --ambiquant-root "${AMBIQUANT_REPO}"
  fi
  fi
fi

# 7) Final subsetting + normalization/clustering in Seurat
"${R_BIN}" "${SEURAT_FINAL_SCRIPT}" \
  "${ROOT_DIR}" \
  "${OUT_QC}" \
  "${SEURAT_FINAL_OUT}" \
  "${SEURAT_FINAL_USE_SOUPX}" \
  "${SEURAT_FINAL_REQUIRE_EMPTYDROPS}" \
  "${SEURAT_FINAL_REQUIRE_DROPLETQC}" \
  "${SEURAT_FINAL_NPCS}" \
  "${SEURAT_FINAL_RESOLUTION}" \
  "${SEURAT_FINAL_SEED}" \
  "${SEURAT_FINAL_MIN_CELLS}" \
  "${SEURAT_FINAL_MIN_FEATURES}" \
  "${SEURAT_FINAL_REGRESS}" \
  "${SEURAT_FINAL_REGRESS_CELL_CYCLE}" \
  "${SEURAT_FINAL_MAX_PERCENT_MT}" \
  "${SEURAT_FINAL_MAX_PERCENT_RIBO}" \
  "${SEURAT_FINAL_BASIC_MIN_FEATURES}"

# 8) Merge (non-integrated) + AddModuleScore + optional integration
if [[ -f "${CELL_MARKER_XLSX}" && -f "${EWING_SIGNATURE_CSV}" ]]; then
  merged_rds="${SEURAT_ANNOTATE_OUT}/merged_annotated.rds"
  integrated_rds="${SEURAT_ANNOTATE_OUT}/integrated_annotated.rds"
  if [[ "${RESUME}" == "true" ]] && [[ -f "${merged_rds}" ]] && { [[ "${SEURAT_DO_INTEGRATION}" != "true" ]] || [[ -f "${integrated_rds}" ]]; }; then
    echo "[INFO] Skipping annotation/integration (resume): outputs exist in ${SEURAT_ANNOTATE_OUT}" >&2
  else
  "${R_BIN}" "${SEURAT_ANNOTATE_SCRIPT}" \
    "${SEURAT_FINAL_OUT}" \
    "${CELL_MARKER_XLSX}" \
    "${EWING_SIGNATURE_CSV}" \
    "${SEURAT_ANNOTATE_OUT}" \
    "${CELL_MARKER_TISSUES}" \
    "${CELL_MARKER_MIN_GENES}" \
    "${CELL_MARKER_CTRL}" \
    "${CELL_MARKER_SEED}" \
    "${SEURAT_DO_INTEGRATION}" \
    "${SEURAT_INTEGRATION_NPCS}" \
    "${SEURAT_INTEGRATION_RESOLUTION}"
  fi
else
  echo "[WARN] CellMarker xlsx or Ewing signature not found; skipping annotation/integration" >&2
fi

echo "QC pipeline complete."
