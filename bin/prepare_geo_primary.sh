#!/bin/bash
set -euo pipefail

# Prepare GEO GSE277083 primary-only Cell Ranger matrices for this pipeline.
# - Downloads GEO matrices + primary cell metadata
# - Splits combined matrices into per-patient folders
# - Writes into pipeline outputs/cellranger (compatible with --skip-cellranger)

# Defaults
GEO="GSE277083"
PIPELINE_ROOT=""
INPUTS_ROOT=""
OUT_ROOT=""
CLEAN="false"

usage() {
  cat <<EOF
Usage: $0 [--geo GSEXXXXXX] [--inputs-root PATH] [--out-root PATH] [--clean]

Defaults:
  --geo          ${GEO}
  --inputs-root  <pipeline_root>/inputs
  --out-root     <pipeline_root>/outputs/cellranger
  --clean        remove existing output root before writing
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --geo)
      GEO="$2"; shift 2 ;;
    --inputs-root)
      INPUTS_ROOT="$2"; shift 2 ;;
    --out-root)
      OUT_ROOT="$2"; shift 2 ;;
    --clean)
      CLEAN="true"; shift 1 ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
 done

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PIPELINE_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)
INPUTS_ROOT=${INPUTS_ROOT:-"${PIPELINE_ROOT}/inputs"}
OUT_ROOT=${OUT_ROOT:-"${PIPELINE_ROOT}/outputs/cellranger"}

mkdir -p "${INPUTS_ROOT}"

SERIES_PREFIX="${GEO:0:6}"  # e.g. GSE277
SERIES_DIR="${SERIES_PREFIX}nnn"  # e.g. GSE277nnn
FTP_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/${SERIES_DIR}/${GEO}/suppl"

DL_ROOT="${INPUTS_ROOT}/${GEO}"
META_ROOT="${INPUTS_ROOT}/${GEO}_metadata"
SPLIT_ROOT="${INPUTS_ROOT}/${GEO}_primary_split"

mkdir -p "${DL_ROOT}/filtered_feature_bc_matrix" "${DL_ROOT}/raw_feature_bc_matrix" "${META_ROOT}" "${SPLIT_ROOT}"

# Download helper
fetch() {
  local url="$1"
  local out="$2"
  if [[ -f "${out}" ]]; then
    return 0
  fi
  echo "Downloading: ${url}"
  curl -L "${url}" -o "${out}"
}

# Encoded filenames
fetch "${FTP_BASE}/${GEO}%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%5Fbarcodes%2Etsv%2Egz" "${DL_ROOT}/filtered_feature_bc_matrix/barcodes.tsv.gz"
fetch "${FTP_BASE}/${GEO}%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%5Ffeatures%2Etsv%2Egz" "${DL_ROOT}/filtered_feature_bc_matrix/features.tsv.gz"
fetch "${FTP_BASE}/${GEO}%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%5Fmatrix%2Emtx%2Egz" "${DL_ROOT}/filtered_feature_bc_matrix/matrix.mtx.gz"

fetch "${FTP_BASE}/${GEO}%5Fraw%5Ffeature%5Fbc%5Fmatrix%5Fbarcodes%2Etsv%2Egz" "${DL_ROOT}/raw_feature_bc_matrix/barcodes.tsv.gz"
fetch "${FTP_BASE}/${GEO}%5Fraw%5Ffeature%5Fbc%5Fmatrix%5Ffeatures%2Etsv%2Egz" "${DL_ROOT}/raw_feature_bc_matrix/features.tsv.gz"
fetch "${FTP_BASE}/${GEO}%5Fraw%5Ffeature%5Fbc%5Fmatrix%5Fmatrix%2Emtx%2Egz" "${DL_ROOT}/raw_feature_bc_matrix/matrix.mtx.gz"

fetch "${FTP_BASE}/${GEO}%5Fcell%5Fmetadata%5Fprimary%2Etxt%2Egz" "${META_ROOT}/${GEO}_cell_metadata_primary.txt.gz"

# Unzip required files (keep .gz)
for d in filtered_feature_bc_matrix raw_feature_bc_matrix; do
  for f in matrix.mtx barcodes.tsv features.tsv; do
    if [[ -f "${DL_ROOT}/${d}/${f}.gz" && ! -f "${DL_ROOT}/${d}/${f}" ]]; then
      gunzip -k "${DL_ROOT}/${d}/${f}.gz"
    fi
  done
 done

if [[ -f "${META_ROOT}/${GEO}_cell_metadata_primary.txt.gz" && ! -f "${META_ROOT}/${GEO}_cell_metadata_primary.txt" ]]; then
  gunzip -k "${META_ROOT}/${GEO}_cell_metadata_primary.txt.gz"
fi

# Export env vars for embedded python
export GEO DL_ROOT META_ROOT SPLIT_ROOT

# Split matrices by primary samples
python3 - <<'PY'
import os
from pathlib import Path
import pandas as pd

GEO = os.environ.get('GEO', 'GSE277083')
DL_ROOT = Path(os.environ['DL_ROOT'])
META_ROOT = Path(os.environ['META_ROOT'])
SPLIT_ROOT = Path(os.environ['SPLIT_ROOT'])

meta_path = META_ROOT / f"{GEO}_cell_metadata_primary.txt"
if not meta_path.exists():
    raise SystemExit(f"Missing metadata: {meta_path}")

meta = pd.read_csv(meta_path, sep='\t', usecols=['Cell','Sample'])
meta['suffix'] = meta['Cell'].str.split('-').str[-1]

# Ensure suffix -> sample is unique
suffix_to_sample = {}
for suffix, grp in meta.groupby('suffix'):
    samples = grp['Sample'].unique().tolist()
    if len(samples) != 1:
        raise SystemExit(f"Suffix {suffix} maps to multiple samples: {samples}")
    suffix_to_sample[suffix] = samples[0]

primary_samples = sorted(set(suffix_to_sample.values()))
print("Primary samples:", primary_samples)

sample_to_suffix = {v:k for k,v in suffix_to_sample.items()}

import numpy as np

def split_matrix(kind):
    print(f"\n=== Splitting {kind} ===")
    in_dir = DL_ROOT / f"{kind}_feature_bc_matrix"
    barcodes_path = in_dir / "barcodes.tsv"
    features_path = in_dir / "features.tsv"
    matrix_path = in_dir / "matrix.mtx"

    if not barcodes_path.exists() or not features_path.exists() or not matrix_path.exists():
        raise SystemExit(f"Missing input in {in_dir}")

    # Prepare output dirs
    for sample in primary_samples:
        out_dir = SPLIT_ROOT / sample / f"{kind}_feature_bc_matrix"
        out_dir.mkdir(parents=True, exist_ok=True)
        feat_out = out_dir / "features.tsv"
        if not feat_out.exists():
            feat_out.write_bytes(features_path.read_bytes())

    # Pass 1: barcodes and index mapping
    sample_idx = {s:i for i,s in enumerate(primary_samples)}
    col_sample = []
    col_newidx = []
    counts = [0]*len(primary_samples)

    barcode_out_files = {}
    for sample in primary_samples:
        out_dir = SPLIT_ROOT / sample / f"{kind}_feature_bc_matrix"
        barcode_out_files[sample] = open(out_dir / "barcodes.tsv", "w", buffering=1024*1024)

    with open(barcodes_path, 'r') as f:
        for line in f:
            bc = line.strip()
            suffix = bc.split('-')[-1]
            sample = suffix_to_sample.get(suffix)
            if sample is None:
                col_sample.append(-1)
                col_newidx.append(0)
                continue
            sid = sample_idx[sample]
            counts[sid] += 1
            col_sample.append(sid)
            col_newidx.append(counts[sid])
            barcode_out_files[sample].write(bc + "\n")

    for fh in barcode_out_files.values():
        fh.close()

    col_sample = np.array(col_sample, dtype=np.int16)
    col_newidx = np.array(col_newidx, dtype=np.int32)

    print("Barcode counts per sample:")
    for s in primary_samples:
        print(f"  {s}: {counts[sample_idx[s]]}")

    # Pass 2: count nnz per sample
    print("Counting nnz per sample (first pass over matrix)...")
    nnz_counts = [0]*len(primary_samples)
    with open(matrix_path, 'r') as f:
        line = f.readline()
        if not line.startswith('%%MatrixMarket'):
            raise SystemExit('Invalid MatrixMarket header')
        # skip comments, read dimensions
        for line in f:
            if line.startswith('%'):
                continue
            parts = line.strip().split()
            if len(parts) == 3:
                nrows = int(parts[0])
                break
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            j = int(parts[1])
            sid = col_sample[j-1]
            if sid >= 0:
                nnz_counts[sid] += 1

    print("NNZ counts per sample:")
    for s in primary_samples:
        print(f"  {s}: {nnz_counts[sample_idx[s]]}")

    # Pass 3: write per-sample matrix
    print("Writing per-sample matrix files (second pass over matrix)...")
    out_files = {}
    for s in primary_samples:
        out_dir = SPLIT_ROOT / s / f"{kind}_feature_bc_matrix"
        out_path = out_dir / "matrix.mtx"
        out_files[s] = open(out_path, 'w', buffering=1024*1024)
        out_files[s].write("%%MatrixMarket matrix coordinate integer general\n")
        out_files[s].write(f"{nrows} {counts[sample_idx[s]]} {nnz_counts[sample_idx[s]]}\n")

    with open(matrix_path, 'r') as f:
        line = f.readline()
        if not line.startswith('%%MatrixMarket'):
            raise SystemExit('Invalid MatrixMarket header')
        for line in f:
            if line.startswith('%'):
                continue
            parts = line.strip().split()
            if len(parts) == 3:
                break
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            i = parts[0]
            j = int(parts[1])
            val = parts[2]
            sid = col_sample[j-1]
            if sid >= 0:
                sname = primary_samples[sid]
                newj = col_newidx[j-1]
                out_files[sname].write(f"{i} {newj} {val}\n")

    for fh in out_files.values():
        fh.close()

    print(f"Done {kind}.")

split_matrix('filtered')
split_matrix('raw')

print("\nAll done. Split root:", SPLIT_ROOT)
PY

# Clean output root if requested
if [[ "${CLEAN}" == "true" ]]; then
  rm -rf "${OUT_ROOT}"
fi

# Refuse to overwrite non-empty output root without --clean
if [[ -d "${OUT_ROOT}" ]] && [[ -n "$(ls -A "${OUT_ROOT}" 2>/dev/null)" ]]; then
  echo "[ERROR] Output root not empty: ${OUT_ROOT}" >&2
  echo "Use --clean to overwrite." >&2
  exit 1
fi

mkdir -p "${OUT_ROOT}"

# Move split samples into output root
if compgen -G "${SPLIT_ROOT}/*" > /dev/null; then
  mv "${SPLIT_ROOT}"/* "${OUT_ROOT}/"
fi

# Ensure .gz versions exist for scanpy.read_10x_mtx
for sample_dir in "${OUT_ROOT}"/*; do
  [[ -d "${sample_dir}" ]] || continue
  for kind in raw_feature_bc_matrix filtered_feature_bc_matrix; do
    dir="${sample_dir}/${kind}"
    [[ -d "${dir}" ]] || continue
    for f in matrix.mtx barcodes.tsv features.tsv; do
      if [[ -f "${dir}/${f}" && ! -f "${dir}/${f}.gz" ]]; then
        gzip -k "${dir}/${f}"
      fi
    done
  done
 done

echo "Prepared primary-only Cell Ranger inputs at: ${OUT_ROOT}"
