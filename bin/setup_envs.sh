#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

usage() {
  cat <<'USAGE'
Usage: ./bin/setup_envs.sh [options]

Options:
  --mamba-bin <path>         Path to micromamba/mamba/conda binary
  --no-lock                  Use env/envs/*.yml instead of env/envs/*.lock.yml
  --skip-r                   Skip creating the R environment
  --skip-py                  Skip creating the Python environment
  --no-github                Skip GitHub-only R packages (DoubletFinder, DropletQC)
  --install-azimuth-refs     Pre-install Azimuth references via SeuratData
  --azimuth-refs <list>      Comma list of Azimuth refs (default set in script)
  --config-in <path>         Input config template (default: ./env/config.env)
  --config-out <path>        Output config (default: ./env/config.local.env)
  -h, --help                 Show this help

Examples:
  ./bin/setup_envs.sh
  ./bin/setup_envs.sh --install-azimuth-refs --azimuth-refs lungref,liverref
USAGE
}

MAMBA_BIN="${MAMBA_BIN:-}"
USE_LOCK=true
SKIP_R=false
SKIP_PY=false
INSTALL_GITHUB=true
INSTALL_AZIMUTH_REFS=false
AZIMUTH_REFERENCES="${AZIMUTH_REFERENCES:-pbmcref,bonemarrowref,lungref,adiposeref,fetusref,liverref}"
CONFIG_IN="${CONFIG_IN:-${ROOT_DIR}/env/config.env}"
CONFIG_OUT="${CONFIG_OUT:-${ROOT_DIR}/env/config.local.env}"
R_INSTALL_SCRIPT="${ROOT_DIR}/r/install_pipeline_packages.R"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mamba-bin)
      MAMBA_BIN="$2"; shift 2 ;;
    --no-lock)
      USE_LOCK=false; shift ;;
    --skip-r)
      SKIP_R=true; shift ;;
    --skip-py)
      SKIP_PY=true; shift ;;
    --no-github)
      INSTALL_GITHUB=false; shift ;;
    --install-azimuth-refs)
      INSTALL_AZIMUTH_REFS=true; shift ;;
    --azimuth-refs)
      AZIMUTH_REFERENCES="$2"; shift 2 ;;
    --config-in)
      CONFIG_IN="$2"; shift 2 ;;
    --config-out)
      CONFIG_OUT="$2"; shift 2 ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "$MAMBA_BIN" ]]; then
  if command -v micromamba >/dev/null 2>&1; then
    MAMBA_BIN="micromamba"
  elif command -v mamba >/dev/null 2>&1; then
    MAMBA_BIN="mamba"
  elif command -v conda >/dev/null 2>&1; then
    MAMBA_BIN="conda"
  else
    echo "No micromamba/mamba/conda found in PATH." >&2
    exit 1
  fi
fi

R_ENV_NAME="ewing-scrna-r"
PY_ENV_NAME="ewing-scrna-py"
OS_NAME="$(uname -s)"
OS_ARCH="$(uname -m)"

is_macos_arm() {
  [[ "${OS_NAME}" == "Darwin" && "${OS_ARCH}" == "arm64" ]]
}

choose_env_file() {
  local base="$1"
  if [[ "${OS_NAME}" == "Darwin" ]]; then
    if [[ "$base" == "r" ]]; then
      if is_macos_arm && [[ -f "${ROOT_DIR}/env/envs/r.macos.yml" ]]; then
        echo "${ROOT_DIR}/env/envs/r.macos.yml"
        return
      fi
      if [[ -f "${ROOT_DIR}/env/envs/r.yml" ]]; then
        echo "${ROOT_DIR}/env/envs/r.yml"
        return
      fi
    fi
    if [[ "$base" == "python" && -f "${ROOT_DIR}/env/envs/python.yml" ]]; then
      echo "${ROOT_DIR}/env/envs/python.yml"
      return
    fi
  fi
  if $USE_LOCK && [[ -f "${ROOT_DIR}/env/envs/${base}.lock.yml" ]]; then
    echo "${ROOT_DIR}/env/envs/${base}.lock.yml"
  else
    echo "${ROOT_DIR}/env/envs/${base}.yml"
  fi
}

env_exists() {
  local env_name="$1"
  "$MAMBA_BIN" env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -Fxq "$env_name"
}

env_prefix() {
  local env_name="$1"
  "$MAMBA_BIN" env list | awk -v env="$env_name" 'NF && $1 !~ /^#/ {if ($1==env) print $NF}' | tail -n 1
}

create_or_update() {
  local env_name="$1"
  local env_file="$2"
  if env_exists "$env_name"; then
    "$MAMBA_BIN" env update -n "$env_name" -f "$env_file"
  else
    "$MAMBA_BIN" env create -n "$env_name" -f "$env_file"
  fi
}

if ! $SKIP_PY; then
  PY_ENV_FILE="$(choose_env_file python)"
  echo "[setup] Creating/updating Python env from: ${PY_ENV_FILE}"
  create_or_update "$PY_ENV_NAME" "$PY_ENV_FILE"
fi

if ! $SKIP_R; then
  R_ENV_FILE="$(choose_env_file r)"
  echo "[setup] Creating/updating R env from: ${R_ENV_FILE}"
  create_or_update "$R_ENV_NAME" "$R_ENV_FILE"
fi

R_PREFIX=""
PY_PREFIX=""
if ! $SKIP_R; then
  R_PREFIX="$(env_prefix "$R_ENV_NAME")"
fi
if ! $SKIP_PY; then
  PY_PREFIX="$(env_prefix "$PY_ENV_NAME")"
fi

if [[ -n "$R_PREFIX" ]]; then
  R_BIN="${R_PREFIX}/bin/Rscript"
fi
if [[ -n "$PY_PREFIX" ]]; then
  PY_BIN="${PY_PREFIX}/bin/python"
fi

if ! $SKIP_R; then
  if [[ ! -f "${R_INSTALL_SCRIPT}" ]]; then
    echo "[setup] R install script not found: ${R_INSTALL_SCRIPT}" >&2
    exit 1
  fi
  R_ENV_BIN_DIR="$(dirname "$R_BIN")"
  R_INSTALL_ARGS=()
  if ! $INSTALL_GITHUB; then
    R_INSTALL_ARGS+=(--skip-github)
  fi
  if $INSTALL_AZIMUTH_REFS; then
    AZIMUTH_REFERENCES_CLEAN="$(echo "$AZIMUTH_REFERENCES" | tr -d ' ')"
    echo "[setup] Installing Azimuth references: ${AZIMUTH_REFERENCES_CLEAN}"
    R_INSTALL_ARGS+=(--install-azimuth-refs --azimuth-refs "${AZIMUTH_REFERENCES_CLEAN}")
  fi
  echo "[setup] Installing required R packages for pipeline"
  if [[ ${#R_INSTALL_ARGS[@]} -gt 0 ]]; then
    PATH="${R_ENV_BIN_DIR}:${PATH}" env -u R_LIBS R_LIBS_USER="NULL" "$R_BIN" "${R_INSTALL_SCRIPT}" "${R_INSTALL_ARGS[@]}"
  else
    PATH="${R_ENV_BIN_DIR}:${PATH}" env -u R_LIBS R_LIBS_USER="NULL" "$R_BIN" "${R_INSTALL_SCRIPT}"
  fi
fi

if [[ -f "$CONFIG_IN" ]]; then
  echo "[setup] Writing config to: ${CONFIG_OUT}"
  cp "$CONFIG_IN" "$CONFIG_OUT"
  if [[ -n "${PY_BIN:-}" ]]; then
    if grep -q '^PYTHON_BIN=' "$CONFIG_OUT"; then
      sed -i.bak "s|^PYTHON_BIN=.*|PYTHON_BIN=\"${PY_BIN}\"|" "$CONFIG_OUT" && rm -f "${CONFIG_OUT}.bak"
    else
      echo "PYTHON_BIN=\"${PY_BIN}\"" >> "$CONFIG_OUT"
    fi
  fi
  if [[ -n "${R_BIN:-}" ]]; then
    if grep -q '^R_BIN=' "$CONFIG_OUT"; then
      sed -i.bak "s|^R_BIN=.*|R_BIN=\"${R_BIN}\"|" "$CONFIG_OUT" && rm -f "${CONFIG_OUT}.bak"
    else
      echo "R_BIN=\"${R_BIN}\"" >> "$CONFIG_OUT"
    fi
  fi
  if grep -q '^MICROMAMBA_BIN=' "$CONFIG_OUT"; then
    sed -i.bak "s|^MICROMAMBA_BIN=.*|MICROMAMBA_BIN=\"${MAMBA_BIN}\"|" "$CONFIG_OUT" && rm -f "${CONFIG_OUT}.bak"
  else
    echo "MICROMAMBA_BIN=\"${MAMBA_BIN}\"" >> "$CONFIG_OUT"
  fi
  if grep -q '^R_LIBS_USER=' "$CONFIG_OUT"; then
    sed -i.bak "s|^R_LIBS_USER=.*|R_LIBS_USER=\"NULL\"|" "$CONFIG_OUT" && rm -f "${CONFIG_OUT}.bak"
  else
    echo "R_LIBS_USER=\"NULL\"" >> "$CONFIG_OUT"
  fi
else
  echo "[setup] Config template not found: ${CONFIG_IN} (skipping config write)"
fi

echo "[setup] Done."
