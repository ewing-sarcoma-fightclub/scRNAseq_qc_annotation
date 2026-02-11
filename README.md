# Ewing Sarcoma scRNA-seq Pipeline

## Quick Install (recommended)
1. Clone the repo:
```bash
git clone https://github.com/ewing-sarcoma-fightclub/scRNAseq_qc_annotation.git
cd scRNAseq_qc_annotation
```
2. Install `micromamba` (or `mamba`/`conda`), then run:
```bash
./bin/setup_envs.sh
```
This will:
- Create/update Python and R environments
- Install required R packages in the correct order
- Write `env/config.local.env` with pinned `R_BIN`/`PYTHON_BIN`
- Pin `R_LIBS_USER=NULL` so pipeline runs use the env library, not user-level R libs

3. Edit config (`setup_envs.sh` already creates `env/config.local.env`):
```bash
[[ -f env/config.local.env ]] || cp env/config.env env/config.local.env
```
Set at minimum:
- `FASTQ_ROOT` (if using FASTQs)
- `AMBIQUANT_REPO` (optional, for ambient contamination tracking)
- `SEURAT_FINAL_REQUIRE_DROPLETQC=false` (if no BAM files)

4. Run pipeline:
```bash
./bin/run_all.sh
```

Matrix-only input (no FASTQs, no BAM):
```bash
./bin/prepare_geo_primary.sh --clean
./bin/run_all.sh --skip-cellranger --cellranger-root ./outputs/cellranger --qc-out ./outputs/qc --config ./env/config.local.env
```

## R Dependency Order
R dependencies are installed by `r/install_pipeline_packages.R` in this order:
1. CRAN baseline: `future`, `mclust`, `remotes`, `statmod`
2. Bioconductor: `AUCell`, `DropletUtils`
3. CRAN stack needed by DropletQC on some systems: `nloptr`, `lme4`, `pbkrtest`, `car`, `rstatix`, `ggpubr`
4. GitHub-only required packages:
   - `chris-mcginnis-ucsf/DoubletFinder`
   - `powellgenomicslab/DropletQC`
5. Optional Azimuth references via `SeuratData::InstallData()`

The installer is idempotent and ends with a required-package verification table.

## Manual Environment Setup (if not using `setup_envs.sh`)
### Linux
```bash
micromamba create -n ewing-scrna-py -f env/envs/python.lock.yml
micromamba create -n ewing-scrna-r  -f env/envs/r.lock.yml
micromamba run -n ewing-scrna-r env -u R_LIBS R_LIBS_USER=NULL Rscript r/install_pipeline_packages.R
```

### macOS (Apple Silicon)
```bash
micromamba create -n ewing-scrna-py -f env/envs/python.yml
micromamba create -n ewing-scrna-r  -f env/envs/r.macos.yml
micromamba run -n ewing-scrna-r env -u R_LIBS R_LIBS_USER=NULL Rscript r/install_pipeline_packages.R
```

### macOS (Intel)
```bash
micromamba create -n ewing-scrna-py -f env/envs/python.yml
micromamba create -n ewing-scrna-r  -f env/envs/r.yml
micromamba run -n ewing-scrna-r env -u R_LIBS R_LIBS_USER=NULL Rscript r/install_pipeline_packages.R
```

If Bioconductor builds fail on macOS with compiler errors:
```bash
xcode-select --install
```
Then set `~/.R/Makevars`:
```make
CC=clang
CXX=clang++
CC17=clang
CXX17=clang++
```

## Optional Azimuth Reference Preload
```bash
./bin/setup_envs.sh --install-azimuth-refs --azimuth-refs lungref,liverref
```

## Quick Start Commands
Run full pipeline (Cell Ranger + QC):
```bash
./bin/run_all.sh
```

Run QC only from existing Cell Ranger outputs:
```bash
./bin/run_all.sh --skip-cellranger --cellranger-root /path/to/cellranger_out --qc-out ./outputs/qc --config ./env/config.local.env
```

Run QC entrypoint directly:
```bash
./bin/pipeline_QC_after_cellranger.sh --root /path/to/cellranger_out --out ./outputs/qc --config ./env/config.local.env
```

## Output Layout
- `outputs/refdata/` Cell Ranger reference
- `outputs/cellranger/` Cell Ranger per-sample outputs
- `outputs/qc/` QC outputs (EmptyDrops, SoupX, DoubletFinder, DropletQC, Seurat QC, optional AmbiQuant)
- `outputs/annotation_integration/` merged annotation/integration outputs

## Repo Layout
- `bin/` orchestration and setup scripts
- `scripts/` per-tool loop wrappers
- `r/` R analysis/install scripts
- `python/` Python utilities
- `resources/` static resources and marker tables
- `env/` env definitions and config templates

## Notes
- Cell Ranger is proprietary and must be installed separately.
- `env/config.local.env` is preferred over `env/config.env`.
- If `AMBIQUANT_REPO` is empty, AmbiQuant steps are skipped.
