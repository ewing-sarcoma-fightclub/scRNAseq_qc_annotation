# Ewing Sarcoma scRNA-seq Pipeline

**Quick Install**
1. Get the code:
```bash
git clone https://github.com/ewing-sarcoma-fightclub/scRNAseq_qc_annotation.git
cd scRNAseq_qc_annotation
```
2. Install micromamba or conda/mamba, then create environments:
```bash
micromamba create -n ewing-scrna-py -f env/envs/python.lock.yml
micromamba create -n ewing-scrna-r  -f env/envs/r.lock.yml
```
3. Install required R GitHub packages (inside the R env):
```bash
micromamba activate ewing-scrna-r
Rscript -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='https://cloud.r-project.org'); remotes::install_github('chris-mcginnis-ucsf/DoubletFinder'); remotes::install_github('powellgenomicslab/DropletQC')"
```
If not using lockfiles, ensure these are installed too:
```bash
Rscript -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install('AUCell')"
Rscript -e "if (!requireNamespace('mclust', quietly=TRUE)) install.packages('mclust', repos='https://cloud.r-project.org')"
```
4. Create local config:
```bash
cp env/config.env env/config.local.env
```
Edit `env/config.local.env` and set:
- `FASTQ_ROOT` (if using FASTQs)
- `AMBIQUANT_REPO` (to track ambient contamination)
- `PYTHON_BIN` / `R_BIN` (paths to env binaries)
- `SEURAT_FINAL_REQUIRE_DROPLETQC=false` (if no BAMs)
5. Run:
```bash
./bin/run_all.sh
```
Matrix-only inputs (no FASTQs and BAMs):
```bash
./bin/prepare_geo_primary.sh --clean
./bin/run_all.sh --skip-cellranger --cellranger-root ./outputs/cellranger --qc-out ./outputs/qc --config ./env/config.local.env
```
Notes:
- Cell Ranger is proprietary and must be installed separately.
- Outputs are written under `./outputs/`.
- Azimuth references are large; pre-install them if you plan to use label transfer:
```bash
./bin/setup_envs.sh --install-azimuth-refs --azimuth-refs lungref,liverref
```


## Quick start
1) Edit `env/config.env` and set `FASTQ_ROOT` and `AMBIQUANT_REPO`.
2) Create environments (see `env/envs/`).
3) Run full pipeline (requires Cell Ranger installed and licensed):

```bash
./bin/run_all.sh
```

To skip Cell Ranger and run QC on an existing Cell Ranger output root:

```bash
./bin/run_all.sh --skip-cellranger --cellranger-root /path/to/cellranger_out --qc-out ./outputs/qc --config ./env/config.local.env
```

To prepare GEO GSE277083 primary-only matrix inputs for this pipeline (no FASTQs/BAMs):

```bash
./bin/prepare_geo_primary.sh --clean
./bin/run_all.sh --skip-cellranger --cellranger-root ./outputs/cellranger --qc-out ./outputs/qc --config ./env/config.local.env
```

Or to run only QC directly (recommended to use `env/config.local.env` if present):

```bash
./bin/pipeline_QC_after_cellranger.sh --root /path/to/cellranger_out --out ./outputs/qc --config ./env/config.local.env
```

If `env/config.local.env` exists, `run_all.sh` and `pipeline_QC_after_cellranger.sh`
will use it automatically. Prefer it over `env/config.env` because it pins the correct
`R_BIN`/`PYTHON_BIN` for your environment.
For final runs, use `env/config.production.env` (a locked snapshot of `env/config.env`
+ `env/config.local.env`) by passing `--config`.

## Repo layout
- `bin/` entrypoints and orchestration scripts
- `scripts/` per-tool loop wrappers
- `r/` R analysis scripts
- `python/` Python utilities
- `resources/` static inputs (CellMarker, signatures, grouping rules, doublet tables)
- `env/` configs + conda/mamba environment files
- `docs/` tutorials and walkthroughs (see `docs/pipeline_tutorial_minimal.md`)

## One-command setup 
This creates/updates the Python + R environments, installs GitHub-only R packages,
and writes `env/config.local.env` with the correct `R_BIN`/`PYTHON_BIN` paths.

```bash
./bin/setup_envs.sh
```

To pre-install Azimuth references (large downloads), run:

```bash
./bin/setup_envs.sh --install-azimuth-refs --azimuth-refs lungref,liverref
```

## Output layout 
- `outputs/refdata/` Cell Ranger reference
- `outputs/cellranger/` Cell Ranger per-sample outputs
- `outputs/qc/` QC outputs (EmptyDrops, Seurat metadata, SoupX, DoubletFinder, DropletQC, AmbiQuant, Seurat QC)
- `outputs/annotation_integration/` Annotation + integration outputs

## Matrix-only inputs (no FASTQs/BAMs)
If you start from matrices only (e.g., GEO), you can skip Cell Ranger and DropletQC.
Example dataset used in this project: Ewing sarcoma scRNA-seq from GEO `GSE277083` (see the associated publication for study details).
Use the helper script to download/clean the GEO inputs into the expected `outputs/cellranger/` layout:

```bash
./bin/prepare_geo_primary.sh --clean
```

Example for GSE277083 primary-only inputs:

```bash
./bin/prepare_geo_primary.sh --clean
./bin/run_all.sh --skip-cellranger --cellranger-root ./outputs/cellranger --qc-out ./outputs/qc --config ./env/config.local.env
```

Also set in `env/config.local.env`:
- `SEURAT_FINAL_REQUIRE_DROPLETQC=false`

## Environments
This pipeline uses separate Python and R environments, and a separate AmbiQuant repo/venv.

### Python (AmbiQuant + QC utilities)
Create the Python environment from `env/envs/python.yml`:

```bash
mamba env create -f env/envs/python.yml
mamba activate ewing-scrna-py
```

Note for macOS (Apple Silicon): `dropkick` is not available on `osx-arm64` via conda.
It is optional for this pipeline and is omitted from `env/envs/python.yml`. If you want it,
install via pip after creating the env:

```bash
pip install dropkick
```

### R (SoupX, Seurat, DropletQC, DoubletFinder)
Create the R environment from `env/envs/r.yml`:

```bash
mamba env create -f env/envs/r.yml
mamba activate ewing-scrna-r
```

If any R packages are missing after env creation, install them from R (Bioconductor or GitHub).

Note for macOS (Apple Silicon): `bioconductor-dropletutils` is not available via conda on `osx-arm64`.
Use the mac-specific env file, then install DropletUtils from Bioconductor:

```bash
mamba env create -f env/envs/r.macos.yml
mamba activate ewing-scrna-r
Rscript -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install('DropletUtils')"
```

#### GitHub-only R packages (required)
These are not guaranteed to be available via conda/Bioconductor. Install them in the R env:

```bash
mamba activate ewing-scrna-r
Rscript -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='https://cloud.r-project.org'); remotes::install_github('chris-mcginnis-ucsf/DoubletFinder'); remotes::install_github('powellgenomicslab/DropletQC')"
```

### AmbiQuant
Clone the AmbiQuant repo and set `AMBIQUANT_REPO` in `env/config.env`. If you use a micromamba/conda env,
set `AMBIQUANT_ENV` and `MICROMAMBA_BIN` as well. If `AMBIQUANT_REPO` is empty, AmbiQuant steps are skipped.


## Reproducible environments 
We provide pinned exports in `env/envs/python.lock.yml` and `env/envs/r.lock.yml` generated from a working install.
These are OS-specific but useful for reproducing exact versions.

```bash
micromamba create -n ewing-scrna-py -f env/envs/python.lock.yml
micromamba create -n ewing-scrna-r -f env/envs/r.lock.yml
```
