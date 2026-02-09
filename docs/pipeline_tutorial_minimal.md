# Minimal scRNA-seq pipeline (Cell Ranger onward)

Assumptions (kept simple for beginners):
- All tools and R packages are already installed and working.
- All input/output directories already exist.
- Each sample is a folder under `FASTQ_ROOT`, and Cell Ranger outputs go under `CELLRANGER_OUT/<sample>/outs/`.

If you start from matrix-only inputs (no FASTQs/BAMs), skip Cell Ranger and DropletQC.
For GSE277083 primary-only inputs in this repo, use `bin/prepare_geo_primary.sh` and set `require_dropletqc=false` in the final step.

Example (matrix-only, no BAMs):
```bash
./bin/prepare_geo_primary.sh --clean
./bin/run_all.sh --skip-cellranger --cellranger-root ./outputs/cellranger --qc-out ./outputs/qc --config ./env/config.local.env
```
Also set:
- `REQUIRE_DROPLETQC="false"` (or `SEURAT_FINAL_REQUIRE_DROPLETQC=false` in `env/config.local.env`)


## Required R packages
- Seurat, SeuratData, Azimuth
- SoupX
- DropletUtils, DropletQC
- DoubletFinder
- Matrix, ggplot2, readxl, HGNChelper

### Version constraints (recommended)
Pinned versions used in this pipeline (from `env/envs/r.lock.yml`):
- Seurat 5.4.0
- Azimuth 0.5.0
- SoupX 1.6.2
- ggplot2 4.0.1
- HGNChelper 0.8.15
- Matrix 1.7_4

GitHub-only packages should be pinned to a commit:
- DoubletFinder (install with a commit hash)
- DropletQC (install with a commit hash)

Below is a single copy/paste script. It avoids `args` inside R and uses explicit variables instead.
Comments explain **why** each step exists.

```bash
#!/bin/bash
set -euo pipefail

# set paths once so every step uses the same locations.
FASTQ_ROOT="/path/to/fastqs"
REF_DIR="/path/to/refdata-gex-GRCh38-2024-A"
OUTPUTS_ROOT="/path/to/outputs"
CELLRANGER_OUT="${OUTPUTS_ROOT}/cellranger"
QC_OUT="${OUTPUTS_ROOT}/qc"
USE_SOUPX="true"
REQUIRE_EMPTYDROPS="true"
REQUIRE_DROPLETQC="true"  # set false if no BAMs (possorted_genome_bam.bam) are available
INCLUDE_INTRONS="true"  # true for snRNA-seq, false for scRNA-seq

# ---------------------------
# 1) Cell Ranger (per sample)
# ---------------------------
# Run Cell Ranger to generate raw/filtered matrices.
for sample_dir in "${FASTQ_ROOT}"/*; do
  sample_name=$(basename "${sample_dir}")
  (
    cd "${CELLRANGER_OUT}"
    cellranger count \
      --id="${sample_name}" \
      --transcriptome="${REF_DIR}" \
      --fastqs="${sample_dir}" \
      --include-introns="${INCLUDE_INTRONS}" \
      --create-bam true
  )
done

# ---------------------------
# 2) EmptyDrops (per sample)
# ---------------------------
# Run EmptyDrops to label real barcodes.
for sample_dir in "${CELLRANGER_OUT}"/*; do
  sample_name=$(basename "${sample_dir}")
  out_dir="${QC_OUT}/emptydrops/${sample_name}"

  Rscript - <<RSCRIPT
  suppressPackageStartupMessages({
    library(DropletUtils)
    library(Matrix)
  })

  # hard-code inputs for clarity (no args).
  sample_dir <- "${sample_dir}/outs"
  sample_name <- "${sample_name}"
  out_dir <- "${out_dir}"

# Read the raw matrix (EmptyDrops needs raw droplets).
  raw_dir <- file.path(sample_dir, "raw_feature_bc_matrix")
  sce <- read10xCounts(raw_dir)

# Run EmptyDrops to label empty vs real droplets.
  set.seed(1)
  ed <- emptyDrops(counts(sce))

  # save results so later steps can filter by these calls.
  out_csv <- file.path(out_dir, paste0(sample_name, "_emptydrops_results.csv"))
  write.csv(as.data.frame(ed), out_csv)

  called <- which(!is.na(ed$FDR) & ed$FDR <= 0.01)
  called_barcodes <- colnames(sce)[called]
  called_file <- file.path(out_dir, paste0(sample_name, "_emptydrops_cells.txt"))
  writeLines(called_barcodes, called_file)
RSCRIPT

done

# ---------------------------
# 3) SoupX (per sample)
# ---------------------------
# Run SoupX to correct ambient RNA contamination.
for sample_dir in "${CELLRANGER_OUT}"/*; do
  sample_name=$(basename "${sample_dir}")
  metadata_csv="${QC_OUT}/seurat_metadata/${sample_name}/metadata.csv"
  out_dir="${QC_OUT}/soupx/${sample_name}"

  Rscript - <<RSCRIPT
  library(SoupX)
  library(Seurat)
  library(Matrix)

  # hard-code inputs for clarity (no args).
  sample_dir <- "${sample_dir}/outs"
  sample_name <- "${sample_name}"
  metadata_path <- "${metadata_csv}"
  out_path <- "${out_dir}"
  npcs <- 30
  resolution <- 0.5
  nfeatures <- 2000
  seed <- 1
  if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

# Load filtered + raw matrices (SoupX needs both).
  filtered_dir <- file.path(sample_dir, "filtered_feature_bc_matrix")
  raw_dir <- file.path(sample_dir, "raw_feature_bc_matrix")
  toc <- CreateSeuratObject(Read10X(filtered_dir))
  tod <- Read10X(raw_dir)

  # Build basic clusters (SoupX uses clusters to estimate contamination).
  build_metadata <- function(res) {
    set.seed(seed)
    toc_tmp <- CreateSeuratObject(Read10X(filtered_dir))
    toc_tmp[["percent.mt"]] <- PercentageFeatureSet(toc_tmp, pattern = "^MT-")
    toc_tmp <- NormalizeData(toc_tmp)
    toc_tmp <- FindVariableFeatures(toc_tmp, selection.method = "vst", nfeatures = nfeatures)
    toc_tmp <- ScaleData(toc_tmp)
    toc_tmp <- RunPCA(toc_tmp, npcs = npcs)
    toc_tmp <- FindNeighbors(toc_tmp, dims = 1:npcs)
    toc_tmp <- FindClusters(toc_tmp, resolution = res, algorithm = 4)

    md <- toc_tmp@meta.data
    md$leiden <- md$seurat_clusters
    md$n_genes <- md$nFeature_RNA
    md$n_genes_by_counts <- md$nFeature_RNA
    md$total_counts <- md$nCount_RNA
    md$pct_counts_mt <- md$percent.mt
    list(meta = md)
  }

  if (!file.exists(metadata_path)) {
    metadata_file <- build_metadata(resolution)$meta
    dir.create(dirname(metadata_path), recursive = TRUE, showWarnings = FALSE)
    write.csv(metadata_file, metadata_path)
  } else {
    metadata_file <- read.csv(metadata_path, row.names = 1)
  }

  if (!("leiden" %in% colnames(metadata_file)) && ("seurat_clusters" %in% colnames(metadata_file))) {
    metadata_file$leiden <- metadata_file$seurat_clusters
  }

  if ("leiden" %in% colnames(metadata_file)) {
    n_clusters <- length(unique(metadata_file$leiden))
    if (n_clusters < 2) {
      message("[WARN] Only one cluster detected at resolution ", resolution, "; attempting higher-resolution re-clustering")
      res_try <- unique(c(max(1.0, resolution * 2), max(2.0, resolution * 4)))
      best_meta <- metadata_file
      best_n <- n_clusters
      for (res in res_try) {
        md_new <- build_metadata(res)$meta
        n_new <- length(unique(md_new$leiden))
        message("[INFO] Re-cluster attempt at resolution ", res, " -> ", n_new, " clusters")
        if (n_new > best_n) {
          best_n <- n_new
          best_meta <- md_new
        }
        if (n_new >= 2) break
      }
      metadata_file <- best_meta
      if (best_n >= 2) {
        dir.create(dirname(metadata_path), recursive = TRUE, showWarnings = FALSE)
        write.csv(metadata_file, metadata_path)
      } else {
        message("[WARN] Re-clustering still produced a single cluster; will fall back during SoupX")
      }
    }
  }

  if (!("leiden" %in% colnames(metadata_file))) {
    stop("metadata.csv missing required 'leiden' column")
  }
  if (!all(colnames(toc) %in% rownames(metadata_file))) {
    stop("metadata.csv missing barcodes present in filtered matrix; do not pre-filter before SoupX")
  }

  metadata_file$clusters <- metadata_file$leiden
  toc <- AddMetaData(object = toc, metadata = metadata_file)
  toc$soupgroup <- toc@meta.data[["leiden"]]

  # correct ambient RNA with SoupX.
  get_counts <- function(obj) {
    tryCatch({
      GetAssayData(obj, assay = "RNA", layer = "counts")
    }, error = function(e) {
      GetAssayData(obj, assay = "RNA", slot = "counts")
    })
  }
  toc_counts <- get_counts(toc)
  sc <- SoupChannel(tod, toc_counts, metaData = metadata_file)
  sc <- setClusters(sc, toc$soupgroup)
  n_clusters <- length(unique(toc$soupgroup))
  if (n_clusters < 2) {
    message("[WARN] Only one cluster detected; skipping autoEstCont and using fallback rho=0.05")
    sc <- setContaminationFraction(sc, 0.05, forceAccept = TRUE)
  } else {
    sc <- tryCatch({
      autoEstCont(sc, doPlot = FALSE, forceAccept = TRUE)
    }, error = function(e) {
      message("[WARN] autoEstCont failed: ", conditionMessage(e))
      message("[WARN] Retrying with relaxed tfidfMin=0.5 and soupQuantile=0.1")
      tryCatch({
        autoEstCont(sc, tfidfMin = 0.5, soupQuantile = 0.1, doPlot = FALSE, forceAccept = TRUE)
      }, error = function(e2) {
        message("[WARN] autoEstCont retry failed: ", conditionMessage(e2))
        message("[WARN] Falling back to rho=0.05")
        setContaminationFraction(sc, 0.05, forceAccept = TRUE)
      })
    })
  }
  out <- adjustCounts(sc, roundToInt = TRUE)

  # save the corrected matrix for downstream QC and clustering.
  writeMM(out, file.path(out_path, paste0(sample_name, "_soupx.mtx")))

  summary_path <- file.path(out_path, paste0(sample_name, "_soupx_summary.txt"))
  summary_lines <- c(
    paste0("Sum of counts before: ", sum(toc_counts)),
    paste0("Sum of counts after: ", sum(out)),
    paste0("Pct of cells left: ", round(sum(out) / sum(toc_counts), digits = 3))
  )
  writeLines(summary_lines, summary_path)
RSCRIPT

done

# ---------------------------
# 4) DoubletFinder (per sample)
# ---------------------------
# Run DoubletFinder to call likely doublets.
for sample_dir in "${CELLRANGER_OUT}"/*; do
  sample_name=$(basename "${sample_dir}")
  soupx_mtx="${QC_OUT}/soupx/${sample_name}/${sample_name}_soupx.mtx"
  metadata_csv="${QC_OUT}/seurat_metadata/${sample_name}/metadata.csv"
  out_dir="${QC_OUT}/doubletfinder/${sample_name}"

  Rscript - <<RSCRIPT
  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(ggplot2)
  })

  # hard-code inputs for clarity (no args).
  sample_dir <- "${sample_dir}/outs"
  sample_name <- "${sample_name}"
  soupx_mtx <- "${soupx_mtx}"
  out_dir <- "${out_dir}"
  doublet_rate <- 0.075  # approximate; see pipeline for auto rate from metrics_summary.csv
  # If you know your chemistry, use the matching table in the full pipeline:
  # DOUBLETFINDER_CHEMISTRY=3p_v3p1 (or 3p_v4, 5p_v2, 3p_ht_v3p1)
  npcs <- 30
  metadata_path <- "${metadata_csv}"
  min_features <- 200
  min_cells <- 3
  max_percent_mt <- Inf
  seed <- 1

# Load counts and metadata (DoubletFinder runs on a Seurat object).
  filtered_dir <- file.path(sample_dir, "filtered_feature_bc_matrix")
  raw_filtered <- Read10X(filtered_dir)
  if (is.list(raw_filtered)) raw_filtered <- raw_filtered[[1]]

  soupx_counts <- readMM(soupx_mtx)
  rownames(soupx_counts) <- rownames(raw_filtered)
  colnames(soupx_counts) <- colnames(raw_filtered)

  meta <- read.csv(metadata_path, row.names = 1)

  seu <- CreateSeuratObject(
    counts = soupx_counts,
    project = sample_name,
    meta.data = meta,
    min.cells = min_cells,
    min.features = min_features
  )

  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  if (is.finite(max_percent_mt)) {
    seu <- subset(seu, subset = percent.mt <= max_percent_mt)
  }

# Compute neighbors/clusters (required by DoubletFinder).
  set.seed(seed)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, npcs = npcs)
  seu <- FindNeighbors(seu, dims = 1:npcs)
  seu <- FindClusters(seu, resolution = 0.5)
  seu <- RunUMAP(seu, dims = 1:npcs)

  # estimate the best pK and call doublets.
  df_exports <- getNamespaceExports("DoubletFinder")
  param_sweep <- if ("paramSweep_v3" %in% df_exports) DoubletFinder::paramSweep_v3 else DoubletFinder::paramSweep
  df_call <- if ("doubletFinder_v3" %in% df_exports) DoubletFinder::doubletFinder_v3 else DoubletFinder::doubletFinder

  sweep.res.list <- param_sweep(seu, PCs = 1:npcs, sct = FALSE)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)

  bc <- bcmvn$BCmetric
  pk <- bcmvn$pK
  if (is.data.frame(bc)) bc <- bc[[1]]
  if (is.data.frame(pk)) pk <- pk[[1]]
  bc <- suppressWarnings(as.numeric(bc))
  pk <- suppressWarnings(as.numeric(as.character(pk)))
  if (length(bc) == 0 || all(is.na(bc)) || length(pk) == 0) {
    best_pK <- 0.05
  } else {
    best_pK <- pk[which.max(bc)]
    if (is.na(best_pK)) best_pK <- 0.05
  }

  nExp_poi <- round(doublet_rate * ncol(seu))
  homotypic.prop <- DoubletFinder::modelHomotypic(seu@meta.data$seurat_clusters)
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

  seu <- df_call(
    seu,
    PCs = 1:npcs,
    pN = 0.25,
    pK = best_pK,
    nExp = nExp_poi.adj,
    reuse.pANN = NULL,
    sct = FALSE
  )

  # save DoubletFinder calls for later filtering.
  meta_out <- file.path(out_dir, paste0(sample_name, "_metadata_with_DoubletFinder.csv"))
  write.csv(seu@meta.data, meta_out)

  df_class_col <- grep("DF.classifications", colnames(seu@meta.data), value = TRUE)

  p1 <- DimPlot(seu, group.by = "seurat_clusters", label = TRUE) + ggtitle("Seurat clusters")
  if (length(df_class_col) == 1) {
    p2 <- DimPlot(seu, group.by = df_class_col) + ggtitle("DoubletFinder")
    ggsave(file.path(out_dir, paste0(sample_name, "_DoubletFinder_UMAP.png")), p2, width = 6, height = 5)
  }
  ggsave(file.path(out_dir, paste0(sample_name, "_Clusters_UMAP.png")), p1, width = 6, height = 5)
RSCRIPT

done

# ---------------------------
# 5) DropletQC (per sample)
# ---------------------------
# Run DropletQC to flag damaged cells using nuclear fraction.
for sample_dir in "${CELLRANGER_OUT}"/*; do
  sample_name=$(basename "${sample_dir}")
  bam_path="${sample_dir}/outs/possorted_genome_bam.bam"
  if [[ ! -f "${bam_path}" ]]; then
    echo "[WARN] DropletQC skipped for ${sample_name}: possorted_genome_bam.bam not found."
    continue
  fi
  out_dir="${QC_OUT}/dropletqc/${sample_name}"
  metadata_csv="${QC_OUT}/seurat_metadata/${sample_name}/metadata.csv"

  Rscript - <<RSCRIPT
  suppressPackageStartupMessages({
    library(DropletQC)
    library(DropletUtils)
    library(Matrix)
  })

  # hard-code inputs for clarity (no args).
  outs_dir <- "${sample_dir}/outs"
  sample_name <- "${sample_name}"
  out_dir <- "${out_dir}"
  tiles <- 100
  cores <- parallel::detectCores() - 1
  metadata_path <- "${metadata_csv}"
  metadata_col <- "leiden"

# Compute nuclear fraction (DropletQC needs it).
  nf <- nuclear_fraction_tags(outs = outs_dir, tiles = tiles, cores = cores, verbose = FALSE)

# Align UMI counts (DropletQC compares counts vs nuclear fraction).
  filtered_dir <- file.path(outs_dir, "filtered_feature_bc_matrix")
  sce <- read10xCounts(filtered_dir)
  umi <- Matrix::colSums(counts(sce))

  common <- intersect(rownames(nf), colnames(sce))
  nf <- nf[common, , drop = FALSE]
  umi <- umi[common]

  nf_umi <- data.frame(
    nf = as.numeric(nf$nuclear_fraction),
    umi = as.numeric(umi),
    row.names = common
  )

  # call empty droplets and damaged cells.
  ed <- identify_empty_drops(nf_umi = nf_umi)

  cell_type <- rep("all", nrow(ed))
  names(cell_type) <- rownames(ed)
  meta <- read.csv(metadata_path, row.names = 1)
  meta_vals <- meta[rownames(ed), metadata_col, drop = TRUE]
  meta_vals <- as.character(meta_vals)
  cell_type[!is.na(meta_vals)] <- meta_vals[!is.na(meta_vals)]
  ed$cell_type <- cell_type

  res <- identify_damaged_cells(ed, verbose = FALSE, output_plots = FALSE)

  qc_df <- res[[1]]
  qc_df$barcode <- rownames(qc_df)

  # save DropletQC calls for filtering in the final step.
  out_csv <- file.path(out_dir, paste0(sample_name, "_dropletqc_results.csv"))
  write.csv(qc_df, out_csv, row.names = FALSE)
RSCRIPT

done

# ---------------------------
# 6) Final Seurat filtering + clustering (all samples)
# ---------------------------
# Combine QC calls to produce the final clean dataset.
Rscript - <<RSCRIPT
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# hard-code inputs for clarity (no args).
root_dir <- "${CELLRANGER_OUT}"
qc_root <- "${QC_OUT}"
out_root <- "${QC_OUT}/seurat_qc"

parse_bool <- function(x, default = TRUE) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) return(default)
  val <- tolower(as.character(x))
  if (val %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (val %in% c("false", "f", "0", "no", "n")) return(FALSE)
  return(default)
}

use_soupx <- parse_bool("${USE_SOUPX}", TRUE)
require_emptydrops <- parse_bool("${REQUIRE_EMPTYDROPS}", TRUE)
require_dropletqc <- parse_bool("${REQUIRE_DROPLETQC}", TRUE)
npcs <- 30
resolution <- 0.5
seed <- 1
min_cells <- 3
min_features <- 200
do_regress <- FALSE
do_regress_cell_cycle <- FALSE
max_percent_mt <- Inf
max_percent_ribo <- Inf
basic_min_features <- min_features

safe_cell_cycle_scoring <- function(obj, s_genes, g2m_genes, label = "dataset") {
  if (length(s_genes) == 0 || length(g2m_genes) == 0) return(obj)
  for (nbin in c(24, 10, 5)) {
    res <- tryCatch(
      CellCycleScoring(obj, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE, nbin = nbin),
      error = function(e) e
    )
    if (!inherits(res, "error")) return(res)
  }
  return(obj)
}

sample_dirs <- list.dirs(root_dir, full.names = TRUE, recursive = FALSE)

for (sample_dir in sample_dirs) {
  sample <- basename(sample_dir)

  filtered_dir <- file.path(sample_dir, "outs", "filtered_feature_bc_matrix")
  raw_filtered <- Read10X(filtered_dir)
  if (is.list(raw_filtered)) raw_filtered <- raw_filtered[[1]]

  barcodes <- colnames(raw_filtered)

  counts <- raw_filtered
  soupx_path <- file.path(qc_root, "soupx", sample, paste0(sample, "_soupx.mtx"))
  if (use_soupx) {
    soupx_counts <- readMM(soupx_path)
    rownames(soupx_counts) <- rownames(raw_filtered)
    colnames(soupx_counts) <- colnames(raw_filtered)
    counts <- soupx_counts
  }

  # assemble metadata from all QC steps.
  meta <- data.frame(row.names = barcodes)

  meta_path <- file.path(qc_root, "seurat_metadata", sample, "metadata.csv")
  meta_in <- read.csv(meta_path, row.names = 1)
  meta <- cbind(meta, meta_in[barcodes, , drop = FALSE])
  meta$orig.ident <- sample

  ed_path <- file.path(qc_root, "emptydrops", sample, paste0(sample, "_emptydrops_results.csv"))
  ed <- read.csv(ed_path, row.names = 1)
  ed_call <- (!is.na(ed$FDR)) & (ed$FDR <= 0.01)
  names(ed_call) <- rownames(ed)
  meta$emptydrops_call <- ed_call[barcodes]

  df_path <- file.path(qc_root, "doubletfinder", sample, paste0(sample, "_metadata_with_DoubletFinder.csv"))
  df <- read.csv(df_path, row.names = 1)
  df_cols <- grep("^DF.classifications", colnames(df), value = TRUE)
  if (length(df_cols) == 1) {
    df_call <- tolower(as.character(df[[df_cols[1]]])) == "doublet"
    names(df_call) <- rownames(df)
    meta$doubletfinder_call <- df_call[barcodes]
  } else {
    meta$doubletfinder_call <- NA
  }

  dqc_path <- file.path(qc_root, "dropletqc", sample, paste0(sample, "_dropletqc_results.csv"))
  if (file.exists(dqc_path)) {
    dqc <- read.csv(dqc_path)
    if ("barcode" %in% colnames(dqc)) {
      rownames(dqc) <- dqc$barcode
    }
    if ("cell_status" %in% colnames(dqc)) {
      dqc_pass <- tolower(as.character(dqc$cell_status)) == "cell"
      names(dqc_pass) <- rownames(dqc)
      meta$dropletqc_pass <- dqc_pass[barcodes]
      meta$dropletqc_status <- dqc$cell_status[barcodes]
    } else {
      meta$dropletqc_pass <- NA
    }
  } else {
    if (require_dropletqc) {
      stop(paste0(
        "DropletQC results missing for sample ", sample, ": ", dqc_path,
        ". Run DropletQC (requires possorted_genome_bam.bam) or set require_dropletqc=false."
      ), call.=FALSE)
    }
    meta$dropletqc_pass <- NA
  }

  # Compute basic QC metrics for consistent filtering.
  tmp <- CreateSeuratObject(
    counts = counts,
    meta.data = meta,
    min.cells = 0,
    min.features = 0,
    project = sample
  )
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
  tmp[["percent.rps"]] <- PercentageFeatureSet(tmp, pattern = "^RPS")
  tmp[["percent.rpl"]] <- PercentageFeatureSet(tmp, pattern = "^RPL")
  tmp$percent.ribo <- tmp$percent.rps + tmp$percent.rpl

  meta$nFeature_RNA <- tmp$nFeature_RNA
  meta$nCount_RNA <- tmp$nCount_RNA
  meta$percent.mt <- tmp$percent.mt
  meta$percent.rps <- tmp$percent.rps
  meta$percent.rpl <- tmp$percent.rpl
  meta$percent.ribo <- tmp$percent.ribo

  # apply QC filters: basic gates + EmptyDrops + DropletQC + DoubletFinder.
  basic_qc_keep <- rep(TRUE, length(barcodes))
  names(basic_qc_keep) <- barcodes
  if (is.finite(basic_min_features)) {
    basic_qc_keep <- basic_qc_keep & (meta$nFeature_RNA >= basic_min_features)
  }
  if (is.finite(max_percent_mt)) {
    basic_qc_keep <- basic_qc_keep & (meta$percent.mt <= max_percent_mt)
  }
  if (is.finite(max_percent_ribo)) {
    basic_qc_keep <- basic_qc_keep & (meta$percent.ribo <= max_percent_ribo)
  }

  df_call_full <- rep(FALSE, length(barcodes))
  names(df_call_full) <- barcodes
  if (!all(is.na(meta$doubletfinder_call))) {
    df_call_full[barcodes] <- as.logical(meta$doubletfinder_call)
    df_call_full[is.na(df_call_full)] <- FALSE
  }

  keep <- basic_qc_keep
  if (require_emptydrops) keep <- keep & !is.na(meta$emptydrops_call) & as.logical(meta$emptydrops_call)
  if (require_dropletqc) keep <- keep & !is.na(meta$dropletqc_pass) & as.logical(meta$dropletqc_pass)
  keep <- keep & !df_call_full

  # create the final filtered Seurat object and run clustering.
  seu <- CreateSeuratObject(
    counts = counts,
    meta.data = meta,
    min.cells = min_cells,
    min.features = min_features,
    project = sample
  )
  seu <- subset(seu, cells = barcodes[keep])

  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu[["percent.rps"]] <- PercentageFeatureSet(seu, pattern = "^RPS")
  seu[["percent.rpl"]] <- PercentageFeatureSet(seu, pattern = "^RPL")

  s_genes <- cc.genes.updated.2019$s.genes
  g2m_genes <- cc.genes.updated.2019$g2m.genes
  s_genes <- intersect(s_genes, rownames(seu))
  g2m_genes <- intersect(g2m_genes, rownames(seu))

  do_regress_cell_cycle_sample <- do_regress_cell_cycle
  if (do_regress_cell_cycle_sample) {
    seu <- NormalizeData(seu)
    seu <- safe_cell_cycle_scoring(seu, s_genes, g2m_genes, label = paste0(sample, " (RNA)"))
  }

  set.seed(seed)
  vars_to_regress <- c()
  if (do_regress) {
    vars_to_regress <- c(vars_to_regress, "percent.mt", "percent.rps", "percent.rpl")
  }
  if (do_regress_cell_cycle_sample) {
    vars_to_regress <- c(vars_to_regress, "S.Score", "G2M.Score")
  }
  if (length(vars_to_regress) > 0) {
    seu_proc <- SCTransform(seu, vars.to.regress = vars_to_regress, verbose = TRUE)
  } else {
    seu_proc <- SCTransform(seu, verbose = FALSE)
  }

  if (!("S.Score" %in% colnames(seu_proc@meta.data)) && length(s_genes) > 0 && length(g2m_genes) > 0) {
    DefaultAssay(seu_proc) <- "SCT"
    seu_proc <- safe_cell_cycle_scoring(seu_proc, s_genes, g2m_genes, label = paste0(sample, " (SCT)"))
  }

  s_genes <- cc.genes.updated.2019$s.genes
  g2m_genes <- cc.genes.updated.2019$g2m.genes
  s_genes <- intersect(s_genes, rownames(seu_proc))
  g2m_genes <- intersect(g2m_genes, rownames(seu_proc))
  if (length(s_genes) > 0 && length(g2m_genes) > 0) {
    seu_proc <- safe_cell_cycle_scoring(seu_proc, s_genes, g2m_genes, label = paste0(sample, " (SCT)"))
  }

  seu_proc <- RunPCA(seu_proc, verbose = FALSE)
  seu_proc <- RunUMAP(seu_proc, dims = 1:npcs, verbose = FALSE)
  seu_proc <- FindNeighbors(seu_proc, dims = 1:npcs, verbose = FALSE)
  seu_proc <- FindClusters(seu_proc, resolution = resolution, verbose = FALSE)

  sample_out <- file.path(out_root, sample)
  if (do_regress) {
    saveRDS(seu_proc, file.path(sample_out, paste0(sample, "_sct_regressed_norm.rds")))
  } else {
    saveRDS(seu_proc, file.path(sample_out, paste0(sample, "_sct_norm.rds")))
  }

  meta_out <- cbind(meta, keep_cell = keep)
  write.csv(meta_out, file.path(sample_out, "final_metadata.csv"))
}
RSCRIPT
```

# ---------------------------
# 7) Annotation + optional integration (all samples)
# ---------------------------
# score CellMarker modules + Ewing signature, run Azimuth per sample,
# and assign cluster-level cell types using per-signature percentile evidence.
# Tumor assignment is driven by high Ewing signature only.
# Immune subtypes require strong module + Azimuth support; non-immune lineages
# can be Azimuth-supported when module evidence is weak/absent (with confidence thresholds).
# Incompatible signals become "Unknown". Multiple strong lineages become "Mixed/Doublet".
# You can override Azimuth refs via AZIMUTH_REFERENCES (comma-separated).
AZIMUTH_REFERENCES="pbmcref,bonemarrowref,lungref,adiposeref,fetusref,liverref" \
Rscript /path/to/pipeline/r/Seurat_merge_annotate_integrate.R \
  "${QC_OUT}/seurat_qc" \
  "/path/to/pipeline/resources/Cell_marker_Seq_human.xlsx" \
  "/path/to/pipeline/resources/Aynaud.csv" \
  "${QC_OUT}/../annotation_integration/seurat_annotation_integration" \
  "blood,bone,embryonic brain,dorsal root ganglion,embryonic stem cell,epithelium,lung,bone marrow,cartilage,skeletal muscle,adipose tissue,blood vessel,skin,embryo" \
  5 20 1 TRUE 30 0.5
# Module signatures used (post tissue filter + drop list):
# - Activated langerhans cell
# - Activated memory B cell
# - Activated T cell
# - Adipocyte
# - Adipocyte progenitor cell
# - Alveolar macrophage
# - B cell
# - Basal cell
# - Basal epithelial cell
# - Cancer cell
# - Capillary cell
# - CAR-T cell
# - Cardiomyocyte
# - Cartilage progenitor cell
# - CD16+ monocyte
# - CD4 T cell
# - CD4+ recently activated effector memory or effector T cell (CTL)
# - CD4+ T cell
# - CD8 T cell
# - CD8+ mucosal-associated invariant T cell(CD8+ MAIT)
# - CD8+ recently activated effector memory or effector T cell (CTL)
# - CD8+ T cell
# - Central memory CD4+ T cell
# - Central memory CD8+ T cell
# - Central memory T cell
# - Conventional dendritic cell 1(cDC1)
# - Conventional dendritic cell 2(cDC2)
# - Conventional dendritic cell(cDC)
# - Cranial neural crest cell
# - Dendritic cell lineage
# - Effector CD8 T cell
# - Effector memory CD8+ T cell
# - Effector memory T cell
# - Endothelial cell
# - Epiblast cell
# - Epidermal cell
# - Epithelial cell
# - Exhausted T(Tex) cell
# - Fibrocartilage chondrocyte
# - Glial cell
# - Goblet cell
# - Granulocyte
# - Hematopoietic cell
# - Hematopoietic progenitor cell
# - Hematopoietic stem cell
# - Hypoblast cell
# - Intestinal cell
# - Keratinocyte
# - Lymphoblastoid cell
# - Lymphocyte
# - Lymphoid-primed multipotent progenitor cell(LMPP)
# - M1 macrophage
# - M2 macrophage
# - Macrophage
# - Mast cell
# - Megakaryocyte
# - Melanocyte
# - Memory B cell
# - Memory T cell
# - Merkel cell
# - Mesenchymal cell
# - Mesothelial cell
# - Migratory langerhans cell
# - Myofibroblast
# - Naive B cell
# - Naive CD4 T cell
# - Naive CD8+ T cell
# - Naive T(Th0) cell
# - Natural killer T(NKT) cell
# - Natural regulatory T cell
# - Neuron
# - Neutrophil
# - Normal cell
# - Pan-endothelial cell
# - Pan-trophectoderm cell
# - Panvascular cell
# - Pericyte
# - Plasmacytoid dendritic cell(pDC)
# - Pluripotent stem cell
# - Primitive endoderm cell
# - Pro-neutrophil
# - Regulatory T(Treg) cell
# - Schwann cell
# - Secretory cell
# - Sensory neuron
# - Smooth muscle cell
# - Stem cell
# - Stromal cell
# - T helper 17(Th17) cell
# - T helper 2(Th2) cell
# - T helper(Th) cell
# - Tissue resident memory T(TRM) cell
# - Trophectoderm cell
# - Vascular cell

```



## Full R scripts (execution order)

### EmptyDrops_per_sample.R

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DropletUtils)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: EmptyDrops_per_sample.R <sample_dir> <sample_name> [out_dir]", call.=FALSE)
}

sample_dir <- args[[1]]
sample_name <- args[[2]]
out_dir <- ifelse(length(args) >= 3, args[[3]], file.path(sample_dir, "EmptyDrops"))

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

raw_dir <- file.path(sample_dir, "raw_feature_bc_matrix")
if (!dir.exists(raw_dir)) {
  raw_dir <- file.path(sample_dir, "outs", "raw_feature_bc_matrix")
}
if (!dir.exists(raw_dir)) stop(paste("raw_feature_bc_matrix not found for", sample_name))

# Read raw matrix
sce <- read10xCounts(raw_dir)

# Ensure barcodes are present for downstream outputs
if (is.null(colnames(sce)) && "Barcode" %in% colnames(colData(sce))) {
  colnames(sce) <- as.character(colData(sce)$Barcode)
}

# Run EmptyDrops
set.seed(1)
ed <- emptyDrops(counts(sce))

# Save full results
out_csv <- file.path(out_dir, paste0(sample_name, "_emptydrops_results.csv"))
write.csv(as.data.frame(ed), out_csv)

# Save called cells (FDR <= 0.01) list
called <- which(!is.na(ed$FDR) & ed$FDR <= 0.01)
called_barcodes <- colnames(sce)[called]
called_file <- file.path(out_dir, paste0(sample_name, "_emptydrops_cells.txt"))
writeLines(called_barcodes, called_file)

cat("EmptyDrops done for", sample_name, "\n")
```

### SoupX_per_sample.R

```r
library(SoupX)
library(Seurat)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: SoupX_per_sample.R <sample_dir> <sample_name> [metadata_csv] [out_dir] [npcs] [resolution] [nfeatures] [seed]", call.=FALSE)
}

sample_dir <- args[[1]]
sample_name <- args[[2]]
metadata_path <- ifelse(length(args) >= 3 && nzchar(args[[3]]), args[[3]], "")
out_path <- ifelse(length(args) >= 4 && nzchar(args[[4]]), args[[4]], file.path(sample_dir, "SoupX"))
npcs <- ifelse(length(args) >= 5 && nzchar(args[[5]]), as.integer(args[[5]]), 30)
resolution <- ifelse(length(args) >= 6 && nzchar(args[[6]]), as.numeric(args[[6]]), 0.5)
nfeatures <- ifelse(length(args) >= 7 && nzchar(args[[7]]), as.integer(args[[7]]), 2000)
seed <- ifelse(length(args) >= 8 && nzchar(args[[8]]), as.integer(args[[8]]), 1)

if (!nzchar(metadata_path)) {
  metadata_path <- file.path(out_path, "metadata.csv")
}

# Resolve filtered/raw directories
filtered_dir <- file.path(sample_dir, "filtered_feature_bc_matrix")
if (!dir.exists(filtered_dir)) {
  filtered_dir <- file.path(sample_dir, "outs", "filtered_feature_bc_matrix")
}
raw_dir <- file.path(sample_dir, "raw_feature_bc_matrix")
if (!dir.exists(raw_dir)) {
  raw_dir <- file.path(sample_dir, "outs", "raw_feature_bc_matrix")
}

if (!dir.exists(filtered_dir)) stop(paste("filtered_feature_bc_matrix not found for", sample_name))
if (!dir.exists(raw_dir)) stop(paste("raw_feature_bc_matrix not found for", sample_name))

if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

# 1. Read data
# filtered counts
toc <- CreateSeuratObject(Read10X(filtered_dir))
# raw counts
# For 10x outputs, Read10X can read the directory directly
# (barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz)
tod <- Read10X(raw_dir)

# 2. Metadata (clusters, filtered cells, QC)
build_metadata <- function(res) {
  set.seed(seed)
  toc_tmp <- CreateSeuratObject(Read10X(filtered_dir))
  toc_tmp[["percent.mt"]] <- PercentageFeatureSet(toc_tmp, pattern = "^MT-")
  toc_tmp <- NormalizeData(toc_tmp)
  toc_tmp <- FindVariableFeatures(toc_tmp, selection.method = "vst", nfeatures = nfeatures)
  toc_tmp <- ScaleData(toc_tmp)
  toc_tmp <- RunPCA(toc_tmp, npcs = npcs)
  toc_tmp <- FindNeighbors(toc_tmp, dims = 1:npcs)
  toc_tmp <- FindClusters(toc_tmp, resolution = res, algorithm = 4)

  md <- toc_tmp@meta.data
  md$leiden <- md$seurat_clusters
  md$n_genes <- md$nFeature_RNA
  md$n_genes_by_counts <- md$nFeature_RNA
  md$total_counts <- md$nCount_RNA
  md$pct_counts_mt <- md$percent.mt
  list(meta = md)
}

if (!file.exists(metadata_path)) {
  metadata_file <- build_metadata(resolution)$meta
  dir.create(dirname(metadata_path), recursive = TRUE, showWarnings = FALSE)
  write.csv(metadata_file, metadata_path)
} else {
  metadata_file <- read.csv(metadata_path, row.names = 1)
}

if (!("leiden" %in% colnames(metadata_file)) && ("seurat_clusters" %in% colnames(metadata_file))) {
  metadata_file$leiden <- metadata_file$seurat_clusters
}

if ("leiden" %in% colnames(metadata_file)) {
  n_clusters <- length(unique(metadata_file$leiden))
  if (n_clusters < 2) {
    message("[WARN] Only one cluster detected at resolution ", resolution, "; attempting higher-resolution re-clustering")
    res_try <- unique(c(max(1.0, resolution * 2), max(2.0, resolution * 4)))
    best_meta <- metadata_file
    best_n <- n_clusters
    for (res in res_try) {
      md_new <- build_metadata(res)$meta
      n_new <- length(unique(md_new$leiden))
      message("[INFO] Re-cluster attempt at resolution ", res, " -> ", n_new, " clusters")
      if (n_new > best_n) {
        best_n <- n_new
        best_meta <- md_new
      }
      if (n_new >= 2) break
    }
    metadata_file <- best_meta
    if (best_n >= 2) {
      dir.create(dirname(metadata_path), recursive = TRUE, showWarnings = FALSE)
      write.csv(metadata_file, metadata_path)
    } else {
      message("[WARN] Re-clustering still produced a single cluster; will fall back during SoupX")
    }
  }
}

if (!("leiden" %in% colnames(metadata_file))) {
  stop("metadata.csv missing required 'leiden' column")
}
if (!all(colnames(toc) %in% rownames(metadata_file))) {
  stop("metadata.csv missing barcodes present in filtered matrix; do not pre-filter before SoupX")
}
metadata_file$clusters <- metadata_file$leiden

# Add metadata
toc <- AddMetaData(object = toc, metadata = metadata_file)
toc$soupgroup <- toc@meta.data[["leiden"]]

# Seurat v4/v5 compatibility for counts
get_counts <- function(obj) {
  tryCatch({
    GetAssayData(obj, assay = "RNA", layer = "counts")
  }, error = function(e) {
    GetAssayData(obj, assay = "RNA", slot = "counts")
  })
}

toc_counts <- get_counts(toc)

# 3. Apply SoupX
sc <- SoupChannel(tod, toc_counts, metaData = metadata_file)
sc <- setClusters(sc, toc$soupgroup)
n_clusters <- length(unique(toc$soupgroup))
if (n_clusters < 2) {
  message("[WARN] Only one cluster detected; skipping autoEstCont and using fallback rho=0.05")
  sc <- setContaminationFraction(sc, 0.05, forceAccept = TRUE)
} else {
  sc <- tryCatch({
    autoEstCont(sc, doPlot = FALSE, forceAccept = TRUE)
  }, error = function(e) {
    message("[WARN] autoEstCont failed: ", conditionMessage(e))
    message("[WARN] Retrying with relaxed tfidfMin=0.5 and soupQuantile=0.1")
    tryCatch({
      autoEstCont(sc, tfidfMin = 0.5, soupQuantile = 0.1, doPlot = FALSE, forceAccept = TRUE)
    }, error = function(e2) {
      message("[WARN] autoEstCont retry failed: ", conditionMessage(e2))
      message("[WARN] Falling back to rho=0.05")
      setContaminationFraction(sc, 0.05, forceAccept = TRUE)
    })
  })
}
out <- adjustCounts(sc, roundToInt = TRUE)

print(paste0('Sum of counts before: ', sum(toc_counts)))
print(paste0('Sum of counts after: ', sum(out)))
print(paste0('Pct of cells left: ', round(sum(out) / sum(toc_counts), digits = 3)))

# 4. Save output (genes x cells to match 10x/Seurat expectations)
writeMM(out, file.path(out_path, paste0(sample_name, '_soupx.mtx')))

# Save SoupX summary stats
summary_path <- file.path(out_path, paste0(sample_name, "_soupx_summary.txt"))
summary_lines <- c(
  paste0("Sum of counts before: ", sum(toc_counts)),
  paste0("Sum of counts after: ", sum(out)),
  paste0("Pct of cells left: ", round(sum(out) / sum(toc_counts), digits = 3))
)
writeLines(summary_lines, summary_path)
```

### DoubletFinder_after_SoupX.R

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: DoubletFinder_after_SoupX.R <sample_dir> <sample_name> [soupx_mtx] [out_dir] [doublet_rate] [npcs] [metadata_csv] [min_features] [min_cells] [max_percent_mt] [seed]", call. = FALSE)
}

sample_dir <- args[[1]]
sample_name <- args[[2]]
soupx_mtx <- ifelse(length(args) >= 3 && nzchar(args[[3]]), args[[3]], file.path(sample_dir, "SoupX", paste0(sample_name, "_soupx.mtx")))
out_dir <- ifelse(length(args) >= 4 && nzchar(args[[4]]), args[[4]], file.path(sample_dir, "DoubletFinder"))
doublet_rate_raw <- ifelse(length(args) >= 5 && nzchar(args[[5]]), args[[5]], "")
npcs <- ifelse(length(args) >= 6 && nzchar(args[[6]]), as.integer(args[[6]]), 30)
metadata_path <- ifelse(length(args) >= 7 && nzchar(args[[7]]), args[[7]], file.path(sample_dir, "metadata.csv"))
min_features <- ifelse(length(args) >= 8 && nzchar(args[[8]]), as.integer(args[[8]]), 200)
min_cells <- ifelse(length(args) >= 9 && nzchar(args[[9]]), as.integer(args[[9]]), 3)
max_percent_mt <- ifelse(length(args) >= 10 && nzchar(args[[10]]), as.numeric(args[[10]]), Inf)
seed <- ifelse(length(args) >= 11 && nzchar(args[[11]]), as.integer(args[[11]]), 1)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

get_script_dir <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  return(getwd())
}

normalize_col <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

get_metrics_values <- function(sample_dir) {
  metrics_path <- file.path(sample_dir, "outs", "metrics_summary.csv")
  if (!file.exists(metrics_path)) {
    metrics_path <- file.path(sample_dir, "metrics_summary.csv")
  }
  if (!file.exists(metrics_path)) return(list())
  df <- tryCatch(read.csv(metrics_path, check.names = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(list())
  cols_norm <- normalize_col(names(df))
  get_col <- function(patterns) {
    idx <- which(cols_norm %in% patterns)
    if (length(idx) == 0) return(NULL)
    df[[idx[1]]][1]
  }
  n_cells <- get_col(c("estimatednumberofcells", "estimatedcells", "numberofcells"))
  mult_rate <- get_col(c("estimatedmultipletrate", "estimatedmultipletbrate", "multipletrate", "multipletbrate"))
  chemistry <- get_col(c("chemistry", "chemistrydescription", "librarychemistry"))
  list(n_cells = n_cells, multiplet_rate = mult_rate, chemistry = chemistry)
}

interp_rate <- function(n_cells, table_df) {
  table_df <- table_df[order(table_df$n_cells), , drop = FALSE]
  if (n_cells <= min(table_df$n_cells)) return(table_df$rate[which.min(table_df$n_cells)])
  if (n_cells >= max(table_df$n_cells)) return(table_df$rate[which.max(table_df$n_cells)])
  hi_idx <- which(table_df$n_cells >= n_cells)[1]
  lo_idx <- hi_idx - 1
  x0 <- table_df$n_cells[lo_idx]
  x1 <- table_df$n_cells[hi_idx]
  y0 <- table_df$rate[lo_idx]
  y1 <- table_df$rate[hi_idx]
  y0 + (y1 - y0) * (n_cells - x0) / (x1 - x0)
}

normalize_chemistry <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "", x)
  x
}

select_rate_table <- function(chemistry, script_dir) {
  chem <- normalize_chemistry(chemistry)
  resources_dir <- file.path(script_dir, "..", "resources")
  table_map <- list(
    "3pv3" = "doublet_rate_table_v3.csv",
    "3pv31" = "doublet_rate_table_v3p1.csv",
    "3pv3p1" = "doublet_rate_table_v3p1.csv",
    "nextgem3pv31" = "doublet_rate_table_v3p1.csv",
    "3pv4" = "doublet_rate_table_v4.csv",
    "3pgemx" = "doublet_rate_table_v4.csv",
    "gemx3p" = "doublet_rate_table_v4.csv",
    "5pv2" = "doublet_rate_table_5p_v2.csv",
    "5pgexv2" = "doublet_rate_table_5p_v2.csv",
    "3phtv31" = "doublet_rate_table_3p_ht_v3p1.csv",
    "3phtv3p1" = "doublet_rate_table_3p_ht_v3p1.csv",
    "3pht" = "doublet_rate_table_3p_ht_v3p1.csv"
  )
  if (length(chem) != 1 || !nzchar(chem)) {
    return(file.path(resources_dir, "doublet_rate_table_v3.csv"))
  }
  if (chem %in% names(table_map)) {
    return(file.path(resources_dir, table_map[[chem]]))
  }
  return(file.path(resources_dir, "doublet_rate_table_v3.csv"))
}

resolve_doublet_rate <- function(raw_value, sample_dir) {
  if (nzchar(raw_value)) {
    val_num <- suppressWarnings(as.numeric(raw_value))
    if (is.finite(val_num) && val_num > 0) return(val_num)
    if (tolower(raw_value) != "auto") {
      warning("Invalid doublet_rate '", raw_value, "'. Falling back to auto.")
    }
  }

  metrics <- get_metrics_values(sample_dir)
  n_cells <- suppressWarnings(as.numeric(metrics$n_cells))
  mult_rate <- suppressWarnings(as.numeric(metrics$multiplet_rate))
  is_valid_scalar <- function(x) length(x) == 1 && is.finite(x)
  if (is_valid_scalar(mult_rate)) {
    if (mult_rate > 1) mult_rate <- mult_rate / 100
    if (mult_rate > 0) return(mult_rate)
  }

  table_path <- Sys.getenv("DOUBLETFINDER_RATE_TABLE", "")
  if (!nzchar(table_path)) {
    chem_env <- Sys.getenv("DOUBLETFINDER_CHEMISTRY", "")
    chem_use <- if (nzchar(chem_env)) chem_env else metrics$chemistry
    table_path <- select_rate_table(chem_use, get_script_dir())
  }
  if (file.exists(table_path) && is_valid_scalar(n_cells) && n_cells > 0) {
    tab <- tryCatch(read.csv(table_path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(tab) && all(c("n_cells", "rate") %in% colnames(tab))) {
      tab$rate <- suppressWarnings(as.numeric(tab$rate))
      tab$n_cells <- suppressWarnings(as.numeric(tab$n_cells))
      tab <- tab[is.finite(tab$n_cells) & is.finite(tab$rate), , drop = FALSE]
      if (nrow(tab) > 0) {
        if (any(tab$rate > 1)) tab$rate <- tab$rate / 100
        return(interp_rate(n_cells, tab))
      }
    }
  }

  warning("Doublet rate could not be derived; defaulting to 0.075")
  0.075
}

doublet_rate <- resolve_doublet_rate(doublet_rate_raw, sample_dir)
message("[INFO] DoubletFinder rate for ", sample_name, ": ", signif(doublet_rate, 3))

# ---- input ----
filtered_dir <- file.path(sample_dir, "filtered_feature_bc_matrix")
if (!dir.exists(filtered_dir)) {
  filtered_dir <- file.path(sample_dir, "outs", "filtered_feature_bc_matrix")
}
if (!file.exists(soupx_mtx)) stop(paste("SoupX matrix not found:", soupx_mtx))
if (!dir.exists(filtered_dir)) stop(paste("filtered_feature_bc_matrix not found:", filtered_dir))

# Read filtered 10x to get gene + barcode names
raw_filtered <- Read10X(filtered_dir)
if (is.list(raw_filtered)) {
  raw_filtered <- raw_filtered[[1]]
}

soupx_counts <- readMM(soupx_mtx)
if (!all(dim(soupx_counts) == dim(raw_filtered))) {
  stop("SoupX matrix dimensions do not match filtered_feature_bc_matrix.")
}

rownames(soupx_counts) <- rownames(raw_filtered)
colnames(soupx_counts) <- colnames(raw_filtered)

meta <- NULL
if (file.exists(metadata_path)) {
  meta <- read.csv(metadata_path, row.names = 1)
}

# ---- Seurat object ----
seu <- CreateSeuratObject(
  counts = soupx_counts,
  project = sample_name,
  meta.data = meta,
  min.cells = min_cells,
  min.features = min_features
)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
if (is.finite(max_percent_mt)) {
  seu <- subset(seu, subset = percent.mt <= max_percent_mt)
}
if (ncol(seu) == 0) {
  stop("No cells remaining after basic QC filtering. Adjust min_features/min_cells/max_percent_mt.", call. = FALSE)
}

# Basic processing (no SCT)
set.seed(seed)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = npcs)
seu <- FindNeighbors(seu, dims = 1:npcs)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:npcs)

# ---- DoubletFinder ----
if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  stop("DoubletFinder not installed. Install via: remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')")
}

df_exports <- getNamespaceExports("DoubletFinder")
param_sweep <- if ("paramSweep_v3" %in% df_exports) DoubletFinder::paramSweep_v3 else DoubletFinder::paramSweep
df_call <- if ("doubletFinder_v3" %in% df_exports) DoubletFinder::doubletFinder_v3 else DoubletFinder::doubletFinder

sweep.res.list <- param_sweep(seu, PCs = 1:npcs, sct = FALSE)
sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- DoubletFinder::find.pK(sweep.stats)

# pick pK with max BCmetric (robust to list/data.frame columns)
bc <- bcmvn$BCmetric
pk <- bcmvn$pK
if (is.data.frame(bc)) bc <- bc[[1]]
if (is.data.frame(pk)) pk <- pk[[1]]
bc <- suppressWarnings(as.numeric(bc))
pk <- suppressWarnings(as.numeric(as.character(pk)))
if (length(bc) == 0 || all(is.na(bc)) || length(pk) == 0) {
  warning("BCmetric/pK not available; falling back to pK=0.05")
  best_pK <- 0.05
} else {
  best_pK <- pk[which.max(bc)]
  if (is.na(best_pK)) {
    warning("best_pK NA; falling back to pK=0.05")
    best_pK <- 0.05
  }
}

nExp_poi <- round(doublet_rate * ncol(seu))
homotypic.prop <- DoubletFinder::modelHomotypic(seu@meta.data$seurat_clusters)
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

seu <- df_call(
  seu,
  PCs = 1:npcs,
  pN = 0.25,
  pK = best_pK,
  nExp = nExp_poi.adj,
  reuse.pANN = NULL,
  sct = FALSE
)

# ---- outputs ----
# Save metadata with DoubletFinder columns
meta_out <- file.path(out_dir, paste0(sample_name, "_metadata_with_DoubletFinder.csv"))
write.csv(seu@meta.data, meta_out)

# Summary table comparing Scrublet vs DoubletFinder if Scrublet exists
scrub_col <- NULL
for (c in c("predicted_doublets", "predicted_doublet")) {
  if (c %in% colnames(seu@meta.data)) scrub_col <- c
}

# Find DF classification column
df_class_col <- grep("DF.classifications", colnames(seu@meta.data), value = TRUE)

if (!is.null(scrub_col) && length(df_class_col) == 1) {
  tab <- table(seu@meta.data[[scrub_col]], seu@meta.data[[df_class_col]])
  tab_out <- file.path(out_dir, paste0(sample_name, "_scrublet_vs_doubletfinder.csv"))
  write.csv(as.data.frame.matrix(tab), tab_out)
}

# UMAP plots
p1 <- DimPlot(seu, group.by = "seurat_clusters", label = TRUE) + ggtitle("Seurat clusters")
if (length(df_class_col) == 1) {
  p2 <- DimPlot(seu, group.by = df_class_col) + ggtitle("DoubletFinder")
  ggsave(file.path(out_dir, paste0(sample_name, "_DoubletFinder_UMAP.png")), p2, width = 6, height = 5)
}
if (!is.null(scrub_col)) {
  p3 <- DimPlot(seu, group.by = scrub_col) + ggtitle("Scrublet")
  ggsave(file.path(out_dir, paste0(sample_name, "_Scrublet_UMAP.png")), p3, width = 6, height = 5)
}

# Save cluster UMAP
ggsave(file.path(out_dir, paste0(sample_name, "_Clusters_UMAP.png")), p1, width = 6, height = 5)

cat("Done. Outputs in:", out_dir, "\n")
```

### DropletQC_per_sample.R

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DropletQC)
  library(DropletUtils)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: DropletQC_per_sample.R <outs_dir> <sample_name> [out_dir] [tiles] [cores] [metadata_csv] [metadata_col]", call.=FALSE)
}

outs_dir <- args[[1]]
sample_name <- args[[2]]
out_dir <- ifelse(length(args) >= 3, args[[3]], file.path(outs_dir, "DropletQC"))
tiles <- ifelse(length(args) >= 4, as.integer(args[[4]]), 100)
cores <- ifelse(length(args) >= 5 && nzchar(args[[5]]), as.integer(args[[5]]), max(1, parallel::detectCores() - 1))
metadata_path <- ifelse(length(args) >= 6 && nzchar(args[[6]]), args[[6]], "")
metadata_col <- ifelse(length(args) >= 7 && nzchar(args[[7]]), args[[7]], "")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Resolve outs directory if a sample folder contains outs/
if (!dir.exists(file.path(outs_dir, "filtered_feature_bc_matrix")) &&
    dir.exists(file.path(outs_dir, "outs", "filtered_feature_bc_matrix"))) {
  outs_dir <- file.path(outs_dir, "outs")
  out_dir <- ifelse(length(args) >= 3, args[[3]], file.path(outs_dir, "DropletQC"))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
}

# Compute nuclear fraction from Cell Ranger outs
nf <- nuclear_fraction_tags(outs = outs_dir, tiles = tiles, cores = cores, verbose = FALSE)

# UMI counts from filtered matrix
filtered_dir <- file.path(outs_dir, "filtered_feature_bc_matrix")
if (!dir.exists(filtered_dir)) stop(paste("filtered_feature_bc_matrix not found:", filtered_dir))

sce <- read10xCounts(filtered_dir)
# Ensure barcodes are present on colnames (some DropletUtils versions store them in colData)
if (is.null(colnames(sce)) && "Barcode" %in% colnames(colData(sce))) {
  colnames(sce) <- as.character(colData(sce)$Barcode)
}
umi <- Matrix::colSums(counts(sce))

# Align barcodes
common <- intersect(rownames(nf), colnames(sce))
if (length(common) == 0) {
  # Try to harmonize barcode suffixes (e.g., add/remove "-1") by matching base barcodes.
  nf_base <- sub("-[0-9]+$", "", rownames(nf))
  sce_base <- sub("-[0-9]+$", "", colnames(sce))
  base_common <- intersect(nf_base, sce_base)
  if (length(base_common) == 0) {
    message("[WARN] No overlapping barcodes between nuclear fraction and filtered matrix. Skipping DropletQC for ", sample_name)
    quit(save = "no", status = 0)
  }
  # Map nuclear-fraction barcodes onto filtered barcodes using base IDs.
  sce_map <- setNames(colnames(sce), sce_base)
  nf_keep <- nf_base %in% base_common
  nf <- nf[nf_keep, , drop = FALSE]
  nf_base <- nf_base[nf_keep]
  mapped <- sce_map[nf_base]
  rownames(nf) <- mapped
  umi <- umi[mapped]
  common <- mapped
} else {
  nf <- nf[common, , drop = FALSE]
  umi <- umi[common]
}

nf_umi <- data.frame(
  nf = as.numeric(nf$nuclear_fraction),
  umi = as.numeric(umi),
  row.names = common
)

# Identify empty droplets
ed <- identify_empty_drops(nf_umi = nf_umi)

# Add cell_type (prefer metadata Leiden/seurat_clusters if available)
cell_type <- rep("all", nrow(ed))
names(cell_type) <- rownames(ed)
if (nzchar(metadata_path) && file.exists(metadata_path)) {
  meta <- read.csv(metadata_path, row.names = 1)
  col <- metadata_col
  if (!nzchar(col)) {
    if ("leiden" %in% colnames(meta)) {
      col <- "leiden"
    } else if ("seurat_clusters" %in% colnames(meta)) {
      col <- "seurat_clusters"
    }
  }
  if (nzchar(col) && col %in% colnames(meta)) {
    meta_vals <- meta[rownames(ed), col, drop = TRUE]
    meta_vals <- as.character(meta_vals)
    cell_type[!is.na(meta_vals)] <- meta_vals[!is.na(meta_vals)]
  } else {
    message("[WARN] DropletQC metadata column not found; using cell_type='all'")
  }
}
ed$cell_type <- cell_type

# Identify damaged cells
res <- identify_damaged_cells(ed, verbose = FALSE, output_plots = FALSE)

# res[[1]] is data frame with cell_status
qc_df <- res[[1]]
qc_df$barcode <- rownames(qc_df)

out_csv <- file.path(out_dir, paste0(sample_name, "_dropletqc_results.csv"))
write.csv(qc_df, out_csv, row.names = FALSE)

cat("DropletQC done for", sample_name, "\n")
```

### Seurat_final_subsetting.R

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# Avoid future globals size errors during SCTransform.
options(future.globals.maxSize = 2 * 1024^3)
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}

safe_cell_cycle_scoring <- function(obj, s_genes, g2m_genes, label = "dataset") {
  if (length(s_genes) == 0 || length(g2m_genes) == 0) {
    message("[WARN] Cell cycle gene lists not found in ", label, "; skipping CellCycleScoring")
    return(obj)
  }
  for (nbin in c(24, 10, 5)) {
    res <- tryCatch(
      CellCycleScoring(obj, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE, nbin = nbin),
      error = function(e) e
    )
    if (!inherits(res, "error")) return(res)
    message("[WARN] CellCycleScoring failed (nbin=", nbin, ") for ", label, ": ", conditionMessage(res))
  }
  message("[WARN] CellCycleScoring skipped for ", label)
  return(obj)
}

parse_bool <- function(x, default = TRUE) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) return(default)
  val <- tolower(as.character(x))
  if (val %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (val %in% c("false", "f", "0", "no", "n")) return(FALSE)
  return(default)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Seurat_final_subsetting.R <root_dir> <qc_out_root> <out_root> [use_soupx] [require_emptydrops] [require_dropletqc] [npcs] [resolution] [seed] [min_cells] [min_features] [do_regress] [do_regress_cell_cycle] [max_percent_mt] [max_percent_ribo] [basic_min_features]", call.=FALSE)
}

root_dir <- args[[1]]
qc_root <- args[[2]]
out_root <- args[[3]]
use_soupx <- parse_bool(ifelse(length(args) >= 4, args[[4]], NA), TRUE)
require_emptydrops <- parse_bool(ifelse(length(args) >= 5, args[[5]], NA), TRUE)
require_dropletqc <- parse_bool(ifelse(length(args) >= 6, args[[6]], NA), TRUE)
npcs <- ifelse(length(args) >= 7 && nzchar(args[[7]]), as.integer(args[[7]]), 30)
resolution <- ifelse(length(args) >= 8 && nzchar(args[[8]]), as.numeric(args[[8]]), 0.5)
seed <- ifelse(length(args) >= 9 && nzchar(args[[9]]), as.integer(args[[9]]), 1)
min_cells <- ifelse(length(args) >= 10 && nzchar(args[[10]]), as.integer(args[[10]]), 3)
min_features <- ifelse(length(args) >= 11 && nzchar(args[[11]]), as.integer(args[[11]]), 200)
do_regress <- parse_bool(ifelse(length(args) >= 12, args[[12]], NA), FALSE)
do_regress_cell_cycle <- parse_bool(ifelse(length(args) >= 13, args[[13]], NA), FALSE)
max_percent_mt <- ifelse(length(args) >= 14 && nzchar(args[[14]]), as.numeric(args[[14]]), Inf)
max_percent_ribo <- ifelse(length(args) >= 15 && nzchar(args[[15]]), as.numeric(args[[15]]), Inf)
basic_min_features <- ifelse(length(args) >= 16 && nzchar(args[[16]]), as.integer(args[[16]]), min_features)

root_dir <- normalizePath(root_dir, mustWork = TRUE)
qc_root <- normalizePath(qc_root, mustWork = TRUE)
out_root <- normalizePath(out_root, mustWork = FALSE)

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

sample_dirs <- list.dirs(root_dir, full.names = TRUE, recursive = FALSE)
if (length(sample_dirs) == 0) {
  stop("No sample directories found in root_dir", call.=FALSE)
}

for (sample_dir in sample_dirs) {
  sample <- basename(sample_dir)
  sample_out <- file.path(out_root, sample)
  out_rds <- if (do_regress) {
    file.path(sample_out, paste0(sample, "_sct_regressed_norm.rds"))
  } else {
    file.path(sample_out, paste0(sample, "_sct_norm.rds"))
  }
  meta_out_path <- file.path(sample_out, "final_metadata.csv")
  if (file.exists(out_rds) && file.exists(meta_out_path)) {
    message("[INFO] Skipping Seurat final (resume): outputs exist for ", sample)
    next
  }

  filtered_dir <- file.path(sample_dir, "filtered_feature_bc_matrix")
  if (!dir.exists(filtered_dir)) {
    alt <- file.path(sample_dir, "outs", "filtered_feature_bc_matrix")
    if (dir.exists(alt)) filtered_dir <- alt
  }
  if (!dir.exists(filtered_dir)) {
    message("[WARN] filtered_feature_bc_matrix not found for ", sample)
    next
  }

  raw_filtered <- Read10X(filtered_dir)
  if (is.list(raw_filtered)) {
    raw_filtered <- raw_filtered[[1]]
  }

  barcodes <- colnames(raw_filtered)

  counts <- raw_filtered
  soupx_path <- file.path(qc_root, "soupx", sample, paste0(sample, "_soupx.mtx"))
  if (use_soupx && file.exists(soupx_path)) {
    soupx_counts <- readMM(soupx_path)
    if (!all(dim(soupx_counts) == dim(raw_filtered))) {
      stop(paste("SoupX matrix dimensions do not match filtered matrix for", sample))
    }
    rownames(soupx_counts) <- rownames(raw_filtered)
    colnames(soupx_counts) <- colnames(raw_filtered)
    counts <- soupx_counts
  } else if (use_soupx) {
    message("[WARN] SoupX matrix missing for ", sample, "; using filtered counts")
  }

  meta <- data.frame(row.names = barcodes)

  meta_path <- file.path(qc_root, "seurat_metadata", sample, "metadata.csv")
  if (file.exists(meta_path)) {
    meta_in <- read.csv(meta_path, row.names = 1)
    missing_barcodes <- setdiff(barcodes, rownames(meta_in))
    if (length(missing_barcodes) > 0) {
      stop(paste("metadata.csv missing", length(missing_barcodes), "barcodes for", sample, "(no pre-filtering before SoupX)."))
    }
    meta <- cbind(meta, meta_in[barcodes, , drop = FALSE])
  }
  # Ensure orig.ident is the sample name (metadata.csv may carry a generic value)
  meta$orig.ident <- sample

  ed_path <- file.path(qc_root, "emptydrops", sample, paste0(sample, "_emptydrops_results.csv"))
  if (file.exists(ed_path)) {
    ed <- read.csv(ed_path, row.names = 1)
    ed_call <- (!is.na(ed$FDR)) & (ed$FDR <= 0.01)
    names(ed_call) <- rownames(ed)
    meta$emptydrops_call <- ed_call[barcodes]
  } else {
    meta$emptydrops_call <- NA
  }

  df_path <- file.path(qc_root, "doubletfinder", sample, paste0(sample, "_metadata_with_DoubletFinder.csv"))
  if (file.exists(df_path)) {
    df <- read.csv(df_path, row.names = 1)
    df_cols <- grep("^DF.classifications", colnames(df), value = TRUE)
    if (length(df_cols) == 1) {
      df_call <- tolower(as.character(df[[df_cols[1]]])) == "doublet"
      names(df_call) <- rownames(df)
      meta$doubletfinder_call <- df_call[barcodes]
    } else {
      meta$doubletfinder_call <- NA
    }
  } else {
    meta$doubletfinder_call <- NA
  }

  dqc_path <- file.path(qc_root, "dropletqc", sample, paste0(sample, "_dropletqc_results.csv"))
  if (file.exists(dqc_path)) {
    dqc <- read.csv(dqc_path)
    if ("barcode" %in% colnames(dqc)) {
      rownames(dqc) <- dqc$barcode
    }
    if ("cell_status" %in% colnames(dqc)) {
      dqc_pass <- tolower(as.character(dqc$cell_status)) == "cell"
      names(dqc_pass) <- rownames(dqc)
      meta$dropletqc_pass <- dqc_pass[barcodes]
      meta$dropletqc_status <- dqc$cell_status[barcodes]
    } else {
      meta$dropletqc_pass <- NA
    }
  } else {
    if (require_dropletqc) {
      stop(paste0(
        "DropletQC results missing for sample ", sample, ": ", dqc_path,
        ". Run DropletQC (requires possorted_genome_bam.bam) or set require_dropletqc=false."
      ), call.=FALSE)
    }
    meta$dropletqc_pass <- NA
  }

  # Compute basic QC metrics on all cells used for downstream filtering
  tmp <- CreateSeuratObject(
    counts = counts,
    meta.data = meta,
    min.cells = 0,
    min.features = 0,
    project = sample
  )
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
  tmp[["percent.rps"]] <- PercentageFeatureSet(tmp, pattern = "^RPS")
  tmp[["percent.rpl"]] <- PercentageFeatureSet(tmp, pattern = "^RPL")
  tmp$percent.ribo <- tmp$percent.rps + tmp$percent.rpl

  meta$nFeature_RNA <- tmp$nFeature_RNA
  meta$nCount_RNA <- tmp$nCount_RNA
  meta$percent.mt <- tmp$percent.mt
  meta$percent.rps <- tmp$percent.rps
  meta$percent.rpl <- tmp$percent.rpl
  meta$percent.ribo <- tmp$percent.ribo

  # Basic QC gate aligned to DoubletFinder prefilter (plus optional ribosomal filter)
  basic_qc_keep <- rep(TRUE, length(barcodes))
  names(basic_qc_keep) <- barcodes
  if (is.finite(basic_min_features)) {
    basic_qc_keep <- basic_qc_keep & (meta$nFeature_RNA >= basic_min_features)
  }
  if (is.finite(max_percent_mt)) {
    basic_qc_keep <- basic_qc_keep & (meta$percent.mt <= max_percent_mt)
  }
  if (is.finite(max_percent_ribo)) {
    basic_qc_keep <- basic_qc_keep & (meta$percent.ribo <= max_percent_ribo)
  }

  # DoubletFinder calls (missing -> FALSE; should be outside basic_qc_keep if thresholds match)
  df_call_full <- rep(FALSE, length(barcodes))
  names(df_call_full) <- barcodes
  if (!all(is.na(meta$doubletfinder_call))) {
    df_call_full[barcodes] <- as.logical(meta$doubletfinder_call)
    df_call_full[is.na(df_call_full)] <- FALSE
  }

  if (!all(is.na(meta$doubletfinder_call))) {
    missing_df <- barcodes[basic_qc_keep & is.na(meta$doubletfinder_call)]
    if (length(missing_df) > 0) {
      message("[WARN] DoubletFinder calls missing for ", length(missing_df), " cells passing basic QC in ", sample, ". Check alignment of thresholds.")
    }
  }

  keep <- basic_qc_keep
  if (require_emptydrops) {
    keep <- keep & !is.na(meta$emptydrops_call) & as.logical(meta$emptydrops_call)
  }
  if (require_dropletqc) {
    keep <- keep & !is.na(meta$dropletqc_pass) & as.logical(meta$dropletqc_pass)
  }
  keep <- keep & !df_call_full

  if (!any(keep, na.rm = TRUE)) {
    message("[WARN] No cells kept after QC filters for ", sample, "; skipping.")
    next
  }

  seu <- CreateSeuratObject(
    counts = counts,
    meta.data = meta,
    min.cells = min_cells,
    min.features = min_features,
    project = sample
  )
  seu <- subset(seu, cells = barcodes[keep])

  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu[["percent.rps"]] <- PercentageFeatureSet(seu, pattern = "^RPS")
  seu[["percent.rpl"]] <- PercentageFeatureSet(seu, pattern = "^RPL")

  s_genes <- cc.genes.updated.2019$s.genes
  g2m_genes <- cc.genes.updated.2019$g2m.genes
  s_genes <- intersect(s_genes, rownames(seu))
  g2m_genes <- intersect(g2m_genes, rownames(seu))

  # Always score cell cycle if gene lists are available; regression remains optional.
  if (length(s_genes) > 0 && length(g2m_genes) > 0) {
    seu <- NormalizeData(seu)
    seu <- safe_cell_cycle_scoring(seu, s_genes, g2m_genes, label = paste0(sample, " (RNA)"))
  } else {
    message("[WARN] Cell cycle gene lists not found in dataset for ", sample, "; skipping scoring")
  }

  do_regress_cell_cycle_sample <- do_regress_cell_cycle

  set.seed(seed)
  vars_to_regress <- c()
  if (do_regress) {
    vars_to_regress <- c(vars_to_regress, "percent.mt", "percent.rps", "percent.rpl")
  }
  if (do_regress_cell_cycle_sample) {
    vars_to_regress <- c(vars_to_regress, "S.Score", "G2M.Score")
  }
  if (length(vars_to_regress) > 0) {
    seu_proc <- SCTransform(seu, vars.to.regress = vars_to_regress, verbose = TRUE)
  } else {
    seu_proc <- SCTransform(seu, verbose = FALSE)
  }

  # Cell-cycle scores already computed on RNA and carried forward; avoid re-scoring after SCT.
  seu_proc <- RunPCA(seu_proc, verbose = FALSE)
  seu_proc <- RunUMAP(seu_proc, dims = 1:npcs, verbose = FALSE)
  seu_proc <- FindNeighbors(seu_proc, dims = 1:npcs, verbose = FALSE)
  seu_proc <- FindClusters(seu_proc, resolution = resolution, verbose = FALSE)

  sample_out <- file.path(out_root, sample)
  dir.create(sample_out, recursive = TRUE, showWarnings = FALSE)

  if (do_regress) {
    saveRDS(seu_proc, file.path(sample_out, paste0(sample, "_sct_regressed_norm.rds")))
  } else {
    saveRDS(seu_proc, file.path(sample_out, paste0(sample, "_sct_norm.rds")))
  }

  meta_out <- cbind(meta, keep_cell = keep)
  write.csv(meta_out, file.path(sample_out, "final_metadata.csv"))
}
```

### Seurat_merge_annotate_integrate.R

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
  library(Azimuth)
})

options(future.globals.maxSize = 2 * 1024^3)
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}

if (!requireNamespace("SeuratData", quietly = TRUE)) {
  stop("SeuratData is required for Azimuth reference installation. Install r-seurat-data.", call.=FALSE)
}

map_genes_to_features <- function(genes, features, label = "genes") {
  genes <- unique(trimws(genes))
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0) return(character())

  present <- genes[genes %in% features]
  missing <- setdiff(genes, features)
  if (length(missing) == 0) return(present)

  if (!requireNamespace("HGNChelper", quietly = TRUE)) {
    stop("HGNChelper is required for gene symbol synonym mapping. Install with conda: r-hgnchelper", call.=FALSE)
  }

  checked <- suppressMessages(HGNChelper::checkGeneSymbols(missing, unmapped.as.na = TRUE))
  suggested <- checked$Suggested.Symbol
  suggested <- suggested[!is.na(suggested)]
  mapped <- suggested[suggested %in% features]

  unique(c(present, mapped))
}

safe_add_module_score <- function(obj, features, name, ctrl, label = "modules") {
  for (nbin in c(24, 10, 5)) {
    res <- tryCatch(
      AddModuleScore(obj, features = features, name = name, ctrl = ctrl, nbin = nbin),
      error = function(e) e
    )
    if (!inherits(res, "error")) return(res)
    message("[WARN] AddModuleScore failed (nbin=", nbin, ") for ", label, ": ", conditionMessage(res))
  }
  stop(paste0("AddModuleScore failed for ", label), call.=FALSE)
}

parse_bool <- function(x, default = TRUE) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) return(default)
  val <- tolower(as.character(x))
  if (val %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (val %in% c("false", "f", "0", "no", "n")) return(FALSE)
  return(default)
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  return(getwd())
}

parse_grouping_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[lines != "" & !startsWith(lines, "#")]
  groups <- list()
  current <- list()
  flush_current <- function() {
    if (length(current) > 0 && !is.null(current$group)) {
      groups[[current$group]] <<- current
    }
  }
  for (ln in lines) {
    if (startsWith(ln, "group:")) {
      flush_current()
      current <- list()
      current$group <- trimws(sub("^group:", "", ln))
    } else if (startsWith(ln, "azimuth_refs:")) {
      val <- trimws(sub("^azimuth_refs:", "", ln))
      current$azimuth_refs <- trimws(unlist(strsplit(val, ",")))
      current$azimuth_refs <- current$azimuth_refs[nzchar(current$azimuth_refs)]
    } else if (startsWith(ln, "azimuth_label_patterns:")) {
      val <- trimws(sub("^azimuth_label_patterns:", "", ln))
      current$azimuth_label_patterns <- trimws(unlist(strsplit(val, ",")))
      current$azimuth_label_patterns <- current$azimuth_label_patterns[nzchar(current$azimuth_label_patterns)]
    } else if (startsWith(ln, "cellmarker_name_patterns:")) {
      val <- trimws(sub("^cellmarker_name_patterns:", "", ln))
      current$cellmarker_name_patterns <- trimws(unlist(strsplit(val, ",")))
      current$cellmarker_name_patterns <- current$cellmarker_name_patterns[nzchar(current$cellmarker_name_patterns)]
    } else if (startsWith(ln, "lineage:")) {
      current$lineage <- trimws(sub("^lineage:", "", ln))
    } else if (startsWith(ln, "notes:")) {
      current$notes <- trimws(sub("^notes:", "", ln))
    }
  }
  flush_current()
  groups
}

normalize_label <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

match_any_pattern <- function(x, patterns, normalize = FALSE) {
  if (length(patterns) == 0 || length(x) == 0) return(rep(FALSE, length(x)))
  if (normalize) {
    x <- normalize_label(x)
    patterns <- normalize_label(patterns)
  }
  res <- rep(FALSE, length(x))
  for (p in patterns) {
    if (!nzchar(p)) next
    res <- res | grepl(p, x, ignore.case = !normalize, fixed = TRUE)
  }
  res
}

ensure_umap <- function(obj, npcs, assay_prefer = c("SCT", "RNA")) {
  assay <- assay_prefer[assay_prefer %in% Assays(obj)][1]
  if (is.na(assay) || !nzchar(assay)) {
    assay <- DefaultAssay(obj)
  }
  DefaultAssay(obj) <- assay
  if (!("pca" %in% Reductions(obj))) {
    if (assay == "RNA") {
      obj <- NormalizeData(obj, verbose = FALSE)
      obj <- FindVariableFeatures(obj, verbose = FALSE)
      obj <- ScaleData(obj, verbose = FALSE)
    }
    feats <- VariableFeatures(obj)
    if (length(feats) == 0) {
      feats <- rownames(obj)
    }
    obj <- RunPCA(obj, npcs = npcs, features = feats, verbose = FALSE)
  }
  if (!("umap" %in% Reductions(obj))) {
    obj <- RunUMAP(obj, dims = 1:npcs, verbose = FALSE)
  }
  obj
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Seurat_merge_annotate_integrate.R <seurat_final_dir> <cellmarker_xlsx> <ewing_signature_csv> <out_dir> [tissues] [min_genes] [ctrl] [seed] [do_integration] [npcs] [resolution]", call.=FALSE)
}

seurat_final_dir <- args[[1]]
cellmarker_xlsx <- args[[2]]
ewing_sig_csv <- args[[3]]
out_dir <- args[[4]]
tissues_arg <- ifelse(length(args) >= 5 && nzchar(args[[5]]), args[[5]], "")
min_genes <- ifelse(length(args) >= 6 && nzchar(args[[6]]), as.integer(args[[6]]), 5)
ctrl <- ifelse(length(args) >= 7 && nzchar(args[[7]]), as.integer(args[[7]]), 20)
seed <- ifelse(length(args) >= 8 && nzchar(args[[8]]), as.integer(args[[8]]), 1)
do_integration <- parse_bool(ifelse(length(args) >= 9, args[[9]], NA), TRUE)
npcs <- ifelse(length(args) >= 10 && nzchar(args[[10]]), as.integer(args[[10]]), 30)
resolution <- ifelse(length(args) >= 11 && nzchar(args[[11]]), as.numeric(args[[11]]), 0.5)
resume <- parse_bool(Sys.getenv("RESUME", "true"), TRUE)

if (!file.exists(cellmarker_xlsx)) stop("CellMarker xlsx not found")
if (!file.exists(ewing_sig_csv)) stop("Ewing sarcoma signature file not found")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Resume: skip if outputs already exist
merged_out <- file.path(out_dir, "merged_annotated.rds")
integrated_out <- file.path(out_dir, "integrated_annotated.rds")
if (resume && file.exists(merged_out) && (!do_integration || file.exists(integrated_out))) {
  message("[INFO] Skipping annotation/integration (resume): outputs already exist in ", out_dir)
  quit(save = "no", status = 0)
}

# Load CellMarker and build modules
x <- read_excel(cellmarker_xlsx)
if (!all(c("species", "cancer_type", "tissue_type", "cell_name", "Symbol") %in% colnames(x))) {
  stop("Required columns not found in CellMarker file")
}
x <- x[x$species == "Human" & x$cancer_type == "Normal", ]
x$tissue_type[is.na(x$tissue_type) | x$tissue_type == ""] <- "Unknown"

if (nzchar(tissues_arg)) {
  if (file.exists(tissues_arg)) {
    tissues <- readLines(tissues_arg, warn = FALSE)
  } else {
    tissues <- unlist(strsplit(tissues_arg, ","))
  }
  tissues <- trimws(tissues)
  tissues <- tissues[nzchar(tissues)]
  if (length(tissues) > 0) {
    tissue_lc <- tolower(x$tissue_type)
    keep_lc <- tolower(tissues)
    x <- x[tissue_lc %in% keep_lc, ]
  }
}

x$Symbol <- trimws(x$Symbol)
x <- x[!is.na(x$Symbol) & x$Symbol != "", ]

cell_names <- sort(unique(x$cell_name))
modules <- list()
for (cn in cell_names) {
  genes <- unique(x$Symbol[x$cell_name == cn])
  genes <- genes[nzchar(genes)]
  if (length(genes) >= min_genes) {
    modules[[cn]] <- genes
  }
}
if (length(modules) == 0) stop("No CellMarker modules passed min_genes after filtering")

# Drop low-information modules (keep CD4/CD8 T cell as requested)
drop_modules <- c(
  "Alveolar cell Type 1",
  "Alveolar type II (ATII) cell",
  "CD14+ monocyte",
  "CD4 recently activated effector memory or effector T cell (CTL)",
  "CD8 recently activated effector memory or effector T cell (CTL)",
  "Ciliated cell",
  "Classical monocyte",
  "Cytotoxic T cell",
  "Dendritic cell",
  "Effector CD4+ T cell",
  "Effector CD8+ T cell",
  "Finally highly effector (TEMRA) memory T cell",
  "Gamma delta(γδ) T cell",
  "Immune cell",
  "Interstitial macrophage",
  "Memory CD4+ T cell",
  "Memory CD8+ T cell",
  "Mesenchymal stromal cell",
  "Monocyte",
  "Myeloid cell",
  "Naive CD4+ T cell",
  "Natural killer cell",
  "Plasma cell",
  "Plasmablast",
  "Progenitor cell",
  "Red blood cell (erythrocyte)",
  "Regulatory CD4+ T cell",
  "T cell",
  "Fibroblast",
  "Reticular fibroblast",
  "Papillary fibroblast",
  "Endothelial progenitor cell",
  "Vascular endothelial cell",
  "Neural progenitor cell",
  "Microglial cell"
)
drop_env <- Sys.getenv("CELL_MARKER_DROP_MODULES", "")
if (nzchar(drop_env)) {
  extra <- trimws(unlist(strsplit(drop_env, ",")))
  extra <- extra[nzchar(extra)]
  drop_modules <- unique(c(drop_modules, extra))
}
modules <- modules[!(names(modules) %in% drop_modules)]

ewing_genes <- tryCatch({
  readLines(ewing_sig_csv, warn = FALSE)
}, error = function(e) {
  as.character(read.csv(ewing_sig_csv, header = FALSE, stringsAsFactors = FALSE)[[1]])
})
ewing_genes <- trimws(ewing_genes)
ewing_genes <- ewing_genes[nzchar(ewing_genes)]

# Load per-sample Seurat objects
sample_dirs <- list.dirs(seurat_final_dir, full.names = TRUE, recursive = FALSE)
if (length(sample_dirs) == 0) stop("No sample directories found in seurat_final_dir")

objs <- list()
for (sd in sample_dirs) {
  sample <- basename(sd)
  rds_reg <- file.path(sd, paste0(sample, "_sct_regressed_norm.rds"))
  rds_noreg <- file.path(sd, paste0(sample, "_sct_norm.rds"))
  rds_in <- if (file.exists(rds_reg)) rds_reg else rds_noreg
  if (!file.exists(rds_in)) {
    message("[WARN] No Seurat RDS found for ", sample)
    next
  }
  obj <- readRDS(rds_in)
  obj <- RenameCells(obj, add.cell.id = sample)
  obj$sample <- sample
  objs[[sample]] <- obj
}

if (length(objs) == 0) stop("No Seurat objects loaded for merge/integration")

# Map gene symbols (including synonyms) to the union of features across samples
features_all <- sort(unique(unlist(lapply(objs, rownames))))
if (mean(grepl("^ENSG", features_all)) > 0.5) {
  message("[WARN] Feature names look like Ensembl IDs; gene symbol synonym mapping may not match.")
}

modules_mapped <- list()
for (cn in names(modules)) {
  genes_mapped <- map_genes_to_features(modules[[cn]], features_all, label = cn)
  if (length(genes_mapped) >= min_genes) {
    modules_mapped[[cn]] <- genes_mapped
  }
}
if (length(modules_mapped) == 0) stop("No CellMarker modules left after synonym mapping")

ewing_genes <- map_genes_to_features(ewing_genes, features_all, label = "Ewing signature")
if (length(ewing_genes) == 0) stop("Ewing sarcoma signature has no genes present after synonym mapping")

module_names_safe <- gsub("[^A-Za-z0-9]+", "_", names(modules_mapped))
module_names_safe <- gsub("^_+|_+$", "", module_names_safe)
module_names_safe <- make.unique(module_names_safe, sep = "_")

module_map <- data.frame(
  module_index = seq_along(modules_mapped),
  module_name = names(modules_mapped),
  module_name_safe = module_names_safe,
  n_genes = vapply(modules_mapped, length, integer(1)),
  stringsAsFactors = FALSE
)

write.csv(module_map, file.path(out_dir, "cellmarker_modules_map.csv"), row.names = FALSE)

# Per-sample annotation (plots generated later on merged UMAP)
score_cols <- paste0("CM_", module_names_safe)
plot_features <- c(score_cols, "EWING_1")
annotated_samples_dir <- file.path(out_dir, "annotated_samples")
dir.create(annotated_samples_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Azimuth references (run per-sample before merging) ----
azimuth_refs_env <- Sys.getenv("AZIMUTH_REFERENCES", "")
if (nzchar(azimuth_refs_env)) {
  azimuth_refs <- strsplit(azimuth_refs_env, ",")[[1]]
  azimuth_refs <- trimws(azimuth_refs)
  azimuth_refs <- azimuth_refs[nzchar(azimuth_refs)]
} else {
  azimuth_refs <- c(
    "pbmcref",       # immune (PBMC)
    "bonemarrowref", # immune (bone marrow)
    "lungref",       # lung
    "adiposeref",    # adipose
    "fetusref",      # fetal development
    "liverref"       # liver
  )
}

# ---- Azimuth homolog table handling (avoid network lookups) ----
azimuth_homologs <- Sys.getenv("AZIMUTH_HOMOLOGS", "")
if (nzchar(azimuth_homologs) && !file.exists(azimuth_homologs)) {
  warning("[Azimuth] AZIMUTH_HOMOLOGS set but file not found: ", azimuth_homologs, ". Skipping gene-name conversion.")
  azimuth_homologs <- ""
}
if (requireNamespace("Azimuth", quietly = TRUE)) {
  orig_convert <- getFromNamespace("ConvertGeneNames", "Azimuth")
  assignInNamespace("ConvertGeneNames", local({
    homologs_path <- azimuth_homologs
    function(object, reference.names, homolog.table = homologs_path) {
      if (nzchar(homologs_path) && file.exists(homologs_path)) {
        return(tryCatch(
          orig_convert(object, reference.names, homolog.table = homologs_path),
          error = function(e) {
            warning("[Azimuth] ConvertGeneNames failed: ", conditionMessage(e))
            object
          }
        ))
      }
      # Skip conversion when offline or no homolog table provided.
      object
    }
  }), ns = "Azimuth")
}

# Install references once if missing
for (ref_name in azimuth_refs) {
  options(timeout = max(getOption("timeout"), 3600))
  installed <- tryCatch(SeuratData::InstalledData(), error = function(e) NULL)
  needs_install <- TRUE
  if (!is.null(installed) && "Dataset" %in% colnames(installed)) {
    needs_install <- !(ref_name %in% installed$Dataset)
  }
  if (needs_install) {
    ok <- tryCatch({
      SeuratData::InstallData(ref_name)
      TRUE
    }, error = function(e) {
      warning("[Azimuth] Failed to install reference ", ref_name, ": ", conditionMessage(e))
      FALSE
    })
    if (!ok) next
  }
  if (!requireNamespace(paste0(ref_name, ".SeuratData"), quietly = TRUE)) {
    warning("[Azimuth] Reference package not available after install: ", ref_name, ". Skipping.")
  }
}

for (sample in names(objs)) {
  annotated_path <- file.path(annotated_samples_dir, paste0(sample, "_annotated.rds"))
  if (resume && file.exists(annotated_path)) {
    obj <- readRDS(annotated_path)
  } else {
    obj <- objs[[sample]]
    if ("SCT" %in% Assays(obj)) {
      DefaultAssay(obj) <- "SCT"
    } else {
      DefaultAssay(obj) <- "RNA"
    }
    genes_in_obj <- rownames(obj)
    module_present_counts <- vapply(modules_mapped, function(g) sum(g %in% genes_in_obj), integer(1))
    modules_present <- modules_mapped[module_present_counts > 0]

    if (length(modules_present) > 0) {
      set.seed(seed)
      obj <- safe_add_module_score(obj, features = modules_present, name = "CM_", ctrl = ctrl, label = paste0("CellMarker:", sample))
      score_cols_raw <- paste0("CM_", seq_along(modules_present))
      present_names_safe <- module_names_safe[match(names(modules_present), names(modules_mapped))]
      score_cols_present <- paste0("CM_", present_names_safe)
      colnames(obj@meta.data)[match(score_cols_raw, colnames(obj@meta.data))] <- score_cols_present
    } else {
      message("[WARN] No CellMarker features present in ", sample, "; CM_ scores set to NA.")
    }

    for (col in score_cols) {
      if (!(col %in% colnames(obj@meta.data))) {
        obj@meta.data[[col]] <- NA_real_
      }
    }
    set.seed(seed)
    obj <- safe_add_module_score(obj, features = list(ewing_sarcoma = ewing_genes), name = "EWING_", ctrl = ctrl, label = paste0("Ewing signature:", sample))
  }

  # Ensure module score columns exist even when resuming
  for (col in score_cols) {
    if (!(col %in% colnames(obj@meta.data))) {
      obj@meta.data[[col]] <- NA_real_
    }
  }

  # Azimuth mapping per-sample (avoid multi-SCT merge issues)
  if (length(azimuth_refs) > 0) {
    if ("SCT" %in% Assays(obj)) {
      DefaultAssay(obj) <- "SCT"
    } else if ("RNA" %in% Assays(obj)) {
      DefaultAssay(obj) <- "RNA"
    }
    for (ref_name in azimuth_refs) {
      if (!requireNamespace(paste0(ref_name, ".SeuratData"), quietly = TRUE)) {
        next
      }
      if (resume) {
        existing <- grepl(paste0("^azimuth_", ref_name, "_"), colnames(obj@meta.data))
        if (any(existing)) next
      }
      message("[Azimuth] Running reference (sample ", sample, "): ", ref_name)
      orig_cols <- colnames(obj@meta.data)
      az_obj <- tryCatch(
        RunAzimuth(obj, reference = ref_name),
        error = function(e) {
          warning("[Azimuth] RunAzimuth failed for ", ref_name, " (sample ", sample, "): ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(az_obj)) next
      new_cols <- setdiff(colnames(az_obj@meta.data), orig_cols)
      if (length(new_cols) > 0) {
        prefix <- paste0("azimuth_", ref_name, "_")
        for (col in new_cols) {
          obj@meta.data[[paste0(prefix, col)]] <- az_obj@meta.data[[col]]
        }
      }
    }
  }

  # Save per-sample annotated object
  saveRDS(obj, annotated_path)

  objs[[sample]] <- obj
}

# Merge non-integrated (after per-sample scoring)
merged <- objs[[1]]
if (length(objs) > 1) {
  merged <- merge(merged, y = objs[2:length(objs)])
}

# ---- Cluster-level group assignment (CellMarker + Azimuth) ----
grouping_file <- Sys.getenv("CELLTYPE_GROUPING_FILE", "")
if (!nzchar(grouping_file)) {
  grouping_file <- file.path(get_script_dir(), "..", "resources", "celltype_grouping.txt")
}
percentile_cutoff <- suppressWarnings(as.numeric(Sys.getenv("CELLTYPE_GROUPING_PERCENTILE", "0.85")))
if (!is.finite(percentile_cutoff) || percentile_cutoff <= 0 || percentile_cutoff >= 1) {
  percentile_cutoff <- 0.85
}
min_mod_cells <- suppressWarnings(as.integer(Sys.getenv("CELLTYPE_GROUPING_MIN_MOD_CELLS", "10")))
min_az_cells <- suppressWarnings(as.integer(Sys.getenv("CELLTYPE_GROUPING_MIN_AZ_CELLS", "10")))
min_mod_frac <- suppressWarnings(as.numeric(Sys.getenv("CELLTYPE_GROUPING_MIN_MOD_FRAC", "0.10")))
min_az_frac <- suppressWarnings(as.numeric(Sys.getenv("CELLTYPE_GROUPING_MIN_AZ_FRAC", "0.10")))
module_only_min_cells <- suppressWarnings(as.integer(Sys.getenv("CELLTYPE_GROUPING_MODULE_ONLY_MIN_CELLS", "20")))
module_only_min_frac <- suppressWarnings(as.numeric(Sys.getenv("CELLTYPE_GROUPING_MODULE_ONLY_MIN_FRAC", "0.25")))
module_only_min_pct <- suppressWarnings(as.numeric(Sys.getenv("CELLTYPE_GROUPING_MODULE_ONLY_MIN_PCT", "90")))
az_nonimmune_pct <- suppressWarnings(as.numeric(Sys.getenv("CELLTYPE_GROUPING_AZ_PCT_NONIMM", "85")))
az_lineage_delta <- suppressWarnings(as.numeric(Sys.getenv("CELLTYPE_GROUPING_AZ_LINEAGE_DELTA", "0.10")))
az_lineage_min <- suppressWarnings(as.numeric(Sys.getenv("CELLTYPE_GROUPING_AZ_LINEAGE_MIN", "0.20")))
ewing_percentile <- suppressWarnings(as.numeric(Sys.getenv("CELLTYPE_GROUPING_EWING_PERCENTILE", "0.75")))
ewing_min_cells <- suppressWarnings(as.integer(Sys.getenv("CELLTYPE_GROUPING_EWING_MIN_CELLS", "20")))
ewing_min_frac <- suppressWarnings(as.numeric(Sys.getenv("CELLTYPE_GROUPING_EWING_MIN_FRAC", "0.20")))
if (!is.finite(min_mod_cells) || min_mod_cells < 1) min_mod_cells <- 10
if (!is.finite(min_az_cells) || min_az_cells < 1) min_az_cells <- 10
if (!is.finite(min_mod_frac) || min_mod_frac < 0 || min_mod_frac > 1) min_mod_frac <- 0.10
if (!is.finite(min_az_frac) || min_az_frac < 0 || min_az_frac > 1) min_az_frac <- 0.10
if (!is.finite(module_only_min_cells) || module_only_min_cells < 1) module_only_min_cells <- 20
if (!is.finite(module_only_min_frac) || module_only_min_frac < 0 || module_only_min_frac > 1) module_only_min_frac <- 0.25
if (!is.finite(module_only_min_pct) || module_only_min_pct < 0 || module_only_min_pct > 100) module_only_min_pct <- NA_real_
if (!is.finite(az_nonimmune_pct) || az_nonimmune_pct < 0 || az_nonimmune_pct > 100) az_nonimmune_pct <- 85
if (!is.finite(az_lineage_delta) || az_lineage_delta < 0) az_lineage_delta <- 0.10
if (!is.finite(az_lineage_min) || az_lineage_min < 0) az_lineage_min <- 0.20
if (!is.finite(ewing_percentile) || ewing_percentile <= 0 || ewing_percentile >= 1) ewing_percentile <- 0.75
if (!is.finite(ewing_min_cells) || ewing_min_cells < 1) ewing_min_cells <- 20
if (!is.finite(ewing_min_frac) || ewing_min_frac < 0 || ewing_min_frac > 1) ewing_min_frac <- 0.20

if (file.exists(grouping_file)) {
  groups <- parse_grouping_file(grouping_file)
  if (length(groups) > 0) {
    group_names_for_plots <- names(groups)
    for (gname in names(groups)) {
      if (is.null(groups[[gname]]$lineage) || !nzchar(groups[[gname]]$lineage)) {
        groups[[gname]]$lineage <- gname
      }
    }
    group_lineage <- vapply(groups, function(g) ifelse(nzchar(g$lineage), g$lineage, g$group), character(1))
    lineage_groups <- split(names(groups), group_lineage)
    immune_lineages <- unique(group_lineage[grepl("^immune_", names(groups))])
    meta <- merged@meta.data

    cluster_raw <- NULL
    if ("seurat_clusters" %in% colnames(meta)) {
      cluster_raw <- meta$seurat_clusters
    } else if (!is.null(Idents(merged))) {
      cluster_raw <- Idents(merged)
    } else {
      cluster_raw <- rep("all", nrow(meta))
      message("[WARN] No cluster column found; assigning all cells to a single cluster.")
    }
    cluster_raw <- as.character(cluster_raw)
    if ("sample" %in% colnames(meta)) {
      cluster_id <- paste(meta$sample, cluster_raw, sep = "_")
    } else {
      cluster_id <- cluster_raw
    }
    meta$cluster_id_for_assignment <- cluster_id

    # Map CellMarker modules to groups using module_map
    module_cols_by_group <- list()
    if (exists("module_map")) {
      for (gname in names(groups)) {
        g <- groups[[gname]]
        patterns <- g$cellmarker_name_patterns
        if (length(patterns) == 0) next
        matches <- match_any_pattern(module_map$module_name, patterns)
        if (any(matches)) {
          cols <- module_map$module_name_safe[matches]
          module_cols_by_group[[gname]] <- cols
        }
      }
    }

    # Azimuth label columns to use (stable levels)
    azimuth_levels <- list(
      pbmcref = "celltype.l2",
      bonemarrowref = "celltype.l2",
      lungref = "ann_level_3",
      adiposeref = "celltype.l1",
      fetusref = "annotation.l1",
      liverref = "celltype.l2"
    )

    # Precompute per-cell group azimuth scores (max across refs)
    group_azimuth_score <- list()
    for (gname in names(groups)) {
      g <- groups[[gname]]
      az_scores <- rep(NA_real_, nrow(meta))
      if (!is.null(g$azimuth_refs) && length(g$azimuth_refs) > 0 && length(g$azimuth_label_patterns) > 0) {
        for (ref in g$azimuth_refs) {
          level <- azimuth_levels[[ref]]
          if (is.null(level)) next
          label_col <- paste0("azimuth_", ref, "_predicted.", level)
          score_col <- paste0("azimuth_", ref, "_predicted.", level, ".score")
          if (!(label_col %in% colnames(meta))) next
          if (!(score_col %in% colnames(meta))) {
            map_col <- paste0("azimuth_", ref, "_mapping.score")
            if (map_col %in% colnames(meta)) {
              score_col <- map_col
            } else {
              next
            }
          }
          labels <- as.character(meta[[label_col]])
          scores <- suppressWarnings(as.numeric(meta[[score_col]]))
          match <- match_any_pattern(labels, g$azimuth_label_patterns, normalize = TRUE)
          scores[!match] <- NA_real_
          if (all(is.na(scores))) next
          if (all(is.na(az_scores))) {
            az_scores <- scores
          } else {
            az_scores <- pmax(az_scores, scores, na.rm = TRUE)
          }
        }
      }
      group_azimuth_score[[gname]] <- az_scores
    }

    # Cluster-level medians
    cluster_ids <- sort(unique(cluster_id))
    cluster_n <- vapply(cluster_ids, function(cid) sum(cluster_id == cid), integer(1))
    cluster_summary <- data.frame(
      cluster_id = cluster_ids,
      n_cells = cluster_n,
      stringsAsFactors = FALSE
    )

    # EWING signature medians and percentiles
    ewing_med <- tapply(meta$EWING_1, cluster_id, median, na.rm = TRUE)
    ewing_med <- ewing_med[cluster_ids]
    ewing_pct <- rep(NA_real_, length(cluster_ids))
    if (sum(!is.na(ewing_med)) >= 1) {
      ewing_pct <- ecdf(ewing_med[!is.na(ewing_med)])(ewing_med) * 100
    }
    ewing_thresh <- NA_real_
    ewing_frac <- rep(NA_real_, length(cluster_ids))
    ewing_n <- tapply(!is.na(meta$EWING_1), cluster_id, sum)
    ewing_n <- ewing_n[cluster_ids]
    if (sum(!is.na(meta$EWING_1)) >= 1) {
      ewing_thresh <- suppressWarnings(as.numeric(quantile(meta$EWING_1, probs = ewing_percentile, na.rm = TRUE, names = FALSE)))
      if (is.finite(ewing_thresh)) {
        ewing_frac <- tapply(meta$EWING_1 >= ewing_thresh, cluster_id, mean, na.rm = TRUE)
        ewing_frac <- ewing_frac[cluster_ids]
      }
    }
    cluster_summary$EWING_1_median <- ewing_med
    cluster_summary$EWING_1_percentile <- ewing_pct
    cluster_summary$EWING_1_frac <- ewing_frac
    cluster_summary$EWING_1_n <- ewing_n
    cluster_summary$EWING_1_threshold <- ewing_thresh

    # Module-level thresholds (per signature)
    module_cols <- character()
    if (exists("module_map")) {
      module_cols <- paste0("CM_", module_map$module_name_safe)
      module_cols <- module_cols[module_cols %in% colnames(meta)]
    }
    module_thresh <- list()
    module_frac <- list()
    module_n <- list()
    module_med <- list()
    module_pct <- list()
    module_strong <- list()
    module_strict <- list()
    for (col in module_cols) {
      name_safe <- sub("^CM_", "", col)
      scores <- meta[[col]]
      mod_thresh <- NA_real_
      if (sum(!is.na(scores)) >= 1) {
        mod_thresh <- suppressWarnings(as.numeric(quantile(scores, probs = percentile_cutoff, na.rm = TRUE, names = FALSE)))
      }
      mod_n <- tapply(!is.na(scores), cluster_id, sum)
      mod_n <- mod_n[cluster_ids]
      mod_frac <- rep(NA_real_, length(cluster_ids))
      if (is.finite(mod_thresh)) {
        mod_frac <- tapply(scores >= mod_thresh, cluster_id, mean, na.rm = TRUE)
        mod_frac <- mod_frac[cluster_ids]
      }
      mod_med <- tapply(scores, cluster_id, median, na.rm = TRUE)
      mod_med <- mod_med[cluster_ids]
      mod_pct <- rep(NA_real_, length(cluster_ids))
      if (sum(!is.na(mod_med)) >= 1) {
        mod_pct <- ecdf(mod_med[!is.na(mod_med)])(mod_med) * 100
      }
      module_thresh[[name_safe]] <- mod_thresh
      module_n[[name_safe]] <- mod_n
      module_frac[[name_safe]] <- mod_frac
      module_med[[name_safe]] <- mod_med
      module_pct[[name_safe]] <- mod_pct
      module_strong[[name_safe]] <- !is.na(mod_n) & !is.na(mod_frac) & mod_n >= min_mod_cells & mod_frac >= min_mod_frac
      module_strict[[name_safe]] <- !is.na(mod_n) & !is.na(mod_frac) & mod_n >= module_only_min_cells & mod_frac >= module_only_min_frac
      if (!is.na(module_only_min_pct)) {
        module_strict[[name_safe]] <- module_strict[[name_safe]] & !is.na(mod_pct) & mod_pct >= module_only_min_pct
      }
    }

    # Group azimuth medians/percentiles + thresholds
    group_names <- names(groups)
    group_az_medians <- list()
    group_az_pct <- list()
    group_az_thresh <- list()
    group_az_frac <- list()
    group_az_n <- list()
    group_az_strong <- list()
    group_az_support <- list()
    group_module_strong_n <- list()
    group_module_strong_any <- list()
    group_module_strict_n <- list()
    group_module_strict_any <- list()

    for (gname in group_names) {
      az_score <- group_azimuth_score[[gname]]

      az_med <- tapply(az_score, cluster_id, median, na.rm = TRUE)
      az_med <- az_med[cluster_ids]

      az_pct <- rep(NA_real_, length(cluster_ids))
      if (sum(!is.na(az_med)) >= 1) {
        az_pct <- ecdf(az_med[!is.na(az_med)])(az_med) * 100
      }

      group_az_medians[[gname]] <- az_med
      group_az_pct[[gname]] <- az_pct

      # Distribution-aware thresholds (global across cells)
      az_thresh <- NA_real_
      if (sum(!is.na(az_score)) >= 1) {
        az_thresh <- suppressWarnings(as.numeric(quantile(az_score, probs = percentile_cutoff, na.rm = TRUE, names = FALSE)))
      }
      group_az_thresh[[gname]] <- az_thresh

      az_n <- tapply(!is.na(az_score), cluster_id, sum)
      az_n <- az_n[cluster_ids]

      az_frac <- rep(NA_real_, length(cluster_ids))
      if (is.finite(az_thresh)) {
        az_frac <- tapply(az_score >= az_thresh, cluster_id, mean, na.rm = TRUE)
        az_frac <- az_frac[cluster_ids]
      }
      group_az_n[[gname]] <- az_n
      group_az_frac[[gname]] <- az_frac

      group_az_strong[[gname]] <- !is.na(az_n) & !is.na(az_frac) & az_n >= min_az_cells & az_frac >= min_az_frac
      az_support <- group_az_strong[[gname]]
      if (!(group_lineage[[gname]] %in% immune_lineages)) {
        az_support <- az_support | (!is.na(az_pct) & az_pct >= az_nonimmune_pct & !is.na(az_n) & az_n >= min_az_cells)
      }
      group_az_support[[gname]] <- az_support

      mods <- module_cols_by_group[[gname]]
      if (is.null(mods) || length(mods) == 0) {
        strong_n <- rep(0L, length(cluster_ids))
        strict_n <- rep(0L, length(cluster_ids))
      } else {
        mods <- mods[mods %in% names(module_strong)]
        if (length(mods) == 0) {
          strong_n <- rep(0L, length(cluster_ids))
          strict_n <- rep(0L, length(cluster_ids))
        } else {
          strong_mat <- vapply(mods, function(m) module_strong[[m]], logical(length(cluster_ids)))
          strict_mat <- vapply(mods, function(m) module_strict[[m]], logical(length(cluster_ids)))
          if (is.null(dim(strong_mat))) {
            strong_mat <- matrix(strong_mat, ncol = 1)
          }
          if (is.null(dim(strict_mat))) {
            strict_mat <- matrix(strict_mat, ncol = 1)
          }
          strong_n <- rowSums(strong_mat)
          strict_n <- rowSums(strict_mat)
        }
      }
      group_module_strong_n[[gname]] <- strong_n
      group_module_strong_any[[gname]] <- strong_n > 0
      group_module_strict_n[[gname]] <- strict_n
      group_module_strict_any[[gname]] <- strict_n > 0
    }

    # Assign groups by cluster (signature-level evidence + Azimuth compatibility)
    cluster_cells <- split(seq_len(nrow(meta)), cluster_id)
    get_best_azimuth_label <- function(idx, gname) {
      g <- groups[[gname]]
      if (is.null(g) || length(idx) == 0) return(NA_character_)
      az_thresh <- group_az_thresh[[gname]]
      labels_all <- character()
      scores_all <- numeric()
      if (!is.null(g$azimuth_refs) && length(g$azimuth_refs) > 0 && length(g$azimuth_label_patterns) > 0) {
        for (ref in g$azimuth_refs) {
          level <- azimuth_levels[[ref]]
          if (is.null(level)) next
          label_col <- paste0("azimuth_", ref, "_predicted.", level)
          score_col <- paste0("azimuth_", ref, "_predicted.", level, ".score")
          if (!(label_col %in% colnames(meta))) next
          if (!(score_col %in% colnames(meta))) {
            map_col <- paste0("azimuth_", ref, "_mapping.score")
            if (map_col %in% colnames(meta)) {
              score_col <- map_col
            } else {
              next
            }
          }
          labels <- as.character(meta[[label_col]][idx])
          scores <- suppressWarnings(as.numeric(meta[[score_col]][idx]))
          match <- match_any_pattern(labels, g$azimuth_label_patterns, normalize = TRUE)
          labels <- labels[match]
          scores <- scores[match]
          if (length(labels) == 0) next
          if (is.finite(az_thresh)) {
            keep <- !is.na(scores) & scores >= az_thresh
            labels <- labels[keep]
            scores <- scores[keep]
          }
          if (length(labels) == 0) next
          labels_all <- c(labels_all, labels)
          scores_all <- c(scores_all, scores)
        }
      }
      if (length(labels_all) == 0) return(NA_character_)
      labs <- unique(labels_all)
      counts <- vapply(labs, function(l) sum(labels_all == l), integer(1))
      med_scores <- vapply(labs, function(l) median(scores_all[labels_all == l], na.rm = TRUE), numeric(1))
      best <- labs[order(-counts, -med_scores)][1]
      best
    }
    get_lineage_az_score <- function(i, lineage) {
      groups_lin <- lineage_groups[[lineage]]
      if (length(groups_lin) == 0) return(NA_real_)
      scores <- vapply(groups_lin, function(g) {
        s <- group_az_medians[[g]][i]
        if (!is.finite(s)) {
          p <- group_az_pct[[g]][i]
          if (is.finite(p)) s <- p / 100
        }
        s
      }, numeric(1))
      scores <- scores[is.finite(scores)]
      if (length(scores) == 0) return(NA_real_)
      max(scores, na.rm = TRUE)
    }

    assigned_group <- rep("Unknown", length(cluster_ids))
    assigned_lineage <- rep(NA_character_, length(cluster_ids))
    assigned_subtype <- rep(NA_character_, length(cluster_ids))
    for (i in seq_along(cluster_ids)) {
      lineage_mods <- list()
      for (lin in names(lineage_groups)) {
        groups_lin <- lineage_groups[[lin]]
        mods <- unique(unlist(lapply(groups_lin, function(g) module_cols_by_group[[g]])))
        mods <- mods[mods %in% names(module_strong)]
        strong_mods <- mods[vapply(mods, function(m) isTRUE(module_strong[[m]][i]), logical(1))]
        has_mod <- length(strong_mods) > 0
        has_az <- any(vapply(groups_lin, function(g) isTRUE(group_az_support[[g]][i]), logical(1)))
        if (lin %in% immune_lineages) {
          if (has_mod && has_az) lineage_mods[[lin]] <- strong_mods
        } else {
          if (has_mod || has_az) lineage_mods[[lin]] <- strong_mods
        }
      }
      lineage_set <- names(lineage_mods)

      az_lineage_set <- character()
      for (gname in group_names) {
        if (isTRUE(group_az_support[[gname]][i])) {
          az_lineage_set <- union(az_lineage_set, group_lineage[[gname]])
        }
      }
      if (length(az_lineage_set) > 1) {
        lin_scores <- vapply(az_lineage_set, function(lin) get_lineage_az_score(i, lin), numeric(1))
        ord <- order(lin_scores, decreasing = TRUE)
        top_score <- lin_scores[ord[1]]
        second_score <- if (length(lin_scores) >= 2) lin_scores[ord[2]] else -Inf
        if (is.finite(top_score) && top_score >= az_lineage_min && (top_score - second_score) >= az_lineage_delta) {
          az_lineage_set <- az_lineage_set[ord[1]]
        }
      }

      ewing_strong <- (!is.na(ewing_frac[i]) && !is.na(ewing_n[i]) &&
                       ewing_n[i] >= ewing_min_cells && ewing_frac[i] >= ewing_min_frac) ||
                      (!is.na(ewing_pct[i]) && ewing_pct[i] >= percentile_cutoff * 100)

      if (ewing_strong) {
        assigned_group[i] <- "Tumor"
        assigned_lineage[i] <- "Tumor"
        next
      }

      if (length(lineage_set) > 1) {
        assigned_group[i] <- "Mixed/Doublet"
        assigned_lineage[i] <- "Mixed/Doublet"
        next
      }

      if (length(lineage_set) == 0) {
        if (length(az_lineage_set) == 1) {
          candidate_lineage <- az_lineage_set[[1]]
          if (!(candidate_lineage %in% immune_lineages)) {
            assigned_group[i] <- candidate_lineage
            assigned_lineage[i] <- candidate_lineage
          }
        }
        next
      }

      candidate_lineage <- lineage_set[[1]]
      assigned_lineage[i] <- candidate_lineage
      conflict <- length(setdiff(az_lineage_set, candidate_lineage)) > 0
      if (conflict) {
        assigned_group[i] <- "Unknown"
        next
      }

      if (candidate_lineage %in% immune_lineages) {
        best_label <- NA_character_
        best_score <- -Inf
        for (gname in lineage_groups[[candidate_lineage]]) {
          if (isTRUE(group_module_strong_any[[gname]][i]) && isTRUE(group_az_strong[[gname]][i])) {
            label <- get_best_azimuth_label(cluster_cells[[cluster_ids[i]]], gname)
            if (!is.na(label) && nzchar(label)) {
              score <- group_az_medians[[gname]][i]
              if (!is.finite(score)) score <- group_az_frac[[gname]][i]
              if (!is.finite(score)) score <- -Inf
              if (score > best_score) {
                best_score <- score
                best_label <- label
              }
            }
          }
        }
        if (!is.na(best_label)) {
          assigned_group[i] <- best_label
          assigned_subtype[i] <- best_label
          next
        }
      }

      # Module-only assignment requires strict thresholds when Azimuth is weak/absent
      strict_ok <- FALSE
      for (gname in lineage_groups[[candidate_lineage]]) {
        if (isTRUE(group_module_strict_any[[gname]][i])) {
          strict_ok <- TRUE
          break
        }
      }
      if (!strict_ok && !(candidate_lineage %in% immune_lineages) && (candidate_lineage %in% az_lineage_set)) {
        assigned_group[i] <- candidate_lineage
      } else if (strict_ok) {
        assigned_group[i] <- candidate_lineage
      } else {
        assigned_group[i] <- "Unknown"
      }
    }

    cluster_summary$assigned_group <- assigned_group
    cluster_summary$assigned_lineage <- assigned_lineage
    cluster_summary$assigned_subtype <- assigned_subtype

    # Add per-group summary to cluster table
    for (gname in group_names) {
      cluster_summary[[paste0(gname, "_azimuth_median")]] <- group_az_medians[[gname]]
      cluster_summary[[paste0(gname, "_azimuth_pct")]] <- group_az_pct[[gname]]
      cluster_summary[[paste0(gname, "_azimuth_frac")]] <- group_az_frac[[gname]]
      cluster_summary[[paste0(gname, "_azimuth_n")]] <- group_az_n[[gname]]
      cluster_summary[[paste0(gname, "_azimuth_strong")]] <- group_az_strong[[gname]]
      cluster_summary[[paste0(gname, "_module_strong_n")]] <- group_module_strong_n[[gname]]
      cluster_summary[[paste0(gname, "_module_strong_any")]] <- group_module_strong_any[[gname]]
      cluster_summary[[paste0(gname, "_module_strict_n")]] <- group_module_strict_n[[gname]]
      cluster_summary[[paste0(gname, "_module_strict_any")]] <- group_module_strict_any[[gname]]
    }

    # Propagate cluster assignment back to cells
    group_map <- setNames(assigned_group, cluster_ids)
    meta$celltype_group <- group_map[cluster_id]
    merged@meta.data <- meta

    write.csv(cluster_summary, file.path(out_dir, "celltype_assignment_by_cluster.csv"), row.names = FALSE)
  } else {
    message("[WARN] No groups found in grouping file: ", grouping_file)
    group_names_for_plots <- character()
  }
} else {
  message("[WARN] Grouping file not found: ", grouping_file)
  group_names_for_plots <- character()
}

merged <- ensure_umap(merged, npcs = npcs)
plots_dir_umap <- file.path(out_dir, "umap_pre")
dir.create(plots_dir_umap, recursive = TRUE, showWarnings = FALSE)
sample_col <- if ("sample" %in% colnames(merged@meta.data)) "sample" else if ("orig.ident" %in% colnames(merged@meta.data)) "orig.ident" else NULL
if (!is.null(sample_col)) {
  ggplot2::ggsave(
    file.path(plots_dir_umap, "umap_pre_by_sample.png"),
    DimPlot(merged, group.by = sample_col) + ggplot2::ggtitle("Merged (pre-integration) by sample"),
    width = 6, height = 5
  )
}
if ("celltype_group" %in% colnames(merged@meta.data)) {
  ggplot2::ggsave(
    file.path(plots_dir_umap, "umap_pre_by_celltype.png"),
    DimPlot(merged, group.by = "celltype_group", label = TRUE) + ggplot2::ggtitle("Merged (pre-integration) by cell type"),
    width = 6, height = 5
  )
}

# Per-sample feature plots on merged UMAP (consistent embedding)
if (!is.null(sample_col)) {
  plots_dir_samples <- file.path(out_dir, "featureplots_samples")
  dir.create(plots_dir_samples, recursive = TRUE, showWarnings = FALSE)
  for (sample in unique(merged@meta.data[[sample_col]])) {
    sample_plot_dir <- file.path(plots_dir_samples, sample)
    dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)
    keep_cells <- rownames(merged@meta.data)[merged@meta.data[[sample_col]] == sample]
    obj_s <- subset(merged, cells = keep_cells)
    for (feat in plot_features) {
      if (feat %in% colnames(obj_s@meta.data)) {
        vals <- obj_s@meta.data[[feat]]
        vals <- vals[!is.na(vals)]
        if (length(vals) == 0 || length(unique(vals)) < 2) {
          next
        }
        ggplot2::ggsave(
          file.path(sample_plot_dir, paste0(feat, ".png")),
          FeaturePlot(obj_s, features = feat, reduction = "umap") + ggplot2::ggtitle(paste0(feat, " (", sample, ")")),
          width = 6, height = 5
        )
      }
    }
  }
}

# Per-sample celltype plots (post-annotation)
if ("celltype_group" %in% colnames(merged@meta.data) && !is.null(sample_col)) {
  plots_dir_celltypes <- file.path(out_dir, "featureplots_celltypes_samples")
  dir.create(plots_dir_celltypes, recursive = TRUE, showWarnings = FALSE)
  safe_filename <- function(x) {
    x <- gsub("[/\\\\]+", "_", x)
    x <- gsub("[[:space:]]+", "_", x)
    x <- gsub("[^A-Za-z0-9_.-]", "_", x)
    x
  }
  plot_groups <- sort(unique(merged@meta.data$celltype_group))
  plot_groups <- plot_groups[nzchar(plot_groups) & plot_groups != "Unknown"]
  for (sample in unique(merged@meta.data[[sample_col]])) {
    sample_dir <- file.path(plots_dir_celltypes, sample)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    keep_cells <- rownames(merged@meta.data)[merged@meta.data[[sample_col]] == sample]
    obj_s <- subset(merged, cells = keep_cells)
    for (g in plot_groups) {
      cells_g <- rownames(obj_s@meta.data)[obj_s@meta.data$celltype_group == g]
      if (length(cells_g) == 0) next
      p <- DimPlot(obj_s, cells.highlight = cells_g, cols = c("grey85", "#D55E00")) +
        ggplot2::ggtitle(paste0(sample, " - ", g))
      ggplot2::ggsave(
        file.path(sample_dir, paste0("celltype_", safe_filename(g), ".png")),
        p,
        width = 6, height = 5
      )
    }
  }
}

saveRDS(merged, file.path(out_dir, "merged_annotated.rds"))
azimuth_cols <- grep("^azimuth_", colnames(merged@meta.data), value = TRUE)
merged_csv_cols <- unique(c(score_cols, "EWING_1", azimuth_cols, "celltype_group"))
write.csv(
  cbind(
    data.frame(barcode = rownames(merged@meta.data)),
    merged@meta.data[, merged_csv_cols, drop = FALSE]
  ),
  file.path(out_dir, "merged_annotations.csv"),
  row.names = FALSE
)

if (do_integration) {
  if (!all(vapply(objs, function(o) "SCT" %in% Assays(o), logical(1)))) {
    stop("SCT assay missing in one or more objects; integration requires SCT")
  }

  features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)
  objs <- PrepSCTIntegration(object.list = objs, anchor.features = features)
  anchors <- FindIntegrationAnchors(object.list = objs, normalization.method = "SCT", anchor.features = features)
  integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

  DefaultAssay(integrated) <- "integrated"
  integrated <- RunPCA(integrated, npcs = npcs, verbose = FALSE)
  integrated <- RunUMAP(integrated, dims = 1:npcs, verbose = FALSE)
  integrated <- FindNeighbors(integrated, dims = 1:npcs, verbose = FALSE)
  integrated <- FindClusters(integrated, resolution = resolution, verbose = FALSE)

  # Transfer annotation fields from merged (by cell name)
  azimuth_cols <- grep("^azimuth_", colnames(merged@meta.data), value = TRUE)
  transfer_cols <- unique(c(score_cols, "EWING_1", azimuth_cols, "celltype_group"))
  for (col in transfer_cols) {
    if (col %in% colnames(merged@meta.data)) {
      integrated@meta.data[[col]] <- merged@meta.data[rownames(integrated@meta.data), col]
    }
  }
  saveRDS(integrated, file.path(out_dir, "integrated_annotated.rds"))

  ggplot2::ggsave(
    file.path(out_dir, "integrated_umap_clusters.png"),
    DimPlot(integrated, group.by = "seurat_clusters", label = TRUE) + ggplot2::ggtitle("Integrated"),
    width = 6, height = 5
  )
  ggplot2::ggsave(
    file.path(out_dir, "integrated_umap_by_sample.png"),
    DimPlot(integrated, group.by = sample_col) + ggplot2::ggtitle("Integrated by sample"),
    width = 6, height = 5
  )
  if ("celltype_group" %in% colnames(integrated@meta.data)) {
    ggplot2::ggsave(
      file.path(out_dir, "integrated_umap_by_celltype.png"),
      DimPlot(integrated, group.by = "celltype_group", label = TRUE) + ggplot2::ggtitle("Integrated by cell type"),
      width = 6, height = 5
    )
  }
  plots_dir_integrated <- file.path(out_dir, "featureplots_integrated")
  dir.create(plots_dir_integrated, recursive = TRUE, showWarnings = FALSE)
  for (feat in plot_features) {
    if (feat %in% colnames(integrated@meta.data)) {
      vals <- integrated@meta.data[[feat]]
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0 || length(unique(vals)) < 2) {
        next
      }
      ggplot2::ggsave(
        file.path(plots_dir_integrated, paste0(feat, ".png")),
        FeaturePlot(integrated, features = feat) + ggplot2::ggtitle(paste0(feat, " (integrated)")),
        width = 6, height = 5
      )
    }
  }
}
```

### Azimuth_label_catalog.R

```r
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(utils)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Azimuth_label_catalog.R <out_file> [refs_csv]", call. = FALSE)
}

out_file <- args[[1]]
refs_arg <- ifelse(length(args) >= 2, args[[2]], "")

get_refs <- function(refs_arg) {
  refs <- character()
  if (nzchar(refs_arg)) {
    refs <- strsplit(refs_arg, ",", fixed = TRUE)[[1]]
  } else if (nzchar(Sys.getenv("AZIMUTH_REFERENCES", ""))) {
    refs <- strsplit(Sys.getenv("AZIMUTH_REFERENCES"), ",", fixed = TRUE)[[1]]
  } else {
    if (requireNamespace("SeuratData", quietly = TRUE)) {
      inst <- tryCatch(SeuratData::InstalledData(), error = function(e) NULL)
      if (!is.null(inst) && "Dataset" %in% colnames(inst)) {
        refs <- inst$Dataset
      }
    }
  }
  refs <- trimws(refs)
  refs[nzchar(refs)]
}

load_reference_object <- function(ref) {
  pkg <- paste0(ref, ".SeuratData")
  if (!requireNamespace(pkg, quietly = TRUE)) return(NULL)

  items <- character()
  items <- tryCatch({
    utils::data(package = pkg)$results[, "Item"]
  }, error = function(e) character())
  if (length(items) == 0) return(NULL)

  obj_name <- if (ref %in% items) ref else items[[1]]
  env <- new.env(parent = emptyenv())
  ok <- tryCatch({
    utils::data(list = obj_name, package = pkg, envir = env)
    TRUE
  }, error = function(e) FALSE)
  if (!ok || !exists(obj_name, envir = env)) return(NULL)

  get(obj_name, envir = env)
}

get_meta_label_cols <- function(meta) {
  is_text <- vapply(meta, function(x) is.character(x) || is.factor(x), logical(1))
  name_match <- grepl("celltype|annotation|ann_level", names(meta), ignore.case = TRUE)
  names(meta)[is_text & name_match]
}

refs <- get_refs(refs_arg)
if (length(refs) == 0) {
  stop("No references found. Provide refs_csv or set AZIMUTH_REFERENCES.", call. = FALSE)
}

con <- file(out_file, open = "w")
on.exit(close(con), add = TRUE)

for (ref in refs) {
  writeLines(paste0("===  ", ref, "  ==="), con)
  obj <- load_reference_object(ref)
  if (is.null(obj)) {
    writeLines("Missing reference path", con)
    writeLines("", con)
    next
  }
  meta <- obj@meta.data
  cols <- get_meta_label_cols(meta)
  if (length(cols) == 0) {
    writeLines("Meta cols:  (none)", con)
    writeLines("", con)
    next
  }
  writeLines(paste0("Meta cols:  ", paste(cols, collapse = ", ")), con)
  for (col in cols) {
    vals <- unique(as.character(meta[[col]]))
    vals <- vals[!is.na(vals) & nzchar(vals)]
    vals <- sort(vals)
    writeLines(paste0("--- ", col, " (", length(vals), " labels)"), con)
    if (length(vals) == 0) {
      writeLines("", con)
    } else {
      writeLines(paste(vals, collapse = "|"), con)
    }
    writeLines("", con)
  }
}

message("Wrote: ", out_file)
```

### build_celltype_grouping.R

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(utils)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: build_celltype_grouping.R <rules_csv> [out_file]", call. = FALSE)
}

rules_csv <- args[[1]]
out_file <- ifelse(length(args) >= 2, args[[2]], "")

get_script_dir <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  return(getwd())
}

if (!nzchar(out_file)) {
  out_file <- file.path(get_script_dir(), "..", "resources", "celltype_grouping.txt")
}

if (!file.exists(rules_csv)) {
  stop("Rules CSV not found: ", rules_csv, call. = FALSE)
}

rules <- read.csv(rules_csv, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("group", "lineage", "azimuth_refs", "azimuth_label_patterns", "cellmarker_name_patterns", "notes")
missing <- setdiff(required_cols, colnames(rules))
if (length(missing) > 0) {
  stop("Rules CSV missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

con <- file(out_file, open = "w")
on.exit(close(con), add = TRUE)

writeLines(c(
  "# Cell type grouping rules (generated)",
  "#",
  "# Sources used to define Azimuth label strings:",
  "# - azimuth_label_catalog.txt (generated from installed references)",
  "#",
  "# Format per group:",
  "# group: <GROUP_ID>",
  "# azimuth_refs: comma-separated references",
  "# azimuth_label_patterns: comma-separated substrings (case-insensitive)",
  "# cellmarker_name_patterns: comma-separated substrings (case-insensitive)",
  "# lineage: coarse lineage label (used for final assignment)",
  "# notes: free text",
  ""
), con)

for (i in seq_len(nrow(rules))) {
  row <- rules[i, , drop = FALSE]
  writeLines(paste0("group: ", row$group), con)
  writeLines(paste0("azimuth_refs: ", row$azimuth_refs), con)
  writeLines(paste0("azimuth_label_patterns: ", row$azimuth_label_patterns), con)
  writeLines(paste0("cellmarker_name_patterns: ", row$cellmarker_name_patterns), con)
  writeLines(paste0("lineage: ", row$lineage), con)
  writeLines(paste0("notes: ", row$notes), con)
  writeLines("", con)
}

message("Wrote: ", out_file)
```

### Optional: regenerate label catalog + grouping rules
If you update Azimuth references or edit the grouping rules, regenerate these files:

```
Rscript /path/to/pipeline/r/Azimuth_label_catalog.R /path/to/pipeline/resources/azimuth_label_catalog.txt pbmcref,bonemarrowref,lungref,adiposeref,fetusref,liverref
Rscript /path/to/pipeline/r/build_celltype_grouping.R /path/to/pipeline/resources/celltype_grouping_rules.csv
```
## Annotation outputs (what to expect)
All outputs are written under `${QC_OUT}/../annotation_integration/seurat_annotation_integration/`:
- `merged_annotated.rds`: merged (non-integrated) object with module scores + Azimuth metadata.
- `merged_annotations.csv`: per-cell module scores and Azimuth columns.
- `celltype_assignment_by_cluster.csv`: cluster-level assignment summary (percentiles + evidence) with
  `assigned_group`, `assigned_lineage`, and `assigned_subtype`.
- `cellmarker_modules_map.csv`: module index → name map.
- `umap_pre/umap_pre_by_sample.png` and `umap_pre/umap_pre_by_celltype.png`.
- `featureplots_samples/<sample>/`: CellMarker + EWING feature plots on the merged UMAP.
- `featureplots_celltypes_samples/<sample>/`: per-celltype highlight plots on the merged UMAP.
- If integration is enabled, `integrated_annotated.rds` plus integrated UMAPs/feature plots.
