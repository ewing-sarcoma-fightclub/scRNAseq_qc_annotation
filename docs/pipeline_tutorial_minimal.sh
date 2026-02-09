#!/bin/bash
set -euo pipefail

# Minimal pipeline: Cell Ranger -> QC (essential steps only)

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PIPELINE_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

# Edit these paths for your run
FASTQ_ROOT="/path/to/fastqs"
REF_DIR="/path/to/refdata-gex-GRCh38-2024-A"
CELLRANGER_OUT="${PIPELINE_ROOT}/outputs/cellranger"
QC_OUT="${PIPELINE_ROOT}/outputs/qc"
USE_SOUPX="true"
REQUIRE_EMPTYDROPS="true"
REQUIRE_DROPLETQC="true"  # set false if no BAMs (possorted_genome_bam.bam) are available

TMP_R_DIR="$(mktemp -d)"

cat <<'RSCRIPT' > "${TMP_R_DIR}/EmptyDrops_per_sample.R"
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DropletUtils)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)

sample_dir <- args[[1]]
sample_name <- args[[2]]
out_dir <- ifelse(length(args) >= 3, args[[3]], file.path(sample_dir, "EmptyDrops"))

raw_dir <- file.path(sample_dir, "outs", "raw_feature_bc_matrix")

# Read raw matrix
sce <- read10xCounts(raw_dir)

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
RSCRIPT

cat <<'RSCRIPT' > "${TMP_R_DIR}/SoupX_per_sample.R"
library(SoupX)
library(Seurat)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)

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

filtered_dir <- file.path(sample_dir, "outs", "filtered_feature_bc_matrix")
raw_dir <- file.path(sample_dir, "outs", "raw_feature_bc_matrix")
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

# 4. Save output (genes x cells to match 10x/Seurat expectations)
writeMM(out, file.path(out_path, paste0(sample_name, "_soupx.mtx")))

# Save SoupX summary stats
summary_path <- file.path(out_path, paste0(sample_name, "_soupx_summary.txt"))
summary_lines <- c(
  paste0("Sum of counts before: ", sum(toc_counts)),
  paste0("Sum of counts after: ", sum(out)),
  paste0("Pct of cells left: ", round(sum(out) / sum(toc_counts), digits = 3))
)
writeLines(summary_lines, summary_path)
RSCRIPT

cat <<'RSCRIPT' > "${TMP_R_DIR}/DoubletFinder_after_SoupX.R"
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)

sample_dir <- args[[1]]
sample_name <- args[[2]]
soupx_mtx <- ifelse(length(args) >= 3 && nzchar(args[[3]]), args[[3]], file.path(sample_dir, "SoupX", paste0(sample_name, "_soupx.mtx")))
out_dir <- ifelse(length(args) >= 4 && nzchar(args[[4]]), args[[4]], file.path(sample_dir, "DoubletFinder"))
doublet_rate <- ifelse(length(args) >= 5 && nzchar(args[[5]]), as.numeric(args[[5]]), 0.075)
npcs <- ifelse(length(args) >= 6 && nzchar(args[[6]]), as.integer(args[[6]]), 30)
metadata_path <- ifelse(length(args) >= 7 && nzchar(args[[7]]), args[[7]], file.path(sample_dir, "metadata.csv"))
min_features <- ifelse(length(args) >= 8 && nzchar(args[[8]]), as.integer(args[[8]]), 200)
min_cells <- ifelse(length(args) >= 9 && nzchar(args[[9]]), as.integer(args[[9]]), 3)
max_percent_mt <- ifelse(length(args) >= 10 && nzchar(args[[10]]), as.numeric(args[[10]]), Inf)
seed <- ifelse(length(args) >= 11 && nzchar(args[[11]]), as.integer(args[[11]]), 1)

# ---- input ----
filtered_dir <- file.path(sample_dir, "outs", "filtered_feature_bc_matrix")

# Read filtered 10x to get gene + barcode names
raw_filtered <- Read10X(filtered_dir)
if (is.list(raw_filtered)) {
  raw_filtered <- raw_filtered[[1]]
}

soupx_counts <- readMM(soupx_mtx)

rownames(soupx_counts) <- rownames(raw_filtered)
colnames(soupx_counts) <- colnames(raw_filtered)

meta <- read.csv(metadata_path, row.names = 1)

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
  best_pK <- 0.05
} else {
  best_pK <- pk[which.max(bc)]
  if (is.na(best_pK)) {
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
RSCRIPT

cat <<'RSCRIPT' > "${TMP_R_DIR}/DropletQC_per_sample.R"
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DropletQC)
  library(DropletUtils)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)

outs_dir <- args[[1]]
sample_name <- args[[2]]
out_dir <- ifelse(length(args) >= 3, args[[3]], file.path(outs_dir, "DropletQC"))
tiles <- ifelse(length(args) >= 4, as.integer(args[[4]]), 100)
cores <- ifelse(length(args) >= 5 && nzchar(args[[5]]), as.integer(args[[5]]), max(1, parallel::detectCores() - 1))
metadata_path <- ifelse(length(args) >= 6 && nzchar(args[[6]]), args[[6]], "")
metadata_col <- ifelse(length(args) >= 7 && nzchar(args[[7]]), args[[7]], "")

# Compute nuclear fraction from Cell Ranger outs
nf <- nuclear_fraction_tags(outs = outs_dir, tiles = tiles, cores = cores, verbose = FALSE)

# UMI counts from filtered matrix
filtered_dir <- file.path(outs_dir, "filtered_feature_bc_matrix")

sce <- read10xCounts(filtered_dir)
umi <- Matrix::colSums(counts(sce))

# Align barcodes
common <- intersect(rownames(nf), colnames(sce))

nf <- nf[common, , drop = FALSE]
umi <- umi[common]

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
}
ed$cell_type <- cell_type

# Identify damaged cells
res <- identify_damaged_cells(ed, verbose = FALSE, output_plots = FALSE)

# res[[1]] is data frame with cell_status
qc_df <- res[[1]]
qc_df$barcode <- rownames(qc_df)

out_csv <- file.path(out_dir, paste0(sample_name, "_dropletqc_results.csv"))
write.csv(qc_df, out_csv, row.names = FALSE)
RSCRIPT

cat <<'RSCRIPT' > "${TMP_R_DIR}/Seurat_final_subsetting.R"
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

safe_cell_cycle_scoring <- function(obj, s_genes, g2m_genes, label = "dataset") {
  if (length(s_genes) == 0 || length(g2m_genes) == 0) {
    return(obj)
  }
  for (nbin in c(24, 10, 5)) {
    res <- tryCatch(
      CellCycleScoring(obj, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE, nbin = nbin),
      error = function(e) e
    )
    if (!inherits(res, "error")) return(res)
  }
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

sample_dirs <- list.dirs(root_dir, full.names = TRUE, recursive = FALSE)

for (sample_dir in sample_dirs) {
  sample <- basename(sample_dir)

  filtered_dir <- file.path(sample_dir, "outs", "filtered_feature_bc_matrix")

  raw_filtered <- Read10X(filtered_dir)
  if (is.list(raw_filtered)) {
    raw_filtered <- raw_filtered[[1]]
  }

  barcodes <- colnames(raw_filtered)

  counts <- raw_filtered
  soupx_path <- file.path(qc_root, "soupx", sample, paste0(sample, "_soupx.mtx"))
  if (use_soupx) {
    soupx_counts <- readMM(soupx_path)
    rownames(soupx_counts) <- rownames(raw_filtered)
    colnames(soupx_counts) <- colnames(raw_filtered)
    counts <- soupx_counts
  }

  meta <- data.frame(row.names = barcodes)

  meta_path <- file.path(qc_root, "seurat_metadata", sample, "metadata.csv")
  meta_in <- read.csv(meta_path, row.names = 1)
  meta <- cbind(meta, meta_in[barcodes, , drop = FALSE])

  # Ensure orig.ident is the sample name (metadata.csv may carry a generic value)
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

  keep <- basic_qc_keep
  if (require_emptydrops) {
    keep <- keep & !is.na(meta$emptydrops_call) & as.logical(meta$emptydrops_call)
  }
  if (require_dropletqc) {
    keep <- keep & !is.na(meta$dropletqc_pass) & as.logical(meta$dropletqc_pass)
  }
  keep <- keep & !df_call_full

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
    if (length(s_genes) > 0 && length(g2m_genes) > 0) {
      seu <- NormalizeData(seu)
      seu <- safe_cell_cycle_scoring(seu, s_genes, g2m_genes, label = paste0(sample, " (RNA)"))
    } else {
      do_regress_cell_cycle_sample <- FALSE
    }
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

# ---- Cell Ranger ----
while IFS= read -r -d '' sample_dir; do
  sample_name=$(basename "${sample_dir}")
  (
    cd "${CELLRANGER_OUT}"
    cellranger count \
      --id="${sample_name}" \
      --transcriptome="${REF_DIR}" \
      --fastqs="${sample_dir}" \
      --include-introns=true \
      --create-bam true
  )
done < <(find "${FASTQ_ROOT}" -mindepth 1 -maxdepth 1 -type d -print0)

# ---- QC after Cell Ranger ----
EMPTYDROPS_OUT="${QC_OUT}/emptydrops"
SEURAT_META_OUT="${QC_OUT}/seurat_metadata"
SOUPX_OUT="${QC_OUT}/soupx"
DOUBLETFINDER_OUT="${QC_OUT}/doubletfinder"
DROPLETQC_OUT="${QC_OUT}/dropletqc"
SEURAT_FINAL_OUT="${QC_OUT}/seurat_qc"

# 1) EmptyDrops
while IFS= read -r -d '' sample_dir; do
  sample_name=$(basename "${sample_dir}")
  out_dir="${EMPTYDROPS_OUT}/${sample_name}"
  Rscript "${TMP_R_DIR}/EmptyDrops_per_sample.R" "${sample_dir}" "${sample_name}" "${out_dir}"
done < <(find "${CELLRANGER_OUT}" -mindepth 1 -maxdepth 1 -type d -print0)

# 2) SoupX per sample
while IFS= read -r -d '' sample_dir; do
  sample_name=$(basename "${sample_dir}")
  metadata_csv="${SEURAT_META_OUT}/${sample_name}/metadata.csv"
  out_dir="${SOUPX_OUT}/${sample_name}"
  Rscript "${TMP_R_DIR}/SoupX_per_sample.R" "${sample_dir}" "${sample_name}" "${metadata_csv}" "${out_dir}"
done < <(find "${CELLRANGER_OUT}" -mindepth 1 -maxdepth 1 -type d -print0)

# 3) DoubletFinder
while IFS= read -r -d '' sample_dir; do
  sample_name=$(basename "${sample_dir}")
  soupx_mtx="${SOUPX_OUT}/${sample_name}/${sample_name}_soupx.mtx"
  metadata_csv="${SEURAT_META_OUT}/${sample_name}/metadata.csv"
  out_dir="${DOUBLETFINDER_OUT}/${sample_name}"
  Rscript "${TMP_R_DIR}/DoubletFinder_after_SoupX.R" "${sample_dir}" "${sample_name}" "${soupx_mtx}" "${out_dir}" \
    "0.075" "30" "${metadata_csv}"
done < <(find "${CELLRANGER_OUT}" -mindepth 1 -maxdepth 1 -type d -print0)

# 4) DropletQC
while IFS= read -r -d '' sample_dir; do
  sample_name=$(basename "${sample_dir}")
  bam_path="${sample_dir}/outs/possorted_genome_bam.bam"
  if [[ ! -f "${bam_path}" ]]; then
    echo "[WARN] DropletQC skipped for ${sample_name}: possorted_genome_bam.bam not found."
    continue
  fi
  out_dir="${DROPLETQC_OUT}/${sample_name}"
  metadata_csv="${SEURAT_META_OUT}/${sample_name}/metadata.csv"
  Rscript "${TMP_R_DIR}/DropletQC_per_sample.R" "${sample_dir}/outs" "${sample_name}" "${out_dir}" "100" "" "${metadata_csv}"
done < <(find "${CELLRANGER_OUT}" -mindepth 1 -maxdepth 1 -type d -print0)

# 5) Final Seurat subsetting + clustering
Rscript "${TMP_R_DIR}/Seurat_final_subsetting.R" "${CELLRANGER_OUT}" "${QC_OUT}" "${SEURAT_FINAL_OUT}" \
  "${USE_SOUPX}" "${REQUIRE_EMPTYDROPS}" "${REQUIRE_DROPLETQC}"
