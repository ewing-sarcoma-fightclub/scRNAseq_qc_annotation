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
