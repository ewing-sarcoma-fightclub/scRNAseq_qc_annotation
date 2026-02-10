#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

script_dir <- {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    getwd()
  }
}

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
    table_path <- select_rate_table(chem_use, script_dir)
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
