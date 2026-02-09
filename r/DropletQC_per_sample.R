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
