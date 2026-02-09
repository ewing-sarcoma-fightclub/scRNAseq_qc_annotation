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
