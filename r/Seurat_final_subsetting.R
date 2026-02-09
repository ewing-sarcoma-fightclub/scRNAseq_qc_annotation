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
