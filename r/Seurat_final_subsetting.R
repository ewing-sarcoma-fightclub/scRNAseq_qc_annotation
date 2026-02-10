#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

script_dir <- {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    getwd()
  }
}

source(file.path(script_dir, "utils_common.R"))
source(file.path(script_dir, "utils_doublet.R"))

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

add_qc_metrics <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.rps"]] <- PercentageFeatureSet(obj, pattern = "^RPS")
  obj[["percent.rpl"]] <- PercentageFeatureSet(obj, pattern = "^RPL")
  obj$percent.ribo <- obj$percent.rps + obj$percent.rpl
  obj
}

resolve_filtered_dir <- function(sample_dir) {
  filtered_dir <- file.path(sample_dir, "filtered_feature_bc_matrix")
  if (!dir.exists(filtered_dir)) {
    alt <- file.path(sample_dir, "outs", "filtered_feature_bc_matrix")
    if (dir.exists(alt)) filtered_dir <- alt
  }
  if (!dir.exists(filtered_dir)) return(NULL)
  filtered_dir
}

apply_soupx_if_available <- function(raw_filtered, qc_root, sample, use_soupx) {
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
  counts
}

load_metadata <- function(qc_root, sample, barcodes) {
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
  meta$orig.ident <- sample
  meta
}

load_emptydrops_calls <- function(qc_root, sample, barcodes, require_emptydrops) {
  ed_path <- file.path(qc_root, "emptydrops", sample, paste0(sample, "_emptydrops_results.csv"))
  if (file.exists(ed_path)) {
    ed <- read.csv(ed_path, row.names = 1)
    ed_call <- (!is.na(ed$FDR)) & (ed$FDR <= 0.01)
    names(ed_call) <- rownames(ed)
    out <- ed_call[barcodes]
  } else {
    if (require_emptydrops) {
      stop("EmptyDrops required but missing for ", sample, ": ", ed_path, call.=FALSE)
    }
    out <- rep(NA, length(barcodes))
  }
  names(out) <- barcodes
  out
}

load_doubletfinder_calls <- function(qc_root, sample, barcodes) {
  df_path <- file.path(qc_root, "doubletfinder", sample, paste0(sample, "_metadata_with_DoubletFinder.csv"))
  if (file.exists(df_path)) {
    df <- read.csv(df_path, row.names = 1)
    df_cols <- grep("^DF\\.classifications", colnames(df), value = TRUE)
    if (length(df_cols) >= 1) {
      if (length(df_cols) > 1) {
        message("[WARN] Multiple DF.classifications columns for ", sample,
                "; using last: ", df_cols[[length(df_cols)]])
      }
      df_col <- df_cols[[length(df_cols)]]
      df_call <- tolower(as.character(df[[df_col]])) == "doublet"
      names(df_call) <- rownames(df)
      out <- df_call[barcodes]
    } else {
      out <- rep(NA, length(barcodes))
    }
  } else {
    out <- rep(NA, length(barcodes))
  }
  names(out) <- barcodes
  out
}

load_dropletqc_calls <- function(qc_root, sample, barcodes, require_dropletqc) {
  dqc_path <- file.path(qc_root, "dropletqc", sample, paste0(sample, "_dropletqc_results.csv"))
  if (file.exists(dqc_path)) {
    dqc <- read.csv(dqc_path)
    if ("barcode" %in% colnames(dqc)) {
      rownames(dqc) <- dqc$barcode
    }
    if ("cell_status" %in% colnames(dqc)) {
      dqc_pass <- tolower(as.character(dqc$cell_status)) == "cell"
      names(dqc_pass) <- rownames(dqc)
      pass <- dqc_pass[barcodes]
      status <- dqc$cell_status[barcodes]
    } else {
      pass <- rep(NA, length(barcodes))
      status <- rep(NA, length(barcodes))
    }
  } else {
    if (require_dropletqc) {
      stop("DropletQC required but missing for ", sample, ": ", dqc_path, call.=FALSE)
    }
    pass <- rep(NA, length(barcodes))
    status <- rep(NA, length(barcodes))
  }
  names(pass) <- barcodes
  names(status) <- barcodes
  list(pass = pass, status = status)
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

doublet_filter_mode <- tolower(Sys.getenv("DOUBLET_FILTER_MODE", "doubletfinder"))
if (!(doublet_filter_mode %in% c("doubletfinder", "union", "intersection", "scrublet"))) {
  warning("[WARN] Unknown DOUBLET_FILTER_MODE='", doublet_filter_mode, "'. Using 'doubletfinder'.")
  doublet_filter_mode <- "doubletfinder"
}

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

  filtered_dir <- resolve_filtered_dir(sample_dir)
  if (is.null(filtered_dir)) {
    message("[WARN] filtered_feature_bc_matrix not found for ", sample)
    next
  }

  raw_filtered <- Read10X(filtered_dir)
  if (is.list(raw_filtered)) {
    raw_filtered <- raw_filtered[[1]]
  }

  barcodes <- colnames(raw_filtered)

  counts <- apply_soupx_if_available(raw_filtered, qc_root, sample, use_soupx)
  meta <- load_metadata(qc_root, sample, barcodes)
  meta$emptydrops_call <- load_emptydrops_calls(qc_root, sample, barcodes, require_emptydrops)
  meta$doubletfinder_call <- load_doubletfinder_calls(qc_root, sample, barcodes)

  scrublet_col <- NULL
  for (c in c("predicted_doublets", "predicted_doublet", "scrublet_doublet", "scrublet_call")) {
    if (c %in% colnames(meta)) {
      scrublet_col <- c
      break
    }
  }
  if (!is.null(scrublet_col)) {
    scrub_call <- coerce_logical_vec(meta[[scrublet_col]])
    meta$scrublet_call <- scrub_call
  } else {
    meta$scrublet_call <- NA
  }

  dqc_calls <- load_dropletqc_calls(qc_root, sample, barcodes, require_dropletqc)
  meta$dropletqc_pass <- dqc_calls$pass
  meta$dropletqc_status <- dqc_calls$status

  # Compute basic QC metrics on all cells used for downstream filtering
  tmp <- CreateSeuratObject(
    counts = counts,
    meta.data = meta,
    min.cells = 0,
    min.features = 0,
    project = sample
  )
  tmp <- add_qc_metrics(tmp)

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
  has_df <- !all(is.na(meta$doubletfinder_call))
  if (has_df) {
    df_call_full[barcodes] <- as.logical(meta$doubletfinder_call)
    df_call_full[is.na(df_call_full)] <- FALSE
  }

  # Scrublet calls (missing -> FALSE)
  scrub_call_full <- rep(FALSE, length(barcodes))
  names(scrub_call_full) <- barcodes
  has_scrublet <- !all(is.na(meta$scrublet_call))
  if (has_scrublet) {
    scrub_call_full[barcodes] <- as.logical(meta$scrublet_call)
    scrub_call_full[is.na(scrub_call_full)] <- FALSE
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

  doublet_fail <- df_call_full
  if (doublet_filter_mode == "scrublet") {
    if (has_scrublet) {
      doublet_fail <- scrub_call_full
    } else {
      warning("[WARN] DOUBLET_FILTER_MODE=scrublet but no scrublet calls found for ", sample, "; falling back to DoubletFinder.")
    }
  } else if (doublet_filter_mode == "union") {
    if (has_df && has_scrublet) {
      doublet_fail <- df_call_full | scrub_call_full
    } else if (has_scrublet) {
      warning("[WARN] DOUBLET_FILTER_MODE=union but DoubletFinder calls missing for ", sample, "; using scrublet only.")
      doublet_fail <- scrub_call_full
    } else if (!has_df) {
      warning("[WARN] DOUBLET_FILTER_MODE=union but no doublet calls found for ", sample, "; skipping doublet filter.")
      doublet_fail <- rep(FALSE, length(barcodes))
    }
  } else if (doublet_filter_mode == "intersection") {
    if (has_df && has_scrublet) {
      doublet_fail <- df_call_full & scrub_call_full
    } else if (has_df) {
      warning("[WARN] DOUBLET_FILTER_MODE=intersection but scrublet calls missing for ", sample, "; using DoubletFinder only.")
      doublet_fail <- df_call_full
    } else if (has_scrublet) {
      warning("[WARN] DOUBLET_FILTER_MODE=intersection but DoubletFinder calls missing for ", sample, "; using scrublet only.")
      doublet_fail <- scrub_call_full
    } else {
      warning("[WARN] DOUBLET_FILTER_MODE=intersection but no doublet calls found for ", sample, "; skipping doublet filter.")
      doublet_fail <- rep(FALSE, length(barcodes))
    }
  }
  keep <- keep & !doublet_fail

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

  seu <- add_qc_metrics(seu)

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
