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
  "Monocyte",
  "Myeloid cell",
  "Naive CD4+ T cell",
  "Progenitor cell",
  "Red blood cell (erythrocyte)",
  "Regulatory CD4+ T cell",
  "T cell",
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
