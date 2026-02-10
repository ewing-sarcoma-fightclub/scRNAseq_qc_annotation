#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
  library(Azimuth)
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

read_gene_id_map <- function(path) {
  if (!file.exists(path)) {
    stop("Gene ID map not found: ", path, call.=FALSE)
  }
  sep <- if (grepl("\\.(tsv|txt|tab)$", path, ignore.case = TRUE)) "\t" else ","
  map <- tryCatch(
    read.csv(path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(map) || ncol(map) < 2) {
    map <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  cols_lc <- tolower(colnames(map))
  ensg_idx <- which(cols_lc %in% c("ensembl", "ensembl_id", "ensembl_gene_id", "gene_id") |
                      grepl("ensembl", cols_lc) | grepl("^ens", cols_lc))[1]
  sym_idx <- which(grepl("symbol", cols_lc) | grepl("gene_name", cols_lc) | grepl("hgnc", cols_lc))[1]
  if (is.na(ensg_idx) || is.na(sym_idx)) {
    if (ncol(map) == 2) {
      ensg_idx <- 1
      sym_idx <- 2
    } else {
      stop("Gene ID map must include Ensembl and symbol columns.", call.=FALSE)
    }
  }
  ensg <- sub("\\..*$", "", as.character(map[[ensg_idx]]))
  sym <- trimws(as.character(map[[sym_idx]]))
  keep <- nzchar(ensg) & nzchar(sym)
  setNames(sym[keep], ensg[keep])
}

aggregate_matrix_by_rows <- function(mat, new_names, fun = "sum") {
  if (length(new_names) == 0 || nrow(mat) == 0) return(mat)
  if (nrow(mat) != length(new_names)) {
    stop("aggregate_matrix_by_rows: length mismatch between matrix rows and names.", call.=FALSE)
  }
  grp <- factor(new_names, levels = unique(new_names))
  if (nlevels(grp) == length(grp)) {
    rownames(mat) <- new_names
    return(mat)
  }
  G <- Matrix::sparseMatrix(
    i = as.integer(grp),
    j = seq_along(grp),
    x = 1,
    dims = c(nlevels(grp), length(grp))
  )
  agg <- G %*% mat
  if (fun == "mean") {
    counts <- as.numeric(table(grp))
    agg <- Matrix::Diagonal(x = 1 / counts) %*% agg
  }
  rownames(agg) <- levels(grp)
  if (!is.null(colnames(mat))) colnames(agg) <- colnames(mat)
  agg
}

read_gene_list_file <- function(path) {
  if (!file.exists(path)) stop("Gene list file not found: ", path, call.=FALSE)
  lines <- readLines(path, warn = FALSE)
  tokens <- unlist(strsplit(lines, "[,;\\t ]+"))
  tokens <- trimws(tokens)
  tokens <- tokens[nzchar(tokens)]

  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv", "tsv", "txt", "tab")) {
    sep <- if (ext == "csv") "," else "\t"
    df <- tryCatch(read.csv(path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE),
                   error = function(e) NULL)
    if (is.null(df)) {
      df <- tryCatch(read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
                     error = function(e) NULL)
    }
    if (!is.null(df) && ncol(df) >= 1) {
      vals <- unlist(lapply(df, as.character), use.names = FALSE)
      vals <- trimws(vals)
      vals <- vals[nzchar(vals)]
      tokens <- c(tokens, vals)
    }
  }

  tokens <- unique(tokens)
  tokens <- tokens[grepl("[A-Za-z]", tokens)]
  drop <- tolower(tokens) %in% c(
    "gene", "genes", "symbol", "symbols", "gene_symbol", "gene_symbols",
    "gene_name", "gene_names", "feature", "features"
  )
  tokens <- tokens[!drop]
  tokens
}

get_assay_layers <- function(obj, assay_name) {
  if (!(assay_name %in% Assays(obj))) return(character())
  assay <- obj[[assay_name]]

  has_slot <- function(x, s) {
    s %in% methods::slotNames(x)
  }

  layers <- character()
  if (requireNamespace("SeuratObject", quietly = TRUE) &&
      exists("Layers", where = asNamespace("SeuratObject"), inherits = FALSE)) {
    layers <- tryCatch(SeuratObject::Layers(assay), error = function(e) character())
  } else if (exists("Layers", where = asNamespace("Seurat"), inherits = FALSE)) {
    layers <- tryCatch(Seurat::Layers(assay), error = function(e) character())
  }

  layers <- unique(as.character(layers))
  layers <- layers[nzchar(layers)]

  if (!("counts" %in% layers) && has_slot(assay, "counts")) layers <- c(layers, "counts")
  if (!("data" %in% layers) && has_slot(assay, "data")) layers <- c(layers, "data")
  if (!("scale.data" %in% layers) && has_slot(assay, "scale.data")) layers <- c(layers, "scale.data")

  unique(layers)
}

get_layer_data <- function(obj, assay_name, layer) {
  tryCatch(
    Seurat::GetAssayData(obj, assay = assay_name, layer = layer),
    error = function(e) {
      tryCatch(
        Seurat::GetAssayData(obj, assay = assay_name, slot = layer),
        error = function(e2) NULL
      )
    }
  )
}

set_layer_data <- function(obj, assay_name, layer, mat) {
  tryCatch(
    Seurat::SetAssayData(obj, assay = assay_name, layer = layer, new.data = mat),
    error = function(e) Seurat::SetAssayData(obj, assay = assay_name, slot = layer, new.data = mat)
  )
}

rename_features_in_object <- function(obj, map_vec) {
  assay_names <- names(obj@assays)
  for (assay_name in assay_names) {
    layers <- unique(c(get_assay_layers(obj, assay_name), "counts", "data", "scale.data"))
    new_features <- character()
    for (layer in layers) {
      mat <- get_layer_data(obj, assay_name, layer)
      if (is.null(mat) || nrow(mat) == 0) next
      old <- rownames(mat)
      base <- sub("\\..*$", "", old)
      mapped <- unname(map_vec[base])
      new <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, old)
      fun <- if (layer == "counts") "sum" else "mean"
      mat2 <- aggregate_matrix_by_rows(mat, new, fun = fun)
      obj <- set_layer_data(obj, assay_name, layer, mat2)
      if (length(new_features) == 0) new_features <- rownames(mat2)
    }

    if (length(new_features) == 0) {
      assay <- obj[[assay_name]]
      new_features <- tryCatch(rownames(assay), error = function(e) character())
    }
    if (length(new_features) == 0) next

    assay <- obj[[assay_name]]
    tryCatch({
      assay@meta.features <- data.frame(row.names = new_features)
    }, error = function(e) {})
    obj[[assay_name]] <- assay

    vf <- tryCatch(VariableFeatures(obj, assay = assay_name), error = function(e) character())
    if (length(vf) > 0) {
      base_vf <- sub("\\..*$", "", vf)
      mapped_vf <- unname(map_vec[base_vf])
      vf_new <- ifelse(!is.na(mapped_vf) & nzchar(mapped_vf), mapped_vf, vf)
      vf_new <- unique(vf_new)
      vf_new <- vf_new[vf_new %in% new_features]
      VariableFeatures(obj, assay = assay_name) <- vf_new
    }
  }
  obj
}

maybe_map_ensembl_to_symbols <- function(objs, map_path, allow_ensg = FALSE) {
  features_all <- sort(unique(unlist(lapply(objs, rownames))))
  ensg_frac <- mean(grepl("^ENSG", features_all))
  if (!is.finite(ensg_frac) || ensg_frac < 0.5) {
    return(list(objs = objs, features_all = features_all, mapped = FALSE))
  }
  if (!nzchar(map_path)) {
    if (allow_ensg) {
      message("[WARN] Feature names look like Ensembl IDs; proceeding without mapping because ALLOW_ENSG_IDS=true.")
      return(list(objs = objs, features_all = features_all, mapped = FALSE))
    }
    stop(
      "Feature names look like Ensembl IDs. Provide ENSEMBL_TO_SYMBOL_MAP to map ENSG→symbol, ",
      "or set ALLOW_ENSG_IDS=true to proceed without mapping.",
      call.=FALSE
    )
  }
  map_vec <- read_gene_id_map(map_path)
  if (length(map_vec) == 0) stop("Gene ID map is empty: ", map_path, call.=FALSE)

  total_features <- 0L
  mapped_features <- 0L
  duplicate_features <- 0L
  for (nm in names(objs)) {
    obj <- objs[[nm]]
    old <- rownames(obj)
    base <- sub("\\..*$", "", old)
    mapped <- unname(map_vec[base])
    total_features <- total_features + length(old)
    mapped_features <- mapped_features + sum(!is.na(mapped) & nzchar(mapped))
    duplicate_features <- duplicate_features + sum(duplicated(ifelse(!is.na(mapped) & nzchar(mapped), mapped, old)))
    objs[[nm]] <- rename_features_in_object(obj, map_vec)
  }
  message("[INFO] Ensembl→symbol mapping: mapped ", mapped_features, "/", total_features,
          " features; duplicates aggregated: ", duplicate_features, ".")
  dup_max <- get_env_num("ENSEMBL_SYMBOL_DUP_MAX", -1)
  if (dup_max >= 0 && duplicate_features > dup_max) {
    stop(
      "Ensembl→symbol mapping produced ", duplicate_features,
      " duplicate symbols (aggregated). Set ENSEMBL_SYMBOL_DUP_MAX to allow, or raise threshold.",
      call.=FALSE
    )
  }
  features_all <- sort(unique(unlist(lapply(objs, rownames))))
  list(objs = objs, features_all = features_all, mapped = TRUE)
}

get_expr_for_aucell <- function(obj, assay = NULL, slot = "data") {
  if (is.null(assay)) assay <- DefaultAssay(obj)
  m <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = slot),
    error = function(e) Seurat::GetAssayData(obj, assay = assay, slot = slot)
  )
  if (min(m) < 0) {
    stop(sprintf("AUCell input has negative values (assay=%s slot=%s). Use slot='counts' or non-negative 'data'.", assay, slot))
  }
  m
}

filter_genesets_present <- function(geneSets, universe, min_genes = 10) {
  gs2 <- lapply(geneSets, function(g) intersect(unique(g), universe))
  keep <- vapply(gs2, function(g) length(g) >= min_genes, logical(1))
  gs2[keep]
}

run_aucell <- function(expr, geneSets, out_rds_prefix,
                       aucMaxRank = NULL, nCores = 1, plotStats = FALSE) {
  if (!requireNamespace("AUCell", quietly = TRUE)) {
    stop("AUCell not installed. Install via BiocManager::install('AUCell').")
  }
  if (is.null(aucMaxRank)) aucMaxRank <- ceiling(nrow(expr) * 0.05)

  rankings_rds <- paste0(out_rds_prefix, ".rankings.rds")
  auc_rds <- paste0(out_rds_prefix, ".auc.rds")

  if (file.exists(rankings_rds)) {
    rankings <- readRDS(rankings_rds)
  } else {
    rankings <- AUCell::AUCell_buildRankings(expr, nCores = nCores, plotStats = plotStats)
    saveRDS(rankings, rankings_rds)
  }

  if (file.exists(auc_rds)) {
    cellsAUC <- readRDS(auc_rds)
  } else {
    cellsAUC <- AUCell::AUCell_calcAUC(geneSets, rankings, aucMaxRank = aucMaxRank, nCores = nCores)
    saveRDS(cellsAUC, auc_rds)
  }

  list(rankings = rankings, cellsAUC = cellsAUC, aucMaxRank = aucMaxRank)
}

aucell_threshold_and_assign <- function(cellsAUC, nCores = 1,
                                       thrP = 0.01, smallestPopPercent = 0.25,
                                       densAdjust = 2, nBreaks = 100,
                                       plotHist = FALSE) {
  if (!requireNamespace("AUCell", quietly = TRUE)) {
    stop("AUCell not installed. Install via BiocManager::install('AUCell').")
  }
  thrObj <- AUCell::AUCell_exploreThresholds(
    cellsAUC,
    thrP = thrP,
    smallestPopPercent = smallestPopPercent,
    densAdjust = densAdjust,
    nBreaks = nBreaks,
    nCores = nCores,
    plotHist = plotHist,
    assignCells = TRUE
  )

  thresholds <- AUCell::getThresholdSelected(thrObj)
  assignments <- AUCell::getAssignments(thrObj)

  list(thrObj = thrObj, thresholds = thresholds, assignments = assignments)
}

attach_aucell_to_seurat <- function(obj, cellsAUC, thresholds, prefix = "AUC_") {
  auc_mat <- AUCell::getAUC(cellsAUC)

  for (gs in rownames(auc_mat)) {
    col <- paste0(prefix, gs)
    obj[[col]] <- as.numeric(auc_mat[gs, colnames(obj)])
  }

  for (gs in names(thresholds)) {
    thr <- thresholds[[gs]]
    col_pos <- paste0(prefix, gs, "_pos")
    scores <- as.numeric(auc_mat[gs, colnames(obj)])
    obj[[col_pos]] <- scores >= thr
  }

  obj
}

azimuth_dominance_qc <- function(obj, sample, refs, out_dir) {
  hi <- get_env_num("AZIMUTH_HIGHCONF", 0.75)
  dom_thr <- get_env_num("AZIMUTH_DOMINANCE_FRAC", 0.80)

  score_cols <- vapply(refs, function(r) paste0("azimuth_", r, "_mapping.score"), character(1))
  present <- score_cols[score_cols %in% colnames(obj@meta.data)]
  if (length(present) < 2) return(invisible(NULL))

  refs_present <- sub("^azimuth_(.*)_mapping\\.score$", "\\1", present)
  mat <- do.call(cbind, lapply(present, function(cn) suppressWarnings(as.numeric(obj@meta.data[[cn]]))))
  colnames(mat) <- refs_present

  max_score <- apply(mat, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  max_ref <- apply(mat, 1, function(x) {
    if (all(is.na(x))) return(NA_character_)
    refs_present[which.max(replace(x, is.na(x), -Inf))][1]
  })

  hi_idx <- which(!is.na(max_score) & max_score >= hi)
  if (length(hi_idx) == 0) return(invisible(NULL))

  tab <- sort(table(max_ref[hi_idx]), decreasing = TRUE)
  frac <- as.numeric(tab) / sum(tab)

  qc <- data.frame(
    sample = sample,
    ref = names(tab),
    n_hi = as.integer(tab),
    frac_hi = frac,
    highconf_threshold = hi,
    stringsAsFactors = FALSE
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(qc, file.path(out_dir, paste0(sample, "_azimuth_dominance.csv")), row.names = FALSE)

  if (length(frac) >= 1 && frac[1] >= dom_thr) {
    warning("[Azimuth QC] Sample ", sample, ": reference '", names(tab)[1],
            "' dominates high-confidence calls (", sprintf("%.1f%%", 100 * frac[1]),
            "). Check reference suitability / gene naming / tissue mismatch.")
  }
  invisible(qc)
}

ensure_aucell_columns <- function(obj, score_cols) {
  for (col in score_cols) {
    if (!(col %in% colnames(obj@meta.data))) {
      obj@meta.data[[col]] <- NA_real_
    }
  }
  for (col in paste0(score_cols, "_pos")) {
    if (!(col %in% colnames(obj@meta.data))) {
      obj@meta.data[[col]] <- NA
    }
  }
  if (!("EWING_1" %in% colnames(obj@meta.data))) {
    obj@meta.data[["EWING_1"]] <- NA_real_
  }
  if (!("EWING_1_pos" %in% colnames(obj@meta.data))) {
    obj@meta.data[["EWING_1_pos"]] <- NA
  }
  obj
}

write_run_manifest <- function(out_dir, args, objs, merged = NULL, integrated = NULL) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
  script_path <- {
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", cmd_args, value = TRUE)
    if (length(file_arg) > 0) {
      normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)
    } else {
      NA_character_
    }
  }

  info_df <- data.frame(
    section = "run",
    key = c("timestamp", "script", "args", "n_samples", "n_cells_merged", "n_cells_integrated"),
    value = c(
      ts,
      script_path,
      paste(args, collapse = " "),
      length(objs),
      if (!is.null(merged)) ncol(merged) else NA,
      if (!is.null(integrated)) ncol(integrated) else NA
    ),
    stringsAsFactors = FALSE
  )

  pkg_list <- c("Seurat", "SeuratObject", "Azimuth", "AUCell", "mclust", "HGNChelper")
  pkg_versions <- vapply(pkg_list, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "not_installed"
  }, character(1))
  pkg_df <- data.frame(
    section = "package",
    key = paste0("pkg_", names(pkg_versions)),
    value = unname(pkg_versions),
    stringsAsFactors = FALSE
  )

  env <- Sys.getenv()
  keep_keys <- names(env)[grepl("^(AUCELL|CELL_|SEURAT_|AZIMUTH_|ANNOTATION_|DOUBLET_|ENSEMBL_|ALLOW_ENSG_IDS|RESUME|PIPELINE_ROOT|RARE_IMMUNE)", names(env))]
  keep_keys <- sort(unique(keep_keys))
  env_df <- data.frame(
    section = "env",
    key = keep_keys,
    value = env[keep_keys],
    stringsAsFactors = FALSE
  )

  manifest_df <- rbind(info_df, pkg_df, env_df)
  write.table(
    manifest_df,
    file.path(out_dir, "run_manifest.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
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

get_grouping_file_path <- function(out_dir) {
  gf <- Sys.getenv("CELLTYPE_GROUPING_FILE", "")
  if (!nzchar(gf)) {
    gf <- file.path(script_dir, "..", "resources", "celltype_grouping.txt")
  }
  gf
}

collect_cellmarker_patterns <- function(grouping_file) {
  if (!file.exists(grouping_file)) return(character())
  groups <- parse_grouping_file(grouping_file)
  pats <- unique(unlist(lapply(groups, function(g) g$cellmarker_name_patterns)))
  pats <- pats[!is.na(pats) & nzchar(pats)]
  pats
}

jaccard <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(1)
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(0)
  inter / uni
}

prune_modules <- function(modules, grouping_file, max_jaccard = 0.7) {
  report <- data.frame(
    module = character(),
    status = character(),
    stage = character(),
    reason = character(),
    n_genes = integer(),
    compared_to = character(),
    jaccard = double(),
    stringsAsFactors = FALSE
  )

  keep_by_grouping <- tolower(Sys.getenv("CELL_MARKER_KEEP_BY_GROUPING", "true"))
  keep_by_grouping <- keep_by_grouping %in% c("true", "t", "1", "yes", "y")
  if (keep_by_grouping) {
    pats <- collect_cellmarker_patterns(grouping_file)
    if (length(pats) > 0) {
      keep <- match_any_pattern(names(modules), pats)
      dropped <- names(modules)[!keep]
      if (length(dropped) > 0) {
        report <- rbind(
          report,
          data.frame(
            module = dropped,
            status = "dropped",
            stage = "grouping",
            reason = "not_in_grouping",
            n_genes = vapply(modules[dropped], length, integer(1)),
            compared_to = "",
            jaccard = NA_real_,
            stringsAsFactors = FALSE
          )
        )
      }
      modules <- modules[keep]
    }
  }

  if (length(modules) <= 1) {
    if (length(modules) == 1) {
      report <- rbind(
        report,
        data.frame(
          module = names(modules),
          status = "kept",
          stage = "final",
          reason = "kept",
          n_genes = vapply(modules, length, integer(1)),
          compared_to = "",
          jaccard = NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }
    return(list(modules = modules, report = report))
  }
  ord <- order(vapply(modules, length, integer(1)), decreasing = TRUE)
  mods <- modules[ord]

  kept <- list()
  for (nm in names(mods)) {
    g <- mods[[nm]]
    redundant <- FALSE
    redundant_against <- ""
    redundant_score <- NA_real_
    if (length(kept) > 0) {
      for (knm in names(kept)) {
        jac <- jaccard(g, kept[[knm]])
        if (jac >= max_jaccard) {
          redundant <- TRUE
          redundant_against <- knm
          redundant_score <- jac
          break
        }
      }
    }
    if (redundant) {
      report <- rbind(
        report,
        data.frame(
          module = nm,
          status = "dropped",
          stage = "jaccard",
          reason = paste0("jaccard>=", max_jaccard),
          n_genes = length(g),
          compared_to = redundant_against,
          jaccard = redundant_score,
          stringsAsFactors = FALSE
        )
      )
    } else {
      kept[[nm]] <- g
      report <- rbind(
        report,
        data.frame(
          module = nm,
          status = "kept",
          stage = "final",
          reason = "kept",
          n_genes = length(g),
          compared_to = "",
          jaccard = NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  list(modules = kept, report = report)
}

get_azimuth_score_columns <- function(meta, ref, level) {
  label_col <- paste0("azimuth_", ref, "_predicted.", level)
  if (!(label_col %in% colnames(meta))) {
    return(list(label_col = label_col, score_col = NA_character_))
  }
  score_col <- paste0("azimuth_", ref, "_predicted.", level, ".score")
  if (!(score_col %in% colnames(meta))) {
    map_col <- paste0("azimuth_", ref, "_mapping.score")
    if (map_col %in% colnames(meta)) {
      score_col <- map_col
    } else {
      score_col <- NA_character_
    }
  }
  list(label_col = label_col, score_col = score_col)
}

calibrate_azimuth_thresholds <- function(meta, azimuth_levels, refs = NULL, method = "mixture",
                                         posterior = 0.95, min_n = 200, fallback_pct = 0.85) {
  if (is.null(refs)) refs <- names(azimuth_levels)
  refs <- unique(refs)
  if (!is.finite(posterior) || posterior <= 0 || posterior >= 1) posterior <- 0.95
  if (!is.finite(min_n) || min_n < 10) min_n <- 10
  if (!is.finite(fallback_pct) || fallback_pct <= 0 || fallback_pct >= 1) fallback_pct <- 0.85
  method <- tolower(method)
  mclust_ok <- TRUE
  if (method == "mixture") {
    mclust_ok <- requireNamespace("mclust", quietly = TRUE)
    if (mclust_ok) {
      suppressPackageStartupMessages(library(mclust))
    } else {
      warning("[Azimuth] mclust not installed; falling back to quantile thresholds.")
    }
  }

  thresholds <- list()
  rows <- list()
  for (ref in refs) {
    level <- azimuth_levels[[ref]]
    if (is.null(level)) next
    cols <- get_azimuth_score_columns(meta, ref, level)
    if (is.na(cols$score_col) || !(cols$score_col %in% colnames(meta))) next

    scores <- suppressWarnings(as.numeric(meta[[cols$score_col]]))
    scores <- scores[is.finite(scores)]
    n_scores <- length(scores)
    thr <- NA_real_
    method_used <- NA_character_

    if (n_scores == 0) {
      method_used <- "no_scores"
    } else if (method == "mixture" && n_scores >= min_n) {
      if (!mclust_ok) {
        method_used <- "quantile_no_mclust"
      } else {
        eps <- 1e-6
        use_logit <- all(scores > 0 & scores < 1)
        x <- if (use_logit) qlogis(pmin(pmax(scores, eps), 1 - eps)) else scores
        fit <- tryCatch(mclust::Mclust(x, G = 2, verbose = FALSE), error = function(e) NULL)
        if (is.null(fit) || is.null(fit$z) || length(fit$parameters$mean) < 2) {
          method_used <- "quantile_fit_failed"
        } else {
          hi_comp <- which.max(fit$parameters$mean)
          post <- fit$z[, hi_comp]
          keep <- which(is.finite(post) & post >= posterior)
          if (length(keep) == 0) {
            method_used <- "quantile_no_posterior"
          } else {
            thr <- suppressWarnings(min(scores[keep], na.rm = TRUE))
            method_used <- "mixture_post"
          }
        }
      }
    } else if (method == "mixture" && n_scores < min_n) {
      method_used <- "quantile_min_n"
    } else {
      method_used <- "quantile_manual"
    }

    if (!is.finite(thr) && n_scores > 0) {
      thr <- suppressWarnings(as.numeric(quantile(scores, probs = fallback_pct, na.rm = TRUE, names = FALSE)))
      if (is.na(method_used) || !nzchar(method_used)) method_used <- "quantile_fallback"
    }

    thresholds[[ref]] <- thr
    rows[[length(rows) + 1]] <- data.frame(
      reference = ref,
      score_col = cols$score_col,
      n_scores = n_scores,
      method = method_used,
      threshold = thr,
      posterior = posterior,
      min_n = min_n,
      fallback_pct = fallback_pct,
      stringsAsFactors = FALSE
    )
  }
  out <- if (length(rows) > 0) do.call(rbind, rows) else data.frame()
  list(thresholds = thresholds, summary = out)
}

summarize_azimuth_references <- function(meta, azimuth_levels, refs = NULL, thresholds = NULL, highconf_pct = 0.9) {
  if (is.null(refs)) refs <- names(azimuth_levels)
  refs <- unique(refs)
  rows <- list()
  for (ref in refs) {
    level <- azimuth_levels[[ref]]
    if (is.null(level)) next
    cols <- get_azimuth_score_columns(meta, ref, level)
    label_col <- cols$label_col
    score_col <- cols$score_col
    if (is.na(score_col) || !(label_col %in% colnames(meta)) || !(score_col %in% colnames(meta))) next
    labels <- as.character(meta[[label_col]])
    scores <- suppressWarnings(as.numeric(meta[[score_col]]))
    has_label <- !is.na(labels) & nzchar(labels)
    n_cells <- sum(has_label)
    score_vals <- scores[has_label & is.finite(scores)]
    thresh <- NA_real_
    if (!is.null(thresholds) && !is.null(thresholds[[ref]]) &&
        length(thresholds[[ref]]) == 1 && is.finite(thresholds[[ref]])) {
      thresh <- thresholds[[ref]]
    } else if (length(score_vals) > 0) {
      thresh <- suppressWarnings(as.numeric(quantile(score_vals, probs = highconf_pct, na.rm = TRUE, names = FALSE)))
    }
    high_conf <- has_label & is.finite(scores) & is.finite(thresh) & scores >= thresh
    n_high <- sum(high_conf)
    high_frac <- if (n_cells > 0) n_high / n_cells else NA_real_
    top_label <- NA_character_
    top_label_frac <- NA_real_
    if (n_high > 0) {
      tab <- sort(table(labels[high_conf]), decreasing = TRUE)
      top_label <- names(tab)[1]
      top_label_frac <- as.numeric(tab[1]) / n_high
    }
    rows[[length(rows) + 1]] <- data.frame(
      reference = ref,
      label_level = level,
      label_col = label_col,
      score_col = score_col,
      n_cells = n_cells,
      n_high_conf = n_high,
      high_conf_frac = high_frac,
      score_threshold = thresh,
      top_label = top_label,
      top_label_frac = top_label_frac,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(data.frame())
  out <- do.call(rbind, rows)
  total_high <- sum(out$n_high_conf, na.rm = TRUE)
  if (is.finite(total_high) && total_high > 0) {
    out$high_conf_share <- out$n_high_conf / total_high
  } else {
    out$high_conf_share <- NA_real_
  }
  out
}

normalize_label <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

normalize_assigned_label <- function(x) {
  x <- as.character(x)
  x[x == "T cell"] <- "T cell lineage"
  x
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

add_group_summary_columns <- function(cluster_summary, group_names, summary_lists) {
  for (suffix in names(summary_lists)) {
    vals <- summary_lists[[suffix]]
    for (gname in group_names) {
      cluster_summary[[paste0(gname, "_", suffix)]] <- vals[[gname]]
    }
  }
  cluster_summary
}

select_de_assay <- function(obj, assay_env = "") {
  if (nzchar(assay_env) && assay_env %in% Assays(obj)) return(assay_env)
  for (a in c("RNA", "SCT")) {
    if (a %in% Assays(obj)) return(a)
  }
  DefaultAssay(obj)
}

downsample_by_cluster <- function(meta, cluster_col, max_per_cluster, seed = 1) {
  if (!(cluster_col %in% colnames(meta))) return(NULL)
  if (!is.finite(max_per_cluster) || max_per_cluster < 1) return(NULL)
  set.seed(seed)
  keep <- character()
  clusts <- unique(meta[[cluster_col]])
  for (cl in clusts) {
    cells <- rownames(meta)[meta[[cluster_col]] == cl]
    if (length(cells) <= max_per_cluster) {
      keep <- c(keep, cells)
    } else {
      keep <- c(keep, sample(cells, max_per_cluster))
    }
  }
  keep
}

downsample_total_cells <- function(meta, cluster_col, cells, max_total, seed = 1) {
  if (!is.finite(max_total) || max_total < 1) return(cells)
  if (length(cells) <= max_total) return(cells)
  set.seed(seed)
  if (!(cluster_col %in% colnames(meta))) {
    return(sample(cells, max_total))
  }
  cl <- meta[cells, cluster_col]
  cl <- as.character(cl)
  uniq <- unique(cl)
  if (max_total < length(uniq)) {
    return(sample(cells, max_total))
  }
  # Keep at least one cell per cluster if possible
  keep <- character()
  for (u in uniq) {
    idx <- which(cl == u)
    keep <- c(keep, cells[sample(idx, 1)])
  }
  remaining <- setdiff(cells, keep)
  if (length(keep) < max_total && length(remaining) > 0) {
    keep <- c(keep, sample(remaining, min(length(remaining), max_total - length(keep))))
  }
  keep
}

ensure_join_layers <- function(obj, assay = NULL) {
  if (is.null(assay)) assay <- DefaultAssay(obj)
  if (!(assay %in% Assays(obj))) return(obj)
  # JoinLayers is required for DE on Seurat v5 objects with multiple layers
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    obj <- tryCatch(
      SeuratObject::JoinLayers(obj, assay = assay),
      error = function(e) obj
    )
  } else if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
    obj <- tryCatch(
      Seurat::JoinLayers(obj, assay = assay),
      error = function(e) obj
    )
  }
  obj
}

ensure_data_layer <- function(obj, assay = NULL) {
  if (is.null(assay)) assay <- DefaultAssay(obj)
  if (!(assay %in% Assays(obj))) return(obj)
  data_mat <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) tryCatch(
      Seurat::GetAssayData(obj, assay = assay, slot = "data"),
      error = function(e2) NULL
    )
  )
  if (is.null(data_mat) || nrow(data_mat) == 0) {
    obj <- tryCatch(
      Seurat::NormalizeData(obj, assay = assay, verbose = FALSE),
      error = function(e) {
        warning("[Validation] NormalizeData failed for assay ", assay, ": ", conditionMessage(e))
        obj
      }
    )
  }
  obj
}

build_group_marker_sets <- function(module_cols_by_group, modules_mapped, module_names_safe) {
  mods_by_safe <- setNames(modules_mapped, module_names_safe)
  out <- list()
  for (gname in names(module_cols_by_group)) {
    mods <- module_cols_by_group[[gname]]
    if (is.null(mods) || length(mods) == 0) next
    genes <- unique(unlist(mods_by_safe[mods]))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    out[[gname]] <- genes
  }
  out
}

build_top_markers <- function(de_markers, fc_col = NULL, top_n = 50) {
  if (is.null(de_markers) || nrow(de_markers) == 0 || !("cluster" %in% colnames(de_markers))) {
    return(list())
  }
  top <- list()
  clusters <- unique(as.character(de_markers$cluster))
  for (cl in clusters) {
    mk <- de_markers[as.character(de_markers$cluster) == cl, , drop = FALSE]
    if (nrow(mk) == 0) next
    if (!is.null(fc_col) && fc_col %in% colnames(mk)) {
      mk <- mk[order(-mk[[fc_col]]), , drop = FALSE]
    }
    genes <- head(as.character(mk$gene), top_n)
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) == 0) next
    top[[cl]] <- list(
      genes = genes,
      genes_norm = unique(toupper(genes))
    )
  }
  top
}

save_feature_plots <- function(obj, features, out_dir, title_suffix = "", reduction = "umap") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (feat in features) {
    if (!(feat %in% colnames(obj@meta.data))) next
    vals <- obj@meta.data[[feat]]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0 || length(unique(vals)) < 2) next
    p <- if (is.null(reduction)) {
      FeaturePlot(obj, features = feat)
    } else {
      FeaturePlot(obj, features = feat, reduction = reduction)
    }
    ggplot2::ggsave(
      file.path(out_dir, paste0(feat, ".png")),
      p + ggplot2::ggtitle(paste0(feat, title_suffix)),
      width = 6, height = 5
    )
  }
}

ensure_umap <- function(obj, npcs, seed.use = 1, assay_prefer = c("SCT", "RNA")) {
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
    } else {
      if (length(VariableFeatures(obj)) == 0) {
        VariableFeatures(obj) <- rownames(obj)
      }
    }

    feats <- VariableFeatures(obj)
    if (length(feats) == 0) {
      feats <- rownames(obj)
    }
    obj <- RunPCA(obj, npcs = npcs, features = feats, verbose = FALSE)
  }

  pca_emb <- tryCatch(Embeddings(obj, "pca"), error = function(e) NULL)
  npcs_avail <- if (!is.null(pca_emb)) min(npcs, ncol(pca_emb)) else npcs
  if (!is.finite(npcs_avail) || npcs_avail < 2) npcs_avail <- 2

  if (!("umap" %in% Reductions(obj))) {
    obj <- tryCatch(
      RunUMAP(obj, dims = 1:npcs_avail, seed.use = seed.use, verbose = FALSE),
      error = function(e) RunUMAP(obj, dims = 1:npcs_avail, verbose = FALSE)
    )
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

aucell_slot <- Sys.getenv("AUCELL_SLOT", "data")
aucell_min_genes <- get_env_int("AUCELL_MIN_GENES", 10, min = 1)
aucell_max_rank_env <- get_env_int("AUCELL_MAX_RANK", 0, min = 0)
aucell_max_rank <- if (is.finite(aucell_max_rank_env) && aucell_max_rank_env > 0) aucell_max_rank_env else NULL
aucell_ncores <- get_env_int("AUCELL_NCORES", 1, min = 1)

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

drop_env <- Sys.getenv("CELL_MARKER_DROP_MODULES", "")
if (nzchar(drop_env)) {
  extra <- trimws(unlist(strsplit(drop_env, ",")))
  extra <- extra[nzchar(extra)]
  drop_modules <- unique(extra)
  modules <- modules[!(names(modules) %in% drop_modules)]
}
grouping_file_early <- get_grouping_file_path(out_dir)
max_j <- get_env_num("CELL_MARKER_MAX_JACCARD", 0.7)
if (!is.finite(max_j) || max_j <= 0 || max_j >= 1) max_j <- 0.7
prune_res <- prune_modules(modules, grouping_file = grouping_file_early, max_jaccard = max_j)
modules <- prune_res$modules
prune_report <- prune_res$report
if (length(modules) == 0) {
  stop("All CellMarker modules were pruned; check grouping patterns / thresholds.")
}

ewing_genes <- tryCatch({
  read_gene_list_file(ewing_sig_csv)
}, error = function(e) {
  as.character(read.csv(ewing_sig_csv, header = FALSE, stringsAsFactors = FALSE)[[1]])
})
ewing_genes <- trimws(ewing_genes)
ewing_genes <- ewing_genes[nzchar(ewing_genes)]
ewing_genes <- ewing_genes[grepl("[A-Za-z]", ewing_genes)]
if (length(ewing_genes) == 0) {
  ewing_genes <- as.character(read.csv(ewing_sig_csv, header = FALSE, stringsAsFactors = FALSE)[[1]])
  ewing_genes <- trimws(ewing_genes)
  ewing_genes <- ewing_genes[nzchar(ewing_genes)]
}

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

# Map Ensembl IDs to gene symbols when needed (avoid silently degraded module mapping)
ens_map_path <- Sys.getenv("ENSEMBL_TO_SYMBOL_MAP", "")
allow_ensg <- parse_bool(Sys.getenv("ALLOW_ENSG_IDS", "false"), FALSE)
map_res <- maybe_map_ensembl_to_symbols(objs, ens_map_path, allow_ensg = allow_ensg)
objs <- map_res$objs
features_all <- map_res$features_all

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

annotation_helpers_dir <- file.path(out_dir, "annotation_helper_files")
dir.create(annotation_helpers_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(module_map, file.path(annotation_helpers_dir, "cell_marker_modules_map.csv"), row.names = FALSE)
if (exists("prune_report") && is.data.frame(prune_report)) {
  write.csv(prune_report, file.path(annotation_helpers_dir, "cellmarker_module_pruning.csv"), row.names = FALSE)
}

# Per-sample annotation (plots generated later on merged UMAP)
score_cols <- paste0("CM_", module_names_safe)
plot_features_all <- c(score_cols, "EWING_1")
annotated_samples_dir <- file.path(out_dir, "annotated_sampled_rds")
dir.create(annotated_samples_dir, recursive = TRUE, showWarnings = FALSE)
doublet_summary_rows <- list()
annotation_doublet_filter <- parse_bool(Sys.getenv("ANNOTATION_DOUBLET_FILTER", "false"), FALSE)
annotation_doublet_require <- parse_bool(Sys.getenv("ANNOTATION_DOUBLET_REQUIRE", "false"), FALSE)
annotation_doublet_mode <- tolower(Sys.getenv(
  "ANNOTATION_DOUBLET_MODE",
  Sys.getenv("DOUBLET_FILTER_MODE", "doubletfinder")
))

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

# Azimuth label columns to use (stable levels)
azimuth_levels <- list(
  pbmcref = "celltype.l2",
  bonemarrowref = "celltype.l2",
  lungref = "ann_level_3",
  adiposeref = "celltype.l1",
  fetusref = "annotation.l1",
  liverref = "celltype.l2"
)

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
  skip_azimuth <- FALSE
  if (resume && file.exists(annotated_path)) {
    obj <- readRDS(annotated_path)
    skip_azimuth <- TRUE
  } else {
    obj <- objs[[sample]]
  }

  dbl_res <- resolve_doublet_calls(obj@meta.data)
  dbl_comb <- combine_doublet_calls(dbl_res$df_call, dbl_res$scrub_call, mode = annotation_doublet_mode)
  n_cells_pre <- ncol(obj)
  n_df <- if (dbl_comb$has_df) sum(dbl_res$df_call == TRUE, na.rm = TRUE) else NA_integer_
  n_scrub <- if (dbl_comb$has_scrub) sum(dbl_res$scrub_call == TRUE, na.rm = TRUE) else NA_integer_
  n_union <- sum(dbl_comb$union_call, na.rm = TRUE)
  n_intersection <- sum(dbl_comb$intersection_call, na.rm = TRUE)
  n_mode <- sum(dbl_comb$doublet_call, na.rm = TRUE)

  if (annotation_doublet_filter) {
    if (!(dbl_comb$has_df || dbl_comb$has_scrub)) {
      msg <- paste0("[Doublet] No doublet calls found for ", sample, "; cannot enforce doublet filtering.")
      if (annotation_doublet_require) {
        stop(msg, call.=FALSE)
      } else {
        warning(msg)
      }
    } else {
      obj$doublet_call <- dbl_comb$doublet_call
      keep_cells <- rownames(obj@meta.data)[!dbl_comb$doublet_call]
      obj <- subset(obj, cells = keep_cells)
    }
  }

  if (!skip_azimuth) {
    if ("SCT" %in% Assays(obj)) {
      DefaultAssay(obj) <- "SCT"
    } else {
      DefaultAssay(obj) <- "RNA"
    }
    genes_in_obj <- rownames(obj)
    module_present_counts <- vapply(modules_mapped, function(g) sum(g %in% genes_in_obj), integer(1))
    modules_present <- modules_mapped[module_present_counts > 0]

    present_names_safe <- module_names_safe[match(names(modules_present), names(modules_mapped))]
    cm_sets <- setNames(modules_present, paste0("CM_", present_names_safe))
    all_sets <- c(cm_sets, list(EWING_1 = ewing_genes))

    expr <- get_expr_for_aucell(obj, assay = DefaultAssay(obj), slot = aucell_slot)
    geneSets <- filter_genesets_present(all_sets, universe = rownames(expr), min_genes = aucell_min_genes)

    if (length(geneSets) > 0) {
      aucell_dir <- file.path(annotation_helpers_dir, "aucell")
      dir.create(aucell_dir, recursive = TRUE, showWarnings = FALSE)

      au_out <- run_aucell(
        expr = expr,
        geneSets = geneSets,
        out_rds_prefix = file.path(aucell_dir, paste0(sample, ".aucell")),
        aucMaxRank = aucell_max_rank,
        nCores = aucell_ncores,
        plotStats = FALSE
      )

      thr_out <- aucell_threshold_and_assign(
        au_out$cellsAUC,
        nCores = aucell_ncores,
        plotHist = FALSE
      )
      thr_df <- data.frame(
        signature = names(thr_out$thresholds),
        threshold = as.numeric(unlist(thr_out$thresholds)),
        stringsAsFactors = FALSE
      )
      write.csv(
        thr_df,
        file.path(aucell_dir, paste0(sample, ".aucell_thresholds.csv")),
        row.names = FALSE
      )

      obj <- attach_aucell_to_seurat(obj, au_out$cellsAUC, thr_out$thresholds, prefix = "")
    } else {
      message("[WARN] No gene sets available for AUCell in ", sample, "; CM_/EWING scores set to NA.")
    }
  }

  dbl_res_post <- resolve_doublet_calls(obj@meta.data)
  dbl_comb_post <- combine_doublet_calls(dbl_res_post$df_call, dbl_res_post$scrub_call, mode = annotation_doublet_mode)
  n_cells_post <- ncol(obj)
  doublet_summary_rows[[length(doublet_summary_rows) + 1]] <- data.frame(
    sample = sample,
    n_cells = n_cells_post,
    n_df_doublet_pre = n_df,
    n_scrublet_doublet_pre = n_scrub,
    n_doublet_union_pre = n_union,
    n_doublet_intersection_pre = n_intersection,
    n_doublet_mode_pre = n_mode,
    n_df_doublet = if (dbl_comb_post$has_df) sum(dbl_res_post$df_call == TRUE, na.rm = TRUE) else NA_integer_,
    n_scrublet_doublet = if (dbl_comb_post$has_scrub) sum(dbl_res_post$scrub_call == TRUE, na.rm = TRUE) else NA_integer_,
    n_doublet_union = sum(dbl_comb_post$union_call, na.rm = TRUE),
    n_doublet_intersection = sum(dbl_comb_post$intersection_call, na.rm = TRUE),
    n_doublet_mode = sum(dbl_comb_post$doublet_call, na.rm = TRUE),
    doublet_mode = annotation_doublet_mode,
    filtered = annotation_doublet_filter,
    n_cells_pre = n_cells_pre,
    n_removed = if (annotation_doublet_filter) (n_cells_pre - n_cells_post) else NA_integer_,
    stringsAsFactors = FALSE
  )

  # Ensure module score columns exist even when resuming
  obj <- ensure_aucell_columns(obj, score_cols)

  # Azimuth mapping per-sample (avoid multi-SCT merge issues)
  if (skip_azimuth) {
    message("[INFO] Skipping Azimuth for ", sample, " (resume: annotated file exists)")
  } else if (length(azimuth_refs) > 0) {
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
  az_qc_dir <- file.path(out_dir, "azimuth_qc")
  azimuth_dominance_qc(obj, sample, azimuth_refs, az_qc_dir)

  # Save per-sample annotated object
  saveRDS(obj, annotated_path)

  objs[[sample]] <- obj
}

if (length(doublet_summary_rows) > 0) {
  doublet_summary <- do.call(rbind, doublet_summary_rows)
  write.csv(doublet_summary, file.path(annotation_helpers_dir, "doublet_summary.csv"), row.names = FALSE)
}

# Merge non-integrated (after per-sample scoring)
merged <- objs[[1]]
if (length(objs) > 1) {
  merged <- merge(merged, y = objs[2:length(objs)])
}

# ---- Azimuth threshold calibration (per reference) ----
az_thresh_method <- tolower(Sys.getenv("AZIMUTH_THRESHOLD_METHOD", "mixture"))
az_thresh_posterior <- get_env_num("AZIMUTH_THRESHOLD_POSTERIOR", 0.95)
az_thresh_min_n <- get_env_int("AZIMUTH_THRESHOLD_MIN_N", 200)
az_fallback_pct <- get_env_num_fallback(
  "AZIMUTH_THRESHOLD_FALLBACK_PCT",
  "CELLTYPE_GROUPING_PERCENTILE",
  0.85
)
az_thresh_out <- calibrate_azimuth_thresholds(
  merged@meta.data,
  azimuth_levels,
  refs = azimuth_refs,
  method = az_thresh_method,
  posterior = az_thresh_posterior,
  min_n = az_thresh_min_n,
  fallback_pct = az_fallback_pct
)
az_ref_thresholds <- az_thresh_out$thresholds
if (nrow(az_thresh_out$summary) > 0) {
  write.csv(
    az_thresh_out$summary,
    file.path(annotation_helpers_dir, "azimuth_thresholds.csv"),
    row.names = FALSE
  )
}

# ---- Azimuth reference sanity summary ----
az_highconf_pct <- get_env_num("AZIMUTH_REF_HIGHCONF_PCT", 0.9)
if (!is.finite(az_highconf_pct) || az_highconf_pct <= 0 || az_highconf_pct >= 1) {
  az_highconf_pct <- 0.9
}
az_dom_share <- get_env_num("AZIMUTH_REF_DOMINANCE_SHARE", 0.8)
if (!is.finite(az_dom_share) || az_dom_share <= 0 || az_dom_share > 1) {
  az_dom_share <- 0.8
}
az_dom_min <- get_env_int("AZIMUTH_REF_DOMINANCE_MIN", 200, min = 1)
az_refs_for_summary <- intersect(azimuth_refs, names(azimuth_levels))
az_summary <- summarize_azimuth_references(
  merged@meta.data,
  azimuth_levels,
  refs = az_refs_for_summary,
  thresholds = az_ref_thresholds,
  highconf_pct = az_highconf_pct
)
if (nrow(az_summary) > 0) {
  write.csv(az_summary, file.path(annotation_helpers_dir, "azimuth_reference_summary.csv"), row.names = FALSE)
  total_high <- sum(az_summary$n_high_conf, na.rm = TRUE)
  if (is.finite(total_high) && total_high >= az_dom_min) {
    idx <- which.max(az_summary$high_conf_share)
    dom_share <- az_summary$high_conf_share[idx]
    if (is.finite(dom_share) && dom_share >= az_dom_share) {
      warning(
        "[Azimuth] High-confidence labels are dominated by reference '",
        az_summary$reference[idx], "' (share=", sprintf("%.2f", dom_share),
        ", n=", az_summary$n_high_conf[idx], "). Validate reference suitability."
      )
    }
  }
}

# ---- Cluster-level group assignment (CellMarker + Azimuth) ----
grouping_file <- get_grouping_file_path(out_dir)
min_mod_cells <- get_env_int("CELLTYPE_GROUPING_MIN_MOD_CELLS", 10, min = 1)
min_az_cells <- get_env_int("CELLTYPE_GROUPING_MIN_AZ_CELLS", 10, min = 1)
min_mod_frac <- get_env_num("CELLTYPE_GROUPING_MIN_MOD_FRAC", 0.10, min = 0, max = 1)
min_az_frac <- get_env_num("CELLTYPE_GROUPING_MIN_AZ_FRAC", 0.10, min = 0, max = 1)
module_only_min_cells <- get_env_int("CELLTYPE_GROUPING_MODULE_ONLY_MIN_CELLS", 20, min = 1)
module_only_min_frac <- get_env_num("CELLTYPE_GROUPING_MODULE_ONLY_MIN_FRAC", 0.25, min = 0, max = 1)
module_only_min_pct <- get_env_num(
  "CELLTYPE_GROUPING_MODULE_ONLY_MIN_PCT",
  90,
  min = 0,
  max = 100,
  invalid_value = NA_real_
)
module_active_min_cells <- get_env_int("CELL_MARKER_MIN_ACTIVE_CELLS", 0, min = 0)
module_active_min_frac <- get_env_num("CELL_MARKER_MIN_ACTIVE_FRAC", 0, min = 0, max = 1)
rare_immune_enable <- parse_bool(Sys.getenv("RARE_IMMUNE_ENABLE", "false"), FALSE)
rare_immune_max_cells <- get_env_int(
  "RARE_IMMUNE_MAX_CELLS",
  50,
  min = 1,
  invalid_value = NA_integer_
)
rare_immune_max_frac <- get_env_num(
  "RARE_IMMUNE_MAX_FRAC",
  0.02,
  min = 0,
  max = 1,
  invalid_value = NA_real_
)
if (!is.na(rare_immune_max_frac) && (rare_immune_max_frac <= 0 || rare_immune_max_frac >= 1)) {
  rare_immune_max_frac <- NA_real_
}
rare_immune_min_az_cells <- get_env_int("RARE_IMMUNE_MIN_AZ_CELLS", 3, min = 1)
rare_immune_min_az_frac <- get_env_num("RARE_IMMUNE_MIN_AZ_FRAC", 0.02, min = 0, max = 1)
az_nonimmune_pct <- get_env_num("CELLTYPE_GROUPING_AZ_PCT_NONIMM", 85, min = 0, max = 100)
az_lineage_delta <- get_env_num("CELLTYPE_GROUPING_AZ_LINEAGE_DELTA", 0.10, min = 0)
az_lineage_min <- get_env_num("CELLTYPE_GROUPING_AZ_LINEAGE_MIN", 0.20, min = 0)
ewing_min_cells <- get_env_int("CELLTYPE_GROUPING_EWING_MIN_CELLS", 20, min = 1)
ewing_min_frac <- get_env_num("CELLTYPE_GROUPING_EWING_MIN_FRAC", 0.20, min = 0, max = 1)

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

    module_keep <- NULL
    if (exists("module_map") && nrow(module_map) > 0) {
      pos_cols <- paste0("CM_", module_map$module_name_safe, "_pos")
      active_n <- vapply(pos_cols, function(pc) {
        if (pc %in% colnames(meta)) {
          sum(as.logical(meta[[pc]]), na.rm = TRUE)
        } else {
          0
        }
      }, numeric(1))
      active_frac <- vapply(pos_cols, function(pc) {
        if (pc %in% colnames(meta)) {
          mean(as.logical(meta[[pc]]), na.rm = TRUE)
        } else {
          0
        }
      }, numeric(1))
      module_activity <- data.frame(
        module_name = module_map$module_name,
        module_name_safe = module_map$module_name_safe,
        active_n = as.integer(active_n),
        active_frac = active_frac,
        stringsAsFactors = FALSE
      )
      if (module_active_min_cells > 0 || module_active_min_frac > 0) {
        keep <- rep(TRUE, nrow(module_activity))
        if (module_active_min_cells > 0) {
          keep <- keep & module_activity$active_n >= module_active_min_cells
        }
        if (module_active_min_frac > 0) {
          keep <- keep & module_activity$active_frac >= module_active_min_frac
        }
        module_activity$keep <- keep
        if (sum(keep) == 0) {
          warning("[Grouping] All modules failed activity thresholds; keeping all modules.")
          module_keep <- module_activity$module_name_safe
          module_activity$keep <- TRUE
        } else {
          module_keep <- module_activity$module_name_safe[keep]
        }
      } else {
        module_activity$keep <- TRUE
        module_keep <- module_activity$module_name_safe
      }
      write.csv(
        module_activity,
        file.path(annotation_helpers_dir, "cellmarker_module_activity.csv"),
        row.names = FALSE
      )
      if (length(module_keep) > 0) {
        module_cols_by_group <- lapply(module_cols_by_group, function(cols) cols[cols %in% module_keep])
      }
    }

    # Precompute per-cell group azimuth scores (max across refs)
    group_azimuth_score <- list()
    for (gname in names(groups)) {
      g <- groups[[gname]]
      az_scores <- rep(NA_real_, nrow(meta))
      if (!is.null(g$azimuth_refs) && length(g$azimuth_refs) > 0 && length(g$azimuth_label_patterns) > 0) {
        for (ref in g$azimuth_refs) {
          level <- azimuth_levels[[ref]]
          if (is.null(level)) next
          cols <- get_azimuth_score_columns(meta, ref, level)
          label_col <- cols$label_col
          score_col <- cols$score_col
          if (is.na(score_col) || !(label_col %in% colnames(meta)) || !(score_col %in% colnames(meta))) next
          labels <- as.character(meta[[label_col]])
          scores <- suppressWarnings(as.numeric(meta[[score_col]]))
          match <- match_any_pattern(labels, g$azimuth_label_patterns, normalize = TRUE)
          scores[!match] <- NA_real_
          ref_thresh <- az_ref_thresholds[[ref]]
          if (length(ref_thresh) == 1 && is.finite(ref_thresh)) {
            scores[scores < ref_thresh] <- NA_real_
          }
          if (all(is.na(scores))) next
          if (all(is.na(az_scores))) {
            az_scores <- scores
          } else {
            both_na <- is.na(az_scores) & is.na(scores)
            az_scores <- pmax(az_scores, scores, na.rm = TRUE)
            az_scores[both_na] <- NA_real_
          }
        }
      }
      group_azimuth_score[[gname]] <- az_scores
    }

    # Cluster-level medians
    cluster_ids <- sort(unique(cluster_id))
    cluster_n <- vapply(cluster_ids, function(cid) sum(cluster_id == cid), integer(1))
    total_cells <- nrow(meta)
    cluster_summary <- data.frame(
      cluster_id = cluster_ids,
      n_cells = cluster_n,
      stringsAsFactors = FALSE
    )

    # EWING signature medians + AUCell-positive fractions
    ewing_med <- tapply(meta$EWING_1, cluster_id, median, na.rm = TRUE)
    ewing_med <- ewing_med[cluster_ids]
    ewing_frac <- rep(NA_real_, length(cluster_ids))
    ewing_n <- tapply(!is.na(meta$EWING_1_pos), cluster_id, sum)
    ewing_n <- ewing_n[cluster_ids]
    if ("EWING_1_pos" %in% colnames(meta)) {
      ewing_frac <- tapply(as.logical(meta$EWING_1_pos), cluster_id, mean, na.rm = TRUE)
      ewing_frac <- ewing_frac[cluster_ids]
    }
    cluster_summary$EWING_1_median <- ewing_med
    cluster_summary$EWING_1_percentile <- NA_real_
    cluster_summary$EWING_1_frac <- ewing_frac
    cluster_summary$EWING_1_n <- ewing_n
    cluster_summary$EWING_1_threshold <- NA_real_

    # Module-level thresholds (per signature; AUCell-positive calls)
    module_cols <- character()
    if (exists("module_map")) {
      if (!is.null(module_keep) && length(module_keep) > 0) {
        module_cols <- paste0("CM_", module_keep)
      } else {
        module_cols <- paste0("CM_", module_map$module_name_safe)
      }
      module_cols <- module_cols[module_cols %in% colnames(meta)]
    }
    module_frac <- list()
    module_n <- list()
    module_med <- list()
    module_pct <- list()
    module_strong <- list()
    module_strict <- list()
    for (col in module_cols) {
      name_safe <- sub("^CM_", "", col)
      scores <- meta[[col]]
      pos_col <- paste0(col, "_pos")
      pos_vals <- if (pos_col %in% colnames(meta)) as.logical(meta[[pos_col]]) else rep(NA, nrow(meta))
      mod_n <- tapply(!is.na(pos_vals), cluster_id, sum)
      mod_n <- mod_n[cluster_ids]
      mod_frac <- tapply(pos_vals, cluster_id, mean, na.rm = TRUE)
      mod_frac <- mod_frac[cluster_ids]
      mod_med <- tapply(scores, cluster_id, median, na.rm = TRUE)
      mod_med <- mod_med[cluster_ids]
      mod_pct <- mod_frac * 100
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

    # Group azimuth medians/percentiles + support
    group_names <- names(groups)
    group_az_medians <- list()
    group_az_pct <- list()
    group_az_frac <- list()
    group_az_n <- list()
    group_az_strong <- list()
    group_az_support <- list()
    group_az_support_rare <- list()
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

      az_n <- tapply(!is.na(az_score), cluster_id, sum)
      az_n <- az_n[cluster_ids]

      az_frac <- tapply(!is.na(az_score), cluster_id, mean, na.rm = TRUE)
      az_frac <- az_frac[cluster_ids]
      group_az_n[[gname]] <- az_n
      group_az_frac[[gname]] <- az_frac

      group_az_strong[[gname]] <- !is.na(az_n) & !is.na(az_frac) & az_n >= min_az_cells & az_frac >= min_az_frac
      az_support <- group_az_strong[[gname]]
      if (!(group_lineage[[gname]] %in% immune_lineages)) {
        az_support <- az_support | (!is.na(az_pct) & az_pct >= az_nonimmune_pct & !is.na(az_n) & az_n >= min_az_cells)
      }
      group_az_support[[gname]] <- az_support
      group_az_support_rare[[gname]] <- !is.na(az_n) & !is.na(az_frac) &
        az_n >= rare_immune_min_az_cells & az_frac >= rare_immune_min_az_frac

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
      labels_all <- character()
      scores_all <- numeric()
      if (!is.null(g$azimuth_refs) && length(g$azimuth_refs) > 0 && length(g$azimuth_label_patterns) > 0) {
        for (ref in g$azimuth_refs) {
          level <- azimuth_levels[[ref]]
          if (is.null(level)) next
          cols <- get_azimuth_score_columns(meta, ref, level)
          label_col <- cols$label_col
          score_col <- cols$score_col
          if (is.na(score_col) || !(label_col %in% colnames(meta)) || !(score_col %in% colnames(meta))) next
          labels <- as.character(meta[[label_col]][idx])
          scores <- suppressWarnings(as.numeric(meta[[score_col]][idx]))
          match <- match_any_pattern(labels, g$azimuth_label_patterns, normalize = TRUE)
          labels <- labels[match]
          scores <- scores[match]
          if (length(labels) == 0) next
          ref_thresh <- az_ref_thresholds[[ref]]
          if (length(ref_thresh) == 1 && is.finite(ref_thresh)) {
            keep <- !is.na(scores) & scores >= ref_thresh
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
      rare_cluster <- FALSE
      if (rare_immune_enable) {
        rare_by_n <- if (!is.na(rare_immune_max_cells)) cluster_n[i] <= rare_immune_max_cells else FALSE
        rare_by_frac <- if (!is.na(rare_immune_max_frac) && total_cells > 0) (cluster_n[i] / total_cells) <= rare_immune_max_frac else FALSE
        rare_cluster <- rare_by_n || rare_by_frac
      }
      for (lin in names(lineage_groups)) {
        groups_lin <- lineage_groups[[lin]]
        mods <- unique(unlist(lapply(groups_lin, function(g) module_cols_by_group[[g]])))
        mods <- mods[mods %in% names(module_strong)]
        strong_mods <- mods[vapply(mods, function(m) isTRUE(module_strong[[m]][i]), logical(1))]
        has_mod <- length(strong_mods) > 0
        has_az <- any(vapply(groups_lin, function(g) isTRUE(group_az_support[[g]][i]), logical(1)))
        has_az_rare <- any(vapply(groups_lin, function(g) isTRUE(group_az_support_rare[[g]][i]), logical(1)))
        if (lin %in% immune_lineages) {
          if (has_mod && has_az) {
            lineage_mods[[lin]] <- strong_mods
          } else if (rare_cluster && has_az_rare) {
            lineage_mods[[lin]] <- character(0)
          }
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
                       ewing_n[i] >= ewing_min_cells && ewing_frac[i] >= ewing_min_frac)

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

    assigned_group_primary <- assigned_group
    assigned_lineage_primary <- assigned_lineage
    assigned_subtype_primary <- assigned_subtype
    cluster_summary$assigned_group_primary <- assigned_group_primary
    cluster_summary$assigned_lineage_primary <- assigned_lineage_primary
    cluster_summary$assigned_subtype_primary <- assigned_subtype_primary

    # Add per-group summary to cluster table
    cluster_summary <- add_group_summary_columns(
      cluster_summary,
      group_names,
      list(
        azimuth_median = group_az_medians,
        azimuth_pct = group_az_pct,
        azimuth_frac = group_az_frac,
        azimuth_n = group_az_n,
        azimuth_strong = group_az_strong,
        module_strong_n = group_module_strong_n,
        module_strong_any = group_module_strong_any,
        module_strict_n = group_module_strict_n,
        module_strict_any = group_module_strict_any
      )
    )

    # Make cluster_id_for_assignment available for DE before triage updates
    merged@meta.data <- meta

    # ---- Orthogonal validation: cluster markers + marker overlap ----
    de_assay <- select_de_assay(merged, Sys.getenv("ANNOTATION_DE_ASSAY", ""))
    de_min_pct <- get_env_num("ANNOTATION_DE_MIN_PCT", 0.25, min = 0)
    de_logfc <- get_env_num("ANNOTATION_DE_LOGFC", 0.25, min = 0)
    de_max_per_cluster <- get_env_int(
      "ANNOTATION_DE_MAX_CELLS_PER_CLUSTER",
      2000,
      min = 1,
      invalid_value = NA_integer_
    )
    de_max_total <- get_env_int(
      "ANNOTATION_DE_MAX_TOTAL_CELLS",
      50000,
      min = 1,
      invalid_value = NA_integer_
    )
    de_top_n <- get_env_int("ANNOTATION_DE_TOP_N", 50, min = 1)
    de_seed <- get_env_int("ANNOTATION_DE_SEED", 1)

    de_group_col <- if ("cluster_id_for_assignment" %in% colnames(merged@meta.data)) {
      "cluster_id_for_assignment"
    } else if ("seurat_clusters" %in% colnames(merged@meta.data)) {
      "seurat_clusters"
    } else {
      NULL
    }

    if (!is.null(de_group_col)) {
      de_obj <- merged
      keep_cells <- downsample_by_cluster(merged@meta.data, de_group_col, de_max_per_cluster, seed = de_seed)
      if (!is.null(keep_cells)) {
        keep_cells <- downsample_total_cells(merged@meta.data, de_group_col, keep_cells, de_max_total, seed = de_seed)
        de_obj <- subset(merged, cells = keep_cells)
      }
      DefaultAssay(de_obj) <- de_assay
      de_obj <- ensure_join_layers(de_obj, assay = de_assay)
      if (identical(de_assay, "SCT")) {
        if ("PrepSCTFindMarkers" %in% getNamespaceExports("Seurat")) {
          de_obj <- tryCatch(
            Seurat::PrepSCTFindMarkers(de_obj),
            error = function(e) {
              warning("[Validation] PrepSCTFindMarkers failed: ", conditionMessage(e))
              de_obj
            }
          )
        } else {
          warning("[Validation] PrepSCTFindMarkers not available in Seurat; SCT DE may fail.")
        }
      } else {
        de_obj <- ensure_data_layer(de_obj, assay = de_assay)
      }

      de_markers <- tryCatch(
        Seurat::FindAllMarkers(
          de_obj,
          only.pos = TRUE,
          min.pct = de_min_pct,
          logfc.threshold = de_logfc,
          group.by = de_group_col
        ),
        error = function(e) {
          warning("[Validation] FindAllMarkers failed: ", conditionMessage(e))
          NULL
        }
      )

      if (!is.null(de_markers) && nrow(de_markers) > 0) {
        if (!("gene" %in% colnames(de_markers))) {
          de_markers$gene <- rownames(de_markers)
        }
        write.csv(
          de_markers,
          file.path(annotation_helpers_dir, "cluster_markers.csv"),
          row.names = FALSE
        )

        group_marker_sets <- build_group_marker_sets(module_cols_by_group, modules_mapped, module_names_safe)
        fc_col <- if ("avg_log2FC" %in% colnames(de_markers)) {
          "avg_log2FC"
        } else if ("avg_logFC" %in% colnames(de_markers)) {
          "avg_logFC"
        } else {
          NULL
        }

        top_markers <- build_top_markers(de_markers, fc_col = fc_col, top_n = de_top_n)

        enrich_rows <- list()
        for (i in seq_len(nrow(cluster_summary))) {
          cl_id <- as.character(cluster_summary$cluster_id[i])
          grp <- as.character(cluster_summary$assigned_group_primary[i])
          top_info <- top_markers[[cl_id]]
          if (is.null(top_info)) next
          top_genes <- top_info$genes
          top_norm <- top_info$genes_norm

          overlap_genes <- character()
          target_groups <- character()
          if (!is.na(grp) && nzchar(grp) && grp %in% names(group_marker_sets)) {
            target_groups <- grp
          } else if (!is.na(cluster_summary$assigned_lineage_primary[i]) &&
                     nzchar(cluster_summary$assigned_lineage_primary[i]) &&
                     exists("lineage_groups") &&
                     cluster_summary$assigned_lineage_primary[i] %in% names(lineage_groups)) {
            target_groups <- lineage_groups[[cluster_summary$assigned_lineage_primary[i]]]
          }

          if (length(target_groups) > 0) {
            set_genes <- unique(unlist(group_marker_sets[target_groups]))
            set_genes <- set_genes[!is.na(set_genes) & nzchar(set_genes)]
            if (length(set_genes) > 0) {
              set_norm <- toupper(set_genes)
              overlap_idx <- top_norm %in% set_norm
              overlap_genes <- top_genes[overlap_idx]
            }
          }

          overlap_n <- length(overlap_genes)
          overlap_frac <- overlap_n / length(top_genes)
          enrich_rows[[length(enrich_rows) + 1]] <- data.frame(
            cluster_id = cl_id,
            assigned_group = grp,
            n_top_markers = length(top_genes),
            overlap_n = overlap_n,
            overlap_frac = overlap_frac,
            overlap_genes = paste(overlap_genes, collapse = ";"),
            stringsAsFactors = FALSE
          )
        }

        if (length(enrich_rows) > 0) {
          enrich_df <- do.call(rbind, enrich_rows)
          write.csv(
            enrich_df,
            file.path(annotation_helpers_dir, "cluster_marker_enrichment.csv"),
            row.names = FALSE
          )
        }

        # Fisher enrichment across all groups (orthogonal validation)
        fisher_rows <- list()
        universe_genes <- rownames(de_obj[[de_assay]])
        universe_norm <- unique(toupper(universe_genes))
        for (i in seq_len(nrow(cluster_summary))) {
          cl_id <- as.character(cluster_summary$cluster_id[i])
          grp <- as.character(cluster_summary$assigned_group_primary[i])
          lin <- as.character(cluster_summary$assigned_lineage_primary[i])
          top_info <- top_markers[[cl_id]]
          if (is.null(top_info)) next
          top_genes <- top_info$genes
          top_norm <- top_info$genes_norm
          top_norm <- intersect(top_norm, universe_norm)
          if (length(top_norm) == 0) next

          for (gname in names(group_marker_sets)) {
            set_genes <- group_marker_sets[[gname]]
            set_genes <- set_genes[!is.na(set_genes) & nzchar(set_genes)]
            if (length(set_genes) == 0) next
            set_norm <- unique(toupper(set_genes))
            set_norm <- intersect(set_norm, universe_norm)
            if (length(set_norm) == 0) next

            a <- length(intersect(top_norm, set_norm))
            b <- length(setdiff(top_norm, set_norm))
            c <- length(setdiff(set_norm, top_norm))
            d <- length(setdiff(universe_norm, union(top_norm, set_norm)))
            if ((a + b + c + d) <= 0) next
            ft <- tryCatch(
              fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater"),
              error = function(e) NULL
            )
            if (is.null(ft)) next
            fisher_rows[[length(fisher_rows) + 1]] <- data.frame(
              cluster_id = cl_id,
              assigned_group = grp,
              assigned_lineage = lin,
              group = gname,
              n_top_markers = length(top_norm),
              overlap_n = a,
              p_value = ft$p.value,
              odds_ratio = unname(ft$estimate),
              stringsAsFactors = FALSE
            )
          }
        }

        summary_df <- NULL
        if (length(fisher_rows) > 0) {
          fisher_df <- do.call(rbind, fisher_rows)
          fisher_df$padj <- ave(
            fisher_df$p_value,
            fisher_df$cluster_id,
            FUN = function(p) p.adjust(p, method = "BH")
          )
          if (exists("group_lineage")) {
            fisher_df$group_lineage <- group_lineage[fisher_df$group]
          } else {
            fisher_df$group_lineage <- NA_character_
          }
          write.csv(
            fisher_df,
            file.path(annotation_helpers_dir, "cluster_marker_enrichment_fisher.csv"),
            row.names = FALSE
          )

          summary_rows <- list()
          split_df <- split(fisher_df, fisher_df$cluster_id)
          for (cid in names(split_df)) {
            sub <- split_df[[cid]]
            sub <- sub[order(sub$padj, sub$p_value, -sub$overlap_n), , drop = FALSE]
            best <- sub[1, , drop = FALSE]
            assigned_group_cur <- as.character(best$assigned_group)
            assigned_lineage_cur <- as.character(best$assigned_lineage)
            best_group <- as.character(best$group)
            best_lineage <- as.character(best$group_lineage)
            match_group <- !is.na(best_group) && !is.na(assigned_group_cur) && best_group == assigned_group_cur
            match_lineage <- !is.na(best_lineage) && !is.na(assigned_lineage_cur) && best_lineage == assigned_lineage_cur
            summary_rows[[length(summary_rows) + 1]] <- data.frame(
              cluster_id = cid,
              assigned_group = assigned_group_cur,
              assigned_lineage = assigned_lineage_cur,
              best_group = best_group,
              best_group_lineage = best_lineage,
              best_overlap_n = best$overlap_n,
              best_p_value = best$p_value,
              best_padj = best$padj,
              match_assigned_group = match_group,
              match_assigned_lineage = match_lineage,
              stringsAsFactors = FALSE
            )
          }
          if (length(summary_rows) > 0) {
            summary_df <- do.call(rbind, summary_rows)
            write.csv(
              summary_df,
              file.path(annotation_helpers_dir, "cluster_marker_enrichment_summary.csv"),
              row.names = FALSE
            )
          }
        }
      }
    }

    triage_enable <- parse_bool(Sys.getenv("ANNOTATION_TRIAGE_ENABLE", "false"), FALSE)
    triage_apply <- parse_bool(Sys.getenv("ANNOTATION_TRIAGE_APPLY", "false"), FALSE)
    triage_padj <- get_env_num("ANNOTATION_TRIAGE_PADJ", 1e-4)
    triage_min_overlap <- get_env_int("ANNOTATION_TRIAGE_MIN_OVERLAP", 5, min = 1)
    triage_unknown_padj <- get_env_num("ANNOTATION_TRIAGE_UNKNOWN_PADJ", 1e-2, invalid_value = NA_real_)
    triage_unknown_min_overlap <- get_env_int(
      "ANNOTATION_TRIAGE_UNKNOWN_MIN_OVERLAP",
      2,
      min = 1,
      invalid_value = NA_integer_
    )
    triage_use_lineage <- parse_bool(Sys.getenv("ANNOTATION_TRIAGE_USE_LINEAGE", "true"), TRUE)
    if (!is.finite(triage_padj) || triage_padj <= 0) triage_padj <- 1e-4
    if (!is.finite(triage_unknown_padj) || triage_unknown_padj <= 0) triage_unknown_padj <- triage_padj
    if (!is.finite(triage_unknown_min_overlap) || triage_unknown_min_overlap < 1) {
      triage_unknown_min_overlap <- triage_min_overlap
    }

    cluster_summary$triage_applied <- FALSE
    cluster_summary$triage_best_group <- NA_character_
    cluster_summary$triage_best_lineage <- NA_character_
    cluster_summary$triage_best_padj <- NA_real_
    cluster_summary$triage_best_overlap <- NA_integer_
    cluster_summary$triage_apply <- triage_apply
    cluster_summary$triage_unknown_padj <- triage_unknown_padj
    cluster_summary$triage_unknown_min_overlap <- triage_unknown_min_overlap

    assigned_group_probable <- assigned_group_primary
    assigned_lineage_probable <- assigned_lineage_primary
    assigned_subtype_probable <- assigned_subtype_primary

    if (triage_enable && exists("summary_df") && !is.null(summary_df) && nrow(summary_df) > 0) {
      triage_map <- split(summary_df, summary_df$cluster_id)
      for (i in seq_along(cluster_ids)) {
        if (!is.na(assigned_group_primary[i]) && assigned_group_primary[i] != "Unknown") next
        cid <- cluster_ids[[i]]
        if (!(cid %in% names(triage_map))) next
        best <- triage_map[[cid]][1, , drop = FALSE]
        if (!is.finite(best$best_padj) || best$best_padj > triage_unknown_padj) next
        if (!is.finite(best$best_overlap_n) || best$best_overlap_n < triage_unknown_min_overlap) next
        best_group <- as.character(best$best_group)
        best_lineage <- as.character(best$best_group_lineage)
        if (!nzchar(best_group) && !nzchar(best_lineage)) next

        new_group <- if (triage_use_lineage && nzchar(best_lineage)) best_lineage else best_group
        assigned_group_probable[i] <- new_group
        if (nzchar(best_lineage)) {
          assigned_lineage_probable[i] <- best_lineage
        }
        if (nzchar(best_group)) {
          assigned_subtype_probable[i] <- best_group
        }
        cluster_summary$triage_applied[i] <- TRUE
        cluster_summary$triage_best_group[i] <- best_group
        cluster_summary$triage_best_lineage[i] <- best_lineage
        cluster_summary$triage_best_padj[i] <- best$best_padj
        cluster_summary$triage_best_overlap[i] <- as.integer(best$best_overlap_n)
      }
    }

    assigned_group_final <- assigned_group_primary
    assigned_lineage_final <- assigned_lineage_primary
    assigned_subtype_final <- assigned_subtype_primary
    if (triage_apply) {
      assigned_group_final <- ifelse(cluster_summary$triage_applied, assigned_group_probable, assigned_group_primary)
      assigned_lineage_final <- ifelse(cluster_summary$triage_applied, assigned_lineage_probable, assigned_lineage_primary)
      assigned_subtype_final <- ifelse(cluster_summary$triage_applied, assigned_subtype_probable, assigned_subtype_primary)
    }

    assigned_group_primary <- normalize_assigned_label(assigned_group_primary)
    assigned_lineage_primary <- normalize_assigned_label(assigned_lineage_primary)
    assigned_group_probable <- normalize_assigned_label(assigned_group_probable)
    assigned_lineage_probable <- normalize_assigned_label(assigned_lineage_probable)
    assigned_group_final <- normalize_assigned_label(assigned_group_final)
    assigned_lineage_final <- normalize_assigned_label(assigned_lineage_final)
    cluster_summary$triage_best_lineage <- normalize_assigned_label(cluster_summary$triage_best_lineage)

    cluster_summary$assigned_group_primary <- assigned_group_primary
    cluster_summary$assigned_lineage_primary <- assigned_lineage_primary
    cluster_summary$assigned_group_probable <- assigned_group_probable
    cluster_summary$assigned_lineage_probable <- assigned_lineage_probable
    cluster_summary$assigned_subtype_probable <- assigned_subtype_probable
    cluster_summary$assigned_group_final <- assigned_group_final
    cluster_summary$assigned_lineage_final <- assigned_lineage_final
    cluster_summary$assigned_subtype_final <- assigned_subtype_final

    group_map_primary <- setNames(assigned_group_primary, cluster_ids)
    group_map_probable <- setNames(assigned_group_probable, cluster_ids)
    group_map_final <- setNames(assigned_group_final, cluster_ids)
    lineage_map_primary <- setNames(assigned_lineage_primary, cluster_ids)
    lineage_map_probable <- setNames(assigned_lineage_probable, cluster_ids)
    lineage_map_final <- setNames(assigned_lineage_final, cluster_ids)
    meta$celltype_group_primary <- group_map_primary[cluster_id]
    meta$celltype_group_probable <- group_map_probable[cluster_id]
    meta$celltype_group_final <- group_map_final[cluster_id]
    meta$celltype_lineage_primary <- lineage_map_primary[cluster_id]
    meta$celltype_lineage_probable <- lineage_map_probable[cluster_id]
    meta$celltype_lineage_final <- lineage_map_final[cluster_id]
    meta$celltype_group <- meta$celltype_group_final
    merged@meta.data <- meta

    write.csv(cluster_summary, file.path(annotation_helpers_dir, "celltype_assignment_by_cluster.csv"), row.names = FALSE)
  } else {
    message("[WARN] No groups found in grouping file: ", grouping_file)
    group_names_for_plots <- character()
  }
} else {
  message("[WARN] Grouping file not found: ", grouping_file)
  group_names_for_plots <- character()
}

plot_features <- NULL
if (exists("module_keep") && length(module_keep) > 0) {
  plot_features <- c(paste0("CM_", module_keep), "EWING_1")
  plot_features <- plot_features[plot_features %in% colnames(merged@meta.data)]
}
if (is.null(plot_features) || length(plot_features) == 0) {
  warning("[Plots] No kept module features found in metadata; falling back to all score columns.")
  plot_features <- plot_features_all
}

merged <- ensure_umap(merged, npcs = npcs, seed.use = seed)
plots_dir_umap <- file.path(out_dir, "umap_pre_integration")
dir.create(plots_dir_umap, recursive = TRUE, showWarnings = FALSE)
sample_col <- if ("sample" %in% colnames(merged@meta.data)) "sample" else if ("orig.ident" %in% colnames(merged@meta.data)) "orig.ident" else NULL
if (!is.null(sample_col)) {
  ggplot2::ggsave(
    file.path(plots_dir_umap, "umap_pre_by_sample.png"),
    DimPlot(merged, group.by = sample_col) + ggplot2::ggtitle("Merged (pre-integration) by sample"),
    width = 6, height = 5
  )
}
## Intentionally omit pre-integration celltype UMAP to reduce redundant plots.

# Per-sample feature plots on merged UMAP (consistent embedding)
if (!is.null(sample_col)) {
  plots_dir_samples <- file.path(out_dir, "featureplots_samples")
  dir.create(plots_dir_samples, recursive = TRUE, showWarnings = FALSE)
  for (sample in unique(merged@meta.data[[sample_col]])) {
    sample_plot_dir <- file.path(plots_dir_samples, sample)
    keep_cells <- rownames(merged@meta.data)[merged@meta.data[[sample_col]] == sample]
    obj_s <- subset(merged, cells = keep_cells)
    save_feature_plots(obj_s, plot_features, sample_plot_dir, title_suffix = paste0(" (", sample, ")"), reduction = "umap")
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
merged_csv_cols <- unique(c(
  score_cols,
  "EWING_1",
  azimuth_cols,
  "celltype_group",
  "celltype_group_primary",
  "celltype_group_probable",
  "celltype_group_final",
  "celltype_lineage_primary",
  "celltype_lineage_probable",
  "celltype_lineage_final"
))
merged_csv_cols <- intersect(merged_csv_cols, colnames(merged@meta.data))
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
  k_anchor <- get_env_int("SEURAT_INTEGRATION_KANCHOR", 5, min = 1)
  anchors <- FindIntegrationAnchors(
    object.list = objs,
    normalization.method = "SCT",
    anchor.features = features,
    k.anchor = k_anchor
  )

  k_weight_env <- get_env_int("SEURAT_INTEGRATION_KWEIGHT", 80, min = 5)
  min_cells <- min(vapply(objs, ncol, integer(1)))
  k_weight <- min(k_weight_env, max(5, floor(min_cells / 4)))
  if (!is.finite(k_weight) || k_weight < 5) k_weight <- 5

  integrated <- NULL
  for (attempt in seq_len(4)) {
    message("[INFO] IntegrateData attempt ", attempt, " with k.weight=", k_weight, " (k.anchor=", k_anchor, ")")
    integrated <- tryCatch(
      IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = k_weight),
      error = function(e) e
    )
    if (!inherits(integrated, "error")) break
    if (grepl("k.weight", conditionMessage(integrated)) && k_weight > 5) {
      k_weight <- max(5, floor(k_weight / 2))
      next
    }
    stop(integrated)
  }

  DefaultAssay(integrated) <- "integrated"
  set.seed(seed)
  integrated <- RunPCA(integrated, npcs = npcs, verbose = FALSE)

  pca_emb_i <- tryCatch(Embeddings(integrated, "pca"), error = function(e) NULL)
  dims_use <- if (!is.null(pca_emb_i)) 1:min(npcs, ncol(pca_emb_i)) else 1:npcs
  if (length(dims_use) < 2) dims_use <- 1:2

  integrated <- tryCatch(
    RunUMAP(integrated, dims = dims_use, seed.use = seed, verbose = FALSE),
    error = function(e) RunUMAP(integrated, dims = dims_use, verbose = FALSE)
  )
  integrated <- FindNeighbors(integrated, dims = dims_use, verbose = FALSE)
  integrated <- tryCatch(
    FindClusters(integrated, resolution = resolution, random.seed = seed, verbose = FALSE),
    error = function(e) FindClusters(integrated, resolution = resolution, verbose = FALSE)
  )

  # Transfer annotation fields from merged (by cell name)
  azimuth_cols <- grep("^azimuth_", colnames(merged@meta.data), value = TRUE)
  transfer_cols <- unique(c(
    score_cols,
    "EWING_1",
    azimuth_cols,
    "celltype_group",
    "celltype_group_primary",
    "celltype_group_probable",
    "celltype_group_final",
    "celltype_lineage_primary",
    "celltype_lineage_probable",
    "celltype_lineage_final"
  ))
  for (col in transfer_cols) {
    if (col %in% colnames(merged@meta.data)) {
      integrated@meta.data[[col]] <- merged@meta.data[rownames(integrated@meta.data), col]
    }
  }
  saveRDS(integrated, file.path(out_dir, "integrated_annotated.rds"))

  plots_dir_umap_post <- file.path(out_dir, "umap_post_integration")
  dir.create(plots_dir_umap_post, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    file.path(plots_dir_umap_post, "integrated_umap_clusters.png"),
    DimPlot(integrated, group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggplot2::ggtitle("Integrated"),
    width = 6, height = 5
  )
  if (!is.null(sample_col) && sample_col %in% colnames(integrated@meta.data)) {
    ggplot2::ggsave(
      file.path(plots_dir_umap_post, "integrated_umap_by_sample.png"),
      DimPlot(integrated, group.by = sample_col) + ggplot2::ggtitle("Integrated by sample"),
      width = 6, height = 5
    )
  }
  if ("celltype_group" %in% colnames(integrated@meta.data)) {
    ggplot2::ggsave(
      file.path(plots_dir_umap_post, "integrated_umap_by_celltype.png"),
      DimPlot(integrated, group.by = "celltype_group", label = TRUE, repel = TRUE) + ggplot2::ggtitle("Integrated by cell type"),
      width = 6, height = 5
    )
  }
  plots_dir_integrated <- file.path(out_dir, "featureplots_integrated")
  save_feature_plots(integrated, plot_features, plots_dir_integrated, title_suffix = " (integrated)", reduction = NULL)
}

write_run_manifest(
  out_dir = out_dir,
  args = args,
  objs = objs,
  merged = merged,
  integrated = if (exists("integrated")) integrated else NULL
)
