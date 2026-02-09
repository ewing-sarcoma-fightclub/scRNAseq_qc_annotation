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
