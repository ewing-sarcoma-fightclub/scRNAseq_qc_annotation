#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(utils)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: build_celltype_grouping.R <rules_csv> [out_file]", call. = FALSE)
}

rules_csv <- args[[1]]
out_file <- ifelse(length(args) >= 2, args[[2]], "")

get_script_dir <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  return(getwd())
}

if (!nzchar(out_file)) {
  out_file <- file.path(get_script_dir(), "..", "resources", "celltype_grouping.txt")
}

if (!file.exists(rules_csv)) {
  stop("Rules CSV not found: ", rules_csv, call. = FALSE)
}

rules <- read.csv(rules_csv, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("group", "lineage", "azimuth_refs", "azimuth_label_patterns", "cellmarker_name_patterns", "notes")
missing <- setdiff(required_cols, colnames(rules))
if (length(missing) > 0) {
  stop("Rules CSV missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

con <- file(out_file, open = "w")
on.exit(close(con), add = TRUE)

writeLines(c(
  "# Cell type grouping rules (generated)",
  "#",
  "# Sources used to define Azimuth label strings:",
  "# - azimuth_label_catalog.txt (generated from installed references)",
  "#",
  "# Format per group:",
  "# group: <GROUP_ID>",
  "# azimuth_refs: comma-separated references",
  "# azimuth_label_patterns: comma-separated substrings (case-insensitive)",
  "# cellmarker_name_patterns: comma-separated substrings (case-insensitive)",
  "# lineage: coarse lineage label (used for final assignment)",
  "# notes: free text",
  ""
), con)

for (i in seq_len(nrow(rules))) {
  row <- rules[i, , drop = FALSE]
  writeLines(paste0("group: ", row$group), con)
  writeLines(paste0("azimuth_refs: ", row$azimuth_refs), con)
  writeLines(paste0("azimuth_label_patterns: ", row$azimuth_label_patterns), con)
  writeLines(paste0("cellmarker_name_patterns: ", row$cellmarker_name_patterns), con)
  writeLines(paste0("lineage: ", row$lineage), con)
  writeLines(paste0("notes: ", row$notes), con)
  writeLines("", con)
}

message("Wrote: ", out_file)
