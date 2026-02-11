#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Force installs to the active R environment library, not user/site libs.
Sys.unsetenv(c("R_LIBS", "R_LIBS_USER"))
.libPaths(.Library)

# On macOS, normalize compiler vars to current R toolchain config.
if (identical(Sys.info()[["sysname"]], "Darwin")) {
  r_bin <- file.path(R.home("bin"), "R")
  r_cfg <- function(key) {
    out <- tryCatch(system2(r_bin, c("CMD", "config", key), stdout = TRUE, stderr = FALSE), error = function(e) "")
    trimws(paste(out, collapse = " "))
  }
  cfg <- c(
    CC = r_cfg("CC"),
    CXX = r_cfg("CXX"),
    CXX11 = r_cfg("CXX11"),
    CXX14 = r_cfg("CXX14"),
    CXX17 = r_cfg("CXX17"),
    CXX20 = r_cfg("CXX20")
  )
  cfg <- cfg[nzchar(cfg)]
  if (length(cfg) > 0) do.call(Sys.setenv, as.list(cfg))
}

usage <- function() {
  cat(
    "Usage: Rscript r/install_pipeline_packages.R [options]\n",
    "\n",
    "Options:\n",
    "  --skip-github                Skip GitHub package installs (DoubletFinder, DropletQC)\n",
    "  --install-azimuth-refs       Install Azimuth SeuratData references\n",
    "  --azimuth-refs <csv>         Comma-separated refs (default: pbmcref,bonemarrowref,lungref,adiposeref,fetusref,liverref)\n",
    "  --ncpus <int>                Number of CPUs for package install (default: 1)\n",
    "  -h, --help                   Show this help\n",
    sep = ""
  )
}

args <- commandArgs(trailingOnly = TRUE)
skip_github <- FALSE
install_azimuth_refs <- FALSE
azimuth_refs <- c("pbmcref", "bonemarrowref", "lungref", "adiposeref", "fetusref", "liverref")
ncpus <- 1L

i <- 1L
while (i <= length(args)) {
  arg <- args[[i]]
  if (arg == "--skip-github") {
    skip_github <- TRUE
    i <- i + 1L
  } else if (arg == "--install-azimuth-refs") {
    install_azimuth_refs <- TRUE
    i <- i + 1L
  } else if (arg == "--azimuth-refs") {
    if (i + 1L > length(args)) stop("--azimuth-refs requires a value", call. = FALSE)
    azimuth_refs <- strsplit(args[[i + 1L]], ",", fixed = TRUE)[[1]]
    azimuth_refs <- trimws(azimuth_refs)
    azimuth_refs <- azimuth_refs[nzchar(azimuth_refs)]
    i <- i + 2L
  } else if (arg == "--ncpus") {
    if (i + 1L > length(args)) stop("--ncpus requires a value", call. = FALSE)
    ncpus <- as.integer(args[[i + 1L]])
    if (!is.finite(ncpus) || ncpus < 1) stop("--ncpus must be a positive integer", call. = FALSE)
    i <- i + 2L
  } else if (arg %in% c("-h", "--help")) {
    usage()
    quit(save = "no", status = 0)
  } else {
    stop(paste("Unknown option:", arg), call. = FALSE)
  }
}

is_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

install_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, is_installed, logical(1))]
  if (length(missing) == 0) return(invisible(NULL))
  message("[install] CRAN: ", paste(missing, collapse = ", "))
  install.packages(missing, dependencies = TRUE, Ncpus = ncpus, lib = .Library)
}

install_bioc <- function(pkgs) {
  install_cran("BiocManager")
  missing <- pkgs[!vapply(pkgs, is_installed, logical(1))]
  if (length(missing) == 0) return(invisible(NULL))
  message("[install] Bioconductor: ", paste(missing, collapse = ", "))
  BiocManager::install(
    missing,
    ask = FALSE,
    update = FALSE,
    Ncpus = ncpus,
    lib = .Library,
    dependencies = TRUE,
    force = TRUE
  )
}

install_github <- function(repo, pkg) {
  if (is_installed(pkg)) return(invisible(NULL))
  install_cran("remotes")
  message("[install] GitHub: ", repo)
  remotes::install_github(repo, upgrade = "never", dependencies = TRUE, lib = .Library)
}

message("[install] Starting pipeline R package installation")

# 1) CRAN baseline packages used directly by pipeline scripts or required by Bioc deps.
install_cran(c("future", "mclust", "remotes", "statmod"))

# 2) Bioconductor packages used directly in pipeline scripts.
install_bioc(c("AUCell", "DropletUtils"))

# 3) Extra CRAN stack needed by DropletQC GitHub install on some systems.
if (!skip_github) {
  install_cran(c("nloptr", "lme4", "pbkrtest", "car", "rstatix", "ggpubr"))
}

# 4) GitHub-only packages.
if (!skip_github) {
  install_github("chris-mcginnis-ucsf/DoubletFinder", "DoubletFinder")
  install_github("powellgenomicslab/DropletQC", "DropletQC")
}

# 5) Optional: Azimuth references via SeuratData.
if (install_azimuth_refs) {
  if (!is_installed("SeuratData")) {
    stop("SeuratData not installed. Install conda package r-seurat-data first.", call. = FALSE)
  }
  for (ref in azimuth_refs) {
    message("[install] Azimuth reference: ", ref)
    SeuratData::InstallData(ref, force = FALSE)
  }
}

required <- c(
  "Seurat", "readxl", "Azimuth", "SeuratData", "SoupX",
  "DropletUtils", "Matrix", "HGNChelper", "future",
  "mclust", "AUCell", "ggplot2"
)
if (!skip_github) {
  required <- c(required, "DoubletFinder", "DropletQC")
}

versions <- vapply(
  required,
  function(pkg) if (is_installed(pkg)) as.character(utils::packageVersion(pkg)) else "MISSING",
  character(1)
)

status_df <- data.frame(package = required, version = versions, row.names = NULL, stringsAsFactors = FALSE)
message("[install] Package status:")
print(status_df, row.names = FALSE)

missing_required <- status_df$package[status_df$version == "MISSING"]
if (length(missing_required) > 0) {
  stop("Missing required R packages: ", paste(missing_required, collapse = ", "), call. = FALSE)
}

message("[install] R package installation complete.")
