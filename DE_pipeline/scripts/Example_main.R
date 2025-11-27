#!/usr/bin/env Rscript

## Quantification_main.R
## Entry point for Kallisto-based quantification + downstream analysis.
## Uses config/config.yml and the KallistoBootstrap layout you showed.

rm(list = ls())
graphics.off()

## ----------------------------------------------------------------------
## 1) Collect input arguments
## ----------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
# project_root <- normalizePath(args[1], mustWork = TRUE)

## ----------------------------------------------------------------------
## 2) Set working directory and check project structure
## ----------------------------------------------------------------------

project_root <- "/hdd/home/roman/FETseq/Differential_expression"
message("Setting project root to: ", project_root)
pipeline_root <- "/hdd/home/roman/GIT_FET/Differential_expression_pipeline"
message("Pipeline root: ", pipeline_root)

# 1) Check that project_root exists and is a directory
if (!dir.exists(project_root)) {
  stop("Project root does not exist: ", project_root,
       "\nPlease create it or correct the path.")
}
if (!dir.exists(pipeline_root)) {
  stop("Project root does not exist: ", pipeline_root,
       "\nPlease create it or correct the path.")
}

# 2) Check required subdirectories inside project_root
required_root_dirs <- c("cluster", "config", "logs", "output")
required_pipeline_dirs <- c("resources", "scripts")

dir_project_paths <- file.path(project_root, required_root_dirs)
missing_project_dirs   <- required_root_dirs[!dir.exists(dir_project_paths)]

dir_pipeline_paths <- file.path(pipeline_root, required_pipeline_dirs)
missing_pipeline_dirs   <- required_pipeline_dirs[!dir.exists(dir_pipeline_paths)]

if (length(missing_project_dirs) > 0L) {
  # Create the missing directories one by one
  for (d in missing_project_dirs) {
    dir_to_create <- file.path(project_root, d)
    dir.create(dir_to_create, recursive = TRUE, showWarnings = FALSE)
  }

  stop(
    "Project root '", project_root,
    "' was missing required subdirectories. ",
    "The following directories were created: ",
    paste(missing_project_dirs, collapse = ", "),
    "\nPlease check them (especially 'config' and 'resources') and rerun."
  )
}

if (length(missing_pipeline_dirs) > 0L) {
  # Create the missing directories one by one
  for (d in missing_pipeline_dirs) {
    dir_to_create <- file.path(pipeline_root, d)
    dir.create(dir_to_create, recursive = TRUE, showWarnings = FALSE)
  }

  stop(
    "Project root '", pipeline_root,
    "' was missing required subdirectories. ",
    "The following directories were created: ",
    paste(missing_pipeline_dirs, collapse = ", "),
    "\nPlease check them (especially 'config' and 'resources') and rerun."
  )
}

# If we reach here, everything is OK
setwd(project_root)
message("Project structure OK. Working directory set to: ", getwd())


## ----------------------------------------------------------------------
## 3) Load required libraries
## ----------------------------------------------------------------------

## Activate renv from project root
# source(file.path(project_root, "renv", "activate.R"))

suppressPackageStartupMessages({
  library(yaml)
  library(DESeq2)
  library(dplyr)
  library(org.Hs.eg.db)
  library(writexl)
  library(fgsea)
  library(msigdbr)
  library(readr)
  library(tibble)
  library(purrr)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(biomaRt)
  # library(UpSetR)
  library(ComplexUpset)
  library(tidyr)
  library(openxlsx)
  library(ggplot2)
  library(scales)
  library(grid)
  library(grDevices)
  library(VennDiagram)
  library(ComplexHeatmap)
  library(circlize)
  library(forcats)
})

## ----------------------------------------------------------------------
## 4) Read config
## ----------------------------------------------------------------------

config_path <- if (length(args) >= 1) args[1] else "config/config.yml"

if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path)
}

message("Using config: ", config_path)
cfg <- yaml::read_yaml(config_path)

## Small helper: make path absolute if it is relative
make_path <- function(root, ..., path = NULL) {
  # Combine ... and path into one relative string
  parts <- c(list(...), path)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  rel   <- do.call(file.path, parts)
  
  if (is.null(rel) || rel == "") return(NULL)
  
  # Absolute if starts with '/' or 'C:'-style drive letter
  if (grepl("^(/|[A-Za-z]:)", rel)) {
    normalizePath(rel, mustWork = FALSE)
  } else {
    normalizePath(file.path(root, rel), mustWork = FALSE)
  }
}


## ----------------------------------------------------------------------
## 2) Resolve project paths
## ----------------------------------------------------------------------

# Kallisto base directory (your KallistoBootstrap)
kallisto_dir <- normalizePath(cfg$kallisto_dir, mustWork = TRUE)
# message("Kallisto output directory : ", kallisto_dir)

# Output directory (relative to project root unless absolute in YAML)
output_dir <- file.path(project_root, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# message("Output dir   : ", output_dir)

# Resources / caches
tx2gene_rds   <- make_path(pipeline_root, cfg$rds_dir, cfg$cache$tx2gene)
msigdb_c5_rds <- make_path(pipeline_root, cfg$rds_dir, cfg$msigdb_c5)
# message("tx2gene RDS  : ", tx2gene_rds)
# message("MSigDB C5 RDS: ", msigdb_c5_rds)

## ----------------------------------------------------------------------
## 3) Global analysis parameters from YAML
## ----------------------------------------------------------------------

DE_LFC   <- cfg$de$lfc
DE_PADJ  <- cfg$de$padj
DE_QVAL  <- cfg$de$q_value

GO_ONTOLOGY         <- cfg$go$ontology
GO_SIMPLIFY_CUTOFF  <- cfg$go$simplify_cutoff
GO_COUNT_MIN        <- cfg$go$count_min
GO_COUNT_MAX        <- cfg$go$count_max

COLOR_KO <- cfg$colors  
GO_TOP_N <- if (!is.null(cfg$go$top_n)) cfg$go$top_n else 5L

FGSEA_NPERM        <- cfg$fgsea$nperm
FGSEA_TOP_N        <- cfg$fgsea$top_n
FGSEA_PADJ_MAX     <- cfg$fgsea$padj_max
FGSEA_ABS_NES_MIN  <- cfg$fgsea$abs_nes_min
FGSEA_MAX_PATHWAYS <- cfg$fgsea$max_pathways

## ----------------------------------------------------------------------
## 4) Samples and mapping to Kallisto folders
## ----------------------------------------------------------------------

## From YAML:
## sets:
##   EWSR1: [EWSR1KO_1, EWSR1KO_2, EWSR1KO_3]
##   FUS: [FUSKO_1, FUSKO_2, FUSKO_3]
##   TAF15: [TAF15KO_1, TAF15KO_2, TAF15KO_3]
##   WT: [WT_1, WT_2, WT_3]
##
## And the folders look like:
##   EWSR1KO_1_ZKRN250022434-1A_... / abundance.tsv
##   EWSR1KO_2_ZKRN250022435-1A_... / abundance.tsv
##   FUSKO_1_ZKRN250022438-1A_... / abundance.tsv
##   FUSKO_2_ZKRN250022439-1A_... / abundance.tsv
##   TAF15KO_1_ZKRN250022440-1A_... / abundance.tsv
##   TAF15KO_2_ZKRN250022441-1A_... / abundance.tsv
##   WT_1_ZKRN250022442-1A_... / abundance.tsv
##   WT_2_ZKRN250022443-1A_... / abundance.tsv
##   etc.
##
## We:
##   1) derive sample names from cfg$sets
##   2) for each sample (e.g. "WT_1") find the unique directory whose name
##      starts with "WT_1_"
##   3) build a named vector of abundance.tsv paths

sets <- cfg$sets          # list: names = conditions, values = character vectors
conditions <- names(sets)
if (is.null(conditions) || length(conditions) == 0) {
  stop("No 'sets' defined in config.yml.")
}

samples <- unlist(sets, use.names = FALSE)
sample_conditions <- rep(conditions, lengths(sets))
sample_info <- data.frame(
  sample    = samples,
  condition = sample_conditions,
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$sample

# List Kallisto subdirectories
all_dirs <- list.dirs(kallisto_dir, full.names = FALSE, recursive = FALSE)
if (length(all_dirs) == 0) {
  stop("No subdirectories found in Kallisto dir: ", kallisto_dir)
}

sample_dirs <- setNames(character(length(samples)), samples)

for (s in samples) {
  prefix <- paste0(s, "_")
  matches <- all_dirs[startsWith(all_dirs, prefix)]
  if (length(matches) == 0L) {
    stop("No Kallisto folder found for sample '", s,
         "' (expected something like '", s, "_XXXX').")
  }
  if (length(matches) > 1L) {
    stop("Multiple Kallisto folders match sample '", s, "': ",
         paste(matches, collapse = ", "))
  }
  sample_dirs[s] <- matches
}

kallisto_paths <- file.path(kallisto_dir, sample_dirs, "abundance.tsv")

names(kallisto_paths) <- names(sample_dirs)

# Simple check that all files exist
missing_tsv <- kallisto_paths[!file.exists(kallisto_paths)]
if (length(missing_tsv) > 0L) {
  stop("These Kallisto abundance.tsv files are missing:\n",
       paste(missing_tsv, collapse = "\n"))
}

# Store sample info in a data frame for downstream use (DESeq2 colData)
SAMPLE_INFO <- data.frame(
  sample    = samples,
  condition = sample_conditions,
  row.names = samples,
  stringsAsFactors = FALSE
)

## These objects are now global and can be used by the helper scripts:
##   kallisto_paths : named vector sample -> abundance.tsv path
##   SAMPLE_INFO    : sample / condition table
##   output_dir     : where to write all results
##   tx2gene_rds    : cached tx2gene if you use it
##   msigdb_c5_rds  : cached MSigDB if you use it
##   DE_LFC, DE_PADJ, GO_* : thresholds
##   COLOR_KO       : FET colors

## ----------------------------------------------------------------------
## 5) Start R log file
## ----------------------------------------------------------------------

log_dir <- file.path(project_root, "logs")

log_file <- file.path(
  log_dir,
  paste0("Quantification_R_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)
message("R log file   : ", log_file)

zz <- file(log_file, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

message("Kallisto output directory : ", kallisto_dir)
message("Output dir   : ", output_dir)
message("tx2gene RDS  : ", tx2gene_rds)
message("MSigDB C5 RDS: ", msigdb_c5_rds)
cat("\n")  
message("Kallisto paths for each of the samples:")
for (s in names(kallisto_paths)) {
  message("  ", s, " -> ", kallisto_paths[[s]])
}
cat("\n")
## ----------------------------------------------------------------------
## 6) Source helper scripts and run pipeline
## ----------------------------------------------------------------------

## Split the huge script into:
##   scripts/Quantification_R/01_counts_from_kallisto.R
##   scripts/Quantification_R/02_deseq.R
##   scripts/Quantification_R/03_go_enrichment.R
##   ...
## Each helper can assume the globals above exist.

helper_dir <- file.path(pipeline_root, "scripts", "Differential_expression_scripts")
helper_files <- list.files(helper_dir, pattern = "\\.[Rr]$", full.names = TRUE)
helper_files <- sort(helper_files)

if (length(helper_files) == 0) {
  stop("No helper scripts found in ", helper_dir)
  cat("\n")
}

## Ensure 99_utils.R is sourced first so helpers are available everywhere
utils_idx  <- grepl("99_utils\\.R$", basename(helper_files))
utils_file <- helper_files[utils_idx]
other_files <- sort(helper_files[!utils_idx])

helper_files_ordered <- c(utils_file, other_files)

for (f in helper_files_ordered) {
  message("Sourcing: ", f)
  cat("\n")
  source(f, chdir = FALSE)
}

message("Quantification pipeline completed successfully.")
# Close log file
sink(type = "message")
sink(type = "output")
close(zz)
