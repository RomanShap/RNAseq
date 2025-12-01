#!/usr/bin/env Rscript
# AS_main.R
# Main orchestrator for the alternative splicing (AS) pipeline.
#
# This script reads a YAML configuration file, sets up paths and logging,
# and sequentially sources numbered step scripts located in
# `scripts/Alternative_splicing_scripts`.  It expects the current working
# directory to be the project root where `config/`, `logs/` and
# `output/` reside.  The pipeline root is inferred from the location of this
# file.

## Load required packages ------------------------------------------------------
suppressMessages({
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required but not installed")
  }
})

## Determine configuration file ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript AS_main.R /path/to/config.yml")
}
config_path <- normalizePath(args[[1]], mustWork = TRUE)

## Determine project root and pipeline root -----------------------------------
# The project root is the current working directory.  All outputs will be
# written under this directory.
project_root <- normalizePath(getwd(), mustWork = TRUE)

# Determine the directory of this script by parsing the --file argument
get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 0) {
    return(NULL)
  }
  return(normalizePath(sub("^--file=", "", file_arg)))
}
this_file <- get_script_path()
if (!is.null(this_file)) {
  pipeline_root <- normalizePath(dirname(this_file), mustWork = TRUE)
} else {
  # Fallback: assume AS_main.R resides in the working directory
  pipeline_root <- normalizePath(file.path(project_root, "AS_pipeline_v2"), mustWork = FALSE)
}

## Create required directories -------------------------------------------------
dirs_to_create <- c("cluster", "config", "logs", "output")
for (d in dirs_to_create) {
  dir_path <- file.path(project_root, d)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message(sprintf("Created directory: %s", dir_path))
  }
}

## Read configuration ----------------------------------------------------------
cfg <- yaml::read_yaml(config_path)
message(sprintf("Loaded configuration from %s", config_path))

## Set up global environment ---------------------------------------------------
# Make the configuration accessible to all sourced scripts via a global
# variable.  The scripts will reference `cfg` directly.
assign("cfg", cfg, envir = .GlobalEnv)

# Define project and pipeline roots globally as well
assign("PROJECT_ROOT", project_root, envir = .GlobalEnv)
assign("PIPELINE_ROOT", pipeline_root, envir = .GlobalEnv)

## Open log file ---------------------------------------------------------------
log_file <- file.path(project_root, "logs", sprintf("AS_pipeline_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE, split = TRUE)
message("Alternative Splicing Pipeline v2 started")
message(sprintf("Project root: %s", project_root))
message(sprintf("Pipeline root: %s", pipeline_root))

## Source helper utilities -----------------------------------------------------
# Load common utilities from the shared repository.  This file must exist
# relative to the pipeline root (e.g. ../shared/99_utils.R).  You may need
# to adjust this path depending on your repository structure.
shared_utils_path <- file.path(dirname(pipeline_root), "shared", "99_utils.R")
if (file.exists(shared_utils_path)) {
  source(shared_utils_path)
  message(sprintf("Sourced shared utilities from %s", shared_utils_path))
} else {
  warning(sprintf("Could not find shared/99_utils.R at %s", shared_utils_path))
}

# Load splicing‑specific utilities
splicing_utils_path <- file.path(pipeline_root, "scripts", "Alternative_splicing_scripts", "99_utils_as.R")
if (file.exists(splicing_utils_path)) {
  source(splicing_utils_path)
  message(sprintf("Sourced splicing utilities from %s", splicing_utils_path))
} else {
  stop(sprintf("Required file not found: %s", splicing_utils_path))
}

## Source and run step scripts -------------------------------------------------
# Identify all step scripts in the Alternative_splicing_scripts folder.  Files
# starting with two digits followed by an underscore are considered steps and
# are executed in lexicographical order.  The utility file 99_utils_as.R is
# skipped.
steps_dir <- file.path(pipeline_root, "scripts", "Alternative_splicing_scripts")
step_files <- list.files(steps_dir, pattern = "^[0-9]{2}_.+\\.R$", full.names = TRUE)
step_files <- sort(step_files)

for (step_file in step_files) {
  step_name <- basename(step_file)
  message(sprintf("\n--- Running step: %s ---", step_name))
  source(step_file)
  message(sprintf("Completed step: %s", step_name))
}

message("\nAll steps completed successfully")
sink(type = "message")
sink(type = "output")
close(log_con)