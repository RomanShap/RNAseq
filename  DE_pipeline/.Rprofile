# 1. Make conda R library visible
conda_lib <- "/hdd/home/roman/miniconda3/envs/bioenv/lib/R/library"
if (dir.exists(conda_lib)) {
  .libPaths(c(conda_lib, .libPaths()))
}

source("renv/activate.R")