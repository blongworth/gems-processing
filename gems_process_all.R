#!/usr/local/bin/Rscript
# read and parse GEMS data to parquet

library(gemstools)
library(tictoc)

out_dir <- "./data/processed/lander/"

tic()
gems_process_data(
  file_dir = "./data/2025_Hadleys_Harbor/GEMS_LANDER",
  out_dir = out_dir,
  dedupe = FALSE
)
toc()
