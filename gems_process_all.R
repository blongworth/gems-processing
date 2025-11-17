#!/usr/local/bin/Rscript
# read and parse GEMS data to parquet

library(gemstools)
library(tictoc)

out_dir <- "./data/processed/lander/"

# ADV data
# Extract S and D lines
# perform time alignment (check that works at 8Hz)
# rotate z plane to minimize mean z for each lander position

# RGA data
# Extract R lines 
# pivot longer to average time for masses
# normalize using prooceanus and seaphox data

# Valve data
# Align valve timestamps to RGA/ADV timebase
# Add "high/low" boolean

# pre-flux dataset
# extract high mean, low mean and mean velocity for each 15 minute timestep
# means from 00:02 to 00:07 of each 7.5 min block
tic()
gems_process_data(
  file_dir = "./data/2025_Hadleys_Harbor/GEMS_LANDER",
  out_dir = out_dir,
  dedupe = FALSE
)
toc()
