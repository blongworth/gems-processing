# process GEMS ADV data into files per lander position
# for MATLAB eddyflux processing

library(tidyverse)
library(arrow)
library(furrr)

adv_file <- "data/processed/lander/adv_data_.parquet"
data_dir <- "data/processed/lander/"

df <- open_dataset(adv_file)
names(df)

df_grp <- df |> 
  collect() |> 
  mutate(is_gap = is.na(u) & !is.na(lag(u, default = NA))) |>
  mutate(group = cumsum(is_gap) + 1) |>
  select(-is_gap) |>
  drop_na(pressure)

# divide up big groups
# 15 min blocks within groups

df_grp <- df_grp  |>
  group_by(group) |>
  mutate(block_time = floor_date(timestamp, unit = "15min"),
         block_start = block_time != lag(block_time),
         block_start = if_else(is.na(block_start), TRUE, block_start)) |>
  ungroup() |>
  mutate(block = cumsum(block_start) + 1) |>
  select(-c(block_time, block_start))

#Flag lander moves
lander_moves <- read_csv("data/2025_Hadleys_Harbor/gems_lander_move_times.csv") |> 
  rename(change_timestamp = timestamp) |> 
  mutate(lander_position = row_number())

df_grp_mv <- df_grp |> 
  arrange(timestamp) |> 
  collect() |> 
  left_join(lander_moves,
            by = join_by(closest(timestamp >= change_timestamp))) |> 
  fill(lander_position, .direction = "down") |> 
  mutate(lander_position = replace_na(lander_position, 1)) |> 
  select(-change_timestamp)

# write out data
df_grp_mv |> 
  write_dataset(file.path(data_dir, "gems_adv_grp_15_mv.parquet"))

ds <- open_dataset(file.path(data_dir, "gems_adv_grp_15_mv.parquet"))

#' Export LECS data in MATLAB eddyfluxpHO2v_41 format
#'
#' @param data A resampled, imputed LECS data_frame, grouped by gap-free runs
#' @param min_time Remove chunks shorter than `min_time` hours
#' @param data_file Output data filename
#' @param timestamp_file Output gap-free timestamp file
#'
#' @return A list containing the output and timestamp data frames
#' @export
lecs_to_ml <- function(data,
                       min_time = .25, 
                       data_file = NULL, 
                       timestamp_file = NULL) {
  df <- data |>
    mutate(time = as.numeric(difftime(timestamp, min(timestamp),
                                      units = "hours")))
  df_ml <- df |>
    mutate(o2 = 0,  # placeholder for O2
           ph = 0) |>  # placeholder for pH
    select(time, vx = u, vy = v, vz = w,
           o2, ph, pres = pressure,
           SNR = amp3, corr = corr3) |> 
    # NA's crash Matlab. These should be taken care of better and earlier
    drop_na()
  
  df_ml_ts <- df %>%
    group_by(block) %>%
    summarise(start = min(time),
              end = max(time)) %>%
    ungroup()
  
  # remove last row of times to make shift indexing work
  df_ml_ts <- df_ml_ts[-nrow(df_ml_ts), ]
  
  if (min_time > 0) {
    df_ml_ts <- df_ml_ts |>
      filter(end - start >= min_time)
  }
  if (!is.null(data_file)) {
    write_delim(df_ml, data_file)
  }
  if (!is.null(timestamp_file)) {
    df_ml_ts %>%
      select(start, end) %>%
      write_delim(timestamp_file,
                  col_names = FALSE)
  }
  list(df_ml, df_ml_ts)
}

# create a file for a given lander position
group_ml <- function(pos, ...) {
  data_path = "data/processed/matlab_input"
  data_file <- file.path(data_path, paste0("lecs_ml_data_", pos, ".dat"))
  timestamp_file <- file.path(data_path, paste0("lecs_ml_timestamps_", pos, ".dat"))
  
  ds |> 
    filter(lander_position == pos) |> 
    arrange(timestamp) |> 
    collect() |> 
    lecs_to_ml(min_time = 0.24, data_file = data_file, timestamp_file = timestamp_file)
}

pos <- ds |> 
  group_by(lander_position) |> 
  summarize(t0 = min(timestamp),
            n = n()) |> 
  arrange(lander_position) |> 
  collect()

# generate files for all lander positions
positions <- pos$lander_position
future_map(positions, group_ml)

# Run all fluxes using matlab/eddyflux_batch.m scripted using run_eddyflux.sh
# system("cd matlab;
#         cp ../data/processed/matlab_input/* .; 
#         ./run_eddyflux.sh")

# get t0 for all files and convert from decimal hours to UTC timestamps
# combine flux data

flux_files <- list.files("matlab", pattern = "outfile2Rot_lecs_ml*", full.names = TRUE)
flux <- read_delim(flux_files, id = "file") |> 
  mutate(lander_position = as.integer(str_extract(file, "(\\d{1,2})(?=\\.dat$)"))) |> 
  left_join(select(pos, lander_position, t0)) |> 
  group_by(lander_position) |> 
  mutate(timestamp = as.POSIXct(timemean * 3600, origin = t0, tz = "GMT"),
         phi2 = na_if(phi2, NaN),
         Cd = na_if(Cd, NaN)) |> 
  select(-c(file, timemean, ...38, t0)) |> 
  arrange(timestamp) |> 
  ungroup()

write_dataset(flux, file.path(data_dir, "lecs_flux_4hz_rot2_fl5_filt.parquet"))