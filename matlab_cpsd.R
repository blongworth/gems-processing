# process GEMS ADV data into files per lander position
# for MATLAB eddyflux processing

library(tidyverse)
library(arrow)
library(furrr)
library(imputeTS)

# adv_file <- "data/processed/lander/adv_data_.parquet"
# data_dir <- "data/processed/lander/"

#' Read and process adv data to matlab flux input format
process_adv_to_ml_input <- function(
  adv_file_name,
  moves_file_name,
  min_correlation = NULL
) {
  adv_df <- read_adv_data(adv_file_name)
  proc_adv_df <- adv_df |>
    scale_adv_velocity() |>
    impute_adv_data() |>
    group_adv_data() |>
    flag_adv_lander_moves(moves_file_name)
  proc_adv_df |>
    generate_matlab_input_files()
  proc_adv_df
}

read_adv_data <- function(adv_file_path, min_correlation = NULL) {
  ds <- open_dataset(adv_file_path)
  if (!is.null(min_correlation)) {
    ds <- ds |>
      filter(
        corr1 > min_correlation,
        corr2 > min_correlation,
        corr3 > min_correlation
      )
  }
  ds |>
    select(timestamp, pressure, u, v, w, amp3, corr3) |>
    collect()
}

scale_adv_velocity <- function(adv_data, scale_factor = 10) {
  adv_data |>
    mutate(
      u = u * scale_factor,
      v = v * scale_factor,
      w = w * scale_factor
    )
}

# Impute by interpolation with a max gap of 1 s
impute_adv_data <- function(adv_data) {
  st <- min(adv_data$timestamp)
  et <- max(adv_data$timestamp)
  ts_8hz <- seq(from = st, to = et, by = 0.125)
  df_ts_8hz <- tibble(timestamp = ts_8hz)

  df_8hz_all <- df_ts_8hz |>
    left_join(adv_data)

  df_8hz_all |>
    mutate(across(-timestamp, \(x) {
      na_interpolation(x, maxgap = 8)
    }))
}


group_adv_data <- function(adv_df) {
  adv_df |>
    mutate(is_gap = is.na(u) & !is.na(lag(u, default = NA))) |>
    mutate(group = cumsum(is_gap) + 1) |>
    select(-is_gap) |>
    drop_na(pressure) |>
    group_by(group) |>
    mutate(
      # divide up big groups into 15 min blocks
      block_time = floor_date(timestamp, unit = "15min"),
      block_start = block_time != lag(block_time),
      block_start = if_else(is.na(block_start), TRUE, block_start)
    ) |>
    ungroup() |>
    mutate(block = cumsum(block_start) + 1) |>
    select(-c(block_time, block_start))
}

#Flag lander moves
flag_adv_lander_moves <- function(grouped_adv_df, moves_file_name) {
  lander_moves <- read_csv(moves_file_name) |>
    rename(change_timestamp = timestamp) |>
    mutate(lander_position = row_number())

  df_grp_mv <- grouped_adv_df |>
    arrange(timestamp) |>
    left_join(
      lander_moves,
      by = join_by(closest(timestamp >= change_timestamp))
    ) |>
    fill(lander_position, .direction = "down") |>
    mutate(lander_position = replace_na(lander_position, 1)) |>
    select(-change_timestamp)

  return(df_grp_mv)
}

# write out data
write_grouped_adv_data <- function(df_grp_mv, grouped_adv_filename) {
  df_grp_mv |>
    write_dataset(grouped_adv_filename)
}

#' Export LECS data in MATLAB eddyfluxpHO2v_41 format
#'
#' @param data A resampled, imputed LECS data_frame, grouped by gap-free runs
#' @param min_time Remove chunks shorter than `min_time` hours
#' @param data_file Output data filename
#' @param timestamp_file Output gap-free timestamp file
#'
#' @return A list containing the output and timestamp data frames
#' @export
lecs_to_ml <- function(
  data,
  min_time = .25,
  data_file = NULL,
  timestamp_file = NULL
) {
  df <- data |>
    mutate(
      time = as.numeric(difftime(timestamp, min(timestamp), units = "hours"))
    )
  df_ml <- df |>
    mutate(
      o2 = 0, # placeholder for O2
      ph = 0
    ) |> # placeholder for pH
    select(
      time,
      vx = u,
      vy = v,
      vz = w,
      o2,
      ph,
      pres = pressure,
      SNR = amp3,
      corr = corr3
    ) |>
    # NA's crash Matlab. These should be taken care of better and earlier
    drop_na()

  df_ml_ts <- df %>%
    group_by(block) %>%
    summarise(start = min(time), end = max(time)) %>%
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
      write_delim(timestamp_file, col_names = FALSE)
  }
  list(df_ml, df_ml_ts)
}

# create a file for a given lander position
group_ml <- function(pos, proc_adv_df) {
  data_path <- "data/processed/matlab_input"
  data_file <- file.path(data_path, paste0("lecs_ml_data_", pos, ".dat"))
  timestamp_file <- file.path(
    data_path,
    paste0("lecs_ml_timestamps_", pos, ".dat")
  )

  proc_adv_df |>
    filter(lander_position == pos) |>
    arrange(timestamp) |>
    collect() |>
    lecs_to_ml(
      min_time = 0.24,
      data_file = data_file,
      timestamp_file = timestamp_file
    )
}

get_lander_positions <- function(proc_adv_df) {
  proc_adv_df |>
    group_by(lander_position) |>
    summarize(t0 = min(timestamp), n = n()) |>
    arrange(lander_position) |>
    collect()
}

# generate files for all lander positions
generate_matlab_input_files <- function(proc_adv_df) {
  positions <- get_lander_positions(proc_adv_df)$lander_position
  purrr::map(positions, \(x) group_ml(x, proc_adv_df))
  # future_map(positions, \(x) group_ml(x, proc_adv_df))
}

copy_matlab_input_files <- function(
  source_dir = "data/processed/matlab_input",
  dest_dir = "matlab"
) {
  files <- list.files(source_dir, full.names = TRUE)
  file.copy(files, dest_dir, overwrite = TRUE)
}

# Run all fluxes using matlab/eddyflux_batch.m scripted using run_eddyflux.sh
run_matlab_eddyflux <- function() {
  copy_matlab_input_files()
  setwd("matlab")
  result <- system2(
    command = "bash",
    args = c(
      "./run_eddyflux.sh"
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  setwd("..")

  if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
    stop("MATLAB processing failed with status ", attr(result, "status"))
  }

  # return file paths of output files
  read_flux_files()
}

# get t0 for all files and convert from decimal hours to UTC timestamps
# combine flux data

read_flux_files <- function(
  matlab_dir = "matlab",
  pattern = "outfile2Rot_lecs_ml*"
) {
  list.files(
    matlab_dir,
    pattern = pattern,
    full.names = TRUE
  )
}

process_flux_data <- function(
  matlab_dir = "matlab",
  pattern = "outfile2Rot_lecs_ml*",
  pos
) {
  read_delim(read_flux_files(matlab_dir, pattern), id = "file") |>
    mutate(
      lander_position = as.integer(str_extract(file, "(\\d{1,2})(?=\\.dat$)"))
    ) |>
    left_join(select(pos, lander_position, t0)) |>
    group_by(lander_position) |>
    mutate(
      timestamp = as.POSIXct(timemean * 3600, origin = t0, tz = "GMT"),
      phi2 = na_if(phi2, NaN),
      Cd = na_if(Cd, NaN)
    ) |>
    select(-c(file, timemean, ...38, t0)) |>
    arrange(timestamp) |>
    ungroup()
}

write_flux_file <- function(flux_df, data_dir = "data/processed/lander/") {
  write_dataset(flux_df, file.path(data_dir, "gems_flux_8hz_rot2.parquet"))
}
