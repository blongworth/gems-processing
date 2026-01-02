# This is a targets workflow for processing RGA and ADV data
# Run the pipeline with: targets::tar_make()

library(targets)
library(tarchetypes)

# Set target options:
tar_option_set(
  packages = c(
    "tidyverse",
    "arrow",
    "plotly",
    "gemstools",
    "suntools",
    "patchwork",
    "data.table",
    "janitor",
    "lubridate",
    "dplyr",
    "imputeTS"
  ),
)

# Source helper functions
source("bin_timeseries.R")
source("functions.R")
source("matlab_cpsd.R")

# config variables
data_dir_path <- "data/processed/lander"
flux_output_file <- file.path(data_dir_path, "gems_flux.parquet")

# Define the pipeline
list(
  ### input files and parameters ###

  tar_target(data_dir, data_dir_path),
  tar_target(
    moves_file_name,
    "data/2025_Hadleys_Harbor/gems_lander_move_times.csv",
    format = "file"
  ),
  tar_target(
    bad_times_file,
    "data/2025_Hadleys_Harbor/gems_2025_bad_data_periods.csv",
    format = "file"
  ),

  ### intermediate file targets ###

  # combined, cleaned RGA and ADV data
  tar_target(
    rga_file,
    "data/processed/lander/gems_rga.parquet",
    format = "file"
  ),
  tar_target(
    raw_adv_path,
    file.path(data_dir, "gems_adv.parquet"),
    format = "file"
  ),
  tar_target(
    status_file,
    file.path(data_dir, "gems_status.parquet"),
    format = "file"
  ),

  tar_target(
    adv_file,
    file.path(data_dir, "adv_bin_rot.parquet"),
    format = "file"
  ),
  tar_target(
    bad_times,
    read_csv(bad_times_file)
  ),
  tar_target(
    seaphox_file,
    "data/2025_Hadleys_Harbor/seaphox/Shallow SeapHox2-0000453-Data-20250716T200119.csv",
    format = "file"
  ),
  tar_target(
    prooceanus_file,
    "data/2025_Hadleys_Harbor/prooceanus/prooceanus_co2_gems_hadley_2025-08-19.txt",
    format = "file"
  ),
  # location for calculating PAR
  tar_target(crds, matrix(c(-70.7003, 41.51875), nrow = 1)),

  # length scale in meters for flux calculation
  tar_target(length_scale, 0.3),

  ### Initial processing from raw files ###

  # Uses functions from gemstools package
  tar_target(
    gems_raw,
    gems_process_data(
      file_dir = "./data/2025_Hadleys_Harbor/GEMS_LANDER",
      out_dir = data_dir,
      dedupe = FALSE
    ),
    format = "file"
  ),

  ### RGA data processing and calibration ###

  # Load and process RGA data
  tar_target(rga_wide, load_and_widen_rga(rga_file)),
  tar_target(rga_clean, remove_bad_rga_periods(rga_wide, bad_times)),
  tar_target(rga_normalized, normalize_rga_by_argon(rga_clean)),

  # Bin timeseries data
  tar_target(
    rga_binned,
    bin_timeseries(
      rga_normalized,
      datetime_col = "timestamp",
      value_cols = c(
        "mass_15_40",
        "mass_28_40",
        "mass_32_40",
        "mass_44_40"
      )
    )
  ),

  # Make one row for each high-low pair
  tar_target(rga_binned_wide, pair_gradient_measurements(rga_binned)),

  tar_target(rga_with_par, add_par(rga_binned_wide, crds)),

  # Add temperature and remove bad times
  tar_target(
    rga_temp,
    add_status_temp(rga_with_par, status_file, bad_times)
  ),

  # Add oxygen calibration
  tar_target(
    rga_with_oxygen,
    add_oxygen(
      rga_temp,
      seaphox_file,
      "2025-07-12 00:00:00",
      "2025-07-16 14:00:00"
    )
  ),

  tar_target(
    seaphox_df_jul,
    load_seaphox_oxygen(
      seaphox_file,
      "2025-07-12 00:00:00",
      "2025-07-16 14:00:00"
    )
  ),

  tar_target(
    proco2_df_jul,
    load_prooceanus_co2(
      prooceanus_file,
      "2025-07-20 00:00:00",
      "2025-07-29 00:00:00"
    )
  ),

  # Add calibrated CO2 to RGA data
  tar_target(
    rga_calibrated,
    add_co2(
      rga_with_oxygen,
      proco2_df_jul
    )
  ),

  ### ADV data processing ###

  # Load and bin ADV data
  tar_target(
    adv_raw,
    open_dataset(adv_file) |> collect()
  ),

  # Bin velocity data to 15 min bins
  tar_target(
    adv_binned,
    adv_raw |>
      group_by(grp = cumsum(inlet == "high")) |>
      summarize(
        mean_timestamp = mean(bin_time),
        timestamp = lubridate::round_date(mean_timestamp, unit = "15 minutes"),
        across(c(pressure, cur_speed, cur_dir, u_rot, v_rot, w_rot), \(x) {
          mean(x)
        })
      ) |>
      dplyr::relocate(timestamp, .before = dplyr::everything())
  ),

  # Select ADV columns
  tar_target(
    adv_select,
    adv_binned |>
      select(
        timestamp,
        pressure,
        cur_speed,
        cur_dir,
        u = u_rot,
        v = v_rot,
        w = w_rot
      )
  ),

  # Join RGA with ADV
  tar_target(
    rga_adv_joined,
    rga_calibrated |>
      left_join(adv_select, by = join_by(timestamp)) |>
      drop_na()
  ),

  ### Flux calculation ###

  # ADV data and flux calculation
  tar_target(
    adv_matlab_input,
    process_adv_to_ml_input(raw_adv_path, moves_file_name)
  ),

  tar_target(pos_df, get_lander_positions(adv_matlab_input)),

  tar_target(matlab_eddyflux, run_matlab_eddyflux(), format = "file"),

  tar_target(flux_dataset, process_flux_data(pos = pos_df)),

  # Join with Ustar and calculate flux
  tar_target(
    rga_adv_flux,
    add_grad_flux(rga_adv_joined, flux_dataset, length_scale)
  ),

  # Write RGA+ADV flux dataset
  tar_target(
    rga_adv_flux_file,
    write_parquet(rga_adv_flux, flux_output_file),
  ),

  # Calculate hourly statistics
  tar_target(
    hourly_stats,
    calculate_hourly_statistics(rga_adv_flux, crds)
  )
)
