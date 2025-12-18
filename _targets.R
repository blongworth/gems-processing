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
    "dplyr"
  ),
)

# Source helper functions
source("bin_timeseries.R")
source("functions.R")


# Define the pipeline
list(
  # input files and parameters
  tar_target(data_dir, "data/processed/lander/"),
  tar_target(
    rga_file,
    "data/processed/landergems_rga/rga.parquet",
    format = "file"
  ),
  tar_target(
    adv_file,
    file.path(data_dir, "adv_bin_rot.parquet"),
    format = "file"
  ),
  tar_target(
    flux_file,
    file.path(data_dir, "gems_flux_8hz_rot2.parquet"),
    format = "file"
  ),
  tar_target(
    bad_times_file,
    "data/2025_Hadleys_Harbor/gems_2025_bad_data_periods.csv",
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
  tar_target(
    status_file,
    file.path(data_dir, "gems_status.parquet"),
    format = "file"
  ),
  tar_target(crds, matrix(c(-70.7003, 41.51875), nrow = 1)),

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

  # Load and process SeaPhox data (July)
  tar_target(
    seaphox_jul,
    load_seaphox_oxygen(
      seaphox_file,
      "2025-07-12 00:00:00",
      "2025-07-16 14:00:00"
    )
  ),

  # Aggregate SeaPhox to 15 min and join with RGA
  tar_target(
    oxygen_calib_data,
    prepare_oxygen_calibration(seaphox_jul, rga_with_par)
  ),

  # Fit oxygen model
  tar_target(
    oxygen_model,
    lm(seaphox_oxy ~ mass_32_40_high_mean, data = oxygen_calib_data)
  ),

  # Add calibrated oxygen to RGA data
  tar_target(
    rga_with_oxygen,
    {
      ox_i <- coef(oxygen_model)[1]
      ox_m <- coef(oxygen_model)[2]
      apply_oxygen_calibration(rga_with_par, ox_i, ox_m)
    }
  ),

  # Load and process ProOceanus CO2 data
  tar_target(
    co2_calib_raw,
    load_prooceanus_co2(
      prooceanus_file,
      "2025-07-17 20:50:00",
      "2025-08-03"
    )
  ),

  # Join CO2 with RGA for calibration
  tar_target(
    co2_calib_data,
    prepare_co2_calibration(co2_calib_raw, rga_with_oxygen)
  ),

  # Fit CO2 model
  tar_target(
    co2_model,
    lm(prooceanus_co2 ~ mass_44_40_high_mean, data = co2_calib_data)
  ),

  # Add calibrated CO2 to RGA data
  tar_target(
    rga_calibrated,
    {
      co2_i <- coef(co2_model)[1]
      co2_m <- coef(co2_model)[2]
      apply_co2_calibration(rga_with_oxygen, co2_i, co2_m)
    }
  ),

  # Write RGA calibrated data
  tar_target(
    rga_cal_file,
    write_parquet(
      rga_calibrated,
      "data/processed/lander/rga_calibrated.parquet"
    )
  ),

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

  # Load and aggregate status temperature
  tar_target(
    status_temp,
    {
      open_dataset(status_file) |>
        select(timestamp, temp) |>
        collect() |>
        mutate(
          timestamp = lubridate::floor_date(timestamp, unit = "15 minutes")
        ) |>
        group_by(timestamp) |>
        summarize(temp = mean(temp, na.rm = TRUE))
    }
  ),

  # Add temperature and remove bad times
  tar_target(
    rga_adv_complete,
    rga_adv_joined |>
      left_join(status_temp, by = join_by(timestamp)) |>
      anti_join(
        bad_times,
        by = join_by(timestamp >= start_time, timestamp <= end_time)
      )
  ),

  # Convert oxygen to umol/l and calculate mean and gradient
  tar_target(
    rga_adv_processed,
    calculate_oxygen_metrics(rga_adv_complete)
  ),

  # Join with Ustar and calculate flux
  tar_target(
    rga_adv_flux,
    calculate_grad_flux(rga_adv_processed, flux_file)
  ),

  # Write RGA+ADV flux dataset
  tar_target(
    rga_adv_flux_file,
    {
      write_parquet(
        rga_adv_flux,
        "data/processed/lander/rga_adv_oxygen_flux.parquet"
      )
      "data/processed/lander/rga_adv_oxygen_flux.parquet"
    },
    format = "file"
  ),

  # Calculate hourly statistics
  tar_target(
    hourly_stats,
    calculate_hourly_statistics(rga_adv_flux, crds)
  )
)
