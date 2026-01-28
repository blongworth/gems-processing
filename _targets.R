# This is a targets workflow for processing RGA and ADV data
# Run the pipeline with: targets::tar_make()

library(targets)
library(tarchetypes)

# Set target options:
tar_option_set(
  packages = c(
    "arrow",
    "plotly",
    "gemstools",
    "suntools",
    "patchwork",
    "data.table",
    "janitor",
    "lubridate",
    "dplyr",
    "imputeTS",
    "readr",
    "quarto",
    "tidyr",
    "purrr",
    "stringr"
  ),
)

# Source helper functions
source("R/binning.R")
source("R/rga.R")
source("R/adv.R")
source("R/calibration.R")
source("R/par.R")
source("R/flux.R")
source("R/analyse.R")
source("R/eelgrass.R")

# config variables
data_dir_path <- "data/processed/lander"
flux_output_file <- file.path(data_dir_path, "gems_2025_flux.parquet")

# define the pipeline
list(
  ### input files and parameters ###

  tar_target(data_dir, data_dir_path),
  tar_target(
    moves_file,
    "data/2025_Hadleys_Harbor/gems_lander_move_times.csv",
    format = "file"
  ),
  tar_target(
    bad_times_file,
    "data/2025_Hadleys_Harbor/gems_2025_bad_data_periods.csv",
    format = "file"
  ),
  tar_target(
    eelgrass_file,
    "data/eelgrass/Naushon_eelgrass_sampling.xlsx",
    format = "file"
  ),

  ### intermediate file targets ###

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

  # minimum correlation for ADV data filtering
  tar_target(min_correlation, 50),

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

  # confirm expected output filenames are present in gems_raw
  tar_target(
    rga_file,
    {
      files <- gems_raw
      expected <- file.path(data_dir, "gems_rga.parquet")
      stopifnot(expected %in% files)
      expected
    },
    format = "file"
  ),
  tar_target(
    adv_file,
    {
      files <- gems_raw
      expected <- file.path(data_dir, "gems_adv.parquet")
      stopifnot(expected %in% files)
      expected
    },
    format = "file"
  ),
  tar_target(
    status_file,
    {
      files <- gems_raw
      expected <- file.path(data_dir, "gems_status.parquet")
      stopifnot(expected %in% files)
      expected
    },
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
    )
  ),
  tar_target(
    rga_binned_summarized,
    summarize_binned_timeseries(
      rga_binned,
      value_cols = c(
        "mass_15_40",
        "mass_28_40",
        "mass_32_40",
        "mass_44_40"
      )
    )
  ),

  # Make one row for each high-low pair and widen
  tar_target(
    rga_binned_wide,
    pair_gradient_measurements(rga_binned_summarized)
  ),
  tar_target(rga_with_par, add_par(rga_binned_wide, crds)),

  # Add temperature and remove bad times
  tar_target(
    rga_temp,
    add_status_temp(rga_with_par, status_file, bad_times)
  ),

  # Add oxygen calibration
  tar_target(
    seaphox_df_jul,
    load_seaphox_oxygen(
      seaphox_file,
      "2025-06-26 00:00:00",
      "2025-07-16 14:00:00"
    )
  ),
  tar_target(
    ox_cal_df,
    make_ox_cal_df(rga_binned, seaphox_df_jul)
  ),
  tar_target(
    ox_model,
    fit_oxygen(ox_cal_df)
  ),
  tar_target(
    rga_oxygen,
    add_oxygen(
      rga_temp,
      ox_model,
      sensor_separation = 1.02
    )
  ),

  # Add CO2 calibration
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
      rga_oxygen,
      proco2_df_jul
    )
  ),

  ### ADV data processing ###

  # Load and bin ADV data
  tar_target(
    adv_bin_rot_df,
    load_and_bin_adv(adv_file, moves_file, min_correlation)
  ),
  tar_target(rga_adv_joined, add_adv(rga_calibrated, adv_bin_rot_df)),

  ### Flux calculation ###

  # ADV data and flux calculation
  tar_target(
    adv_matlab_input,
    process_adv_to_ml_input(adv_file, moves_file, min_correlation)
  ),
  tar_target(pos_df, get_lander_positions(adv_matlab_input)),
  tar_target(
    matlab_eddyflux,
    run_matlab_eddyflux(adv_matlab_input),
    format = "file"
  ),
  tar_target(flux_dataset, process_flux_data(matlab_eddyflux, pos = pos_df)),

  # Join with Ustar and calculate flux
  tar_target(
    rga_adv_flux,
    add_grad_flux(rga_adv_joined, flux_dataset, length_scale)
  ),
  tar_target(
    hourly_flux,
    calc_hourly_flux(rga_adv_flux)
  ),

  # Write RGA+ADV flux dataset
  tar_target(
    rga_adv_flux_files,
    {
      write_parquet(hourly_flux, paste0(flux_output_file, ".parquet"))
      write_csv(hourly_flux, paste0(flux_output_file, ".csv"))
    }
  ),

  # Filtered hourly flux to remove June and October
  tar_target(
    filtered_hourly_flux,
    filter(
      hourly_flux,
      timestamp >= as.POSIXct("2025-07-01") &
        timestamp < as.POSIXct("2025-09-24")
    )
  ),

  # Calculate statistics
  tar_target(
    hourly_stats,
    calculate_hourly_statistics(filtered_hourly_flux)
  ),
  tar_target(
    monthly_stats,
    calculate_monthly_statistics(filtered_hourly_flux)
  ),

  # load eelgrass data
  tar_target(
    eelgrass,
    get_eelgrass_data(eelgrass_file)
  ),

  # Visualizations
  tar_quarto(co2_report, "reports/gems_co2_issue.qmd"),
  tar_quarto(eelgrass_report, "reports/eelgrass.qmd"),
  tar_quarto(flux_report, "reports/gems_flux_report.qmd")
)
