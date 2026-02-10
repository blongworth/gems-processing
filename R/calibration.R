# calibration for Oxygen and CO2

#' Load and process SeaPhox oxygen data
#'
#' @param file_path Path to SeaPhox CSV file
#' @param start_time Start timestamp for filtering
#' @param end_time End timestamp for filtering
#'
#' @return Data frame with timestamp and oxygen columns
#'
#' @export
load_seaphox_oxygen <- function(seaphox_path, start_time, end_time) {
  data.table::fread(seaphox_path) |>
    janitor::clean_names() |>
    mutate(timestamp = lubridate::mdy_hms(date_time_utc_00_00, tz = "UTC")) |>
    select(
      timestamp,
      seaphox_pH = internal_p_h_p_h,
      seaphox_temp_c = p_h_temperature_celsius,
      seaphox_pressure_db = pressure_decibar,
      seaphox_salinity_psu = salinity_psu,
      seaphox_oxygen_ml_l = oxygen_ml_l
    ) |>
    filter(timestamp >= start_time, timestamp <= end_time)
}

# TODO: fit to interpolated 15 minute O2
# TODO: try fitting to mean of values in each high O2 period

#' Fit oxygen data to RGA data
#'
#' @param seaphox_df Data frame with timestamp and oxygen columns
#' @param rga_df RGA data to join with
#'
#' @return A linear model fitting oxygen to RGA mass ratio
#'
#' @export
make_ox_cal_df <- function(
  rga_df,
  seaphox_df
) {
  seaphox_df <- seaphox_df |>
    assign_inlets() |>
    filter_inlet_window(window_start = 30, window_end = 420) |>
    calculate_period_means() |>
    filter(inlet == "high") |>
    select(timestamp, seaphox_oxygen_ml_l)
  # TODO: calculate umol/l here with seaphox temp

  rga_df <- rga_df |>
    calculate_period_means() |>
    filter(inlet == "high") |>
    select(timestamp, mass_32_40)

  ox_cal_df <- dplyr::left_join(
    seaphox_df,
    rga_df,
    by = dplyr::join_by(timestamp)
  )
  ox_cal_df
}

fit_oxygen <- function(ox_cal_df) {
  ox_model <- lm(seaphox_oxygen_ml_l ~ mass_32_40, data = ox_cal_df)
  return(ox_model)
}

#' Add oxygen data to RGA data
#'
#' @param seaphox_df Data frame with timestamp and oxygen columns
#' @param rga_df RGA data to join with
#'
#' @return Joined data frame ready for linear regression
#'
#' @export
add_oxygen <- function(
  rga_df,
  ox_model,
  sensor_separation = 1.02
) {
  ox_i <- coef(ox_model)[1]
  ox_m <- coef(ox_model)[2]

  rga_df |>
    mutate(
      # TODO: should umol conversion be done here or in seaphox data?
      oxygen_high = ox_i + ox_m * mass_32_40_high,
      oxygen_low = ox_i + ox_m * mass_32_40_low,
      ox_high_umol_l = o2_ml_l_to_umol_l(oxygen_high, adv_temp),
      ox_low_umol_l = o2_ml_l_to_umol_l(oxygen_low, adv_temp),
      ox_mean_umol_l = (ox_low_umol_l + ox_high_umol_l) / 2,
      ox_gradient_umol_l_m = (ox_high_umol_l - ox_low_umol_l) /
        sensor_separation
    )
}

### CO2 CALIBRATION ###

#' Load and process ProOceanus CO2 data
#'
#' @param file_path Path to ProOceanus data file
#' @param start_time Start timestamp for filtering
#' @param end_time End timestamp for filtering
#'
#' @return Data frame with timestamp and co2 columns
#'
#' @export
load_prooceanus_co2 <- function(file_path, start_time, end_time) {
  read_prooceanus(file_path) |>
    filter(ts >= as.Date(start_time), ts <= as.Date(end_time)) |>
    select(timestamp = ts, prooceanus_co2_ppm = co2, cell_pressure)
}

#' Make CO2 calibration data frame
#'
#' @param prooceanus_df Data frame with timestamp and CO2 columns
#' @param rga_df RGA data to join with
#'
#' @return A combined data frame ready for linear regression
#'
#' @export
make_co2_cal_df <- function(
  rga_df,
  prooceanus_df,
  status_file
) {
  prooceanus_df <- prooceanus_df |>
    assign_inlets() |>
    filter_inlet_window(window_start = 30, window_end = 420) |>
    calculate_period_means() |>
    filter(inlet == "high") |>
    select(timestamp, prooceanus_co2_ppm, cell_pressure)

  rga_df <- rga_df |>
    calculate_period_means() |>
    filter(inlet == "high") |>
    select(timestamp, mass_44_40)

  status_temp_df <- open_dataset(status_file) |>
    select(timestamp, adv_temp = temp) |>
    collect() |>
    assign_inlets() |>
    filter_inlet_window(window_start = 30, window_end = 420) |>
    calculate_period_means() |>
    filter(inlet == "high") |>
    select(timestamp, adv_temp)

  co2_cal_df <- dplyr::left_join(
    prooceanus_df,
    rga_df,
    by = dplyr::join_by(timestamp)
  ) |>
    left_join(status_temp_df) |>
    mutate(
      prooceanus_co2_umol_l = co2_ppm_to_umol_per_l(
        xco2_ppm = prooceanus_co2_ppm,
        temp_c = adv_temp, # status temp
        sal_psu = 31.425, # mean seaphox salinity
        pressure_mbar = cell_pressure
      )
    )

  co2_cal_df
}

#' Fit CO2 data to RGA data
fit_co2 <- function(co2_cal_df) {
  co2_model <- lm(prooceanus_co2_umol_l ~ mass_44_40, data = co2_cal_df)
  return(co2_model)
}

#' Add oxygen data to RGA data
#'
#' @param seaphox_df Data frame with timestamp and oxygen columns
#' @param rga_df RGA data to join with
#'
#' @return Joined data frame ready for linear regression
#'
#' @export
add_co2 <- function(
  rga_df,
  co2_model,
  sensor_separation = 1.02
) {
  co2_i <- coef(co2_model)[1]
  co2_m <- coef(co2_model)[2]

  rga_df |>
    mutate(
      co2_high_umol_l = co2_i + co2_m * mass_44_40_high,
      co2_low_umol_l = co2_i + co2_m * mass_44_40_low,
      co2_mean_umol_l = (co2_low_umol_l + co2_high_umol_l) / 2,
      co2_gradient_umol_l_m = (co2_high_umol_l - co2_low_umol_l) /
        sensor_separation
    )
}

#' Calculate oxygen units conversion and gradients
#'
#' @param rga_adv_data Data frame with oxygen and temperature
#' @param sensor_separation Vertical separation between sensors in meters (default: 1.02)
#'
#' @return Data frame with converted units and calculated gradients
#'
#' @export
calculate_oxygen_metrics <- function(rga_adv_data, sensor_separation = 1.02) {
  rga_adv_data |>
    mutate(
      ox_high_umol_l = o2_ml_l_to_umol_l(oxygen_high, adv_temp),
      ox_low_umol_l = o2_ml_l_to_umol_l(oxygen_low, adv_temp),
      ox_mean_umol_l = (ox_low_umol_l + ox_high_umol_l) / 2,
      ox_gradient_umol_l_m = (ox_low_umol_l - ox_high_umol_l) /
        sensor_separation
    )
}
