# Processing functions for RGA and ADV workflow
library(tidyverse)
library(arrow)
library(suntools)
library(gemstools)

load_and_widen_rga <- function(rga_file) {
  open_dataset(rga_file) |>
    collect() |>
    process_rga_to_wide()
}

remove_bad_rga_periods <- function(rga_data, bad_times) {
  rga_data |>
    anti_join(
      bad_times,
      by = join_by(timestamp >= start_time, timestamp <= end_time)
    )
}

normalize_rga_by_argon <- function(rga_data) {
  rga_data |>
    mutate(
      across(starts_with("mass_"), ~ .x / mass_40, .names = "{.col}_40")
    )
}

pair_gradient_measurements <- function(rga_binned) {
  rga_binned |>
    group_by(grp = cumsum(inlet == "high")) |>
    mutate(mean_timestamp = mean(bin_time)) |>
    ungroup() |>
    pivot_wider(
      id_cols = c(grp, mean_timestamp),
      names_from = inlet,
      values_from = starts_with("mass_"),
      names_sep = "_"
    ) |>
    select(-grp) |>
    mutate(
      timestamp = lubridate::round_date(mean_timestamp, unit = "15 minutes"),
      across(
        ends_with("_high"),
        \(x) (x + get(sub("_high$", "_low", cur_column()))) / 2,
        .names = "{.col}_mean"
      ),
      across(
        ends_with("_high"),
        \(x) x - get(sub("_high$", "_low", cur_column())),
        .names = "{.col}_diff"
      ),
    ) |>
    dplyr::relocate(timestamp, .before = dplyr::everything())
}

#' Calculate PAR from solar position
#'
#' @param crds Coordinate matrix (longitude, latitude)
#' @param timestamp Timestamp vector
#'
#' @return Numeric vector of PAR values
#'
#' @details
#' Calculates predicted PAR from solar altitude angle
#'
#' @export
calculate_par <- function(crds, timestamp) {
  altitude <- solarpos(crds, timestamp)[, 2]
  altitude <- replace(altitude, altitude <= 0, 0)
  insolation <- sin(altitude * pi / 180)
  par <- 0.45 * insolation * 4.6 * 1000
  par
}

#' Add PAR to binned RGA data
#'' @param rga_binned Binned RGA data frame
#' @param crds Coordinate matrix (longitude, latitude)
#' @return Data frame with PAR column added
#' @export
add_par <- function(rga_df, crds) {
  rga_df |>
    mutate(par = calculate_par(crds, timestamp))
}

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
  seaphox_df,
  sensor_separation = 1.02
) {
  seaphox_15min <- seaphox_df |>
    mutate(timestamp = lubridate::floor_date(timestamp, "15 minute")) |>
    dplyr::group_by(timestamp) |>
    summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

  ox_cal_df <- dplyr::left_join(
    seaphox_15min,
    rga_df,
    by = dplyr::join_by(timestamp)
  )

  ox_model <- lm(seaphox_oxygen_ml_l ~ mass_32_40_high_mean, data = ox_cal_df)
  ox_i <- coef(ox_model)[1]
  ox_m <- coef(ox_model)[2]

  rga_df |>
    mutate(
      oxygen_high = ox_i + ox_m * mass_32_40_high,
      oxygen_low = ox_i + ox_m * mass_32_40_low,
      ox_high_umol_l = o2_ml_l_to_umol_l(oxygen_high, adv_temp),
      ox_low_umol_l = o2_ml_l_to_umol_l(oxygen_low, adv_temp),
      ox_mean_umol_l = (ox_low_umol_l + ox_high_umol_l) / 2,
      ox_gradient_umol_l_m = (ox_low_umol_l - ox_high_umol_l) /
        sensor_separation
    )
}

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
    mutate(ts = floor_date(ts, unit = "15 mins")) |>
    group_by(ts) |>
    summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))
}


#' Add CO2 data to RGA data
#'
#' @param co2_raw Raw CO2 data from ProOceanus
#' @param rga_data RGA data with oxygen
#'
#' @return Joined data frame ready for linear regression
#'
#' @export
add_co2 <- function(
  rga_df,
  co2_raw,
  sensor_separation = 1.02
) {
  co2_df <- co2_raw |>
    select(timestamp = ts, prooceanus_co2_ppm = co2, cell_pressure)
  co2_cal_df <- dplyr::inner_join(
    co2_df,
    rga_df,
    by = dplyr::join_by(timestamp)
  ) |>
    mutate(
      prooceanus_co2_umol_l = co2_ppm_to_umol_per_l(
        xco2_ppm = prooceanus_co2_ppm,
        temp_c = adv_temp, # status temp
        sal_psu = 31.425, # mean seaphox salinity
        pressure_mbar = cell_pressure
      )
    )

  co2_model <- lm(
    prooceanus_co2_umol_l ~ mass_44_40_high_mean,
    data = co2_cal_df
  )
  co2_i <- coef(co2_model)[1]
  co2_m <- coef(co2_model)[2]

  rga_df |>
    mutate(
      co2_high_umol_l = co2_i + co2_m * mass_44_40_high,
      co2_low_umol_l = co2_i + co2_m * mass_44_40_low,
      co2_mean_umol_l = (co2_low_umol_l + co2_high_umol_l) / 2,
      co2_gradient_umol_l_m = (co2_low_umol_l - co2_high_umol_l) /
        sensor_separation
    )
}


#' Apply CO2 calibration coefficients to RGA data
#'
#' @param rga_data RGA data with raw mass ratios
#' @param intercept Calibration intercept
#' @param slope Calibration slope
#'
#' @return Data frame with co2_high and co2_low columns
#'
#' @export
apply_co2_calibration <- function(rga_data, intercept, slope) {
  rga_data |>
    mutate(
      co2_high = intercept + slope * mass_44_40_high,
      co2_low = intercept + slope * mass_44_40_low
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


#' Calculate flux from concentration gradient and velocity
#'
#' @param data Data frame with oxygen gradient and Ustar values
#' @param length_scale Turbulent length scale parameter (default: 0.3)
#'
#' @return Data frame with ox_flux column added
#'
#' @export
calculate_grad_flux <- function(data, grad_var, length_scale = 0.3) {
  data |>
    mutate(
      lscale = length_scale,
      ox_flux = -1 * Ustar * lscale * {{ grad_var }}
    )
}

get_ustar <- function(flux_dataset) {
  flux_dataset |>
    select(timestamp, Ustar) |>
    mutate(
      timestamp = lubridate::floor_date(timestamp, unit = "15 minutes")
    ) |>
    group_by(timestamp) |>
    summarize(Ustar = mean(Ustar, na.rm = TRUE))
}

# Join with Ustar and calculate flux
add_grad_flux <- function(rga_adv_processed, flux_dataset, length_scale) {
  ustar_data <- get_ustar(flux_dataset)
  rga_adv_processed |>
    left_join(ustar_data, by = join_by(timestamp)) |>
    mutate(
      lscale = length_scale,
      ox_flux = -1 * Ustar * lscale * ox_gradient_umol_l_m,
      co2_flux = -1 * Ustar * lscale * co2_gradient_umol_l_m
    )
}

#' Calculate hourly statistics for oxygen metrics and flux
#'
#' @param flux_data Data frame with flux and concentration data
#' @param coordinates Coordinate matrix for solar noon calculation
#'
#' @return Data frame with hourly statistics grouped by solar hour
#'
#' @export
calculate_hourly_statistics <- function(flux_data, coordinates) {
  flux_data |>
    mutate(
      hour = hour(timestamp),
      diff_noon = timestamp - 16800,
      solar_hour = hour(diff_noon)
    ) |>
    group_by(solar_hour) |>
    summarize(
      ox_mean_umol_l_hourly = mean(ox_mean_umol_l, na.rm = TRUE),
      ox_mean_umol_l_hourly_sd = sd(ox_mean_umol_l, na.rm = TRUE),
      ox_gradient_umol_l_m_hourly = mean(ox_gradient_umol_l_m, na.rm = TRUE),
      ox_gradient_umol_l_m_hourly_sd = sd(ox_gradient_umol_l_m, na.rm = TRUE),
      ox_flux_hourly = mean(ox_flux, na.rm = TRUE),
      ox_flux_hourly_sd = sd(ox_flux, na.rm = TRUE),
      .groups = "drop"
    )
}

# add status temperature
add_status_temp <- function(rga_df, status_file, bad_times) {
  status_temp_df <- open_dataset(status_file) |>
    select(timestamp, adv_temp = temp) |>
    collect() |>
    mutate(
      timestamp = lubridate::floor_date(timestamp, unit = "15 minutes")
    ) |>
    group_by(timestamp) |>
    summarize(adv_temp = mean(adv_temp, na.rm = TRUE))

  rga_df |>
    left_join(status_temp_df, by = join_by(timestamp)) |>
    anti_join(
      bad_times,
      by = join_by(timestamp >= start_time, timestamp <= end_time)
    )
}
