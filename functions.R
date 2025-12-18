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
load_seaphox_oxygen <- function(file_path, start_time, end_time) {
  data.table::fread(file_path) |>
    janitor::clean_names() |>
    mutate(timestamp = lubridate::mdy_hms(date_time_utc_00_00, tz = "UTC")) |>
    select(
      timestamp,
      pH = internal_p_h_p_h,
      temp = p_h_temperature_celsius,
      pressure = pressure_decibar,
      sal = salinity_psu,
      oxygen = oxygen_ml_l
    ) |>
    filter(timestamp >= start_time, timestamp <= end_time)
}


#' Aggregate calibration data to 15-minute intervals
#'
#' @param calib_data Data frame with timestamp and oxygen columns
#' @param rga_data RGA data to join with
#'
#' @return Joined data frame ready for linear regression
#'
#' @export
prepare_oxygen_calibration <- function(calib_data, rga_data) {
  seaphox_15min <- calib_data |>
    mutate(timestamp = lubridate::floor_date(timestamp, "15 minute")) |>
    dplyr::group_by(timestamp) |>
    dplyr::summarise(seaphox_oxy = mean(oxygen))

  dplyr::inner_join(seaphox_15min, rga_data, by = dplyr::join_by(timestamp))
}


#' Apply oxygen calibration coefficients to RGA data
#'
#' @param rga_data RGA data with raw mass ratios
#' @param intercept Calibration intercept
#' @param slope Calibration slope
#'
#' @return Data frame with oxygen_high and oxygen_low columns
#'
#' @export
apply_oxygen_calibration <- function(rga_data, intercept, slope) {
  rga_data |>
    mutate(
      oxygen_high = intercept + slope * mass_32_40_high,
      oxygen_low = intercept + slope * mass_32_40_low
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
    summarize(co2 = mean(co2, na.rm = TRUE))
}


#' Prepare CO2 calibration data
#'
#' @param co2_raw Raw CO2 data from ProOceanus
#' @param rga_data RGA data with oxygen
#'
#' @return Joined data frame ready for linear regression
#'
#' @export
prepare_co2_calibration <- function(co2_raw, rga_data) {
  co2_df <- co2_raw |>
    rename(timestamp = ts, prooceanus_co2 = co2)

  dplyr::inner_join(co2_df, rga_data, by = dplyr::join_by(timestamp))
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
      ox_high_umol_l = o2_ml_l_to_umol_l(oxygen_high, temp),
      ox_low_umol_l = o2_ml_l_to_umol_l(oxygen_low, temp),
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
calculate_oxygen_flux <- function(data, length_scale = 0.3) {
  data |>
    mutate(
      lscale = length_scale,
      ox_flux = -1 * Ustar * lscale * (ox_gradient_umol_l_m)
    )
}

get_ustar <- function(flux_file) {
  open_dataset(file.path(flux_file)) |>
    collect() |>
    select(timestamp, Ustar) |>
    mutate(
      timestamp = lubridate::floor_date(timestamp, unit = "15 minutes")
    ) |>
    group_by(timestamp) |>
    summarize(Ustar = mean(Ustar, na.rm = TRUE))
}


# Join with Ustar and calculate flux
calculate_grad_flux <- function(rga_adv_processed, flux_file) {
  ustar_data <- get_ustar(flux_file)
  rga_adv_processed |>
    left_join(ustar_data, by = join_by(timestamp)) |>
    calculate_oxygen_flux()
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
