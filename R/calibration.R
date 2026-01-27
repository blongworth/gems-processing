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

#' Fit oxygen data to RGA data
#'
#' @param seaphox_df Data frame with timestamp and oxygen columns
#' @param rga_df RGA data to join with
#'
#' @return A linear model fitting oxygen to RGA mass ratio
#'
#' @export
fit_oxygen <- function(
  rga_df,
  seaphox_df
) {
  seaphox_df <- seaphox_df |>
    mutate(timestamp = lubridate::floor_date(timestamp, "minute"))

  rga_df <- rga_df |>
    filter(inlet == "high") |>
    mutate(timestamp = lubridate::floor_date(timestamp, "minute")) |>
    group_by(timestamp) |>
    summarize(mass_32_40 = mean(mass_32_40, na.rm = TRUE))

  ox_cal_df <- dplyr::left_join(
    seaphox_df,
    rga_df,
    by = dplyr::join_by(timestamp)
  )

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
      oxygen_high = ox_i + ox_m * mass_32_40_high,
      oxygen_low = ox_i + ox_m * mass_32_40_low,
      ox_high_umol_l = o2_ml_l_to_umol_l(oxygen_high, adv_temp),
      ox_low_umol_l = o2_ml_l_to_umol_l(oxygen_low, adv_temp),
      ox_mean_umol_l = (ox_low_umol_l + ox_high_umol_l) / 2,
      ox_gradient_umol_l_m = (ox_high_umol_l - ox_low_umol_l) /
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
    filter(ts >= as.Date(start_time), ts <= as.Date(end_time))
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
    mutate(ts = floor_date(ts, unit = "15 mins")) |>
    group_by(ts) |>
    select(timestamp = ts, prooceanus_co2_ppm = co2, cell_pressure) |>
    summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))
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
      co2_gradient_umol_l_m = (co2_high_umol_l - co2_low_umol_l) /
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
