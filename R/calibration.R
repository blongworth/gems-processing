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
    mutate(
      timestamp = lubridate::mdy_hms(date_time_utc_00_00, tz = "UTC"),
      seaphox_oxygen_umol_l = o2_ml_l_to_umol_l(
        oxygen_ml_l,
        p_h_temperature_celsius
      )
    ) |>
    select(
      timestamp,
      seaphox_pH = internal_p_h_p_h,
      seaphox_temp_c = p_h_temperature_celsius,
      seaphox_pressure_db = pressure_decibar,
      seaphox_salinity_psu = salinity_psu,
      seaphox_oxygen_umol_l
    ) |>
    filter(timestamp >= start_time, timestamp <= end_time)
}

#' Build period-level reference data for alternating inlet sensors
#'
#' @param df Input time series with timestamp column
#' @param value_cols Columns to retain after aggregation
#' @param inlet Which inlet period to keep
#' @param window_start Averaging window start in seconds
#' @param window_end Averaging window end in seconds
#'
#' @return Period-mean data for one inlet
#' @export
build_reference_periods <- function(
  df,
  value_cols,
  inlet = "high",
  window_start = 30,
  window_end = 420
) {
  df |>
    assign_inlets() |>
    filter_inlet_window(window_start = window_start, window_end = window_end) |>
    calculate_period_means() |>
    filter(inlet == .env$inlet) |>
    select(timestamp, mean_time, all_of(value_cols))
}

#' Add mean, high-low gradient columns for a paired concentration
#'
#' @param df Wide period-level data
#' @param high_col High inlet concentration column
#' @param low_col Low inlet concentration column
#' @param mean_name Output mean column name
#' @param gradient_name Output gradient column name
#' @param sensor_separation Vertical separation in meters
#'
#' @return Data frame with mean and gradient columns added
#' @export
add_gradient_metrics <- function(
  df,
  high_col,
  low_col,
  mean_name,
  gradient_name,
  sensor_separation = 1.02
) {
  high_col <- rlang::ensym(high_col)
  low_col <- rlang::ensym(low_col)

  df |>
    mutate(
      !!mean_name := ((!!high_col) + (!!low_col)) / 2,
      !!gradient_name := ((!!high_col) - (!!low_col)) / sensor_separation
    )
}

#' Fit oxygen data to RGA data
#'
#' @param seaphox_df Data frame with timestamp and oxygen columns
#' @param rga_df RGA data to join with
#'
#' @return A combined data frame ready for linear regression
#'
#' @export
make_oxygen_calibration_df <- function(
  rga_df,
  seaphox_df
) {
  seaphox_periods <- build_reference_periods(
    seaphox_df,
    value_cols = c("seaphox_oxygen_umol_l", "seaphox_temp_c")
  ) |>
    select(timestamp, mean_time, seaphox_oxygen_umol_l)

  rga_df |>
    calculate_period_means() |>
    filter(inlet == "high") |>
    transmute(timestamp, mean_time, mass_32_40) |>
    left_join(seaphox_periods, by = join_by(timestamp))
}

# Backward-compatible wrapper
make_ox_cal_df <- function(rga_df, seaphox_df) {
  make_oxygen_calibration_df(rga_df, seaphox_df)
}

fit_oxygen <- function(ox_cal_df) {
  lm(seaphox_oxygen_umol_l ~ mass_32_40, data = ox_cal_df)
}

#' Add oxygen data to period-level RGA data
#'
#' @param rga_df Period-level RGA data with high/low mass ratios
#' @param ox_model Linear model predicting oxygen in umol/L
#' @param sensor_separation Vertical separation between sensors in meters
#'
#' @return Period-level dataset with oxygen concentration and gradient columns
#'
#' @export
add_oxygen <- function(
  rga_df,
  ox_model,
  sensor_separation = 1.02
) {
  ox_umol_i <- coef(ox_model)[1]
  ox_umol_m <- coef(ox_model)[2]

  rga_df |>
    mutate(
      ox_high_umol_l = ox_umol_i + ox_umol_m * mass_32_40_high,
      ox_low_umol_l = ox_umol_i + ox_umol_m * mass_32_40_low
    ) |>
    add_gradient_metrics(
      high_col = ox_high_umol_l,
      low_col = ox_low_umol_l,
      mean_name = "ox_mean_umol_l",
      gradient_name = "ox_gradient_umol_l_m",
      sensor_separation = sensor_separation
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
make_co2_calibration_df <- function(
  rga_df,
  prooceanus_df,
  status_file = NULL
) {
  prooceanus_periods <- build_reference_periods(
    prooceanus_df,
    value_cols = c("prooceanus_co2_ppm", "cell_pressure")
  )

  status_temp_df <- NULL
  if (!is.null(status_file)) {
    status_temp_df <- open_dataset(status_file) |>
      select(timestamp, adv_temp = temp) |>
      collect() |>
      build_reference_periods(value_cols = "adv_temp") |>
      select(timestamp, adv_temp)
  }

  co2_cal_df <- rga_df |>
    calculate_period_means() |>
    filter(inlet == "high") |>
    transmute(timestamp, mean_time, mass_44_40) |>
    left_join(prooceanus_periods, by = join_by(timestamp))

  if (!is.null(status_temp_df)) {
    co2_cal_df <- co2_cal_df |>
      left_join(status_temp_df, by = join_by(timestamp))
  }

  co2_cal_df |>
    filter(
      !is.na(prooceanus_co2_ppm),
      is.finite(prooceanus_co2_ppm),
      prooceanus_co2_ppm >= 0,
      !is.na(cell_pressure),
      is.finite(cell_pressure),
      !is.na(adv_temp),
      is.finite(adv_temp)
    ) |>
    mutate(
      prooceanus_co2_umol_l = co2_ppm_to_umol_per_l(
        xco2_ppm = prooceanus_co2_ppm,
        temp_c = adv_temp,
        sal_psu = 31.425,
        pressure_mbar = cell_pressure
      )
    )
}

# Backward-compatible wrapper
make_co2_cal_df <- function(
  rga_df,
  prooceanus_df,
  status_file = NULL
) {
  make_co2_calibration_df(rga_df, prooceanus_df)
}

#' Fit CO2 data to RGA data
fit_co2 <- function(co2_cal_df) {
  lm(prooceanus_co2_umol_l ~ mass_44_40, data = co2_cal_df)
}

#' Add CO2 data to period-level RGA data
#'
#' @param rga_df Period-level RGA data with high/low mass ratios
#' @param co2_model Linear model predicting CO2 in umol/L
#' @param sensor_separation Vertical separation between sensors in meters
#'
#' @return Period-level dataset with CO2 concentration and gradient columns
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
      co2_low_umol_l = co2_i + co2_m * mass_44_40_low
    ) |>
    add_gradient_metrics(
      high_col = co2_high_umol_l,
      low_col = co2_low_umol_l,
      mean_name = "co2_mean_umol_l",
      gradient_name = "co2_gradient_umol_l_m",
      sensor_separation = sensor_separation
    )
}

#' Recompute oxygen gradients from period-level oxygen concentrations
#'
#' @param rga_adv_data Data frame with oxygen and temperature
#' @param sensor_separation Vertical separation between sensors in meters
#'
#' @return Data frame with oxygen concentration and gradient columns updated
#'
#' @export
calculate_oxygen_metrics <- function(rga_adv_data, sensor_separation = 1.02) {
  rga_adv_data |>
    add_gradient_metrics(
      high_col = ox_high_umol_l,
      low_col = ox_low_umol_l,
      mean_name = "ox_mean_umol_l",
      gradient_name = "ox_gradient_umol_l_m",
      sensor_separation = sensor_separation
    )
}
