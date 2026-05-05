# analysis functions

#' Calculate hourly statistics for oxygen metrics and flux
#'
#' @param flux_data Data frame with flux and concentration data
#' @param coordinates Coordinate matrix for solar noon calculation
#'
#' @return Data frame with hourly statistics grouped by solar hour
#'
#' @export
calculate_hourly_statistics <- function(flux_data) {
  flux_data |>
    mutate(
      diff_noon = timestamp - 16800,
      solar_hour = hour(diff_noon)
    ) |>
    select(-timestamp, -diff_noon) |>
    group_by(solar_hour) |>
    summarise(
      across(
        everything(),
        c(
          mean = \(x) mean(x, na.rm = TRUE),
          sd = \(x) sd(x, na.rm = TRUE),
          se = \(x) sd(x, na.rm = TRUE) / sqrt(length(x))
        )
      ),
      .groups = "drop"
    )
}

calculate_monthly_statistics <- function(flux_data) {
  flux_data |>
    mutate(
      diff_noon = timestamp - 16800,
      solar_hour = hour(diff_noon),
      month = month(timestamp)
    ) |>
    select(-timestamp, -diff_noon) |>
    group_by(solar_hour, month) |>
    summarise(
      across(
        everything(),
        c(
          mean = \(x) mean(x, na.rm = TRUE),
          sd = \(x) sd(x, na.rm = TRUE),
          se = \(x) sd(x, na.rm = TRUE) / sqrt(length(x))
        )
      ),
      .groups = "drop"
    )
}

calculate_seasonal_statistics <- function(flux_data, coordinates) {
  flux_data |>
    mutate(
      diff_noon = timestamp - 16800,
      solar_hour = hour(diff_noon),
      season = case_when(
        month %in% 10:12 ~ "Fall",
        month %in% 1:3 ~ "Winter",
        month %in% 4:6 ~ "Spring",
        TRUE ~ "Summer"
      )
    ) |>
    select(-timestamp, -diff_noon) |>
    group_by(solar_hour, season) |>
    summarise(
      across(
        everything(),
        c(
          mean = \(x) mean(x, na.rm = TRUE),
          sd = \(x) sd(x, na.rm = TRUE),
          se = \(x) sd(x, na.rm = TRUE) / sqrt(length(x))
        )
      ),
      .groups = "drop"
    )
}

#' @title Calculate daily Net Ecosystem Metabolism (NEM)
calculate_nem_daily <- function(hourly_flux) {
  hourly_flux |>
    mutate(day = as.Date(timestamp)) |>
    group_by(day) |>
    summarise(
      n_hours = n(),
      n_missing = sum(is.na(ox_flux)),
      nem_mmol_m2_h = if_else(
        all(is.na(ox_flux)),
        NA_real_,
        sum(ox_flux, na.rm = TRUE) / 24
      ),
      .groups = "drop"
    )
}

#' @title Calculate Monthly Net Ecosystem Metabolism (NEM)
calculate_nem_monthly <- function(monthly_stats) {
  monthly_stats |>
    group_by(month) |>
    summarise(
      n_hours = n(),
      n_missing = sum(is.na(ox_flux_mean)),
      nem_mmol_m2_day = if_else(
        all(is.na(ox_flux_mean)),
        NA_real_,
        sum(ox_flux_mean, na.rm = TRUE)
      ),
      .groups = "drop"
    )
}

#' Build 15-minute SeapHOx chemistry periods for carbonate calculations
#'
#' @param seaphox_df Data frame with SeapHOx data
#' @param start_time Start time for filtering
#' @param end_time End time for filtering
#'
#' @return Data frame with 15-minute SeapHOx chemistry periods
#' @export
build_seaphox_carb_periods <- function(
  seaphox_df,
  start_time = as.POSIXct("2025-06-26", tz = "UTC"),
  end_time = as.POSIXct("2025-07-16 14:00:00", tz = "UTC")
) {
  seaphox_df |>
    filter(timestamp >= start_time, timestamp <= end_time) |>
    mutate(timestamp = floor_date(timestamp, unit = "15 min")) |>
    group_by(timestamp) |>
    summarise(
      seaphox_temp_c = mean(seaphox_temp_c, na.rm = TRUE),
      seaphox_salinity_psu = mean(seaphox_salinity_psu, na.rm = TRUE),
      seaphox_pressure_db = mean(seaphox_pressure_db, na.rm = TRUE),
      seaphox_pH = mean(seaphox_pH, na.rm = TRUE),
      .groups = "drop"
    )
}

umol_l_to_mol_kg <- function(x_umol_l, density_kg_m3) {
  x_umol_l / density_kg_m3 / 1e3
}

mol_kg_to_umol_l <- function(x_mol_kg, density_kg_m3) {
  x_mol_kg * density_kg_m3 * 1e3
}

calculate_dic_umol_l <- function(
  ph,
  co2_umol_l,
  salinity_psu,
  temp_c,
  pressure_db = 0
) {
  result <- rep(NA_real_, length(ph))
  valid <- is.finite(ph) &
    is.finite(co2_umol_l) &
    is.finite(salinity_psu) &
    is.finite(temp_c) &
    is.finite(pressure_db)

  if (!any(valid)) {
    return(result)
  }

  pressure_bar <- pressure_db[valid] / 10
  density_kg_m3 <- seacarb::rho(
    S = salinity_psu[valid],
    T = temp_c[valid],
    P = pressure_bar
  )

  co2_mol_kg <- umol_l_to_mol_kg(co2_umol_l[valid], density_kg_m3)
  dic_mol_kg <- seacarb::carb(
    flag = 1,
    var1 = ph[valid],
    var2 = co2_mol_kg,
    S = salinity_psu[valid],
    T = temp_c[valid],
    P = pressure_bar
  )$DIC

  result[valid] <- mol_kg_to_umol_l(dic_mol_kg, density_kg_m3)
  result
}

calculate_alkalinity_umol_l <- function(
  ph,
  co2_umol_l,
  salinity_psu,
  temp_c,
  pressure_db = 0
) {
  result <- rep(NA_real_, length(ph))
  valid <- is.finite(ph) &
    is.finite(co2_umol_l) &
    is.finite(salinity_psu) &
    is.finite(temp_c) &
    is.finite(pressure_db)

  if (!any(valid)) {
    return(result)
  }

  pressure_bar <- pressure_db[valid] / 10
  density_kg_m3 <- seacarb::rho(
    S = salinity_psu[valid],
    T = temp_c[valid],
    P = pressure_bar
  )

  co2_mol_kg <- umol_l_to_mol_kg(co2_umol_l[valid], density_kg_m3)
  alk_mol_kg <- seacarb::carb(
    flag = 1,
    var1 = ph[valid],
    var2 = co2_mol_kg,
    S = salinity_psu[valid],
    T = temp_c[valid],
    P = pressure_bar
  )$ALK

  result[valid] <- mol_kg_to_umol_l(alk_mol_kg, density_kg_m3)
  result
}

calculate_dic_from_alkalinity <- function(
  co2_umol_l,
  alkalinity_umol_l,
  salinity_psu,
  temp_c,
  pressure_db = 0
) {
  result <- rep(NA_real_, length(co2_umol_l))
  valid <- is.finite(co2_umol_l) &
    is.finite(alkalinity_umol_l) &
    is.finite(salinity_psu) &
    is.finite(temp_c) &
    is.finite(pressure_db)

  if (!any(valid)) {
    return(result)
  }

  pressure_bar <- pressure_db[valid] / 10
  density_kg_m3 <- seacarb::rho(
    S = salinity_psu[valid],
    T = temp_c[valid],
    P = pressure_bar
  )

  co2_mol_kg <- umol_l_to_mol_kg(co2_umol_l[valid], density_kg_m3)
  alk_mol_kg <- umol_l_to_mol_kg(alkalinity_umol_l[valid], density_kg_m3)
  dic_mol_kg <- seacarb::carb(
    flag = 4,
    var1 = co2_mol_kg,
    var2 = alk_mol_kg,
    S = salinity_psu[valid],
    T = temp_c[valid],
    P = pressure_bar
  )$DIC

  result[valid] <- mol_kg_to_umol_l(dic_mol_kg, density_kg_m3)
  result
}

add_dic <- function(rga_calibrated, seaphox_carb_periods, sensor_separation) {
  rga_calibrated |>
    left_join(seaphox_carb_periods, by = "timestamp") |>
    mutate(
      dic_high_umol_l = calculate_dic_umol_l(
        ph = seaphox_pH,
        co2_umol_l = co2_high_umol_l,
        salinity_psu = seaphox_salinity_psu,
        temp_c = seaphox_temp_c,
        pressure_db = seaphox_pressure_db
      ),
      dic_low_umol_l = calculate_dic_umol_l(
        ph = seaphox_pH,
        co2_umol_l = co2_low_umol_l,
        salinity_psu = seaphox_salinity_psu,
        temp_c = seaphox_temp_c,
        pressure_db = seaphox_pressure_db
      ),
      alk_high_umol_l = calculate_alkalinity_umol_l(
        ph = seaphox_pH,
        co2_umol_l = co2_high_umol_l,
        salinity_psu = seaphox_salinity_psu,
        temp_c = seaphox_temp_c,
        pressure_db = seaphox_pressure_db
      ),
      dic_low_predicted_umol_l = calculate_dic_from_alkalinity(
        co2_umol_l = co2_low_umol_l,
        alkalinity_umol_l = alk_high_umol_l,
        salinity_psu = seaphox_salinity_psu,
        temp_c = seaphox_temp_c,
        pressure_db = seaphox_pressure_db
      )
    ) |>
    add_gradient_metrics(
      high_col = dic_high_umol_l,
      low_col = dic_low_umol_l,
      mean_name = "dic_mean_umol_l",
      gradient_name = "dic_gradient_umol_l_m",
      sensor_separation = sensor_separation
    ) |>
    mutate(
      dic_gradient_predicted_umol_l_m = (dic_high_umol_l - dic_low_predicted_umol_l) /
        sensor_separation
    )
}

add_scalar_flux <- function(data, grad_var, flux_name, length_scale) {
  von_karman <- 0.41

  data |>
    mutate(
      lscale = length_scale,
      !!rlang::sym(flux_name) := -1 *
        Ustar *
        von_karman *
        lscale *
        {{ grad_var }} *
        3600
    )
}

#' Create dataset with DIC and DIC flux from SeapHOx and calibrated CO2
#'
#' @param rga_calibrated Data frame with calibrated gas concentrations
#' @param seaphox_df Data frame with SeapHOx data
#' @param flux_dataset Data frame with ADV-derived friction velocity
#' @param rga_adv_flux Data frame with oxygen and CO2 fluxes
#' @param length_scale Turbulent length scale in meters
#' @param sensor_separation Vertical sensor separation in meters
#'
#' @return Data frame with DIC, DIC gradient, and flux columns
#' @export
carbonate_calculations <- function(
  rga_calibrated,
  seaphox_df,
  flux_dataset,
  rga_adv_flux,
  length_scale,
  sensor_separation
) {
  seaphox_carb_periods <- build_seaphox_carb_periods(seaphox_df)
  ustar_data <- get_ustar(flux_dataset)

  add_dic(
    rga_calibrated = rga_calibrated,
    seaphox_carb_periods = seaphox_carb_periods,
    sensor_separation = sensor_separation
  ) |>
    left_join(ustar_data, by = "timestamp") |>
    add_scalar_flux(
      grad_var = dic_gradient_umol_l_m,
      flux_name = "dic_flux",
      length_scale = length_scale
    ) |>
    add_scalar_flux(
      grad_var = dic_gradient_predicted_umol_l_m,
      flux_name = "dic_flux_corrected",
      length_scale = length_scale
    ) |>
    left_join(
      select(
        rga_adv_flux,
        timestamp,
        ox_flux,
        co2_flux
      ),
      by = "timestamp"
    ) |>
    select(
      timestamp,
      par,
      seaphox_temp_c,
      seaphox_salinity_psu,
      seaphox_pH,
      adv_temp,
      alk_high_umol_l,
      dic_high_umol_l,
      dic_low_umol_l,
      dic_low_predicted_umol_l,
      dic_mean_umol_l,
      dic_gradient_umol_l_m,
      dic_gradient_predicted_umol_l_m,
      dic_flux,
      dic_flux_corrected,
      ox_flux,
      co2_flux
    )
}
