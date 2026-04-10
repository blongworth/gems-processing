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
      nem_mmol_m2_h = sum(ox_flux) / 24
    )
}

#' @title Calculate Monthly Net Ecosystem Metabolism (NEM)
calculate_nem_monthly <- function(monthly_stats) {
  monthly_stats |>
    group_by(month) |>
    summarise(
      nem_mmol_m2_day = sum(ox_flux_mean)
    )
}

#' Create dataset with pH from Seaphox
#'
#' @param seaphox_df Data frame with Seaphox data
#' @param start_time Start time for filtering
#' @param end_time End time for filtering
#' @return Data frame with timestamp and pH columns
#' @export
carbonate_calculations <- function(
  rga_calibrated,
  seaphox_df,
  flux_dataset,
  rga_adv_flux,
  length_scale,
  sensor_separation
) {
  df <- seaphox_df |>
    filter(
      timestamp >= as.POSIXct("2025-06-26"),
      timestamp <= as.POSIXct("2025-07-16 14:00:00")
    ) |>
    mutate(timestamp = floor_date(timestamp, unit = "15 min")) |>
    group_by(timestamp) |>
    summarise(
      seaphox_temp_c = mean(seaphox_temp_c),
      seaphox_salinity_psu = mean(seaphox_salinity_psu),
      seaphox_pH = mean(seaphox_pH, na.rm = TRUE)
    ) |>
    left_join(rga_calibrated, by = "timestamp") |>
    mutate(
      dic_high_umol_l = seacarb::carb(
        flag = 1,
        var1 = seaphox_pH,
        var2 = co2_high_umol_l,
        S = seaphox_salinity_psu,
        T = seaphox_temp_c
      )$DIC,
      dic_low_umol_l = seacarb::carb(
        flag = 1,
        var1 = seaphox_pH,
        var2 = co2_low_umol_l,
        S = seaphox_salinity_psu,
        T = seaphox_temp_c
      )$DIC,
      dic_gradient_umol_l_m = (dic_high_umol_l - dic_low_umol_l) /
        sensor_separation
    )

  df |>
    left_join(
      flux_dataset |>
        select(timestamp, Ustar) |>
        mutate(timestamp = floor_date(timestamp, unit = "15 min")),
      by = "timestamp"
    ) |>
    mutate(
      lscale = length_scale,
      dic_flux = -1 * Ustar * 0.41 * lscale * dic_gradient_umol_l_m * 3600
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
      dic_flux,
      ox_flux,
      co2_flux
    )
}
