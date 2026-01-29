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
      hour = hour(timestamp),
      diff_noon = timestamp - 16800,
      solar_hour = hour(diff_noon),
      month = as.integer(month(timestamp)),
      season = case_when(
        month %in% 10:12 ~ "Fall",
        month %in% 1:3 ~ "Winter",
        month %in% 4:6 ~ "Spring",
        TRUE ~ "Summer"
      )
    ) |>
    group_by(solar_hour) |>
    summarise(
      across(
        c(-timestamp),
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
      hour = hour(timestamp),
      diff_noon = timestamp - 16800,
      solar_hour = hour(diff_noon),
      month = month(timestamp),
      season = case_when(
        month %in% 10:12 ~ "Fall",
        month %in% 1:3 ~ "Winter",
        month %in% 4:6 ~ "Spring",
        TRUE ~ "Summer"
      )
    ) |>
    group_by(solar_hour, month) |>
    summarise(
      across(
        c(-timestamp),
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
      hour = hour(timestamp),
      diff_noon = timestamp - 16800,
      solar_hour = hour(diff_noon),
      month = month(timestamp),
      season = case_when(
        month %in% 10:12 ~ "Fall",
        month %in% 1:3 ~ "Winter",
        month %in% 4:6 ~ "Spring",
        TRUE ~ "Summer"
      )
    ) |>
    group_by(solar_hour, season) |>
    summarise(
      across(
        c(-timestamp),
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
