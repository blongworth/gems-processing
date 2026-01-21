# analysis functions

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
      co2_gradient_umol_l_m_hourly = mean(co2_gradient_umol_l_m, na.rm = TRUE),
      co2_gradient_umol_l_m_hourly_sd = sd(co2_gradient_umol_l_m, na.rm = TRUE),
      ox_flux_hourly = mean(ox_flux, na.rm = TRUE),
      ox_flux_hourly_sd = sd(ox_flux, na.rm = TRUE),
      .groups = "drop"
    )
}
