# PAR estimation

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

#' Generate DLI for a range of dates
make_dli_df <- function(start_date, end_date, crds) {
  date_seq <- seq.Date(as.Date(start_date), as.Date(end_date), by = "day")
  dli_df <- data.frame(date = date_seq) |>
    rowwise() |>
    mutate(
      par_values = list(
        calculate_par(
          crds,
          seq.POSIXt(
            as.POSIXct(date),
            as.POSIXct(date + 1) - 3600,
            by = "hour"
          )
        )
      ),
      dli_mol_m2_day = sum(unlist(par_values) * 0.0036) # 3600 / 1e6
    ) |>
    select(date, dli_mol_m2_day)
  return(dli_df)
}
