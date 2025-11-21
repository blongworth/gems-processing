library(dplyr)
library(lubridate)

utils::globalVariables(c("minute_of_hour", "bin_index", "bin_start_minute", "date_part", "hour", "bin_id", "bin_time", "inlet"))

#' Bin timeseries into 7.5-minute chunks by clock time
#'
#' Groups a timeseries into 7.5-minute chunks starting at the top of each hour,
#' summarizes by averaging a specified range of minutes in each chunk, and adds alternating inlet labels.
#'
#' @param data A tibble or data frame with a datetime column and numeric value columns to summarize
#' @param datetime_col Name of the datetime column (as string, e.g., "timestamp")
#' @param value_cols Names of numeric columns to average (character vector)
#' @param avg_start Starting minute within each chunk to begin averaging (default: 1)
#' @param avg_end Ending minute within each chunk to stop averaging (default: 7)
#'
#' @return A tibble with binned data, averaged values, bin times, and inlet labels
#'
#' @examples
#' \dontrun{
#' library(lubridate)
#' df <- tibble(
#'   timestamp = seq(ymd_hms("2023-01-01 00:00:00"), by = "1 min", length.out = 60),
#'   value1 = rnorm(60),
#'   value2 = rnorm(60)
#' )
#' result <- bin_timeseries(df, "timestamp", c("value1", "value2"))
#' result <- bin_timeseries(df, "timestamp", c("value1", "value2"), avg_start = 2, avg_end = 6)
#' }
#'
#' @export
bin_timeseries <- function(data, datetime_col, value_cols, 
                           avg_start = 2, avg_end = 7) {
  
  data %>%
    # Calculate minutes within the hour
    mutate(
      minute_of_hour = minute(!!sym(datetime_col)),
      hour = hour(!!sym(datetime_col)),
      date_part = as_date(!!sym(datetime_col))
    ) %>%
    # Assign to 7.5-minute bins based on clock time
    # Bins: 0:00-7:29, 7:30-14:59, 15:00-22:29, 22:30-29:59, 30:00-37:29, etc.
    mutate(
      bin_index = floor(minute_of_hour / 7.5),
      bin_start_minute = bin_index * 7.5
    ) %>%
    # Filter to only include minutes avg_start through avg_end of each chunk
    filter(minute_of_hour >= bin_start_minute + avg_start &
           minute_of_hour < bin_start_minute + avg_end) %>%
    # Create bin identifier
    mutate(
      bin_id = paste0(date_part, "_", hour, "_", bin_index)
    ) %>%
    # Group and summarize
    group_by(bin_id, date_part, hour, bin_index) %>%
    summarise(
      across(all_of(value_cols), \(x) mean(x, na.rm = TRUE)),
      bin_start_minute = first(bin_start_minute),
      .groups = "drop"
    ) %>%
    # Arrange chronologically
    arrange(date_part, hour, bin_index) %>%
    # Add alternating inlet labels
    mutate(
      inlet = ifelse(bin_index %% 2 == 0, "low", "high")
    ) %>%
    # Reconstruct approximate bin start time (top of hour + bin offset)
    mutate(
      bin_time = as_datetime(paste0(date_part, " ", 
                             sprintf("%02d", hour), ":", 
                             sprintf("%02.0f", bin_start_minute), ":00"))
    ) %>%
    # Clean up temporary columns and reorder
    select(bin_time, inlet, all_of(value_cols)) %>%
    as_tibble()
}

#' Process RGA data into a wide format
#'
#' This function takes raw RGA data in a long format and transforms it into a
#' wide format. It identifies measurement cycles based on a specific mass
#' (e.g., mass 18 for water), calculates a mean timestamp for each cycle,
#' and then pivots the data so that each mass has its own column.
#'
#' @param data A tibble or data frame with RGA data. Expected columns are
#'   `send`, `timestamp`, `mass`, and `pressure`.
#' @param cycle_start_mass The mass value that indicates the start of a new
#'   measurement cycle. Defaults to 18.
#'
#' @return A tibble in wide format, 
#' with a row for each cycle and columns for each mass.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'rga_data' is a tibble with the required columns
#' rga_wide <- process_rga_to_wide(rga_data)
#' }
process_rga_to_wide <- function(data, cycle_start_mass = 18) {
    data |>
        mutate(cycle = cumsum(mass == cycle_start_mass)) %>%
        group_by(cycle) %>%
        mutate(cycle_ts = mean(timestamp)) %>%
        ungroup() %>%
        select(timestamp = cycle_ts, mass, pressure) %>%
        tidyr::pivot_wider(
            names_from = mass, names_prefix = "mass_",
            values_from = pressure,
            values_fn = mean
        )
}

#' Rotate coordinate axes to minimize mean z-velocity
#'
#' Takes velocity components (u, v, w) and finds the optimal rotation angles
#' that minimize the mean z-velocity component. Uses a two-step rotation approach:
#' first rotating around the z-axis, then around the rotated y-axis.
#'
#' @param data A tibble or data frame with velocity columns
#' @param u_col Name of x-velocity column (as string)
#' @param v_col Name of y-velocity column (as string)
#' @param w_col Name of z-velocity column (as string)
#'
#' @return A list containing:
#'   - rotated_data: tibble with rotated velocity components (u_rot, v_rot, w_rot)
#'   - theta: rotation angle around z-axis (radians)
#'   - phi: rotation angle around y-axis (radians)
#'   - mean_w_rotated: mean z-velocity after rotation
#'
#' @examples
#' \dontrun{
#' result <- rotate_to_minimize_z(data, "u", "v", "w")
#' rotated_df <- result$rotated_data
#' }
#'
#' @export
rotate_to_minimize_z <- function(data, u_col, v_col, w_col) {
  
  # Extract velocity vectors
  u <- data[[u_col]]
  v <- data[[v_col]]
  w <- data[[w_col]]
  
  # Step 1: Find theta to align mean horizontal velocity with x-axis
  # Minimize by rotating in xy-plane
  mean_u <- mean(u, na.rm = TRUE)
  mean_v <- mean(v, na.rm = TRUE)
  theta <- atan2(mean_v, mean_u)
  
  # Step 2: Rotate to find phi that minimizes mean w
  # First rotate u,v by -theta
  u_temp <- u * cos(theta) + v * sin(theta)
  v_temp <- -u * sin(theta) + v * cos(theta)
  
  # Then find phi to align horizontal velocity with x-axis in xz-plane
  mean_u_temp <- mean(u_temp, na.rm = TRUE)
  mean_w <- mean(w, na.rm = TRUE)
  phi <- atan2(mean_w, mean_u_temp)
  
  # Apply full rotation
  u_rot <- u_temp * cos(phi) + w * sin(phi)
  v_rot <- v_temp
  w_rot <- -u_temp * sin(phi) + w * cos(phi)
  
  # Return rotated data and angles
  result <- data %>%
    mutate(
      u_rot = u_rot,
      v_rot = v_rot,
      w_rot = w_rot
    ) %>%
    select(all_of(c(setdiff(names(data), c(u_col, v_col, w_col)))), 
           u_rot, v_rot, w_rot)
  
  # list(
  #   rotated_data = result,
  #   theta = theta,
  #   phi = phi,
  #   mean_w_rotated = mean(w_rot, na.rm = TRUE)
  # )
    
    result
}