# functions for adv data

load_and_bin_adv <- function(adv_raw_file, moves_file, min_correlation = NULL) {
  ds <- open_dataset(adv_raw_file)
  if (!is.null(min_correlation)) {
    ds <- ds |>
      filter(
        corr1 > min_correlation,
        corr2 > min_correlation,
        corr3 > min_correlation
      )
  }
  df_binned <- ds |>
    select(timestamp, pressure, u, v, w) |>
    collect() |>
    bin_timeseries(
      datetime_col = "timestamp",
      value_cols = c("pressure", "u", "v", "w")
    ) |>
    mutate(
      cur_speed = sqrt(v^2 + u^2),
      cur_dir = atan2(u, v) * (180 / pi)
    )

  lander_moves <- read_csv(moves_file) |>
    rename(change_timestamp = timestamp) |>
    mutate(lander_position = row_number())

  df_binned_mv <- df_binned |>
    arrange(bin_time) |>
    left_join(
      lander_moves,
      by = join_by(closest(bin_time > change_timestamp))
    ) |>
    fill(lander_position, .direction = "down") |>
    mutate(lander_position = replace_na(lander_position, 1)) |>
    select(-change_timestamp)

  df_bin_rot <- df_binned_mv |>
    group_by(lander_position) |>
    nest() %>%
    mutate(
      data = map(data, \(x) {
        rotate_to_minimize_z(x, "u", "v", "w")$rotated_data
      })
    ) %>%
    unnest(data) %>%
    ungroup()

  df_bin_rot
}

add_adv <- function(rga_data, adv_bin_rot_df) {
  adv_binned <- adv_bin_rot_df |>
    group_by(grp = cumsum(inlet == "high")) |>
    summarize(
      mean_timestamp = mean(bin_time),
      timestamp = lubridate::round_date(
        mean_timestamp,
        unit = "15 minutes"
      ),
      across(c(pressure, cur_speed, cur_dir, u_rot, v_rot, w_rot), \(x) {
        mean(x)
      })
    ) |>
    dplyr::relocate(timestamp, .before = dplyr::everything())

  adv_select <- adv_binned |>
    select(
      timestamp,
      pressure,
      cur_speed,
      cur_dir,
      u = u_rot,
      v = v_rot,
      w = w_rot
    )

  rga_data |>
    left_join(adv_select, by = join_by(timestamp)) |>
    drop_na()
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
    select(
      all_of(c(setdiff(names(data), c(u_col, v_col, w_col)))),
      u_rot,
      v_rot,
      w_rot
    )

  list(
    rotated_data = result,
    theta = theta,
    phi = phi,
    mean_w_rotated = mean(w_rot, na.rm = TRUE)
  )

  #  result
}
