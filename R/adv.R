# functions for adv data
library(arrow)
library(dplyr)
library(readr)
library(tidyr)

load_and_bin_adv <- function(adv_raw_file, moves_file, min_correlation) {
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
      datetime_col = "timestamp"
    ) |>
    summarize_binned_timeseries(
      value_cols = c("pressure", "u", "v", "w")
    ) |>
    mutate(
      cur_speed = sqrt(v^2 + u^2),
      cur_dir = atan2(u, v) * (180 / pi)
    )

  lander_moves <- read_csv(moves_file) |>
    dplyr::rename(change_timestamp = timestamp) |>
    dplyr::mutate(lander_position = dplyr::row_number())

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
    nest() |>
    mutate(
      data = map(data, \(x) {
        rotate_to_minimize_z(x, "u", "v", "w")$rotated_data
      })
    ) |>
    unnest(data) |>
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
#'   - theta: rotation angle around z-axis (degrees)
#'   - phi: rotation angle around y-axis (degrees)
#'   - theta_rad: rotation angle around z-axis (radians)
#'   - phi_rad: rotation angle around y-axis (radians)
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
  theta_rad <- atan2(mean_v, mean_u)

  # Step 2: Rotate to find phi that minimizes mean w
  # First rotate u,v by -theta
  u_temp <- u * cos(theta_rad) + v * sin(theta_rad)
  v_temp <- -u * sin(theta_rad) + v * cos(theta_rad)

  # Then find phi to align horizontal velocity with x-axis in xz-plane
  mean_u_temp <- mean(u_temp, na.rm = TRUE)
  mean_w <- mean(w, na.rm = TRUE)
  phi_rad <- atan2(mean_w, mean_u_temp)

  # Apply full rotation
  u_rot <- u_temp * cos(phi_rad) + w * sin(phi_rad)
  v_rot <- v_temp
  w_rot <- -u_temp * sin(phi_rad) + w * cos(phi_rad)

  theta <- theta_rad * 180 / pi
  phi <- phi_rad * 180 / pi

  # Return rotated data and angles
  result <- data |>
    mutate(
      u_rot = u_rot,
      v_rot = v_rot,
      w_rot = w_rot
    ) |>
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
    theta_rad = theta_rad,
    phi_rad = phi_rad,
    mean_w_rotated = mean(w_rot, na.rm = TRUE)
  )

  #  result
}

calculate_adv_rotations <- function(
  adv_data,
  min_velocity = 0.045,
  max_velocity = 0.3
) {
  rotation_counts <- adv_data |>
    group_by(lander_position) |>
    summarize(
      n_total = n(),
      .groups = "drop"
    )

  adv_data |>
    mutate(horizontal_velocity = sqrt(u^2 + v^2)) |>
    filter(
      horizontal_velocity >= min_velocity,
      horizontal_velocity < max_velocity
    ) |>
    group_by(lander_position) |>
    summarize(
      n = n(),
      mean_horizontal_velocity = mean(horizontal_velocity, na.rm = TRUE),
      mean_u = mean(u, na.rm = TRUE),
      mean_v = mean(v, na.rm = TRUE),
      mean_w = mean(w, na.rm = TRUE),
      .groups = "drop"
    ) |>
    right_join(rotation_counts, by = join_by(lander_position)) |>
    mutate(
      n = replace_na(n, 0L),
      n_removed = n_total - n,
      min_velocity = min_velocity,
      max_velocity = max_velocity
    ) |>
    mutate(
      theta_rad = atan2(mean_v, mean_u),
      mean_u_rot1 = mean_u * cos(theta_rad) + mean_v * sin(theta_rad),
      mean_v_rot1 = -mean_u * sin(theta_rad) + mean_v * cos(theta_rad),
      phi_rad = atan2(mean_w, mean_u_rot1),
      theta_deg = theta_rad * 180 / pi,
      phi_deg = phi_rad * 180 / pi,
      mean_u_rot = mean_u_rot1 * cos(phi_rad) + mean_w * sin(phi_rad),
      mean_v_rot = mean_v_rot1,
      mean_w_rot = -mean_u_rot1 * sin(phi_rad) + mean_w * cos(phi_rad)
    ) |>
    select(
      lander_position,
      min_velocity,
      max_velocity,
      n_total,
      n,
      n_removed,
      mean_horizontal_velocity,
      theta_deg,
      phi_deg,
      theta_rad,
      phi_rad,
      mean_u,
      mean_v,
      mean_w,
      mean_u_rot,
      mean_v_rot,
      mean_w_rot
    )
}

calculate_adv_rotations_with_rotate <- function(
  adv_data,
  min_velocity = 0.045,
  max_velocity = 0.3
) {
  rotation_counts <- adv_data |>
    group_by(lander_position) |>
    summarize(
      n_total = n(),
      .groups = "drop"
    )

  rotations <- adv_data |>
    mutate(horizontal_velocity = sqrt(u^2 + v^2)) |>
    filter(
      horizontal_velocity >= min_velocity,
      horizontal_velocity < max_velocity
    ) |>
    group_by(lander_position) |>
    group_modify(\(.x, .y) {
      rotation <- rotate_to_minimize_z(.x, "u", "v", "w")
      rotated_data <- rotation$rotated_data

      tibble::tibble(
        n = nrow(.x),
        mean_horizontal_velocity = mean(.x$horizontal_velocity, na.rm = TRUE),
        theta_deg = rotation$theta,
        phi_deg = rotation$phi,
        theta_rad = rotation$theta_rad,
        phi_rad = rotation$phi_rad,
        mean_u = mean(.x$u, na.rm = TRUE),
        mean_v = mean(.x$v, na.rm = TRUE),
        mean_w = mean(.x$w, na.rm = TRUE),
        mean_u_rot = mean(rotated_data$u_rot, na.rm = TRUE),
        mean_v_rot = mean(rotated_data$v_rot, na.rm = TRUE),
        mean_w_rot = mean(rotated_data$w_rot, na.rm = TRUE)
      )
    }) |>
    ungroup()

  rotations |>
    right_join(rotation_counts, by = join_by(lander_position)) |>
    mutate(
      n = replace_na(n, 0L),
      n_removed = n_total - n,
      min_velocity = min_velocity,
      max_velocity = max_velocity
    ) |>
    select(
      lander_position,
      min_velocity,
      max_velocity,
      n_total,
      n,
      n_removed,
      mean_horizontal_velocity,
      theta_deg,
      phi_deg,
      theta_rad,
      phi_rad,
      mean_u,
      mean_v,
      mean_w,
      mean_u_rot,
      mean_v_rot,
      mean_w_rot
    )
}
