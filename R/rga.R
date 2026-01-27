# RGA data processing functions

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
process_rga_to_wide <- function(data, cycle_start_mass = 18) {
  data |>
    mutate(cycle = cumsum(mass == cycle_start_mass)) %>%
    group_by(cycle) %>%
    mutate(cycle_ts = mean(timestamp)) %>%
    ungroup() %>%
    select(timestamp = cycle_ts, mass, pressure) %>%
    tidyr::pivot_wider(
      names_from = mass,
      names_prefix = "mass_",
      values_from = pressure,
      values_fn = mean
    )
}

load_and_widen_rga <- function(rga_file) {
  open_dataset(rga_file) |>
    collect() |>
    process_rga_to_wide()
}

remove_bad_rga_periods <- function(rga_data, bad_times) {
  rga_data |>
    anti_join(
      bad_times,
      by = join_by(timestamp >= start_time, timestamp <= end_time)
    )
}

normalize_rga_by_argon <- function(rga_data) {
  rga_data |>
    mutate(
      across(starts_with("mass_"), ~ .x / mass_40, .names = "{.col}_40")
    )
}

pair_gradient_measurements <- function(rga_binned) {
  rga_binned |>
    group_by(grp = cumsum(inlet == "high")) |>
    mutate(mean_timestamp = mean(bin_time)) |>
    ungroup() |>
    pivot_wider(
      id_cols = c(grp, mean_timestamp),
      names_from = inlet,
      values_from = starts_with("mass_"),
      names_sep = "_"
    ) |>
    select(-grp) |>
    mutate(
      timestamp = lubridate::round_date(mean_timestamp, unit = "15 minutes"),
      # TODO: rename "high_mean" to just "mean"
      across(
        ends_with("_high"),
        \(x) (x + get(sub("_high$", "_low", cur_column()))) / 2,
        .names = "{.col}_mean"
      ),
      across(
        ends_with("_high"),
        \(x) x - get(sub("_high$", "_low", cur_column())),
        .names = "{.col}_diff"
      ),
    ) |>
    dplyr::relocate(timestamp, .before = dplyr::everything())
}

# add status temperature
add_status_temp <- function(rga_df, status_file, bad_times) {
  status_temp_df <- open_dataset(status_file) |>
    select(timestamp, adv_temp = temp) |>
    collect() |>
    mutate(
      timestamp = lubridate::floor_date(timestamp, unit = "15 minutes")
    ) |>
    group_by(timestamp) |>
    summarize(adv_temp = mean(adv_temp, na.rm = TRUE))

  rga_df |>
    left_join(status_temp_df, by = join_by(timestamp)) |>
    anti_join(
      bad_times,
      by = join_by(timestamp >= start_time, timestamp <= end_time)
    )
}
