# resampling and binning functions

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
#' @export
bin_timeseries <- function(
  data,
  datetime_col,
  value_cols,
  avg_start = 60,
  avg_end = 360
) {
  data %>%
    # Calculate seconds within the hour
    mutate(
      second_of_hour = minute(!!sym(datetime_col)) *
        60 +
        second(!!sym(datetime_col)),
      hour = hour(!!sym(datetime_col)),
      date_part = as_date(!!sym(datetime_col))
    ) %>%
    # Assign to 7.5-minute bins based on clock time
    # Bins: 0-449s, 450-899s, 900-1349s, etc. (7.5 min = 450 seconds)
    mutate(
      bin_index = floor(second_of_hour / 450),
      bin_start_second = bin_index * 450,
      inlet = ifelse(bin_index %% 2 == 0, "low", "high")
    ) %>%
    # Filter to only include seconds avg_start through avg_end of each chunk
    filter(
      second_of_hour >= bin_start_second + avg_start &
        second_of_hour < bin_start_second + avg_end
    ) %>%
    # Create bin identifier
    mutate(
      bin_id = paste0(date_part, "_", hour, "_", bin_index)
    )
}

#' Summarize binned timeseries by bin
#'
#' summarizes by averaging a specified range of minutes in each chunk, and adds alternating inlet labels.
#'
#' @param data A tibble with bin id fields from `bin_timeseries()`
#' @param value_cols Names of numeric columns to average (character vector)
#'
#' @return A tibble with binned data, averaged values, bin times, and inlet labels
#'
#' @export
summarize_binned_timeseries <- function(
  data,
  value_cols
) {
  data %>%
    # Group and summarize
    group_by(bin_id, date_part, hour, bin_index) %>%
    summarise(
      across(all_of(value_cols), \(x) mean(x, na.rm = TRUE)),
      inlet = first(inlet),
      bin_start_second = first(bin_start_second),
      .groups = "drop"
    ) %>%
    # Arrange chronologically
    arrange(date_part, hour, bin_index) %>%
    # Reconstruct approximate bin start time (top of hour + bin offset)
    mutate(
      bin_time = as_datetime(paste0(
        date_part,
        " ",
        sprintf("%02d", hour),
        ":",
        sprintf("%02.0f", bin_start_second %/% 60),
        ":",
        sprintf("%02.0f", bin_start_second %% 60)
      ))
    ) %>%
    # Clean up temporary columns and reorder
    select(bin_time, inlet, all_of(value_cols)) %>%
    as_tibble()
}
