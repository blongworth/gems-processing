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


#' Interpolate cyclic inlet data on a 15-minute grid
#'
#' @param df Data frame with a POSIXct-compatible \code{timestamp} column.
#' @param cols_to_interp Character vector of numeric column names to interpolate;
#'   if \code{NULL}, all numeric columns except \code{timestamp} and \code{inlet} are used.
#' @param block_start_sec Numeric. Start (in seconds from block start) of the
#'   sub-window within each 7.5-minute block used for interpolation.
#' @param block_end_sec Numeric. End (in seconds from block start) of the
#'   sub-window within each 7.5-minute block used for interpolation.
#'
#' @return A data frame with one row per 15-minute grid time and columns:
#'   \code{time_grid}, \code{timestamp_low}, \code{<var>_low},
#'   \code{timestamp_high}, \code{<var>_high}.
#'
#' @details
#'   Assumes a repeating 15-minute cycle starting at the top of the hour, with
#'   7.5 minutes at inlet \code{"low"} followed by 7.5 minutes at inlet \code{"high"}.
#'   Data are first restricted to the given sub-window within each block, then
#'   piecewise-linear interpolation is performed at the midpoints of each
#'   low/high block for each 15-minute cycle.
#'
#' @examples
#' \dontrun{
#' res <- interpolate_inlets(
#'   df,
#'   cols_to_interp = c("temp", "pressure"),
#'   block_start_sec = 60,
#'   block_end_sec = 420
#' )
#' }
interpolate_inlets <- function(
  df,
  cols_to_interp = NULL,
  block_start_sec = 60,
  block_end_sec = 420
) {
  cycle_length_min <- 15
  half_cycle_min <- 7.5
  cycle_length_sec <- cycle_length_min * 60
  half_cycle_sec <- half_cycle_min * 60

  df <- df %>%
    mutate(
      timestamp = as.POSIXct(timestamp, tz = "UTC")
    ) %>%
    arrange(timestamp)

  # auto-detect columns to interpolate if not provided
  if (is.null(cols_to_interp)) {
    cols_to_interp <- df %>%
      select(-timestamp, -any_of("inlet")) %>%
      select(where(is.numeric)) %>%
      names()
  }

  # origin time: floor to hour
  origin_time <- floor_date(min(df$timestamp), unit = "hour")

  # assign inlet (low first half, high second half) and block_index/start
  df <- df %>%
    mutate(
      offset_sec = as.numeric(difftime(timestamp, origin_time, units = "secs")),
      cycle_pos = offset_sec %% cycle_length_sec,
      inlet = if_else(cycle_pos < half_cycle_sec, "low", "high"),
      block_index = floor(offset_sec / half_cycle_sec),
      block_start = origin_time + block_index * half_cycle_sec,
      time_in_block_sec = as.numeric(difftime(
        timestamp,
        block_start,
        units = "secs"
      ))
    )

  # enforce sensible bounds on block window
  if (block_start_sec < 0) {
    warning("block_start_sec < 0; clamping to 0.")
    block_start_sec <- 0
  }
  if (block_end_sec > half_cycle_sec) {
    warning("block_end_sec > half block length; clamping to half_cycle_sec.")
    block_end_sec <- half_cycle_sec
  }
  if (block_start_sec >= block_end_sec) {
    stop("block_start_sec must be < block_end_sec.")
  }

  # keep only the middle part of each block as defined by [block_start_sec, block_end_sec]
  df_mid <- df %>%
    filter(
      time_in_block_sec >= block_start_sec,
      time_in_block_sec <= block_end_sec
    )

  if (nrow(df_mid) == 0) {
    stop(
      "No data remaining after applying block_start_sec and block_end_sec filters."
    )
  }

  # 15-minute grid aligned with origin_time
  grid_15min <- seq(
    from = origin_time,
    to = ceiling_date(max(df$timestamp), unit = "hour"),
    by = "15 min"
  )

  # midpoint times inside low and high windows (same as before)
  half_window_sec <- half_cycle_sec / 2 # 3.75 min
  times_low <- grid_15min + half_window_sec
  times_high <- grid_15min + (half_cycle_sec + half_window_sec)

  # helper: interpolate many columns at specified target_times using df_mid
  interp_many_cols <- function(df_sub, target_times, cols) {
    x <- as.numeric(df_sub$timestamp)
    xout <- as.numeric(target_times)

    interpolated_list <- map(cols, function(col_name) {
      y <- df_sub[[col_name]]

      out <- approx(
        x = x,
        y = y,
        xout = xout,
        method = "linear",
        rule = 2
      )

      tibble(
        timestamp = as.POSIXct(
          out$x,
          origin = "1970-01-01",
          tz = attr(df_sub$timestamp, "tzone")
        ),
        !!col_name := out$y
      )
    })

    reduce(interpolated_list, full_join, by = "timestamp") %>%
      arrange(timestamp)
  }

  # split the middle-of-block data by inlet
  df_low_mid <- df_mid %>% filter(inlet == "low")
  df_high_mid <- df_mid %>% filter(inlet == "high")

  if (nrow(df_low_mid) == 0) {
    stop("No 'low' inlet data remaining after block filtering.")
  }
  if (nrow(df_high_mid) == 0) {
    stop("No 'high' inlet data remaining after block filtering.")
  }

  interp_low <- interp_many_cols(df_low_mid, times_low, cols_to_interp) %>%
    rename(timestamp_low = timestamp)

  interp_high <- interp_many_cols(df_high_mid, times_high, cols_to_interp) %>%
    rename(timestamp_high = timestamp)

  result <- tibble(time_grid = grid_15min) %>%
    bind_cols(
      tibble(timestamp_low = interp_low$timestamp_low),
      interp_low %>%
        select(-timestamp_low) %>%
        rename_with(~ paste0(.x, "_low"), all_of(cols_to_interp)),
      tibble(timestamp_high = interp_high$timestamp_high),
      interp_high %>%
        select(-timestamp_high) %>%
        rename_with(~ paste0(.x, "_high"), all_of(cols_to_interp))
    ) %>%
    select(
      time_grid,
      timestamp_low,
      ends_with("_low"),
      timestamp_high,
      ends_with("_high")
    )

  result
}
