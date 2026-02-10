# resampling and binning functions
# Functions for processing alternating inlet timeseries data
# Data alternates between two inlets every 7.5 minutes starting at the top of the hour

#' Assign inlet labels to timeseries data
#'
#' @param df Data frame with a POSIXct "timestamp" column
#' @return Data frame with added 'inlet' and 'period_start' columns
#' @export
assign_inlets <- function(df) {
  df %>%
    mutate(
      seconds_in_hour = as.numeric(difftime(
        timestamp,
        floor_date(timestamp, "hour"),
        units = "secs"
      )),
      period_index = floor(seconds_in_hour / 450),
      inlet = if_else(period_index %% 2 == 0, "low", "high"),
      period_start = floor_date(timestamp, "hour") + seconds(period_index * 450)
    ) %>%
    select(-seconds_in_hour, -period_index)
}

filter_inlet_window <- function(df, window_start = 0, window_end = 450) {
  df %>%
    mutate(
      seconds_in_period = as.numeric(difftime(
        timestamp,
        period_start,
        units = "secs"
      ))
    ) %>%
    filter(
      seconds_in_period >= window_start,
      seconds_in_period < window_end
    )
}

#' Calculate mean within a window for each period
#'
#' @param df Data frame with inlet assignments (output from assign_inlets)
#' @param window_start Start of window in seconds from period start (default: 0)
#' @param window_end End of window in seconds from period start (default: 450)
#' @return Data frame with period means and mean time for all numeric columns
#' @export
calculate_period_means <- function(df) {
  # Identify numeric columns to average
  exclude_cols <- c("timestamp", "inlet", "period_start")
  numeric_cols <- names(df)[
    sapply(df, is.numeric) & !names(df) %in% exclude_cols
  ]

  if (length(numeric_cols) == 0) {
    stop("No numeric columns found to calculate means")
  }

  df %>%
    group_by(inlet, period_start) %>%
    summarise(
      mean_time = mean(timestamp),
      n_samples = n(),
      across(
        all_of(numeric_cols),
        \(.x) mean(.x, na.rm = TRUE),
        # .names = "{.col}_mean"
        .names = "{.col}"
      ),
      .groups = "drop"
    ) %>%
    mutate(
      mean_time = as.POSIXct(
        mean_time,
        origin = "1970-01-01",
        tz = attr(df$timestamp, "tzone")
      )
    ) |>
    rename(timestamp = period_start)
}


#' Interpolate inlet data to regular 15-minute grid
#'
#' @param period_means Data frame with period means (output from calculate_period_means)
#' @return Data frame with interpolated values on 15-minute grid for each inlet
#' @export
interpolate_to_grid <- function(period_means) {
  # Create 15-minute grid from data extent
  time_grid <- seq(
    from = floor_date(min(period_means$mean_time), "15 min"),
    to = ceiling_date(max(period_means$mean_time), "15 min"),
    by = "15 min"
  )

  # Identify columns to interpolate
  mean_cols <- names(period_means)[
    !names(period_means) %in%
      c("inlet", "timestamp", "mean_time", "n_samples", "seconds_in_period")
  ]

  # Interpolate each inlet separately
  period_means %>%
    select(-timestamp) %>%
    group_by(inlet) %>%
    arrange(mean_time) %>%
    group_modify(\(.x, .y) {
      # Create data frame with interpolated values for all columns
      result <- tibble(timestamp = time_grid)

      for (col in mean_cols) {
        interp <- approx(
          x = .x$mean_time,
          y = .x[[col]],
          xout = time_grid,
          method = "linear",
          rule = 2
        )
        result[[col]] <- interp$y
      }
      result
    }) %>%
    ungroup()
}


#' Complete pipeline: process raw data to interpolated grid
#'
#' @param df Data frame with raw timeseries data (must have "timestamp" column)
#' @param window_start Start of averaging window in seconds (default: 0)
#' @param window_end End of averaging window in seconds (default: 450)
#' @param interp_method Interpolation method: "linear" or "nearest" (default: "linear")
#' @return Data frame with interpolated values on 15-minute grid for all numeric columns
#' @export
interpolate_rga <- function(
  df,
  value_cols = NULL,
  window_start = 0,
  window_end = 450
) {
  df <- df %>%
    assign_inlets() %>%
    filter_inlet_window(window_start, window_end) %>%
    calculate_period_means() %>%
    interpolate_to_grid()

  if (!is.null(value_cols)) {
    df <- df %>%
      select(timestamp, inlet, all_of(value_cols))
  }

  df
}

widen_binned_rga <- function(rga_binned, value_cols = NULL) {
  if (is.null(value_cols)) {
    value_cols <- names(rga_binned)[
      !names(rga_binned) %in% c("timestamp", "inlet")
    ]
  }
  rga_binned |>
    pivot_wider(
      names_from = inlet,
      values_from = all_of(value_cols),
      names_glue = "{.value}_{inlet}"
    ) %>%
    arrange(timestamp)
}
