# Functions for reading and parsing raw LECS data

# TODO: add option to correct future status dates using met data before adv cor

#' LECS raw file parser
#'
#' Parse raw files to a df with row_num and type
#' Remove type chars from beginning of lines
#'
#' @param files a vector of file paths
#'
#' @return a dataframe with row_num, type, and unparsed data
#' @export
#' @importFrom stringr str_sub
lecs_read_file <- function(file) {
  l <- readLines(file)
  type <- str_sub(l, 1, 1)
  data <- str_sub(l, start = 3)
  data.frame(row_num = 1:length(l), type = type, data = data)
}

#' Parse LECS metadata
#'
#' Add row numbers, line type, number and line of send (web data),
#'
#' @param df_raw a raw LECS dataframe
#'
#' @return a dataframe with added metadata
#' @export
lecs_add_metadata <- function(df_raw) {
  df_raw |>
    mutate(row_num_i = dplyr::row_number(),
           type = stringr::str_match(X1, "\\[\\d+\\]([RGDMS$!]):?")[,2],
           send = cumsum(type == "$"),
           line = as.integer(stringr::str_match(X1, "\\[(\\d+)\\][RGDMS$!]:?")[,2]),
           data = stringr::str_remove(X1, "\\[\\d+\\][RGDMS$!]:?")) |>
    select(row_num_i, send, type, line, data) |>
    group_by(send) |>
    arrange(line, .by_group = TRUE) |>
    ungroup() |>
    mutate(row_num = dplyr::row_number())
}

#' parse LECS post times
#'
#' @param df a LECS dataframe with added metadata
#'
#' @return a dataframe of post times
#' @export
lecs_post_times <- function(df) {
  df |>
    dplyr::filter(type == "$") |>
    mutate(timestamp = lubridate::ymd_hms(data),
           row_count = row_num - lag(row_num)) |>
    select(timestamp, send, row_count)
}

#' parse RGA data
#'
#' @param df a dataframe with added metadata
#'
#' @return a dataframe of RGA data
#' @export
gems_rga_data <- function(df) {
  df |>
    dplyr::filter(type == "R") |>
    tidyr::separate(data,
                    into = c('timestamp', 'mass', 'current'),
                    sep = ',') |>
    mutate(across(c('mass', 'current'), as.integer),
           timestamp = lubridate::ymd_hms(timestamp),
    )
}

#' parse LECS met data
#'
#' @param df a LECS dataframe with added metadata
#'
#' @return a dataframe of met data
#' @export
lecs_met_data <- function(df) {
  df |>
    dplyr::filter(type == "M") |>
    tidyr::separate(data,
             into = c('hour', 'min', 'sec',
                      'day', 'month', 'year',
                      'PAR', 'wind_speed', 'wind_dir'),
             sep = ',') |>
    mutate(wind_dir = stringr::str_sub(wind_dir, 1, 6),
           across(c('hour', 'min', 'sec',
                      'day', 'month', 'year'), as.integer),
           across(c('PAR', 'wind_speed', 'wind_dir'), as.numeric),
           wind_speed = ifelse(wind_speed < 99, wind_speed, NA),
           wind_dir = ifelse(wind_dir < 360, wind_dir, NA),
           timestamp = lubridate::make_datetime(year, month, day,
                                     hour, min, sec, tz = "UTC"),
           # Convert to UTC. Handle EST/EDT
           timestamp = dplyr::if_else(timestamp > "2023-11-15 11:00:00" &
                                 timestamp < "2024-03-13 00:00:00",
                               timestamp + 5 * 3600,
                               timestamp + 4 * 3600
                               )
                              )|>
    select(-row_num, -type, -any_of("line"))
}

#' Clean met data
#'
#' Final output columns are selected in lecs_parse_file
#'
#' @param met
#'
#' @export
lecs_clean_met <- function(met) {
  met |>
    filter(timestamp > "2023-03-22",
           timestamp < "2024-09-28")
}

#' parse GEMS turbo status data
#'
#' @param df a GEMS dataframe with added metadata
#'
#' @return a dataframe of status data
#' @export
gems_turbo_status <- function(df) {
  df |>
    dplyr::filter(type == "!") |>
    tidyr::separate(data,
                    into = c('timestamp',
                             'status', 'speed', 'power',
                             'voltage', 'e_temp', 'p_temp',
                             'm_temp', 'filament'),
                    sep = ',') |>
    dplyr::mutate(timestamp = lubridate::ymd_hms(timestamp),
                  status = as.integer(status),
                  speed = as.integer(speed), # Hz
                  power = as.integer(power), # Watts
                  voltage = as.integer(voltage) / 100, # Volts
                  e_temp = as.integer(e_temp), # Electronics C
                  p_temp = as.integer(p_temp), # Pump C
                  m_temp = as.integer(m_temp), # Motor C
                  filament = as.numeric(filament)) # Amps
}

#' parse LECS status data
#'
#' @param df a LECS dataframe with added metadata
#'
#' @return a dataframe of status data
#' @export
lecs_status_data <- function(df) {
  df |>
    dplyr::filter(type == "S") |>
    tidyr::separate(data,
                    into = c('timestamp',
                             'adv_min', 'adv_sec', 'adv_day',
                             'adv_hour', 'adv_year', 'adv_month',
                             'bat', 'soundspeed', 'heading', 'pitch',
                             'roll', 'temp'),
                    sep = ',') |>
    mutate(dplyr::across(c('adv_min', 'adv_day',
                           'adv_hour', 'adv_year', 'adv_month'),
                         as.integer),
           dplyr::across(c('adv_sec',
                           'bat', 'soundspeed', 'heading', 'pitch',
                           'roll', 'temp',), as.numeric),
           dplyr::across(c('bat', 'soundspeed',
                           'heading', 'pitch', 'roll') , ~(.x) * .1),
           temp = temp * 0.01,
           timestamp = lubridate::ymd_hms(timestamp),
           adv_timestamp = lubridate::make_datetime(adv_year + 2000, adv_month, adv_day,
                                                    adv_hour, adv_min, adv_sec, tz = "UTC"))
}

#' Clean status data
#'
#' Final output columns are selected in lecs_parse_file
#'
#' @param status
#'
#' @export
lecs_clean_status <- function(status) {
  status |>
    filter(timestamp > "2024-03-22",
           timestamp < "2025-09-28",
           temp < 30,
           temp > 0,
           bat > 9,
           bat < 20,
           soundspeed > 1000,
           soundspeed < 2000,
           heading == 0,
           pitch > -180,
           pitch < 180,
           roll > -180,
           roll < 180,
           adv_day > 0, adv_day < 32,
           adv_month > 0, adv_month < 13,
           adv_min >= 0, adv_min < 61,
           adv_hour >= 0, adv_hour < 24,
           adv_year > 0, adv_year < 100) |>
    mutate(orig_timestamp = timestamp,
           timestamp = correct_status_timestamp_jitter(orig_timestamp, adv_timestamp),
           adv_timestamp_cor = correct_status_timestamp_adv(orig_timestamp, adv_timestamp))
  # Final output columns are selected in lecs_parse_file
}

#' parse LECS ADV data
#'
#' Run time alignment after this
#' Final output columns are selected in lecs_parse_file
#'
#' @param df a LECS dataframe with added metadata
#'
#' @return a dataframe of ADV data
#' @export
#' @importFrom magrittr %>%
lecs_adv_data <- function(df) {
df |>
    dplyr::filter(type == "D") |>
    tidyr::separate(data,
                    into = c('count', 'pressure',
                             'u', 'v', 'w',
                             'amp1', 'amp2', 'amp3',
                             'corr1', 'corr2', 'corr3',
                             'ana_in', 'ana_in2', 'ph_counts'),
                    sep = ',') |>
    mutate(pressure = as.numeric(pressure),
           dplyr::across(c('count',
                             'u', 'v', 'w',
                             'amp1', 'amp2', 'amp3',
                             'corr1', 'corr2', 'corr3',
                           'ana_in', 'ana_in2', 'ph_counts'),
                         as.integer),
           # scale velocity to m/s. 1 count = 0.001 m/s for high range
           # 1 count = 0.0001 m/s for low range
           # should define this elsewhere!
           dplyr::across(c('u', 'v', 'w'), ~(.x) * .0001),
           pressure = pressure / 1000,
    )
}

#' Clean adv data
#'
#' Final output columns are selected in lecs_parse_file
#'
#' @param adv
#'
#' @export
lecs_clean_adv_data <- function(adv) {
  adv |>
    #tidyr::drop_na() |>
    filter(
           !is.na(count),
           count >= 0,
           count < 256,
           pressure > 0.4,
           pressure < 2.6,
           temp > 1.6,
           temp < 26,
           corr1 <= 100,
           corr2 <= 100,
           corr3 <= 100,
           !is.na(u),
           !is.na(v),
           !is.na(w),
           !is.na(corr1),
           !is.na(corr2),
           !is.na(corr3),
           !is.na(amp1),
           !is.na(amp2),
           !is.na(amp3),
           # filter velocities to remove impossible velocities
           u > -20, u < 20,
           v > -20, v < 20,
           w > -4, w < 4) |>
    #Replace outliers with NA with rolling Hampel filter
    # use purrr::possibly to catch findOutliers errors
    # and replace with original vector
    # TODO: change this to flag outliers
    mutate(across(c(pressure, u, v, w, amp1, amp2, amp3,
                  corr1, corr2, corr3, ph_counts, temp, DO, DO_percent, ph_counts, ox_umol_l, pH),
                  \(x) replace(x, purrr::possibly(seismicRoll::findOutliers, NULL)(x), NA)))
  # Final output columns are selected in lecs_parse_file
}

#' Calculate number of missing lines
#'
#' @param count Vector of ADV count
#' @param line vector indicating line of send
#'
#' @return A vector with a guess at the number of missing lines from
#' previous to current
#' @export
lecs_missing <- function(count, line = NULL) {
  missing <-  dplyr::case_when(count > 255 | dplyr::lag(count) > 255 ~
                                 NA_integer_,
                               count > dplyr::lag(count) ~
                                 count - 1L - dplyr::lag(count),
                               TRUE ~
                                 255L + count - dplyr::lag(count))
  if (!is.null(line)) {
    missing = replace(missing, missing > 0 & lag(line) != line - 1L, NA)
  }
  missing
}

#' Read seaphox data
#'
#' @param file seaphox file path
#'
#' @return a dataframe of seaphox data
#' @export
read_seaphox <- function(file) {
  data.table::fread(file) |>
    janitor::clean_names() |>
    mutate(timestamp = lubridate::mdy_hms(date_time_utc_04_00,
                          tz = "UTC")) |>
    select(timestamp,
           pH = internal_p_h_p_h,
           temp = p_h_temperature_celsius,
           pressure = pressure_decibar,
           sal = salinity_psu,
           oxygen = oxygen_ml_l)
}

#' Read ProOceanus CO2 files
#'
#' @param file
#'
#' @return a data frame of CO2 data
#' @export
read_prooceanus <- function(file) {
  df_names <- c("type", "year", "month", "day", "hour", "minute", "second",
                "zero_ad", "cur_ad", "co2", "irga_temp", "humidity", "hum_sens_temp",
                "cell_pressure", "bat_v", "board_temp", "ana_in_1", "ana_in_2",
                "d_in_1", "d_in_2")
  df_spec <- readr::cols(
    type = 'c',
    year = 'i',
    month = 'i',
    day = 'i',
    hour = 'i',
    minute = 'i',
    second = 'i',
    zero_ad = 'i',
    cur_ad = 'i',
    co2 = 'd',
    irga_temp = 'd',
    humidity = 'd',
    hum_sens_temp = 'd',
    cell_pressure = 'd',
    bat_v = 'd',
    board_temp = 'd',
    ana_in_1 = 'd',
    ana_in_2 = 'd',
    d_in_1 = 'd',
    d_in_2 = 'c'
  )
  #find line with "File Contents:", next line is header
  skipped_lines <- readLines(file) |>
    stringr::str_detect("File Contents:") |>
    which() + 1
  readr::read_csv(file, col_names = df_names, col_types = df_spec, skip = skipped_lines, comment = "%") |>
    mutate(dplyr::across(c(month, day, hour, minute, second), as.integer),
           ts = lubridate::make_datetime(year, month, day, hour, minute, second))
}

#' Read Underwater PAR Odyssey files
#'
#' @param file
#'
#' @return a data frame of PAR data
#' @export
lecs_read_par_odyssey <- function(file) {
  header <- c("scan", "date", "time", "raw_par", "cal_par")
  read_csv(file, skip = 9, col_names = header) |>
    mutate(date = as.Date(date, "%d/%m/%Y"),
           timestamp = lubridate::parse_date_time(paste(date, time),
                                                  "Ymd HMS",
                                                  tz = "UTC")) |>
    select(scan, timestamp, raw_par, cal_par)
}
