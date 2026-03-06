# Plots for GEMS Data

# RGA mass plot for July data
plot_rga_masses <- function(rga_binned) {
  jul_rga <- rga_binned |>
    filter(
      timestamp > as.POSIXct("2025-07-15 00:00:00"),
      timestamp < as.POSIXct("2025-07-20 00:00:00")
    )

  mep <- jul_rga |>
    ggplot(aes(timestamp, mass_15_40, color = inlet)) +
    geom_line(
      alpha = 0.6
    ) +
    labs(
      title = "Methane (mass 15)",
      x = NULL,
      y = "Mass 15:40"
    )

  nip <- jul_rga |>
    ggplot(aes(timestamp, mass_28_40, color = inlet)) +
    geom_line(
      alpha = 0.6
    ) +
    labs(
      title = "Nitrogen (mass 28)",
      x = NULL,
      y = "Mass 28:40"
    ) +
    theme(legend.position = "none")

  sup <- jul_rga |>
    ggplot(aes(timestamp, mass_34_40, color = inlet)) +
    geom_line(
      alpha = 0.6
    ) +
    labs(
      title = "Hydrogen Sulfide (mass 34)",
      x = NULL,
      y = "Mass 34:40"
    ) +
    theme(legend.position = "none")

  cop <- jul_rga |>
    ggplot(aes(timestamp, mass_44_40, color = inlet)) +
    geom_line(
      alpha = 0.6
    ) +
    labs(
      title = "Carbon Dioxide (mass 44)",
      x = NULL,
      y = "Mass 44:40"
    ) +
    theme(legend.position = "none")

  oxp <- rga_binned |>
    filter(
      timestamp > as.POSIXct("2025-07-15 00:00:00"),
      timestamp < as.POSIXct("2025-07-20 00:00:00")
    ) |>
    ggplot(aes(timestamp, mass_32_40, color = inlet)) +
    geom_line(
      alpha = 0.6
    ) +
    labs(
      x = NULL,
      y = "Mass 32:40"
    )

  op <- oxp +
    ggtitle("Oxygen (mass 32)") +
    theme(legend.position = "none")

  mep / nip / op / sup / cop
}

# QMS system noise comparison plot
plot_qms_noise_comparison <- function(rga_binned, gems_2022_file) {
  # Jeff's processed 2022 data. Unclear what inlet this is from.
  mat <- readMat(here::here(gems_2022_file))

  # Adjust these names to match your .mat structure:
  # e.g. mat$date, mat$data_D, etc.
  date_matlab <- as.vector(mat$date) # numeric MATLAB datenums
  data_D <- as.matrix(mat$data.D) # numeric matrix

  # Helper: convert MATLAB datenum -> POSIXct (UTC)
  matlab_to_posix <- function(datenum) {
    # MATLAB datenum origin: 0000-01-01, R: 1970-01-01
    # Difference in days:
    # as.numeric(as.Date("1970-01-01") - as.Date("0000-01-01")) == 719529
    origin_shift <- 719529
    as.POSIXct(
      (datenum - origin_shift) * 86400,
      origin = "1970-01-01",
      tz = "UTC"
    )
  }

  date <- matlab_to_posix(date_matlab)

  # Apply the MATLAB filters:

  # 1) id=find(data_D(:,10)<1.2*10^7);
  keep_idx1 <- data_D[, 10] >= 1.2 * 10^7
  data_D <- data_D[keep_idx1, ]
  date <- date[keep_idx1]

  # 2) id=find(date<datenum(2022,8,10,14,30,0) | date>datenum(2022,8,31,7,0,0));

  start_matlab <- as.numeric(as.POSIXct("2022-08-10 14:30:00", tz = "UTC")) /
    86400 +
    719529
  end_matlab <- as.numeric(as.POSIXct("2022-08-31 07:00:00", tz = "UTC")) /
    86400 +
    719529

  start_r <- matlab_to_posix(start_matlab)
  end_r <- matlab_to_posix(end_matlab)

  keep_idx2 <- date >= start_r & date <= end_r
  data_D <- data_D[keep_idx2, ]
  date <- date[keep_idx2]

  # Build a data.frame for ggplot
  df <- data.frame(
    date = date,
    ratio_O2_Ar = data_D[, 3] / data_D[, 5] # data_D(:,3)./data_D(:,5)
  )

  jp <- df |>
    filter(date > "2022-08-14", date < "2022-08-18") |>
    ggplot(aes(date, ratio_O2_Ar)) +
    geom_line(color = cb_print_4[1], linewidth = 0.2, alpha = 0.6) +
    labs(x = "Date (2022)", y = "Mass 32:40")

  bp <- rga_binned |>
    filter(
      inlet == "high",
      timestamp > "2025-08-14",
      timestamp < "2025-08-18"
    ) |>
    ggplot(aes(timestamp, mass_32_40)) +
    geom_line(color = cb_print_4[1], linewidth = 0.2, alpha = 0.6) +
    labs(x = "Date (2025)", y = "Mass 32:40")

  (jp / bp) + plot_annotation(tag_levels = 'A')
}

# Argon normalization plot
plot_argon_normalization <- function(rga_binned) {
  sep_rga <- rga_binned |>
    filter(
      timestamp > as.POSIXct("2025-09-04 00:00:00"),
      timestamp < as.POSIXct("2025-09-14 00:00:00")
    )

  nup <- sep_rga |>
    ggplot(aes(timestamp, mass_28, color = inlet)) +
    geom_line(
      alpha = 0.6
    ) +
    labs(
      x = NULL,
      y = "Mass 28 [Torr]"
    )

  nnp <- sep_rga |>
    ggplot(aes(timestamp, mass_28_40, color = inlet)) +
    geom_line(
      alpha = 0.6
    ) +
    labs(
      x = NULL,
      y = "Mass 28:40"
    ) +
    theme(legend.position = "none")

  (nup / nnp) + plot_annotation(tag_levels = 'A')
}

# Oxygen timeseries
plot_oxygen_timeseries <- function(rga_calibrated) {
  rga_calibrated |>
    arrange(timestamp) |>
    select(timestamp, mass_32_40_high, mass_32_40_low) |>
    mutate(
      time_diff = as.numeric(difftime(
        timestamp,
        lag(timestamp),
        units = "hours"
      )),
      mass_32_40_high = ifelse(time_diff > 1, NA, mass_32_40_high),
      mass_32_40_low = ifelse(time_diff > 1, NA, mass_32_40_low)
    ) |>
    pivot_longer(
      cols = c(mass_32_40_low, mass_32_40_high),
      names_to = "inlet",
      names_prefix = "mass_32_40_",
      values_to = "mass_32_40"
    ) |>
    ggplot(aes(timestamp, mass_32_40, color = inlet)) +
    geom_line(na.rm = TRUE, alpha = 0.6) +
    labs(x = NULL, y = "Mass 32:40")
}

# Short oxygen timeseries for July
plot_oxygen_timeseries_july <- function(rga_calibrated) {
  rga_calibrated |>
    filter(
      timestamp > as.POSIXct("2025-07-15 00:00:00"),
      timestamp < as.POSIXct("2025-07-20 00:00:00")
    ) |>
    arrange(timestamp) |>
    select(timestamp, mass_32_40_high, mass_32_40_low) |>
    mutate(
      time_diff = as.numeric(difftime(
        timestamp,
        lag(timestamp),
        units = "hours"
      )),
      mass_32_40_high = ifelse(time_diff > 1, NA, mass_32_40_high),
      mass_32_40_low = ifelse(time_diff > 1, NA, mass_32_40_low)
    ) |>
    pivot_longer(
      cols = c(mass_32_40_low, mass_32_40_high),
      names_to = "inlet",
      names_prefix = "mass_32_40_",
      values_to = "mass_32_40"
    ) |>
    ggplot(aes(timestamp, mass_32_40, color = inlet)) +
    geom_line(na.rm = TRUE, alpha = 0.6) +
    labs(x = NULL, y = "Mass 32:40")
}

# Optode vs normalized mass 32
plot_oxygen_calibration <- function(ox_cal_df, ox_model) {
  om <- ox_model$coefficients[2]
  oi <- ox_model$coefficients[1]
  r2 <- summary(ox_model)$r.squared

  cal_plot_df <- ox_cal_df |>
    filter(timestamp > "2025-07-11 10:00:00") |>
    mutate(
      rga_ox = mass_32_40 * om + oi
    )

  eq_label <- paste0(
    "y = ",
    round(oi, 2),
    " + ",
    round(om, 2),
    "x\n",
    "R² = ",
    round(r2, 3)
  )

  ggplot(cal_plot_df, aes(mass_32_40, seaphox_oxygen_ml_l)) +
    geom_smooth(method = "lm", se = FALSE, color = "darkgrey") +
    annotate(
      "text",
      size = 3,
      x = Inf,
      y = -Inf,
      label = eq_label,
      hjust = 1.1,
      vjust = -1
    ) +
    geom_point(
      shape = 1,
      color = "blue",
      size = 2,
      stroke = 0.5,
      alpha = 0.2
    ) +
    #xlim(3, 9) +
    #ylim(3, 9) +
    labs(
      y = "Optode oxygen [ml/l]",
      x = "GEMS mass 32:40"
    )
}

# Fitted Oxygen
plot_cal_timeseries <- function(cal_plot_df, ox_model) {
  om <- ox_model$coefficients[2]
  oi <- ox_model$coefficients[1]
  r2 <- summary(ox_model)$r.squared
  cal_plot_df <- ox_cal_df |>
    filter(timestamp > "2025-07-11 10:00:00") |>
    mutate(
      rga_ox = mass_32_40 * om + oi
    )

  cal_plot_df |>
    pivot_longer(
      cols = c(seaphox_oxygen_ml_l, rga_ox),
      names_to = "source",
      values_to = "oxygen_ml_l"
    ) |>
    mutate(source = ifelse(source == "rga_ox", "GEMS", "Optode")) |>
    ggplot(aes(timestamp, oxygen_ml_l, color = source)) +
    geom_line(
      alpha = 0.6
    ) +
    labs(x = NULL, y = "Oxygen [ml/l]")
}

# Diel Flux
plot_diel_flux <- function(hourly_flux) {
  gop <- hourly_flux |>
    mutate(
      time_diff = as.numeric(difftime(
        timestamp,
        lag(timestamp),
        units = "hours"
      )),
      ox_flux = ifelse(time_diff > 1, NA, ox_flux),
    ) |>
    ggplot(aes(timestamp, ox_flux)) +
    geom_hline(yintercept = 0) +
    geom_line(color = cb_print_4[1], na.rm = TRUE) +
    ylim(-15, 15) +
    labs(
      x = NULL,
      y = expression("Oxygen flux [" * mmol ~ m^-2 ~ h^-1 * "]")
    )

  gcp = ggplot(hourly_flux, aes(timestamp, co2_flux)) +
    geom_line()

  tp <- hourly_flux |>
    mutate(
      time_diff = as.numeric(difftime(
        timestamp,
        lag(timestamp),
        units = "hours"
      )),
      adv_temp_mean = ifelse(time_diff > 1, NA, adv_temp_mean),
    ) |>
    ggplot(aes(timestamp, adv_temp_mean)) +
    geom_line(color = cb_print_4[1]) +
    labs(x = NULL, y = "Temp [C]")

  gpp = ggplot(hourly_flux, aes(timestamp, par_mean)) +
    geom_line() +
    labs(
      x = NULL,
      y = expression(
        PAR ~ "[" *
          mu *
          "mol " *
          m^{
            -2
          } *
          " " *
          s^{
            -1
          } *
          "]"
      )
    )
  # DLI daily PAR
  min_flux_ts <- as.Date(min(hourly_flux[["timestamp"]]))
  max_flux_ts <- as.Date(max(hourly_flux[["timestamp"]]))

  par_model_df <- par_model_df |>
    filter(date >= min_flux_ts, date <= max_flux_ts)

  dpp <- ggplot(par_model_df, aes(date, dli_mol_m2_day)) +
    geom_line(color = cb_print_4[1]) +
    labs(
      x = NULL,
      y = expression(
        DLI ~ "[mol " *
          m^{
            -2
          } *
          " " *
          {
            day
          }^{
            -1
          } *
          "]"
      )
    )

  (gop / dpp / tp) +
    plot_layout(heights = c(3, 2, 2)) +
    plot_annotation(tag_levels = 'A')
}

# Representative daily flux
plot_rep_daily_flux <- function(hourly_flux) {
  jul_flux <- hourly_flux |>
    filter(
      timestamp > as.POSIXct("2025-07-15 00:00:00"),
      timestamp < as.POSIXct("2025-07-19 12:00:00")
    )
  jfp <- jul_flux |>
    ggplot(aes(timestamp, ox_flux)) +
    geom_hline(yintercept = 0) +
    geom_line(color = cb_print_4[1]) +
    labs(
      x = NULL,
      y = expression("Oxygen flux [" * mmol ~ m^-2 ~ h^-1 * "]")
    )

  jpp <- jul_flux |>
    ggplot(aes(timestamp, par_mean)) +
    geom_line(color = cb_print_4[1]) +
    labs(
      x = NULL,
      y = expression(
        PAR ~ "[" *
          mu *
          "mol " *
          m^{
            -2
          } *
          " " *
          s^{
            -1
          } *
          "]"
      )
    )
  (jfp / jpp) +
    plot_layout(heights = c(3, 1)) +
    plot_annotation(tag_levels = 'A')
}

# Diel Concentration and gradient
plot_grad_diel <- function(hourly_stats) {
  hmp <- hourly_stats |>
    ggplot(aes(solar_hour, ox_mean_umol_l_mean_mean)) +
    geom_line() +
    geom_pointrange(
      aes(
        ymin = ox_mean_umol_l_mean_mean - ox_mean_umol_l_mean_se,
        ymax = ox_mean_umol_l_mean_mean + ox_mean_umol_l_mean_se
      ),
      size = 0.2,
      color = cb_print_4[1]
    ) +
    labs(x = NULL, y = expression("Oxygen [" * mmol ~ m^-3 * "]"))

  hgp <- hourly_stats |>
    ggplot(aes(solar_hour, ox_gradient_umol_l_m_mean_mean)) +
    geom_line() +
    geom_pointrange(
      aes(
        ymin = ox_gradient_umol_l_m_mean_mean - ox_gradient_umol_l_m_mean_se,
        ymax = ox_gradient_umol_l_m_mean_mean + ox_gradient_umol_l_m_mean_se
      ),
      size = 0.2,
      color = cb_print_4[1]
    ) +
    geom_hline(yintercept = 0) +
    labs(
      x = "Hour of Day",
      y = expression("Oxygen Gradient [" * mmol ~ m^-4 * "]")
    )

  (hmp / hgp) + plot_annotation(tag_levels = 'A')
}

# Oxygen flux and predicted par
plot_flux_par <- function(hourly_stats) {
  flp <- hourly_stats |>
    ggplot(aes(solar_hour, ox_flux_mean)) +
    geom_line() +
    geom_pointrange(
      aes(
        ymin = ox_flux_mean - ox_flux_se,
        ymax = ox_flux_mean + ox_flux_se
      ),
      size = 0.2,
      color = cb_print_4[1]
    ) +
    geom_hline(yintercept = 0) +
    labs(x = NULL, y = expression("Oxygen flux [" * mmol ~ m^-2 ~ h^-1 * "]"))

  parp <- hourly_stats |>
    ggplot(aes(solar_hour, par_mean_mean)) +
    geom_line(color = cb_print_4[1]) +
    labs(
      x = "Hour of Day",
      y = expression(
        PAR ~ "[" *
          mu *
          "mol " *
          m^{
            -2
          } *
          " " *
          s^{
            -1
          } *
          "]"
      )
    )

  (flp / parp) +
    plot_layout(heights = c(3, 1)) +
    plot_annotation(tag_levels = 'A')
}

# CO2 concentration, gradient and flux
plot_co2_flux <- function(hourly_stats) {
  hmp <- hourly_stats |>
    ggplot(aes(solar_hour, co2_mean_umol_l_mean_mean)) +
    geom_line() +
    geom_pointrange(
      aes(
        ymin = co2_mean_umol_l_mean_mean - co2_mean_umol_l_mean_se,
        ymax = co2_mean_umol_l_mean_mean + co2_mean_umol_l_mean_se
      ),
      size = 0.2,
      color = cb_print_4[1]
    ) +
    labs(x = NULL, y = expression("CO"[2] ~ "[" * mmol ~ m^-3 * "]"))

  hgp <- hourly_stats |>
    ggplot(aes(solar_hour, co2_gradient_umol_l_m_mean_mean)) +
    geom_line() +
    geom_pointrange(
      aes(
        ymin = co2_gradient_umol_l_m_mean_mean - co2_gradient_umol_l_m_mean_se,
        ymax = co2_gradient_umol_l_m_mean_mean + co2_gradient_umol_l_m_mean_se
      ),
      size = 0.2,
      color = cb_print_4[1]
    ) +
    geom_hline(yintercept = 0) +
    labs(
      x = NULL,
      y = expression("CO"[2] ~ "Gradient [" * mmol ~ m^-4 * "]")
    )

  flp <- hourly_stats |>
    ggplot(aes(solar_hour, co2_flux_mean)) +
    geom_line() +
    geom_pointrange(
      aes(
        ymin = co2_flux_mean - co2_flux_se,
        ymax = co2_flux_mean + co2_flux_se
      ),
      size = 0.2,
      color = cb_print_4[1]
    ) +
    geom_hline(yintercept = 0) +
    labs(
      x = NULL,
      y = expression("CO"[2] ~ "flux [" * mmol ~ m^-2 ~ h^-1 * "]")
    )

  parp <- hourly_stats |>
    ggplot(aes(solar_hour, par_mean_mean)) +
    geom_line(color = cb_print_4[1]) +
    labs(
      x = "Hour of Day",
      y = expression(
        PAR ~ "[" *
          mu *
          "mol " *
          m^{
            -2
          } *
          " " *
          s^{
            -1
          } *
          "]"
      )
    )

  (hmp / hgp / flp / parp) +
    plot_layout(heights = c(3, 3, 3, 1)) +
    plot_annotation(tag_levels = 'A')
}

# CO2 vs O2 gradients and fluxes
plot_co2_vs_o2 <- function(rga_calibrated, hourly_flux) {
  # Function to map hour (0-24) to color: dark blue at night, yellow in day
  hour_to_color <- function(hours) {
    # Normalize hours into [0,24)
    hours <- hours %% 24

    night_col <- grDevices::col2rgb("#1a4eb5ff") # dark blue
    day_col <- grDevices::col2rgb("#fce23aff") # yellow

    # Initialize RGB matrix
    rgb_mat <- matrix(NA, nrow = 3, ncol = length(hours))

    for (i in seq_along(hours)) {
      h <- hours[i]

      if (h <= 6) {
        # Deep night: constant dark blue
        mix <- night_col
      } else if (h > 6 && h < 8) {
        # Dawn: interpolate from blue -> yellow between 5 and 7
        t <- (h - 6) / (8 - 6) # 0 at 5, 1 at 7
        mix <- (1 - t) * night_col + t * day_col
      } else if (h >= 8 && h <= 17) {
        # Full day: constant yellow
        mix <- day_col
      } else if (h > 17 && h < 19) {
        # Dusk: interpolate from yellow -> blue between 17 and 19
        t <- (h - 17) / (19 - 17) # 0 at 17, 1 at 19
        mix <- (1 - t) * day_col + t * night_col
      } else {
        # Night again: constant dark blue
        mix <- night_col
      }

      rgb_mat[, i] <- mix
    }

    grDevices::rgb(
      red = rgb_mat[1, ] / 255,
      green = rgb_mat[2, ] / 255,
      blue = rgb_mat[3, ] / 255
    )
  }

  # A convenience palette generator for gradientn
  hour_palette <- function(n = 256) {
    # Sample hours evenly across 0-24
    hours_seq <- seq(0, 24, length.out = n)
    hour_to_color(hours_seq)
  }

  gp <- rga_calibrated |>
    mutate(hour = (hour(timestamp) - 4) %% 24) |>
    ggplot(aes(ox_gradient_umol_l_m, co2_gradient_umol_l_m, color = hour)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(size = 1, alpha = .6) +
    scale_x_reverse() +
    scale_y_reverse() +
    scale_colour_gradientn(
      colours = hour_palette(256),
      limits = c(0, 24),
      breaks = c(0, 6, 12, 18, 24),
      labels = c("00:00", "06:00", "12:00", "18:00", "24:00"),
      name = "Hour of day"
    ) +
    labs(
      x = expression("O"[2] ~ "Gradient [" * mmol ~ m^-4 * "]"),
      y = expression("CO"[2] ~ "Gradient [" * mmol ~ m^-4 * "]"),
      color = "PAR"
    ) +
    theme(legend.position = "none")

  fp <- hourly_flux |>
    filter(ox_flux < 30) |>
    mutate(hour = (hour(timestamp) - 4) %% 24) |>
    ggplot(aes(ox_flux, co2_flux, color = hour)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(size = 1, alpha = .5) +
    scale_colour_gradientn(
      colours = hour_palette(256),
      limits = c(0, 24),
      breaks = c(0, 6, 12, 18, 24),
      labels = c("00:00", "06:00", "12:00", "18:00", "24:00"),
      name = "Hour of day"
    ) +
    labs(
      x = expression("O"[2] ~ "flux [" * mmol ~ m^-2 ~ h^-1 * "]"),
      y = expression("CO"[2] ~ "flux [" * mmol ~ m^-2 ~ h^-1 * "]"),
      color = "Hour of Day"
    )

  (gp + fp) + plot_annotation(tag_levels = 'A')
}

# diel fluxes by month
plot_diel_monthly <- function(monthly_stats) {
  monthly_stats <- mutate(monthly_stats, month = as.factor(month))
  flp <- monthly_stats |>
    ggplot(aes(solar_hour, ox_flux_mean, color = month)) +
    geom_line() +
    geom_pointrange(
      aes(
        ymin = ox_flux_mean - ox_flux_se,
        ymax = ox_flux_mean + ox_flux_se
      ),
      size = 0.2,
      # position = position_dodge(width = 0.5)
    ) +
    geom_hline(yintercept = 0) +
    labs(x = NULL, y = expression("Oxygen flux [" * mmol ~ m^-2 ~ h^-1 * "]"))

  parp <- monthly_stats |>
    ggplot(aes(solar_hour, par_mean_mean, color = month)) +
    geom_line() +
    labs(
      x = "Hour of Day",
      y = expression(
        PAR ~ "[" *
          mu *
          "mol " *
          m^{
            -2
          } *
          " " *
          s^{
            -1
          } *
          "]"
      )
    )

  (flp / parp) +
    plot_layout(heights = c(3, 1)) +
    plot_annotation(tag_levels = 'A')
}

# Net ecosystem metabolism
plot_nem <- function(monthly_nem, daily_nem) {
  monthly_nem$month_factor <- factor(
    monthly_nem$month,
    levels = 1:12,
    labels = month.abb # or month.abb
  )
  mp <- ggplot(monthly_nem, aes(month_factor, nem_mmol_m2_day)) +
    geom_col(fill = cb_print_4[1]) +
    geom_hline(yintercept = 0) +
    labs(
      x = NULL,
      y = expression("NEM [" * mmol ~ m^-2 ~ d^-1 * "]")
    )

  dp <- ggplot(daily_nem, aes(day, nem_mmol_m2_h)) +
    geom_col(fill = cb_print_4[1]) +
    geom_hline(yintercept = 0) +

    labs(
      y = expression("NEM [" * mmol ~ m^-2 ~ h^-1 * "]")
    )

  (mp / dp) +
    plot_annotation(tag_levels = "A")
}

# Eelgrass growth
plot_eelgrass <- function(eelgrass) {
  bmp <- ggplot(eelgrass, aes(day_of_year, density_g_m2, color = year)) +
    geom_smooth(
      inherit.aes = FALSE,
      aes(day_of_year, density_g_m2),
      se = FALSE
    ) +
    geom_line() +
    geom_pointrange(
      aes(
        ymin = density_g_m2 - density_g_m2_sd,
        ymax = density_g_m2 + density_g_m2_sd
      ),
      size = 0.2
    ) +
    scale_x_continuous(
      breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
      labels = month.abb
    ) +
    labs(
      x = NULL,
      y = expression("Density [" * g ~ m^-2 * "]")
    )

  blp <- ggplot(eelgrass, aes(day_of_year, length_cm, color = year)) +
    geom_smooth(inherit.aes = FALSE, aes(day_of_year, length_cm), se = FALSE) +
    geom_line() +
    geom_line() +
    geom_pointrange(
      aes(
        ymin = length_cm - length_cm_sd,
        ymax = length_cm + length_cm_sd
      ),
      size = 0.2
    ) +
    scale_x_continuous(
      breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
      labels = month.abb
    ) +
    labs(x = NULL, y = "Length [cm]")

  (bmp / blp) + plot_annotation(tag_levels = 'A')
}
