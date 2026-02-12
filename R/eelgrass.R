# Functions for working with eelgrass quadrat data

read_eelgrass_density <- function(eelgrass_file) {
    df <- readxl::read_excel(
        eelgrass_file,
        sheet = 1,
        col_names = c("date", "type", "density_g_m2", "density_g_m2_sd"),
        col_types = c("text", "text", "numeric", "numeric"),
        skip = 1
    ) |>
        dplyr::mutate(date = as.Date(date, format = "%m/%d/%Y"))
}

read_eelgrass_length <- function(eelgrass_file) {
    df <- readxl::read_excel(
        eelgrass_file,
        sheet = 2,
        col_names = c(
            "date",
            "length_cm",
            "length_cm_sd",
            "shoot_density_m2",
            "shoot_density_m2_sd"
        ),
        col_types = c("date", "numeric", "numeric", "numeric", "numeric"),
        skip = 1
    )
}

combine_eelgrass_data <- function(density_df, length_df) {
    combined_df <- density_df |>
        dplyr::filter(type == "Grass") |>
        dplyr::select(-type) |>
        dplyr::inner_join(length_df, by = "date") |>
        dplyr::mutate(
            day_of_year = yday(date),
            year = as.factor(lubridate::year(date)),
            biomass_per_shoot_g = density_g_m2 / shoot_density_m2,
            biomass_per_shoot_g_sd = sqrt(
                (density_g_m2_sd / shoot_density_m2)^2 +
                    (density_g_m2 * length_cm_sd / shoot_density_m2^2)^2
            )
        )
    return(combined_df)
}

read_2021_eelgrass <- function(sg_length_2021_file, sg_weight_2021_file) {
    sg_length_2021 <- read_csv(sg_length_2021_file) |>
        janitor::clean_names() |>
        mutate(date = as.POSIXct(date, format = "%m/%d/%Y")) |>
        group_by(date) |>
        summarize(length_cm = mean(l), length_cm_sd = sd(l))
    sg_weight_2021 <- read_csv(sg_weight_2021_file) |>
        janitor::clean_names() |>
        mutate(
            date = as.POSIXct(date, format = "%m/%d/%Y"),
            density_g_m2 = sea_grass_weight / 0.04
        ) |>
        group_by(date) |>
        summarize(
            density_g_m2_sd = sd(density_g_m2),
            density_g_m2 = mean(density_g_m2)
        )
    combined_2021 <- sg_length_2021 |>
        inner_join(sg_weight_2021, by = "date") |>
        mutate(
            day_of_year = yday(date),
            year = as.factor(lubridate::year(date)),
            biomass_per_shoot_g = NA_real_,
            biomass_per_shoot_g_sd = NA_real_,
            shoot_density_m2 = NA_real_,
            shoot_density_m2_sd = NA_real_
        ) |>
        select(
            date,
            density_g_m2,
            density_g_m2_sd,
            length_cm,
            length_cm_sd,
            shoot_density_m2,
            shoot_density_m2_sd,
            day_of_year,
            year,
            biomass_per_shoot_g,
            biomass_per_shoot_g_sd
        )
    return(combined_2021)
}


get_eelgrass_data <- function(
    eelgrass_file,
    sg_length_2021_file,
    sg_weight_2021_file
) {
    density_df <- read_eelgrass_density(eelgrass_file)
    length_df <- read_eelgrass_length(eelgrass_file)
    combined_df <- combine_eelgrass_data(density_df, length_df)
    eelgrass_2021 <- read_2021_eelgrass(
        sg_length_2021_file,
        sg_weight_2021_file
    )
    combined_df <- bind_rows(combined_df, eelgrass_2021)
    combined_df$year <- factor(
        combined_df$year,
        levels = sort(unique(combined_df$year))
    )
    return(combined_df)
}
