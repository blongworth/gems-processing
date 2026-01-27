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
            biomass_per_shoot_g = density_g_m2 / shoot_density_m2,
            biomass_per_shoot_g_sd = sqrt(
                (density_g_m2_sd / shoot_density_m2)^2 +
                    (density_g_m2 * length_cm_sd / shoot_density_m2^2)^2
            )
        )
    return(combined_df)
}

get_eelgrass_data <- function(eelgrass_file) {
    density_df <- read_eelgrass_density(eelgrass_file)
    length_df <- read_eelgrass_length(eelgrass_file)
    combined_df <- combine_eelgrass_data(density_df, length_df)
    return(combined_df)
}
