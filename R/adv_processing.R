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
