# Script to read and combine CSV files from a directory
# Uses tidyverse for data manipulation

library(tidyverse)
library(arrow)

# Define the directory containing CSV files
csv_dir <- "data/raw/GEMS_VALVE/"

# Get list of all CSV files in the directory
csv_files <- list.files(
  path = csv_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

# Function to process each file
process_file <- function(file_path) {
  file_path
  # Read all lines from the file
  all_lines <- readLines(file_path)

  # Extract lines starting with "Rebooted"
  reboot_lines <- all_lines[grepl("^Rebooted", all_lines)]

  # Get lines that don't start with "Rebooted" (CSV data)
  csv_lines <- all_lines[!grepl("^Rebooted", all_lines)]

  # Parse CSV data
  csv_data <- read_csv(
    I(csv_lines),
    col_types = cols(
      timestamp = col_datetime(),
      voltage = col_double(),
      current = col_double(),
      valve_position = col_double()
    )
  ) %>%
    mutate(
      source_file = basename(file_path),
      reboot_info = paste(reboot_lines, collapse = "; ")
    )

  return(csv_data)
}

# Read and combine all CSV files
combined_data <- csv_files %>%
  map_dfr(process_file)

# View the combined data
glimpse(combined_data)

# Optional: Write the combined data to a new CSV file
write_csv(combined_data, "data/processed/valve_combined.csv")

write_parquet(combined_data, "data/processed/valve_combined.parquet")

# Print summary
cat("Successfully combined", length(csv_files), "CSV files\n")
cat("Total rows:", nrow(combined_data), "\n")
cat("Total columns:", ncol(combined_data), "\n")
