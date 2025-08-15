## code to prepare `DATASET` dataset goes here

# Create package data objects for ZooMSS
# This script creates the .rda files in the data/ directory

# Load GroupInputs from CSV
GroupInputs <- readr::read_csv("data-raw/GroupInputs.csv", show_col_types = FALSE)


# Save the data
usethis::use_data(GroupInputs, overwrite = TRUE)

