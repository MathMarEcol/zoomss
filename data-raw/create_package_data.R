# Create package data objects for ZooMSS
# This script creates the .rda files in the data/ directory

# Load TestGroups from CSV
TestGroups <- readr::read_csv("TestGroups.csv", show_col_types = FALSE)

# Ensure data/ directory exists
if (!dir.exists("../data")) {
  dir.create("../data")
}

# Save as package data
save(TestGroups, file = "../data/TestGroups.rda")

cat("✓ Created TestGroups.rda with", nrow(TestGroups), "functional groups\n")
cat("✓ Package data ready for use\n")
