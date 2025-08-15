## Test script for reproduction mechanism in ZooMSS
## This script tests the new reproduction functionality

# Clear workspace
rm(list = ls())

# Set working directory to ZooMSS folder
setwd("c:/Users/kjmurphy/OneDrive - University of Tasmania/Documents/GitHub/ZooMSS")

# Source required functions
source("fZooMSS_Model.R")

# Load test group parameters
Groups <- read.csv("TestGroups.csv")

# Test 1: Run model WITHOUT reproduction (traditional pinning)
cat("Testing model WITHOUT reproduction (traditional pinning)...\n")
input_params_no_repro <- list(
  tmax = 10,           # 10 years
  dt = 0.01,           # time step
  phyto_int = -3,      # phytoplankton intercept
  phyto_slope = -1,    # phytoplankton slope
  phyto_max = -8,      # maximum phytoplankton size
  sst = 15             # sea surface temperature
)

# Ensure reproduction is off for all groups
Groups_no_repro <- Groups
Groups_no_repro$Reproduction <- 0

# Run model without reproduction
results_no_repro <- fZooMSS_Model(input_params_no_repro, Groups_no_repro, SaveTimeSteps = TRUE)

cat("Model without reproduction completed successfully!\n")
cat("Final fish abundances (smallest size classes):\n")
fish_groups <- which(Groups_no_repro$Type == "Fish")
for (i in fish_groups) {
  min_size_idx <- which(round(log10(results_no_repro$model$param$w), digits = 2) == Groups_no_repro$W0[i])[1]
  final_abundance <- results_no_repro$abundances[i, min_size_idx]
  cat(sprintf("  %s: %.2e\n", Groups_no_repro$Species[i], final_abundance))
}

# Test 2: Run model WITH reproduction
cat("\nTesting model WITH reproduction...\n")
input_params_repro <- input_params_no_repro

# Ensure reproduction is on for fish groups only
Groups_repro <- Groups
Groups_repro$Reproduction[Groups_repro$Type == "Fish"] <- 1
Groups_repro$Reproduction[Groups_repro$Type == "Zooplankton"] <- 0

# Run model with reproduction
results_repro <- fZooMSS_Model(input_params_repro, Groups_repro, SaveTimeSteps = TRUE)

cat("Model with reproduction completed successfully!\n")
cat("Final fish abundances (smallest size classes):\n")
for (i in fish_groups) {
  min_size_idx <- which(round(log10(results_repro$model$param$w), digits = 2) == Groups_repro$W0[i])[1]
  final_abundance <- results_repro$abundances[i, min_size_idx]
  cat(sprintf("  %s: %.2e\n", Groups_repro$Species[i], final_abundance))
}

# Test 3: Compare results
cat("\nComparison of results:\n")
cat("Reproduction affects fish abundance dynamics: ")
if (exists("results_repro$model$reproduction_rate")) {
  cat("YES - Reproduction rates saved successfully\n")
} else {
  cat("NO - Check reproduction implementation\n")
}

cat("Fish pinning status:\n")
cat(sprintf("  Without reproduction: Pinning %s\n", 
            ifelse(results_no_repro$model$param$reproduction_on, "OFF", "ON")))
cat(sprintf("  With reproduction: Pinning %s\n", 
            ifelse(results_repro$model$param$reproduction_on, "OFF", "ON")))

cat("\nTest completed! Check results for biological realism.\n")
