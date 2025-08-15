## Detailed diagnostic script for reproduction mechanism
## This script examines reproduction rates and energy allocation

# Clear workspace
rm(list = ls())

# Set working directory
setwd("c:/Users/kjmurphy/OneDrive - University of Tasmania/Documents/GitHub/ZooMSS")

# Load required functions
source("fZooMSS_Model.R")

# Load groups with reproduction enabled
Groups <- read.csv("TestGroups.csv")
Groups$Reproduction[Groups$Type == "Fish"] <- 1

# Set up parameters for longer run to see reproduction dynamics
input_params <- list(
  tmax = 50,           # 50 years - longer to see reproduction effects
  dt = 0.01,           # smaller time step
  phyto_int = -3,      
  phyto_slope = -1,    
  phyto_max = -8,      
  sst = 15             
)

cat("Running model with reproduction for", input_params$tmax, "years...\n")

# Run model with reproduction and save time steps
results <- fZooMSS_Model(input_params, Groups, SaveTimeSteps = TRUE)

cat("Model completed successfully!\n")

# Analyze results
fish_groups <- which(Groups$Type == "Fish")
param <- results$model$param

cat("\n=== REPRODUCTION DIAGNOSTICS ===\n")

# Check final abundances
cat("\nFinal fish abundances across all size classes:\n")
for (i in fish_groups) {
  total_abundance <- sum(results$abundances[i,])
  cat(sprintf("  %s: %.2e (total)\n", Groups$Species[i], total_abundance))
  
  # Check specific size classes
  min_size_idx <- which(round(log10(param$w), digits = 2) == Groups$W0[i])[1]
  max_size_idx <- which(round(log10(param$w), digits = 2) == Groups$Wmax[i])[1]
  
  cat(sprintf("    Smallest size class (%.1f): %.2e\n", Groups$W0[i], results$abundances[i, min_size_idx]))
  
  if (max_size_idx > min_size_idx) {
    mid_size_idx <- round((min_size_idx + max_size_idx) / 2)
    cat(sprintf("    Mid size class (%.1f): %.2e\n", log10(param$w[mid_size_idx]), results$abundances[i, mid_size_idx]))
  }
  
  cat(sprintf("    Largest size class (%.1f): %.2e\n", Groups$Wmax[i], results$abundances[i, max_size_idx]))
}

# Check if reproduction rates were saved
if (!is.null(results$model$reproduction_rate)) {
  cat("\nReproduction rates successfully saved.\n")
  
  # Check reproduction in mature size classes
  cat("\nReproduction rates in final time step:\n")
  final_time_idx <- dim(results$model$reproduction_rate)[1]
  
  for (i in fish_groups) {
    # Find mature size classes (> maturity size)
    maturity_size <- Groups$Wmat[i]
    mature_sizes <- which(log10(param$w) > maturity_size)
    
    if (length(mature_sizes) > 0) {
      avg_repro_rate <- mean(results$model$reproduction_rate[final_time_idx, i, mature_sizes])
      cat(sprintf("  %s (mature sizes): %.2e\n", Groups$Species[i], avg_repro_rate))
    }
  }
  
  # Check recruitment
  if (!is.null(results$model$recruitment)) {
    cat("\nRecruitment rates in final time step:\n")
    for (i in fish_groups) {
      recruitment_rate <- results$model$recruitment[final_time_idx, i]
      cat(sprintf("  %s: %.2e\n", Groups$Species[i], recruitment_rate))
    }
  }
} else {
  cat("\nWARNING: Reproduction rates not saved - check implementation.\n")
}

# Check energy allocation parameters
cat("\n=== ENERGY ALLOCATION ===\n")
cat("Growth efficiency (first fish group, first 5 size classes):\n")
first_fish <- fish_groups[1]
for (j in 1:min(5, length(param$w))) {
  growth_eff <- results$model$assim_eff[first_fish, j]
  repro_eff <- results$model$reproduction_eff[first_fish, j]
  cat(sprintf("  Size %.1f: Growth=%.3f, Reproduction=%.3f\n", 
              log10(param$w[j]), growth_eff, repro_eff))
}

# Check maturity sizes
cat("\nMaturity sizes for fish:\n")
for (i in fish_groups) {
  cat(sprintf("  %s: %.1f log10(g)\n", Groups$Species[i], Groups$Wmat[i]))
}

cat("\n=== RECOMMENDATIONS ===\n")
cat("If fish abundances are still very low:\n")
cat("1. Increase reproduction efficiency (currently 30% at maturity)\n")
cat("2. Reduce maturity size to allow earlier reproduction\n")
cat("3. Increase initial fish abundance\n")
cat("4. Run model for longer time to allow population buildup\n")

# Save detailed results for further analysis
save(results, Groups, input_params, file = "reproduction_diagnostics.RData")
cat("\nResults saved to 'reproduction_diagnostics.RData'\n")
