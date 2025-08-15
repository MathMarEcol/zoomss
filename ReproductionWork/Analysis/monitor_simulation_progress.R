## Monitor Simulation Progress and Provide Interim Results
## This script checks for completion and provides preliminary analysis

# Function to check if simulation is complete
check_simulation_complete <- function() {
  return(file.exists("comprehensive_reproduction_comparison_results.RData"))
}

# Function to provide interim analysis if results are available
provide_interim_analysis <- function() {
  if (check_simulation_complete()) {
    cat("=== SIMULATION COMPLETED! ===\n\n")
    
    # Load results
    load("comprehensive_reproduction_comparison_results.RData")
    
    cat("PRELIMINARY RESULTS:\n")
    print(comparison_table)
    
    cat("\n=== KEY FINDINGS ===\n")
    
    # Quick analysis of major differences
    scenarios <- unique(comparison_table$Scenario)
    
    for (scenario in scenarios) {
      cat(sprintf("\n%s Scenario:\n", gsub("_", " ", toupper(scenario))))
      
      scenario_data <- comparison_table[comparison_table$Scenario == scenario, ]
      
      old_fish <- scenario_data$Fish_Biomass[scenario_data$Model_Version == "Old_ZooMSS"]
      new_no_repro <- scenario_data$Fish_Biomass[scenario_data$Model_Version == "New_ZooMSS_No_Repro"]
      new_with_repro <- scenario_data$Fish_Biomass[scenario_data$Model_Version == "New_ZooMSS_With_Repro"]
      
      cat(sprintf("  Fish Biomass:\n"))
      cat(sprintf("    Old ZooMSS: %.2e\n", old_fish))
      cat(sprintf("    New (no repro): %.2e (%.1f%% vs old)\n", 
                  new_no_repro, ((new_no_repro - old_fish) / old_fish) * 100))
      cat(sprintf("    New (with repro): %.2e (%.1f%% vs old)\n", 
                  new_with_repro, ((new_with_repro - old_fish) / old_fish) * 100))
      
      # Check for major differences
      if (abs((new_with_repro - old_fish) / old_fish) > 0.1) {
        cat("    *** SIGNIFICANT DIFFERENCE with reproduction! ***\n")
      }
      
      if (abs((new_no_repro - old_fish) / old_fish) > 0.05) {
        cat("    *** Framework itself affects results ***\n")
      }
    }
    
    # Overall patterns
    cat("\n=== OVERALL PATTERNS ===\n")
    
    # Chlorophyll response
    old_data <- comparison_table[comparison_table$Model_Version == "Old_ZooMSS", ]
    old_data <- old_data[order(old_data$Scenario), ]
    
    repro_data <- comparison_table[comparison_table$Model_Version == "New_ZooMSS_With_Repro", ]
    repro_data <- repro_data[order(repro_data$Scenario), ]
    
    cat("Chlorophyll Response (Fish Biomass):\n")
    cat(sprintf("  Old ZooMSS - Low to High ratio: %.2f\n", 
                old_data$Fish_Biomass[old_data$Scenario == "high"] / 
                old_data$Fish_Biomass[old_data$Scenario == "low"]))
    cat(sprintf("  New (with repro) - Low to High ratio: %.2f\n", 
                repro_data$Fish_Biomass[repro_data$Scenario == "high"] / 
                repro_data$Fish_Biomass[repro_data$Scenario == "low"]))
    
    cat("\n=== NEXT STEPS ===\n")
    cat("✓ Run analyze_reproduction_results.R for detailed analysis\n")
    cat("✓ Check generated plots and time series\n")
    cat("✓ Review model_version_summary.md for interpretation\n")
    
    return(TRUE)
  } else {
    cat("Simulation still running...\n")
    return(FALSE)
  }
}

# Main execution
cat("=== SIMULATION PROGRESS MONITOR ===\n")
cat("Checking for completed results...\n\n")

if (provide_interim_analysis()) {
  cat("\nSimulation analysis complete!\n")
} else {
  cat("Simulation in progress. Run this script again when complete.\n")
  cat("Expected completion time: ~30-60 minutes for 9 simulations of 120 years each.\n")
}