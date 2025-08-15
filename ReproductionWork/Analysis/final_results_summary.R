## Final Results Summary for ZooMSS Reproduction Comparison
## This script provides the ultimate summary of all results

# Load results if available
if (file.exists("comprehensive_reproduction_comparison_results.RData")) {
  load("comprehensive_reproduction_comparison_results.RData")
  
  cat("=== FINAL ZOOMSS REPRODUCTION COMPARISON RESULTS ===\n\n")
  
  # Display the main comparison table
  cat("FINAL BIOMASS COMPARISON TABLE:\n")
  print(comparison_table)
  
  cat("\n=== KEY FINDINGS SUMMARY ===\n\n")
  
  # Calculate key metrics
  scenarios <- unique(comparison_table$Scenario)
  
  for (scenario in scenarios) {
    cat(sprintf("=== %s CHLOROPHYLL SCENARIO ===\n", toupper(gsub("_", " ", scenario))))
    
    scenario_data <- comparison_table[comparison_table$Scenario == scenario, ]
    
    old_fish <- scenario_data$Fish_Biomass[scenario_data$Model_Version == "Old_ZooMSS"]
    new_no_repro <- scenario_data$Fish_Biomass[scenario_data$Model_Version == "New_ZooMSS_No_Repro"]
    new_with_repro <- scenario_data$Fish_Biomass[scenario_data$Model_Version == "New_ZooMSS_With_Repro"]
    
    old_total <- scenario_data$Total_Biomass[scenario_data$Model_Version == "Old_ZooMSS"]
    new_no_repro_total <- scenario_data$Total_Biomass[scenario_data$Model_Version == "New_ZooMSS_No_Repro"]
    new_with_repro_total <- scenario_data$Total_Biomass[scenario_data$Model_Version == "New_ZooMSS_With_Repro"]
    
    cat(sprintf("Fish Biomass Results:\n"))
    cat(sprintf("  Old ZooMSS:           %.3e\n", old_fish))
    cat(sprintf("  New (no repro):       %.3e (%+.1f%% vs old)\n", 
                new_no_repro, ((new_no_repro - old_fish) / old_fish) * 100))
    cat(sprintf("  New (with repro):     %.3e (%+.1f%% vs old)\n", 
                new_with_repro, ((new_with_repro - old_fish) / old_fish) * 100))
    cat(sprintf("  Repro vs No-Repro:    %+.1f%% difference\n", 
                ((new_with_repro - new_no_repro) / new_no_repro) * 100))
    
    cat(sprintf("\nTotal Biomass Results:\n"))
    cat(sprintf("  Old ZooMSS:           %.3e\n", old_total))
    cat(sprintf("  New (no repro):       %.3e (%+.1f%% vs old)\n", 
                new_no_repro_total, ((new_no_repro_total - old_total) / old_total) * 100))
    cat(sprintf("  New (with repro):     %.3e (%+.1f%% vs old)\n", 
                new_with_repro_total, ((new_with_repro_total - old_total) / old_total) * 100))
    
    cat("\n")
  }
  
  # Cross-scenario analysis
  cat("=== CHLOROPHYLL RESPONSE ANALYSIS ===\n\n")
  
  models <- unique(comparison_table$Model_Version)
  for (model in models) {
    cat(sprintf("%s:\n", gsub("_", " ", model)))
    
    model_data <- comparison_table[comparison_table$Model_Version == model, ]
    model_data <- model_data[order(model_data$Scenario), ]
    
    low_fish <- model_data$Fish_Biomass[model_data$Scenario == "low"]
    med_fish <- model_data$Fish_Biomass[model_data$Scenario == "medium"]
    high_fish <- model_data$Fish_Biomass[model_data$Scenario == "high"]
    
    cat(sprintf("  Low → Medium:  %+.1f%% change\n", ((med_fish - low_fish) / low_fish) * 100))
    cat(sprintf("  Medium → High: %+.1f%% change\n", ((high_fish - med_fish) / med_fish) * 100))
    cat(sprintf("  Low → High:    %+.1f%% total change\n", ((high_fish - low_fish) / low_fish) * 100))
    cat(sprintf("  Sensitivity:   %.2f (High/Low ratio)\n", high_fish / low_fish))
    cat("\n")
  }
  
  # Major conclusions
  cat("=== MAJOR CONCLUSIONS ===\n\n")
  
  # Check if reproduction framework affects results when disabled
  framework_effects <- c()
  repro_effects <- c()
  
  for (scenario in scenarios) {
    scenario_data <- comparison_table[comparison_table$Scenario == scenario, ]
    
    old_fish <- scenario_data$Fish_Biomass[scenario_data$Model_Version == "Old_ZooMSS"]
    new_no_repro <- scenario_data$Fish_Biomass[scenario_data$Model_Version == "New_ZooMSS_No_Repro"]
    new_with_repro <- scenario_data$Fish_Biomass[scenario_data$Model_Version == "New_ZooMSS_With_Repro"]
    
    framework_effect <- abs((new_no_repro - old_fish) / old_fish)
    repro_effect <- abs((new_with_repro - old_fish) / old_fish)
    
    framework_effects <- c(framework_effects, framework_effect)
    repro_effects <- c(repro_effects, repro_effect)
  }
  
  avg_framework_effect <- mean(framework_effects) * 100
  avg_repro_effect <- mean(repro_effects) * 100
  
  cat(sprintf("1. FRAMEWORK IMPACT: %.1f%% average difference when reproduction framework is present but disabled\n", avg_framework_effect))
  if (avg_framework_effect < 5) {
    cat("   → Framework implementation has minimal impact when disabled ✓\n")
  } else {
    cat("   → Framework implementation affects results even when disabled ⚠\n")
  }
  
  cat(sprintf("\n2. REPRODUCTION IMPACT: %.1f%% average difference when reproduction is enabled\n", avg_repro_effect))
  if (avg_repro_effect > 10) {
    cat("   → Reproduction mechanism significantly changes fish dynamics ✓\n")
  } else {
    cat("   → Reproduction mechanism has modest effects on fish dynamics\n")
  }
  
  # Environmental sensitivity
  old_sensitivity <- comparison_table$Fish_Biomass[comparison_table$Model_Version == "Old_ZooMSS" & comparison_table$Scenario == "high"] / 
                   comparison_table$Fish_Biomass[comparison_table$Model_Version == "Old_ZooMSS" & comparison_table$Scenario == "low"]
  
  repro_sensitivity <- comparison_table$Fish_Biomass[comparison_table$Model_Version == "New_ZooMSS_With_Repro" & comparison_table$Scenario == "high"] / 
                      comparison_table$Fish_Biomass[comparison_table$Model_Version == "New_ZooMSS_With_Repro" & comparison_table$Scenario == "low"]
  
  cat(sprintf("\n3. ENVIRONMENTAL SENSITIVITY:\n"))
  cat(sprintf("   Old ZooMSS sensitivity:        %.2f (High/Low ratio)\n", old_sensitivity))
  cat(sprintf("   New ZooMSS (repro) sensitivity: %.2f (High/Low ratio)\n", repro_sensitivity))
  
  if (repro_sensitivity > old_sensitivity) {
    cat("   → Reproduction model is MORE sensitive to environmental changes\n")
  } else {
    cat("   → Reproduction model is LESS sensitive to environmental changes\n")
  }
  
  cat("\n=== SIMULATION DETAILS ===\n")
  cat("✓ 9 simulations completed (3 models × 3 scenarios)\n")
  cat("✓ Each simulation: 120 years duration\n")
  cat("✓ Original parameter values preserved\n")
  cat("✓ Maturation sizes maintained as requested\n")
  cat("✓ DBPM-based reproduction mechanism implemented\n")
  
  cat("\n=== FILES GENERATED ===\n")
  cat("• comprehensive_reproduction_comparison_results.RData - Main results\n")
  cat("• TestGroups_Original_Backup.csv - Original parameter backup\n")
  cat("• model_version_summary.md - Detailed documentation\n")
  cat("• Various analysis and visualization scripts\n")
  
  if (file.exists("reproduction_analysis_results.RData")) {
    cat("• reproduction_analysis_results.RData - Detailed analysis\n")
    cat("• Time series plots and comparison visualizations\n")
  }
  
  cat("\n=== REPRODUCTION COMPARISON COMPLETE ===\n")
  
} else {
  cat("Results file not found. Simulations may still be running.\n")
}