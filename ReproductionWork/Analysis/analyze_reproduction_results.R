## Analysis and Visualization of Reproduction Comparison Results
## This script analyzes the results from comprehensive_reproduction_comparison.R

# Load required libraries
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(reshape2)) install.packages("reshape2")
if (!require(gridExtra)) install.packages("gridExtra")

library(ggplot2)
library(reshape2)
library(gridExtra)

# Load results
if (file.exists("comprehensive_reproduction_comparison_results.RData")) {
  load("comprehensive_reproduction_comparison_results.RData")
  cat("Results loaded successfully!\n")
} else {
  stop("Results file not found. Please run comprehensive_reproduction_comparison.R first.")
}

# Function to extract time series data
extract_time_series <- function(result) {
  if (is.null(result$model$Biomass)) {
    # Calculate biomass from abundances if not available
    n_times <- dim(result$model$N)[1]
    n_groups <- dim(result$model$N)[2]
    weights <- result$model$param$w
    
    biomass_ts <- array(NA, dim = c(n_times, n_groups))
    for (t in 1:n_times) {
      for (g in 1:n_groups) {
        biomass_ts[t, g] <- sum(result$model$N[t, g, ] * weights)
      }
    }
  } else {
    biomass_ts <- result$model$Biomass
  }
  
  return(biomass_ts)
}

# Create time series plots
create_time_series_plots <- function() {
  cat("Creating time series plots...\n")
  
  # Extract time series for each simulation
  ts_data <- list()
  
  for (result_name in names(results)) {
    result <- results[[result_name]]
    biomass_ts <- extract_time_series(result)
    
    # Get fish and zooplankton groups
    fish_groups <- which(result$groups_used$Type == "Fish")
    zoo_groups <- which(result$groups_used$Type == "Zooplankton")
    
    # Calculate total fish and zoo biomass over time
    fish_biomass_ts <- rowSums(biomass_ts[, fish_groups, drop = FALSE])
    zoo_biomass_ts <- rowSums(biomass_ts[, zoo_groups, drop = FALSE])
    
    # Time vector (assuming isave = 100 and dt = 0.01)
    time_years <- seq(0, result$model$param$tmax, length.out = length(fish_biomass_ts))
    
    ts_data[[result_name]] <- data.frame(
      Time = time_years,
      Fish_Biomass = fish_biomass_ts,
      Zoo_Biomass = zoo_biomass_ts,
      Total_Biomass = fish_biomass_ts + zoo_biomass_ts,
      Model = result$version,
      Scenario = result$scenario
    )
  }
  
  # Combine all time series data
  all_ts_data <- do.call(rbind, ts_data)
  
  # Create plots for each scenario
  scenarios <- unique(all_ts_data$Scenario)
  
  plots <- list()
  
  for (scenario in scenarios) {
    scenario_data <- all_ts_data[all_ts_data$Scenario == scenario, ]
    
    # Fish biomass plot
    p1 <- ggplot(scenario_data, aes(x = Time, y = Fish_Biomass, color = Model)) +
      geom_line(size = 1) +
      labs(title = paste("Fish Biomass -", gsub("_", " ", scenario)),
           x = "Time (years)", y = "Fish Biomass") +
      theme_minimal() +
      scale_color_manual(values = c("Old_ZooMSS" = "blue", 
                                   "New_ZooMSS_No_Repro" = "green", 
                                   "New_ZooMSS_With_Repro" = "red"))
    
    # Total biomass plot
    p2 <- ggplot(scenario_data, aes(x = Time, y = Total_Biomass, color = Model)) +
      geom_line(size = 1) +
      labs(title = paste("Total Biomass -", gsub("_", " ", scenario)),
           x = "Time (years)", y = "Total Biomass") +
      theme_minimal() +
      scale_color_manual(values = c("Old_ZooMSS" = "blue", 
                                   "New_ZooMSS_No_Repro" = "green", 
                                   "New_ZooMSS_With_Repro" = "red"))
    
    plots[[paste(scenario, "fish", sep = "_")]] <- p1
    plots[[paste(scenario, "total", sep = "_")]] <- p2
  }
  
  return(plots)
}

# Create comparison bar plots
create_comparison_plots <- function() {
  cat("Creating comparison bar plots...\n")
  
  # Reshape comparison table for plotting
  comparison_melted <- melt(comparison_table, 
                           id.vars = c("Model_Version", "Scenario"),
                           measure.vars = c("Fish_Biomass", "Zoo_Biomass", "Total_Biomass"))
  
  # Fish biomass comparison
  fish_plot <- ggplot(comparison_table, aes(x = Scenario, y = Fish_Biomass, fill = Model_Version)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Final Fish Biomass Comparison",
         x = "Chlorophyll Scenario", y = "Fish Biomass") +
    theme_minimal() +
    scale_fill_manual(values = c("Old_ZooMSS" = "lightblue", 
                                "New_ZooMSS_No_Repro" = "lightgreen", 
                                "New_ZooMSS_With_Repro" = "lightcoral")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Total biomass comparison
  total_plot <- ggplot(comparison_table, aes(x = Scenario, y = Total_Biomass, fill = Model_Version)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Final Total Biomass Comparison",
         x = "Chlorophyll Scenario", y = "Total Biomass") +
    theme_minimal() +
    scale_fill_manual(values = c("Old_ZooMSS" = "lightblue", 
                                "New_ZooMSS_No_Repro" = "lightgreen", 
                                "New_ZooMSS_With_Repro" = "lightcoral")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Fish/Zoo ratio comparison
  ratio_plot <- ggplot(comparison_table, aes(x = Scenario, y = Fish_Zoo_Ratio, fill = Model_Version)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Fish/Zooplankton Biomass Ratio Comparison",
         x = "Chlorophyll Scenario", y = "Fish/Zoo Ratio") +
    theme_minimal() +
    scale_fill_manual(values = c("Old_ZooMSS" = "lightblue", 
                                "New_ZooMSS_No_Repro" = "lightgreen", 
                                "New_ZooMSS_With_Repro" = "lightcoral")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(fish = fish_plot, total = total_plot, ratio = ratio_plot))
}

# Analyze reproduction dynamics (if available)
analyze_reproduction_dynamics <- function() {
  cat("Analyzing reproduction dynamics...\n")
  
  repro_results <- results[grepl("With_Repro", names(results))]
  
  if (length(repro_results) == 0) {
    cat("No reproduction results found.\n")
    return(NULL)
  }
  
  repro_analysis <- list()
  
  for (result_name in names(repro_results)) {
    result <- repro_results[[result_name]]
    
    if (!is.null(result$model$reproduction_rate) && !is.null(result$model$recruitment)) {
      # Calculate average reproduction rates and recruitment over last 20% of simulation
      n_times <- dim(result$model$reproduction_rate)[1]
      last_20_percent <- floor(0.8 * n_times):n_times
      
      avg_repro_rates <- apply(result$model$reproduction_rate[last_20_percent, , ], c(2, 3), mean)
      avg_recruitment <- apply(result$model$recruitment[last_20_percent, ], 2, mean)
      
      fish_groups <- which(result$groups_used$Type == "Fish")
      
      repro_analysis[[result_name]] <- list(
        reproduction_rates = avg_repro_rates[fish_groups, ],
        recruitment = avg_recruitment[fish_groups],
        scenario = result$scenario
      )
    }
  }
  
  return(repro_analysis)
}

# Generate detailed report
generate_detailed_report <- function() {
  cat("\n=== DETAILED ANALYSIS REPORT ===\n\n")
  
  # Model version effects across scenarios
  cat("1. MODEL VERSION EFFECTS:\n")
  for (scenario in unique(comparison_table$Scenario)) {
    cat(sprintf("\n%s Scenario:\n", gsub("_", " ", scenario)))
    
    scenario_data <- comparison_table[comparison_table$Scenario == scenario, ]
    
    # Calculate percentage differences
    old_fish <- scenario_data$Fish_Biomass[scenario_data$Model_Version == "Old_ZooMSS"]
    new_no_repro_fish <- scenario_data$Fish_Biomass[scenario_data$Model_Version == "New_ZooMSS_No_Repro"]
    new_with_repro_fish <- scenario_data$Fish_Biomass[scenario_data$Model_Version == "New_ZooMSS_With_Repro"]
    
    cat(sprintf("  Fish Biomass Changes:\n"))
    cat(sprintf("    New (no repro) vs Old: %+.1f%%\n", 
                ((new_no_repro_fish - old_fish) / old_fish) * 100))
    cat(sprintf("    New (with repro) vs Old: %+.1f%%\n", 
                ((new_with_repro_fish - old_fish) / old_fish) * 100))
    cat(sprintf("    New (with repro) vs New (no repro): %+.1f%%\n", 
                ((new_with_repro_fish - new_no_repro_fish) / new_no_repro_fish) * 100))
  }
  
  # Chlorophyll response analysis
  cat("\n\n2. CHLOROPHYLL RESPONSE ANALYSIS:\n")
  for (model in unique(comparison_table$Model_Version)) {
    cat(sprintf("\n%s:\n", gsub("_", " ", model)))
    
    model_data <- comparison_table[comparison_table$Model_Version == model, ]
    model_data <- model_data[order(model_data$Scenario), ]
    
    low_fish <- model_data$Fish_Biomass[model_data$Scenario == "low"]
    med_fish <- model_data$Fish_Biomass[model_data$Scenario == "medium"]
    high_fish <- model_data$Fish_Biomass[model_data$Scenario == "high"]
    
    cat(sprintf("  Fish Biomass Response to Chlorophyll:\n"))
    cat(sprintf("    Low: %.2e\n", low_fish))
    cat(sprintf("    Medium: %.2e\n", med_fish))
    cat(sprintf("    High: %.2e\n", high_fish))
    cat(sprintf("    High/Low Ratio: %.2f\n", high_fish / low_fish))
  }
  
  # Reproduction analysis
  repro_analysis <- analyze_reproduction_dynamics()
  if (!is.null(repro_analysis)) {
    cat("\n\n3. REPRODUCTION DYNAMICS:\n")
    for (result_name in names(repro_analysis)) {
      analysis <- repro_analysis[[result_name]]
      cat(sprintf("\n%s:\n", gsub("_", " ", result_name)))
      cat(sprintf("  Average Recruitment Rates:\n"))
      for (i in 1:length(analysis$recruitment)) {
        cat(sprintf("    Fish Group %d: %.2e\n", i, analysis$recruitment[i]))
      }
    }
  }
}

# Main analysis execution
cat("=== REPRODUCTION COMPARISON ANALYSIS ===\n\n")

# Create visualizations
time_series_plots <- create_time_series_plots()
comparison_plots <- create_comparison_plots()

# Save plots
cat("Saving plots...\n")
ggsave("fish_biomass_comparison.png", comparison_plots$fish, width = 10, height = 6)
ggsave("total_biomass_comparison.png", comparison_plots$total, width = 10, height = 6)
ggsave("fish_zoo_ratio_comparison.png", comparison_plots$ratio, width = 10, height = 6)

# Save time series plots for each scenario
for (i in seq(1, length(time_series_plots), 2)) {
  scenario_name <- names(time_series_plots)[i]
  scenario_name <- gsub("_fish", "", scenario_name)
  
  combined_plot <- grid.arrange(time_series_plots[[i]], time_series_plots[[i+1]], ncol = 1)
  ggsave(paste0("time_series_", scenario_name, ".png"), combined_plot, width = 12, height = 8)
}

# Generate detailed report
generate_detailed_report()

# Save analysis results
save(comparison_table, time_series_plots, comparison_plots, 
     file = "reproduction_analysis_results.RData")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("✓ Time series plots created and saved\n")
cat("✓ Comparison plots created and saved\n")
cat("✓ Detailed analysis report generated\n")
cat("✓ Results saved to 'reproduction_analysis_results.RData'\n")
cat("\nCheck the generated PNG files for visualizations.\n")