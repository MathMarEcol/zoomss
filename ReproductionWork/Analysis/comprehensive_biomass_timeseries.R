## Comprehensive Biomass Time Series Analysis for ZooMSS Reproduction Comparison
## Focus on biomass dynamics over time for community and functional groups

# Load required libraries
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(reshape2)) install.packages("reshape2")
if (!require(gridExtra)) install.packages("gridExtra")
if (!require(viridis)) install.packages("viridis")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")
if (!require(RColorBrewer)) install.packages("RColorBrewer")

library(ggplot2)
library(reshape2)
library(gridExtra)
library(viridis)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Create output directory
if (!dir.exists("ZooMSS_Biomass_TimeSeries")) {
  dir.create("ZooMSS_Biomass_TimeSeries")
}

# Load results
if (!file.exists("comprehensive_reproduction_comparison_results.RData")) {
  stop("Results file not found. Please run comprehensive_reproduction_comparison.R first.")
}

load("comprehensive_reproduction_comparison_results.RData")
cat("=== COMPREHENSIVE BIOMASS TIME SERIES ANALYSIS ===\n\n")

# Set up color schemes and theme
model_colors <- c("Old_ZooMSS" = "#1f77b4", 
                  "New_ZooMSS_No_Repro" = "#ff7f0e", 
                  "New_ZooMSS_With_Repro" = "#d62728")

# Enhanced color palette for functional groups
group_colors <- c(
  "Flagellates" = "#e41a1c", "Ciliates" = "#377eb8", "Larvaceans" = "#4daf4a",
  "OmniCopepods" = "#984ea3", "CarnCopepods" = "#ff7f00", "Euphausiids" = "#ffff33",
  "Chaetognaths" = "#a65628", "Salps" = "#f781bf", "Jellyfish" = "#999999",
  "Fish_Small" = "#1b9e77", "Fish_Med" = "#d95f02", "Fish_Large" = "#7570b3"
)

custom_theme <- theme_bw() + 
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# Function to extract biomass time series from model results
extract_biomass_timeseries <- function(result) {
  if (is.null(result$model$N)) {
    cat("Warning: No abundance time series found in result\n")
    return(NULL)
  }
  
  n_times <- dim(result$model$N)[1]
  n_groups <- dim(result$model$N)[2]
  n_sizes <- dim(result$model$N)[3]
  
  # Get weights and time vector
  weights <- result$model$param$w
  dx <- result$model$param$dx
  time_years <- seq(0, result$model$param$tmax, length.out = n_times)
  
  # Calculate biomass time series for each group
  biomass_ts <- array(NA, dim = c(n_times, n_groups))
  
  for (t in 1:n_times) {
    for (g in 1:n_groups) {
      # Calculate total biomass for this group at this time
      biomass_ts[t, g] <- sum(result$model$N[t, g, ] * weights * dx, na.rm = TRUE)
    }
  }
  
  # Create data frame
  ts_data <- data.frame()
  for (g in 1:n_groups) {
    temp_data <- data.frame(
      time = time_years,
      biomass = biomass_ts[, g],
      group = result$groups_used$Species[g],
      type = result$groups_used$Type[g],
      model = result$version,
      scenario = result$scenario
    )
    ts_data <- rbind(ts_data, temp_data)
  }
  
  return(ts_data)
}

# Extract time series data for all results
cat("Extracting biomass time series from all model runs...\n")
all_timeseries_data <- data.frame()

for (result_name in names(results)) {
  cat(sprintf("Processing %s...\n", result_name))
  result <- results[[result_name]]
  
  ts_data <- extract_biomass_timeseries(result)
  if (!is.null(ts_data)) {
    all_timeseries_data <- rbind(all_timeseries_data, ts_data)
  }
}

if (nrow(all_timeseries_data) == 0) {
  stop("No time series data could be extracted!")
}

cat(sprintf("Extracted time series data: %d rows\n", nrow(all_timeseries_data)))

# Function to create community-level time series plots
create_community_timeseries <- function(data) {
  plots <- list()
  
  for (scenario in unique(data$scenario)) {
    scenario_data <- data[data$scenario == scenario, ]
    
    # 1. Total community biomass
    community_biomass <- scenario_data %>%
      group_by(time, model) %>%
      summarise(total_biomass = sum(biomass, na.rm = TRUE), .groups = 'drop')
    
    p1 <- ggplot(community_biomass, aes(x = time, y = total_biomass, color = model)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = model_colors) +
      scale_y_log10(labels = scales::scientific) +
      labs(title = paste("Total Community Biomass:", gsub("_", " ", toupper(scenario))),
           x = "Time (years)", 
           y = "Total Biomass (g, log scale)",
           color = "Model Version") +
      custom_theme +
      theme(legend.position = "right")
    
    # 2. Fish vs Zooplankton biomass
    type_biomass <- scenario_data %>%
      group_by(time, model, type) %>%
      summarise(type_biomass = sum(biomass, na.rm = TRUE), .groups = 'drop')
    
    p2 <- ggplot(type_biomass, aes(x = time, y = type_biomass, color = model, linetype = type)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = model_colors) +
      scale_y_log10(labels = scales::scientific) +
      labs(title = paste("Fish vs Zooplankton Biomass:", gsub("_", " ", toupper(scenario))),
           x = "Time (years)", 
           y = "Biomass (g, log scale)",
           color = "Model Version",
           linetype = "Group Type") +
      custom_theme +
      theme(legend.position = "right")
    
    # 3. Fish:Zooplankton ratio over time
    fish_zoo_ratio <- type_biomass %>%
      pivot_wider(names_from = type, values_from = type_biomass) %>%
      mutate(fish_zoo_ratio = Fish / Zooplankton)
    
    p3 <- ggplot(fish_zoo_ratio, aes(x = time, y = fish_zoo_ratio, color = model)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = model_colors) +
      scale_y_log10() +
      labs(title = paste("Fish:Zooplankton Ratio:", gsub("_", " ", toupper(scenario))),
           x = "Time (years)", 
           y = "Fish:Zooplankton Biomass Ratio (log scale)",
           color = "Model Version") +
      custom_theme +
      theme(legend.position = "right")
    
    plots[[paste(scenario, "community", sep = "_")]] <- p1
    plots[[paste(scenario, "fish_zoo", sep = "_")]] <- p2
    plots[[paste(scenario, "ratio", sep = "_")]] <- p3
  }
  
  return(plots)
}

# Function to create functional group time series plots
create_functional_group_timeseries <- function(data) {
  plots <- list()
  
  for (scenario in unique(data$scenario)) {
    scenario_data <- data[data$scenario == scenario, ]
    
    # Fish groups time series
    fish_data <- scenario_data[scenario_data$type == "Fish", ]
    if (nrow(fish_data) > 0) {
      p_fish <- ggplot(fish_data, aes(x = time, y = biomass, color = group)) +
        geom_line(linewidth = 0.9) +
        facet_wrap(~model, scales = "free_y", ncol = 1) +
        scale_color_manual(values = group_colors[names(group_colors) %in% unique(fish_data$group)]) +
        scale_y_log10(labels = scales::scientific) +
        labs(title = paste("Fish Group Biomass Time Series:", gsub("_", " ", toupper(scenario))),
             x = "Time (years)", 
             y = "Biomass (g, log scale)",
             color = "Fish Group") +
        custom_theme +
        theme(legend.position = "bottom",
              strip.text = element_text(size = 10))
      
      plots[[paste(scenario, "fish_groups", sep = "_")]] <- p_fish
    }
    
    # Zooplankton groups time series
    zoo_data <- scenario_data[scenario_data$type == "Zooplankton", ]
    if (nrow(zoo_data) > 0) {
      p_zoo <- ggplot(zoo_data, aes(x = time, y = biomass, color = group)) +
        geom_line(linewidth = 0.9) +
        facet_wrap(~model, scales = "free_y", ncol = 1) +
        scale_color_manual(values = group_colors[names(group_colors) %in% unique(zoo_data$group)]) +
        scale_y_log10(labels = scales::scientific) +
        labs(title = paste("Zooplankton Group Biomass Time Series:", gsub("_", " ", toupper(scenario))),
             x = "Time (years)", 
             y = "Biomass (g, log scale)",
             color = "Zooplankton Group") +
        custom_theme +
        theme(legend.position = "bottom",
              strip.text = element_text(size = 10))
      
      plots[[paste(scenario, "zoo_groups", sep = "_")]] <- p_zoo
    }
    
    # All groups combined (smaller plot for overview)
    p_all <- ggplot(scenario_data, aes(x = time, y = biomass, color = group)) +
      geom_line(linewidth = 0.7, alpha = 0.8) +
      facet_grid(type ~ model, scales = "free_y") +
      scale_color_manual(values = group_colors[names(group_colors) %in% unique(scenario_data$group)]) +
      scale_y_log10(labels = scales::scientific) +
      labs(title = paste("All Functional Groups:", gsub("_", " ", toupper(scenario))),
           x = "Time (years)", 
           y = "Biomass (g, log scale)",
           color = "Functional Group") +
      custom_theme +
      theme(legend.position = "bottom",
            strip.text = element_text(size = 9),
            axis.text = element_text(size = 8))
    
    plots[[paste(scenario, "all_groups", sep = "_")]] <- p_all
  }
  
  return(plots)
}

# Function to create equilibrium analysis plots
create_equilibrium_analysis <- function(data) {
  plots <- list()
  
  # Calculate final 20 years average (equilibrium)
  max_time <- max(data$time)
  equilibrium_data <- data[data$time >= (max_time - 20), ]
  
  equilibrium_summary <- equilibrium_data %>%
    group_by(model, scenario, group, type) %>%
    summarise(mean_biomass = mean(biomass, na.rm = TRUE),
              sd_biomass = sd(biomass, na.rm = TRUE),
              cv_biomass = sd_biomass / mean_biomass,
              .groups = 'drop')
  
  # Community level equilibrium comparison
  community_equilibrium <- equilibrium_summary %>%
    group_by(model, scenario) %>%
    summarise(total_biomass = sum(mean_biomass, na.rm = TRUE), .groups = 'drop')
  
  p1 <- ggplot(community_equilibrium, aes(x = scenario, y = total_biomass, fill = model)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = model_colors) +
    scale_y_log10(labels = scales::scientific) +
    labs(title = "Equilibrium Community Biomass Comparison",
         x = "Chlorophyll Scenario", 
         y = "Total Biomass (g, log scale)",
         fill = "Model Version") +
    custom_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Fish vs Zooplankton equilibrium
  type_equilibrium <- equilibrium_summary %>%
    group_by(model, scenario, type) %>%
    summarise(type_biomass = sum(mean_biomass, na.rm = TRUE), .groups = 'drop')
  
  p2 <- ggplot(type_equilibrium, aes(x = scenario, y = type_biomass, fill = model)) +
    geom_col(position = "dodge", alpha = 0.8) +
    facet_wrap(~type, scales = "free_y") +
    scale_fill_manual(values = model_colors) +
    scale_y_log10(labels = scales::scientific) +
    labs(title = "Equilibrium Fish vs Zooplankton Biomass",
         x = "Chlorophyll Scenario", 
         y = "Biomass (g, log scale)",
         fill = "Model Version") +
    custom_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  plots[["equilibrium_community"]] <- p1
  plots[["equilibrium_fish_zoo"]] <- p2
  
  return(plots)
}

# Generate all plots
cat("Creating community-level time series plots...\n")
community_plots <- create_community_timeseries(all_timeseries_data)

cat("Creating functional group time series plots...\n")
functional_plots <- create_functional_group_timeseries(all_timeseries_data)

cat("Creating equilibrium analysis plots...\n")
equilibrium_plots <- create_equilibrium_analysis(all_timeseries_data)

# Combine all plots
all_plots <- c(community_plots, functional_plots, equilibrium_plots)

# Save all plots
cat("Saving biomass time series plots...\n")
for (plot_name in names(all_plots)) {
  filename <- paste0("ZooMSS_Biomass_TimeSeries/", plot_name, ".png")
  
  # Adjust plot size based on type
  if (grepl("groups", plot_name)) {
    ggsave(filename, all_plots[[plot_name]], width = 14, height = 10, dpi = 300)
  } else if (grepl("all_groups", plot_name)) {
    ggsave(filename, all_plots[[plot_name]], width = 16, height = 12, dpi = 300)
  } else {
    ggsave(filename, all_plots[[plot_name]], width = 12, height = 8, dpi = 300)
  }
  
  cat(sprintf("Saved: %s\n", filename))
}

# Save plot objects and data
save(all_plots, all_timeseries_data, 
     file = "ZooMSS_Biomass_TimeSeries/biomass_timeseries_analysis.RData")

# Create summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")

# Final biomass comparison
final_time <- max(all_timeseries_data$time)
final_data <- all_timeseries_data[all_timeseries_data$time == final_time, ]

summary_stats <- final_data %>%
  group_by(model, scenario, type) %>%
  summarise(total_biomass = sum(biomass, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = type, values_from = total_biomass) %>%
  mutate(fish_zoo_ratio = Fish / Zooplankton,
         total_community = Fish + Zooplankton)

print("Final Biomass Summary (after 120 years):")
print(summary_stats)

# Calculate percentage differences
cat("\n=== REPRODUCTION MECHANISM IMPACT ===\n")
for (scenario in unique(summary_stats$scenario)) {
  cat(sprintf("\n%s Scenario:\n", gsub("_", " ", toupper(scenario))))
  
  scenario_stats <- summary_stats[summary_stats$scenario == scenario, ]
  
  old_fish <- scenario_stats$Fish[scenario_stats$model == "Old_ZooMSS"]
  new_no_repro_fish <- scenario_stats$Fish[scenario_stats$model == "New_ZooMSS_No_Repro"]
  new_with_repro_fish <- scenario_stats$Fish[scenario_stats$model == "New_ZooMSS_With_Repro"]
  
  cat("  Fish Biomass Changes:\n")
  cat(sprintf("    Reproduction ON vs OFF: %+.1f%%\n", 
              100 * (new_with_repro_fish - new_no_repro_fish) / new_no_repro_fish))
  cat(sprintf("    Reproduction ON vs Old: %+.1f%%\n", 
              100 * (new_with_repro_fish - old_fish) / old_fish))
  
  old_total <- scenario_stats$total_community[scenario_stats$model == "Old_ZooMSS"]
  new_with_repro_total <- scenario_stats$total_community[scenario_stats$model == "New_ZooMSS_With_Repro"]
  
  cat(sprintf("  Total Community Biomass Change (Reproduction vs Old): %+.1f%%\n", 
              100 * (new_with_repro_total - old_total) / old_total))
}

cat("\n=== BIOMASS TIME SERIES ANALYSIS COMPLETE ===\n")
cat("Generated comprehensive biomass time series plots:\n")
cat("✓ Total community biomass time series for each scenario\n")
cat("✓ Fish vs zooplankton biomass time series\n")
cat("✓ Fish:zooplankton ratio time series\n")
cat("✓ Individual fish group time series by model version\n")
cat("✓ Individual zooplankton group time series by model version\n")
cat("✓ All functional groups overview plots\n")
cat("✓ Equilibrium biomass comparison plots\n")
cat("✓ Summary statistics and reproduction mechanism impact analysis\n")
cat(sprintf("✓ All plots saved in 'ZooMSS_Biomass_TimeSeries/' directory (%d plots total)\n", length(all_plots)))