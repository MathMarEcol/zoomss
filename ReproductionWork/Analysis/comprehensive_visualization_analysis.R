## Comprehensive Visualization Analysis for ZooMSS Reproduction Comparison
## Creates detailed figures comparing size spectra, abundance, growth, reproduction, and diet

# Load required libraries
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(reshape2)) install.packages("reshape2")
if (!require(gridExtra)) install.packages("gridExtra")
if (!require(viridis)) install.packages("viridis")
if (!require(RColorBrewer)) install.packages("RColorBrewer")

library(ggplot2)
library(reshape2)
library(gridExtra)
library(viridis)
library(RColorBrewer)

# Create figures directory
if (!dir.exists("ZooMSS_Reproduction_Figures")) {
  dir.create("ZooMSS_Reproduction_Figures")
}

# Load results
if (!file.exists("comprehensive_reproduction_comparison_results.RData")) {
  stop("Results file not found. Please run comprehensive_reproduction_comparison.R first.")
}

load("comprehensive_reproduction_comparison_results.RData")
cat("=== COMPREHENSIVE VISUALIZATION ANALYSIS ===\n\n")

# Set up color schemes
model_colors <- c("Old_ZooMSS" = "#2E86AB", 
                  "New_ZooMSS_No_Repro" = "#A23B72", 
                  "New_ZooMSS_With_Repro" = "#F18F01")

scenario_colors <- c("low" = "#440154", "medium" = "#21908C", "high" = "#FDE725")

# Function to create size spectrum plots
create_size_spectrum_plots <- function() {
  cat("Creating size spectrum plots...\n")
  
  plots <- list()
  
  for (scenario_name in names(chlorophyll_scenarios)) {
    scenario_data <- list()
    
    for (result_name in names(results)) {
      if (grepl(scenario_name, result_name)) {
        result <- results[[result_name]]
        
        # Get final abundances and weights
        final_abundances <- result$abundances
        weights <- result$model$param$w
        log_weights <- log10(weights)
        
        # Create data frame for plotting
        for (i in 1:nrow(final_abundances)) {
          group_name <- result$groups_used$Species[i]
          group_type <- result$groups_used$Type[i]
          
          temp_data <- data.frame(
            log_weight = log_weights,
            abundance = final_abundances[i, ],
            biomass = final_abundances[i, ] * weights,
            group = group_name,
            type = group_type,
            model = result$version,
            scenario = scenario_name
          )
          
          scenario_data[[paste(result$version, group_name, sep = "_")]] <- temp_data
        }
      }
    }
    
    # Combine all data for this scenario
    combined_data <- do.call(rbind, scenario_data)
    
    # Size spectrum plot (abundance)
    p1 <- ggplot(combined_data, aes(x = log_weight, y = log10(abundance + 1e-10), 
                                   color = model, linetype = type)) +
      geom_line(size = 1) +
      facet_wrap(~group, scales = "free_y") +
      scale_color_manual(values = model_colors) +
      labs(title = paste("Size Spectrum - Abundance:", gsub("_", " ", toupper(scenario_name))),
           x = "Log10 Body Weight (g)", y = "Log10 Abundance",
           color = "Model Version", linetype = "Group Type") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Size spectrum plot (biomass)
    p2 <- ggplot(combined_data, aes(x = log_weight, y = log10(biomass + 1e-10), 
                                   color = model, linetype = type)) +
      geom_line(size = 1) +
      facet_wrap(~group, scales = "free_y") +
      scale_color_manual(values = model_colors) +
      labs(title = paste("Size Spectrum - Biomass:", gsub("_", " ", toupper(scenario_name))),
           x = "Log10 Body Weight (g)", y = "Log10 Biomass",
           color = "Model Version", linetype = "Group Type") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    plots[[paste(scenario_name, "abundance", sep = "_")]] <- p1
    plots[[paste(scenario_name, "biomass", sep = "_")]] <- p2
  }
  
  return(plots)
}

# Function to create growth rate comparison plots
create_growth_rate_plots <- function() {
  cat("Creating growth rate plots...\n")
  
  plots <- list()
  
  for (scenario_name in names(chlorophyll_scenarios)) {
    growth_data <- list()
    
    for (result_name in names(results)) {
      if (grepl(scenario_name, result_name)) {
        result <- results[[result_name]]
        
        # Get growth rates
        growth_rates <- result$growth
        weights <- result$model$param$w
        log_weights <- log10(weights)
        
        # Create data frame for plotting
        for (i in 1:nrow(growth_rates)) {
          group_name <- result$groups_used$Species[i]
          group_type <- result$groups_used$Type[i]
          
          temp_data <- data.frame(
            log_weight = log_weights,
            growth_rate = growth_rates[i, ],
            group = group_name,
            type = group_type,
            model = result$version,
            scenario = scenario_name
          )
          
          growth_data[[paste(result$version, group_name, sep = "_")]] <- temp_data
        }
      }
    }
    
    # Combine all data for this scenario
    combined_data <- do.call(rbind, growth_data)
    
    # Growth rate plot
    p <- ggplot(combined_data, aes(x = log_weight, y = growth_rate, 
                                  color = model, linetype = type)) +
      geom_line(size = 1) +
      facet_wrap(~group, scales = "free_y") +
      scale_color_manual(values = model_colors) +
      labs(title = paste("Growth Rates:", gsub("_", " ", toupper(scenario_name))),
           x = "Log10 Body Weight (g)", y = "Growth Rate",
           color = "Model Version", linetype = "Group Type") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    plots[[scenario_name]] <- p
  }
  
  return(plots)
}

# Function to create reproduction allocation plots
create_reproduction_plots <- function() {
  cat("Creating reproduction allocation plots...\n")
  
  plots <- list()
  
  # Only analyze reproduction-enabled results
  repro_results <- results[grepl("With_Repro", names(results))]
  
  if (length(repro_results) == 0) {
    cat("No reproduction results found for detailed analysis.\n")
    return(plots)
  }
  
  for (result_name in names(repro_results)) {
    result <- repro_results[[result_name]]
    scenario <- result$scenario
    
    if (!is.null(result$model$reproduction_rate)) {
      # Get reproduction rates and weights
      repro_rates <- apply(result$model$reproduction_rate, c(2, 3), mean, na.rm = TRUE)
      weights <- result$model$param$w
      log_weights <- log10(weights)
      
      # Create data frame for fish groups only
      fish_groups <- which(result$groups_used$Type == "Fish")
      repro_data <- list()
      
      for (i in fish_groups) {
        group_name <- result$groups_used$Species[i]
        
        temp_data <- data.frame(
          log_weight = log_weights,
          reproduction_rate = repro_rates[i, ],
          group = group_name,
          scenario = scenario
        )
        
        repro_data[[group_name]] <- temp_data
      }
      
      # Combine data
      combined_data <- do.call(rbind, repro_data)
      
      # Reproduction allocation plot
      p <- ggplot(combined_data, aes(x = log_weight, y = reproduction_rate, color = group)) +
        geom_line(size = 1.2) +
        scale_color_brewer(type = "qual", palette = "Set1") +
        labs(title = paste("Reproduction Allocation:", gsub("_", " ", toupper(scenario))),
             x = "Log10 Body Weight (g)", y = "Reproduction Rate",
             color = "Fish Group") +
        theme_minimal() +
        theme(legend.position = "bottom")
      
      plots[[scenario]] <- p
    }
  }
  
  return(plots)
}

# Function to create diet comparison plots
create_diet_plots <- function() {
  cat("Creating diet comparison plots...\n")
  
  plots <- list()
  
  for (scenario_name in names(chlorophyll_scenarios)) {
    diet_data <- list()
    
    for (result_name in names(results)) {
      if (grepl(scenario_name, result_name)) {
        result <- results[[result_name]]
        
        # Get diet data
        diets <- result$diets
        
        # Create data frame for plotting
        for (i in 1:dim(diets)[1]) {
          group_name <- result$groups_used$Species[i]
          group_type <- result$groups_used$Type[i]
          
          # Diet composition (assuming columns are: pico, nano, micro, then functional groups)
          diet_composition <- diets[i, ]
          
          temp_data <- data.frame(
            prey_type = c("Pico_Phyto", "Nano_Phyto", "Micro_Phyto", 
                         paste("Group", 4:length(diet_composition))),
            diet_fraction = diet_composition,
            predator = group_name,
            predator_type = group_type,
            model = result$version,
            scenario = scenario_name
          )
          
          diet_data[[paste(result$version, group_name, sep = "_")]] <- temp_data
        }
      }
    }
    
    # Combine all data for this scenario
    combined_data <- do.call(rbind, diet_data)
    
    # Diet composition plot
    p <- ggplot(combined_data, aes(x = predator, y = diet_fraction, 
                                  fill = prey_type)) +
      geom_bar(stat = "identity", position = "stack") +
      facet_wrap(~model, scales = "free_x") +
      scale_fill_viridis_d() +
      labs(title = paste("Diet Composition:", gsub("_", " ", toupper(scenario_name))),
           x = "Predator Group", y = "Diet Fraction",
           fill = "Prey Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom")
    
    plots[[scenario_name]] <- p
  }
  
  return(plots)
}

# Function to create time series plots
create_time_series_plots <- function() {
  cat("Creating time series plots...\n")
  
  plots <- list()
  
  for (scenario_name in names(chlorophyll_scenarios)) {
    ts_data <- list()
    
    for (result_name in names(results)) {
      if (grepl(scenario_name, result_name)) {
        result <- results[[result_name]]
        
        if (!is.null(result$model$N)) {
          # Calculate biomass time series
          n_times <- dim(result$model$N)[1]
          n_groups <- dim(result$model$N)[2]
          weights <- result$model$param$w
          
          # Calculate total biomass by group over time
          biomass_ts <- array(NA, dim = c(n_times, n_groups))
          for (t in 1:n_times) {
            for (g in 1:n_groups) {
              biomass_ts[t, g] <- sum(result$model$N[t, g, ] * weights)
            }
          }
          
          # Time vector
          time_years <- seq(0, result$model$param$tmax, length.out = n_times)
          
          # Separate fish and zooplankton
          fish_groups <- which(result$groups_used$Type == "Fish")
          zoo_groups <- which(result$groups_used$Type == "Zooplankton")
          
          fish_biomass <- rowSums(biomass_ts[, fish_groups, drop = FALSE])
          zoo_biomass <- rowSums(biomass_ts[, zoo_groups, drop = FALSE])
          
          temp_data <- data.frame(
            time = time_years,
            fish_biomass = fish_biomass,
            zoo_biomass = zoo_biomass,
            total_biomass = fish_biomass + zoo_biomass,
            model = result$version,
            scenario = scenario_name
          )
          
          ts_data[[result$version]] <- temp_data
        }
      }
    }
    
    # Combine data
    combined_data <- do.call(rbind, ts_data)
    
    # Time series plot
    p1 <- ggplot(combined_data, aes(x = time, y = fish_biomass, color = model)) +
      geom_line(size = 1) +
      scale_color_manual(values = model_colors) +
      labs(title = paste("Fish Biomass Time Series:", gsub("_", " ", toupper(scenario_name))),
           x = "Time (years)", y = "Fish Biomass",
           color = "Model Version") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    p2 <- ggplot(combined_data, aes(x = time, y = total_biomass, color = model)) +
      geom_line(size = 1) +
      scale_color_manual(values = model_colors) +
      labs(title = paste("Total Biomass Time Series:", gsub("_", " ", toupper(scenario_name))),
           x = "Time (years)", y = "Total Biomass",
           color = "Model Version") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    plots[[paste(scenario_name, "fish", sep = "_")]] <- p1
    plots[[paste(scenario_name, "total", sep = "_")]] <- p2
  }
  
  return(plots)
}

# Function to create summary comparison plots
create_summary_plots <- function() {
  cat("Creating summary comparison plots...\n")
  
  # Final biomass comparison
  p1 <- ggplot(comparison_table, aes(x = Scenario, y = Fish_Biomass, fill = Model_Version)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = model_colors) +
    scale_y_log10() +
    labs(title = "Final Fish Biomass Comparison (Log Scale)",
         x = "Chlorophyll Scenario", y = "Fish Biomass (log10)",
         fill = "Model Version") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  p2 <- ggplot(comparison_table, aes(x = Scenario, y = Total_Biomass, fill = Model_Version)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = model_colors) +
    scale_y_log10() +
    labs(title = "Final Total Biomass Comparison (Log Scale)",
         x = "Chlorophyll Scenario", y = "Total Biomass (log10)",
         fill = "Model Version") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  p3 <- ggplot(comparison_table, aes(x = Scenario, y = Fish_Zoo_Ratio, fill = Model_Version)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = model_colors) +
    labs(title = "Fish/Zooplankton Biomass Ratio Comparison",
         x = "Chlorophyll Scenario", y = "Fish/Zoo Ratio",
         fill = "Model Version") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(list(fish_biomass = p1, total_biomass = p2, fish_zoo_ratio = p3))
}

# Generate all plots
cat("Generating comprehensive visualizations...\n\n")

# Create all plot sets
size_spectrum_plots <- create_size_spectrum_plots()
growth_rate_plots <- create_growth_rate_plots()
reproduction_plots <- create_reproduction_plots()
diet_plots <- create_diet_plots()
time_series_plots <- create_time_series_plots()
summary_plots <- create_summary_plots()

# Save all plots
cat("Saving plots to ZooMSS_Reproduction_Figures/...\n")

# Size spectrum plots
for (plot_name in names(size_spectrum_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures/size_spectrum_", plot_name, ".png")
  ggsave(filename, size_spectrum_plots[[plot_name]], width = 12, height = 8, dpi = 300)
}

# Growth rate plots
for (plot_name in names(growth_rate_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures/growth_rates_", plot_name, ".png")
  ggsave(filename, growth_rate_plots[[plot_name]], width = 12, height = 8, dpi = 300)
}

# Reproduction plots
for (plot_name in names(reproduction_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures/reproduction_allocation_", plot_name, ".png")
  ggsave(filename, reproduction_plots[[plot_name]], width = 10, height = 6, dpi = 300)
}

# Diet plots
for (plot_name in names(diet_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures/diet_composition_", plot_name, ".png")
  ggsave(filename, diet_plots[[plot_name]], width = 12, height = 8, dpi = 300)
}

# Time series plots
for (plot_name in names(time_series_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures/time_series_", plot_name, ".png")
  ggsave(filename, time_series_plots[[plot_name]], width = 12, height = 6, dpi = 300)
}

# Summary plots
for (plot_name in names(summary_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures/summary_", plot_name, ".png")
  ggsave(filename, summary_plots[[plot_name]], width = 10, height = 6, dpi = 300)
}

# Create a combined summary figure
combined_summary <- grid.arrange(
  summary_plots$fish_biomass,
  summary_plots$total_biomass,
  summary_plots$fish_zoo_ratio,
  ncol = 1
)

ggsave("ZooMSS_Reproduction_Figures/combined_summary.png", combined_summary, 
       width = 12, height = 15, dpi = 300)

# Save plot objects for further analysis
save(size_spectrum_plots, growth_rate_plots, reproduction_plots, diet_plots, 
     time_series_plots, summary_plots, 
     file = "ZooMSS_Reproduction_Figures/all_plots.RData")

cat("\n=== VISUALIZATION ANALYSIS COMPLETE ===\n")
cat("Generated figures:\n")
cat("✓ Size spectrum plots (abundance and biomass) for each scenario\n")
cat("✓ Growth rate comparisons across model versions\n")
cat("✓ Reproduction allocation plots (for reproduction-enabled model)\n")
cat("✓ Diet composition comparisons\n")
cat("✓ Time series plots showing 120-year dynamics\n")
cat("✓ Summary comparison plots\n")
cat("✓ Combined summary figure\n")
cat("\nAll figures saved in 'ZooMSS_Reproduction_Figures/' directory\n")
cat("Plot objects saved in 'all_plots.RData' for further analysis\n")