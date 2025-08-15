## Diagnostic and Time Series Analysis for ZooMSS Reproduction Comparison
## Investigates reproduction plot issues and creates comprehensive time series plots

# Load required libraries
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(reshape2)) install.packages("reshape2")
if (!require(gridExtra)) install.packages("gridExtra")
if (!require(viridis)) install.packages("viridis")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")

library(ggplot2)
library(reshape2)
library(gridExtra)
library(viridis)
library(dplyr)
library(tidyr)

# Create diagnostic figures directory
if (!dir.exists("ZooMSS_Diagnostic_TimeSeries")) {
  dir.create("ZooMSS_Diagnostic_TimeSeries")
}

# Load results
if (!file.exists("comprehensive_reproduction_comparison_results.RData")) {
  stop("Results file not found. Please run comprehensive_reproduction_comparison.R first.")
}

load("comprehensive_reproduction_comparison_results.RData")
cat("=== DIAGNOSTIC AND TIME SERIES ANALYSIS ===\n\n")

# Set up color schemes and theme
model_colors <- c("Old_ZooMSS" = "#1f77b4", 
                  "New_ZooMSS_No_Repro" = "#ff7f0e", 
                  "New_ZooMSS_With_Repro" = "#d62728")

group_colors <- RColorBrewer::brewer.pal(12, "Set3")
names(group_colors) <- c("Flagellates", "Ciliates", "Larvaceans", "OmniCopepods", 
                        "CarnCopepods", "Euphausiids", "Chaetognaths", "Salps", 
                        "Jellyfish", "Fish_Small", "Fish_Med", "Fish_Large")

custom_theme <- theme_bw() + 
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10)
  )

# Function to diagnose reproduction data
diagnose_reproduction_data <- function() {
  cat("Diagnosing reproduction data...\n")
  
  repro_results <- results[grepl("With_Repro", names(results))]
  
  for (result_name in names(repro_results)) {
    result <- repro_results[[result_name]]
    cat(sprintf("\n=== %s ===\n", result_name))
    
    if (!is.null(result$model$reproduction_rate)) {
      repro_rates <- result$model$reproduction_rate
      cat(sprintf("Reproduction rate array dimensions: %s\n", paste(dim(repro_rates), collapse = " x ")))
      cat(sprintf("Range of reproduction rates: %.6f to %.6f\n", 
                  min(repro_rates, na.rm = TRUE), max(repro_rates, na.rm = TRUE)))
      
      # Check if there are non-zero values
      non_zero_count <- sum(repro_rates > 1e-10, na.rm = TRUE)
      cat(sprintf("Number of non-zero reproduction rates: %d\n", non_zero_count))
      
      # Check fish groups specifically
      fish_groups <- which(result$groups_used$Type == "Fish")
      for (i in fish_groups) {
        group_name <- result$groups_used$Species[i]
        group_repro <- repro_rates[, i, ]
        cat(sprintf("  %s: range %.6f to %.6f, non-zero: %d\n", 
                    group_name, 
                    min(group_repro, na.rm = TRUE), 
                    max(group_repro, na.rm = TRUE),
                    sum(group_repro > 1e-10, na.rm = TRUE)))
      }
    } else {
      cat("No reproduction rate data found!\n")
    }
    
    if (!is.null(result$model$recruitment)) {
      recruitment <- result$model$recruitment
      cat(sprintf("Recruitment array dimensions: %s\n", paste(dim(recruitment), collapse = " x ")))
      cat(sprintf("Range of recruitment: %.6f to %.6f\n", 
                  min(recruitment, na.rm = TRUE), max(recruitment, na.rm = TRUE)))
    } else {
      cat("No recruitment data found!\n")
    }
  }
}

# Function to create proper reproduction plots using time series data
create_proper_reproduction_plots <- function() {
  cat("Creating proper reproduction plots using time series data...\n")
  
  plots <- list()
  repro_results <- results[grepl("With_Repro", names(results))]
  
  if (length(repro_results) == 0) {
    cat("No reproduction results found.\n")
    return(plots)
  }
  
  for (result_name in names(repro_results)) {
    result <- repro_results[[result_name]]
    scenario <- result$scenario
    
    if (!is.null(result$model$reproduction_rate)) {
      repro_rates <- result$model$reproduction_rate
      weights <- result$model$param$w
      log_weights <- log10(weights)
      
      # Get time vector
      n_times <- dim(repro_rates)[1]
      time_years <- seq(0, result$model$param$tmax, length.out = n_times)
      
      # Fish groups only
      fish_groups <- which(result$groups_used$Type == "Fish")
      
      # Create time series plot of reproduction rates
      repro_ts_data <- data.frame()
      
      for (i in fish_groups) {
        group_name <- result$groups_used$Species[i]
        w0_log <- result$groups_used$W0[i]
        wmax_log <- result$groups_used$Wmax[i]
        wmat_log <- result$groups_used$Wmat[i]
        
        # Find size classes within group range
        size_mask <- log_weights >= w0_log & log_weights <= wmax_log
        
        if (any(size_mask)) {
          # Calculate total reproduction rate for this group over time
          group_repro_ts <- rowSums(repro_rates[, i, size_mask, drop = FALSE], na.rm = TRUE)
          
          temp_data <- data.frame(
            time = time_years,
            reproduction_rate = group_repro_ts,
            group = group_name,
            scenario = scenario,
            maturity_size = wmat_log
          )
          
          repro_ts_data <- rbind(repro_ts_data, temp_data)
        }
      }
      
      if (nrow(repro_ts_data) > 0) {
        # Time series of reproduction rates
        p1 <- ggplot(repro_ts_data, aes(x = time, y = reproduction_rate, color = group)) +
          geom_line(linewidth = 1) +
          scale_color_brewer(type = "qual", palette = "Set1") +
          labs(title = paste("Reproduction Rate Time Series:", gsub("_", " ", toupper(scenario))),
               x = "Time (years)", 
               y = "Total Reproduction Rate",
               color = "Fish Group") +
          custom_theme
        
        # Size-specific reproduction rates (final time step)
        final_repro <- repro_rates[n_times, , ]
        size_repro_data <- data.frame()
        
        for (i in fish_groups) {
          group_name <- result$groups_used$Species[i]
          w0_log <- result$groups_used$W0[i]
          wmax_log <- result$groups_used$Wmax[i]
          wmat_log <- result$groups_used$Wmat[i]
          
          size_mask <- log_weights >= w0_log & log_weights <= wmax_log
          
          if (any(size_mask)) {
            temp_data <- data.frame(
              log_weight = log_weights[size_mask],
              reproduction_rate = final_repro[i, size_mask],
              group = group_name,
              scenario = scenario,
              maturity_size = wmat_log
            )
            
            size_repro_data <- rbind(size_repro_data, temp_data)
          }
        }
        
        if (nrow(size_repro_data) > 0) {
          p2 <- ggplot(size_repro_data, aes(x = log_weight, y = reproduction_rate, color = group)) +
            geom_line(linewidth = 1.2) +
            geom_vline(aes(xintercept = maturity_size, color = group), 
                       linetype = "dashed", alpha = 0.7) +
            scale_color_brewer(type = "qual", palette = "Set1") +
            labs(title = paste("Final Reproduction Allocation:", gsub("_", " ", toupper(scenario))),
                 x = "Log10 Body Weight (g)", 
                 y = "Reproduction Rate",
                 color = "Fish Group",
                 caption = "Dashed lines show maturity sizes") +
            custom_theme
          
          plots[[paste(scenario, "timeseries", sep = "_")]] <- p1
          plots[[paste(scenario, "allocation", sep = "_")]] <- p2
        }
      }
    }
  }
  
  return(plots)
}

# Function to create comprehensive time series plots
create_comprehensive_time_series <- function() {
  cat("Creating comprehensive time series plots...\n")
  
  plots <- list()
  
  for (scenario_name in names(chlorophyll_scenarios)) {
    cat(sprintf("Processing %s scenario...\n", scenario_name))
    
    # Collect time series data for this scenario
    ts_data_list <- list()
    
    for (result_name in names(results)) {
      if (grepl(scenario_name, result_name)) {
        result <- results[[result_name]]
        
        if (!is.null(result$model$N)) {
          n_times <- dim(result$model$N)[1]
          n_groups <- dim(result$model$N)[2]
          weights <- result$model$param$w
          
          # Time vector
          time_years <- seq(0, result$model$param$tmax, length.out = n_times)
          
          # Calculate biomass time series for each group
          group_biomass_ts <- array(NA, dim = c(n_times, n_groups))
          
          for (t in 1:n_times) {
            for (g in 1:n_groups) {
              group_biomass_ts[t, g] <- sum(result$model$N[t, g, ] * weights)
            }
          }
          
          # Create data frame
          ts_data <- data.frame(
            time = rep(time_years, n_groups),
            biomass = as.vector(group_biomass_ts),
            group = rep(result$groups_used$Species, each = n_times),
            type = rep(result$groups_used$Type, each = n_times),
            model = result$version,
            scenario = scenario_name
          )
          
          ts_data_list[[result$version]] <- ts_data
        }
      }
    }
    
    if (length(ts_data_list) > 0) {
      # Combine all data for this scenario
      combined_ts_data <- do.call(rbind, ts_data_list)
      
      # 1. Total community biomass time series
      community_ts <- combined_ts_data %>%
        group_by(time, model) %>%
        summarise(total_biomass = sum(biomass), .groups = 'drop')
      
      p1 <- ggplot(community_ts, aes(x = time, y = total_biomass, color = model)) +
        geom_line(linewidth = 1.2) +
        scale_color_manual(values = model_colors) +
        scale_y_log10() +
        labs(title = paste("Total Community Biomass:", gsub("_", " ", toupper(scenario_name))),
             x = "Time (years)", 
             y = "Total Biomass (g, log scale)",
             color = "Model Version") +
        custom_theme
      
      # 2. Fish vs Zooplankton biomass time series
      type_ts <- combined_ts_data %>%
        group_by(time, model, type) %>%
        summarise(type_biomass = sum(biomass), .groups = 'drop')
      
      p2 <- ggplot(type_ts, aes(x = time, y = type_biomass, color = model, linetype = type)) +
        geom_line(linewidth = 1) +
        scale_color_manual(values = model_colors) +
        scale_y_log10() +
        labs(title = paste("Fish vs Zooplankton Biomass:", gsub("_", " ", toupper(scenario_name))),
             x = "Time (years)", 
             y = "Biomass (g, log scale)",
             color = "Model Version",
             linetype = "Group Type") +
        custom_theme
      
      # 3. Individual functional group time series (fish only for clarity)
      fish_ts <- combined_ts_data %>%
        filter(type == "Fish")
      
      if (nrow(fish_ts) > 0) {
        p3 <- ggplot(fish_ts, aes(x = time, y = biomass, color = group)) +
          geom_line(linewidth = 0.8) +
          facet_wrap(~model, scales = "free_y") +
          scale_color_manual(values = group_colors[names(group_colors) %in% unique(fish_ts$group)]) +
          scale_y_log10() +
          labs(title = paste("Fish Group Biomass Time Series:", gsub("_", " ", toupper(scenario_name))),
               x = "Time (years)", 
               y = "Biomass (g, log scale)",
               color = "Fish Group") +
          custom_theme +
          theme(legend.position = "right")
      }
      
      # 4. Zooplankton functional groups
      zoo_ts <- combined_ts_data %>%
        filter(type == "Zooplankton")
      
      if (nrow(zoo_ts) > 0) {
        p4 <- ggplot(zoo_ts, aes(x = time, y = biomass, color = group)) +
          geom_line(linewidth = 0.8) +
          facet_wrap(~model, scales = "free_y") +
          scale_color_manual(values = group_colors[names(group_colors) %in% unique(zoo_ts$group)]) +
          scale_y_log10() +
          labs(title = paste("Zooplankton Group Biomass Time Series:", gsub("_", " ", toupper(scenario_name))),
               x = "Time (years)", 
               y = "Biomass (g, log scale)",
               color = "Zooplankton Group") +
          custom_theme +
          theme(legend.position = "right")
      }
      
      plots[[paste(scenario_name, "community", sep = "_")]] <- p1
      plots[[paste(scenario_name, "fish_zoo", sep = "_")]] <- p2
      if (exists("p3")) plots[[paste(scenario_name, "fish_groups", sep = "_")]] <- p3
      if (exists("p4")) plots[[paste(scenario_name, "zoo_groups", sep = "_")]] <- p4
    }
  }
  
  return(plots)
}

# Run diagnostics
cat("Running reproduction diagnostics...\n")
diagnose_reproduction_data()

# Generate plots
cat("\nGenerating diagnostic and time series plots...\n")

# Create reproduction plots
reproduction_plots <- create_proper_reproduction_plots()

# Create time series plots
time_series_plots <- create_comprehensive_time_series()

# Save all plots
cat("Saving diagnostic and time series plots...\n")

# Reproduction plots
for (plot_name in names(reproduction_plots)) {
  filename <- paste0("ZooMSS_Diagnostic_TimeSeries/reproduction_", plot_name, ".png")
  ggsave(filename, reproduction_plots[[plot_name]], width = 10, height = 6, dpi = 300)
}

# Time series plots
for (plot_name in names(time_series_plots)) {
  filename <- paste0("ZooMSS_Diagnostic_TimeSeries/timeseries_", plot_name, ".png")
  ggsave(filename, time_series_plots[[plot_name]], width = 12, height = 8, dpi = 300)
}

# Save plot objects
save(reproduction_plots, time_series_plots, 
     file = "ZooMSS_Diagnostic_TimeSeries/diagnostic_timeseries_plots.RData")

cat("\n=== DIAGNOSTIC AND TIME SERIES ANALYSIS COMPLETE ===\n")
cat("Generated plots:\n")
cat("✓ Reproduction diagnostic plots with time series and size allocation\n")
cat("✓ Total community biomass time series\n")
cat("✓ Fish vs zooplankton biomass time series\n")
cat("✓ Individual functional group time series\n")
cat("✓ All plots saved in 'ZooMSS_Diagnostic_TimeSeries/' directory\n")