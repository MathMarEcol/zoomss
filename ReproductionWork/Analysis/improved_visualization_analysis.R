## Improved Visualization Analysis for ZooMSS Reproduction Comparison
## Uses existing ZooMSS plotting functions and fixes visualization issues

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

# Source the ZooMSS extras functions
source("fZooMSS_Xtras.R")

# Create improved figures directory
if (!dir.exists("ZooMSS_Reproduction_Figures_Improved")) {
  dir.create("ZooMSS_Reproduction_Figures_Improved")
}

# Load results
if (!file.exists("comprehensive_reproduction_comparison_results.RData")) {
  stop("Results file not found. Please run comprehensive_reproduction_comparison.R first.")
}

load("comprehensive_reproduction_comparison_results.RData")
cat("=== IMPROVED VISUALIZATION ANALYSIS ===\n\n")

# Set up improved color schemes and theme
model_colors <- c("Old_ZooMSS" = "#1f77b4", 
                  "New_ZooMSS_No_Repro" = "#ff7f0e", 
                  "New_ZooMSS_With_Repro" = "#d62728")

scenario_colors <- c("low" = "#440154", "medium" = "#21908C", "high" = "#FDE725")

# Use theme_bw() as requested
custom_theme <- theme_bw() + 
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
  )

# Function to create proper size spectrum plots using ZooMSS functions
create_proper_size_spectrum_plots <- function() {
  cat("Creating improved size spectrum plots...\n")
  
  plots <- list()
  
  for (scenario_name in names(chlorophyll_scenarios)) {
    scenario_results <- list()
    
    # Get results for this scenario
    for (result_name in names(results)) {
      if (grepl(scenario_name, result_name)) {
        scenario_results[[results[[result_name]]$version]] <- results[[result_name]]
      }
    }
    
    # Create combined size spectrum data
    spectrum_data <- data.frame()
    
    for (model_name in names(scenario_results)) {
      result <- scenario_results[[model_name]]
      
      # Use ZooMSS biomass function
      biomass_by_size <- fZooMSS_SizeBiomass(list(result$abundances), result$model$param$w)[[1]]
      abundance_by_size <- apply(result$abundances, 2, sum)
      
      log_weights <- log10(result$model$param$w)
      
      temp_data <- data.frame(
        log_weight = log_weights,
        abundance = abundance_by_size,
        biomass = biomass_by_size,
        model = model_name,
        scenario = scenario_name
      )
      
      spectrum_data <- rbind(spectrum_data, temp_data)
    }
    
    # Abundance spectrum plot
    p1 <- ggplot(spectrum_data, aes(x = log_weight, y = log10(abundance + 1e-12), 
                                   color = model)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = model_colors) +
      labs(title = paste("Size Spectrum - Abundance:", gsub("_", " ", toupper(scenario_name))),
           x = "Log10 Body Weight (g)", 
           y = "Log10 Abundance (individuals)",
           color = "Model Version") +
      custom_theme
    
    # Biomass spectrum plot
    p2 <- ggplot(spectrum_data, aes(x = log_weight, y = log10(biomass + 1e-12), 
                                   color = model)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = model_colors) +
      labs(title = paste("Size Spectrum - Biomass:", gsub("_", " ", toupper(scenario_name))),
           x = "Log10 Body Weight (g)", 
           y = "Log10 Biomass (g)",
           color = "Model Version") +
      custom_theme
    
    plots[[paste(scenario_name, "abundance", sep = "_")]] <- p1
    plots[[paste(scenario_name, "biomass", sep = "_")]] <- p2
  }
  
  return(plots)
}

# Function to create proper growth rate plots
create_proper_growth_rate_plots <- function() {
  cat("Creating improved growth rate plots...\n")
  
  plots <- list()
  
  for (scenario_name in names(chlorophyll_scenarios)) {
    growth_data <- data.frame()
    
    for (result_name in names(results)) {
      if (grepl(scenario_name, result_name)) {
        result <- results[[result_name]]
        
        # Get growth rates properly
        if (!is.null(result$growth)) {
          growth_rates <- result$growth
          weights <- result$model$param$w
          log_weights <- log10(weights)
          
          # Create data for each group
          for (i in 1:nrow(growth_rates)) {
            group_name <- result$groups_used$Species[i]
            group_type <- result$groups_used$Type[i]
            
            # Only include size classes within the group's range
            w0_idx <- which.min(abs(log_weights - result$groups_used$W0[i]))
            wmax_idx <- which.min(abs(log_weights - result$groups_used$Wmax[i]))
            
            if (w0_idx <= wmax_idx) {
              size_range <- w0_idx:wmax_idx
              
              temp_data <- data.frame(
                log_weight = log_weights[size_range],
                growth_rate = growth_rates[i, size_range],
                group = group_name,
                type = group_type,
                model = result$version,
                scenario = scenario_name
              )
              
              growth_data <- rbind(growth_data, temp_data)
            }
          }
        }
      }
    }
    
    if (nrow(growth_data) > 0) {
      # Growth rate plot
      p <- ggplot(growth_data, aes(x = log_weight, y = growth_rate, 
                                  color = model)) +
        geom_line(size = 1) +
        facet_wrap(~group, scales = "free", ncol = 3) +
        scale_color_manual(values = model_colors) +
        labs(title = paste("Growth Rates:", gsub("_", " ", toupper(scenario_name))),
             x = "Log10 Body Weight (g)", 
             y = "Growth Rate (per year)",
             color = "Model Version") +
        custom_theme
      
      plots[[scenario_name]] <- p
    }
  }
  
  return(plots)
}

# Function to create proper reproduction allocation plots
create_proper_reproduction_plots <- function() {
  cat("Creating improved reproduction allocation plots...\n")
  
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
      # Get average reproduction rates over last 20% of simulation
      n_times <- dim(result$model$reproduction_rate)[1]
      last_20_percent <- floor(0.8 * n_times):n_times
      
      avg_repro_rates <- apply(result$model$reproduction_rate[last_20_percent, , ], c(2, 3), mean, na.rm = TRUE)
      weights <- result$model$param$w
      log_weights <- log10(weights)
      
      # Create data frame for fish groups only
      fish_groups <- which(result$groups_used$Type == "Fish")
      repro_data <- data.frame()
      
      for (i in fish_groups) {
        group_name <- result$groups_used$Species[i]
        
        # Only include size classes within the group's range
        w0_idx <- which.min(abs(log_weights - result$groups_used$W0[i]))
        wmax_idx <- which.min(abs(log_weights - result$groups_used$Wmax[i]))
        
        if (w0_idx <= wmax_idx) {
          size_range <- w0_idx:wmax_idx
          
          temp_data <- data.frame(
            log_weight = log_weights[size_range],
            reproduction_rate = avg_repro_rates[i, size_range],
            group = group_name,
            scenario = scenario,
            maturity_size = result$groups_used$Wmat[i]
          )
          
          repro_data <- rbind(repro_data, temp_data)
        }
      }
      
      if (nrow(repro_data) > 0) {
        # Reproduction allocation plot with maturity lines
        p <- ggplot(repro_data, aes(x = log_weight, y = reproduction_rate, color = group)) +
          geom_line(size = 1.2) +
          geom_vline(aes(xintercept = maturity_size, color = group), 
                     linetype = "dashed", alpha = 0.7) +
          scale_color_brewer(type = "qual", palette = "Set1") +
          labs(title = paste("Reproduction Allocation:", gsub("_", " ", toupper(scenario))),
               x = "Log10 Body Weight (g)", 
               y = "Reproduction Rate",
               color = "Fish Group",
               caption = "Dashed lines show maturity sizes") +
          custom_theme
        
        plots[[scenario]] <- p
      }
    }
  }
  
  return(plots)
}

# Function to create proper diet plots with correct prey categories
create_proper_diet_plots <- function() {
  cat("Creating improved diet composition plots...\n")
  
  plots <- list()
  
  for (scenario_name in names(chlorophyll_scenarios)) {
    diet_data <- data.frame()
    
    for (result_name in names(results)) {
      if (grepl(scenario_name, result_name)) {
        result <- results[[result_name]]
        
        # Get diet data - should be predator x prey matrix
        if (!is.null(result$diets)) {
          diets <- result$diets
          
          # Create proper diet tibble using ZooMSS function
          diet_tibble <- fZooMSS_MakeDietTibble(diets, result$model)
          diet_tibble$model <- result$version
          diet_tibble$scenario <- scenario_name
          
          diet_data <- rbind(diet_data, diet_tibble)
        }
      }
    }
    
    if (nrow(diet_data) > 0) {
      # Filter out very small diet fractions for clarity
      diet_data_filtered <- diet_data %>%
        filter(Diet > 0.01) %>%  # Only show prey that make up >1% of diet
        group_by(Predator, model) %>%
        mutate(Diet_prop = Diet / sum(Diet)) %>%  # Renormalize
        ungroup()
      
      # Diet composition plot
      p <- ggplot(diet_data_filtered, aes(x = Predator, y = Diet_prop, fill = Prey)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_wrap(~model, scales = "free_x") +
        scale_fill_viridis_d(option = "plasma") +
        labs(title = paste("Diet Composition:", gsub("_", " ", toupper(scenario_name))),
             x = "Predator Group", 
             y = "Diet Proportion",
             fill = "Prey Type") +
        custom_theme +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      plots[[scenario_name]] <- p
    }
  }
  
  return(plots)
}

# Function to create proper time series plots
create_proper_time_series_plots <- function() {
  cat("Creating improved time series plots...\n")
  
  plots <- list()
  
  for (scenario_name in names(chlorophyll_scenarios)) {
    ts_data <- data.frame()
    
    for (result_name in names(results)) {
      if (grepl(scenario_name, result_name)) {
        result <- results[[result_name]]
        
        if (!is.null(result$model$N)) {
          # Use ZooMSS functions to calculate biomass time series
          n_times <- dim(result$model$N)[1]
          
          # Convert to list format for ZooMSS functions
          N_list <- list()
          for (t in 1:n_times) {
            N_list[[t]] <- result$model$N[t, , ]
          }
          
          # Calculate species biomass time series
          species_biomass_ts <- fZooMSS_SpeciesBiomass(N_list, result$model)
          
          # Convert to data frame
          biomass_matrix <- do.call(rbind, species_biomass_ts)
          
          # Time vector
          time_years <- seq(0, result$model$param$tmax, length.out = n_times)
          
          # Separate fish and zooplankton
          fish_groups <- which(result$groups_used$Type == "Fish")
          zoo_groups <- which(result$groups_used$Type == "Zooplankton")
          
          fish_biomass <- rowSums(biomass_matrix[, fish_groups, drop = FALSE])
          zoo_biomass <- rowSums(biomass_matrix[, zoo_groups, drop = FALSE])
          
          temp_data <- data.frame(
            time = time_years,
            fish_biomass = fish_biomass,
            zoo_biomass = zoo_biomass,
            total_biomass = fish_biomass + zoo_biomass,
            model = result$version,
            scenario = scenario_name
          )
          
          ts_data <- rbind(ts_data, temp_data)
        }
      }
    }
    
    if (nrow(ts_data) > 0) {
      # Fish biomass time series
      p1 <- ggplot(ts_data, aes(x = time, y = fish_biomass, color = model)) +
        geom_line(size = 1) +
        scale_color_manual(values = model_colors) +
        scale_y_log10() +
        labs(title = paste("Fish Biomass Time Series:", gsub("_", " ", toupper(scenario_name))),
             x = "Time (years)", 
             y = "Fish Biomass (g, log scale)",
             color = "Model Version") +
        custom_theme
      
      # Total biomass time series
      p2 <- ggplot(ts_data, aes(x = time, y = total_biomass, color = model)) +
        geom_line(size = 1) +
        scale_color_manual(values = model_colors) +
        scale_y_log10() +
        labs(title = paste("Total Biomass Time Series:", gsub("_", " ", toupper(scenario_name))),
             x = "Time (years)", 
             y = "Total Biomass (g, log scale)",
             color = "Model Version") +
        custom_theme
      
      plots[[paste(scenario_name, "fish", sep = "_")]] <- p1
      plots[[paste(scenario_name, "total", sep = "_")]] <- p2
    }
  }
  
  return(plots)
}

# Function to create improved summary plots
create_improved_summary_plots <- function() {
  cat("Creating improved summary comparison plots...\n")
  
  # Final biomass comparison
  p1 <- ggplot(comparison_table, aes(x = Scenario, y = Fish_Biomass, fill = Model_Version)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.3) +
    scale_fill_manual(values = model_colors) +
    scale_y_log10() +
    labs(title = "Final Fish Biomass Comparison",
         x = "Chlorophyll Scenario", 
         y = "Fish Biomass (g, log scale)",
         fill = "Model Version") +
    custom_theme
  
  p2 <- ggplot(comparison_table, aes(x = Scenario, y = Total_Biomass, fill = Model_Version)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.3) +
    scale_fill_manual(values = model_colors) +
    scale_y_log10() +
    labs(title = "Final Total Biomass Comparison",
         x = "Chlorophyll Scenario", 
         y = "Total Biomass (g, log scale)",
         fill = "Model Version") +
    custom_theme
  
  p3 <- ggplot(comparison_table, aes(x = Scenario, y = Fish_Zoo_Ratio, fill = Model_Version)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.3) +
    scale_fill_manual(values = model_colors) +
    labs(title = "Fish/Zooplankton Biomass Ratio Comparison",
         x = "Chlorophyll Scenario", 
         y = "Fish/Zooplankton Ratio",
         fill = "Model Version") +
    custom_theme
  
  return(list(fish_biomass = p1, total_biomass = p2, fish_zoo_ratio = p3))
}

# Generate all improved plots
cat("Generating improved visualizations...\n\n")

# Create all plot sets
size_spectrum_plots <- create_proper_size_spectrum_plots()
growth_rate_plots <- create_proper_growth_rate_plots()
reproduction_plots <- create_proper_reproduction_plots()
diet_plots <- create_proper_diet_plots()
time_series_plots <- create_proper_time_series_plots()
summary_plots <- create_improved_summary_plots()

# Save all plots with high quality
cat("Saving improved plots to ZooMSS_Reproduction_Figures_Improved/...\n")

# Size spectrum plots
for (plot_name in names(size_spectrum_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures_Improved/size_spectrum_", plot_name, ".png")
  ggsave(filename, size_spectrum_plots[[plot_name]], width = 10, height = 6, dpi = 300)
}

# Growth rate plots
for (plot_name in names(growth_rate_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures_Improved/growth_rates_", plot_name, ".png")
  ggsave(filename, growth_rate_plots[[plot_name]], width = 12, height = 8, dpi = 300)
}

# Reproduction plots
for (plot_name in names(reproduction_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures_Improved/reproduction_allocation_", plot_name, ".png")
  ggsave(filename, reproduction_plots[[plot_name]], width = 10, height = 6, dpi = 300)
}

# Diet plots
for (plot_name in names(diet_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures_Improved/diet_composition_", plot_name, ".png")
  ggsave(filename, diet_plots[[plot_name]], width = 12, height = 8, dpi = 300)
}

# Time series plots
for (plot_name in names(time_series_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures_Improved/time_series_", plot_name, ".png")
  ggsave(filename, time_series_plots[[plot_name]], width = 10, height = 6, dpi = 300)
}

# Summary plots
for (plot_name in names(summary_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures_Improved/summary_", plot_name, ".png")
  ggsave(filename, summary_plots[[plot_name]], width = 10, height = 6, dpi = 300)
}

# Create a combined summary figure
combined_summary <- grid.arrange(
  summary_plots$fish_biomass,
  summary_plots$total_biomass,
  summary_plots$fish_zoo_ratio,
  ncol = 1
)

ggsave("ZooMSS_Reproduction_Figures_Improved/combined_summary.png", combined_summary, 
       width = 10, height = 12, dpi = 300)

# Save plot objects for further analysis
save(size_spectrum_plots, growth_rate_plots, reproduction_plots, diet_plots, 
     time_series_plots, summary_plots, 
     file = "ZooMSS_Reproduction_Figures_Improved/all_improved_plots.RData")

cat("\n=== IMPROVED VISUALIZATION ANALYSIS COMPLETE ===\n")
cat("Generated improved figures:\n")
cat("✓ Size spectrum plots using proper ZooMSS functions\n")
cat("✓ Growth rate plots with correct size ranges for each group\n")
cat("✓ Reproduction allocation plots with maturity size indicators\n")
cat("✓ Diet composition plots using ZooMSS diet functions\n")
cat("✓ Time series plots with proper biomass calculations\n")
cat("✓ Summary comparison plots with theme_bw()\n")
cat("✓ Combined summary figure\n")
cat("\nAll improved figures saved in 'ZooMSS_Reproduction_Figures_Improved/' directory\n")
cat("Plot objects saved in 'all_improved_plots.RData' for further analysis\n")