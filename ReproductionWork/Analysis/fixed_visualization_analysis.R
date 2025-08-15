## Fixed Visualization Analysis for ZooMSS Reproduction Comparison
## Addresses library issues and creates proper visualizations

# Load required libraries
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(reshape2)) install.packages("reshape2")
if (!require(gridExtra)) install.packages("gridExtra")
if (!require(viridis)) install.packages("viridis")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")
if (!require(purrr)) install.packages("purrr")

library(ggplot2)
library(reshape2)
library(gridExtra)
library(viridis)
library(dplyr)
library(tidyr)
library(purrr)

# Create fixed figures directory
if (!dir.exists("ZooMSS_Reproduction_Figures_Fixed")) {
  dir.create("ZooMSS_Reproduction_Figures_Fixed")
}

# Load results
if (!file.exists("comprehensive_reproduction_comparison_results.RData")) {
  stop("Results file not found. Please run comprehensive_reproduction_comparison.R first.")
}

load("comprehensive_reproduction_comparison_results.RData")
cat("=== FIXED VISUALIZATION ANALYSIS ===\n\n")

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
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  )

# Function to create proper size spectrum plots
create_fixed_size_spectrum_plots <- function() {
  cat("Creating fixed size spectrum plots...\n")
  
  plots <- list()
  
  for (scenario_name in names(chlorophyll_scenarios)) {
    spectrum_data <- data.frame()
    
    # Get results for this scenario
    for (result_name in names(results)) {
      if (grepl(scenario_name, result_name)) {
        result <- results[[result_name]]
        
        # Calculate total abundance and biomass by size class
        total_abundance_by_size <- colSums(result$abundances)
        total_biomass_by_size <- colSums(sweep(result$abundances, 2, result$model$param$w, "*"))
        
        log_weights <- log10(result$model$param$w)
        
        temp_data <- data.frame(
          log_weight = log_weights,
          abundance = total_abundance_by_size,
          biomass = total_biomass_by_size,
          model = result$version,
          scenario = scenario_name
        )
        
        spectrum_data <- rbind(spectrum_data, temp_data)
      }
    }
    
    # Remove zero values for log plotting
    spectrum_data$abundance[spectrum_data$abundance <= 0] <- 1e-12
    spectrum_data$biomass[spectrum_data$biomass <= 0] <- 1e-12
    
    # Abundance spectrum plot
    p1 <- ggplot(spectrum_data, aes(x = log_weight, y = log10(abundance), color = model)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = model_colors) +
      labs(title = paste("Size Spectrum - Abundance:", gsub("_", " ", toupper(scenario_name))),
           x = "Log10 Body Weight (g)", 
           y = "Log10 Abundance (individuals)",
           color = "Model Version") +
      custom_theme
    
    # Biomass spectrum plot
    p2 <- ggplot(spectrum_data, aes(x = log_weight, y = log10(biomass), color = model)) +
      geom_line(linewidth = 1.2) +
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

# Function to create fixed growth rate plots
create_fixed_growth_rate_plots <- function() {
  cat("Creating fixed growth rate plots...\n")
  
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
            
            # Find size range for this group
            w0_log <- result$groups_used$W0[i]
            wmax_log <- result$groups_used$Wmax[i]
            
            # Only include size classes within the group's range
            size_mask <- log_weights >= w0_log & log_weights <= wmax_log
            
            if (any(size_mask)) {
              temp_data <- data.frame(
                log_weight = log_weights[size_mask],
                growth_rate = growth_rates[i, size_mask],
                group = group_name,
                type = group_type,
                model = result$version,
                scenario = scenario_name
              )
              
              # Remove any infinite or NA values
              temp_data <- temp_data[is.finite(temp_data$growth_rate), ]
              
              if (nrow(temp_data) > 0) {
                growth_data <- rbind(growth_data, temp_data)
              }
            }
          }
        }
      }
    }
    
    if (nrow(growth_data) > 0) {
      # Growth rate plot
      p <- ggplot(growth_data, aes(x = log_weight, y = growth_rate, color = model)) +
        geom_line(linewidth = 1) +
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

# Function to create fixed reproduction allocation plots
create_fixed_reproduction_plots <- function() {
  cat("Creating fixed reproduction allocation plots...\n")
  
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
        w0_log <- result$groups_used$W0[i]
        wmax_log <- result$groups_used$Wmax[i]
        wmat_log <- result$groups_used$Wmat[i]
        
        # Only include size classes within the group's range
        size_mask <- log_weights >= w0_log & log_weights <= wmax_log
        
        if (any(size_mask)) {
          temp_data <- data.frame(
            log_weight = log_weights[size_mask],
            reproduction_rate = avg_repro_rates[i, size_mask],
            group = group_name,
            scenario = scenario,
            maturity_size = wmat_log
          )
          
          # Remove any infinite or NA values
          temp_data <- temp_data[is.finite(temp_data$reproduction_rate), ]
          
          if (nrow(temp_data) > 0) {
            repro_data <- rbind(repro_data, temp_data)
          }
        }
      }
      
      if (nrow(repro_data) > 0) {
        # Reproduction allocation plot with maturity lines
        p <- ggplot(repro_data, aes(x = log_weight, y = reproduction_rate, color = group)) +
          geom_line(linewidth = 1.2) +
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

# Function to create fixed diet plots
create_fixed_diet_plots <- function() {
  cat("Creating fixed diet composition plots...\n")
  
  plots <- list()
  
  for (scenario_name in names(chlorophyll_scenarios)) {
    diet_data <- data.frame()
    
    for (result_name in names(results)) {
      if (grepl(scenario_name, result_name)) {
        result <- results[[result_name]]
        
        # Get diet data - should be predator x prey matrix
        if (!is.null(result$diets)) {
          diets <- result$diets
          
          # Create diet data frame manually
          predator_names <- result$groups_used$Species
          prey_names <- c("Pico_Phyto", "Nano_Phyto", "Micro_Phyto", predator_names)
          
          # Ensure we don't exceed the number of columns in diets
          n_prey <- min(ncol(diets), length(prey_names))
          
          for (i in 1:nrow(diets)) {
            for (j in 1:n_prey) {
              if (diets[i, j] > 0.01) {  # Only include significant diet components
                temp_data <- data.frame(
                  Predator = predator_names[i],
                  Prey = prey_names[j],
                  Diet = diets[i, j],
                  model = result$version,
                  scenario = scenario_name
                )
                diet_data <- rbind(diet_data, temp_data)
              }
            }
          }
        }
      }
    }
    
    if (nrow(diet_data) > 0) {
      # Normalize diet proportions within each predator-model combination
      diet_data <- diet_data %>%
        group_by(Predator, model) %>%
        mutate(Diet_prop = Diet / sum(Diet)) %>%
        ungroup()
      
      # Diet composition plot
      p <- ggplot(diet_data, aes(x = Predator, y = Diet_prop, fill = Prey)) +
        geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
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

# Function to create fixed summary plots
create_fixed_summary_plots <- function() {
  cat("Creating fixed summary comparison plots...\n")
  
  # Final biomass comparison
  p1 <- ggplot(comparison_table, aes(x = Scenario, y = Fish_Biomass, fill = Model_Version)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = model_colors) +
    scale_y_log10() +
    labs(title = "Final Fish Biomass Comparison",
         x = "Chlorophyll Scenario", 
         y = "Fish Biomass (g, log scale)",
         fill = "Model Version") +
    custom_theme
  
  p2 <- ggplot(comparison_table, aes(x = Scenario, y = Total_Biomass, fill = Model_Version)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = model_colors) +
    scale_y_log10() +
    labs(title = "Final Total Biomass Comparison",
         x = "Chlorophyll Scenario", 
         y = "Total Biomass (g, log scale)",
         fill = "Model Version") +
    custom_theme
  
  p3 <- ggplot(comparison_table, aes(x = Scenario, y = Fish_Zoo_Ratio, fill = Model_Version)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = model_colors) +
    labs(title = "Fish/Zooplankton Biomass Ratio Comparison",
         x = "Chlorophyll Scenario", 
         y = "Fish/Zooplankton Ratio",
         fill = "Model Version") +
    custom_theme
  
  return(list(fish_biomass = p1, total_biomass = p2, fish_zoo_ratio = p3))
}

# Generate all fixed plots
cat("Generating fixed visualizations...\n\n")

# Create all plot sets
size_spectrum_plots <- create_fixed_size_spectrum_plots()
growth_rate_plots <- create_fixed_growth_rate_plots()
reproduction_plots <- create_fixed_reproduction_plots()
diet_plots <- create_fixed_diet_plots()
summary_plots <- create_fixed_summary_plots()

# Save all plots with high quality
cat("Saving fixed plots to ZooMSS_Reproduction_Figures_Fixed/...\n")

# Size spectrum plots
for (plot_name in names(size_spectrum_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures_Fixed/size_spectrum_", plot_name, ".png")
  ggsave(filename, size_spectrum_plots[[plot_name]], width = 10, height = 6, dpi = 300)
}

# Growth rate plots
for (plot_name in names(growth_rate_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures_Fixed/growth_rates_", plot_name, ".png")
  ggsave(filename, growth_rate_plots[[plot_name]], width = 12, height = 8, dpi = 300)
}

# Reproduction plots
for (plot_name in names(reproduction_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures_Fixed/reproduction_allocation_", plot_name, ".png")
  ggsave(filename, reproduction_plots[[plot_name]], width = 10, height = 6, dpi = 300)
}

# Diet plots
for (plot_name in names(diet_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures_Fixed/diet_composition_", plot_name, ".png")
  ggsave(filename, diet_plots[[plot_name]], width = 12, height = 8, dpi = 300)
}

# Summary plots
for (plot_name in names(summary_plots)) {
  filename <- paste0("ZooMSS_Reproduction_Figures_Fixed/summary_", plot_name, ".png")
  ggsave(filename, summary_plots[[plot_name]], width = 10, height = 6, dpi = 300)
}

# Create a combined summary figure
combined_summary <- grid.arrange(
  summary_plots$fish_biomass,
  summary_plots$total_biomass,
  summary_plots$fish_zoo_ratio,
  ncol = 1
)

ggsave("ZooMSS_Reproduction_Figures_Fixed/combined_summary.png", combined_summary, 
       width = 10, height = 12, dpi = 300)

# Save plot objects for further analysis
save(size_spectrum_plots, growth_rate_plots, reproduction_plots, diet_plots, 
     summary_plots, 
     file = "ZooMSS_Reproduction_Figures_Fixed/all_fixed_plots.RData")

cat("\n=== FIXED VISUALIZATION ANALYSIS COMPLETE ===\n")
cat("Generated fixed figures:\n")
cat("✓ Size spectrum plots with proper abundance and biomass calculations\n")
cat("✓ Growth rate plots with correct size ranges for each group\n")
cat("✓ Reproduction allocation plots with maturity size indicators\n")
cat("✓ Diet composition plots with simplified prey categories\n")
cat("✓ Summary comparison plots with theme_bw()\n")
cat("✓ Combined summary figure\n")
cat("\nAll fixed figures saved in 'ZooMSS_Reproduction_Figures_Fixed/' directory\n")
cat("Plot objects saved in 'all_fixed_plots.RData' for further analysis\n")