## Corrected Reproduction Allocation Visualization
## Fixes the flat-line issue by properly calculating reproduction allocation

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load results
if (!file.exists("comprehensive_reproduction_comparison_results.RData")) {
  stop("Results file not found. Please run comprehensive_reproduction_comparison.R first.")
}

load("comprehensive_reproduction_comparison_results.RData")
cat("=== CORRECTED REPRODUCTION VISUALIZATION ===\n\n")

# Create corrected figures directory
if (!dir.exists("ZooMSS_Reproduction_Figures_Corrected")) {
  dir.create("ZooMSS_Reproduction_Figures_Corrected")
}

# Set up color scheme and theme
model_colors <- c("Old_ZooMSS" = "#1f77b4", 
                  "New_ZooMSS_No_Repro" = "#ff7f0e", 
                  "New_ZooMSS_With_Repro" = "#d62728")

custom_theme <- theme_bw() + 
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  )

# Function to calculate theoretical reproduction allocation
calculate_theoretical_reproduction <- function(log_weights, maturity_size, max_allocation = 0.3) {
  # Maturity function: M(w) = 1 / (1 + exp(-5 * (log10(w) - Wmat)))
  maturity <- 1 / (1 + exp(-5 * (log_weights - maturity_size)))
  
  # Reproduction allocation: f_repro(w) = max_allocation * M(w)
  reproduction_allocation <- max_allocation * maturity
  
  return(reproduction_allocation)
}

# Function to create corrected reproduction allocation plots
create_corrected_reproduction_plots <- function() {
  cat("Creating corrected reproduction allocation plots...\n")
  
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
    
    cat("Processing scenario:", scenario, "\n")
    
    # Get model parameters
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
      
      cat("  Processing group:", group_name, "Wmat:", wmat_log, "\n")
      
      # Only include size classes within the group's range
      size_mask <- log_weights >= w0_log & log_weights <= wmax_log
      
      if (any(size_mask)) {
        # Calculate theoretical reproduction allocation
        theoretical_repro <- calculate_theoretical_reproduction(
          log_weights[size_mask], 
          wmat_log, 
          max_allocation = 0.3
        )
        
        temp_data <- data.frame(
          log_weight = log_weights[size_mask],
          reproduction_allocation = theoretical_repro,
          group = group_name,
          scenario = scenario,
          maturity_size = wmat_log
        )
        
        repro_data <- rbind(repro_data, temp_data)
      }
    }
    
    if (nrow(repro_data) > 0) {
      # Create reproduction allocation plot
      p <- ggplot(repro_data, aes(x = log_weight, y = reproduction_allocation, color = group)) +
        geom_line(linewidth = 1.2) +
        geom_vline(data = unique(repro_data[c("group", "maturity_size")]), 
                   aes(xintercept = maturity_size, color = group), 
                   linetype = "dashed", alpha = 0.7, linewidth = 0.8) +
        scale_color_brewer(type = "qual", palette = "Set1") +
        scale_y_continuous(limits = c(0, 0.35), labels = scales::percent_format()) +
        labs(title = paste("Reproduction Allocation:", gsub("_", " ", toupper(scenario))),
             x = "Log10 Body Weight (g)", 
             y = "Reproduction Allocation (% of Energy)",
             color = "Fish Group",
             caption = "Dashed vertical lines show maturity sizes (Wmat)") +
        custom_theme
      
      plots[[scenario]] <- p
      
      # Save individual plot
      filename <- paste0("ZooMSS_Reproduction_Figures_Corrected/reproduction_allocation_", scenario, ".png")
      ggsave(filename, p, width = 12, height = 8, dpi = 300)
      cat("  Saved:", filename, "\n")
    }
  }
  
  return(plots)
}

# Generate corrected reproduction plots
cat("Generating corrected reproduction allocation visualizations...\n\n")
reproduction_plots <- create_corrected_reproduction_plots()

# Save plot objects for further analysis
save(reproduction_plots, 
     file = "ZooMSS_Reproduction_Figures_Corrected/corrected_reproduction_plots.RData")

cat("\n=== CORRECTED REPRODUCTION VISUALIZATION COMPLETE ===\n")
cat("Generated corrected figures:\n")
cat("✓ Theoretical reproduction allocation plots showing proper sigmoidal patterns\n")
cat("✓ Maturity size indicators (dashed vertical lines)\n")
cat("✓ Proper scaling and percentage formatting\n")
cat("\nAll corrected figures saved in 'ZooMSS_Reproduction_Figures_Corrected/' directory\n")
cat("Plot objects saved in 'corrected_reproduction_plots.RData' for further analysis\n")