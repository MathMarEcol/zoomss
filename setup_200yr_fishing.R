## ZooMSS 200-year simulation with fishing mortality for fish groups
## This script runs a ZooMSS simulation for 200 years using default values
## except with added fishing mortality for the three fish groups
##
## Code written by Dr Jason Everett (UQ/UNSW/CSIRO), Dr Ryan Heneghan (QUT) and Mr Patrick Sykes (UQ)
## Modified for 200-year fishing simulation

# library(Rcpp) # Only needed if we are running with Rcpp code.
source("fZooMSS_Model.R") #source the model code
source("fZooMSS_Xtras.R")


# Setup user defined parameters -------------------------------------------

# Create environmental data using default values
# Using typical oceanic conditions: SST = 15°C, Chlorophyll = 0.5 mg/m³
enviro_data <- fZooMSS_CalculatePhytoParam(data.frame(cellID = 1,
                                                      sst = 15,
                                                      chlo = 0.5))

# Add delta time (years) - using smaller timestep for longer simulation
enviro_data$dt <- 0.01

# Set length of simulation to 200 years
enviro_data$tmax <- 200

# Setup jobname
jobname <- "20250806_ZooMSS_200yr_Fishing"  # Job name for this 200-year fishing simulation
enviro_row <- 1 # Which row of the environmental data to use

HPC <- FALSE # Running locally, not on HPC
SaveTimeSteps <- TRUE # Save all time steps for analysis

# Load functional group information with fishing mortality for fish groups
Groups <- read.csv("TestGroups_WithFishing.csv", stringsAsFactors = FALSE)

# Print fishing mortality values for fish groups to confirm setup
cat("Fishing mortality rates applied:\n")
fish_groups <- Groups[Groups$Type == "Fish", ]
for(i in 1:nrow(fish_groups)){
  cat(sprintf("%s: F = %0.2f yr^-1 (sizes %s to %s log10(g))\n", 
              fish_groups$Species[i], 
              fish_groups$Fmort[i],
              fish_groups$Fmort_W0[i],
              fish_groups$Fmort_Wmax[i]))
}

# No need to change anything below here. ----------------------------------

if (HPC == TRUE){
  ID <- as.integer(Sys.getenv('PBS_ARRAY_INDEX')) # Get the array run number on HPC
  } else {
    ID <- enviro_row
  }

ID_char <- sprintf("%04d",ID) # Set the ID as a 4 digit character so it will sort properly

input_params <- enviro_data[ID,]
rm(enviro_data)

# Create RawOutput directory if it doesn't exist
if (!dir.exists("RawOutput")) {
  dir.create("RawOutput")
}

cat("Starting ZooMSS 200-year simulation with fishing mortality...\n")
cat("This may take several minutes to complete.\n")

# Run the model and time it
model_runtime <- system.time(
  out <- fZooMSS_Model(input_params, Groups, SaveTimeSteps)
)

out$model$model_runtime <- model_runtime

cat(sprintf("Model completed in %0.2f seconds\n", model_runtime[3]))

# Save the output
output_file <- paste0("RawOutput/", jobname, "_", ID_char,".RDS")
saveRDS(out, file = output_file)
cat(sprintf("Output saved to: %s\n", output_file))

# Plotting ----------------------------------------------------------------

source("fZooMSS_Plot.R")

cat("Generating plots...\n")

# Generate standard plots
(ggPPMR <- fZooMSS_Plot_PPMR(out))
(ggSizeSpec <- fZooMSS_Plot_SizeSpectra(out))

## Since we saved timesteps, plot the timeseries
(ggAbundTS <- fZooMSS_Plot_AbundTimeSeries(out))
(ggGrowthTS <- fZooMSS_Plot_GrowthTimeSeries(out))
(ggPredTS <- fZooMSS_Plot_PredTimeSeries(out))

# Save plots
ggsave("ZooMSS_200yr_Fishing_PPMR.png", ggPPMR, width = 10, height = 8, dpi = 300)
ggsave("ZooMSS_200yr_Fishing_SizeSpectra.png", ggSizeSpec, width = 10, height = 8, dpi = 300)
ggsave("ZooMSS_200yr_Fishing_AbundanceTimeSeries.png", ggAbundTS, width = 12, height = 8, dpi = 300)
ggsave("ZooMSS_200yr_Fishing_GrowthTimeSeries.png", ggGrowthTS, width = 12, height = 8, dpi = 300)
ggsave("ZooMSS_200yr_Fishing_PredationTimeSeries.png", ggPredTS, width = 12, height = 8, dpi = 300)

cat("Plots saved as PNG files.\n")
cat("Simulation complete!\n")

# Print summary of final state
cat("\n=== SIMULATION SUMMARY ===\n")
cat(sprintf("Duration: %d years\n", input_params$tmax))
cat(sprintf("Time step: %0.3f years\n", input_params$dt))
cat(sprintf("Total runtime: %0.2f seconds\n", model_runtime[3]))
cat(sprintf("Output file: %s\n", output_file))
cat("Fishing mortality applied to:\n")
for(i in 1:nrow(fish_groups)){
  cat(sprintf("  %s: F = %0.2f yr^-1\n", fish_groups$Species[i], fish_groups$Fmort[i]))
}
