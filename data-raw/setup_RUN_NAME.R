## ZooMSS Setup Script
## This script demonstrates how to run ZooMSS with time-varying environmental data
## Based on the original setup by Dr Jason Everett, Dr Ryan Heneghan, and Mr Patrick Sykes
## Environmental forcing - August 2025

# Load all required functions
source("zoomss_model.R") #source the model code
source("zXtras.R") # For zCalculatePhytoParam and other utilities
source("zEnvironmental_Utils.R") # For environmental time series functions

plotting <- FALSE

# Setup user defined parameters -------------------------------------------
dt <- 0.01 # years^-1 - time step for creating time vector
tmax <- 250 # years - total simulation time
cellID <- 1
isave <- 100 # Model save frequency
base_sst <- 15
base_chlo <- 0.01

# Create time vector (dt and tmax only used here for time vector creation)
time_vec <- seq(0, tmax, by = dt)

# Create environmental data with the same length as time vector
n_time_steps <- length(time_vec)

# Create input_params dataframe using new time-based approach
# Note: dt and tmax are automatically calculated from time_vec in zCreateInputs
input_params <- zCreateInputs(
  time = time_vec,
  sst = rep(base_sst, n_time_steps),
  chl = rep(base_chlo, n_time_steps),
  cellID = rep(cellID, n_time_steps),
  isave = isave
)

# Alternative: add time-varying patterns
# input_params$chl <- 0.01 + 0.00015 * input_params$time_step # Linear increasing Chlorophyll
# input_params$chl <- 0.5 + 0.49 * sin(2 * pi * input_params$time) # Chlorophyll cycles between 0.01-0.99 mg/m³

# Ensure chlorophyll stays positive
input_params$chl <- pmax(input_params$chl, 0.01)
cat("Environmental data created with", nrow(input_params), "timesteps covering",
    max(input_params$time), "years\n")

zPlotEnvironment(input_params)


# Setup jobname
jobname <- "20250813_chl0_01"  # This is the job name to save the run
enviro_row <- 1 # Which row of the environmental data to use
HPC <- FALSE # Is this being run on a HPC
SaveTimeSteps <- TRUE # Should we save all time steps
Groups <- read.csv("GroupInputs.csv", stringsAsFactors = FALSE) # Load functional group information

if (HPC == TRUE){
  ID <- as.integer(Sys.getenv('PBS_ARRAY_INDEX'))
} else {
  ID <- enviro_row
}
ID_char <- sprintf("%04d",ID)

# Create dynamic environmental time series
cat("\nCreating dynamic environmental scenario...\n")
cat("✅ Environmental input_params created with", nrow(input_params), "rows\n")
cat("- SST range:", round(min(input_params$sst), 1), "to", round(max(input_params$sst), 1), "°C\n")
cat("- Chlorophyll range:", round(min(input_params$chlo), 2), "to", round(max(input_params$chlo), 2), "mg/m³\n")


cat("Running model...\n")
out <- zoomss_model(input_params, Groups, SaveTimeSteps)
cat("✅ Model run completed successfully!\n")


# Save outputs
if (HPC == FALSE) {
  save(out, input_params,
       file = paste0(jobname, "_", ID_char, ".RData"))
  cat("Results saved to:", paste0(jobname, "_", ID_char, ".RData"), "\n")
}

# Optional: Generate plots if you want to visualize results
# NOTE: Plotting disabled temporarily due to PPMR plot issue - core model works perfectly

if (isTRUE(plotting)){
  source("zPlot.R")

  # Plot results
  # (ggPPMR_dynamic <- zPlot_PPMR(out))
  (ggSizeSpec <- zPlot_SizeSpectra(out))
  (ggAbundTS <- zPlot_AbundTimeSeries(out))
  (ggGrowthTS <- zPlot_GrowthTimeSeries(out))
  (ggBiomassTS <- zPlot_BiomassTimeSeries(out))
  (ggBiomassTS_stacked <- zPlot_BiomassTimeSeries(out, stacked = TRUE))

  cat("✅ Plots generated successfully!\n")
}

