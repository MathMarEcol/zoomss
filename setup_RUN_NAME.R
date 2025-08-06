## Dynamic ZooMSS Setup Script
## This script demonstrates how to run ZooMSS with time-varying environmental data
## Based on the original setup by Dr Jason Everett, Dr Ryan Heneghan, and Mr Patrick Sykes
## Dynamic environmental forcing only - August 2025

# Load all required functions
source("fZooMSS_Model.R") #source the dynamic model code
source("fZooMSS_Xtras.R") # For fZooMSS_CalculatePhytoParam and other utilities
source("fZooMSS_Environmental_Utils.R") # For environmental time series functions

# Setup user defined parameters -------------------------------------------
dt <- 0.01 # years^-1
tmax <- 250 # years
cellID <- 1
isave <- 100 # Match static model's save frequency
base_sst <- 15
base_chlo <- 0.01

# Create expanded input_params dataframe - one row per timestep
input_params <- data.frame(
  cellID = cellID,
  dt = dt,
  tmax = tmax,
  isave = isave,
  time_step = seq(1, 1/dt * tmax, 1),
  sst = base_sst,
  chlo = base_chlo
)

# input_params$chlo <- 0.1 + 0.0002 * input_params$time_step # Linear increasing Chlorophyll
# input_params$chlo <- 0.5 + 0.49 * sin(2 * pi * (input_params$time_step * dt)) # Chlorophyll cycles between 0.01-0.99 mg/m³

# Ensure chlorophyll stays positive
input_params$chlo <- pmax(input_params$chlo, 0.01)
cat("Environmental data created with", nrow(input_params), "timesteps covering",
    tmax, "years\n")

fZooMSS_PlotEnvironment(input_params)


# Setup jobname
jobname <- "20250806_chl0_01"  # This is the job name to save the run
enviro_row <- 1 # Which row of the environmental data to use
HPC <- FALSE # Is this being run on a HPC
SaveTimeSteps <- TRUE # Should we save all time steps
Groups <- read.csv("TestGroups.csv", stringsAsFactors = FALSE) # Load functional group information

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


cat("Running dynamic model...\n")
out <- fZooMSS_Model(input_params, Groups, SaveTimeSteps)
cat("✅ Dynamic model run completed successfully!\n")


# Save outputs
if (HPC == FALSE) {
  save(out, input_params,
       file = paste0(jobname, "_", ID_char, ".RData"))
  cat("Results saved to:", paste0(jobname, "_", ID_char, ".RData"), "\n")
}

# Optional: Generate plots if you want to visualize results
# NOTE: Plotting disabled temporarily due to PPMR plot issue - core model works perfectly

source("fZooMSS_Plot.R")

# Plot results
# (ggPPMR_dynamic <- fZooMSS_Plot_PPMR(out))
(ggSizeSpec <- fZooMSS_Plot_SizeSpectra(out))
(ggAbundTS <- fZooMSS_Plot_AbundTimeSeries(out))
(ggGrowthTS <- fZooMSS_Plot_GrowthTimeSeries(out))
(ggBiomassTS <- fZooMSS_Plot_BiomassTimeSeries(out))
(ggBiomassTS_stacked <- fZooMSS_Plot_BiomassTimeSeries(out, stacked = TRUE))

cat("✅ Plots generated successfully!\n")


