## ZooMSS 100-year simulation with HIGH fishing mortality for fish groups
## This script runs a ZooMSS simulation for 100 years using default values
## except with HIGH fishing mortality for the three fish groups
##
## Code written by Dr Jason Everett (UQ/UNSW/CSIRO), Dr Ryan Heneghan (QUT) and Mr Patrick Sykes (UQ)
## Modified for 100-year high fishing simulation

# library(Rcpp) # Only needed if we are running with Rcpp code.
source("fZooMSS_Model.R") #source the model code
source("fZooMSS_Xtras.R")


# Setup user defined parameters -------------------------------------------

# Create environmental data using default values
# Using typical oceanic conditions: SST = 15°C, Chlorophyll = 0.5 mg/m³
enviro_data <- fZooMSS_CalculatePhytoParam(data.frame(cellID = 1,
                                                      sst = 15,
                                                      chlo = 0.5))

# Add delta time (years) - using smaller timestep for stability
enviro_data$dt <- 0.01

# Set length of simulation to 100 years
enviro_data$tmax <- 100

# Setup jobname
jobname <- "20250806_ZooMSS_100yr_HighFishing"  # Job name for this 100-year high fishing simulation
enviro_row <- 1 # Which row of the environmental data to use

HPC <- FALSE # Running locally, not on HPC
SaveTimeSteps <- TRUE # Save all time steps for analysis

# Load functional group information with HIGH fishing mortality for fish groups
Groups <- read.csv("TestGroups_WithFishing.csv", stringsAsFactors = FALSE)

# Print fishing mortality values for fish groups to confirm setup
cat("HIGH fishing mortality rates applied:\n")
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

cat("Starting ZooMSS 100-year simulation with HIGH fishing mortality...\n")
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

# Simple plotting without the problematic functions
cat("Generating basic summary...\n")

# Calculate some basic metrics at the final time step
final_abundances <- out$model$N[dim(out$model$N)[1],,]
total_biomass_by_group <- apply(final_abundances * rep(out$model$param$w, each = 12), 1, sum)

cat("\n=== FINAL STATE SUMMARY ===\n")
cat("Total biomass by functional group at end of simulation:\n")
for(i in 1:length(Groups$Species)){
  cat(sprintf("%s: %0.2e g/m³\n", Groups$Species[i], total_biomass_by_group[i]))
}

# Print summary of simulation
cat("\n=== SIMULATION SUMMARY ===\n")
cat(sprintf("Duration: %d years\n", input_params$tmax))
cat(sprintf("Time step: %0.3f years\n", input_params$dt))
cat(sprintf("Total runtime: %0.2f seconds\n", model_runtime[3]))
cat(sprintf("Output file: %s\n", output_file))
cat("HIGH fishing mortality applied to:\n")
for(i in 1:nrow(fish_groups)){
  cat(sprintf("  %s: F = %0.2f yr^-1\n", fish_groups$Species[i], fish_groups$Fmort[i]))
}

cat("\nSimulation complete! The output is compatible with the ZooMSS Dashboard.\n")
