## This is the development version of the model published in Heneghan et al., (2020):
## This version models multiple zooplankton functional groups, and three fish groups
##
## If you save timesteps, you can examine the model dynamics at
##   https://jaseeverett.shinyapps.io/ZooMSS_Dashboard/
##  Alternatively you can clone the ZooMSS_Dashboard from Github
##   https://github.com/MathMarEcol/ZooMSS_Dashboard
##
##
## Code written by Dr Jason Everett (UQ/UNSW/CSIRO), Dr Ryan Heneghan (QUT) and Mr Patrick Sykes (UQ)
## Last updated Friday 31st March 2023

# library(Rcpp) # Only needed if we are running with Rcpp code.
source("fZooMSS_Model.R") #source the model code
source("fZooMSS_Xtras.R")


# Setup user defined parameters -------------------------------------------

# enviro_data <- readRDS("enviroData_oneDeg_20210728.rds") # Load environmental data at 1 degree resolution
# enviro_data <- readRDS("enviroData_fiveDeg_20200317.rds") # Load environmental data at 5 degree resolution

# # You can also create your own environmental data using the below.
enviro_data <- fZooMSS_CalculatePhytoParam(data.frame(cellID = 1,
                                                      sst = 15,
                                                      chlo = 0.5))

# Add delta time (years)
enviro_data$dt <- 0.01

# Set length of simulation (years)
enviro_data$tmax <- 50

# Setup jobname
jobname <- "DATE_JOBNAME"  # This is the job name used on the HPC queue, and also to save the run: Recommend: YYYYMMDD_AbbrevExperimentName.
enviro_row <- 1 # Which row of the environmental data do you want to run if HPC=FALSE.

HPC <- FALSE # Is this being run on a HPC for all cells or will we manually choose the row of the enviro_data to be used.
SaveTimeSteps <- TRUE # Should we save all time steps. This can be very large if tmax is large


Groups <- read.csv("TestGroups.csv", stringsAsFactors = FALSE) # Load in functional group information. This can be edited directly.


# No need to change anything below here. ----------------------------------

if (HPC == TRUE){
  ID <- as.integer(Sys.getenv('PBS_ARRAY_INDEX')) # Get the array run number on HPC
  } else {
    ID <- enviro_row
  }

ID_char <- sprintf("%04d",ID) # Set the ID as a 4 digit character so it will sort properly

input_params <- enviro_data[ID,]
rm(enviro_data)

out$model$model_runtime <- system.time(
  out <- fZooMSS_Model(input_params, Groups, SaveTimeSteps)
)

# Save the output if you want
saveRDS(out, file = paste0("RawOutput/", jobname, "_", ID_char,".RDS"))


# Plotting ----------------------------------------------------------------

source("fZooMSS_Plot.R")

(ggPPMR <- fZooMSS_Plot_PPMR(out))
(ggSizeSpec <- fZooMSS_Plot_SizeSpectra(out))

## If you have saved the timesteps you can plot the timeseries
(ggAbundTS <- fZooMSS_Plot_AbundTimeSeries(out))
(ggGrowthTS <- fZooMSS_Plot_GrowthTimeSeries(out))
(ggPredTS <- fZooMSS_Plot_PredTimeSeries(out))

