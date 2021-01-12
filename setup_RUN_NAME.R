## This is the development version of the model published in Heneghan et al., (2020):
## This version models multiple zooplankton functional groups, and three fish groups
##
## If you save timesteps, you can examine the model dynamics at
##   https://jaseeverett.shinyapps.io/ZooMSS_Dashboard/
##  Alternatively you can close the ZooMSS_Dashboard from Github
##   https://github.com/MathMarEcol/ZooMSS_Dashboard
##
## Last updated 28th April 2020
## Additional commenting 12th January 2021

# library(Rcpp) # Only needed if we are running with Rcpp code.

source("fZooMSS_Model.R") #source the model code
source("fZooMSS_CalculatePhytoParam.R")

# enviro_data <- readRDS("envirofull_20200317.RDS") # Load environmental data.

# You can also create your own environmental data using the below.
enviro_data <- fZooMSS_CalculatePhytoParam(data.frame(cellID = 1,
                                                      sst = 5,
                                                      chlo = 2,
                                                      dt = 0.01))

enviro_data$tmax <- 1000 # Set length of simulation (years)

# jobname <- "DATE_JOBNAME"  # This is the job name used on the HPC queue, and also to save the run: Recommend: YYYYMMDD_AbbrevExperimentName.
jobname <- paste0("YiwenExp_Chl_", enviro_data$chlo, "_SST_", enviro_data$sst)
enviro_row <- 1 # Which row of the environmental data do you want to run if HPC=FALSE.

HPC <- FALSE # Is this being run on a HPC for all cells or will we manually choose the row of the enviro_data to be used.
SaveTimeSteps <- TRUE # Should we save all time steps. This can be very large if tmax is large

Groups <- read.csv("TestGroups.csv", stringsAsFactors = FALSE) # Load in functional group information. This can be edited directly.

### No need to change anything below here.
if (HPC == TRUE){
  ID <- as.integer(Sys.getenv('PBS_ARRAY_INDEX')) # Get the array run number on HPC
  } else {
    ID <- enviro_row
  }
ID_char <- sprintf("%04d",ID) # Set the ID as a 4 digit character so it will sort properly

input_params <- enviro_data[ID,]

out$model$model_runtime <- system.time(
  out <- fZooMSS_Model(input_params, Groups, SaveTimeSteps)
)

saveRDS(out, file = paste0("RawOutput/", jobname, "_", ID_char,".RDS"))
