## An extension of the model published in Heneghan et al., (2016):
## Models multiple zooplankton functional groups, and three fish groups
## This code is to run the model across multiple cores
##
## Updated Friday 31st January 2020
## Updated Friday 13th March 2020
## Updated Tuesday 17th March 2020

source("fZooMSS_Model.R") #source the model code

enviro_data <- readRDS("envirofull_20200317.RDS") # Load environmental data
enviro_data$tmax <- 1000 # Set length of simulation (years)
Groups <- read.csv("TestGroups_20200330.csv") # Load in functional group information

jobname <- 'DATE_JOBNAME' #job name used on queue: Recommend: YYYYMMDD_AbbrevExperimentName.
enviro_row <- 10 # Which row of the environmental data do you want to run if HPC=FALSE

HPC <- FALSE # Is this being run on a HPC or will we choose the row
SaveTimeSteps <- TRUE # Should we save all time steps

### No need to change anything below here.
if (HPC == TRUE){
  ID <- as.integer(Sys.getenv('PBS_ARRAY_INDEX')) # Get the array run number on HPC
  } else { 
    ID <- enviro_row
  } 
ID_char <- sprintf("%04d",ID) # Set the ID as a 4 digit character so it will sort properly

input_params <- enviro_data[ID,]

out <- fZooMSS_Model(input_params, SaveTimeSteps)

saveRDS(out, file = paste0("RawOutput/", jobname, "_", ID_char,".RDS"))
