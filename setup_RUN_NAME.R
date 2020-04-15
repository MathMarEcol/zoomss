## An extension of the model published in Heneghan et al., (2016):
## Models multiple zooplankton functional groups, and three fish groups
## This code is to run the model across multiple cores
##
## Updated Friday 31st January 2020
## Updated Friday 13th March 2020
## Updated Tuesday 17th March 2020
## Updated Tuesday 31st March 2020
## Updated Monday 13th April 2020

source("fZooMSS_Model.R") #source the model code

enviro_data <- readRDS("envirofull_20200317.RDS") # Load environmental data
enviro_data$tmax <- 25 # Set length of simulation (years)
Groups <- read.csv("TestGroups.csv", stringsAsFactors = FALSE) # Load in functional group information

jobname <- 'DATE_JOBNAME' #job name used on queue: Recommend: YYYYMMDD_AbbrevExperimentName.
enviro_row <- 1 # Which row of the environmental data do you want to run if HPC=FALSE

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

out$model$model_runtime <- system.time(
  out <- fZooMSS_Model(input_params, Groups, SaveTimeSteps)
)

saveRDS(out, file = paste0("RawOutput/", jobname, "_", ID_char,".RDS"))
#
# source("fZooMSS_CheckIdent.R")
# out_old <- readRDS("RawOutput/DATE_JOBNAME_0001_20200414.RDS")
# # out_old <- readRDS("RawOutput/DATE_JOBNAME_0350_20200401.RDS")
# # out_old <- readRDS("RawOutput/DATE_JOBNAME_0611_20200401.RDS")
#
# fZooMSS_CheckIdent(out, out_old)
#


# 25 years - after updates
# user  system elapsed
# 58.126  11.076  65.322
# 57.617   8.869  67.118

# 55.924   8.022  63.699 # Removed all list accessing
# 56.792  10.589  63.379

# 59.514   8.850  68.369 # Added list again
# 54.737   8.760  63.391 # Remove list accessing again
# 61.679  11.131  68.390
# 58.708   9.891  68.621

# 57.719  12.591  69.833 # 15 April

# 50.791   8.354  58.669 $ With Rcpp (but wrong answer)
