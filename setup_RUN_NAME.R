## An extension of the model published in Heneghan et al., (2016):
## Models multiple zooplankton functional groups, and three fish groups
## This code is to run the model across multiple cores
##
##
## Updated Friday 31st January 2020
## Updated Friday 13th March 2020
## Updated Tuesday 17th March 2020

source("fZooMSS_Model.R") #source the model code

save_all <- 1 # Should we save all time steps
HPC <- 0 # Is this being run on a HPC or will we choose the cell

enviro_data <- readRDS("envirofull_20200312.RDS") # Load environmental data

enviro_data$tmaxx <- 500 # Set length of simulation (years)

#job specifics
Groups <- read.csv("TestGroups.csv") # Load in functional group information

# Should the filter feeder groups (salps/larv) have fixed PPMR? 
# TRUE(default - for most runs) = YES, 
# FALSE = NO (Use for PPMR=100 and PPMR=1000 runs).
fixed_filterPPMR <- TRUE

jobname <- 'DATE_JOBNAME' #job name used on queue: Recommend: YYYYMMDD_AbbrevExperimentName.

if (HPC == 1){
  ID <- as.integer(Sys.getenv('PBS_ARRAY_INDEX'))} else { # Get the array run number
    ID <- 100 # Test Run
    } 

ID_char <- sprintf("%04d",ID) # Set the ID as a 4 digit character so it will sort properly

system.time(
  out <- fZooMSS_Model(enviro_data$sst[ID], enviro_data$chlo[ID], 
                enviro_data$a[ID], enviro_data$b[ID], 
                enviro_data$phyto_max[ID], enviro_data$dt[ID], 
                enviro_data$tmaxx[ID], fixed_filterPPMR, save_all)
)

saveRDS(out, file = paste0("RawOutput/", jobname, "_", ID_char,".RDS"))
