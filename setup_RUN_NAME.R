## An extension of the model published in Heneghan et al., (2016):
## Models multiple zooplankton functional groups, and three fish groups
## This code is to run the model across multiple cores
##
## Jason's version for UNSWs Katana
##
## Updated Friday 31st January 2020

source("fZooMSS_Model.R") #source the model code

# Choose environmental data to use
# enviro_data <- readRDS("enviro5d.RDS") # full set of 5x5 degree grids
enviro_data <- readRDS("envirofull_20200312.RDS") # 200 square subset used for testing
enviro_data$tmaxx <- 5 # Set length of simulation (years)

#job specifics
Groups <- read.csv("TestGroups.csv") # Load in functional group information

# Should the filter feeder groups (salps/larv) have fixed PPMR? 
# TRUE(default - for most runs) = YES, 
# FALSE = NO (Use for PPMR=100 and PPMR=1000 runs).
fixed_filterPPMR <- TRUE

jobname <- 'DATE_JOB.NAME' #job name used on queue

ID <- as.integer(Sys.getenv('PBS_ARRAY_INDEX')) # Get the array run number

ID <- 1 # Test Run

ID_char <- sprintf("%04d",ID) # Set the ID as a 4 digit character so it will sort properly

system.time(
  out_3d <- ZooMSS(enviro_data$sst[ID], enviro_data$chlo[ID], 
                enviro_data$a[ID], enviro_data$b[ID], 
                enviro_data$phyto_max[ID], enviro_data$dt[ID], 
                enviro_data$tmaxx[ID], fixed_filterPPMR)
)

# saveRDS(out, file = paste0("RawOutput/", jobname, "_", ID_char,".RDS"))


# system elapsed
# 6.508   1.289   7.713 # 3D with colSums(aperm())
# 11.707   1.222  11.706 # 3D
# 21.917   3.778  25.774 # 4D


