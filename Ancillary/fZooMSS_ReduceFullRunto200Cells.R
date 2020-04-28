
fZooMSS_ReduceFullRunto200Cells = function(){
  library(tidyverse)

  envirofull <- read_rds("/Users/jason/Dropbox/Multi-Zoo Size Spectrum Model/_LatestModel/envirofull_20200317.RDS")
  enviro200 <- read_rds("/Users/jason/Dropbox/Multi-Zoo Size Spectrum Model/_LatestModel/enviro200_20200317.RDS")

  fi <- NA
  for (f in 1:nrow(enviro200)){
    fii <- which(envirofull$chlo==enviro200$chlo[f] & envirofull$sst==enviro200$sst[f])
    if(length(fii)>1){
      print(fii)
    }
    fi[f] <- fii
  }

  return(fi)
}