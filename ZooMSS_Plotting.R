
library(tidyverse)
res <- readRDS("Output/res_Control_20200131.RDS")
Groups <- read_csv("TestGroups.csv")

enviro_data <- readRDS("enviro200.RDS")

### ACTUAL BIOMASS
biom_act <- function(N, groups, cut_point1, cut_point2, enviro_data, en){

  w <- 10^(seq(from = -10.7, to =  7, 0.1)) # Size bins of whole model

  zoo_groups = which(is.na(groups$prop) == FALSE)[-c(1,2)] # Which rows are zooplankton
  num_zoo = length(zoo_groups)

  abund_mat <- matrix(0, nrow = dim(enviro_data)[1], ncol = num_zoo)


  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    B_ave = N_ave*matrix(w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
    weight_cut = which(round(log10(w),2) >= cut_point1 & round(log10(w),2) <= cut_point2)
    zoo_abunds = rowSums(B_ave[zoo_groups, weight_cut])
    zoo_props = zoo_abunds
    abund_mat[i,] <- zoo_props
  }

  if(en == "chlo"){
    par(mfrow = c(4,2))
    for( i in 1:num_zoo){
      plot(log10(enviro_data$chlo), abund_mat[,i], main = groups$species[i+2], xlab = "log10(Chlo)", ylab = expression(paste("Biomass (g m"^-3, ")")))
    }
  }

  if(en == "sst"){
    par(mfrow = c(4,2))
    for( i in 1:num_zoo){
      plot((enviro_data$sst), abund_mat[,i], main = groups$species[i+2], xlab = "SST", ylab = expression(paste("Biomass (# m"^-3, ")")))
    }
  }
  par(mfrow = c(1,1))


}

####### PROPORTIONS OF BIOMASS
biom_props <- function(N, groups, cut_point1, cut_point2, enviro_data, en){

  w <- 10^(seq(from = -10.7, to =  7, 0.1)) # Size bins of whole model
  zoo_groups = which(is.na(groups$prop) == FALSE)[-c(1,2)] # Which rows are zooplankton
  num_zoo = length(zoo_groups)

  abund_mat <- matrix(0, nrow = dim(enviro_data)[1], ncol = num_zoo)


  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    B_ave = N_ave*matrix(w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
    weight_cut = which(round(log10(w),2) >= cut_point1 & round(log10(w),2) <= cut_point2)
    zoo_abunds = rowSums(B_ave[zoo_groups, weight_cut])
    zoo_props = zoo_abunds/sum(zoo_abunds)
    abund_mat[i,] <- zoo_props
  }

  if(en == "chlo"){
    par(mfrow = c(4,2))
    for( i in 1:num_zoo){
      plot(log10(enviro_data$chlo), abund_mat[,i], main = groups$species[i+2], xlab = "log10(Chlo)", ylab = "Biom Prop")
    }
  }

  if(en == "sst"){
    par(mfrow = c(4,2))
    for( i in 1:num_zoo){
      plot((enviro_data$sst), abund_mat[,i], main = groups$species[i+2], xlab = "SST", ylab = "Biom Prop")
    }
  }
  par(mfrow = c(1,1))


}


filename = "BiomProps_Chlo.pdf"
pdf(filename, width = 5.8, height = 8.3)
biom_props(res, groups = Groups, -10.7, 3, enviro_data, en = "chlo")
dev.off()

filename = "BiomProps_SST.pdf"
pdf(filename, width = 5.8, height = 8.3)
biom_props(res, groups = Groups, -10.7, 3, enviro_data, en = "sst")
dev.off()

filename = "BiomAct_Chlo.pdf"
pdf(filename, width = 5.8, height = 8.3)
biom_act(res, groups = Groups, -10.7, 3, enviro_data, en = "chlo")
dev.off()

filename = "BiomAct_SST.pdf"
pdf(filename, width = 5.8, height = 8.3)
biom_act(res, groups = Groups, -10.7, 3, enviro_data, en = "sst")
dev.off()
