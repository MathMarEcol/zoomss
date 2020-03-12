ZooMSS <- function(sst, chlo, a, b, phyto_max, dt, tmaxx, fixed_filterPPMR){

  source("fZooMSS_Params.R")
  source("fZooMSS_Setup.R")
  source("fZooMSS_Project.R")

  enviro_vector <- data.frame("sst" = sst, "chlo" = chlo, "a" = a, "b" = b,
                              "phyto_max" = phyto_max, "dt" = dt)

  ##################### RUN THE MODEL ################################################
  param <- Params(Groups, enviro_vector, tmax = tmaxx, f_mort = 0) # Set up parameter list
  model <- Setup(param) # Set up model equation stuff
  modelss <- Project(model, fish_on = TRUE) # Run the model
      
  saveRDS(model, "Output/ModelParameters.RDS")
  
  ################### OUTPUT ABUNDANCES ##############################################
  ave_abundances = colMeans(modelss$N[(ceiling(0.5*dim(modelss$N)[1])):(dim(modelss$N)[1]),,], dim = 1)
  ave_diets = colMeans(modelss$diet[(ceiling(0.5*dim(modelss$diet)[1])):(dim(modelss$diet)[1]),,], dim = 1)
  ave_growth = colMeans(modelss$gg[(ceiling(0.5*dim(modelss$gg)[1])):(dim(modelss$gg)[1]),,], dim = 1)
  ave_pred = colMeans(modelss$Z[(ceiling(0.5*dim(modelss$Z)[1])):(dim(modelss$Z)[1]),,], dim = 1)
  
  results = list("abundances" = ave_abundances, # Save mean abundance
                 "diets" = ave_diets,  # Save mean diets
                 "growth" = ave_growth,  # Save mean growth
                 "predation" = ave_pred, # Save mean predation
                 "model" = modelss) # Save whole model
  return(results)
}
