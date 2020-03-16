fZooMSS_Model <- function(sst, chlo, a, b, phyto_max, dt, tmaxx, fixed_filterPPMR, save_all){

  source("fZooMSS_Params.R")
  source("fZooMSS_Setup.R")
  source("fZooMSS_Project.R")

  enviro_vector <- data.frame("sst" = sst, "chlo" = chlo, "a" = a, "b" = b,
                              "phyto_max" = phyto_max, "dt" = dt)

  ##################### RUN THE MODEL ################################################
  param <- fZooMSS_Params(Groups, enviro_vector, tmax = tmaxx, f_mort = 0) # Set up parameter list
  model <- fZooMSS_Setup(param) # Set up model equation stuff
  model_output <- fZooMSS_Project(model, fish_on = TRUE) # Run the model
      
  saveRDS(model, "Output/ModelParameters.RDS")
  
  ################### OUTPUT ABUNDANCES ##############################################
  ave_abundances = colMeans(model_output$N[(ceiling(0.5*dim(model_output$N)[1])):(dim(model_output$N)[1]),,], dim = 1)
  ave_diets = colMeans(model_output$diet[(ceiling(0.5*dim(model_output$diet)[1])):(dim(model_output$diet)[1]),,], dim = 1)
  ave_growth = colMeans(model_output$gg[(ceiling(0.5*dim(model_output$gg)[1])):(dim(model_output$gg)[1]),,], dim = 1)
  ave_pred = colMeans(model_output$M2[(ceiling(0.5*dim(model_output$M2)[1])):(dim(model_output$M2)[1]),,], dim = 1)
  
  if (save_all == 1){
  results = list("abundances" = ave_abundances, # Save mean abundance
                 "diets" = ave_diets,  # Save mean diets
                 "growth" = ave_growth,  # Save mean growth
                 "predation" = ave_pred, # Save mean predation
                 "model" = model_output) # Save whole model
  }
  if (save_all == 0){
    results = list("abundances" = ave_abundances, # Save mean abundance
                   "diets" = ave_diets,  # Save mean diets
                   "growth" = ave_growth,  # Save mean growth
                   "predation" = ave_pred) # Save mean predation
  }
  
  return(results)
}
