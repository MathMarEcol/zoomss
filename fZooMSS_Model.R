fZooMSS_Model <- function(input_params, save_all){

  source("fZooMSS_Params.R")
  source("fZooMSS_Setup.R")
  source("fZooMSS_Project.R")

  ##################### RUN THE MODEL ################################################
  param <- fZooMSS_Params(Groups, input_params) # Set up parameter list
  model <- fZooMSS_Setup(param) # Set up model equation stuff
  model_output <- fZooMSS_Project(model) # Run the model
  
  ################### OUTPUT ABUNDANCES ##############################################
  ave_abundances = colMeans(model_output$N[(ceiling(0.5*dim(model_output$N)[1])):(dim(model_output$N)[1]),,], dim = 1)
  ave_diets = colMeans(model_output$diet[(ceiling(0.5*dim(model_output$diet)[1])):(dim(model_output$diet)[1]),,], dim = 1)
  ave_growth = colMeans(model_output$gg[(ceiling(0.5*dim(model_output$gg)[1])):(dim(model_output$gg)[1]),,], dim = 1)
  ave_pred = colMeans(model_output$M2[(ceiling(0.5*dim(model_output$M2)[1])):(dim(model_output$M2)[1]),,], dim = 1)
  
  if (SaveTimeSteps == TRUE){
  results = list("abundances" = ave_abundances, # Save mean abundance
                 "diets" = ave_diets,  # Save mean diets
                 "growth" = ave_growth,  # Save mean growth
                 "predation" = ave_pred, # Save mean predation
                 "param" = param, # Save parameters
                 "model" = model_output) # Save whole model
  }
  if (SaveTimeSteps == FALSE){
    results = list("abundances" = ave_abundances, # Save mean abundance
                   "diets" = ave_diets,  # Save mean diets
                   "growth" = ave_growth,  # Save mean growth
                   "predation" = ave_pred, # Save parameters
                   "param" = param) # Save mean predation
  }
  
  return(results)
}
