fZooMSS_Model <- function(input_params, Groups, save_all){

  source("fZooMSS_Params.R")
  source("fZooMSS_Setup.R")
  source("fZooMSS_Run.R")

  ##################### RUN THE MODEL ################################################
  param <- fZooMSS_Params(Groups, input_params) # Set up parameter list
  model <- fZooMSS_Setup(param) # Set up model equation stuff
  # cat(tracemem(model$N), "\n")
  model_output <- fZooMSS_Run(model) # Run the model

  ################### OUTPUT ABUNDANCES ##############################################
  ave_abundances <- colMeans(model_output$N[(ceiling(0.5*dim(model_output$N)[1])):(dim(model_output$N)[1]),,], dims = 1)
  ave_diets <- colMeans(model_output$diet[(ceiling(0.5*dim(model_output$diet)[1])):(dim(model_output$diet)[1]),,], dims = 1)
  ave_growth <- colMeans(model_output$gg[(ceiling(0.5*dim(model_output$gg)[1])):(dim(model_output$gg)[1]),,], dims = 1)
  ave_mort <- colMeans(model_output$Z[(ceiling(0.5*dim(model_output$Z)[1])):(dim(model_output$Z)[1]),,], dims = 1)

  if (SaveTimeSteps == TRUE){
    model$Abundance <- rowSums(model$N) ## Save Total Abundance
    # model$Biomass <- rowSums(model$N*matrix(model$w, nrow = model$param$ngrps, ncol = model$ngrid, byrow = TRUE)) ## Save biomass

  results <- list("abundances" = ave_abundances, # Save mean abundance
                 "diets" = ave_diets,  # Save mean diets
                 "growth" = ave_growth,  # Save mean growth
                 "mortality" = ave_mort, # Save mean predation
                 "model" = model_output) # Save whole model
  }
  if (SaveTimeSteps == FALSE){
    reduce_output <- list(param = model_output$param) # Create a new list so we can have identical structure of model$param
    results <- list("abundances" = ave_abundances, # Save mean abundance
                   "diets" = ave_diets,  # Save mean diets
                   "growth" = ave_growth,  # Save mean growth
                   "mortality" = ave_mort, # Save mean predation
                   "model" = reduce_output) # Save parameters only
  }

  return(results)
}
