#' Run Complete ZooMSS Model Simulation
#'
#' @title Main ZooMSS model function for complete simulations
#' @description This is the main wrapper function that orchestrates a complete ZooMSS 
#'   model simulation from parameter setup through model execution to output processing.
#' @details This function coordinates the entire ZooMSS modeling workflow:
#'   1. Validates that environmental time series data is provided
#'   2. Sets up model parameters using the Groups data and input parameters
#'   3. Initializes the model structure and feeding kernels
#'   4. Runs the model forward in time with dynamic environmental forcing
#'   5. Processes outputs by averaging the final 50% of the simulation
#'   6. Returns organized results including abundances, diets, growth, and mortality
#'   
#'   The function can optionally save full time series data or just averaged results
#'   depending on the SaveTimeSteps parameter. This is the primary entry point for
#'   running ZooMSS simulations with environmental forcing.
#'
#' @param input_params Data frame containing model parameters and environmental time series.
#'   Must include columns: tmax (simulation length), dt (time step), isave (save frequency),
#'   and environmental data (time_step, sst, chlo) if running with environmental forcing.
#' @param Groups Data frame defining functional groups with their biological parameters.
#'   Must include columns defining species characteristics, size ranges, and feeding parameters.
#' @param SaveTimeSteps Logical indicating whether to save full time series (TRUE) or 
#'   just model parameters (FALSE). TRUE saves complete model output including time series
#'   of abundance and biomass.
#'
#' @return List containing:
#'   \itemize{
#'     \item abundances: Average abundances across the final 50% of simulation
#'     \item diets: Average diet compositions across functional groups
#'     \item growth: Average growth rates for each group and size class
#'     \item mortality: Average mortality rates
#'     \item model: Complete model output (if SaveTimeSteps=TRUE) or just parameters (if FALSE)
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # Load functional groups data
#' data(Groups)
#' 
#' # Create environmental time series
#' env_data <- fZooMSS_CreateSimpleTimeSeries(1000, 0.01)
#' 
#' # Set up input parameters
#' input_params <- data.frame(
#'   tmax = 10,
#'   dt = 0.01, 
#'   isave = 10,
#'   time_step = 1:1000,
#'   sst = env_data$sst,
#'   chlo = env_data$chlo
#' )
#' 
#' # Run the model
#' results <- fZooMSS_Model(input_params, Groups, SaveTimeSteps = FALSE)
#' }
#'
fZooMSS_Model <- function(input_params, Groups, SaveTimeSteps){

  input_params <- untibble(input_params)
  
  # Validate that input_params has the required environmental data
  if (nrow(input_params) > 1 && all(c("time_step", "sst", "chlo") %in% names(input_params))) {
    # Environmental time series found - proceed with model
  } else {
    stop("No environmental time series provided and input_params doesn't contain expanded time series data")
  }

  ################### RUN THE MODEL ###################
  param <- fZooMSS_Params(Groups, input_params) # Set up parameter list
  model <- fZooMSS_Setup(param) # Set up model equation stuff
  model_output <- fZooMSS_Run(model) # Run the model

  ################### AVERAGE THE LAST 50 % OF THE MODEL RUN ###################

  ave_abundances <- fZooMSS_AveOutput(model_output$N)
  ave_diets <- fZooMSS_AveOutput(model_output$diet)
  ave_growth <- fZooMSS_AveOutput(model_output$gg)
  ave_mort <- fZooMSS_AveOutput(model_output$Z)

  if (SaveTimeSteps == TRUE){
    model_output$Abundance <- rowSums(model$N, dims = 2) ## Save Total Abundance
    model_output$Biomass <- colSums(aperm(sweep(model$N, 3, model$param$w, "*"), c(3,1,2)))

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
