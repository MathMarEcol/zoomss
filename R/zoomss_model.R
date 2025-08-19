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
#'   Must include columns: time (time vector in years), sst (sea surface temperature), 
#'   chl (chlorophyll), and isave (save frequency). Can optionally include cellID for 
#'   spatial data. The time step (dt) and maximum time (tmax) are automatically calculated 
#'   from the time vector. Can be created using zCreateInputs().
#' @param Groups Data frame defining functional groups with their biological parameters.
#'   Must include columns defining species characteristics, size ranges, and feeding parameters.
#'   If NULL, uses default ZooMSS functional groups. Can be obtained/customized using
#'   zGetGroups().
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
#' # Basic usage with default groups
#' env_data <- zCreateSimpleTimeSeries(1000, 0.01)
#' input_params <- zCreateInputs(env_data$time, env_data$sst, env_data$chl, isave = 50)
#' results <- zoomss_model(input_params, SaveTimeSteps = FALSE)
#'
#' # Using custom groups
#' Groups <- zGetGroups()  # Get default groups
#' Groups$W0[1] <- -12.5          # Modify a parameter
#' results <- zoomss_model(input_params, Groups, SaveTimeSteps = FALSE)
#'
#' # Loading groups from file
#' custom_groups <- zGetGroups(source = "file", file = "my_groups.csv")
#' results <- zoomss_model(input_params, custom_groups, SaveTimeSteps = TRUE)
#' }
#'
zoomss_model <- function(input_params, Groups = NULL, SaveTimeSteps = TRUE){

  # Handle default Groups parameter
  if (is.null(Groups)) {
    Groups <- zGetGroups(source = "default")
    message("Using default ZooMSS functional groups. Use zGetGroups() to customize.")
  } else {
    # Validate user-provided Groups
    zValidateGroups(Groups)
  }

  input_params <- untibble(input_params)

  # Validate that input_params has the required environmental data
  if (nrow(input_params) > 1 && all(c("time_step", "sst", "chl") %in% names(input_params))) {
    # Environmental time series found - proceed with model
  } else {
    stop("No environmental time series provided and input_params doesn't contain expanded time series data")
  }

  ################### RUN THE MODEL ###################
  param <- zoomss_params(Groups, input_params) # Set up parameter list
  model <- zoomss_setup(param) # Set up model equation stuff
  model_output <- zoomss_run(model) # Run the model

  ################### AVERAGE THE LAST 50 % OF THE MODEL RUN ###################

  # TODO Remove averaging of the last 50 %.
  ave_abundances <- zAveOutput(model_output$N)
  ave_diets <- zAveOutput(model_output$diet)
  ave_growth <- zAveOutput(model_output$gg)
  ave_mort <- zAveOutput(model_output$Z)

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
    # Create a new list so we can have identical structure of model$param
    reduce_output <- list(param = model_output$param)
    results <- list("abundances" = ave_abundances, # Save mean abundance
                   "diets" = ave_diets,  # Save mean diets
                   "growth" = ave_growth,  # Save mean growth
                   "mortality" = ave_mort, # Save mean predation
                   "model" = reduce_output) # Save parameters only
  }

  return(results)
}
