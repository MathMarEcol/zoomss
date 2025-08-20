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
#'   This is the primary entry point for
#'   running ZooMSS simulations with environmental forcing.
#'
#' @param input_params Data frame containing model parameters and environmental time series.
#'   Must include columns: time (time vector in years), sst (sea surface temperature),
#'   and chl (chlorophyll). Can optionally include cellID for spatial data. The time step (dt)
#'   and maximum time (tmax) are automatically calculated from the time vector. Can be created using createInputParams().
#' @param Groups Data frame defining functional groups with their biological parameters.
#'   Must include columns defining species characteristics, size ranges, and feeding parameters.
#'   If NULL, uses default ZooMSS functional groups. Can be obtained/customized using
#'   getGroups().
#' @param isave Save frequency in time steps (default: 50)
#'
#' @return List containing:
#'   \itemize{
#'     \item abundances: Average abundances across the final 50% of simulation
#'     \item diets: Average diet compositions across functional groups
#'     \item growth: Average growth rates for each group and size class
#'     \item mortality: Average mortality rates
#'     \item model: Complete model output including time series data and parameters
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with default groups
#' env_data <- createEnviroData(10, 0.01)
#' input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)
#' results <- zoomss_model(input_params, isave = 50)
#'
#' # Using custom groups
#' Groups <- getGroups()  # Get default groups
#' Groups$W0[1] <- -12.5          # Modify a parameter
#' results <- zoomss_model(input_params, Groups, isave = 100)
#'
#' # Loading groups from file
#' custom_groups <- getGroups(source = "file", file = "my_groups.csv")
#' results <- zoomss_model(input_params, custom_groups)
#' }
#'
zoomss_model <- function(input_params, Groups = NULL, isave = 10){

  # Handle default Groups parameter
  if (is.null(Groups)) {
    Groups <- getGroups(source = "default")
  } else {
    # Validate user-provided Groups
    validateGroups(Groups)
  }

  input_params <- untibble(input_params)

  # Validate that input_params has the required environmental data
  if (nrow(input_params) > 1 && all(c("time", "sst", "chl") %in% names(input_params))) {
    # Environmental time series found - proceed with model
  } else {
    stop("No environmental time series provided and input_params doesn't contain expanded time series data")
  }

  ################### RUN THE MODEL ###################
  param <- zoomss_params(Groups, input_params, isave) # Set up parameter list
  model <- zoomss_setup(param) # Set up model equation stuff
  model_output <- zoomss_run(model) # Run the model

  # model_output$Abundance <- rowSums(model$N, dims = 2) ## Save Total Abundance
  # model_output$Biomass <- colSums(aperm(sweep(model$N, 3, model$param$w, "*"), c(3,1,2)))


  return(model_output)
}
