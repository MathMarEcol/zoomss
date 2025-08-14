#' Merge Multiple Cell Outputs from ZooMSS Runs
#'
#' @title Combine ZooMSS outputs from multiple spatial cells
#' @description This function merges ZooMSS model outputs from multiple spatial cells 
#'   (grid cells) into single consolidated datasets for analysis and visualization.
#' @details When running ZooMSS across multiple spatial locations (e.g., different 
#'   ocean grid cells), this function combines the individual cell outputs into 
#'   unified datasets. The function:
#'   - Loads all RDS files from the "RawOutput" directory
#'   - Extracts abundances, growth rates, and diet data from each cell
#'   - Combines them into lists where each element represents one spatial cell
#'   - Saves the merged results in various formats (RDS and RData)
#'   - Preserves the model parameters from the first cell (assumed identical across cells)
#'   
#'   This is essential for spatial ZooMSS analyses where the model is run 
#'   independently for multiple locations and the results need to be combined
#'   for regional or global analysis.
#'
#' @return NULL (function saves files directly to disk)
#' @export
#'
#' @note This function assumes:
#'   - All individual cell outputs are stored as RDS files in "RawOutput/" directory
#'   - All cells used identical model parameters
#'   - Output directory "Output/" exists or can be created
#'   - The function is run from the directory containing "RawOutput/" subdirectory
#'
#' @examples
#' \dontrun{
#' # Typically run after completing multiple ZooMSS cell runs
#' fZooMSS_MergeMultipleCells()
#' 
#' # The function will create files like:
#' # "Output/model_[run_name].RDS"
#' # "Output/res_[run_name].RDS" 
#' # "Output/growth_[run_name].RDS"
#' # "Output/diets_[run_name].RDS"
#' # "Output/full_[run_name].RData"
#' }
#'
fZooMSS_MergeMultipleCells <- function(){
  run <- basename(getwd())

  library(tidyverse)
  files <- sort(list.files("RawOutput", full.names = TRUE))

  res <- list()
  diets <- list()
  growth <- list()

  for (f in 1:length(files)){
    out <- readr::read_rds(files[f])

    res[[f]] <- out$abundances
    growth[[f]] <- out$growth
    diets[[f]] <- out$diets

    if (f == 1){
      mdl <- out$model
    }

    rm(out)
  }

  saveRDS(mdl, paste0("Output/model_",run,".RDS"))
  saveRDS(res, paste0("Output/res_",run,".RDS"))
  saveRDS(growth, paste0("Output/growth_",run,".RDS"))
  saveRDS(diets, paste0("Output/diets_",run,".RDS"))

  save(list = c("res", "growth", "diets", "mdl"), file = paste0(paste0("Output/full_",run,".RData")))

  print("Consider running `tar -cf RawOutput.tar RawOutput` to reduce the number of files to be synced")

}
