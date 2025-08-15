#' Get Default ZooMSS Functional Groups
#'
#' @title Load default or custom functional groups for ZooMSS model
#' @description Provides access to the default ZooMSS functional groups or loads custom
#'   groups from a file. This function is the primary way to obtain Groups data for
#'   ZooMSS model runs.
#' @details This function provides flexible access to functional groups data:
#'   - **Default groups**: Returns the standard ZooMSS functional groups (9 groups)
#'   - **Custom file**: Loads and validates groups from a user-provided CSV file
#'   - **Template creation**: Exports default groups to a file for user modification
#'   
#'   The default groups include: Flagellates, Ciliates, Larvaceans, OmniCopepods,
#'   CarnCopepods, Euphausiids, Chaetognaths, Salps, and Jellyfish.
#'   
#'   All groups data is validated to ensure it contains required columns and
#'   reasonable parameter values for successful model runs.
#'
#' @param source Character string specifying data source. Options:
#'   \itemize{
#'     \item "default": Use built-in ZooMSS functional groups
#'     \item "file": Load groups from a CSV file
#'     \item "template": Export default groups to a file for modification
#'   }
#' @param file Path to CSV file when source="file" or source="template"
#'
#' @return Data frame containing functional groups with required columns:
#'   Species, Type, W0, Wmax, and other biological parameters
#' @export
#'
#' @examples
#' \dontrun{
#' # Use default groups
#' Groups <- zGetGroups()
#' 
#' # Create a template file for modification
#' zGetGroups(source = "template", file = "my_groups.csv")
#' 
#' # Load custom groups from file
#' custom_groups <- zGetGroups(source = "file", file = "my_groups.csv")
#' 
#' # Modify default groups programmatically
#' Groups <- zGetGroups()
#' Groups$W0[Groups$Species == "Flagellates"] <- -12.5  # Modify minimum size
#' }
#'
zGetGroups <- function(source = c("default", "file", "template"), file = NULL) {
  
  source <- match.arg(source)
  
  switch(source,
    "default" = {
      # Load from package data - this would reference the built-in TestGroups
      groups <- zLoadDefaultGroups()
      message("Using default ZooMSS functional groups (9 groups)")
      return(groups)
    },
    
    "file" = {
      if (is.null(file)) {
        stop("file path must be specified when source='file'")
      }
      if (!file.exists(file)) {
        stop("File not found: ", file)
      }
      
      groups <- utils::read.csv(file, stringsAsFactors = FALSE)
      message("Loaded ", nrow(groups), " functional groups from: ", file)
      
      # Validate the loaded groups
      zValidateGroups(groups)
      return(groups)
    },
    
    "template" = {
      if (is.null(file)) {
        stop("file path must be specified when source='template'")
      }
      
      # Get default groups and write to file
      default_groups <- zLoadDefaultGroups()
      utils::write.csv(default_groups, file, row.names = FALSE)
      
      message("✅ Template functional groups written to: ", file)
      message("Edit this file to customize groups, then load with:")
      message("Groups <- zGetGroups(source='file', file='", file, "')")
      
      return(default_groups)
    }
  )
}

#' Load Default Functional Groups Data
#'
#' @title Internal function to load default ZooMSS groups
#' @description Loads the default functional groups from the package data or CSV file.
#'   This is an internal function used by zGetGroups().
#' @details This function handles the actual loading of default groups data,
#'   whether from package data (if available) or from the CSV file in data-raw.
#'
#' @return Data frame with default functional groups
#' @keywords internal
#'
zLoadDefaultGroups <- function() {
  
  # First try to load from package data
  tryCatch({
    # This will work when the package is properly installed
    data("TestGroups", package = "zoomss", envir = environment())
    if (exists("TestGroups")) {
      message("ℹ Loaded default functional groups from package data")
      return(TestGroups)
    }
  }, error = function(e) {
    # Package data not available, try other locations
  })
  
  # Try loading from inst/extdata (for installed packages)
  package_file <- system.file("extdata", "TestGroups.csv", package = "zoomss")
  
  if (package_file != "") {
    groups <- readr::read_csv(package_file, show_col_types = FALSE)
    message("ℹ Loaded default functional groups from package extdata")
    return(groups)
  }
  
  # Try loading from data-raw (for development)
  if (file.exists("data-raw/TestGroups.csv")) {
    groups <- readr::read_csv("data-raw/TestGroups.csv", show_col_types = FALSE)
    message("ℹ Loaded default functional groups from data-raw/TestGroups.csv")
    return(groups)
  }
  
  # If we get here, no default groups were found
  stop("❌ Default groups file not found. Please ensure TestGroups data is available.",
       "\n   Try installing the package or ensure TestGroups.csv exists in data-raw/")
}

#' Validate Functional Groups Data
#'
#' @title Validate ZooMSS functional groups data structure and values
#' @description Performs comprehensive validation of functional groups data to ensure
#'   it meets ZooMSS model requirements.
#' @details This function validates:
#'   - Required column names are present
#'   - Data types are correct
#'   - Parameter values are within reasonable ranges
#'   - No missing values in critical columns
#'   - Size ranges are logical (W0 < Wmax)
#'
#' @param groups Data frame containing functional groups data
#'
#' @return TRUE if validation passes (invisibly), otherwise throws an error
#' @export
#'
#' @examples
#' \dontrun{
#' Groups <- zGetGroups()
#' zValidateGroups(Groups)  # Should pass
#' 
#' # This would fail validation:
#' bad_groups <- Groups
#' bad_groups$W0 <- NULL
#' zValidateGroups(bad_groups)  # Error: missing required column
#' }
#'
zValidateGroups <- function(groups) {
  
  # Load assertthat for validation
  if (!requireNamespace("assertthat", quietly = TRUE)) {
    stop("assertthat package required for groups validation")
  }
  
  # Check that groups is a data frame
  assertthat::assert_that(is.data.frame(groups), 
                         msg = "Groups must be a data frame")
  
  # Check required columns exist (based on actual TestGroups.csv structure)
  required_cols <- c("Species", "Type", "FeedType", "Prop", "W0", "Wmax", "Wmat",
                     "SearchCoef", "SearchExp", "PPMRscale", "PPMR", "FeedWidth",
                     "GrossGEscale", "Carbon", "Repro", "Fmort", "Fmort_W0", 
                     "Fmort_Wmax", "PlotColour")
  
  missing_cols <- setdiff(required_cols, names(groups))
  assertthat::assert_that(length(missing_cols) == 0,
                         msg = paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  
  # Check data types and ranges
  assertthat::assert_that(is.character(groups$Species) || is.factor(groups$Species),
                         msg = "Species column must be character or factor")
  
  assertthat::assert_that(all(!is.na(groups$Species)),
                         msg = "Species names cannot be NA")
  
  assertthat::assert_that(all(!duplicated(groups$Species)),
                         msg = "Species names must be unique")
  
  # Check size parameters
  assertthat::assert_that(is.numeric(groups$W0),
                         msg = "W0 (minimum weight) must be numeric")
  
  assertthat::assert_that(is.numeric(groups$Wmax),
                         msg = "Wmax (maximum weight) must be numeric")
  
  assertthat::assert_that(all(groups$W0 < groups$Wmax),
                         msg = "W0 must be less than Wmax for all groups")
  
  # Check reasonable size ranges (log10 weights)
  assertthat::assert_that(all(groups$W0 >= -15 & groups$W0 <= 5),
                         msg = "W0 values should be between -15 and 5 (log10 grams)")
  
  assertthat::assert_that(all(groups$Wmax >= -10 & groups$Wmax <= 10),
                         msg = "Wmax values should be between -10 and 10 (log10 grams)")
  
  # Check maturation weight
  assertthat::assert_that(all(groups$Wmat >= groups$W0 & groups$Wmat <= groups$Wmax),
                         msg = "Wmat must be between W0 and Wmax")
  
  # Check Type values
  valid_types <- c("Zooplankton", "Fish")
  assertthat::assert_that(all(groups$Type %in% valid_types),
                         msg = paste("Type must be one of:", paste(valid_types, collapse = ", ")))
  
  # Check FeedType values  
  valid_feedtypes <- c("Carnivore", "Omnivore", "FilterFeeder", "Heterotroph")
  assertthat::assert_that(all(groups$FeedType %in% valid_feedtypes),
                         msg = paste("FeedType must be one of:", paste(valid_feedtypes, collapse = ", ")))
  
  # Check biological parameters are positive where required
  assertthat::assert_that(all(groups$SearchCoef > 0),
                         msg = "SearchCoef must be positive")
  
  assertthat::assert_that(all(groups$SearchExp > 0),
                         msg = "SearchExp must be positive")
  
  assertthat::assert_that(all(groups$FeedWidth > 0),
                         msg = "FeedWidth must be positive")
  
  assertthat::assert_that(all(groups$GrossGEscale > 0),
                         msg = "GrossGEscale must be positive")
  
  assertthat::assert_that(all(groups$Carbon > 0 & groups$Carbon <= 1),
                         msg = "Carbon content must be between 0 and 1")
  
  # Check fishing mortality is non-negative
  assertthat::assert_that(all(groups$Fmort >= 0),
                         msg = "Fmort (fishing mortality) must be non-negative")
  
  message("✅ Functional groups validation passed")
  return(invisible(TRUE))
}
