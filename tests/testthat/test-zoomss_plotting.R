# Test suite for ZooMSS plotting functions
# Tests cover: plotPPMR, plotSizeSpectra, plotTimeSeries, plotEnvironment
# Following tidyverse and R packages testing standards

# Setup test data --------------------------------------------------------

env_data <- createEnviroData(
  n_years = 20,
  dt = 0.1,
  seasonal = FALSE,
  base_sst = 20,
  base_chl = 1.0
)

# Get groups data
Groups <- getGroups()

# Create mock model results for plotting tests
mdl <- zoomss_model(input_params = env_data, Groups = Groups, isave = 1)


# Tests fosr plotPPMR() ---------------------------------------------------

test_that("plotPPMR returns ggplot object", {
  # Test that function returns a ggplot object
  expect_no_error(plot_obj <- plotPPMR(mdl, idx = 1))
  expect_s3_class(plot_obj, "ggplot")
})

# test_that("plotPPMR handles model with all zooplankton groups", {
#
#   mdl2 <- mdl
#
#   # Remove fish groups (those with NA PPMRscale)
#   fish_rows <- which(is.na(mdl2$param$Groups$PPMRscale))
#   if(length(fish_rows) > 0) {
#     mdl2$param$Groups <- mdl2$param$Groups[-fish_rows, ]
#     mdl2$abundance <- mdl2$abundance[, -fish_rows, ]
#     # Also need to update diet array if it exists
#     if(!is.null(mdl2$diet)) {
#       # Diet array has groups as predators, so remove those rows
#       mdl2$diet <- mdl2$diet[, -fish_rows, ]
#     }
#   }
#
#   # This test may fail with bandwidth error if too few groups remain
#   # In that case, we expect either success or a specific bandwidth error
#   result <- tryCatch({
#     plot_obj <- plotPPMR(mdl2, idx = 1)
#     expect_s3_class(plot_obj, "ggplot")
#     TRUE
#   }, error = function(e) {
#     # Accept bandwidth errors as they're expected with few data points
#     expect_true(grepl("bandwidth|points", e$message, ignore.case = TRUE))
#     TRUE
#   })
#
#   expect_true(result)
#   rm(mdl2)
# })

test_that("plotPPMR handles empty abundance data gracefully", {
  mdl2 <- mdl

  # Set all abundances to zero
  mdl2$abundance[] <- 0

  expect_warning(plot_obj <- plotPPMR(mdl2, idx = 1), "Non-finite or zero total biomass")
  expect_s3_class(plot_obj, "ggplot")

  rm(mdl2)
})

# Tests for plotSizeSpectra() --------------------------------------------

test_that("plotSizeSpectra returns ggplot object", {
  expect_no_error(plot_obj <- plotSizeSpectra(mdl, n_years = 5))
  expect_s3_class(plot_obj, "ggplot")
})

test_that("plotSizeSpectra filters zero abundances", {
  mdl2 <- mdl

  # Set some abundances to zero (time, groups, size)
  mdl2$abundance[1, 1, 1:50] <- 0

  expect_no_error(plot_obj <- plotSizeSpectra(mdl2, n_years = 5))
  expect_s3_class(plot_obj, "ggplot")

  rm(mdl2)
})

test_that("plotSizeSpectra handles all zero abundances", {
  mdl2 <- mdl

  # Set all abundances to zero except one small value (time, groups, size)
  mdl2$abundance[] <- 0
  mdl2$abundance[1, 1, 1] <- 0.001

  expect_no_error(plot_obj <- plotSizeSpectra(mdl2, n_years = 5))
  expect_s3_class(plot_obj, "ggplot")
  rm(mdl2)
})

# Tests for plotTimeSeries() ---------------------------------------------

test_that("plotTimeSeries returns ggplot object for abundance", {
  expect_no_error(plot_obj <- plotTimeSeries(mdl, by = "abundance"))
  expect_s3_class(plot_obj, "ggplot")
})

test_that("plotTimeSeries returns ggplot object for biomass", {
  expect_no_error(plot_obj <- plotTimeSeries(mdl, by = "biomass"))
  expect_s3_class(plot_obj, "ggplot")
})

test_that("plotTimeSeries returns ggplot object for mortality", {
  expect_no_error(plot_obj <- plotTimeSeries(mdl, by = "mortality"))
  expect_s3_class(plot_obj, "ggplot")
})

test_that("plotTimeSeries returns ggplot object for growth", {
  expect_no_error(plot_obj <- plotTimeSeries(mdl, by = "growth"))
  expect_s3_class(plot_obj, "ggplot")
})

test_that("plotTimeSeries works with stacked biomass", {
  expect_no_error(plot_obj <- plotTimeSeries(mdl, by = "biomass", type = "stack"))
  expect_s3_class(plot_obj, "ggplot")
})

test_that("plotTimeSeries works with proportional biomass", {
  expect_no_error(plot_obj <- plotTimeSeries(mdl, by = "biomass", type = "fill"))
  expect_s3_class(plot_obj, "ggplot")
})

test_that("plotTimeSeries works with species filtering", {
  expect_no_error(plot_obj <- plotTimeSeries(mdl, by = "abundance", species = c("Flagellates")))
  expect_s3_class(plot_obj, "ggplot")
})

test_that("plotTimeSeries handles invalid by parameter", {
  expect_error(plotTimeSeries(mdl, by = "invalid"), "must be one of")
})

test_that("plotTimeSeries handles missing time series data", {
  mdl2 <- mdl
  mdl2$abundance <- NULL  # Remove time series data
  expect_error(plotTimeSeries(mdl2, by = "abundance"), "Time series data not available")
})

test_that("plotTimeSeries handles invalid species names", {
  expect_error(plotTimeSeries(mdl, by = "abundance", species = c("InvalidSpecies")),
               "All elements of 'species' must be in")
})

test_that("plotTimeSeries handles all invalid species names", {
  expect_error(plotTimeSeries(mdl, by = "abundance", species = c("InvalidSpecies1", "InvalidSpecies2")),
               "All elements of 'species' must be in")
})

# Tests for plotEnvironment() --------------------------------------------

test_that("plotEnvironment returns ggplot object or list", {
  expect_no_error(plot_obj <- plotEnvironment(env_data))

  # Should return either a ggplot (if patchwork available) or a list
  expect_true(inherits(plot_obj, "ggplot") || is.list(plot_obj))
})

test_that("plotEnvironment handles missing columns gracefully", {

  env_data2 <- dplyr::select(env_data, -chl)  # Remove chl column

  # Missing chl column should cause error in pivot_longer
  expect_error(plotEnvironment(env_data2))

  rm(env_data2)
})

test_that("plotEnvironment works with small datasets", {
  env_data2 <- env_data[1:5,]

  expect_no_error(plot_obj <- plotEnvironment(env_data2))
  expect_true(inherits(plot_obj, "ggplot") || is.list(plot_obj))

  rm(env_data2)
})

