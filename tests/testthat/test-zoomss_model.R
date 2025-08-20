# Test suite for ZooMSS model functions
# Tests cover the main model workflow: zoomss_model, zoomss_params, zoomss_setup, zoomss_run
# Following tidyverse and R packages testing standards

# Setup test data --------------------------------------------------------

# Create minimal test environment data
create_test_env_data <- function(n_years = 20, dt = 0.1) {
  time_vec <- seq(0, n_years - dt, by = dt)

  list(
    time = time_vec,
    sst = 15 + 3 * sin(2 * pi * time_vec),  # Seasonal temperature variation
    chl = 2 + 1 * sin(2 * pi * time_vec + pi/2)  # Seasonal chlorophyll variation
  )
}

# Create minimal test groups data
create_test_groups <- function() {
  data.frame(
    Species = c("TestZoo1", "TestFish1"),
    Type = c("Zooplankton", "Fish"),
    FeedType = c("FilterFeeder", "Carnivore"),
    Prop = c(1.0, NA),
    W0 = c(-12.0, -6.0),      # Use exact values that match w_log10 grid
    Wmax = c(-6.0, 3.0),      # Use exact values that match w_log10 grid
    Wmat = c(-8.0, 0.0),      # Use exact values that match w_log10 grid
    SearchCoef = c(640, 640),
    SearchExp = c(0.8, 0.8),
    PPMRscale = c(1.0, 1.0),
    PPMR = c(NA, 1000),
    FeedWidth = c(2.0, 3.0),
    GrossGEscale = c(2.5, 2.5),
    Carbon = c(0.1, 0.1),
    Repro = c(0.0, 0.0),
    Fmort = c(0.0, 0.1),
    Fmort_W0 = c(-8.0, 0.0),  # Use exact values that match w_log10 grid
    Fmort_Wmax = c(-6.0, 3.0),  # Use exact values that match w_log10 grid
    PlotColour = c("blue", "red"),
    stringsAsFactors = FALSE
  )
}

# Tests for zoomss_model() -----------------------------------------------

test_that("zoomss_model runs with minimal input", {
  skip_if_not_installed("zoomss")

  # Create test data
  env_data <- create_test_env_data(n_years = 20, dt = 0.1)
  input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)
  test_groups <- create_test_groups()

  # Test basic model run
  expect_no_error({
    result <- zoomss_model(input_params, test_groups, isave = 5)
  })

  # Test with default groups
  expect_no_error({
    result <- zoomss_model(input_params, isave = 5)
  })
})

test_that("zoomss_model returns expected structure", {
  skip_if_not_installed("zoomss")

  env_data <- create_test_env_data(n_years = 20, dt = 0.1)
  input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)
  test_groups <- create_test_groups()

  result <- zoomss_model(input_params, test_groups, isave = 5)

  # Check that result is a list
  expect_type(result, "list")

  # Check for expected components
  expect_true("param" %in% names(result))
  expect_true("N" %in% names(result))
  expect_true("gg" %in% names(result))

  # Check dimensions are reasonable
  expect_true(is.array(result$N))
  expect_true(is.array(result$gg))
  expect_equal(length(dim(result$N)), 3)  # time x groups x size_classes
})

test_that("zoomss_model validates input parameters", {
  skip_if_not_installed("zoomss")

  # Test with missing environmental data
  bad_input <- data.frame(time = 1:10)
  test_groups <- create_test_groups()

  expect_error(
    zoomss_model(bad_input, test_groups),
    "No environmental time series provided"
  )

  # Test with invalid groups
  bad_groups <- data.frame(Species = "Test")
  env_data <- create_test_env_data(n_years = 20, dt = 0.1)
  input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)

  expect_error(
    zoomss_model(input_params, bad_groups),
    "Missing required columns"
  )
})

# Tests for zoomss_params() ----------------------------------------------

test_that("zoomss_params creates valid parameter list", {
  skip_if_not_installed("zoomss")

  env_data <- create_test_env_data(n_years = 20, dt = 0.1)
  input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)
  test_groups <- create_test_groups()

  params <- zoomss_params(test_groups, input_params, isave = 5)

  # Check basic structure
  expect_type(params, "list")
  expect_true("Groups" %in% names(params))
  expect_true("ngrps" %in% names(params))
  expect_true("dt" %in% names(params))
  expect_true("tmax" %in% names(params))

  # Check calculated values
  expect_equal(params$ngrps, nrow(test_groups))
  expect_equal(params$dt, 0.1)
  expect_true(params$tmax > 0)

  # Check group indices
  expect_true("zoo_grps" %in% names(params))
  expect_true("fish_grps" %in% names(params))
  expect_equal(length(params$zoo_grps), sum(test_groups$Type == "Zooplankton"))
  expect_equal(length(params$fish_grps), sum(test_groups$Type == "Fish"))
})

test_that("zoomss_params calculates time parameters correctly", {
  skip_if_not_installed("zoomss")

  # Test different time configurations
  env_data1 <- create_test_env_data(n_years = 20, dt = 0.1)
  input_params1 <- createInputParams(env_data1$time, env_data1$sst, env_data1$chl)
  test_groups <- create_test_groups()

  params1 <- zoomss_params(test_groups, input_params1, isave = 5)
  expect_equal(params1$dt, 0.1)
  expect_equal(params1$tmax, max(env_data1$time))

  # Test with different dt
  env_data2 <- create_test_env_data(n_years = 20, dt = 0.05)
  input_params2 <- createInputParams(env_data2$time, env_data2$sst, env_data2$chl)

  params2 <- zoomss_params(test_groups, input_params2, isave = 10)
  expect_equal(params2$dt, 0.05)
  expect_equal(params2$tmax, max(env_data2$time))
})

test_that("zoomss_params validates uniform time steps", {
  skip_if_not_installed("zoomss")

  test_groups <- create_test_groups()

  # Create non-uniform time steps
  bad_time <- c(0, 0.1, 0.25, 0.3, 0.4)  # Non-uniform steps

  expect_error(
    createInputParams(bad_time, rep(15, 5), rep(2, 5)),
    "Time steps are not uniform"
  )
})

# Tests for zoomss_setup() -----------------------------------------------

test_that("zoomss_setup creates model structure", {
  skip_if_not_installed("zoomss")

  env_data <- create_test_env_data(n_years = 20, dt = 0.1)
  input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)
  test_groups <- create_test_groups()

  params <- zoomss_params(test_groups, input_params, isave = 5)
  model <- zoomss_setup(params)

  # Check model structure
  expect_type(model, "list")
  expect_true("param" %in% names(model))
  expect_identical(model$param, params)

  # Check for feeding kernels
  kernel_names <- c("dynam_growthkernel", "dynam_diffkernel",
                   "dynam_dietkernel", "dynam_mortkernel",
                   "phyto_growthkernel", "phyto_diffkernel", "phyto_dietkernel")

  for (kernel in kernel_names) {
    expect_true(kernel %in% names(model), info = paste("Missing kernel:", kernel))
    expect_true(is.array(model[[kernel]]), info = paste("Kernel not array:", kernel))
  }

  # Check abundance array
  expect_true("N" %in% names(model))
  expect_true(is.array(model$N))
  expect_equal(length(dim(model$N)), 3)  # time x groups x size_classes
})

test_that("zoomss_setup initializes mortality and efficiency matrices", {
  skip_if_not_installed("zoomss")

  env_data <- create_test_env_data(n_years = 20, dt = 0.1)
  input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)
  test_groups <- create_test_groups()

  params <- zoomss_params(test_groups, input_params, isave = 5)
  model <- zoomss_setup(params)

  # Check mortality matrices
  expect_true("M_sb_base" %in% names(model))
  expect_true("fish_mort" %in% names(model))
  expect_true(is.matrix(model$M_sb_base))
  expect_true(is.matrix(model$fish_mort))

  # Check efficiency matrix
  expect_true("assim_eff" %in% names(model))
  expect_true(is.matrix(model$assim_eff))

  # Check dimensions match groups and size classes
  expect_equal(nrow(model$M_sb_base), params$ngrps)
  expect_equal(nrow(model$fish_mort), params$ngrps)
  expect_equal(nrow(model$assim_eff), params$ngrps)
})

# Tests for zoomss_run() -------------------------------------------------

test_that("zoomss_run executes model simulation", {
  skip_if_not_installed("zoomss")

  env_data <- create_test_env_data(n_years = 20, dt = 0.1)
  input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)
  test_groups <- create_test_groups()

  params <- zoomss_params(test_groups, input_params, isave = 5)
  model <- zoomss_setup(params)

  expect_no_error({
    result <- zoomss_run(model)
  })

  # Check result structure
  expect_type(result, "list")
  expect_true("N" %in% names(result))
  expect_true("param" %in% names(result))

  # Check that abundances are finite and non-negative
  expect_true(all(is.finite(result$N)))
  expect_true(all(result$N >= 0))
})

test_that("zoomss_run produces time series output", {
  skip_if_not_installed("zoomss")

  env_data <- create_test_env_data(n_years = 20, dt = 0.1)
  input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)
  test_groups <- create_test_groups()

  params <- zoomss_params(test_groups, input_params, isave = 2)
  model <- zoomss_setup(params)
  result <- zoomss_run(model)

  # Check time series dimensions
  expect_true(dim(result$N)[1] > 1)  # Multiple time steps saved
  expect_equal(dim(result$N)[2], nrow(test_groups))  # Correct number of groups
  expect_equal(dim(result$N)[3], params$ngrid)  # Correct number of size classes

  # Check that we have reasonable output structure
  time_steps_saved <- dim(result$N)[1]
  expected_saves <- floor(params$itimemax / params$isave)
  expect_true(time_steps_saved <= expected_saves + 1)  # Allow for initial condition
})

# Integration tests -------------------------------------------------------

test_that("Full model workflow produces consistent results", {
  skip_if_not_installed("zoomss")

  # Run model twice with same inputs
  env_data <- create_test_env_data(n_years = 20, dt = 0.1)
  input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)
  test_groups <- create_test_groups()

  result1 <- zoomss_model(input_params, test_groups, isave = 5)
  result2 <- zoomss_model(input_params, test_groups, isave = 5)

  # Results should be identical (deterministic model)
  expect_equal(result1$N, result2$N)
  expect_equal(result1$Biomass, result2$Biomass)
})

test_that("Model handles different group configurations", {
  skip_if_not_installed("zoomss")

  env_data <- create_test_env_data(n_years = 20, dt = 0.1)
  input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)

  # Test with minimal groups (1 zoo, 1 fish)
  minimal_groups <- create_test_groups()
  result_minimal <- zoomss_model(input_params, minimal_groups, isave = 5)
  expect_equal(dim(result_minimal$N)[2], 2)  # 2 groups

  # Test with default groups
  result_default <- zoomss_model(input_params, isave = 5)
  expect_true(dim(result_default$N)[2] > 2)  # More groups in default
})



# Edge cases and error handling ------------------------------------------

test_that("Model handles edge cases gracefully", {
  skip_if_not_installed("zoomss")

  # Very short simulation
  env_data_short <- create_test_env_data(n_years = 5, dt = 0.1)
  input_params_short <- createInputParams(env_data_short$time, env_data_short$sst, env_data_short$chl)
  test_groups <- create_test_groups()

  expect_no_error({
    result_short <- zoomss_model(input_params_short, test_groups, isave = 1)
  })

  # Large isave parameter (save less frequently)
  env_data <- create_test_env_data(n_years = 20, dt = 0.1)
  input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)

  expect_no_error({
    result_sparse <- zoomss_model(input_params, test_groups, isave = 50)
  })
})

test_that("Parameter validation catches common errors", {
  skip_if_not_installed("zoomss")

  env_data <- create_test_env_data(n_years = 20, dt = 0.1)
  input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)

  # Test with groups missing required columns
  bad_groups <- data.frame(
    Species = "Test",
    Type = "Zooplankton"
    # Missing other required columns
  )

  expect_error(
    zoomss_model(input_params, bad_groups),
    "Missing required columns"
  )
})
