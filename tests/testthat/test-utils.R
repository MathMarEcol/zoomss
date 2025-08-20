library(testthat)
library(zoomss)
library(dplyr)
library(tibble)

# Setup real model data once for all tests  
enviro_data <- createEnviroData(n_years = 20)
groups <- getGroups() %>% 
  dplyr::filter(!Species %in% c("Apex_Fish", "Demersal_Fish", "Small_Pelagic_Fish"))

# Create real model data with isave parameter for full functionality  
mdl <- zoomss_model(
  input_params = enviro_data,
  Groups = groups,
  isave = 50  # This is important for proper test results
)

# Create a second model with zero abundance for edge case testing
mdl2 <- mdl
mdl2$N[] <- 0

# Tests for reduceSize() -------------------------------------------------

test_that("reduceSize reduces size dimension correctly", {
  # Test with real model data
  result <- reduceSize(mdl$N)
  
  # Check dimensions: should reduce last dimension (sizes)
  expect_equal(dim(result), c(dim(mdl$N)[1], dim(mdl$N)[2]))
  
  # Check that values are sums across size dimension
  manual_sum <- apply(mdl$N, c(1, 2), sum)
  expect_equal(result, manual_sum)
})

test_that("reduceSize handles 2D arrays", {
  # Test with 2D array (single time slice from model)
  result <- reduceSize(mdl$N[1, , ])
  
  # For 2D array, reduceSize returns the same array (no size dimension to reduce)
  expect_equal(result, mdl$N[1, , ])
})

# Tests for reduceSpecies() ----------------------------------------------

test_that("reduceSpecies reduces group dimension correctly", {
  # Test with real model data
  result <- reduceSpecies(mdl$N)
  
  # Check dimensions: should reduce middle dimension (groups)
  expect_equal(dim(result), c(dim(mdl$N)[1], dim(mdl$N)[3]))
  
  # Check that values are sums across group dimension
  manual_sum <- apply(mdl$N, c(1, 3), sum)
  expect_equal(result, manual_sum)
})

test_that("reduceSpecies handles 2D arrays", {
  # Test with 2D array (single time slice from model)
  # reduceSpecies is designed for 3D arrays, so 2D should error
  expect_error(reduceSpecies(mdl$N[1, , ]), "'MARGIN' does not match dim")
})

# Tests for reduceAll() --------------------------------------------------

test_that("reduceAll reduces all dimensions correctly", {
  # Test with real model data
  result <- reduceAll(mdl$N)
  
  # Should reduce to vector (time dimension only)
  expect_true(is.vector(result))
  expect_equal(length(result), dim(mdl$N)[1])
  
  # Check that values are sums across all but first dimension
  manual_sum <- apply(mdl$N, 1, sum)
  expect_equal(result, manual_sum)
})

# Tests for getBiomass() --------------------------------------------------

test_that("getBiomass converts abundance to wet weight biomass correctly", {
  # Test with real model data
  result <- getBiomass(mdl, units = "ww")
  
  # Check dimensions match input
  expect_equal(dim(result), dim(mdl$N))
  
  # Check that biomass = abundance * weight (approximately, allowing for rounding)
  # Test a few specific cases
  idx <- which(mdl$N > 0, arr.ind = TRUE)[seq_len(5), ]
  for(i in seq_len(nrow(idx))) {
    t_idx <- idx[i, 1]
    g_idx <- idx[i, 2] 
    s_idx <- idx[i, 3]
    expected_biomass <- mdl$N[t_idx, g_idx, s_idx] * mdl$param$w[s_idx]
    expect_equal(result[t_idx, g_idx, s_idx], expected_biomass, tolerance = 1e-10)
  }
})

# Tests for untibble() ---------------------------------------------------

test_that("untibble converts tibble to data frame", {
  # Create a test tibble
  test_tibble <- tibble::tibble(x = 1:5, y = letters[1:5])
  
  result <- untibble(test_tibble)
  
  # Check that result is a data frame, not a tibble
  expect_true(is.data.frame(result))
  expect_false(tibble::is_tibble(result))
  
  # Check that data is preserved
  expect_equal(result$x, test_tibble$x)
  expect_equal(result$y, test_tibble$y)
})
