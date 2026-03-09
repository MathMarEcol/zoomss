# Generate synthetic environmental data for ZooMSS testing

Creates simple synthetic environmental time series with optional
seasonal variation for testing ZooMSS model runs when real environmental
data is not available.

## Usage

``` r
createEnviroData(
  n_years,
  dt,
  base_sst = 15,
  base_chl = 0.5,
  seasonal = TRUE,
  sst_amplitude = 3,
  chl_amplitude = 0.2
)
```

## Arguments

- n_years:

  Number of years to generate

- dt:

  Time step size in years

- base_sst:

  Base sea surface temperature in deg C (default: 15)

- base_chl:

  Base chlorophyll concentration in mg/m^3 (default: 0.5)

- seasonal:

  Logical, whether to add seasonal variation (default: TRUE)

- sst_amplitude:

  Amplitude of SST seasonal variations in deg C (default: 3)

- chl_amplitude:

  Amplitude of chlorophyll seasonal variations in mg/m^3 (default: 0.2)

## Value

Data frame with columns: time, sst, chl

## Details

Create Environmental Time Series

This function generates synthetic sea surface temperature and
chlorophyll time series that can be used for testing ZooMSS model
behavior. The function can create either static environmental conditions
or seasonal cycles with sinusoidal variation. This is particularly
useful for:

- Testing model sensitivity to environmental forcing

- Creating idealized scenarios for model exploration

- Generating data when real environmental data is unavailable

The seasonal option creates SST and chlorophyll cycles that are out of
phase, mimicking typical ocean patterns where chlorophyll peaks when SST
is lower.

## Examples

``` r
# Create seasonal environmental data
env_data <- createEnviroData(
  n_years = 10,
  dt = 0.01,
  seasonal = TRUE
)

# Create static environmental conditions
static_data <- createEnviroData(
  n_years = 5,
  dt = 0.01,
  seasonal = FALSE,
  base_sst = 20,
  base_chl = 1.0
)
```
