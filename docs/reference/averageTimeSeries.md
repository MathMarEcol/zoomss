# Calculate mean of final portion of ZooMSS time series

Calculates the mean of the final n years of a time series to obtain
equilibrium values after model spin-up period.

## Usage

``` r
averageTimeSeries(mdl, var, n_years = 10)
```

## Arguments

- mdl:

  ZooMSS model results object containing model parameters and output
  arrays

- var:

  Character string specifying which variable to extract and average
  (e.g., "N", "Growth", "Mort")

- n_years:

  Number of years from the end of the time series to average (default:
  10)

## Value

2D array with averaged values (groups x size_classes)

## Details

Calculate Average Output from Model Time Series

This function removes the initial transient period from time series data
and calculates the mean of the final n years, providing representative
steady-state values. Essential for obtaining equilibrium abundances,
growth rates, and other model outputs after the model has reached
dynamic equilibrium.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run ZooMSS model
results <- zoomss_model(input_params, Groups)

# Average final 3 years of abundance data
avg_abundance <- averageTimeSeries(results, "N", n_years = 3)

# Average final 10 years of growth data (default)
avg_growth <- averageTimeSeries(results, "gg")
} # }
```
