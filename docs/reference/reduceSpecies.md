# Aggregate ZooMSS abundances across all species

Aggregates abundance values across all species bins for each functional
group and size class using the specified method.

## Usage

``` r
reduceSpecies(x, method = "sum")
```

## Arguments

- x:

  3D array outptut from ZooMSS model

- method:

  Character string specifying aggregation method: "sum" (default) or
  "mean".

## Value

Array with species dimension reduced using the specified method

## Details

Aggregate ZooMSS abundances across all species

This function collapses the species dimension by applying the specified
method (sum or mean) across all species bins.
