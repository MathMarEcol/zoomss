# Aggregate ZooMSS abundances across all size classes

Sums abundance values across all size classes for each functional group,
providing total abundance per group.

## Usage

``` r
reduceSize(x, method = "sum")
```

## Arguments

- x:

  3D array outptut from ZooMSS model

- method:

  Character string specifying aggregation method: "sum" (default) or
  "mean".

## Value

List of vectors with total abundance per functional group

## Details

Sum ZooMSS Output Across Size Bins

This function collapses the size dimension of ZooMSS output by summing
across all size classes. Useful for analyzing total abundance patterns
without size structure detail.

## Examples

``` r
if (FALSE) { # \dontrun{
# After running ZooMSS model
results <- zoomss_model(input_params, Groups)
total_abundances <- reduceSize(results$abundances)
} # }
```
