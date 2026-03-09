# Aggregate abundances across all groups and size classes

Calculates total abundance across all functional groups and size classes
using the specified method.

## Usage

``` r
reduceAll(x, method = "sum")
```

## Arguments

- x:

  3D array outptut from ZooMSS model

- method:

  Character string specifying aggregation method: "sum" (default) or
  "mean".

## Value

Vector of total abundance values (one per spatial cell)

## Details

Aggregate abundances across all groups and size classes

This function provides the most aggregated view of ZooMSS output by
applying the method across both functional groups and size classes.
