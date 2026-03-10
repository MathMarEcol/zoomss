# Convert tibble to data frame for efficiency

Removes tibble attributes and converts to a plain data frame for
improved speed and memory efficiency in computational workflows.

## Usage

``` r
untibble(tibble)
```

## Arguments

- tibble:

  A tibble or data frame object to convert

## Value

Plain data frame without tibble attributes

## Details

Remove Tibble Attributes

This utility function strips tibble-specific attributes that can slow
down operations in tight computational loops. Used internally by ZooMSS
for performance optimization when working with large datasets.
