# Validate ZooMSS functional groups data structure and values

Performs comprehensive validation of functional groups data to ensure it
meets ZooMSS model requirements.

## Usage

``` r
validateGroups(groups)
```

## Arguments

- groups:

  Data frame containing functional groups data

## Value

TRUE if validation passes (invisibly), otherwise throws an error

## Details

Validate Functional Groups Data

This function validates:

- Required column names are present

- Data types are correct

- Parameter values are within reasonable ranges

- No missing values in critical columns

- Size ranges are logical (W0 \< Wmax)

## Examples

``` r
if (FALSE) { # \dontrun{
Groups <- getGroups()
validateGroups(Groups)  # Should pass

# This would fail validation:
bad_groups <- Groups
bad_groups$W0 <- NULL
validateGroups(bad_groups)  # Error: missing required column
} # }
```
