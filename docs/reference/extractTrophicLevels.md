# Compute trophic levels for functional groups using diet composition

Calculates trophic levels for each functional group based on their diet
composition using an iterative Gauss-Seidel algorithm.

## Usage

``` r
extractTrophicLevels(mdl)
```

## Arguments

- mdl:

  ZooMSS model results object containing 3D diet data (mdl\$diet).
  Dimensions are time, groups, prey_items where columns 1:3 are always
  phytoplankton size classes and remaining columns are zooplankton/fish
  groups.

## Value

Matrix where rows are time steps and columns are functional groups

## Details

Calculate Trophic Levels from Diet Matrix

This function computes trophic levels by:

- Starting with phytoplankton at trophic level 1.0

- Initializing all other groups at trophic level 2.0

- Iteratively updating trophic levels based on weighted diet composition

- Continuing until convergence (difference \< 0.01) or maximum
  iterations (100)

- Processing 3D diet arrays with time series data

Trophic level calculation follows: TL = 1 + sum(diet_fraction_i \*
TL_prey_i)

The function calculates trophic levels for each time step separately and
dynamically determines the number of groups from the diet matrix
dimensions.

This provides a quantitative measure of each group's position in the
food web and is useful for analyzing ecosystem structure and energy
transfer efficiency.

## Examples

``` r
if (FALSE) { # \dontrun{
# After running ZooMSS model with 3D time series
results <- zoomss_model(input_params, Groups)
trophic_levels <- extractTrophicLevels(results)  # Returns matrix (time x groups)

# View trophic levels by group for final time step
final_tl <- trophic_levels[nrow(trophic_levels), ]
names(final_tl) <- results$param$Groups$Species
print(final_tl)
} # }
```
