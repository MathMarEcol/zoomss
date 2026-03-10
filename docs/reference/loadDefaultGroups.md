# Internal function to load default ZooMSS groups

Loads the default functional groups from the package data or CSV file.
This is an internal function used by getGroups().

## Usage

``` r
loadDefaultGroups()
```

## Value

Data frame with default functional groups

## Details

Load Default Functional Groups Data

This function handles the actual loading of default groups data, whether
from package data (if available) or from the CSV file in data-raw.
