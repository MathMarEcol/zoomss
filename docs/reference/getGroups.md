# Load default or custom functional groups for ZooMSS model

Provides access to the default ZooMSS functional groups or loads custom
groups from a file. This function is the primary way to obtain Groups
data for ZooMSS model runs.

## Usage

``` r
getGroups(source = "default", file = NULL)
```

## Arguments

- source:

  Character string specifying data source. Options:

  - "default": Use built-in ZooMSS functional groups (default)

  - "file": Load groups from a CSV file

  - "template": Export default groups to a file for modification

- file:

  Path to CSV file when source="file" or source="template"

## Value

Data frame containing functional groups with required columns: Species,
Type, W0, Wmax, and other biological parameters

## Details

Get Default ZooMSS Functional Groups

This function provides flexible access to functional groups data:

- **Default groups**: Returns the standard ZooMSS functional groups (12
  groups)

- **Custom file**: Loads and validates groups from a user-provided CSV
  file

- **Template creation**: Exports default groups to a file for user
  modification

The default groups include: Flagellates, Ciliates, Larvaceans,
OmniCopepods, CarnCopepods, Euphausiids, Chaetognaths, Salps, Jellyfish,
and three Fish groups (Small, Medium, Large).

All groups data is validated to ensure it contains required columns and
reasonable parameter values for successful model runs.

## Examples

``` r
if (FALSE) { # \dontrun{
# Use default groups
Groups <- getGroups()

# Create a template file for modification
getGroups(source = "template", file = "my_groups.csv")

# Load custom groups from file
custom_groups <- getGroups(source = "file", file = "my_groups.csv")

# Modify default groups programmatically
Groups <- getGroups()
Groups$W0[Groups$Species == "Flagellates"] <- -12.5  # Modify minimum size
} # }
```
