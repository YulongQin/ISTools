# Set IST package options

Convenience function to set multiple IST package options at once.

## Usage

``` r
IST_options(...)
```

## Arguments

- ...:

  Named options to set. Valid options include:

  - `verbose`: Logical, print verbose output

  - `parallel_workers`: Integer, number of parallel workers

  - `project`: Character, default project name

  - `default_assay`: Character, default assay name

  - `check_python`: Logical, check Python configuration

## Value

NULL (invisibly)

## Examples

``` r
if (FALSE) { # \dontrun{
# Set multiple options at once
IST_options(verbose = FALSE, parallel_workers = 8)

# Set project name
IST_options(project = "MyProject")
} # }
```
