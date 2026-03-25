# IST package options

Options that control the behavior of the IST package. These can be set
via [`options`](https://rdrr.io/r/base/options.html).

## Options

- `IST.project`:

  Default project name (default: "IST_Project")

- `IST.verbose`:

  Whether to print verbose output (default: TRUE)

- `IST.check_python`:

  Whether to check Python configuration (default: TRUE)

- `IST.default_assay`:

  Default assay name (default: "Spatial")

- `IST.parallel_workers`:

  Number of parallel workers (default: 4)

- `IST.temp_dir`:

  Temporary directory (default: tempdir())

## Examples

``` r
if (FALSE) { # \dontrun{
options(IST.verbose = FALSE)
options(IST.parallel_workers = 8)
getOption("IST.project")
} # }
```
