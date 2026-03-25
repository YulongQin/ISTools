# Get IST package information

Returns information about the IST package including version, available
functions, and current package options.

## Usage

``` r
IST_info()
```

## Value

A list containing package information with the following components:

- `version`: Package version

- `functions`: List of available functions by category

- `supported_formats`: Supported spatial transcriptomics formats

- `supported_hosts`: Supported host species

- `supported_pathogens`: Supported pathogen types

- `score_methods`: Available scoring methods

- `blur_methods`: Available blurring methods

- `options`: Current package options

## Examples

``` r
if (FALSE) { # \dontrun{
info <- IST_info()
print(info$version)
print(info$supported_formats)
} # }
```
