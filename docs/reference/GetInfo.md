# Get information from IST object

Retrieves specific information from the IST object's IST_info slot.

## Usage

``` r
GetInfo(IST_obj, info_key, ...)

# S3 method for class 'IST'
GetInfo(IST_obj = NULL, info_key = NULL, sub_key = NULL, ...)
```

## Arguments

- IST_obj:

  An IST object

- info_key:

  Character, main information category (e.g., "samp_info", "data_info")

- ...:

  Additional arguments passed to methods

- sub_key:

  Character vector, specific sub-keys to retrieve (optional)

## Value

Requested information from the IST object

## Examples

``` r
if (FALSE) { # \dontrun{
# Get sample column name
samp_colnm <- GetInfo(ist_obj, info_key = "data_info", sub_key = "samp_colnm")

# Get all data information
data_info <- GetInfo(ist_obj, info_key = "data_info")
} # }
```
