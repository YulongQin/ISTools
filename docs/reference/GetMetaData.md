# Get metadata from IST object

Retrieves metadata from the IST object's meta_data_record, optionally
combining multiple metadata sources and adding coordinate information.

## Usage

``` r
GetMetaData(IST_obj, meta_key, ...)

# S3 method for class 'IST'
GetMetaData(IST_obj = NULL, meta_key = NULL, add_coord = TRUE, ...)
```

## Arguments

- IST_obj:

  An IST object

- meta_key:

  Character vector, metadata keys to retrieve

- ...:

  Additional arguments passed to methods

- add_coord:

  Logical, whether to add coordinate columns (default: TRUE)

## Value

List of metadata data frames, one for each meta_key

## Examples

``` r
if (FALSE) { # \dontrun{
# Get raw metadata
raw_meta <- GetMetaData(ist_obj, meta_key = "raw")

# Get multiple metadata sources combined
combined_meta <- GetMetaData(ist_obj,
                             meta_key = list(c("raw", "coord"),
                                            c("M1_SpotDetect_Gene_20240101")))
} # }
```
