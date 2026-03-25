# Remove metadata from IST object

Removes specified metadata from the IST object's meta_data_record.

## Usage

``` r
RemoveMetaData(IST_obj, meta_key, ...)

# S3 method for class 'IST'
RemoveMetaData(IST_obj = NULL, meta_key = NULL, ...)
```

## Arguments

- IST_obj:

  An IST object

- meta_key:

  Character vector, metadata keys to remove

- ...:

  Additional arguments passed to methods

## Value

Modified IST object with removed metadata

## Examples

``` r
if (FALSE) { # \dontrun{
# Remove temporary metadata
IST_obj <- RemoveMetaData(ist_obj, meta_key = "temp_analysis")
} # }
```
