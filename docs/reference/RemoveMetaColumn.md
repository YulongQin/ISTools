# Remove columns from existing metadata

Removes specified columns from an existing metadata data frame in the
IST object.

## Usage

``` r
RemoveMetaColumn(IST_obj, ...)

# S3 method for class 'IST'
RemoveMetaColumn(IST_obj = NULL, meta_key = NULL, remove_colnm = NULL, ...)
```

## Arguments

- IST_obj:

  An IST object

- ...:

  Additional arguments passed to methods

- meta_key:

  Character, metadata key to modify

- remove_colnm:

  Character vector, column names to remove

## Value

Modified IST object with removed metadata columns

## Examples

``` r
if (FALSE) { # \dontrun{
# Remove temporary columns
IST_obj <- RemoveMetaColumn(ist_obj,
                            meta_key = "raw",
                            remove_colnm = c("temp_score", "temp_label"))
} # }
```
