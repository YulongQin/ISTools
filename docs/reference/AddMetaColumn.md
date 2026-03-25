# Add columns to existing metadata

Adds new columns to an existing metadata data frame in the IST object.

## Usage

``` r
AddMetaColumn(IST_obj, ...)

# S3 method for class 'IST'
AddMetaColumn(
  IST_obj = NULL,
  meta_key = NULL,
  add_data = NULL,
  ignore_rownm = FALSE,
  ...
)
```

## Arguments

- IST_obj:

  An IST object

- ...:

  Additional arguments passed to methods

- meta_key:

  Character, metadata key to modify

- add_data:

  Data frame, new columns to add

- ignore_rownm:

  Logical, whether to ignore row name matching (default: FALSE)

## Value

Modified IST object with added metadata columns

## Examples

``` r
if (FALSE) { # \dontrun{
# Add new columns to raw metadata
new_cols <- data.frame(new_score = runif(ncol(ist_obj)))
IST_obj <- AddMetaColumn(ist_obj,
                         meta_key = "raw",
                         add_data = new_cols)
} # }
```
