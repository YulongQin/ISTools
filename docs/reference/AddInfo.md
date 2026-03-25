# Add information to IST object

Appends information to existing slots in the IST_info, primarily for
adding comments and notes.

## Usage

``` r
AddInfo(IST_obj, info_key, sub_key, info_value, ...)

# S3 method for class 'IST'
AddInfo(
  IST_obj = NULL,
  info_key = "comment_info",
  sub_key = NULL,
  info_value = NULL,
  ...
)
```

## Arguments

- IST_obj:

  An IST object

- info_key:

  Character, information category (default: "comment_info")

- sub_key:

  Character, specific sub-key to add to

- info_value:

  Value to append

- ...:

  Additional arguments passed to methods

## Value

Modified IST object

## Examples

``` r
if (FALSE) { # \dontrun{
# Add a comment about analysis parameters
IST_obj <- AddInfo(ist_obj,
                   info_key = "comment_info",
                   sub_key = "parameters",
                   info_value = "Used DBSCAN with eps=5, minPts=10")
} # }
```
