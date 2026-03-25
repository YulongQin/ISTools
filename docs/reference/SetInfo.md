# Set information in IST object

Modifies specific information in the IST object's IST_info slot.

## Usage

``` r
SetInfo(IST_obj, info_key, sub_key, info_value, ...)

# S3 method for class 'IST'
SetInfo(
  IST_obj = NULL,
  info_key = NULL,
  sub_key = NULL,
  info_value = NULL,
  ...
)
```

## Arguments

- IST_obj:

  An IST object

- info_key:

  Character, main information category

- sub_key:

  Character, specific sub-key to set (optional)

- info_value:

  Value to set

- ...:

  Additional arguments passed to methods

## Value

Modified IST object

## Examples

``` r
if (FALSE) { # \dontrun{
# Update project description
IST_obj <- SetInfo(ist_obj,
                   info_key = "project_info",
                   sub_key = "description",
                   info_value = "Updated project description")
} # }
```
