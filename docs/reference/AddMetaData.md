# Add metadata to IST object

Adds a new metadata data frame to the IST object's meta_data_record and
updates the tracking information.

## Usage

``` r
AddMetaData(IST_obj, ...)

# S3 method for class 'IST'
AddMetaData(
  IST_obj = NULL,
  meta_key = NULL,
  add_data = NULL,
  dir_nm = NA,
  grp_nm = NA,
  asso_key = NULL,
  description = NULL,
  ...
)
```

## Arguments

- IST_obj:

  An IST object

- ...:

  Additional arguments passed to methods

- meta_key:

  Character, unique identifier for the metadata

- add_data:

  Data frame, metadata to add

- dir_nm:

  Character, directory name for output (optional)

- grp_nm:

  Character, group name for output (optional)

- asso_key:

  Character, associated metadata key (optional)

- description:

  Character, description of the metadata (optional)

## Value

Modified IST object with added metadata

## Examples

``` r
if (FALSE) { # \dontrun{
# Add custom metadata
IST_obj <- AddMetaData(ist_obj,
                       meta_key = "custom_annotation",
                       add_data = custom_metadata,
                       description = "Manual cell annotations")
} # }
```
