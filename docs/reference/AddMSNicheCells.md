# Add columns to MultiSampNiche cells

Adds new columns to the cell data in MultiSampNiche objects.

## Usage

``` r
AddMSNicheCells(IST_obj, loop_id, meta_key, select_colnm, niche_key, ...)

# S3 method for class 'IST'
AddMSNicheCells(
  IST_obj = NULL,
  loop_id = "LoopAllMulti",
  meta_key = NULL,
  select_colnm = NULL,
  niche_key = NULL,
  ...
)
```

## Arguments

- IST_obj:

  An IST object

- loop_id:

  Character vector, multi-sample analysis identifiers

- meta_key:

  Character, metadata key containing new data

- select_colnm:

  Character vector, columns to add

- niche_key:

  Character, niche key to modify

- ...:

  Additional arguments passed to methods

## Value

Modified IST object with updated multi-sample niche cell data

## Examples

``` r
if (FALSE) { # \dontrun{
# Add new annotations to multi-sample niche cells
IST_obj <- AddMSNicheCells(ist_obj,
                           loop_id = "LoopAllMulti",
                           meta_key = "custom_annotation",
                           select_colnm = c("cell_type", "confidence"),
                           niche_key = "niche_virulence")
} # }
```
