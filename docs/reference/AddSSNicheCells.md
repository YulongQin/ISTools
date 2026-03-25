# Add columns to SingleSampNiche cells

Adds new columns to the cell data in SingleSampNiche objects.

## Usage

``` r
AddSSNicheCells(IST_obj, loop_id, meta_key, select_colnm, niche_key, ...)

# S3 method for class 'IST'
AddSSNicheCells(
  IST_obj = NULL,
  loop_id = "LoopAllSamp",
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

  Character vector, sample identifiers

- meta_key:

  Character, metadata key containing new data

- select_colnm:

  Character vector, columns to add

- niche_key:

  Character, niche key to modify

- ...:

  Additional arguments passed to methods

## Value

Modified IST object with updated niche cell data

## Examples

``` r
if (FALSE) { # \dontrun{
# Add new annotation columns to niche cells
IST_obj <- AddSSNicheCells(ist_obj,
                           loop_id = "LoopAllSamp",
                           meta_key = "custom_annotation",
                           select_colnm = c("cell_type", "confidence"),
                           niche_key = "niche_virulence")
} # }
```
