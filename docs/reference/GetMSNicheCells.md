# Get MultiSampNiche cells

Retrieves cell-level data from MultiSampNiche objects for specified
multi-sample analyses.

## Usage

``` r
GetMSNicheCells(IST_obj, loop_id, niche_key, ...)

# S3 method for class 'IST'
GetMSNicheCells(
  IST_obj = NULL,
  loop_id = "LoopAllMulti",
  niche_key = NULL,
  ...
)
```

## Arguments

- IST_obj:

  An IST object

- loop_id:

  Character vector, multi-sample analysis identifiers

- niche_key:

  Character, niche key to retrieve

- ...:

  Additional arguments passed to methods

## Value

List of cell data frames, one per multi-sample analysis

## Examples

``` r
if (FALSE) { # \dontrun{
multi_cells <- GetMSNicheCells(ist_obj,
                               loop_id = c("comparison1", "comparison2"),
                               niche_key = "niche_virulence")
} # }
```
