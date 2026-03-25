# Get SingleSampNiche cells

Retrieves cell-level data from SingleSampNiche objects for specified
samples.

## Usage

``` r
GetSSNicheCells(IST_obj, loop_id, niche_key, ...)

# S3 method for class 'IST'
GetSSNicheCells(IST_obj = NULL, loop_id = "LoopAllSamp", niche_key = NULL, ...)
```

## Arguments

- IST_obj:

  An IST object

- loop_id:

  Character vector, sample identifiers

- niche_key:

  Character, niche key to retrieve

- ...:

  Additional arguments passed to methods

## Value

List of cell data frames, one per sample

## Examples

``` r
if (FALSE) { # \dontrun{
niche_cells <- GetSSNicheCells(ist_obj,
                               loop_id = c("sample1", "sample2"),
                               niche_key = "niche_virulence")
} # }
```
