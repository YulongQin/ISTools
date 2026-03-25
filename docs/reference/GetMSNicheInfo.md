# Get MultiSampNiche information

Retrieves niche information from MultiSampNiche objects for specified
multi-sample analyses.

## Usage

``` r
GetMSNicheInfo(IST_obj, loop_id, niche_key, ...)

# S3 method for class 'IST'
GetMSNicheInfo(IST_obj = NULL, loop_id = "LoopAllMulti", niche_key = NULL, ...)
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

List of niche information data frames, one per multi-sample analysis

## Examples

``` r
if (FALSE) { # \dontrun{
multi_info <- GetMSNicheInfo(ist_obj,
                             loop_id = c("comparison1", "comparison2"),
                             niche_key = "niche_virulence")
} # }
```
