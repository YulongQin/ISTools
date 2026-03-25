# Create a MultiSampNiche object for cross-sample comparison

Creates and stores a MultiSampNiche object within an IST object,
combining niche analysis results from multiple samples for comparative
analysis.

## Usage

``` r
CreateMultiSampNiche(
  IST_obj = NULL,
  multi_id = NULL,
  loop_id = "LoopAllSamp",
  compare_mode = NULL,
  niche_key = NULL,
  description = NULL
)

# S3 method for class 'IST'
CreateMultiSampNiche(
  IST_obj = NULL,
  multi_id = NULL,
  loop_id = "LoopAllSamp",
  compare_mode = NULL,
  niche_key = NULL,
  description = NULL
)

CreateMSNiche(...)
```

## Arguments

- IST_obj:

  An IST object containing SingleSampNiche objects

- multi_id:

  Character, unique identifier for the multi-sample analysis

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- compare_mode:

  Character, comparison mode - "Comparative" or "Temporal"

- niche_key:

  Character, niche key to combine across samples

- description:

  Character, description of the multi-sample analysis

- ...:

  Additional arguments passed to methods

## Value

Modified IST object with added MultiSampNiche information

## Examples

``` r
if (FALSE) { # \dontrun{
# Create MultiSampNiche object for comparative analysis
IST_obj <- CreateMultiSampNiche(
  IST_obj = ist_obj,
  multi_id = "comparison_infected_vs_control",
  loop_id = c("infected1", "infected2", "control1", "control2"),
  compare_mode = "Comparative",
  niche_key = "niche_virulence"
)
} # }
```
