# Expand niche regions based on distance threshold

Expands existing niche regions by including neighboring spots within a
specified distance threshold. Useful for defining broader influence
zones around detected niches.

## Usage

``` r
NicheExpand(
  IST_obj = NULL,
  meta_key = NULL,
  loop_id = "LoopAllSamp",
  pos_colnm = NULL,
  neg_value = "neg",
  center_colnm = NULL,
  expand_dist = 1,
  description = NULL
)
```

## Arguments

- IST_obj:

  An IST object containing niche detection results

- meta_key:

  Character, metadata key containing niche information

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- pos_colnm:

  Character, column name containing positive spot labels

- neg_value:

  Character, value indicating negative spots (default: "neg")

- center_colnm:

  Character, column name indicating niche center spots (default: NULL)

- expand_dist:

  Numeric, distance threshold for expansion (default: 1)

- description:

  Character, description of the analysis (default: NULL)

## Value

Returns the modified IST object with expanded ROI labels and updated
distance calculations

## Examples

``` r
if (FALSE) { # \dontrun{
# Expand niches by including spots within distance 5
IST_obj <- NicheExpand(
  IST_obj = ist_object,
  meta_key = "M2_NicheDetect_STS_20240101",
  pos_colnm = "ROI_label",
  expand_dist = 5
)
} # }
```
