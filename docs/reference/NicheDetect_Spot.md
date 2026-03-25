# Spot-level niche detection

Performs niche detection at the individual spot level, treating each
positive spot as its own ROI. Calculates distances from all spots to the
nearest positive spot.

## Usage

``` r
NicheDetect_Spot(
  IST_obj = NULL,
  loop_id = "LoopAllSamp",
  pos_colnm = NULL,
  neg_value = "neg",
  meta_key = NULL,
  description = NULL
)
```

## Arguments

- IST_obj:

  An IST object containing spatial transcriptomics data

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- pos_colnm:

  Character, column name containing positive spot labels

- neg_value:

  Character, value indicating negative spots (default: "neg")

- meta_key:

  Character, metadata key containing positive spot information

- description:

  Character, description of the analysis (default: NULL)

## Value

Returns the modified IST object with added metadata containing ROI
labels and distances to nearest positive spots

## Examples

``` r
if (FALSE) { # \dontrun{
# Perform spot-level niche detection
IST_obj <- NicheDetect_Spot(
  IST_obj = ist_object,
  meta_key = "M1_SpotDetect_Gene_20240101",
  pos_colnm = "Label_geneA"
)
} # }
```
