# Create a SingleSampNiche object from analysis results

Creates and stores a SingleSampNiche object within an IST object,
organizing niche detection results for a single sample. This function
extracts relevant cell and gene information based on specified niche
parameters.

## Usage

``` r
CreateSingleSampNiche(
  IST_obj = NULL,
  loop_id = "LoopAllSamp",
  niche_key = NULL,
  meta_key = NULL,
  ROI_type = NULL,
  pos_colnm = NULL,
  neg_value = "neg",
  center_colnm = NULL,
  edge_colnm = NULL,
  all_label_colnm = NULL,
  all_dist_colnm = NULL,
  other_colnm = NULL,
  description = NULL
)

# S3 method for class 'IST'
CreateSingleSampNiche(
  IST_obj = NULL,
  loop_id = "LoopAllSamp",
  niche_key = NULL,
  meta_key = NULL,
  ROI_type = NULL,
  pos_colnm = NULL,
  neg_value = "neg",
  center_colnm = NULL,
  edge_colnm = NULL,
  all_label_colnm = NULL,
  all_dist_colnm = NULL,
  other_colnm = NULL,
  description = NULL
)

CreateSSNiche(...)
```

## Arguments

- IST_obj:

  An IST object containing analysis results

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- niche_key:

  Character, unique identifier for the niche analysis

- meta_key:

  Character, metadata key containing the source data

- ROI_type:

  Character, type of ROI - "ROI", "Spot", or "Unknown"

- pos_colnm:

  Character, column name containing positive spot labels

- neg_value:

  Character, value indicating negative spots (default: "neg")

- center_colnm:

  Character, column name indicating center spots (optional)

- edge_colnm:

  Character, column name indicating edge spots (optional)

- all_label_colnm:

  Character, column name for all spot labels (optional)

- all_dist_colnm:

  Character, column name for distance to center (optional)

- other_colnm:

  Character vector, additional columns to include (optional)

- description:

  Character, description of the niche analysis (optional)

- ...:

  Additional arguments passed to methods

## Value

Modified IST object with added SingleSampNiche information

## Examples

``` r
if (FALSE) { # \dontrun{
# Create SingleSampNiche object
IST_obj <- CreateSingleSampNiche(
  IST_obj = ist_obj,
  niche_key = "niche_virulence",
  meta_key = "M2_NicheDetect_STS_20240101",
  ROI_type = "ROI",
  pos_colnm = "ROI_label",
  center_colnm = "ROI_center",
  edge_colnm = "ROI_edge",
  all_label_colnm = "All_ROI_label",
  all_dist_colnm = "All_Dist2ROIcenter"
)
} # }
```
