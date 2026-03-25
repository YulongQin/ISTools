# Automated ROI detection using density-based spatial clustering

This function automatically detects regions of interest (ROIs) using
density-based clustering (DBSCAN) on spatial transcriptomics data. It
includes preprocessing steps like density filtering and iterative spot
updating.

## Usage

``` r
NicheDetect_STS(
  IST_obj = NULL,
  meta_key = NULL,
  loop_id = "LoopAllSamp",
  method_level = "region",
  method_region = "convex",
  update_spots = TRUE,
  pos_colnm = NULL,
  neg_value = "neg",
  density_thres = 0.9,
  ROI_size = 10,
  minPts = NULL,
  k_kNNdist = NULL,
  description = NULL,
  grp_nm = NULL,
  dir_nm = "M2_NicheDetect_STS"
)
```

## Arguments

- IST_obj:

  An IST object containing spatial transcriptomics data

- meta_key:

  Character, metadata key containing positive spot information

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- method_level:

  Character, detection level - "region" or "spot" (default: "region")

- method_region:

  Character, hull method - "convex" or "concave" (default: "convex")

- update_spots:

  Logical, whether to iteratively update positive spots (default: TRUE)

- pos_colnm:

  Character, column name containing positive spot labels

- neg_value:

  Character, value indicating negative/non-ROI spots (default: "neg")

- density_thres:

  Numeric, density threshold for filtering (0-1, default: 0.9)

- ROI_size:

  Integer, minimum number of spots to form an ROI (default: 10)

- minPts:

  Integer, minimum points parameter for DBSCAN (default: NULL,
  automatically determined based on data format)

- k_kNNdist:

  Integer, k value for kNN distance calculation (default: NULL, uses
  minPts if not specified)

- description:

  Character, description of the analysis (default: NULL)

- grp_nm:

  Character, group name for output organization (default: NULL, uses
  timestamp)

- dir_nm:

  Character, directory name for output (default: "M2_NicheDetect_STS")

## Value

Returns the modified IST object with added metadata containing ROI
labels, distances to ROI centers, edge information, and region
classifications

## Examples

``` r
if (FALSE) { # \dontrun{
# Automatically detect ROIs based on positive spots
IST_obj <- NicheDetect_STS(
  IST_obj = ist_object,
  meta_key = "M1_SpotDetect_Gene_20240101",
  pos_colnm = "Label_geneA",
  density_thres = 0.9,
  ROI_size = 10
)
} # }
```
