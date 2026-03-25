# Calculate Niche Co-localization using MISTyR

Performs spatial co-localization analysis using the MISTyR framework to
identify intra-view, juxtaview, and paraview relationships between cell
types or gene expression patterns within niches.

## Usage

``` r
CalNicheCoLoc(
  IST_obj = NULL,
  loop_id = "LoopAllSamp",
  niche_key = NULL,
  meta_key = NULL,
  group_by = NULL,
  group_use = NULL,
  features = NULL,
  feature_colnm = NULL,
  juxtaview_radius = 15,
  paraview_radius = 10,
  heatmap_cutoff = 0,
  comm_cutoff = 1,
  return_data = TRUE,
  grp_nm = NULL,
  dir_nm = "M3_CalNicheCoLoc"
)
```

## Arguments

- IST_obj:

  An IST object containing niche analysis results

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- niche_key:

  Character, niche key to analyze

- meta_key:

  Character, metadata key for analysis when niche_key is NULL

- group_by:

  Character, column name for cell type grouping

- group_use:

  Character vector, specific groups to include

- features:

  Character vector, gene names to analyze

- feature_colnm:

  Character vector, metadata column names to analyze

- juxtaview_radius:

  Numeric, radius for juxtaview in spatial units (default: 15)

- paraview_radius:

  Numeric, radius for paraview in spatial units (default: 10)

- heatmap_cutoff:

  Numeric, cutoff for interaction heatmap display (default: 0)

- comm_cutoff:

  Numeric, cutoff for community detection (default: 1)

- return_data:

  Logical, whether to return the results list (default: TRUE)

- grp_nm:

  Character, group name for output organization (default: NULL)

- dir_nm:

  Character, directory name for output (default: "M3_CalNicheCoLoc")

## Value

If return_data = TRUE, returns a list of MISTyR results per sample

## Examples

``` r
if (FALSE) { # \dontrun{
# Analyze co-localization of cell types within niches
results <- CalNicheCoLoc(
  IST_obj = ist_obj,
  niche_key = "niche_virulence",
  group_by = "cell_type",
  juxtaview_radius = 20,
  paraview_radius = 15
)
} # }
```
