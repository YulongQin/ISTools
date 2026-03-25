# Calculate Niche Cell Type Composition

Analyzes and visualizes the cell type composition within and outside
niches. Generates stacked barplots and faceted barplots showing cell
type proportions for niche vs. bystander regions.

## Usage

``` r
CalNicheComposition(
  IST_obj = NULL,
  samp_type = "SS",
  loop_id = "LoopAllSamp",
  samp_grp_index = FALSE,
  meta_key = NULL,
  niche_key = NULL,
  group_by = NULL,
  col = COLOR_LIST[["PALETTE_WHITE_BG"]],
  return_data = FALSE
)
```

## Arguments

- IST_obj:

  An IST object containing niche analysis results

- samp_type:

  Character, sample type - "SS" (single sample) or "MS" (multi-sample)
  (default: "SS")

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- samp_grp_index:

  Logical, whether to group by sample groups in MS mode (default: FALSE)

- meta_key:

  Character, metadata key for MS mode when niche_key is NULL

- niche_key:

  Character, niche key to analyze (only one value supported)

- group_by:

  Character, column name for cell type grouping

- col:

  Color palette for visualization (default:
  COLOR_LIST\$PALETTE_WHITE_BG)

- return_data:

  Logical, whether to return the plot data (default: FALSE)

## Value

If return_data = TRUE, returns a list of plots per sample; otherwise
NULL

## Examples

``` r
if (FALSE) { # \dontrun{
# Single-sample niche composition
CalNicheComposition(
  IST_obj = ist_obj,
  samp_type = "SS",
  niche_key = "niche_virulence",
  group_by = "cell_type"
)

# Multi-sample niche composition with sample grouping
CalNicheComposition(
  IST_obj = ist_obj,
  samp_type = "MS",
  loop_id = "LoopAllMulti",
  niche_key = "niche_virulence",
  samp_grp_index = TRUE,
  group_by = "cell_type"
)
} # }
```
