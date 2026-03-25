# Plot Niche Cell Communication Results

Visualizes CellChat analysis results with various plot types including
circle plots, heatmaps, bubble plots, and spatial signaling plots.

## Usage

``` r
Plot_NicheCellComm(
  IST_obj = NULL,
  CellComm_data = NULL,
  samp_type = "SS",
  loop_id = "LoopAllSamp",
  sources.use = NULL,
  targets.use = NULL,
  signaling = NULL,
  pairLR.use = NULL,
  col = NULL
)
```

## Arguments

- IST_obj:

  An IST object containing niche analysis results

- CellComm_data:

  List, CellChat results from CalNicheCellComm

- samp_type:

  Character, sample type - "SS" or "MS" (default: "SS")

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- sources.use:

  Character vector, source cell types to include

- targets.use:

  Character vector, target cell types to include

- signaling:

  Character vector, signaling pathways to plot

- pairLR.use:

  Data frame, specific ligand-receptor pairs to plot

- col:

  Named vector, colors for cell types

## Value

NULL (invisible), generates plots

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot CellChat results
Plot_NicheCellComm(
  IST_obj = ist_obj,
  CellComm_data = cellcomm_results,
  signaling = c("TGFb", "WNT"),
  sources.use = "Epithelial",
  targets.use = "Immune"
)
} # }
```
