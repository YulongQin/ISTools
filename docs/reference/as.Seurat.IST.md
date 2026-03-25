# Convert IST object to Seurat object

Converts an IST object back to a standard Seurat object, preserving all
assay data and basic metadata but removing IST-specific analysis
structures.

## Usage

``` r
as.Seurat.IST(IST_obj)
```

## Arguments

- IST_obj:

  An IST object to be converted

## Value

A Seurat object containing the same assay data and basic metadata

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert IST object back to Seurat
seurat_obj <- as.Seurat.IST(ist_obj)
} # }
```
