# Convert Seurat rds file to AnnData h5ad file

Converts a Seurat rds file to a Python AnnData h5ad file for use with
scanpy and other Python-based tools.

## Usage

``` r
rds2annh5ad(
  seurat_obj = NULL,
  python_path = "D:/APP/anaconda3/envs/scanpy/python.exe",
  data_type = "stRNA",
  convert_mode = "seurat",
  assay_id = NULL,
  X_index = "rawX",
  binsize = 1,
  SCT_index = FALSE,
  reduction_index = FALSE,
  image_index = NULL,
  return_object = FALSE,
  grp_nm = NULL,
  dir_nm = "M0_rds2annh5ad"
)
```

## Arguments

- seurat_obj:

  A Seurat object to convert

- python_path:

  Character, path to Python executable with scanpy installed (default:
  "D:/APP/anaconda3/envs/scanpy/python.exe")

- data_type:

  Character, data type - "stRNA" (spatial) or "scRNA" (single-cell)
  (default: "stRNA")

- convert_mode:

  Character, conversion method - "seurat" (using SeuratDisk) (default:
  "seurat")

- assay_id:

  Character, assay name to convert (default: NULL, auto-detected)

- X_index:

  Character, which matrix to use as counts - "X" or "rawX" (default:
  "rawX")

- binsize:

  Numeric, bin size for spatial coordinate scaling (default: 1)

- SCT_index:

  Logical, whether to include SCT-transformed data (default: FALSE)

- reduction_index:

  Logical, whether to include dimensional reductions (default: FALSE)

- image_index:

  Logical, whether to include spatial image data (default: NULL,
  auto-detected)

- return_object:

  Logical, whether to return the object (default: FALSE)

- grp_nm:

  Character, group name for output organization (default: NULL)

- dir_nm:

  Character, directory name for output (default: "M0_rds2annh5ad")

## Value

NULL (invisible), saves h5ad file to disk

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert Seurat object to h5ad
rds2annh5ad(
  seurat_obj = seurat_object,
  python_path = "~/anaconda3/envs/scanpy/bin/python",
  data_type = "stRNA",
  grp_nm = "sample1"
)
} # }
```
