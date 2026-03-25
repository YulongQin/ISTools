# Add annotations to SingleSampNiche genes

Adds new columns to the gene data in SingleSampNiche objects.

## Usage

``` r
AddSSNicheGenes(IST_obj, ...)

# S3 method for class 'IST'
AddSSNicheGenes(
  IST_obj = NULL,
  gene = NULL,
  label = NULL,
  add_colnm = NULL,
  loop_id = "LoopAllSamp",
  niche_key = NULL,
  ...
)
```

## Arguments

- IST_obj:

  An IST object

- ...:

  Additional arguments passed to methods

- gene:

  Character vector, gene names to annotate

- label:

  Character vector, labels to assign to genes (optional)

- add_colnm:

  Character, name of new column to add

- loop_id:

  Character vector, sample identifiers

- niche_key:

  Character, niche key to modify

## Value

Modified IST object with updated niche gene data

## Examples

``` r
if (FALSE) { # \dontrun{
# Add gene family annotations
IST_obj <- AddSSNicheGenes(ist_obj,
                           gene = c("gene1", "gene2", "gene3"),
                           label = c("toxin", "toxin", "adhesin"),
                           add_colnm = "gene_family",
                           niche_key = "niche_virulence")
} # }
```
