# Summarise `malaria-profiler` outputs

Summarise `malaria-profiler` outputs

## Usage

``` r
summarise_mp_outputs(
  target_dir,
  variants_file,
  metadata_file,
  mutation_type = c("confirmed", "potential"),
  geo_localisation = TRUE
)
```

## Arguments

- target_dir:

  A character with the full path to the folder where the output files
  (.json files) from malaria-profiler are stored

- variants_file:

  A character with the path to the file with the confirmed mutations
  associated to the resistance to the drugs of interest

- metadata_file:

  A character with the full path to file with the sample metadata

- mutation_type:

  A character that determines whether to summarise the outputs for the
  confirmed or potential drug resistance mutations. Possible values are
  `confirmed` or `potential`. Default is `confirmed`

- geo_localisation:

  A logical used to determine whether to summarise the outputs about the
  potential region of origin of the samples. If `TRUE`, a column with
  inferred sample potential regions of origin will be added to the
  sample metadata table. Default is `TRUE`

## Value

A list with the following three elements:

1.  variants: a data frame with the genomic coordinates of the detected
    drug resistance mutations and `x` columns of `0` and `1` used to
    determine whether the sample carries a particular mutation (`1`) or
    no (`0`), where `x` is the number of processed samples.

2.  metadata: a data frame with the metadata table augmented with two
    columns: the sample origins, and the fraction of detected SNP
    markers used for geographic origin prediction

3.  region_probabilities: a data frame with the probabilities obtained
    during the geographic origin prediction. This is empty when
    `geo_localisation = FALSE`

## Examples

``` r
if (FALSE) { # \dontrun{
  # establish the connection to the database
  results <- summarise_mp_outputs(
    target_dir = tempdir(),
    variants_file = file.path(tempdir(), "variant.csv"),
    metadata_file = file.path(tempdir(), "metadata.txt"),
    mutation_type = "confirmed",
    geo_localisation = "TRUE"
  )
} # }
```
