# Extract the prediction region of origin for the specified sample

Extract the prediction region of origin for the specified sample

## Usage

``` r
get_inferred_region(json, metadata, region_proba, sample_name)
```

## Arguments

- json:

  A list object obtained from reading an output file in JSon format

- region_proba:

  A data frame used to store the prediction probabilities for each
  region

- sample_name:

  A character with a sample name

## Value

A list with the updated sample metadata and region probabilities data
frame
