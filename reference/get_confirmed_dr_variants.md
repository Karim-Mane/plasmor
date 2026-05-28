# Summarise the malaria-profiler outputs for confirmed drug resistance mutations

Summarise the malaria-profiler outputs for confirmed drug resistance
mutations

## Usage

``` r
get_confirmed_dr_variants(json, variants, sample_name)
```

## Arguments

- json:

  A list object obtained from reading an output file in JSon format

- variants:

  A data frame that will be used to store the variant profile for each
  sample

- sample_name:

  A character with a sample name

## Value

An updated variant table where the profile of the specified sample is
added
