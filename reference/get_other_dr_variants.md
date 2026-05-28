# Get the sample profiles for the potential drug resistance mutations

Get the sample profiles for the potential drug resistance mutations

## Usage

``` r
get_other_dr_variants(nested_list, sample_name, tmp_others_dr)
```

## Arguments

- nested_list:

  A list with the information on the potential drug resistance mutations
  found on the given sample

- sample_name:

  A character with a sample name

- tmp_others_dr:

  A data frame used to store the sample profiles

## Value

A data frame with the sample profiles augmented with a new column
corresponding to the profile of the specified sample
