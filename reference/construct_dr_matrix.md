# Add a drug resistance profile for the specified sample to the data frame of drug resistance profiles based on potential drug resistance mutations

Add a drug resistance profile for the specified sample to the data frame
of drug resistance profiles based on potential drug resistance mutations

## Usage

``` r
construct_dr_matrix(dat, sample_name, tmp_others_dr)
```

## Arguments

- dat:

  A data frame with the resistance profile of the target sample based on
  the potential drug resistance mutations

- sample_name:

  A character with a sample name

- tmp_others_dr:

  A data frame used to store the sample profiles

## Value

A data frame potential drug resistance profile augmented with one column
corresponding the profile of the target sample
