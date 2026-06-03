
<!-- README.md is generated from README.Rmd. Please edit that file. -->

<!-- The code to render this README is stored in .github/workflows/render-readme.yaml -->

<!-- Variables marked with double curly braces will be transformed beforehand: -->

<!-- - `packagename` is extracted from the DESCRIPTION file -->

<!-- - `gh_repo` is extracted via a special environment variable in GitHub Actions -->

<!-- Having these variables allow us to re-use this the same README template across repos -->

<!-- without ever hardcoding repo-specific elements. -->

# plasmor <img src="man/figures/logo.svg" align="right" width="120" alt="" />

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/license/mit/)
[![R-CMD-check](https://github.com/Karim-Mane/plasmor/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Karim-Mane/plasmor/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/Karim-Mane/plasmor/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Karim-Mane/plasmor?branch=main)
[![lifecycle-concept](https://raw.githubusercontent.com/reconverse/reconverse.github.io/master/images/badge-concept.svg)](https://www.reconverse.org/lifecycle.html#concept)
<!-- badges: end -->

plasmor provides functions to analyse different types of data generated
by a researcher working in the field of parasitology, specially on the
malaria parasite *P. falciparum*.

<!-- This sentence is optional and can be removed -->

plasmor is developed at [IHU Méditerranée
Infection](https://www.mediterranee-infection.com/accueil-2/) as part of
Karim’s PhD at [Aix Marseille Université](https://www.univ-amu.fr/en).

## Installation

You can install the development version of plasmor from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("Karim-Mane/plasmor")
```

## Example

These examples illustrate some of the current functionalities

## Development

### Lifecycle

This package is currently a *concept*, as defined by the [RECON software
lifecycle](https://www.reconverse.org/lifecycle.html). This means that
essential features and mechanisms are still being developed, and the
package is not ready for use outside of the development team.

### Contributions

Contributions are welcome via [pull
requests](https://github.com/Karim-Mane/plasmor/pulls).

### Code of Conduct

Please note that the plasmor project is released with a [Contributor
Code of
Conduct](https://github.com/epiverse-trace/.github/blob/main/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

## Citing this package

``` r
citation("plasmor")
#> To cite package 'plasmor' in publications use:
#> 
#>   Mané K, Delandre O, javelle E (2026). _plasmor:
#>   Analyse Plasmodium falciparum Data_. R package version
#>   0.0.1.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {plasmor: Analyse Plasmodium falciparum Data},
#>     author = {Karim Mané and Océane Delandre and Emilie Javelle},
#>     year = {2026},
#>     note = {R package version 0.0.1},
#>   }
```
