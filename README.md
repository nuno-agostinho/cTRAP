# cTRAP

<!-- badges: start -->
[![GitHub Actions Build status][ghActionsIcon]][ghActions]
[![Coverage status][codecovBadge]][codecov]
<!-- badges: end -->

`cTRAP` is an R package designed to compare differential gene
expression results with those from known cellular perturbations (such as gene 
knock-down, overexpression or small molecules) derived from the 
[Connectivity Map][clue.io] [(Subramanian et al., Cell 2017)][subramanian2017].
Such analyses allow not only to infer the molecular causes of the observed 
difference in gene expression but also to identify small molecules that could 
drive or revert specific transcriptomic alterations.

## Installing

[cTRAP is available in Bioconductor][bioconductor] and can be installed with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("cTRAP")
```

If you prefer, you can install cTRAP from GitHub instead:

``` r
install.packages("remotes")
remotes::install_github("nuno-agostinho/cTRAP")
```

[clue.io]: https://clue.io/
[subramanian2017]: https://doi.org/10.1016/j.cell.2017.10.049
[codecov]: https://codecov.io/github/nuno-agostinho/cTRAP?branch=master
[codecovBadge]: https://img.shields.io/codecov/c/github/nuno-agostinho/cTRAP/master.svg
[ghActions]: https://github.com/nuno-agostinho/cTRAP/actions
[ghActionsIcon]: https://github.com/nuno-agostinho/cTRAP/workflows/R-CMD-check/badge.svg
[bioconductor]: http://bioconductor.org/packages/cTRAP
