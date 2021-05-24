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

### Bioconductor

[cTRAP is available in Bioconductor][bioconductor] and can be installed with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("cTRAP")
```

### GitHub

cTRAP can also be installed from GitHub instead:

``` r
install.packages("remotes")
remotes::install_github("nuno-agostinho/cTRAP")
```

### Docker

The Docker images are based on [Bioconductor Docker][biocDocker] and contain cTRAP and its dependencies.

1. Pull the latest Docker image:

```
docker pull ghcr.io/nuno-agostinho/ctrap:latest
```

2. Start RStudio Web from the Docker image:

```
docker run -e PASSWORD=bioc -p 8787:8787 ghcr.io/nuno-agostinho/ctrap:latest
```

3. Go to RStudio Web via your web browser at https://localhost:8787
4. Login in RStudio Web with user `rstudio` and password `bioc`
5. Load the package in RStudio Web using `library(cTRAP)`

[clue.io]: https://clue.io/
[subramanian2017]: https://doi.org/10.1016/j.cell.2017.10.049
[codecov]: https://codecov.io/github/nuno-agostinho/cTRAP?branch=master
[codecovBadge]: https://img.shields.io/codecov/c/github/nuno-agostinho/cTRAP/master.svg
[ghActions]: https://github.com/nuno-agostinho/cTRAP/actions
[ghActionsIcon]: https://github.com/nuno-agostinho/cTRAP/workflows/R-CMD-check-bioc/badge.svg
[bioconductor]: http://bioconductor.org/packages/cTRAP
[biocDocker]: https://github.com/Bioconductor/bioconductor_docker
