# l1000

`l1000` is an R package designed to compare differential gene
expression results with those from known cellular perturbations (such as gene 
knock-down, overexpression or small molecules) derived from the 
[Connectivity Map][clue.io] (Subramanian et al., Cell 2017). Such analyses allow
not only to infer the molecular causes of the observed difference in gene
expression but also to identify small molecules that could drive or revert
specific transcriptomic alterations.

## Installation

You can install l1000 from GitHub with:

``` r
install.packages("devtools")
devtools::install_github("nuno-agostinho/l1000")
```

[clue.io]: https://clue.io/
