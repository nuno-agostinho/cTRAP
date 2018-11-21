# 1.0.3 (12 November, 2018)

* Include copyright text and full license for source code distributed from cmapR

# 1.0.2 (11 November, 2018)

* Fix comparison against CMap perturbations using gene set enrichment analysis 
(the resulting score was the additive inverse of the real scores)

# 1.0.1 (2 November, 2018)

* Update title, author names, version and README
* Remove biomaRt dependency
* By default, `getL1000conditions` now shows CMap perturbation types except for
controls
* Compare against CMap perturbations (`compareAgainstL1000` function):
    - Remove "_t" from resulting column names (as the t-statistic may or may not
    be used)
    - Select p-value adjustment method when performing correlation analyses
    (Benjamini-Hochberg is set by default)
* Documentation:
    - Fix obsolete function calls in function documentation
    - Hide non-exported functions from reference PDF manual
