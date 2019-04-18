# 1.0.3 (18 April, 2019)

* Update and copy-edit the tutorial and function documentation
* Replace all references to `L1000` with `CMap`, as appropriate, including in
dataset and function names:
    - `getL1000perturbationTypes()` -> `getCMapPerturbationTypes()`
    - `getL1000conditions()`        -> `getCMapConditions()`
    - `downloadL1000data()`         -> `loadCMapData()`
    - `filterL1000metadata()`       -> `filterCMapMetadata()`
    - `loadL1000perturbations()`    -> `loadCMapPerturbations()`
    - `compareAgainstL1000()`       -> `compareAgainstCMap()`
    - `plotL1000comparison()`       -> `plot()`
* Rename `downloadENCODEsamples()` to `loadENCODEsamples()`
* Integrate compound metadata from CMap:
    - Allow to download CMap compound metadata using `loadCMapData()`
* Improve comparison against CMap perturbations (`compareAgainstCMap()`):
    - Correctly set name of perturbagens depending on the perturbation type
    (genes, biological agents or compounds)
    - Improve performance when correlating against multiple cell lines
    - Remove `cellLine` argument
    - By default, calculate mean across cell lines if there is more than one 
    cell line available; mean calculation can now also be avoided if the 
    argument `cellLineMean = FALSE`
    - Automatically order results based on the correlation coefficient (for 
    Spearman's and Pearson's correlation) or the weighted connectivity score 
    (WTCS) score (for Gene Set Enrichment Analysis) using the mean across cell
    lines or the results for the first cell line alone
    - Improve comparison performance when testing for correlation
    - Fix incorrect label of identifiers
    - Report run time and settings used
* Improve plotting of L1000 comparisons (`plot`):
    - Allow to call `plot()` with a `cmapComparison` object
    - When displaying Gene Set Enrichment Analysis (GSEA) plots, automatically
    render results for both top and bottom genes 79by default
    - Improve rendering performance of the GSEA plot
    - Plot comparison results against all perturbations
* Update demo datasets:
    - Update the `cmapPerturbationsSmallMolecules` and 
    `cmapPerturbationsKnockdown` datasets according to new internal changes and
    fix their respective code in the documentation
* Replace instances of -666 in CMap data to show up as missing values
* Include copyright text and full license for source code distributed from cmapR
* Fix error when subsetting a `cmapPerturbations` object with only one row

# 1.0.2 (11 November, 2018)

* Fix comparison against CMap perturbations using gene set enrichment analysis 
(the resulting score was the additive inverse of the real scores)

# 1.0.1 (2 November, 2018)

* Update title, author names, version and README
* Remove biomaRt dependency
* By default, `getL1000conditions()` now shows CMap perturbation types except 
for controls
* Compare against CMap perturbations (`compareAgainstL1000()` function):
    - Remove "_t" from resulting column names (as the t-statistic may or may not
    be used)
    - Select p-value adjustment method when performing correlation analyses
    (Benjamini-Hochberg is set by default)
* Documentation:
    - Fix obsolete function calls in function documentation
    - Hide non-exported functions from reference PDF manual
