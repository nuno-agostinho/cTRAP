# 1.0.3 (5 December, 2018)

* Replace instances of -666 in CMap data to show up as missing values
* Include copyright text and full license for source code distributed from cmapR
* Compare against CMap perturbations (`compareAgainstL1000()` function):
    - Correctly set name of perturbagens depending on the perturbation type
    (genes, biological agents or compounds)
    - By default, calculate mean across cell lines if there is more than one 
    cell line available; mean calculation can now also be avoided if the 
    argument `cellLineMean = FALSE`
    - Automatically order results based on the correlation coefficient (for 
    Spearman's and Pearson's correlation) or the weighted connectivity score 
    (WTCS) score (for Gene Set Enrichment Analysis) using the mean across cell
    lines or the results for the first cell line alone
    - Improve comparison performance when testing for correlation
* Included datasets:
    - Update the `l1000perturbationsSmallMolecules` and 
    `l1000perturbationsKnockdown` datasets according to new internal changes and
    fix their respective code in the documentation

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
