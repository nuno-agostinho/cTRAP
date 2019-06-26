# 1.2.1 (10 June, 2019)

* Update the tutorial and function documentation
* Replace all `L1000` references with `CMap`, as appropriate, including in
function names:
    - `getL1000perturbationTypes()` -> `getCMapPerturbationTypes()`
    - `getL1000conditions()`        -> `getCMapConditions()`
    - `downloadL1000data()`         -> `loadCMapData()`
    - `filterL1000metadata()`       -> `filterCMapMetadata()`
    - `loadL1000perturbations()`    -> `prepareCMapPerturbations()`
    - `compareAgainstL1000()`       -> `compareAgainstCMap()`
    - `plotL1000comparison()`       -> `plot()`
* Improve loading of ENCODE samples (`loadENCODEsamples()`):
    - Rename function from `downloadENCODEsamples()` to `loadENCODEsamples()`
    - Allow to load ENCODE samples regarding multiple cell lines and experiment 
    targets using `loadENCODEsamples()`
* Improve CMap data and metadata retrieval:
    - By default, do not return control pertubation types when using
    `getCMapPerturbationTypes()` (unless if using argument `control = TRUE`)
    - Allow to parse CMap identifiers using `parseCMapID()`
    - Allow to load CMap's compound metadata using `loadCMapData()`
    - Ask to download CMap perturbations z-scores file for differential 
    expression if not found (avoiding downloading a huge file without user 
    consent)
* Improve preparation of CMap perturbations (`prepareCMapPerturbations()`):
    - Allow to load CMap metadata directly from files when using filepaths as
    arguments of `prepareCMapPerturbations()`
    - Significantly decrease memory required to use cTRAP by loading chunks of
    z-scores from CMap perturbations on-demand (a slight decrease in time
    performance is expected), unless `prepareCMapPerturbations()` is run with
    argument `loadZscores = TRUE`
    - Display summary of loaded perturbations after running 
    `prepareCMapPerturbations()`
* Improve comparison against CMap perturbations (`compareAgainstCMap()`):
    - Redesigned output: long (instead of wide) table
    - By default, calculate mean across cell lines if there is more than one 
    cell line available; disabled if argument `cellLineMean = FALSE`
    - Allow to rank (or not) individual cell line perturbations (argument
    `rankIndividualCellLinePerturbations`) when the mean is calculated
    - Allow to perform multiple comparison methods if desired (by providing a 
    vector of supported methods via the `method` argument)
    - Calculate the rank product's rank to assess ranks across multiple methods
* Improve metadata display for perturbations from CMap comparisons 
(`cmapComparison` object, obtained after running `compareAgainstCMap()`):
    - Show further metadata information (including compound data, if available)
    related with a given perturbation by calling `print()` with a
    `cmapComparison` object and a specific pertubation
    - Show a complete table with metadata (and compound information, if 
    available) when calling `as.table()` with a `cmapComparison` object
* Improve plotting of CMap comparisons (`plot()`):
    - Plot comparison results against all perturbations by calling `plot()` with
    a `cmapComparison` object (non-ranked perturbations may also be plotted 
    with `plotNonRankedPerturbations = TRUE`)
    - Plot scatter and Gene Set Enrichment Analysis (GSEA) plots between 
    differential expression results and a single perturbation by calling 
    `plot()` with a `cmapPerturbations` object (if an identifier regarding the
    summary of multiple perturbations scores across cell lines is given, the
    plots are coloured by cell line)
    - When displaying GSEA plots, automatically render results for both top and 
    bottom genes by default

## Bug fixes and minor changes

* CMap metadata minor improvements:
    - Improve list returned by `getCMapConditions()`, including sorting of dose 
    and time points
    - Correctly set instances of `-666` in CMap metadata as missing values and 
    fix specific issues with metadata (such as doses displayed as
    `300 ng|300 ng`)
* CMap perturbation minor improvements:
    - Fix error when subsetting a `cmapPerturbations` object with only one row
    - Improve performance when subsetting `cmapPerturbations` objects
* Minor improvements to `compareAgainstCMap()`:
    - Correctly set name of perturbagens depending on the perturbation type
    (genes, biological agents or compounds)
    - Improve performance when correlating against multiple cell lines
    - Remove `cellLine` argument (please filter conditions with upstream
    functions such as `filterCMapMetadata()`)
    - Fix incorrect label of first column identifiers
    - Report run time and settings used
    - Perform comparisons against perturbations disregarding their cell lines
    (faster runtime)
* Minor improvements to `plot()`:
    - Improve rendering performance of the GSEA plot
    - Fix disproportion between top and bottom enrichment score panels in GSEA
    plots
* Update demo datasets:
    - Update the `cmapPerturbationsCompounds` and `cmapPerturbationsKD` datasets 
    according to new internal changes and fix their respective code in the 
    documentation
* Include license and copyright text for cmapR code

# 1.0.3 (3 December, 2018)

* Add tag ImmunoOncology to BiocViews

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
