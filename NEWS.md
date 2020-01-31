# 1.4.1 (18 November, 2019)

* List available gene expression and drug sensitivity associations
(`listExpressionDrugSensitivityAssociation()`)

## Minor changes

* Copy-edit CMap-related console messages
* Improve tutorial
* Copy-edit function documentation

# 1.4 (25 October, 2019)

## New features

* Predict targeting drugs (`predictTargetingDrug()`):
    - Based on expression and drug sensitivity associations derived from NCI60,
    CTRP and GDSC data (see `loadExpressionDrugSensitivityAssociation()`)
    - Compare user-provided differential expression profile with gene expression
    and drug sensitivity associations to predict targeting drugs and their
    targeted genes
    - Compounds are ranked based on their relative targeting potential
    - Plot candidate targeting drugs against ranked compound perturbations using
    `plotTargetingDrugsVSsimilarPerturbations()`, highlighting compounds that
    selectively select against cells with a similar differential gene expression
    profile
* Analyse drug set enrichment (`performDSEA()`):
    - Prepare drug sets based on a table with compound identifiers and 
    respective 2D and 3D molecular descriptors using `prepareDrugSets()`
    - Test drug set enrichment on results from
    `rankSimilarPerturbations()` (when ranking against compound perturbations)
    and `predictTargetingDrugs()`
* Convert ENSEMBL identifiers to gene symbols using
`convertENSEMBLtoGeneSymbols()`

## Major changes

* Update the tutorial and function documentation
* Remove most `L1000` instances, including in function names:
    - `getL1000perturbationTypes()` -> `getCMapPerturbationTypes()`
    - `getL1000conditions()`        -> `getCMapConditions()`
    - `downloadL1000data()`         -> `loadCMapData()`
    - `filterL1000metadata()`       -> `filterCMapMetadata()`
    - `loadL1000perturbations()`    -> `prepareCMapPerturbations()`
    - `compareAgainstL1000()`       -> `rankSimilarPerturbations()`
    - `plotL1000comparison()`       -> `plot()`
* Improve loading of ENCODE samples (`loadENCODEsamples()`):
    - Rename function from `downloadENCODEsamples()` to `loadENCODEsamples()`
    - Load ENCODE samples regarding multiple cell lines and experiment targets
    using `loadENCODEsamples()`
* Improve CMap data and metadata retrieval:
    - By default, do not return control perturbation types when using
    `getCMapPerturbationTypes()` (unless if using argument `control = TRUE`)
    - Parse CMap identifiers using `parseCMapID()`
    - Load CMap's compound metadata using `loadCMapData()`
    - Ask to download CMap perturbations z-scores file for differential 
    expression if not found (avoiding downloading a huge file without user 
    consent)
* Improve preparation of CMap perturbations (`prepareCMapPerturbations()`):
    - Allow to load CMap metadata directly from files when using file paths as
    arguments of `prepareCMapPerturbations()`
    - Significantly decrease memory required to use cTRAP by loading chunks of
    z-scores from CMap perturbations on-demand (a slight decrease in time
    performance is expected), unless `prepareCMapPerturbations()` is run with
    argument `loadZscores = TRUE`
    - Display summary of loaded perturbations after running 
    `prepareCMapPerturbations()`
* Improve ranking of similar perturbations (`rankSimilarPerturbation()`):
    - Redesigned output: long (instead of wide) table
    - By default, calculate mean across cell lines if there is more than one 
    cell line available; disabled if argument `cellLineMean = FALSE`
    - Allow to rank (or not) individual cell line perturbations (argument
    `rankIndividualCellLinePerturbations`) when the mean is calculated
    - Allow to perform multiple comparison methods if desired (by providing a 
    vector of supported methods via the `method` argument)
    - Calculate the rank product's rank to assess ranks across multiple methods
    - Sort results based on rank product's rank (or the rank of the only
    comparison method performed, otherwise)
    - Include information for calculated means across cell lines in metadata
    - Include run time as an attribute
* Improve metadata display for a `similarPerturbations` object, obtained after
running `rankSimilarPerturbations()`:
    - Show further metadata information (including compound data, if available)
    related with a given perturbation by calling `print()` with a
    `similarPerturbations` object and a specific perturbation identifier
    - Show a complete table with metadata (and compound information, if 
    available) when calling `as.table()` with a `similarPerturbations` object
* Improve plotting (`plot()`):
    - Plot comparison results against all compared data by calling `plot()` with
    the results obtained after running `rankSimilarPerturbations()` or
    `predictTargetingDrugs()`; non-ranked compared data can also be plotted with
    argument `plotNonRankedPerturbations = TRUE`
    - Render scatter and Gene Set Enrichment Analysis (GSEA) plots between
    differential expression results and a single perturbation by calling 
    `plot()` with a `perturbationChanges` object (if an identifier regarding the 
    summary of multiple perturbations scores across cell lines is given, the
    plots are coloured by cell line)
    - When displaying GSEA plots, plot results for most up- and down-regulated
    user-provided differentially expressed genes (by default)
    - Improve GSEA plot style, including rug plot in enrichment score plot
    (replacing the gene hit plot)

## Bug fixes and minor changes

* CMap metadata minor improvements:
    - Improve list returned by `getCMapConditions()`, including sorting of dose 
    and time points
    - Correctly set instances of `-666` in CMap metadata as missing values and 
    fix specific issues with metadata (such as doses displayed as
    `300 ng|300 ng`)
    - In compound metadata, fix missing values showing as literal "NA" values
* CMap perturbation minor improvements:
    - Fix error when subsetting a `perturbationChanges` object with only one row
    - Improve performance when subsetting `perturbationChanges` objects
* Minor improvements to `rankSimilarPerturbations()`:
    - Correctly set name of perturbations depending on their type (genes,
    biological agents or compounds)
    - Improve performance when correlating against multiple cell lines
    - Remove `cellLine` argument (please filter conditions with upstream
    functions such as `filterCMapMetadata()`)
    - Fix incorrect label of first column identifiers
    - Report run time and settings used
    - Perform comparisons against perturbations disregarding their cell lines
    (faster runtime)
    - Fix error when trying to calculate the mean for cell lines with no 
    intersecting conditions available
    - Clearly state to the user when no intersecting genes were found between
    input dataset and CMap data
* Minor improvements to `plot()`:
    - Improve rendering performance of the GSEA plot
    - Fix disproportionate height between top and bottom enrichment score panels
    in GSEA plots
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
