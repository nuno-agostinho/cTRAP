# cTRAP 1.10.1 (7 June, 2021)

* `convertGeneIdentifiers()` replaces `convertENSEMBLtoGeneSymbols()`:
    - Use AnnotationHub to convert to gene symbols (instead of biomaRt that has
    been unstable) and allow to 

# cTRAP 1.10.0 (18 March, 2021)

## Improvements to graphical interface functions:

* New `launchDrugSetEnrichmentAnalysis()` function to analyse drug set
enrichment and visualize respective results
* `launchCMapDataLoader()`:
    - Now allows to load multiple CMap perturbation types simultaneously
    - Keep selected timepoint, dosage and cell line options when selecting
    another perturbation type
    - Add bubble plot of CMap perturbation types
* `launchResultPlotter()`:
    - Now allows to view tables below specific plots and drag-and-select those
    plots to filter data in those same tables
    - When plotting targeting drugs and similar perturbations, update available
    columns and correctly use user-selected column to plot
* `launchMetadataViewer()` now correctly parses values from `Input`
attributes as numeric

## Major changes

* `prepareCMapPerturbations()`: directly set perturbation type, cell line,
timepoint and dosage conditions as arguments
* `rankSimilarPerturbations()` and `predictTargetingDrugs()`:
    - Avoid redundant loading of data chunks, slightly decreasing run time
    - Lower memory footprint when using NCI60's gene expression and drug
    sensitivity association (now available in HDF5 files) by loading and
    processing data in chunks
    - Faster GSEA-based score calculation (up to 4-7 times faster)
    - New `threads` argument allows to set number of parallel threads (not
    supported on Windows)
    - New `chunkGiB` argument allows to set size of data chunks when reading
    from supported HDF5 files (decreases peak RAM usage)
    - New `verbose` argument allows to increase details printed in the console
* `prepareDrugSets()`: allow greater control on the creation of bins based on
numeric columns, including the setting of maximum number of bins per column and
minimum bin size
* `analyseDrugSetEnrichment()` and `plotDrugSetEnrichment()`: allow to select
columns to use when comparing compound identifiers between datasets

## Bug fixes and minor changes

* `filterCMapMetadata()`: allow filtering CMap metadata based on multiple
perturbation types
* `prepareDrugSets()`: fix issues with 3D descriptors containing missing values
* `plot()`:
    - Fix wrong labels when plotting `targetingDrugs` objects
    - Avoid printing "NA" in labels identifying metadata for perturbations
* `plotTargetingDrugsVSsimilarPerturbations()`:
    - Fix highlighting of plot points depending whether drug activity is
    directly proportional to drug sensitivity
    - Include rug plot
* When subsetting a `perturbationChanges` or an
`expressionDrugSensitivityAssociation` object, passing only one argument
extracts its columns as in previous versions of cTRAP (similarly to when
subsetting a `data.frame`)
* `analyseDrugSetEnrichment()`: for the resulting table, the name of the first
column was renamed from `pathway` to `descriptor`

# cTRAP 1.8 (23 October, 2020)

## Interactive functions for loading data and analysing results

* New Shiny-based graphical interface functions:
    - `launchDiffExprLoader()`: load differential expression data
    - `launchCMapDataLoader()`: load CMap data
    - `launchResultPlotter()`: view and plot data results
    - `launchMetadataViewer()`: check metadata of a given object

## Major changes

* `downloadENCODEknockdownMetadata()`: metadata is automatically saved to a file
in order to avoid downloading metadata every time this function is run
* `plotTargetingDrugsVSsimilarPerturbations()`:
    - automatically look for matching compounds in multiple columns of both
    datasets
    - allow to manually select columns on which to merge datasets
* `prepareDrugSets()`: drug sets based on numeric molecular descriptors are now
prepared using evenly-distributed intervals
* Simplify tutorial

# cTRAP 1.6.1 (17 August, 2020)

* `listExpressionDrugSensitivityAssociation()` lists available gene expression
and drug sensitivity associations
* First argument of `rankSimilarPerturbations()` and `predictTargetingDrugs()`
changed name from `diffExprGenes` to `input` and now accepts:
    - `Named numeric vector` containing differential gene expression values
    with gene symbols as names, as before;
    - `Character vector` containing a custom gene set to test for enrichment
    (only to use with GSEA).
* In `rankSimilarPerturbations()` and `predictTargetingDrugs()`, when performing
`gsea` method, allow to set different gene set size for top up- and
down-regulated genes with `geneSize` argument:
    - e.g. `geneSize=c(100, 200)` creates gene sets from the top 100 up-
    and top 200 down-regulated genes
    - using `geneSize=c(150, 150)` or `geneSize=150` is equivalent
* Plotting:
    - `plot()` now supports plotting `predictTargetingDrugs()` results for a
    given drug, e.g. `plot(targetingDrugs, "1425")`
    - `plot()` nows allows to set plot title with argument `title`
    - `plot()` now plots results based on available methods instead of trying
    to plot based on results from `spearman` method only
    - GSEA plots now support two or less gene hits
    - GSEA plots now support plotting of multiple perturbations
    - GESA plots now show the first and last values of ranked genes
    - `plotDrugSetEnrichment()` now returns a list whose names are drug set
    names
* `as.table()` improvements:
    - Return cell identifiers and gene information (if available and as needed)
    - Support `predictTargetingDrugs()` results
    - Return results ordered as found on input

## Bug fixes and minor changes

* `downloadENCODEknockdownMetadata()` now correctly retrieves metadata following
a change in the metadata content from ENCODE
* Fix bugs when rendering GSEA plots due to deprecated functions in `ggplot2`
* Improve tutorial
* Copy-edit CMap-related console messages
* Copy-edit function documentation

# cTRAP 1.4 (25 October, 2019)

## New features

* Predict targeting drugs (`predictTargetingDrugs()`):
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
* Analyse drug set enrichment (`analyseDrugSetEnrichment()`):
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
* Improve ranking of similar perturbations (`rankSimilarPerturbations()`):
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
* Include license and copyright text for `cmapR` code

# cTRAP 1.0.3 (3 December, 2018)

* Add tag ImmunoOncology to BiocViews

# cTRAP 1.0.2 (11 November, 2018)

* Fix comparison against CMap perturbations using gene set enrichment analysis 
(the resulting score was the additive inverse of the real scores)

# cTRAP 1.0.1 (2 November, 2018)

* Update title, author names, version and README
* Remove biomaRt dependency
* By default, `getL1000conditions()` now shows CMap perturbation types except 
for controls
* Compare against CMap perturbations (`compareAgainstL1000()`):
    - Remove "_t" from resulting column names (as the t-statistic may or may not
    be used)
    - Select p-value adjustment method when performing correlation analyses
    (Benjamini-Hochberg is set by default)
* Documentation:
    - Fix obsolete function calls in function documentation
    - Hide non-exported functions from reference PDF manual
