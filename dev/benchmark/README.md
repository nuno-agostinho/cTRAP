# Time and memory profiling of cTRAP

> Nuno Agostinho, 27 November 2020

cTRAP is a multi-threaded R package composed of three modules. These scripts
test the critical module of ranking user-provided differential expression
results against differential expression results from CMap perturbations.

They also benchmark the prediction of targeting drugs (using the NCI60 gene
expression and drug sensitivity association, the most time-consuming option)
and drug set enrichment analysis.

## cTRAP performance milestones (dev versions)

* 1.8.0: release version for reference
* 1.8.1 (b13ee45): faster GSEA-based score calculation
* 1.8.1 (c34566c): avoid redundant loading of chunks from CMap perturbation data
* 1.8.1 (3e3720d): multi-thread support in systems that support forking (e.g.
Linux and macOS, but not Windows) and print times for measurable actions (to
directly compare with memory profile)
* 1.8.1 (9b96229): fix issues with missing values when preparing drug
descriptor sets
* 1.8.1 (9852a1a): improve drug set enrichment analysis (fix bugs and allow to
match compounds as done by `plotTargetingDrugsVSsimilarPerturbations()`
* 1.8.1 (296f9b21): minimise RAM usage when predicting targeting drugs while
using NCI60 gene expression and drug sensitivity correlation matrix

## General instructions

- Run [runRankCMapPerturbations.sh](scripts/runRankCMapPerturbations.sh)
to profile time using `Sys.time()` (no debugger attached)

- Run [runRankCMapPerturbations_heaptrack.sh](scripts/runRankCMapPerturbations_heaptrack.sh)
to profile memory with heaptrack memory profiler (timed with `Sys.time()`)

  - Convert heaptrack output to massif version (so we can plot in R) via
  [convertHeaptrackToMassif.sh](scripts/convertHeaptrackToMassif.sh)
  
  - Plot heap memory profiling with
  [R/memoryConsumptionPlot.R](R/memoryConsumptionPlot.R)

## Ranking CMap perturbations

### Input

- *User-provided data:* named numeric vector containing t-statistics of
differential expression (name corresponds to the gene symbol)
- *CMap perturbations:* publicly available differential expression z-scores;
~21GB file automatically downloaded

### CMap perturbation data loading

CMap perturbation data is first filtered according to available variables (cell 
lines, timepoints, drug dosage, perturbation types). Only the data matching the
user criteria is loaded into memory.

CMap perturbation types tested:
- *knockdown:* consensus signature from shRNAs targeting the same gene
- *overexpression:* cDNA for overexpression of wild-type gene
- *compound*

Given that the CMap perturbation data is too big for usually available RAM,
there are two options of loading CMap perturbation data:
- *On-demand (default):* load ~1GB chunks of filtered z-scores while comparing
data
- *Pre-load*: load all filtered z-scores into memory before comparing data

### Similarity ranking

CMap data is ranked against user-provided differential expression results. The
less similar the data, the higher the final rank value. Similarity is measured
using:
- Spearman's correlation coefficient
- Pearson's correlation coefficient
- GSEA-based score (weighted connectivity score as described in CMap original
article)

The values of these scores are ranked. The ranks themselves are then 
summarised via the rank product's rank (i.e. the final rank).

