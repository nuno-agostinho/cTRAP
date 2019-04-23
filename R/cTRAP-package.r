#' cTRAP package
#'
#' Compare differential gene expression results with those from big datasets
#' (e.g. CMap), allowing to infer which types of perturbations may explain the
#' observed difference in gene expression.
#'
#' \strong{Input:} To use this package, a named vector of differentially
#' expressed gene metric is needed, where its values represent the significance
#' and magnitude of the differentially expressed genes (e.g. t-statistic) and
#' its names are gene symbols.
#'
#' \strong{Workflow:} The differentially expressed genes will be compared
#' against selected perturbation conditions by:
#' \itemize{
#'     \item{Spearman or Pearson correlation with z-scores of differentially
#'     expressed genes after perturbations from CMap. Use function
#'     \code{compareAgainstCMap} with \code{method = "spearman"} or
#'     \code{method = "pearson"}}
#'     \item{Gene set enrichment analysis (GSEA) using the (around) 12 000 genes
#'     from CMap. Use function \code{compareAgainstCMap} with
#'     \code{method = gsea}.}
#' }
#'
#' Available perturbation conditions for CMap include:
#' \itemize{
#'     \item{Cell line(s).}
#'     \item{Perturbation type (gene knockdown, gene upregulation or drug
#'     intake).}
#'     \item{Drug concentration.}
#'     \item{Time points.}
#' }
#'
#' Values for each perturbation type can be listed with
#' \code{getCMapPerturbationTypes()}
#'
#' \strong{Output:} The output includes a data frame of ranked perturbations
#' based on the associated statistical values and respective p-values.
#'
#' @name cTRAP
#' @docType package
NULL

#' Gene expression data sample
#'
#' @description
#' Gene expression data sample obtained by running the following code:
#'
#' \preformatted{
#' gene <- "EIF4G1"
#' cellLine <- "HepG2"
#'
#' ENCODEmetadata <- downloadENCODEknockdownMetadata(cellLine, gene)
#' table(ENCODEmetadata$`Experiment target`)
#' length(unique(ENCODEmetadata$`Experiment target`))
#'
#' ENCODEsamples <- loadENCODEsamples(ENCODEmetadata)[[1]]
#' counts <- prepareENCODEgeneExpression(ENCODEsamples)
#'
#' # Remove low coverage (at least 10 counts shared across two samples)
#' minReads   <- 10
#' minSamples <- 2
#' filter <- rowSums(counts[ , -c(1, 2)] >= minReads) >= minSamples
#' counts <- counts[filter, ]
#' }
#'
#' @name counts
#' @docType data
#' @keywords internal
NULL

#' Differential expression's t-statistics sample
#'
#' @description
#' Differential expression's t-statistics sample obtained by running the
#' following code:
#'
#' \preformatted{
#' gene <- "EIF4G1"
#' cellLine <- "HepG2"
#'
#' ENCODEmetadata <- downloadENCODEknockdownMetadata(cellLine, gene)
#' table(ENCODEmetadata$`Experiment target`)
#' length(unique(ENCODEmetadata$`Experiment target`))
#'
#' ENCODEsamples <- loadENCODEsamples(ENCODEmetadata)[[1]]
#' counts <- prepareENCODEgeneExpression(ENCODEsamples)
#'
#' # Remove low coverage (at least 10 counts shared across two samples)
#' minReads   <- 10
#' minSamples <- 2
#' filter <- rowSums(counts[ , -c(1, 2)] >= minReads) >= minSamples
#' counts <- counts[filter, ]
#'
#' # Perform differential gene expression analysis
#' diffExpr <- performDifferentialExpression(counts)
#'
#' # Get t-statistics of differential expression with respective gene names
#' diffExprStat <- diffExpr$t
#' names(diffExprStat) <- diffExpr$Gene_symbol
#' }
#'
#' @name diffExprStat
#' @docType data
#' @keywords internal
NULL

#' ENCODE metadata sample
#'
#' @description
#' ENCODE metadata sample obtained by running the following code:
#'
#' \preformatted{
#' gene <- "EIF4G1"
#' cellLine <- "HepG2"
#' ENCODEmetadata <- downloadENCODEknockdownMetadata(cellLine, gene)
#' }
#'
#' @name ENCODEmetadata
#' @docType data
#' @keywords internal
NULL

#' CMap metadata
#'
#' @description
#' CMap metadata obtained by running the following code:
#'
#' \preformatted{
#' cmapMetadata <- loadCMapData("cmapMetadata.txt", "metadata")
#' cmapMetadata <- filterCMapMetadata(cmapMetadata, cellLine = "HEPG2",
#'                                    timepoint = "2 h")
#' }
#'
#' @name cmapMetadata
#' @docType data
#' @keywords internal
NULL

#' CMap perturbations sample for knockdown experiments
#'
#' @description
#' CMap perturbations sample for knockdown experiments obtained by running the
#' following code:
#'
#' \preformatted{
#' # Code for loading CMap gene KD HepG2 data
#' cellLine <- "HepG2"
#' cmapMetadata <- loadCMapData("cmapMetadata.txt", "metadata")
#' cmapMetadataKnockdown <- filterCMapMetadata(
#'   cmapMetadata, cellLine=cellLine,
#'   perturbationType="Consensus signature from shRNAs targeting the same gene")
#' cmapZscores  <- loadCMapData("cmapZscores.gctx", "zscores",
#'                              cmapMetadataKnockdown$sig_id)
#' cmapGeneInfo <- loadCMapData("cmapGeneInfo.txt", "geneInfo")
#'
#' cmapPerturbationsKnockdown <- loadCMapPerturbations(
#'   cmapMetadataKnockdown, cmapZscores, cmapGeneInfo)
#'
#' # Select only some perturbations (to reduce file size)
#' data("diffExprStat")
#'
#' compareKnockdown <- list()
#' compareKnockdown$spearman <- compareAgainstCMap(
#'     diffExprStat, cmapPerturbationsKnockdown, cellLine, method="spearman")
#' compareKnockdown$pearson <- compareAgainstCMap(
#'     diffExprStat, cmapPerturbationsKnockdown, cellLine, method="pearson")
#' compareKnockdown$gsea <- compareAgainstCMap(
#'     diffExprStat, cmapPerturbationsKnockdown, cellLine, method="gsea",
#'     geneSize=150)
#'
#' genes  <- lapply(compareKnockdown, "[[", "genes")
#' filter <- c(unlist(lapply(genes, head)), unlist(lapply(genes, tail)))
#' filter <- unique(filter)
#' cmapPerturbationsKnockdown <- cmapPerturbationsKnockdown[ , filter]
#' }
#'
#' @name cmapPerturbationsKnockdown
#' @docType data
#' @keywords internal
NULL

#' CMap perturbations sample for small molecules
#'
#' @description
#' CMap perturbations sample for small molecules obtained by running the
#' following code:
#'
#' \preformatted{
#' cellLine <- c("HepG2", "HUH7")
#' cmapMetadata <- loadCMapData("cmapMetadata.txt", "metadata")
#' cmapMetadataSmallMolecules <- filterCMapMetadata(
#'     cmapMetadata, cellLine=cellLine, timepoint="24 h",
#'     dosage="5 \\U00B5M", # \\U00B5 is the unicode code for the micro symbol
#'     perturbationType="Compound")
#' cmapZscores  <- loadCMapData("cmapZscores.gctx", "zscores",
#'                              cmapMetadataSmallMolecules$sig_id)
#' cmapGeneInfo <- loadCMapData("cmapGeneInfo.txt")
#'
#' cmapPerturbationsSmallMolecules <- loadCMapPerturbations(
#'     cmapMetadataSmallMolecules, cmapZscores, cmapGeneInfo)
#' }
#'
#' @name cmapPerturbationsSmallMolecules
#' @docType data
#' @keywords internal
NULL

#' Sample of ENCODE samples
#'
#' @description
#' Sample of ENCODE samples obtained by running the following code:
#'
#' \preformatted{
#' # Download ENCODE metadata for a specific cell line and gene
#' cellLine <- "HepG2"
#' gene <- "EIF4G1"
#' ENCODEmetadata <- downloadENCODEknockdownMetadata(cellLine, gene)
#'
#' # Load samples based on filtered ENCODE metadata
#' ENCODEsamples <- loadENCODEsamples(ENCODEmetadata)[[1]]
#'
#' # Get small subset of whole dataset
#' filter <- ENCODEsamples[[1]]$expected_count > 0
#' genes <- head(ENCODEsamples[[1]][filter, ])$gene_id
#' for (k in seq(ENCODEsamples)) {
#'     ENCODEsamples[[k]] <- ENCODEsamples[[k]][
#'         ENCODEsamples[[k]]$gene_id %in% genes, ]
#' }
#' }
#'
#' @name ENCODEsamples
#' @docType data
#' @keywords internal
NULL
