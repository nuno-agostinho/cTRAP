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
#'     \code{\link{rankSimilarPerturbations}} with \code{method = "spearman"} or
#'     \code{method = "pearson"}}
#'     \item{Gene set enrichment analysis (GSEA) using the (around) 12 000 genes
#'     from CMap. Use function \code{\link{rankSimilarPerturbations}} with
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

#' ENCODE metadata sample
#'
#' @description
#' ENCODE metadata sample obtained by running the following code:
#'
#' \preformatted{
#' gene <- "EIF4G1"
#' cellLine <- "HepG2"
#' ENCODEmetadata <- downloadENCODEknockdownMetadata(cellLine, gene)
#'
#' table(ENCODEmetadata$`Experiment target`)
#' length(unique(ENCODEmetadata$`Experiment target`))
#' }
#'
#' @name ENCODEmetadata
#' @docType data
#' @keywords internal
NULL

#' Gene expression data sample
#'
#' @description
#' Gene expression data sample obtained by running the following code:
#'
#' \preformatted{
#' data("ENCODEmetadata")
#' ENCODEsamples <- loadENCODEsamples(ENCODEmetadata)[[1]]
#' counts <- prepareENCODEgeneExpression(ENCODEsamples)
#'
#' # Remove low coverage (at least 10 counts shared across two samples)
#' minReads   <- 10
#' minSamples <- 2
#' filter <- rowSums(counts[ , -c(1, 2)] >= minReads) >= minSamples
#' counts <- counts[filter, ]
#'
#' # Convert ENSEMBL identifier to gene symbol
#' counts$gene_id <- convertGeneIdentifiers(counts$gene_id)
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
#' data("counts")
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

#' CMap metadata
#'
#' @description
#' CMap metadata obtained by running the following code:
#'
#' \preformatted{
#' cmapMetadata <- filterCMapMetadata("cmapMetadata.txt", cellLine = "HEPG2",
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
#' cmapMetadataKD <- filterCMapMetadata(
#'   "cmapMetadata.txt", cellLine=cellLine,
#'   perturbationType="Consensus signature from shRNAs targeting the same gene")
#'
#' cmapPerturbationsKD <- prepareCMapPerturbations(
#'   cmapMetadataKD, "cmapZscores.gctx", "cmapGeneInfo.txt",
#'   loadZscores=TRUE)
#'
#' data("diffExprStat")
#' compareKD <- rankSimilarPerturbations(diffExprStat, cmapPerturbationsKD)
#'
#' # Select only some perturbations (to reduce file size)
#' filter <- c(head(order(compareKD$spearman_rank)),
#'             tail(order(compareKD$spearman_rank)),
#'             head(order(compareKD$pearson_rank)),
#'             tail(order(compareKD$pearson_rank)),
#'             head(order(compareKD$gsea_rank)),
#'             tail(order(compareKD$gsea_rank)))
#' filter <- unique(compareKD[[1]][filter])
#' cmapPerturbationsKD <- cmapPerturbationsKD[ , filter]
#'
#' # Remove non-ASCII characters for portability reasons
#' metadata <- attr(cmapPerturbationsKD, "metadata")
#' metadata$pert_idose <- gsub("\u00B5", "micro", metadata$pert_idose)
#' metadata$pert_dose_unit <- gsub("\u00B5", "micro", metadata$pert_dose_unit)
#' attr(cmapPerturbationsKD, "metadata") <- metadata
#' }
#'
#' @name cmapPerturbationsKD
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
#' cmapMetadataCompounds <- filterCMapMetadata(
#'     "cmapMetadata.txt", cellLine=cellLine, timepoint="24 h",
#'     dosage="5 \u00B5M", perturbationType="Compound")
#'
#' cmapPerturbationsCompounds <- prepareCMapPerturbations(
#'     cmapMetadataCompounds, "cmapZscores.gctx", "cmapGeneInfo.txt",
#'     "cmapCompoundInfo_drugs.txt", loadZscores=TRUE)
#'
#' # Remove non-ASCII characters for portability reasons
#' metadata <- attr(cmapPerturbationsCompounds, "metadata")
#' metadata$pert_idose <- gsub("\u00B5", "micro", metadata$pert_idose)
#' metadata$pert_dose_unit <- gsub("\u00B5", "micro", metadata$pert_dose_unit)
#' attr(cmapPerturbationsCompounds, "metadata") <- metadata
#' }
#'
#' @name cmapPerturbationsCompounds
#' @docType data
#' @keywords internal
NULL
